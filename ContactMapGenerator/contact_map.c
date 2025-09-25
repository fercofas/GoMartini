// contact_map.c  — sparse counters version for large assemblies

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <float.h>
#include <assert.h>
#include <math.h>
#include <limits.h>
#include <stdint.h>

#include "contact_map.h"

// GLOBALS
pdb_str        *pdb      = NULL;
atomaux_str    *atomauxs = NULL;
residue_str    *residues = NULL;
surface_str    *surface  = NULL;

// ---------- small helpers ----------

float dist(float3 c1,float3 c2) {
    float dx = c2.x - c1.x, dy = c2.y - c1.y, dz = c2.z - c1.z;
    return sqrtf(dx*dx + dy*dy + dz*dz);
}

static void make_surface(surface_str *surf, float3 c,int fiba,int fibb,float vrad) {
    int phi_aux = 0;
    for(int k=0; k<fibb; k++){
        surface_str * s = surf + k;

        s->i = -1; s->j = -1; s->d = FLT_MAX;

        phi_aux += fiba;
        if(phi_aux > fibb) phi_aux -= fibb;

        float theta = acosf(1.0f - 2.0f*k/(float)fibb);
        float phi   = 2.0f*(float)M_PI*phi_aux/(float)fibb;
        s->x = c.x + vrad*sinf(theta)*cosf(phi);
        s->y = c.y + vrad*sinf(theta)*sinf(phi);
        s->z = c.z + vrad*cosf(theta);
    }
}

// ---------- getters over your structs ----------

float3 SURFACE_COORDINATES(int k) {
    surface_str *s = surface + k;
    float3 r = {s->x,s->y,s->z};
    return r;
}
float3 RESIDUE_CENTER_OF_MASS(int i) {
    residue_str *res = residues+i;
    float3 rc = {res->x,res->y,res->z};
    return rc;
}
float RESIDUE_RADIUS(int i) {
    return residues[i].rad;
}
float3 ATOM_COORDINATES(int i,int j) {
    residue_str *res = residues+i;
    atom_pdb_str *a = res->a + j;
    float3 c = {a->x,a->y,a->z};
    return c;
}
int ATOM_INDEX(int i1,int j1) {
    residue_str * res1 = residues+i1;
    atom_pdb_str * atom1 = res1->a + j1;
    return (int)(atom1 - pdb->atoms);
}
float vdW_RADIUS(int i1,int j1) {
    int k = ATOM_INDEX(i1,j1);
    return atomauxs[k].vrad;
}
int BASES_IN_RESIDUE(int i1) {
    return residues[i1].n;
}

static float DISTANCE_C_ALPHA(int i1,int i2) {
    // Legacy heuristic: CA is the second atom. Works for standard protein PDB/mmCIF.
    int j1=1,j2=1;
    return dist(ATOM_COORDINATES(i1,j1),ATOM_COORDINATES(i2,j2));
}

// ---------- build residues and aux ----------

void init_atomauxs_and_residues() {
    int   resCount  = -1;
    int   lastResSeq = INT_MIN;
    int   lastModel  = INT_MIN;
    char  lastChain[8] = {0};
    char  lastICode = '\0';

    for (int k=0;k<pdb->natoms;k++) {
        atom_pdb_str * atom = pdb->atoms + k;

        bool newResidue =
            (resCount < 0) ||
            (atom->model != lastModel) ||
            (strcmp(atom->chainID, lastChain) != 0) ||
            (atom->resSeq != lastResSeq) ||
            (atom->iCode[0] != lastICode);

        if (newResidue) {
            lastModel  = atom->model;
            lastResSeq = atom->resSeq;
            strncpy(lastChain, atom->chainID, sizeof(lastChain)-1);
            lastICode = atom->iCode[0];

            resCount++;
            residue_str * res = residues + resCount;
            memset(res, 0, sizeof(*res));
            res->a = atom;
            strncpy(res->name,  atom->resName, sizeof(res->name)-1);
            res->seq   = atom->resSeq;
            strncpy(res->chain, atom->chainID, sizeof(res->chain)-1);
            res->model = atom->model;
        }

        atomaux_str * aux = atomauxs + k;
        protein_map(atom,aux);

        residue_str * res = residues + resCount;
        res->n++;
        res->x += atom->x;
        res->y += atom->y;
        res->z += atom->z;
    }
    assert(resCount+1 == pdb->nresidues);

    for (int k=0; k<pdb->nresidues; k++) {
        residue_str * res = residues + k;
        res->x /= res->n;
        res->y /= res->n;
        res->z /= res->n;
    }
}

// ---------- sparse counters for (i1,i2) ----------

typedef struct {
    uint64_t key;  // ((uint64_t)i1 << 32) | (uint64_t)i2
    int over;
    int cont;
    int stab;
    int dest;
} PairCount;

typedef struct {
    size_t cap;   // power of two
    size_t used;
    PairCount *a; // open addressing
} CounterMap;

static const uint64_t EMPTY_KEY = UINT64_MAX;

static uint64_t mix64(uint64_t x){
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
}
static size_t next_pow2(size_t v){
    if (v < 2) return 2;
    v--;
    for (size_t i=1;i<sizeof(size_t)*8;i<<=1) v |= v>>i;
    return v+1;
}
static void cmap_init(CounterMap *m, size_t hint){
    m->cap = next_pow2(hint*2 + 1024);
    m->used = 0;
    m->a = (PairCount*)malloc(m->cap * sizeof(PairCount));
    if (!m->a) { m->cap=0; return; }
    for (size_t i=0;i<m->cap;i++){ m->a[i].key = EMPTY_KEY; m->a[i].over=m->a[i].cont=m->a[i].stab=m->a[i].dest=0; }
}
static void cmap_free(CounterMap *m){
    free(m->a); m->a=NULL; m->cap=m->used=0;
}
static PairCount *cmap_find_slot(CounterMap *m, uint64_t key){
    size_t mask = m->cap - 1;
    size_t i = (size_t)(mix64(key) & mask);
    for (;;){
        if (m->a[i].key == EMPTY_KEY) return &m->a[i];
        if (m->a[i].key == key) return &m->a[i];
        i = (i + 1) & mask;
    }
}
static bool cmap_rehash(CounterMap *m, size_t newcap){
    PairCount *old = m->a;
    size_t oldcap = m->cap;
    m->a = (PairCount*)malloc(newcap * sizeof(PairCount));
    if (!m->a) { m->a = old; return false; }
    m->cap = newcap; m->used = 0;
    for (size_t i=0;i<m->cap;i++){ m->a[i].key = EMPTY_KEY; m->a[i].over=m->a[i].cont=m->a[i].stab=m->a[i].dest=0; }
    for (size_t i=0;i<oldcap;i++){
        if (old[i].key != EMPTY_KEY){
            PairCount *s = cmap_find_slot(m, old[i].key);
            *s = old[i];
            m->used++;
        }
    }
    free(old);
    return true;
}
static PairCount *cmap_get(CounterMap *m, int i1, int i2){
    uint64_t key = ((uint64_t)(uint32_t)i1 << 32) | (uint64_t)(uint32_t)i2;
    if (m->cap == 0) return NULL;
    PairCount *s = cmap_find_slot(m, key);
    if (s->key == EMPTY_KEY){
        s->key = key;
        s->over = s->cont = s->stab = s->dest = 0;
        m->used++;
        if (m->used * 10 > m->cap * 7) {
            cmap_rehash(m, m->cap ? m->cap*2 : 2048);
        }
    }
    return s;
}

// ---------- main ----------

int main (int argc,char **argv) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <pdb|cif>\n", argv[0]);
        return 1;
    }

    printf(
"                         CONTACT MAPS FROM PDB FILES                          \n"
"\n"
" This software is written by:\n"
"       Rodrigo Azevedo Moreira da Silva\n"
"\n"
" Copyright (c) 2020 - IPPT-PAN\n"
"       Institute of Fundamental Techonological Research\n"
"       Polish Academy of Sciences\n"
" MIT LICENSE, check out LICENSE for more informations.\n"
"\n"
    );

    printf("Reading file:    %s\n",argv[1]);

    // read structure into heap block
    pdb = read_structure_alloc(argv[1]);
    if (!pdb) {
        fprintf(stderr, "Could not read structure: %s\n", argv[1]);
        return 1;
    }
    int natoms    = pdb->natoms;
    int nresidues = pdb->nresidues;
    printf("pdb natoms:      %d\n",natoms);
    printf("pdb nresidues:   %d\n",nresidues);

    // AUX ARRAYS
    atomauxs = (atomaux_str*)calloc((size_t)natoms,sizeof(atomaux_str));
    residues = (residue_str*)calloc((size_t)nresidues,sizeof(residue_str));
    if (!atomauxs || !residues) {
        fprintf(stderr, "Out of memory for auxiliaries\n");
        return 2;
    }
    init_atomauxs_and_residues();

    // Fibonacci grid for surface sampling
    int fib = 14, fiba = 0, fibb = 1;
    for (int f=0;f<fib;f++) { int fibc=fiba+fibb; fiba=fibb; fibb=fibc; }

    surface = (surface_str*)calloc((size_t)fibb,sizeof(surface_str));
    if (!surface) { fprintf(stderr,"Out of memory for surface\n"); return 2; }
    printf("Fibonacci grid:  %d\n",fibb);

    // MODEL PARAMETERS
    float alpha        = 1.24f;
    float water_radius = 2.80f;
    printf("ALPHA:        %7.2f\n",alpha);
    printf("WATER_RADIUS: %7.2f\n",water_radius);

    // Sparse counters
    CounterMap cmap;
    // heuristic: ~64 neighbors per residue
    size_t hint_pairs = (size_t)nresidues * 64u;
    cmap_init(&cmap, hint_pairs);
    if (!cmap.a) { fprintf(stderr, "Out of memory for counters\n"); return 2; }

    // CONTACTS
    for(int i1=0; i1 < nresidues; i1++) {
        for(int j1=0; j1 < BASES_IN_RESIDUE(i1); j1++) {
            make_surface(surface, ATOM_COORDINATES(i1,j1), fiba, fibb, vdW_RADIUS(i1,j1)+water_radius);

            for(int i2=0; i2 < nresidues; i2++) {
                if (residues[i1].model != residues[i2].model) continue;
                // coarse filter on residue COM distance
                if (dist(RESIDUE_CENTER_OF_MASS(i1),RESIDUE_CENTER_OF_MASS(i2)) > 14.0f) continue;

                for(int j2=0; j2 < BASES_IN_RESIDUE(i2); j2++) {
                    if (i1==i2 && j1==j2) continue;

                    float distance = dist(ATOM_COORDINATES(i1,j1),ATOM_COORDINATES(i2,j2));

                    PairCount *pc = NULL;

                    // OV criterion
                    if(distance <= ((vdW_RADIUS(i1,j1)+vdW_RADIUS(i2,j2))*alpha)) {
                        if (!pc) pc = cmap_get(&cmap, i1, i2);
                        if (pc) pc->over = 1;
                    }

                    // CSU coverage, plus stabilize/destabilize classification
                    if (distance <= vdW_RADIUS(i1,j1)+vdW_RADIUS(i2,j2)+water_radius) {
                        for (int k=0; k<fibb; k++) {
                            surface_str * s = surface + k;
                            if( dist(SURFACE_COORDINATES(k), ATOM_COORDINATES(i2,j2)) < vdW_RADIUS(i2,j2)+water_radius
                                && distance <= s->d) {
                                s->d = distance; s->i = i2; s->j = j2;
                            }
                        }
                    }

                    // classify bond type if this atom pair gets chosen by any surface point later
                    // we count after the surface loop below
                    (void)pc;
                }
            }

            // accumulate surface hits for this i1,j1
            for(int k=0; k<fibb; k++) {
                surface_str * s = surface + k;
                int i2 = s->i;
                int j2 = s->j;
                if (i2 >= 0 && j2 >= 0) {
                    int at1 = ATOMTYPE(i1,j1);
                    int at2 = ATOMTYPE(i2,j2);
                    if (at1>0 && at2>0) {
                        PairCount *pc = cmap_get(&cmap, i1, i2);
                        if (pc) {
                            pc->cont += 1;
                            int btype = BONDTYPE(at1,at2);
                            if (btype <= 4) pc->stab += 1;
                            if (btype == 5) pc->dest += 1;
                        }
                    }
                }
            }

        } // j1
    }     // i1

    // dump table header
    printf("\n"
"Residue-Residue Contacts\n"
"\n"
"ID       - atom identification\n"
"I1,I2    - serial residue id\n"
"AA       - 3-letter code of aminoacid\n"
"C        - chain\n"
"I(PDB)   - residue number in PDB file\n"
"DCA      - distance between CA\n"
"CMs      - OV , CSU , oCSU , rCSU\n"
"           (CSU does not take into account chemical properties of atoms)\n"
"rCSU     - net contact from rCSU\n"
"Count    - number of contacts between residues\n"
"MODEL    - model number\n"
"\n"
"      ID    I1  AA  C I(PDB)     I2  AA  C I(PDB)        DCA       CMs    rCSU   Count Model\n"
"============================================================================================\n"
    );

    // collect entries
    size_t nitems = 0;
    for (size_t i=0;i<cmap.cap;i++) if (cmap.a[i].key != EMPTY_KEY) nitems++;
    typedef struct { int i1,i2; PairCount *pc; } RowRef;
    RowRef *rows = (RowRef*)malloc(nitems * sizeof(RowRef));
    size_t rpos=0;
    for (size_t i=0;i<cmap.cap;i++){
        if (cmap.a[i].key != EMPTY_KEY){
            int i1 = (int)(cmap.a[i].key >> 32);
            int i2 = (int)(cmap.a[i].key & 0xffffffffu);
            rows[rpos++] = (RowRef){i1,i2,&cmap.a[i]};
        }
    }
    // sort by i1 then i2 for stable output
    int cmp_rowref(const void *a, const void *b){
        const RowRef *ra = (const RowRef*)a, *rb = (const RowRef*)b;
        if (ra->i1 != rb->i1) return (ra->i1 < rb->i1)? -1 : 1;
        if (ra->i2 != rb->i2) return (ra->i2 < rb->i2)? -1 : 1;
        return 0;
    }
    qsort(rows, nitems, sizeof(RowRef), cmp_rowref);

    int count = 0;
    for (size_t t=0; t<nitems; t++) {
        int i1 = rows[t].i1;
        int i2 = rows[t].i2;
        PairCount *pc = rows[t].pc;

        if (i1==i2) continue; // skip self
        if (residues[i1].model != residues[i2].model) continue;

        int over = pc->over;
        int cont = pc->cont;
        int stab = pc->stab;
        int dest = pc->dest;
        int ocsu = stab;
        int rcsu = stab - dest;

        if (over > 0 || cont > 0) {
            count++;
            printf("R %6d ",count);
            printf("%5d %4s %s %4d",i1+1,residues[i1].name,residues[i1].chain,residues[i1].seq);
            printf("    ");
            printf("%5d %4s %s %4d",i2+1,residues[i2].name,residues[i2].chain,residues[i2].seq);
            printf("     %8.4f     ",DISTANCE_C_ALPHA(i1,i2));
            printf("%d %d %d %d", over, cont != 0 ? 1 : 0, ocsu != 0 ? 1 : 0, rcsu >  0 ? 1 : 0);
            printf("%6d  %6d %4d\n", rcsu, cont, residues[i1].model);
        }
    }

    // cleanup
    free(rows);
    cmap_free(&cmap);
    free(pdb);
    free(atomauxs);
    free(residues);
    free(surface);

    return 0;
}

