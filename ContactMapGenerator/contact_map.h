#ifndef CONTACT_MAP_H
#define CONTACT_MAP_H

#include <stddef.h>
#include <stdbool.h>

/* simple 3D float */
typedef struct {
    float x, y, z;
} float3;

/* atom record used across PDB and mmCIF */
typedef struct {
    char  field[7];      /* "ATOM" or "HETATM" */
    char  name[5];       /* atom name (trimmed) */
    char  resName[5];    /* residue name (3 letters) */
    char  altLoc[2];     /* alternate location id or space */
    char  chainID[8];    /* chain identifier, may be >1 char in mmCIF */
    char  iCode[2];      /* insertion code or space */
    char  element[3];    /* element symbol */
    char  charge[3];     /* formal charge string */
    int   serial;        /* atom serial (optional) */
    int   resSeq;        /* residue sequence number */
    float x, y, z;       /* coordinates */
    float occupancy;     /* occupancy */
    float tempFactor;    /* B factor */
    int   model;         /* model number */
} atom_pdb_str;

/* flexible array container for all atoms */
typedef struct {
    int natoms;
    int nresidues;
    atom_pdb_str atoms[1]; /* flexible array */
} pdb_str;

/* residue summary */
typedef struct {
    int n;               /* number of atoms in residue */
    float x, y, z;       /* center of mass */
    float rad;           /* residue radius */
    char name[4];        /* residue 3-letter name */
    int  seq;            /* residue number */
    char chain[8];       /* full chain id */
    int  model;          /* model number */
    atom_pdb_str *a;     /* pointer to first atom of this residue */
} residue_str;

/* per-atom auxiliary (vdW etc.) */
typedef struct {
    int   atype;
    int   nb;
    float vrad;
    int   keyresidue;
    int   keyatom;
} atomaux_str;

/* surface sampling point */
typedef struct {
    float x, y, z;
    float d;
    int   i, j;
} surface_str;

/* protein map and readers */
extern bool     protein_map(atom_pdb_str *atom, atomaux_str *vdw);
extern void     read_pdb(char *name, pdb_str *pdb);

/* allocating readers */
extern pdb_str* read_pdb_alloc(const char *name);
extern pdb_str* read_cif_alloc(const char *name);
extern pdb_str* read_structure_alloc(const char *name);

/* chemistry helpers implemented elsewhere */
extern float ATOMTYPE(unsigned int i, unsigned int j);
extern int   BONDTYPE(int i, int j);

#endif /* CONTACT_MAP_H */

