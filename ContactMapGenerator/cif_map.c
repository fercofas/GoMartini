// cif_map.c
// Robust mmCIF reader that fills a heap-allocated pdb_str using the atom_site loop.
// Exposes read_cif_alloc, read_pdb_alloc (wrapper on your read_pdb), and read_structure_alloc.

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <stdbool.h>

#include "contact_map.h"

#if !defined(_WIN32)
#include <strings.h>  // strcasecmp
#endif

#ifndef CIF_DEBUG
#define CIF_DEBUG 0   // set to 0 when you finish debugging
#endif

// --------------------- small helpers ---------------------

static void rstrip_newline(char *s) {
    size_t n = strlen(s);
    while (n && (s[n-1] == '\n' || s[n-1] == '\r')) s[--n] = '\0';
}

// trim all trailing whitespace (spaces, tabs, CR/LF)
static void rstrip_ws(char *s) {
    size_t n = strlen(s);
    while (n && isspace((unsigned char)s[n-1])) s[--n] = '\0';
}

static char *lstrip(char *s) {
    while (*s && isspace((unsigned char)*s)) s++;
    return s;
}

static int ieq(const char *a, const char *b) {
#ifdef _WIN32
    return _stricmp(a, b) == 0;
#else
    return strcasecmp(a, b) == 0;
#endif
}

static int has_ext(const char *path, const char *ext) {
    const char *dot = strrchr(path, '.');
    return dot && ieq(dot + 1, ext);
}

static int is_missing(const char *s) {
    return (s == NULL || s[0] == '\0' || strcmp(s, ".") == 0 || strcmp(s, "?") == 0);
}

static void copy_trim(char *dst, size_t dstsz, const char *src, size_t maxkeep) {
    size_t j = 0;
    for (size_t i = 0; src && src[i] && j < maxkeep; ++i) {
        if (!isspace((unsigned char)src[i])) dst[j++] = src[i];
    }
    if (dstsz == 0) return;
    if (j >= dstsz) j = dstsz - 1;
    dst[j] = '\0';
}

static void upper_copy(char *dst, size_t dstsz, const char *src, size_t maxkeep) {
    size_t j = 0;
    for (size_t i = 0; src && src[i] && j < maxkeep; ++i) {
        if (!isspace((unsigned char)src[i])) dst[j++] = (char)toupper((unsigned char)src[i]);
    }
    if (dstsz == 0) return;
    if (j >= dstsz) j = dstsz - 1;
    dst[j] = '\0';
}

// Tokenizer for mmCIF rows with simple quoted field support.
static int split_cif_row(char *line, char **tokens, int max_tokens) {
    int nt = 0;
    char *p = line;
    while (*p) {
        while (isspace((unsigned char)*p)) p++;
        if (!*p || *p == '#') break;

        char *start = p;
        char q = 0;
        if (*p == '\'' || *p == '"') {
            q = *p++;
            start = p;
            while (*p && *p != q) p++;
        } else {
            while (*p && !isspace((unsigned char)*p)) p++;
        }
        size_t len = (size_t)(p - start);
        if (q && *p == q) p++;
        if (len) {
            if (nt >= max_tokens) return nt;
            tokens[nt] = (char*)malloc(len + 1);
            if (!tokens[nt]) return -1;
            memcpy(tokens[nt], start, len);
            tokens[nt][len] = '\0';
            nt++;
        }
    }
    return nt;
}

static void free_tokens(char **tokens, int nt) {
    for (int i = 0; i < nt; ++i) free(tokens[i]);
}

// ----------------- atom_site column bookkeeping -----------------

typedef struct {
    int i_group_PDB;

    int i_auth_asym_id,  i_label_asym_id;
    int i_auth_seq_id,   i_label_seq_id;
    int i_auth_comp_id,  i_label_comp_id;
    int i_auth_atom_id,  i_label_atom_id;
    int i_pdb_atom_name;
    int i_pdb_model;

    int i_id, i_altloc, i_inscode, i_element, i_charge;
    int i_x, i_y, i_z, i_occ, i_b;

    int ncols;
    bool ready;
} atom_site_cols;

static void cols_init(atom_site_cols *c) {
    memset(c, 0xff, sizeof(*c)); // set all to -1
    c->ncols = 0;
    c->ready = false;
}

static int has_atom_site_prefix(const char *s) {
    s = lstrip((char*)s);
#ifdef _WIN32
    return _strnicmp(s, "_atom_site.", 11) == 0;
#else
    return strncasecmp(s, "_atom_site.", 11) == 0;
#endif
}

// compare suffix (after '.') ignoring trailing spaces in header
static int find_suffix(char **hdr, int nhdr, const char *suffix) {
    size_t slen = strlen(suffix);
    for (int i = 0; i < nhdr; ++i) {
        const char *p = lstrip(hdr[i]);
        const char *dot = strchr(p, '.');
        const char *s = dot ? dot + 1 : p;
        // take s up to first whitespace
        size_t len = 0;
        while (s[len] && !isspace((unsigned char)s[len])) len++;
        if (len == slen) {
#ifdef _WIN32
            if (_strnicmp(s, suffix, (unsigned int)len) == 0) return i;
#else
            if (strncasecmp(s, suffix, len) == 0) return i;
#endif
        }
    }
    return -1;
}

static void bulk_map_headers(atom_site_cols *c, char **hdr, int nhdr) {
    int grp = find_suffix(hdr, nhdr, "group_PDB");
    if (grp >= 0) c->i_group_PDB = grp;

    c->i_label_asym_id = find_suffix(hdr, nhdr, "label_asym_id");
    c->i_auth_asym_id  = find_suffix(hdr, nhdr, "auth_asym_id");

    c->i_label_seq_id  = find_suffix(hdr, nhdr, "label_seq_id");
    c->i_auth_seq_id   = find_suffix(hdr, nhdr, "auth_seq_id");

    c->i_label_comp_id = find_suffix(hdr, nhdr, "label_comp_id");
    c->i_auth_comp_id  = find_suffix(hdr, nhdr, "auth_comp_id");

    c->i_pdb_atom_name = find_suffix(hdr, nhdr, "pdbx_PDB_atom_name");
    c->i_auth_atom_id  = find_suffix(hdr, nhdr, "auth_atom_id");
    c->i_label_atom_id = find_suffix(hdr, nhdr, "label_atom_id");

    c->i_pdb_model     = find_suffix(hdr, nhdr, "pdbx_PDB_model_num");

    c->i_id            = find_suffix(hdr, nhdr, "id");
    c->i_altloc        = find_suffix(hdr, nhdr, "label_alt_id");
    if (c->i_altloc < 0) c->i_altloc = find_suffix(hdr, nhdr, "pdbx_alt_id");
    c->i_inscode       = find_suffix(hdr, nhdr, "pdbx_PDB_ins_code");
    c->i_element       = find_suffix(hdr, nhdr, "type_symbol");
    c->i_charge        = find_suffix(hdr, nhdr, "pdbx_formal_charge");

    c->i_x             = find_suffix(hdr, nhdr, "Cartn_x");
    c->i_y             = find_suffix(hdr, nhdr, "Cartn_y");
    c->i_z             = find_suffix(hdr, nhdr, "Cartn_z");
}

static bool cols_complete(const atom_site_cols *c) {
    return (c->i_auth_asym_id  >= 0 || c->i_label_asym_id  >= 0) &&
           (c->i_auth_seq_id   >= 0 || c->i_label_seq_id   >= 0) &&
           (c->i_auth_comp_id  >= 0 || c->i_label_comp_id  >= 0) &&
           (c->i_pdb_atom_name >= 0 || c->i_auth_atom_id   >= 0 || c->i_label_atom_id >= 0) &&
           c->i_x >= 0 && c->i_y >= 0 && c->i_z >= 0;
}

#if CIF_DEBUG
static void debug_print_header(const atom_site_cols *c, char **hdr, int nhdr) {
    fprintf(stderr, "CIF atom_site header (%d cols):\n", nhdr);
    for (int i = 0; i < nhdr; ++i) fprintf(stderr, "  %s\n", hdr[i]);
    fprintf(stderr,
        "Mapped indices: group=%d label_atom=%d auth_atom=%d pdb_atom=%d "
        "label_comp=%d auth_comp=%d label_asym=%d auth_asym=%d "
        "label_seq=%d auth_seq=%d x=%d y=%d z=%d occ=%d b=%d id=%d alt=%d ins=%d elem=%d charge=%d model=%d ready=%d\n",
        c->i_group_PDB, c->i_label_atom_id, c->i_auth_atom_id, c->i_pdb_atom_name,
        c->i_label_comp_id, c->i_auth_comp_id, c->i_label_asym_id, c->i_auth_asym_id,
        c->i_label_seq_id, c->i_auth_seq_id, c->i_x, c->i_y, c->i_z, c->i_occ, c->i_b,
        c->i_id, c->i_altloc, c->i_inscode, c->i_element, c->i_charge, c->i_pdb_model,
        (int)c->ready);
}
#endif

// ----------------- atom name normalization -----------------

static void normalize_atom_name(char *name) {
    size_t L = strlen(name);
    if (L == 0) return;
    if (name[0] == 'H') return; // keep hydrogens as-is

    if (L == 1) {
        if (name[0]=='A'||name[0]=='B'||name[0]=='G'||name[0]=='D'||name[0]=='E'||name[0]=='Z'){
            char tmp[5]; tmp[0]='C'; tmp[1]=name[0]; tmp[2]='\0';
            strncpy(name, tmp, 4); name[4]='\0';
        }
        return;
    }
    if (L == 2) {
        if ((name[0]=='D'||name[0]=='E') && (name[1]=='1'||name[1]=='2')) {
            char tmp[5]; tmp[0]='C'; tmp[1]=name[0]; tmp[2]=name[1]; tmp[3]='\0';
            strncpy(name, tmp, 4); name[4]='\0';
        }
    }
}

static const char* choose3(char **tok, int i1, int i2, int i3) {
    if (i1 >= 0 && !is_missing(tok[i1])) return tok[i1];
    if (i2 >= 0 && !is_missing(tok[i2])) return tok[i2];
    if (i3 >= 0 && !is_missing(tok[i3])) return tok[i3];
    return ".";
}

// ------------------------- core fill -------------------------

static void fill_atom_from_tokens(const atom_site_cols *c, char **tok, atom_pdb_str *out, int model_id_default) {
    memset(out, 0, sizeof(*out));
    strncpy(out->field, "ATOM", sizeof(out->field)-1);

    out->serial = (c->i_id >= 0 && !is_missing(tok[c->i_id])) ? atoi(tok[c->i_id]) : 0;

    const char *atom_id = choose3(tok, c->i_pdb_atom_name, c->i_auth_atom_id, c->i_label_atom_id);
    copy_trim(out->name, sizeof(out->name), atom_id, 4);
    normalize_atom_name(out->name);

    if (c->i_altloc >= 0 && !is_missing(tok[c->i_altloc]))
        copy_trim(out->altLoc, sizeof(out->altLoc), tok[c->i_altloc], 1);
    else { out->altLoc[0] = ' '; out->altLoc[1] = '\0'; }

    const char *comp_id = choose3(tok, c->i_label_comp_id, c->i_auth_comp_id, -1);
    upper_copy(out->resName, sizeof(out->resName), comp_id, 4);

    const char *asym_id = choose3(tok, c->i_auth_asym_id, c->i_label_asym_id, -1);
    copy_trim(out->chainID, sizeof(out->chainID), asym_id, (size_t)(sizeof(out->chainID)-1));
    
    const char *seq_id = choose3(tok, c->i_auth_seq_id, c->i_label_seq_id, -1);
    out->resSeq = is_missing(seq_id) ? 0 : atoi(seq_id);

    if (c->i_inscode >= 0 && !is_missing(tok[c->i_inscode]))
        copy_trim(out->iCode, sizeof(out->iCode), tok[c->i_inscode], 1);
    else { out->iCode[0] = ' '; out->iCode[1] = '\0'; }

    if (c->i_element >= 0 && !is_missing(tok[c->i_element]))
        upper_copy(out->element, sizeof(out->element), tok[c->i_element], 2);
    if (c->i_charge >= 0 && !is_missing(tok[c->i_charge]))
        copy_trim(out->charge, sizeof(out->charge), tok[c->i_charge], 2);

    out->x = (c->i_x >= 0 && !is_missing(tok[c->i_x])) ? (float)atof(tok[c->i_x]) : 0.0f;
    out->y = (c->i_y >= 0 && !is_missing(tok[c->i_y])) ? (float)atof(tok[c->i_y]) : 0.0f;
    out->z = (c->i_z >= 0 && !is_missing(tok[c->i_z])) ? (float)atof(tok[c->i_z]) : 0.0f;
    out->occupancy  = (c->i_occ >= 0 && !is_missing(tok[c->i_occ])) ? (float)atof(tok[c->i_occ]) : 1.0f;
    out->tempFactor = (c->i_b   >= 0 && !is_missing(tok[c->i_b]))   ? (float)atof(tok[c->i_b])   : 0.0f;

    if (c->i_pdb_model >= 0 && !is_missing(tok[c->i_pdb_model]))
        out->model = atoi(tok[c->i_pdb_model]);
    else
        out->model = model_id_default;
}

// --------------------- alloc helpers ---------------------

static pdb_str *alloc_pdb_block(int natoms) {
    if (natoms <= 0) return NULL;
    size_t bytes = sizeof(pdb_str) + (size_t)(natoms - 1) * sizeof(atom_pdb_str);
    pdb_str *p = (pdb_str*)malloc(bytes);
    if (p) memset(p, 0, bytes);
    return p;
}

// ------------------- CIF allocating reader -------------------

pdb_str* read_cif_alloc(const char *name) {
    FILE *fp = fopen(name, "r");
    if (!fp) {
        fprintf(stderr, "Could not open CIF file: %s\n", name);
        return NULL;
    }

    char line[8192];
    bool in_atom_loop = false;
    atom_site_cols cols; cols_init(&cols);
    int model_id = 1;

    size_t cap = 0, n = 0;
    atom_pdb_str *arr = NULL;

    enum { MAX_HDR = 1024 };
    char *hdr[MAX_HDR]; int nhdr = 0;

#if CIF_DEBUG
    size_t dbg_rows = 0;
#endif

    while (fgets(line, sizeof(line), fp)) {
        rstrip_newline(line);
        char *sline = lstrip(line);
        if (*sline == '#') continue;

        if (!in_atom_loop) {
            if (strcmp(sline, "loop_") == 0) {
                long pos = ftell(fp);
                nhdr = 0;
                while (fgets(line, sizeof(line), fp)) {
                    rstrip_newline(line);
                    rstrip_ws(line);               // trim trailing spaces in header lines
                    char *p = lstrip(line);
                    if (*p == '_') {
                        if (nhdr < MAX_HDR) hdr[nhdr++] = strdup(p);
                    } else break;
                }
                bool is_atom_site = false;
                for (int i = 0; i < nhdr; ++i) {
                    if (has_atom_site_prefix(hdr[i])) { is_atom_site = true; break; }
                }
                if (!is_atom_site) {
                    for (int i = 0; i < nhdr; ++i) free(hdr[i]);
                    nhdr = 0;
                    fseek(fp, pos, SEEK_SET);
                    continue;
                }

                cols_init(&cols);
                bulk_map_headers(&cols, hdr, nhdr);
                cols.ncols = nhdr;

                // positional fallback for common PDBx layout
                if (cols.i_group_PDB == 0 && cols.ncols >= 18) {
                    if (cols.i_id            < 0) cols.i_id            = 1;
                    if (cols.i_element       < 0) cols.i_element       = 2;
                    if (cols.i_label_atom_id < 0) cols.i_label_atom_id = 3;
                    if (cols.i_altloc        < 0) cols.i_altloc        = 4;
                    if (cols.i_label_comp_id < 0) cols.i_label_comp_id = 5;
                    if (cols.i_label_asym_id < 0) cols.i_label_asym_id = 6;
                    if (cols.i_label_seq_id  < 0) cols.i_label_seq_id  = 7;
                    if (cols.i_inscode       < 0) cols.i_inscode       = 8;
                    if (cols.i_x             < 0) cols.i_x             = 9;
                    if (cols.i_y             < 0) cols.i_y             = 10;
                    if (cols.i_z             < 0) cols.i_z             = 11;
                    if (cols.i_auth_asym_id  < 0) cols.i_auth_asym_id  = 12;
                    if (cols.i_auth_seq_id   < 0) cols.i_auth_seq_id   = 13;
                    if (cols.i_pdb_model     < 0) cols.i_pdb_model     = 14;
                    if (cols.i_occ           < 0) cols.i_occ           = 15;
                    if (cols.i_b             < 0) cols.i_b             = 16;
                    if (cols.i_charge        < 0) cols.i_charge        = 17;
                }

                cols.ready = cols_complete(&cols);

#if CIF_DEBUG
                debug_print_header(&cols, hdr, nhdr);
#endif

                for (int i = 0; i < nhdr; ++i) free(hdr[i]);
                nhdr = 0;

                if (!cols.ready) {
                    fprintf(stderr, "CIF atom_site loop is missing required columns in %s\n", name);
                    continue;
                }

                in_atom_loop = true;
                // fall through with current line as first data row
            } else {
                continue;
            }
        }

        // data rows
        if (*sline == '\0') continue;
        if (*sline == '_' || strcmp(sline, "loop_") == 0 || strncmp(sline, "data_", 5) == 0) {
            in_atom_loop = false;
            fseek(fp, -(long)strlen(line) - 1, SEEK_CUR); // push back this line (approx)
            continue;
        }

        enum { MAXTOK = 2048 };
        char *tok[MAXTOK] = {0};

        char buf[8192];
        strncpy(buf, sline, sizeof(buf) - 1);
        buf[sizeof(buf) - 1] = '\0';
        int nt = split_cif_row(buf, tok, MAXTOK);
        if (nt < 0) { fclose(fp); free(arr); return NULL; }

        // try continuation if short
        while (nt < cols.ncols) {
            long save = ftell(fp);
            if (!fgets(line, sizeof(line), fp)) break;
            rstrip_newline(line);
            char *cont = lstrip(line);
            if (*cont == '_' || strcmp(cont, "loop_") == 0 || strncmp(cont, "data_", 5) == 0 || *cont == '#') {
                fseek(fp, save, SEEK_SET);
                break;
            }
            size_t cur = strlen(buf), add = strlen(cont);
            if (cur + 1 + add + 1 < sizeof(buf)) {
                buf[cur] = ' ';
                memcpy(buf + cur + 1, cont, add + 1);
                free_tokens(tok, nt);
                nt = split_cif_row(buf, tok, MAXTOK);
                if (nt < 0) break;
            } else break;
        }

#if CIF_DEBUG
        if (dbg_rows < 4 && nt > 0) {
            int ia = (cols.i_pdb_atom_name>=0? cols.i_pdb_atom_name:
                     (cols.i_auth_atom_id >=0? cols.i_auth_atom_id:
                     (cols.i_label_atom_id>=0? cols.i_label_atom_id: -1)));
            int ic = (cols.i_label_comp_id>=0? cols.i_label_comp_id:
                     (cols.i_auth_comp_id >=0? cols.i_auth_comp_id: -1));
            int ich = (cols.i_auth_asym_id>=0? cols.i_auth_asym_id:
                      (cols.i_label_asym_id>=0? cols.i_label_asym_id: -1));
            int isq = (cols.i_auth_seq_id>=0? cols.i_auth_seq_id:
                      (cols.i_label_seq_id>=0? cols.i_label_seq_id: -1));
            fprintf(stderr, "ROW[%zu] atom=%s res=%s chain=%s seq=%s\n", dbg_rows,
                (ia>=0? tok[ia]: "."), (ic>=0? tok[ic]: "."),
                (ich>=0? tok[ich]: "."), (isq>=0? tok[isq]: "."));
            dbg_rows++;
        }
#endif

        if (nt >= cols.ncols) {
            if (n == cap) {
                cap = cap ? cap * 2 : 4096;
                arr = (atom_pdb_str*)realloc(arr, cap * sizeof(*arr));
                if (!arr) { fclose(fp); return NULL; }
            }
            fill_atom_from_tokens(&cols, tok, &arr[n], model_id);
            n++;
        }
        free_tokens(tok, nt);
    }

    fclose(fp);

    if (n == 0) {
        fprintf(stderr, "No atoms read from CIF file: %s\n", name);
        free(arr);
        return NULL;
    }

    // pack
    pdb_str *out = alloc_pdb_block((int)n);
    if (!out) { free(arr); return NULL; }
    out->natoms = (int)n;

    // residue count by transitions (model, chain, resSeq, iCode)
    int nres = 0;
    if (n > 0) {
        nres = 1;
        for (size_t i = 1; i < n; ++i) {
            if (arr[i].model != arr[i-1].model ||
                arr[i].chainID[0] != arr[i-1].chainID[0] ||
                arr[i].resSeq     != arr[i-1].resSeq     ||
                arr[i].iCode[0]   != arr[i-1].iCode[0]) {
                nres++;
            }
        }
    }
    out->nresidues = nres;

    memcpy(out->atoms, arr, n * sizeof(atom_pdb_str));
    free(arr);
    return out;
}

// --------------- PDB allocating wrapper on read_pdb ---------------

static int count_pdb_atoms_(const char *path) {
    FILE *f = fopen(path, "r");
    if (!f) return 0;
    char line[256];
    int n = 0;
    while (fgets(line, sizeof(line), f)) {
        if (!strncmp(line, "ATOM  ", 6) || !strncmp(line, "HETATM", 6)) n++;
    }
    fclose(f);
    return n;
}

pdb_str* read_pdb_alloc(const char *name) {
    int nmax = count_pdb_atoms_(name);
    if (nmax <= 0) return NULL;

    pdb_str *blk = alloc_pdb_block(nmax);
    if (!blk) return NULL;

    read_pdb((char*)name, blk); // legacy non-allocating reader

    if (blk->natoms <= 0 || blk->natoms > nmax) {
        free(blk);
        return NULL;
    }
    int n = blk->natoms;
    size_t bytes = sizeof(pdb_str) + (size_t)(n - 1) * sizeof(atom_pdb_str);
    pdb_str *shrunk = (pdb_str*)realloc(blk, bytes);
    return shrunk ? shrunk : blk;
}

// ---------------------- autodetect reader ----------------------

pdb_str* read_structure_alloc(const char *name) {
    if (has_ext(name, "cif") || has_ext(name, "mmcif"))
        return read_cif_alloc(name);
    return read_pdb_alloc(name);
}

