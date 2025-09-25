// pdb_alloc.c
// Wrapper de asignación para leer PDB usando tu read_pdb existente.

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "contact_map.h"

static int count_pdb_atoms(const char *path) {
    FILE *f = fopen(path, "r");
    if (!f) return 0;
    char line[256];
    int n = 0;
    while (fgets(line, sizeof(line), f)) {
        // Cuenta ATOM y HETATM
        if (!strncmp(line, "ATOM  ", 6) || !strncmp(line, "HETATM", 6)) n++;
    }
    fclose(f);
    return n;
}

static pdb_str* alloc_pdb_block(int natoms) {
    if (natoms <= 0) return NULL;
    size_t bytes = sizeof(pdb_str) + (size_t)(natoms - 1) * sizeof(atom_pdb_str);
    pdb_str *p = (pdb_str*)malloc(bytes);
    if (p) memset(p, 0, bytes);
    return p;
}

pdb_str* read_pdb_alloc(const char *name) {
    int nmax = count_pdb_atoms(name);
    if (nmax <= 0) return NULL;

    // reserva un bloque suficientemente grande
    pdb_str *blk = alloc_pdb_block(nmax);
    if (!blk) return NULL;

    // tu lector existente llena el bloque y establece natoms y nresidues
    read_pdb((char*)name, blk);

    // valida
    if (blk->natoms <= 0 || blk->natoms > nmax) {
        free(blk);
        return NULL;
    }

    // reduce al tamaño exacto
    int n = blk->natoms;
    size_t bytes = sizeof(pdb_str) + (size_t)(n - 1) * sizeof(atom_pdb_str);
    pdb_str *shrunk = (pdb_str*)realloc(blk, bytes);
    return shrunk ? shrunk : blk; // si realloc falla, devuelve el bloque grande
}

