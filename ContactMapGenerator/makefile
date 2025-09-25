# Portable makefile
# Usage:
#   make                # build
#   make install        # symlink to ~/bin/contact_map
#   make clean          # remove build artifacts

# Compiler and flags (override from CLI if you want)
CC      ?= gcc
CFLAGS  ?= -Wall -O3 -march=native
LDFLAGS ?=
LDLIBS  ?= -lm

# Install prefix
PREFIX  ?= $(HOME)
BINDIR  := $(PREFIX)/bin

DEPS := contact_map.h
OBJ = chemical_map.o pdb_map.o protein_map.o contact_map.o cif_map.o


contact_map: $(OBJ)
	$(CC) $(LDFLAGS) -o $@ $^ $(LDLIBS)

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

.PHONY: clean install

clean:
	rm -f *.o contact_map

install: contact_map
	mkdir -p "$(BINDIR)"
	ln -sf "$(CURDIR)/contact_map" "$(BINDIR)/contact_map"

