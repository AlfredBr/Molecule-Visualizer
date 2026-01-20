/*
 * MolVis - Molecule Types
 *
 * Platform-agnostic molecule data structures.
 * Shared between Windows (CUDA) and macOS (Metal) implementations.
 */

#ifndef MOLECULE_TYPES_H
#define MOLECULE_TYPES_H

// Molecule limits
#define MAX_ATOMS 200
#define MAX_BONDS 250

// Atom types (CPK coloring)
#define ATOM_H  0   // Hydrogen - white
#define ATOM_C  1   // Carbon - dark gray
#define ATOM_N  2   // Nitrogen - blue
#define ATOM_O  3   // Oxygen - red
#define ATOM_P  4   // Phosphorus - orange
#define ATOM_S  5   // Sulfur - yellow
#define ATOM_CL 6   // Chlorine - green
#define ATOM_BR 7   // Bromine - dark red
#define ATOM_F  8   // Fluorine - light green
#define ATOM_I  9   // Iodine - purple
#define ATOM_NA 10  // Sodium
#define ATOM_SI 11  // Silicon
#define ATOM_B  12  // Boron
#define ATOM_FE 13  // Iron
#define ATOM_CU 14  // Copper
#define ATOM_AL 15  // Aluminum
#define ATOM_TI 16  // Titanium
#define ATOM_PT 17  // Platinum

// Atom structure
struct Atom {
    float x, y, z;
    int type;
    float radius;
};

// Bond structure
struct Bond {
    int atom1, atom2;
    int order;  // 1=single, 2=double, 3=triple
};

// Molecule structure
struct Molecule {
    Atom atoms[MAX_ATOMS];
    Bond bonds[MAX_BONDS];
    int numAtoms;
    int numBonds;
    char name[64];
};

#endif // MOLECULE_TYPES_H
