/*
 * MolVis - Molecule Database Implementation
 *
 * Complete database of 233 molecules organized by category.
 * Molecule builders imported from legacy cuda_molecule.cu
 */

#include "molecule_db.h"
#include <string.h>
#include <stdlib.h>
#include <math.h>

#ifndef PI
#define PI 3.14159265358979323846f
#endif

#ifndef TWO_PI
#define TWO_PI (2.0f * PI)
#endif

// ============== Category Names ==============

static const char* categoryNames[] = {
    "Simple Molecules",       // CAT_SIMPLE
    "Organic Compounds",      // CAT_ORGANIC
    "Amino Acids",           // CAT_AMINO_ACIDS
    "Sugars & Carbohydrates", // CAT_SUGARS
    "Vitamins",              // CAT_VITAMINS
    "Neurotransmitters",     // CAT_NEUROTRANS
    "Hormones",              // CAT_HORMONES
    "Pharmaceuticals",       // CAT_PHARMA
    "Household Chemicals",   // CAT_HOUSEHOLD
    "Acids",                 // CAT_ACIDS
    "Plastic Monomers",      // CAT_PLASTICS
    "Fatty Acids",           // CAT_FATS
    "Metal Compounds",       // CAT_METALS
    "Energy Molecules",      // CAT_ENERGY
    "Flavors & Fragrances",  // CAT_FLAVORS
    "Other"                  // CAT_OTHER
};

const char* molecule_get_category_name(int category) {
    if (category < 0 || category >= CAT_COUNT) return "Unknown";
    return categoryNames[category];
}

// ============== Helper Functions ==============

static void addAtom(Molecule* mol, float x, float y, float z, int type) {
    if (mol->numAtoms >= MAX_ATOMS) return;
    Atom* a = &mol->atoms[mol->numAtoms];
    a->x = x; a->y = y; a->z = z;
    a->type = type;
    float radii[] = {
        0.25f, 0.40f, 0.38f, 0.35f, 0.45f, 0.45f, 0.45f, 0.50f, 0.35f, 0.55f,
        0.55f, 0.48f, 0.42f, 0.55f, 0.50f, 0.50f, 0.52f, 0.55f
    };
    a->radius = radii[type < 18 ? type : 0];
    mol->numAtoms++;
}

static void addBond(Molecule* mol, int a1, int a2, int order) {
    if (mol->numBonds >= MAX_BONDS) return;
    mol->bonds[mol->numBonds].atom1 = a1;
    mol->bonds[mol->numBonds].atom2 = a2;
    mol->bonds[mol->numBonds].order = order;
    mol->numBonds++;
}

static void centerMolecule(Molecule* mol) {
    float cx = 0, cy = 0, cz = 0;
    for (int i = 0; i < mol->numAtoms; i++) {
        cx += mol->atoms[i].x; cy += mol->atoms[i].y; cz += mol->atoms[i].z;
    }
    if (mol->numAtoms > 0) {
        cx /= mol->numAtoms; cy /= mol->numAtoms; cz /= mol->numAtoms;
        for (int i = 0; i < mol->numAtoms; i++) {
            mol->atoms[i].x -= cx; mol->atoms[i].y -= cy; mol->atoms[i].z -= cz;
        }
    }
}

// ============== MOLECULE BUILDERS (233 total) ==============
// Imported from legacy/cuda_molecule.cu

void buildWater(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Water (H2O)");

    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_O);
    addAtom(mol, 0.76f, 0.59f, 0.0f, ATOM_H);
    addAtom(mol, -0.76f, 0.59f, 0.0f, ATOM_H);

    addBond(mol, 0, 1, 1);
    addBond(mol, 0, 2, 1);

    centerMolecule(mol);
}

// Build Methane (CH4)
void buildMethane(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Methane (CH4)");

    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_C);
    addAtom(mol, 0.63f, 0.63f, 0.63f, ATOM_H);
    addAtom(mol, -0.63f, -0.63f, 0.63f, ATOM_H);
    addAtom(mol, -0.63f, 0.63f, -0.63f, ATOM_H);
    addAtom(mol, 0.63f, -0.63f, -0.63f, ATOM_H);

    addBond(mol, 0, 1, 1);
    addBond(mol, 0, 2, 1);
    addBond(mol, 0, 3, 1);
    addBond(mol, 0, 4, 1);

    centerMolecule(mol);
}

// Build Benzene (C6H6)
void buildBenzene(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Benzene (C6H6)");

    float r = 1.4f;
    float rH = 2.2f;

    // Carbon ring
    for (int i = 0; i < 6; i++) {
        float angle = i * PI / 3.0f;
        addAtom(mol, r * cosf(angle), r * sinf(angle), 0.0f, ATOM_C);
    }

    // Hydrogen atoms
    for (int i = 0; i < 6; i++) {
        float angle = i * PI / 3.0f;
        addAtom(mol, rH * cosf(angle), rH * sinf(angle), 0.0f, ATOM_H);
    }

    // C-C bonds (alternating single/double for aromatic)
    for (int i = 0; i < 6; i++) {
        addBond(mol, i, (i + 1) % 6, (i % 2 == 0) ? 2 : 1);
    }

    // C-H bonds
    for (int i = 0; i < 6; i++) {
        addBond(mol, i, i + 6, 1);
    }

    centerMolecule(mol);
}

// Build Ethanol (C2H5OH)
void buildEthanol(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Ethanol (C2H5OH)");

    // Carbons
    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_C);      // C1
    addAtom(mol, 1.5f, 0.0f, 0.0f, ATOM_C);      // C2

    // Oxygen
    addAtom(mol, 2.2f, 1.1f, 0.0f, ATOM_O);      // O

    // Hydrogens on C1
    addAtom(mol, -0.5f, 0.9f, 0.3f, ATOM_H);
    addAtom(mol, -0.5f, -0.5f, 0.8f, ATOM_H);
    addAtom(mol, -0.5f, -0.4f, -0.9f, ATOM_H);

    // Hydrogens on C2
    addAtom(mol, 2.0f, -0.5f, 0.85f, ATOM_H);
    addAtom(mol, 2.0f, -0.5f, -0.85f, ATOM_H);

    // Hydrogen on O
    addAtom(mol, 3.1f, 1.0f, 0.0f, ATOM_H);

    // Bonds
    addBond(mol, 0, 1, 1);  // C-C
    addBond(mol, 1, 2, 1);  // C-O
    addBond(mol, 0, 3, 1);  // C-H
    addBond(mol, 0, 4, 1);
    addBond(mol, 0, 5, 1);
    addBond(mol, 1, 6, 1);
    addBond(mol, 1, 7, 1);
    addBond(mol, 2, 8, 1);  // O-H

    centerMolecule(mol);
}

// Build Caffeine (C8H10N4O2)
void buildCaffeine(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Caffeine (C8H10N4O2)");

    // Purine ring system (approximate coordinates)
    // Imidazole ring
    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_N);       // N1
    addAtom(mol, 1.2f, 0.5f, 0.0f, ATOM_C);       // C2
    addAtom(mol, 1.2f, 1.9f, 0.0f, ATOM_N);       // N3
    addAtom(mol, 0.0f, 2.4f, 0.0f, ATOM_C);       // C4
    addAtom(mol, -0.8f, 1.2f, 0.0f, ATOM_C);      // C5

    // Pyrimidine ring
    addAtom(mol, -0.8f, 3.6f, 0.0f, ATOM_N);      // N6
    addAtom(mol, 0.0f, 4.8f, 0.0f, ATOM_C);       // C7
    addAtom(mol, 1.4f, 4.6f, 0.0f, ATOM_N);       // N8
    addAtom(mol, 1.8f, 3.3f, 0.0f, ATOM_C);       // C9

    // Carbonyl oxygens
    addAtom(mol, 2.3f, -0.2f, 0.0f, ATOM_O);      // O1
    addAtom(mol, -0.4f, 5.9f, 0.0f, ATOM_O);      // O2

    // Methyl groups
    addAtom(mol, -0.5f, -1.3f, 0.0f, ATOM_C);     // CH3 on N1
    addAtom(mol, -2.0f, 3.8f, 0.0f, ATOM_C);      // CH3 on N6
    addAtom(mol, 2.0f, 5.8f, 0.0f, ATOM_C);       // CH3 on N8

    // Hydrogen on imidazole C8 position
    addAtom(mol, 3.0f, 3.2f, 0.0f, ATOM_H);       // H on C9 (C8 position in IUPAC)

    // Hydrogens on methyl groups (simplified - 3 each)
    addAtom(mol, -1.5f, -1.5f, 0.0f, ATOM_H);
    addAtom(mol, -0.1f, -1.8f, 0.8f, ATOM_H);
    addAtom(mol, -0.1f, -1.8f, -0.8f, ATOM_H);

    addAtom(mol, -2.3f, 4.8f, 0.0f, ATOM_H);
    addAtom(mol, -2.5f, 3.3f, 0.8f, ATOM_H);
    addAtom(mol, -2.5f, 3.3f, -0.8f, ATOM_H);

    addAtom(mol, 1.5f, 6.7f, 0.0f, ATOM_H);
    addAtom(mol, 2.5f, 5.8f, 0.9f, ATOM_H);
    addAtom(mol, 2.5f, 5.8f, -0.9f, ATOM_H);

    // Bonds - rings
    addBond(mol, 0, 1, 1);
    addBond(mol, 1, 2, 2);
    addBond(mol, 2, 3, 1);
    addBond(mol, 3, 4, 2);
    addBond(mol, 4, 0, 1);
    addBond(mol, 3, 5, 1);
    addBond(mol, 5, 6, 1);
    addBond(mol, 6, 7, 1);
    addBond(mol, 7, 8, 1);
    addBond(mol, 8, 2, 1);

    // Carbonyl bonds
    addBond(mol, 1, 9, 2);
    addBond(mol, 6, 10, 2);

    // Methyl bonds
    addBond(mol, 0, 11, 1);
    addBond(mol, 5, 12, 1);
    addBond(mol, 7, 13, 1);

    // C-H bond (C9/C8 position)
    addBond(mol, 8, 14, 1);

    // Methyl H bonds
    addBond(mol, 11, 15, 1);
    addBond(mol, 11, 16, 1);
    addBond(mol, 11, 17, 1);
    addBond(mol, 12, 18, 1);
    addBond(mol, 12, 19, 1);
    addBond(mol, 12, 20, 1);
    addBond(mol, 13, 21, 1);
    addBond(mol, 13, 22, 1);
    addBond(mol, 13, 23, 1);

    centerMolecule(mol);
}

// Build Adenine (DNA base)
void buildAdenine(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Adenine (DNA base)");

    // Purine ring
    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_N);       // N1
    addAtom(mol, 1.3f, 0.3f, 0.0f, ATOM_C);       // C2
    addAtom(mol, 1.8f, 1.6f, 0.0f, ATOM_N);       // N3
    addAtom(mol, 0.9f, 2.6f, 0.0f, ATOM_C);       // C4
    addAtom(mol, -0.4f, 2.3f, 0.0f, ATOM_C);      // C5
    addAtom(mol, -0.8f, 1.0f, 0.0f, ATOM_C);      // C6

    // Imidazole ring
    addAtom(mol, 1.2f, 3.9f, 0.0f, ATOM_N);       // N7
    addAtom(mol, 0.0f, 4.5f, 0.0f, ATOM_C);       // C8
    addAtom(mol, -1.0f, 3.5f, 0.0f, ATOM_N);      // N9

    // Amino group
    addAtom(mol, -2.1f, 0.7f, 0.0f, ATOM_N);      // NH2

    // Hydrogens
    addAtom(mol, 2.0f, -0.4f, 0.0f, ATOM_H);      // H on C2
    addAtom(mol, -0.2f, 5.5f, 0.0f, ATOM_H);      // H on C8
    addAtom(mol, -1.9f, 3.8f, 0.0f, ATOM_H);      // H on N9
    addAtom(mol, -2.6f, 1.5f, 0.0f, ATOM_H);      // H on NH2
    addAtom(mol, -2.6f, -0.1f, 0.0f, ATOM_H);     // H on NH2

    // Bonds
    addBond(mol, 0, 1, 2);
    addBond(mol, 1, 2, 1);
    addBond(mol, 2, 3, 2);
    addBond(mol, 3, 4, 1);
    addBond(mol, 4, 5, 2);
    addBond(mol, 5, 0, 1);
    addBond(mol, 3, 6, 1);
    addBond(mol, 6, 7, 2);
    addBond(mol, 7, 8, 1);
    addBond(mol, 8, 4, 1);
    addBond(mol, 5, 9, 1);
    addBond(mol, 1, 10, 1);
    addBond(mol, 7, 11, 1);
    addBond(mol, 8, 12, 1);
    addBond(mol, 9, 13, 1);
    addBond(mol, 9, 14, 1);

    centerMolecule(mol);
}

// Build Glucose (C6H12O6)
void buildGlucose(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Glucose (C6H12O6)");

    // Pyranose ring (chair conformation, simplified)
    float r = 1.4f;

    // Ring carbons and oxygen
    addAtom(mol, r, 0.0f, 0.3f, ATOM_C);          // C1
    addAtom(mol, r * 0.5f, r * 0.866f, -0.3f, ATOM_C);   // C2
    addAtom(mol, -r * 0.5f, r * 0.866f, 0.3f, ATOM_C);   // C3
    addAtom(mol, -r, 0.0f, -0.3f, ATOM_C);        // C4
    addAtom(mol, -r * 0.5f, -r * 0.866f, 0.3f, ATOM_C);  // C5
    addAtom(mol, r * 0.5f, -r * 0.866f, -0.3f, ATOM_O);  // Ring O

    // C6 (CH2OH group)
    addAtom(mol, -r * 0.9f, -r * 1.5f, 0.0f, ATOM_C);    // C6

    // OH groups
    addAtom(mol, r * 1.5f, 0.3f, 1.0f, ATOM_O);          // O on C1
    addAtom(mol, r * 0.9f, r * 1.4f, -1.0f, ATOM_O);     // O on C2
    addAtom(mol, -r * 0.9f, r * 1.4f, 1.0f, ATOM_O);     // O on C3
    addAtom(mol, -r * 1.5f, 0.0f, -1.0f, ATOM_O);        // O on C4
    addAtom(mol, -r * 0.5f, -r * 2.3f, 0.0f, ATOM_O);    // O on C6

    // Hydrogens (simplified - one per carbon/oxygen)
    addAtom(mol, r * 1.3f, -0.5f, -0.5f, ATOM_H);        // H on C1
    addAtom(mol, r * 0.8f, r * 0.5f, 0.5f, ATOM_H);      // H on C2
    addAtom(mol, -r * 0.2f, r * 1.1f, -0.5f, ATOM_H);    // H on C3
    addAtom(mol, -r * 0.7f, -0.3f, 0.5f, ATOM_H);        // H on C4
    addAtom(mol, -r * 0.8f, -r * 0.6f, -0.5f, ATOM_H);   // H on C5
    addAtom(mol, -r * 1.5f, -r * 1.3f, 0.7f, ATOM_H);    // H on C6
    addAtom(mol, -r * 1.3f, -r * 1.5f, -0.8f, ATOM_H);   // H on C6

    // OH hydrogens
    addAtom(mol, r * 2.2f, 0.0f, 1.3f, ATOM_H);
    addAtom(mol, r * 1.5f, r * 1.2f, -1.5f, ATOM_H);
    addAtom(mol, -r * 1.5f, r * 1.2f, 1.5f, ATOM_H);
    addAtom(mol, -r * 2.2f, 0.0f, -1.3f, ATOM_H);
    addAtom(mol, -r * 0.9f, -r * 2.8f, 0.5f, ATOM_H);

    // Ring bonds
    addBond(mol, 0, 1, 1);
    addBond(mol, 1, 2, 1);
    addBond(mol, 2, 3, 1);
    addBond(mol, 3, 4, 1);
    addBond(mol, 4, 5, 1);
    addBond(mol, 5, 0, 1);

    // C6 bond
    addBond(mol, 4, 6, 1);

    // C-O bonds
    addBond(mol, 0, 7, 1);
    addBond(mol, 1, 8, 1);
    addBond(mol, 2, 9, 1);
    addBond(mol, 3, 10, 1);
    addBond(mol, 6, 11, 1);

    // C-H bonds
    addBond(mol, 0, 12, 1);
    addBond(mol, 1, 13, 1);
    addBond(mol, 2, 14, 1);
    addBond(mol, 3, 15, 1);
    addBond(mol, 4, 16, 1);
    addBond(mol, 6, 17, 1);
    addBond(mol, 6, 18, 1);

    // O-H bonds
    addBond(mol, 7, 19, 1);
    addBond(mol, 8, 20, 1);
    addBond(mol, 9, 21, 1);
    addBond(mol, 10, 22, 1);
    addBond(mol, 11, 23, 1);

    centerMolecule(mol);
}

// Build Aspirin (C9H8O4)
void buildAspirin(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Aspirin/Bayer (C9H8O4)");

    float r = 1.4f;

    // Benzene ring
    for (int i = 0; i < 6; i++) {
        float angle = i * PI / 3.0f;
        addAtom(mol, r * cosf(angle), r * sinf(angle), 0.0f, ATOM_C);
    }

    // Carboxylic acid group (COOH)
    addAtom(mol, 2.5f, 0.5f, 0.0f, ATOM_C);       // C7
    addAtom(mol, 3.2f, 1.5f, 0.0f, ATOM_O);       // O1 (=O)
    addAtom(mol, 3.0f, -0.6f, 0.0f, ATOM_O);      // O2 (OH)

    // Acetyl group (OCOCH3)
    addAtom(mol, r * cosf(PI/3.0f) - 0.8f, r * sinf(PI/3.0f) + 0.8f, 0.0f, ATOM_O);  // O3
    addAtom(mol, r * cosf(PI/3.0f) - 1.0f, r * sinf(PI/3.0f) + 2.2f, 0.0f, ATOM_C);  // C8
    addAtom(mol, r * cosf(PI/3.0f) - 2.3f, r * sinf(PI/3.0f) + 2.8f, 0.0f, ATOM_O);  // O4 (=O)
    addAtom(mol, r * cosf(PI/3.0f) + 0.2f, r * sinf(PI/3.0f) + 3.2f, 0.0f, ATOM_C);  // C9 (CH3)

    // Hydrogens on benzene (4 of them, positions 2,3,4,5)
    float rH = 2.4f;
    addAtom(mol, rH * cosf(2*PI/3.0f), rH * sinf(2*PI/3.0f), 0.0f, ATOM_H);
    addAtom(mol, rH * cosf(PI), rH * sinf(PI), 0.0f, ATOM_H);
    addAtom(mol, rH * cosf(4*PI/3.0f), rH * sinf(4*PI/3.0f), 0.0f, ATOM_H);
    addAtom(mol, rH * cosf(5*PI/3.0f), rH * sinf(5*PI/3.0f), 0.0f, ATOM_H);

    // H on COOH
    addAtom(mol, 3.8f, -0.5f, 0.0f, ATOM_H);

    // H on CH3
    addAtom(mol, r * cosf(PI/3.0f) + 0.0f, r * sinf(PI/3.0f) + 4.1f, 0.0f, ATOM_H);
    addAtom(mol, r * cosf(PI/3.0f) + 0.9f, r * sinf(PI/3.0f) + 2.8f, 0.7f, ATOM_H);
    addAtom(mol, r * cosf(PI/3.0f) + 0.9f, r * sinf(PI/3.0f) + 2.8f, -0.7f, ATOM_H);

    // Benzene bonds
    for (int i = 0; i < 6; i++) {
        addBond(mol, i, (i + 1) % 6, (i % 2 == 0) ? 2 : 1);
    }

    // COOH bonds
    addBond(mol, 0, 6, 1);    // C1-C7
    addBond(mol, 6, 7, 2);    // C7=O1
    addBond(mol, 6, 8, 1);    // C7-O2

    // Acetyl bonds
    addBond(mol, 1, 9, 1);    // C2-O3
    addBond(mol, 9, 10, 1);   // O3-C8
    addBond(mol, 10, 11, 2);  // C8=O4
    addBond(mol, 10, 12, 1);  // C8-C9

    // C-H bonds
    addBond(mol, 2, 13, 1);
    addBond(mol, 3, 14, 1);
    addBond(mol, 4, 15, 1);
    addBond(mol, 5, 16, 1);
    addBond(mol, 8, 17, 1);
    addBond(mol, 12, 18, 1);
    addBond(mol, 12, 19, 1);
    addBond(mol, 12, 20, 1);

    centerMolecule(mol);
}

// Build Ammonia (NH3)
void buildAmmonia(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Ammonia (NH3)");

    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_N);
    addAtom(mol, 0.94f, 0.0f, 0.34f, ATOM_H);
    addAtom(mol, -0.47f, 0.81f, 0.34f, ATOM_H);
    addAtom(mol, -0.47f, -0.81f, 0.34f, ATOM_H);

    addBond(mol, 0, 1, 1);
    addBond(mol, 0, 2, 1);
    addBond(mol, 0, 3, 1);

    centerMolecule(mol);
}

// Build Sodium Hypochlorite / Bleach (NaOCl)
void buildBleach(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Bleach/NaOCl");

    // Ionic compound: Na+ and OCl-
    addAtom(mol, -1.5f, 0.0f, 0.0f, ATOM_NA);   // Na+
    addAtom(mol, 0.5f, 0.0f, 0.0f, ATOM_O);     // O
    addAtom(mol, 2.0f, 0.0f, 0.0f, ATOM_CL);    // Cl

    // O-Cl bond (hypochlorite ion)
    addBond(mol, 1, 2, 1);
    // Na-O ionic interaction (shown as single bond for visualization)
    addBond(mol, 0, 1, 1);

    centerMolecule(mol);
}

// Build Sodium Hydroxide / Lye (NaOH)
void buildLye(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Lye/NaOH");

    // Ionic compound: Na+ and OH-
    addAtom(mol, -1.2f, 0.0f, 0.0f, ATOM_NA);   // Na+
    addAtom(mol, 0.5f, 0.0f, 0.0f, ATOM_O);     // O
    addAtom(mol, 1.5f, 0.0f, 0.0f, ATOM_H);     // H

    // O-H bond (hydroxide ion)
    addBond(mol, 1, 2, 1);
    // Na-O ionic interaction (shown as single bond for visualization)
    addBond(mol, 0, 1, 1);

    centerMolecule(mol);
}

// Build Table Salt (NaCl)
void buildTableSalt(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Table Salt (NaCl)");

    addAtom(mol, -1.2f, 0.0f, 0.0f, ATOM_NA);   // Na+
    addAtom(mol, 1.2f, 0.0f, 0.0f, ATOM_CL);    // Cl-

    addBond(mol, 0, 1, 1);  // Ionic bond shown as single

    centerMolecule(mol);
}

// Build Baking Soda (NaHCO3)
void buildBakingSoda(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Baking Soda (NaHCO3)");

    addAtom(mol, -2.0f, 0.0f, 0.0f, ATOM_NA);   // Na+
    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_C);     // C (carbonate center)
    addAtom(mol, 1.2f, 0.0f, 0.0f, ATOM_O);     // O (double bonded)
    addAtom(mol, -0.6f, 1.0f, 0.0f, ATOM_O);    // O- (ionic to Na)
    addAtom(mol, -0.6f, -1.0f, 0.0f, ATOM_O);   // OH
    addAtom(mol, -0.6f, -2.0f, 0.0f, ATOM_H);   // H

    addBond(mol, 1, 2, 2);  // C=O
    addBond(mol, 1, 3, 1);  // C-O-
    addBond(mol, 1, 4, 1);  // C-OH
    addBond(mol, 4, 5, 1);  // O-H
    addBond(mol, 0, 3, 1);  // Na-O ionic

    centerMolecule(mol);
}

// Build Washing Soda (Na2CO3)
void buildWashingSoda(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Washing Soda (Na2CO3)");

    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_C);     // C (carbonate center)
    addAtom(mol, 1.2f, 0.0f, 0.0f, ATOM_O);     // O (double bonded)
    addAtom(mol, -0.6f, 1.0f, 0.0f, ATOM_O);    // O-
    addAtom(mol, -0.6f, -1.0f, 0.0f, ATOM_O);   // O-
    addAtom(mol, -1.8f, 1.8f, 0.0f, ATOM_NA);   // Na+
    addAtom(mol, -1.8f, -1.8f, 0.0f, ATOM_NA);  // Na+

    addBond(mol, 0, 1, 2);  // C=O
    addBond(mol, 0, 2, 1);  // C-O-
    addBond(mol, 0, 3, 1);  // C-O-
    addBond(mol, 2, 4, 1);  // O-Na ionic
    addBond(mol, 3, 5, 1);  // O-Na ionic

    centerMolecule(mol);
}

// Build Sodium Fluoride (NaF) - toothpaste
void buildSodiumFluoride(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Sodium Fluoride (NaF)");

    addAtom(mol, -1.1f, 0.0f, 0.0f, ATOM_NA);   // Na+
    addAtom(mol, 1.1f, 0.0f, 0.0f, ATOM_F);     // F-

    addBond(mol, 0, 1, 1);

    centerMolecule(mol);
}

// Build Sodium Nitrate (NaNO3) - fertilizer/preservative
void buildSodiumNitrate(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Sodium Nitrate (NaNO3)");

    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_N);     // N (center)
    addAtom(mol, 1.2f, 0.0f, 0.0f, ATOM_O);     // O
    addAtom(mol, -0.6f, 1.0f, 0.0f, ATOM_O);    // O
    addAtom(mol, -0.6f, -1.0f, 0.0f, ATOM_O);   // O-
    addAtom(mol, -1.8f, -1.8f, 0.0f, ATOM_NA);  // Na+

    addBond(mol, 0, 1, 2);  // N=O
    addBond(mol, 0, 2, 2);  // N=O
    addBond(mol, 0, 3, 1);  // N-O-
    addBond(mol, 3, 4, 1);  // O-Na ionic

    centerMolecule(mol);
}

// Build Sodium Nitrite (NaNO2) - cured meats
void buildSodiumNitrite(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Sodium Nitrite (NaNO2)");

    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_N);     // N (center)
    addAtom(mol, 1.2f, 0.0f, 0.0f, ATOM_O);     // O
    addAtom(mol, -0.6f, -1.0f, 0.0f, ATOM_O);   // O-
    addAtom(mol, -1.8f, -1.8f, 0.0f, ATOM_NA);  // Na+

    addBond(mol, 0, 1, 2);  // N=O
    addBond(mol, 0, 2, 1);  // N-O-
    addBond(mol, 2, 3, 1);  // O-Na ionic

    centerMolecule(mol);
}

// Build Sodium Sulfate (Na2SO4) - Glauber's salt, detergents
void buildSodiumSulfate(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Sodium Sulfate (Na2SO4)");

    // Sulfate ion SO4 2- with two Na+
    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_S);      // S (center)
    addAtom(mol, 1.4f, 0.0f, 0.0f, ATOM_O);      // O
    addAtom(mol, -1.4f, 0.0f, 0.0f, ATOM_O);     // O
    addAtom(mol, 0.0f, 1.4f, 0.0f, ATOM_O);      // O-
    addAtom(mol, 0.0f, -1.4f, 0.0f, ATOM_O);     // O-
    addAtom(mol, 0.0f, 2.8f, 0.0f, ATOM_NA);     // Na+
    addAtom(mol, 0.0f, -2.8f, 0.0f, ATOM_NA);    // Na+

    addBond(mol, 0, 1, 2);  // S=O
    addBond(mol, 0, 2, 2);  // S=O
    addBond(mol, 0, 3, 1);  // S-O-
    addBond(mol, 0, 4, 1);  // S-O-
    addBond(mol, 3, 5, 1);  // O-Na ionic
    addBond(mol, 4, 6, 1);  // O-Na ionic

    centerMolecule(mol);
}

// Build MSG - Monosodium Glutamate (C5H8NNaO4)
void buildMSG(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "MSG (C5H8NNaO4)");

    // Glutamate backbone
    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_C);      // C1 alpha carbon
    addAtom(mol, -1.3f, 0.5f, 0.0f, ATOM_N);     // NH2 (amino)
    addAtom(mol, 1.3f, 0.7f, 0.0f, ATOM_C);      // C2 carboxyl
    addAtom(mol, 1.3f, 2.0f, 0.0f, ATOM_O);      // O (C=O)
    addAtom(mol, 2.4f, 0.0f, 0.0f, ATOM_O);      // O- (carboxylate)
    addAtom(mol, 0.0f, -1.5f, 0.0f, ATOM_C);     // C3 (CH2)
    addAtom(mol, 0.0f, -3.0f, 0.0f, ATOM_C);     // C4 (CH2)
    addAtom(mol, 0.0f, -4.5f, 0.0f, ATOM_C);     // C5 carboxyl
    addAtom(mol, 1.1f, -5.2f, 0.0f, ATOM_O);     // O (C=O)
    addAtom(mol, -1.1f, -5.2f, 0.0f, ATOM_O);    // O-
    addAtom(mol, 3.5f, 0.0f, 0.0f, ATOM_NA);     // Na+

    // Hydrogens
    addAtom(mol, 0.0f, 0.5f, 0.9f, ATOM_H);      // 11: H on alpha C
    addAtom(mol, -1.4f, 1.5f, 0.0f, ATOM_H);     // 12: H on NH2
    addAtom(mol, -2.1f, 0.0f, 0.0f, ATOM_H);     // 13: H on NH2
    addAtom(mol, 0.9f, -1.5f, 0.5f, ATOM_H);     // 14: H on CH2
    addAtom(mol, -0.9f, -1.5f, 0.5f, ATOM_H);    // 15: H on CH2
    addAtom(mol, 0.9f, -3.0f, 0.5f, ATOM_H);     // 16: H on CH2
    addAtom(mol, -0.9f, -3.0f, 0.5f, ATOM_H);    // 17: H on CH2
    addAtom(mol, -1.8f, -5.8f, 0.0f, ATOM_H);    // 18: H on terminal COOH

    addBond(mol, 0, 1, 1);   // C-NH2
    addBond(mol, 0, 2, 1);   // C-COOH
    addBond(mol, 2, 3, 2);   // C=O
    addBond(mol, 2, 4, 1);   // C-O-
    addBond(mol, 0, 5, 1);   // C-CH2
    addBond(mol, 5, 6, 1);   // CH2-CH2
    addBond(mol, 6, 7, 1);   // CH2-COOH
    addBond(mol, 7, 8, 2);   // C=O
    addBond(mol, 7, 9, 1);   // C-OH (protonated end)
    addBond(mol, 4, 10, 1);  // O-Na+ ionic
    addBond(mol, 0, 11, 1);
    addBond(mol, 1, 12, 1);
    addBond(mol, 1, 13, 1);
    addBond(mol, 5, 14, 1);
    addBond(mol, 5, 15, 1);
    addBond(mol, 6, 16, 1);
    addBond(mol, 6, 17, 1);
    addBond(mol, 9, 18, 1);  // O-H on terminal COOH

    centerMolecule(mol);
}

// Build Citric Acid (C6H8O7) - citrus fruits
void buildCitricAcid(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Citric Acid (C6H8O7)");

    // Central carbon with OH
    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_C);      // C1 central (with OH)
    addAtom(mol, 0.0f, 1.4f, 0.0f, ATOM_O);      // OH on central
    addAtom(mol, 0.0f, 2.4f, 0.0f, ATOM_H);      // H on OH

    // Carboxyl on central carbon
    addAtom(mol, 0.0f, 0.0f, 1.5f, ATOM_C);      // COOH carbon
    addAtom(mol, 1.0f, 0.0f, 2.2f, ATOM_O);      // =O
    addAtom(mol, -1.0f, 0.0f, 2.2f, ATOM_O);     // OH
    addAtom(mol, -1.0f, 0.0f, 3.2f, ATOM_H);     // H

    // CH2-COOH arm 1
    addAtom(mol, 1.5f, 0.0f, -0.3f, ATOM_C);     // CH2
    addAtom(mol, 2.8f, 0.0f, 0.3f, ATOM_C);      // COOH
    addAtom(mol, 3.5f, 1.0f, 0.0f, ATOM_O);      // =O
    addAtom(mol, 3.2f, -1.0f, 1.0f, ATOM_O);     // OH
    addAtom(mol, 4.0f, -1.0f, 1.5f, ATOM_H);     // H

    // CH2-COOH arm 2
    addAtom(mol, -1.5f, 0.0f, -0.3f, ATOM_C);    // CH2
    addAtom(mol, -2.8f, 0.0f, 0.3f, ATOM_C);     // COOH
    addAtom(mol, -3.5f, 1.0f, 0.0f, ATOM_O);     // =O
    addAtom(mol, -3.2f, -1.0f, 1.0f, ATOM_O);    // OH
    addAtom(mol, -4.0f, -1.0f, 1.5f, ATOM_H);    // H

    // CH2 hydrogens
    addAtom(mol, 1.5f, 0.9f, -0.9f, ATOM_H);
    addAtom(mol, 1.5f, -0.9f, -0.9f, ATOM_H);
    addAtom(mol, -1.5f, 0.9f, -0.9f, ATOM_H);
    addAtom(mol, -1.5f, -0.9f, -0.9f, ATOM_H);

    // Bonds
    addBond(mol, 0, 1, 1);   // C-OH
    addBond(mol, 1, 2, 1);   // O-H
    addBond(mol, 0, 3, 1);   // C-COOH
    addBond(mol, 3, 4, 2);   // C=O
    addBond(mol, 3, 5, 1);   // C-OH
    addBond(mol, 5, 6, 1);   // O-H
    addBond(mol, 0, 7, 1);   // C-CH2
    addBond(mol, 7, 8, 1);   // CH2-COOH
    addBond(mol, 8, 9, 2);   // C=O
    addBond(mol, 8, 10, 1);  // C-OH
    addBond(mol, 10, 11, 1); // O-H
    addBond(mol, 0, 12, 1);  // C-CH2
    addBond(mol, 12, 13, 1); // CH2-COOH
    addBond(mol, 13, 14, 2); // C=O
    addBond(mol, 13, 15, 1); // C-OH
    addBond(mol, 15, 16, 1); // O-H
    addBond(mol, 7, 17, 1);
    addBond(mol, 7, 18, 1);
    addBond(mol, 12, 19, 1);
    addBond(mol, 12, 20, 1);

    centerMolecule(mol);
}

// Build Carbonic Acid (H2CO3) - soda/carbonation
void buildCarbonicAcid(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Carbonic Acid (H2CO3)");

    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_C);     // C
    addAtom(mol, 1.2f, 0.0f, 0.0f, ATOM_O);     // =O
    addAtom(mol, -0.6f, 1.0f, 0.0f, ATOM_O);    // OH
    addAtom(mol, -0.6f, -1.0f, 0.0f, ATOM_O);   // OH
    addAtom(mol, -0.6f, 1.9f, 0.0f, ATOM_H);    // H
    addAtom(mol, -0.6f, -1.9f, 0.0f, ATOM_H);   // H

    addBond(mol, 0, 1, 2);  // C=O
    addBond(mol, 0, 2, 1);  // C-OH
    addBond(mol, 0, 3, 1);  // C-OH
    addBond(mol, 2, 4, 1);  // O-H
    addBond(mol, 3, 5, 1);  // O-H

    centerMolecule(mol);
}

// Build Boric Acid (H3BO3) - antiseptic
void buildBoricAcid(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Boric Acid (H3BO3)");

    // Boron center with three OH groups
    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_B);     // B (boron)
    addAtom(mol, 1.3f, 0.0f, 0.0f, ATOM_O);     // OH
    addAtom(mol, -0.65f, 1.1f, 0.0f, ATOM_O);   // OH
    addAtom(mol, -0.65f, -1.1f, 0.0f, ATOM_O);  // OH
    addAtom(mol, 2.1f, 0.0f, 0.0f, ATOM_H);
    addAtom(mol, -0.65f, 1.9f, 0.0f, ATOM_H);
    addAtom(mol, -0.65f, -1.9f, 0.0f, ATOM_H);

    addBond(mol, 0, 1, 1);
    addBond(mol, 0, 2, 1);
    addBond(mol, 0, 3, 1);
    addBond(mol, 1, 4, 1);
    addBond(mol, 2, 5, 1);
    addBond(mol, 3, 6, 1);

    centerMolecule(mol);
}

// Build Oxalic Acid (C2H2O4) - found in spinach/rhubarb
void buildOxalicAcid(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Oxalic Acid (C2H2O4)");

    addAtom(mol, -0.7f, 0.0f, 0.0f, ATOM_C);    // C1
    addAtom(mol, 0.7f, 0.0f, 0.0f, ATOM_C);     // C2
    addAtom(mol, -1.4f, 1.0f, 0.0f, ATOM_O);    // =O
    addAtom(mol, -1.4f, -1.0f, 0.0f, ATOM_O);   // OH
    addAtom(mol, 1.4f, 1.0f, 0.0f, ATOM_O);     // =O
    addAtom(mol, 1.4f, -1.0f, 0.0f, ATOM_O);    // OH
    addAtom(mol, -1.4f, -1.8f, 0.0f, ATOM_H);
    addAtom(mol, 1.4f, -1.8f, 0.0f, ATOM_H);

    addBond(mol, 0, 1, 1);  // C-C
    addBond(mol, 0, 2, 2);  // C=O
    addBond(mol, 0, 3, 1);  // C-OH
    addBond(mol, 1, 4, 2);  // C=O
    addBond(mol, 1, 5, 1);  // C-OH
    addBond(mol, 3, 6, 1);  // O-H
    addBond(mol, 5, 7, 1);  // O-H

    centerMolecule(mol);
}

// Build Tartaric Acid (C4H6O6) - wine/cream of tartar
void buildTartaricAcid(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Tartaric Acid (C4H6O6)");

    // HOOC-CHOH-CHOH-COOH
    addAtom(mol, -2.5f, 0.0f, 0.0f, ATOM_C);    // COOH carbon
    addAtom(mol, -3.2f, 1.0f, 0.0f, ATOM_O);    // =O
    addAtom(mol, -3.2f, -1.0f, 0.0f, ATOM_O);   // OH
    addAtom(mol, -3.2f, -1.8f, 0.0f, ATOM_H);
    addAtom(mol, -1.0f, 0.0f, 0.0f, ATOM_C);    // CHOH
    addAtom(mol, -1.0f, 1.4f, 0.0f, ATOM_O);    // OH
    addAtom(mol, -1.0f, 2.2f, 0.0f, ATOM_H);
    addAtom(mol, -1.0f, -0.5f, 0.9f, ATOM_H);   // H on C
    addAtom(mol, 1.0f, 0.0f, 0.0f, ATOM_C);     // CHOH
    addAtom(mol, 1.0f, 1.4f, 0.0f, ATOM_O);     // OH
    addAtom(mol, 1.0f, 2.2f, 0.0f, ATOM_H);
    addAtom(mol, 1.0f, -0.5f, 0.9f, ATOM_H);    // H on C
    addAtom(mol, 2.5f, 0.0f, 0.0f, ATOM_C);     // COOH carbon
    addAtom(mol, 3.2f, 1.0f, 0.0f, ATOM_O);     // =O
    addAtom(mol, 3.2f, -1.0f, 0.0f, ATOM_O);    // OH
    addAtom(mol, 3.2f, -1.8f, 0.0f, ATOM_H);

    addBond(mol, 0, 1, 2);   // C=O
    addBond(mol, 0, 2, 1);   // C-OH
    addBond(mol, 2, 3, 1);   // O-H
    addBond(mol, 0, 4, 1);   // C-CHOH
    addBond(mol, 4, 5, 1);   // C-OH
    addBond(mol, 5, 6, 1);   // O-H
    addBond(mol, 4, 7, 1);   // C-H
    addBond(mol, 4, 8, 1);   // CHOH-CHOH
    addBond(mol, 8, 9, 1);   // C-OH
    addBond(mol, 9, 10, 1);  // O-H
    addBond(mol, 8, 11, 1);  // C-H
    addBond(mol, 8, 12, 1);  // C-COOH
    addBond(mol, 12, 13, 2); // C=O
    addBond(mol, 12, 14, 1); // C-OH
    addBond(mol, 14, 15, 1); // O-H

    centerMolecule(mol);
}

// Build Malic Acid (C4H6O5) - apples
void buildMalicAcid(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Malic Acid (C4H6O5)");

    // HOOC-CH2-CHOH-COOH
    addAtom(mol, -2.5f, 0.0f, 0.0f, ATOM_C);    // COOH carbon
    addAtom(mol, -3.2f, 1.0f, 0.0f, ATOM_O);    // =O
    addAtom(mol, -3.2f, -1.0f, 0.0f, ATOM_O);   // OH
    addAtom(mol, -3.2f, -1.8f, 0.0f, ATOM_H);
    addAtom(mol, -1.0f, 0.0f, 0.0f, ATOM_C);    // CH2
    addAtom(mol, -1.0f, 0.5f, 0.9f, ATOM_H);
    addAtom(mol, -1.0f, -0.5f, 0.9f, ATOM_H);
    addAtom(mol, 0.5f, 0.0f, 0.0f, ATOM_C);     // CHOH
    addAtom(mol, 0.5f, 1.4f, 0.0f, ATOM_O);     // OH
    addAtom(mol, 0.5f, 2.2f, 0.0f, ATOM_H);
    addAtom(mol, 0.5f, -0.5f, 0.9f, ATOM_H);    // H on C
    addAtom(mol, 2.0f, 0.0f, 0.0f, ATOM_C);     // COOH carbon
    addAtom(mol, 2.7f, 1.0f, 0.0f, ATOM_O);     // =O
    addAtom(mol, 2.7f, -1.0f, 0.0f, ATOM_O);    // OH
    addAtom(mol, 2.7f, -1.8f, 0.0f, ATOM_H);

    addBond(mol, 0, 1, 2);   // C=O
    addBond(mol, 0, 2, 1);   // C-OH
    addBond(mol, 2, 3, 1);   // O-H
    addBond(mol, 0, 4, 1);   // C-CH2
    addBond(mol, 4, 5, 1);   // CH2-H
    addBond(mol, 4, 6, 1);   // CH2-H
    addBond(mol, 4, 7, 1);   // CH2-CHOH
    addBond(mol, 7, 8, 1);   // C-OH
    addBond(mol, 8, 9, 1);   // O-H
    addBond(mol, 7, 10, 1);  // C-H
    addBond(mol, 7, 11, 1);  // C-COOH
    addBond(mol, 11, 12, 2); // C=O
    addBond(mol, 11, 13, 1); // C-OH
    addBond(mol, 13, 14, 1); // O-H

    centerMolecule(mol);
}

// Build Hydrofluoric Acid (HF)
void buildHydrofluoricAcid(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Hydrofluoric Acid (HF)");

    addAtom(mol, -0.5f, 0.0f, 0.0f, ATOM_H);
    addAtom(mol, 0.5f, 0.0f, 0.0f, ATOM_F);

    addBond(mol, 0, 1, 1);

    centerMolecule(mol);
}

// Build Hydrobromic Acid (HBr)
void buildHydrobromicAcid(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Hydrobromic Acid (HBr)");

    addAtom(mol, -0.7f, 0.0f, 0.0f, ATOM_H);
    addAtom(mol, 0.7f, 0.0f, 0.0f, ATOM_BR);

    addBond(mol, 0, 1, 1);

    centerMolecule(mol);
}

// Build Hydroiodic Acid (HI)
void buildHydroiodicAcid(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Hydroiodic Acid (HI)");

    addAtom(mol, -0.8f, 0.0f, 0.0f, ATOM_H);
    addAtom(mol, 0.8f, 0.0f, 0.0f, ATOM_I);

    addBond(mol, 0, 1, 1);

    centerMolecule(mol);
}

// Build Perchloric Acid (HClO4)
void buildPerchloricAcid(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Perchloric Acid (HClO4)");

    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_CL);    // Cl (center)
    addAtom(mol, 1.4f, 0.0f, 0.0f, ATOM_O);     // =O
    addAtom(mol, -1.4f, 0.0f, 0.0f, ATOM_O);    // =O
    addAtom(mol, 0.0f, 1.4f, 0.0f, ATOM_O);     // =O
    addAtom(mol, 0.0f, -1.4f, 0.0f, ATOM_O);    // OH
    addAtom(mol, 0.0f, -2.3f, 0.0f, ATOM_H);

    addBond(mol, 0, 1, 2);
    addBond(mol, 0, 2, 2);
    addBond(mol, 0, 3, 2);
    addBond(mol, 0, 4, 1);
    addBond(mol, 4, 5, 1);

    centerMolecule(mol);
}

// Build Acrylic Acid (C3H4O2) - used in polymers
void buildAcrylicAcid(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Acrylic Acid (C3H4O2)");

    addAtom(mol, -1.2f, 0.0f, 0.0f, ATOM_C);    // CH2=
    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_C);     // =CH
    addAtom(mol, 1.3f, 0.0f, 0.0f, ATOM_C);     // COOH
    addAtom(mol, 2.0f, 1.0f, 0.0f, ATOM_O);     // =O
    addAtom(mol, 2.0f, -1.0f, 0.0f, ATOM_O);    // OH
    addAtom(mol, 2.8f, -1.0f, 0.0f, ATOM_H);
    addAtom(mol, -1.8f, 0.9f, 0.0f, ATOM_H);
    addAtom(mol, -1.8f, -0.9f, 0.0f, ATOM_H);
    addAtom(mol, 0.0f, 1.0f, 0.0f, ATOM_H);

    addBond(mol, 0, 1, 2);  // C=C
    addBond(mol, 1, 2, 1);  // C-COOH
    addBond(mol, 2, 3, 2);  // C=O
    addBond(mol, 2, 4, 1);  // C-OH
    addBond(mol, 4, 5, 1);  // O-H
    addBond(mol, 0, 6, 1);
    addBond(mol, 0, 7, 1);
    addBond(mol, 1, 8, 1);

    centerMolecule(mol);
}

// ============== PLASTICS (Monomers/Oligomers) ==============

// Build Vinyl Chloride (C2H3Cl) - PVC monomer
void buildVinylChloride(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Vinyl Chloride/PVC (C2H3Cl)");

    addAtom(mol, -0.6f, 0.0f, 0.0f, ATOM_C);    // CH2=
    addAtom(mol, 0.6f, 0.0f, 0.0f, ATOM_C);     // =CHCl
    addAtom(mol, 1.5f, 0.0f, 0.0f, ATOM_CL);    // Cl
    addAtom(mol, -1.1f, 0.9f, 0.0f, ATOM_H);
    addAtom(mol, -1.1f, -0.9f, 0.0f, ATOM_H);
    addAtom(mol, 0.6f, 1.0f, 0.0f, ATOM_H);

    addBond(mol, 0, 1, 2);  // C=C
    addBond(mol, 1, 2, 1);  // C-Cl
    addBond(mol, 0, 3, 1);
    addBond(mol, 0, 4, 1);
    addBond(mol, 1, 5, 1);

    centerMolecule(mol);
}

// Build Methyl Methacrylate (C5H8O2) - PMMA/Plexiglas monomer
void buildMethylMethacrylate(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "MMA/Plexiglas (C5H8O2)");

    // CH2=C(CH3)-COOCH3
    addAtom(mol, -1.8f, 0.0f, 0.0f, ATOM_C);    // CH2=
    addAtom(mol, -0.6f, 0.0f, 0.0f, ATOM_C);    // =C(CH3)
    addAtom(mol, -0.6f, 1.4f, 0.0f, ATOM_C);    // CH3
    addAtom(mol, 0.7f, -0.5f, 0.0f, ATOM_C);    // C(=O)
    addAtom(mol, 0.7f, -1.7f, 0.0f, ATOM_O);    // =O
    addAtom(mol, 1.9f, 0.2f, 0.0f, ATOM_O);     // O-
    addAtom(mol, 3.1f, -0.3f, 0.0f, ATOM_C);    // OCH3
    addAtom(mol, -2.3f, 0.9f, 0.0f, ATOM_H);
    addAtom(mol, -2.3f, -0.9f, 0.0f, ATOM_H);
    addAtom(mol, -0.1f, 1.9f, 0.9f, ATOM_H);
    addAtom(mol, -0.1f, 1.9f, -0.9f, ATOM_H);
    addAtom(mol, -1.6f, 1.8f, 0.0f, ATOM_H);
    addAtom(mol, 3.6f, 0.2f, 0.9f, ATOM_H);
    addAtom(mol, 3.6f, 0.2f, -0.9f, ATOM_H);
    addAtom(mol, 3.1f, -1.3f, 0.0f, ATOM_H);

    addBond(mol, 0, 1, 2);
    addBond(mol, 1, 2, 1);
    addBond(mol, 1, 3, 1);
    addBond(mol, 3, 4, 2);
    addBond(mol, 3, 5, 1);
    addBond(mol, 5, 6, 1);
    addBond(mol, 0, 7, 1);
    addBond(mol, 0, 8, 1);
    addBond(mol, 2, 9, 1);
    addBond(mol, 2, 10, 1);
    addBond(mol, 2, 11, 1);
    addBond(mol, 6, 12, 1);
    addBond(mol, 6, 13, 1);
    addBond(mol, 6, 14, 1);

    centerMolecule(mol);
}

// Build Tetrafluoroethylene (C2F4) - PTFE/Teflon monomer
void buildTetrafluoroethylene(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "TFE/Teflon (C2F4)");

    addAtom(mol, -0.6f, 0.0f, 0.0f, ATOM_C);
    addAtom(mol, 0.6f, 0.0f, 0.0f, ATOM_C);
    addAtom(mol, -1.4f, 0.9f, 0.0f, ATOM_F);
    addAtom(mol, -1.4f, -0.9f, 0.0f, ATOM_F);
    addAtom(mol, 1.4f, 0.9f, 0.0f, ATOM_F);
    addAtom(mol, 1.4f, -0.9f, 0.0f, ATOM_F);

    addBond(mol, 0, 1, 2);
    addBond(mol, 0, 2, 1);
    addBond(mol, 0, 3, 1);
    addBond(mol, 1, 4, 1);
    addBond(mol, 1, 5, 1);

    centerMolecule(mol);
}

// Build Caprolactam (C6H11NO) - Nylon-6 precursor
void buildCaprolactam(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Caprolactam/Nylon-6 (C6H11NO)");

    // 7-membered ring with amide
    float r = 1.5f;
    addAtom(mol, r*cosf(0), r*sinf(0), 0.0f, ATOM_C);         // 0 C=O
    addAtom(mol, r*cosf(0.9f), r*sinf(0.9f), 0.0f, ATOM_N);   // 1 NH
    addAtom(mol, r*cosf(1.8f), r*sinf(1.8f), 0.0f, ATOM_C);   // 2
    addAtom(mol, r*cosf(2.7f), r*sinf(2.7f), 0.0f, ATOM_C);   // 3
    addAtom(mol, r*cosf(3.6f), r*sinf(3.6f), 0.0f, ATOM_C);   // 4
    addAtom(mol, r*cosf(4.5f), r*sinf(4.5f), 0.0f, ATOM_C);   // 5
    addAtom(mol, r*cosf(5.4f), r*sinf(5.4f), 0.0f, ATOM_C);   // 6
    addAtom(mol, r*cosf(0)+0.8f, r*sinf(0)+0.8f, 0.0f, ATOM_O); // 7 =O
    addAtom(mol, r*cosf(0.9f)+0.6f, r*sinf(0.9f)+0.6f, 0.0f, ATOM_H); // 8 NH

    addBond(mol, 0, 1, 1);  // C-N
    addBond(mol, 1, 2, 1);  // N-C
    addBond(mol, 2, 3, 1);
    addBond(mol, 3, 4, 1);
    addBond(mol, 4, 5, 1);
    addBond(mol, 5, 6, 1);
    addBond(mol, 6, 0, 1);  // close ring
    addBond(mol, 0, 7, 2);  // C=O
    addBond(mol, 1, 8, 1);  // N-H

    centerMolecule(mol);
}

// Build Acrylonitrile (C3H3N) - ABS component
void buildAcrylonitrile(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Acrylonitrile/ABS (C3H3N)");

    // CH2=CH-CN
    addAtom(mol, -1.2f, 0.0f, 0.0f, ATOM_C);    // CH2=
    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_C);     // =CH
    addAtom(mol, 1.2f, 0.0f, 0.0f, ATOM_C);     // C
    addAtom(mol, 2.3f, 0.0f, 0.0f, ATOM_N);     // N (triple bond)
    addAtom(mol, -1.7f, 0.9f, 0.0f, ATOM_H);
    addAtom(mol, -1.7f, -0.9f, 0.0f, ATOM_H);
    addAtom(mol, 0.0f, 1.0f, 0.0f, ATOM_H);

    addBond(mol, 0, 1, 2);  // C=C
    addBond(mol, 1, 2, 1);  // C-C
    addBond(mol, 2, 3, 3);  // C≡N
    addBond(mol, 0, 4, 1);
    addBond(mol, 0, 5, 1);
    addBond(mol, 1, 6, 1);

    centerMolecule(mol);
}

// Build 1,3-Butadiene (C4H6) - ABS/rubber component
void buildButadiene(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "1,3-Butadiene/ABS (C4H6)");

    // CH2=CH-CH=CH2
    addAtom(mol, -1.8f, 0.0f, 0.0f, ATOM_C);    // CH2=
    addAtom(mol, -0.6f, 0.0f, 0.0f, ATOM_C);    // =CH
    addAtom(mol, 0.6f, 0.0f, 0.0f, ATOM_C);     // CH=
    addAtom(mol, 1.8f, 0.0f, 0.0f, ATOM_C);     // =CH2
    addAtom(mol, -2.3f, 0.9f, 0.0f, ATOM_H);
    addAtom(mol, -2.3f, -0.9f, 0.0f, ATOM_H);
    addAtom(mol, -0.6f, 1.0f, 0.0f, ATOM_H);
    addAtom(mol, 0.6f, 1.0f, 0.0f, ATOM_H);
    addAtom(mol, 2.3f, 0.9f, 0.0f, ATOM_H);
    addAtom(mol, 2.3f, -0.9f, 0.0f, ATOM_H);

    addBond(mol, 0, 1, 2);
    addBond(mol, 1, 2, 1);
    addBond(mol, 2, 3, 2);
    addBond(mol, 0, 4, 1);
    addBond(mol, 0, 5, 1);
    addBond(mol, 1, 6, 1);
    addBond(mol, 2, 7, 1);
    addBond(mol, 3, 8, 1);
    addBond(mol, 3, 9, 1);

    centerMolecule(mol);
}

// Build Bisphenol A (C15H16O2) - Polycarbonate/epoxy component
void buildBisphenolA(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Bisphenol A/BPA (C15H16O2)");

    // Two phenol rings connected by C(CH3)2
    // Left phenol ring
    addAtom(mol, -3.0f, 0.0f, 0.0f, ATOM_C);   // 0
    addAtom(mol, -3.6f, 1.2f, 0.0f, ATOM_C);   // 1
    addAtom(mol, -3.0f, 2.4f, 0.0f, ATOM_C);   // 2
    addAtom(mol, -1.8f, 2.4f, 0.0f, ATOM_C);   // 3
    addAtom(mol, -1.2f, 1.2f, 0.0f, ATOM_C);   // 4
    addAtom(mol, -1.8f, 0.0f, 0.0f, ATOM_C);   // 5
    addAtom(mol, -3.6f, 3.6f, 0.0f, ATOM_O);   // 6 OH
    addAtom(mol, -4.4f, 3.6f, 0.0f, ATOM_H);   // 7

    // Central carbon with 2 CH3
    addAtom(mol, 0.0f, 1.2f, 0.0f, ATOM_C);    // 8 central C
    addAtom(mol, 0.0f, 2.4f, 0.9f, ATOM_C);    // 9 CH3
    addAtom(mol, 0.0f, 2.4f, -0.9f, ATOM_C);   // 10 CH3

    // Right phenol ring
    addAtom(mol, 1.2f, 1.2f, 0.0f, ATOM_C);    // 11
    addAtom(mol, 1.8f, 0.0f, 0.0f, ATOM_C);    // 12
    addAtom(mol, 3.0f, 0.0f, 0.0f, ATOM_C);    // 13
    addAtom(mol, 3.6f, 1.2f, 0.0f, ATOM_C);    // 14
    addAtom(mol, 3.0f, 2.4f, 0.0f, ATOM_C);    // 15
    addAtom(mol, 1.8f, 2.4f, 0.0f, ATOM_C);    // 16
    addAtom(mol, 3.6f, 3.6f, 0.0f, ATOM_O);    // 17 OH
    addAtom(mol, 4.4f, 3.6f, 0.0f, ATOM_H);    // 18

    // Bonds for left ring
    addBond(mol, 0, 1, 2);
    addBond(mol, 1, 2, 1);
    addBond(mol, 2, 3, 2);
    addBond(mol, 3, 4, 1);
    addBond(mol, 4, 5, 2);
    addBond(mol, 5, 0, 1);
    addBond(mol, 2, 6, 1);
    addBond(mol, 6, 7, 1);

    // Central connections
    addBond(mol, 4, 8, 1);
    addBond(mol, 8, 9, 1);
    addBond(mol, 8, 10, 1);
    addBond(mol, 8, 11, 1);

    // Bonds for right ring
    addBond(mol, 11, 12, 2);
    addBond(mol, 12, 13, 1);
    addBond(mol, 13, 14, 2);
    addBond(mol, 14, 15, 1);
    addBond(mol, 15, 16, 2);
    addBond(mol, 16, 11, 1);
    addBond(mol, 15, 17, 1);
    addBond(mol, 17, 18, 1);

    centerMolecule(mol);
}

// Build Terephthalic Acid (C8H6O4) - PET component
void buildTerephthalicAcid(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Terephthalic Acid/PET (C8H6O4)");

    // Benzene ring with -COOH at 1 and 4 positions
    float r = 1.4f;
    for (int i = 0; i < 6; i++) {
        float angle = i * 3.14159f / 3.0f;
        addAtom(mol, r*cosf(angle), r*sinf(angle), 0.0f, ATOM_C);
    }
    // COOH at position 0
    addAtom(mol, 2.5f, 0.0f, 0.0f, ATOM_C);     // 6 COOH
    addAtom(mol, 3.2f, 1.0f, 0.0f, ATOM_O);     // 7 =O
    addAtom(mol, 3.2f, -1.0f, 0.0f, ATOM_O);    // 8 OH
    addAtom(mol, 4.0f, -1.0f, 0.0f, ATOM_H);    // 9
    // COOH at position 3
    addAtom(mol, -2.5f, 0.0f, 0.0f, ATOM_C);    // 10 COOH
    addAtom(mol, -3.2f, 1.0f, 0.0f, ATOM_O);    // 11 =O
    addAtom(mol, -3.2f, -1.0f, 0.0f, ATOM_O);   // 12 OH
    addAtom(mol, -4.0f, -1.0f, 0.0f, ATOM_H);   // 13
    // Ring H
    addAtom(mol, r*cosf(3.14159f/3.0f)+0.6f, r*sinf(3.14159f/3.0f)+0.6f, 0.0f, ATOM_H);
    addAtom(mol, r*cosf(2*3.14159f/3.0f)-0.6f, r*sinf(2*3.14159f/3.0f)+0.6f, 0.0f, ATOM_H);
    addAtom(mol, r*cosf(4*3.14159f/3.0f)-0.6f, r*sinf(4*3.14159f/3.0f)-0.6f, 0.0f, ATOM_H);
    addAtom(mol, r*cosf(5*3.14159f/3.0f)+0.6f, r*sinf(5*3.14159f/3.0f)-0.6f, 0.0f, ATOM_H);

    // Ring bonds (alternating)
    addBond(mol, 0, 1, 2);
    addBond(mol, 1, 2, 1);
    addBond(mol, 2, 3, 2);
    addBond(mol, 3, 4, 1);
    addBond(mol, 4, 5, 2);
    addBond(mol, 5, 0, 1);
    // COOH bonds
    addBond(mol, 0, 6, 1);
    addBond(mol, 6, 7, 2);
    addBond(mol, 6, 8, 1);
    addBond(mol, 8, 9, 1);
    addBond(mol, 3, 10, 1);
    addBond(mol, 10, 11, 2);
    addBond(mol, 10, 12, 1);
    addBond(mol, 12, 13, 1);
    // H bonds
    addBond(mol, 1, 14, 1);
    addBond(mol, 2, 15, 1);
    addBond(mol, 4, 16, 1);
    addBond(mol, 5, 17, 1);

    centerMolecule(mol);
}

// Build Vinyl Acetate (C4H6O2) - EVA component
void buildVinylAcetate(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Vinyl Acetate/EVA (C4H6O2)");

    // CH2=CH-O-CO-CH3
    addAtom(mol, -2.0f, 0.0f, 0.0f, ATOM_C);    // CH2=
    addAtom(mol, -0.8f, 0.0f, 0.0f, ATOM_C);    // =CH
    addAtom(mol, 0.2f, -0.8f, 0.0f, ATOM_O);    // O
    addAtom(mol, 1.4f, -0.3f, 0.0f, ATOM_C);    // C=O
    addAtom(mol, 1.4f, 1.0f, 0.0f, ATOM_O);     // =O
    addAtom(mol, 2.6f, -1.0f, 0.0f, ATOM_C);    // CH3
    addAtom(mol, -2.5f, 0.9f, 0.0f, ATOM_H);
    addAtom(mol, -2.5f, -0.9f, 0.0f, ATOM_H);
    addAtom(mol, -0.8f, 1.0f, 0.0f, ATOM_H);
    addAtom(mol, 3.1f, -0.5f, 0.9f, ATOM_H);
    addAtom(mol, 3.1f, -0.5f, -0.9f, ATOM_H);
    addAtom(mol, 2.6f, -2.0f, 0.0f, ATOM_H);

    addBond(mol, 0, 1, 2);
    addBond(mol, 1, 2, 1);
    addBond(mol, 2, 3, 1);
    addBond(mol, 3, 4, 2);
    addBond(mol, 3, 5, 1);
    addBond(mol, 0, 6, 1);
    addBond(mol, 0, 7, 1);
    addBond(mol, 1, 8, 1);
    addBond(mol, 5, 9, 1);
    addBond(mol, 5, 10, 1);
    addBond(mol, 5, 11, 1);

    centerMolecule(mol);
}

// Build Dimethylsiloxane (C2H6OSi) - Silicone unit
void buildDimethylsiloxane(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "PDMS/Silicone (C2H6OSi)");

    // (CH3)2-Si-O unit
    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_SI);    // Si
    addAtom(mol, 1.8f, 0.0f, 0.0f, ATOM_O);     // O
    addAtom(mol, -1.0f, 1.2f, 0.0f, ATOM_C);    // CH3
    addAtom(mol, -1.0f, -1.2f, 0.0f, ATOM_C);   // CH3
    // H on CH3 groups
    addAtom(mol, -0.5f, 1.7f, 0.9f, ATOM_H);
    addAtom(mol, -0.5f, 1.7f, -0.9f, ATOM_H);
    addAtom(mol, -2.0f, 1.5f, 0.0f, ATOM_H);
    addAtom(mol, -0.5f, -1.7f, 0.9f, ATOM_H);
    addAtom(mol, -0.5f, -1.7f, -0.9f, ATOM_H);
    addAtom(mol, -2.0f, -1.5f, 0.0f, ATOM_H);

    addBond(mol, 0, 1, 1);  // Si-O
    addBond(mol, 0, 2, 1);  // Si-C
    addBond(mol, 0, 3, 1);  // Si-C
    addBond(mol, 2, 4, 1);
    addBond(mol, 2, 5, 1);
    addBond(mol, 2, 6, 1);
    addBond(mol, 3, 7, 1);
    addBond(mol, 3, 8, 1);
    addBond(mol, 3, 9, 1);

    centerMolecule(mol);
}

// Build Adipic Acid (C6H10O4) - Nylon-6,6 component
void buildAdipicAcid(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Adipic Acid/Nylon-6,6 (C6H10O4)");

    // HOOC-(CH2)4-COOH
    addAtom(mol, -3.5f, 0.0f, 0.0f, ATOM_C);    // COOH
    addAtom(mol, -4.2f, 1.0f, 0.0f, ATOM_O);    // =O
    addAtom(mol, -4.2f, -1.0f, 0.0f, ATOM_O);   // OH
    addAtom(mol, -5.0f, -1.0f, 0.0f, ATOM_H);
    addAtom(mol, -2.0f, 0.0f, 0.0f, ATOM_C);    // CH2
    addAtom(mol, -0.7f, 0.0f, 0.0f, ATOM_C);    // CH2
    addAtom(mol, 0.7f, 0.0f, 0.0f, ATOM_C);     // CH2
    addAtom(mol, 2.0f, 0.0f, 0.0f, ATOM_C);     // CH2
    addAtom(mol, 3.5f, 0.0f, 0.0f, ATOM_C);     // COOH
    addAtom(mol, 4.2f, 1.0f, 0.0f, ATOM_O);     // =O
    addAtom(mol, 4.2f, -1.0f, 0.0f, ATOM_O);    // OH
    addAtom(mol, 5.0f, -1.0f, 0.0f, ATOM_H);
    // CH2 hydrogens
    addAtom(mol, -2.0f, 0.5f, 0.9f, ATOM_H);
    addAtom(mol, -2.0f, 0.5f, -0.9f, ATOM_H);
    addAtom(mol, -0.7f, 0.5f, 0.9f, ATOM_H);
    addAtom(mol, -0.7f, 0.5f, -0.9f, ATOM_H);
    addAtom(mol, 0.7f, 0.5f, 0.9f, ATOM_H);
    addAtom(mol, 0.7f, 0.5f, -0.9f, ATOM_H);
    addAtom(mol, 2.0f, 0.5f, 0.9f, ATOM_H);
    addAtom(mol, 2.0f, 0.5f, -0.9f, ATOM_H);

    addBond(mol, 0, 1, 2);
    addBond(mol, 0, 2, 1);
    addBond(mol, 2, 3, 1);
    addBond(mol, 0, 4, 1);
    addBond(mol, 4, 5, 1);
    addBond(mol, 5, 6, 1);
    addBond(mol, 6, 7, 1);
    addBond(mol, 7, 8, 1);
    addBond(mol, 8, 9, 2);
    addBond(mol, 8, 10, 1);
    addBond(mol, 10, 11, 1);
    addBond(mol, 4, 12, 1);
    addBond(mol, 4, 13, 1);
    addBond(mol, 5, 14, 1);
    addBond(mol, 5, 15, 1);
    addBond(mol, 6, 16, 1);
    addBond(mol, 6, 17, 1);
    addBond(mol, 7, 18, 1);
    addBond(mol, 7, 19, 1);

    centerMolecule(mol);
}

// Build Hexamethylenediamine (C6H16N2) - Nylon-6,6 component
void buildHexamethylenediamine(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "HMDA/Nylon-6,6 (C6H16N2)");

    // H2N-(CH2)6-NH2
    addAtom(mol, -4.5f, 0.0f, 0.0f, ATOM_N);    // NH2
    addAtom(mol, -3.0f, 0.0f, 0.0f, ATOM_C);    // CH2
    addAtom(mol, -1.5f, 0.0f, 0.0f, ATOM_C);    // CH2
    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_C);     // CH2
    addAtom(mol, 1.5f, 0.0f, 0.0f, ATOM_C);     // CH2
    addAtom(mol, 3.0f, 0.0f, 0.0f, ATOM_C);     // CH2
    addAtom(mol, 4.5f, 0.0f, 0.0f, ATOM_C);     // CH2
    addAtom(mol, 6.0f, 0.0f, 0.0f, ATOM_N);     // NH2
    // NH2 hydrogens
    addAtom(mol, -5.0f, 0.8f, 0.0f, ATOM_H);
    addAtom(mol, -5.0f, -0.8f, 0.0f, ATOM_H);
    addAtom(mol, 6.5f, 0.8f, 0.0f, ATOM_H);
    addAtom(mol, 6.5f, -0.8f, 0.0f, ATOM_H);

    addBond(mol, 0, 1, 1);
    addBond(mol, 1, 2, 1);
    addBond(mol, 2, 3, 1);
    addBond(mol, 3, 4, 1);
    addBond(mol, 4, 5, 1);
    addBond(mol, 5, 6, 1);
    addBond(mol, 6, 7, 1);
    addBond(mol, 0, 8, 1);
    addBond(mol, 0, 9, 1);
    addBond(mol, 7, 10, 1);
    addBond(mol, 7, 11, 1);

    centerMolecule(mol);
}

// Build Isoprene (C5H8) - Natural rubber monomer
void buildIsoprene(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Isoprene/Rubber (C5H8)");

    // CH2=C(CH3)-CH=CH2
    addAtom(mol, -2.0f, 0.0f, 0.0f, ATOM_C);    // CH2=
    addAtom(mol, -0.7f, 0.0f, 0.0f, ATOM_C);    // =C
    addAtom(mol, -0.7f, 1.4f, 0.0f, ATOM_C);    // CH3
    addAtom(mol, 0.5f, -0.7f, 0.0f, ATOM_C);    // CH=
    addAtom(mol, 1.8f, -0.2f, 0.0f, ATOM_C);    // =CH2
    addAtom(mol, -2.5f, 0.9f, 0.0f, ATOM_H);
    addAtom(mol, -2.5f, -0.9f, 0.0f, ATOM_H);
    addAtom(mol, -0.2f, 1.9f, 0.9f, ATOM_H);
    addAtom(mol, -0.2f, 1.9f, -0.9f, ATOM_H);
    addAtom(mol, -1.7f, 1.8f, 0.0f, ATOM_H);
    addAtom(mol, 0.5f, -1.7f, 0.0f, ATOM_H);
    addAtom(mol, 2.3f, 0.7f, 0.0f, ATOM_H);
    addAtom(mol, 2.3f, -1.1f, 0.0f, ATOM_H);

    addBond(mol, 0, 1, 2);
    addBond(mol, 1, 2, 1);
    addBond(mol, 1, 3, 1);
    addBond(mol, 3, 4, 2);
    addBond(mol, 0, 5, 1);
    addBond(mol, 0, 6, 1);
    addBond(mol, 2, 7, 1);
    addBond(mol, 2, 8, 1);
    addBond(mol, 2, 9, 1);
    addBond(mol, 3, 10, 1);
    addBond(mol, 4, 11, 1);
    addBond(mol, 4, 12, 1);

    centerMolecule(mol);
}

// Build Polyethylene dimer (C4H10) - PE structure
void buildPEDimer(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "PE Dimer/HDPE (C4H10)");

    // Actually butane - shows PE repeat
    addAtom(mol, -1.5f, 0.0f, 0.0f, ATOM_C);
    addAtom(mol, -0.5f, 0.0f, 0.0f, ATOM_C);
    addAtom(mol, 0.5f, 0.0f, 0.0f, ATOM_C);
    addAtom(mol, 1.5f, 0.0f, 0.0f, ATOM_C);
    // H atoms
    addAtom(mol, -2.0f, 0.5f, 0.9f, ATOM_H);
    addAtom(mol, -2.0f, 0.5f, -0.9f, ATOM_H);
    addAtom(mol, -2.0f, -0.9f, 0.0f, ATOM_H);
    addAtom(mol, -0.5f, 0.5f, 0.9f, ATOM_H);
    addAtom(mol, -0.5f, -0.5f, -0.9f, ATOM_H);
    addAtom(mol, 0.5f, 0.5f, 0.9f, ATOM_H);
    addAtom(mol, 0.5f, -0.5f, -0.9f, ATOM_H);
    addAtom(mol, 2.0f, 0.5f, 0.9f, ATOM_H);
    addAtom(mol, 2.0f, 0.5f, -0.9f, ATOM_H);
    addAtom(mol, 2.0f, -0.9f, 0.0f, ATOM_H);

    addBond(mol, 0, 1, 1);
    addBond(mol, 1, 2, 1);
    addBond(mol, 2, 3, 1);
    addBond(mol, 0, 4, 1);
    addBond(mol, 0, 5, 1);
    addBond(mol, 0, 6, 1);
    addBond(mol, 1, 7, 1);
    addBond(mol, 1, 8, 1);
    addBond(mol, 2, 9, 1);
    addBond(mol, 2, 10, 1);
    addBond(mol, 3, 11, 1);
    addBond(mol, 3, 12, 1);
    addBond(mol, 3, 13, 1);

    centerMolecule(mol);
}

// Build Polypropylene dimer (C6H14) - PP structure
void buildPPDimer(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "PP Dimer/Polypropylene (C6H14)");

    // Two propylene units
    addAtom(mol, -1.5f, 0.0f, 0.0f, ATOM_C);    // CH3
    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_C);     // CH
    addAtom(mol, 0.0f, -1.4f, 0.0f, ATOM_C);    // CH3 branch
    addAtom(mol, 1.5f, 0.0f, 0.0f, ATOM_C);     // CH2
    addAtom(mol, 3.0f, 0.0f, 0.0f, ATOM_C);     // CH
    addAtom(mol, 3.0f, -1.4f, 0.0f, ATOM_C);    // CH3 branch

    addBond(mol, 0, 1, 1);
    addBond(mol, 1, 2, 1);
    addBond(mol, 1, 3, 1);
    addBond(mol, 3, 4, 1);
    addBond(mol, 4, 5, 1);

    centerMolecule(mol);
}

// Build PVC Trimer - shows the chlorine pattern
void buildPVCTrimer(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "PVC Trimer (-CHCl-CH2-)3");

    // -CHCl-CH2-CHCl-CH2-CHCl-CH2-
    addAtom(mol, -3.0f, 0.0f, 0.0f, ATOM_C);    // CHCl
    addAtom(mol, -3.0f, 1.5f, 0.0f, ATOM_CL);
    addAtom(mol, -1.5f, 0.0f, 0.0f, ATOM_C);    // CH2
    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_C);     // CHCl
    addAtom(mol, 0.0f, 1.5f, 0.0f, ATOM_CL);
    addAtom(mol, 1.5f, 0.0f, 0.0f, ATOM_C);     // CH2
    addAtom(mol, 3.0f, 0.0f, 0.0f, ATOM_C);     // CHCl
    addAtom(mol, 3.0f, 1.5f, 0.0f, ATOM_CL);

    addBond(mol, 0, 1, 1);
    addBond(mol, 0, 2, 1);
    addBond(mol, 2, 3, 1);
    addBond(mol, 3, 4, 1);
    addBond(mol, 3, 5, 1);
    addBond(mol, 5, 6, 1);
    addBond(mol, 6, 7, 1);

    centerMolecule(mol);
}

// Build Polystyrene dimer - shows benzene pendant groups
void buildPSDimer(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "PS Dimer/Polystyrene");

    // Backbone
    addAtom(mol, -1.5f, 0.0f, 0.0f, ATOM_C);    // CH
    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_C);     // CH2
    addAtom(mol, 1.5f, 0.0f, 0.0f, ATOM_C);     // CH
    addAtom(mol, 3.0f, 0.0f, 0.0f, ATOM_C);     // CH3

    // First benzene ring (pendant on atom 0)
    float r = 1.2f;
    addAtom(mol, -1.5f + r*cosf(3.14159f/2), -r - r*sinf(3.14159f/2), 0.0f, ATOM_C);
    addAtom(mol, -1.5f + r*cosf(3.14159f/6), -r - r*sinf(3.14159f/6), 0.0f, ATOM_C);
    addAtom(mol, -1.5f + r*cosf(-3.14159f/6), -r - r*sinf(-3.14159f/6), 0.0f, ATOM_C);
    addAtom(mol, -1.5f + r*cosf(-3.14159f/2), -r - r*sinf(-3.14159f/2), 0.0f, ATOM_C);
    addAtom(mol, -1.5f + r*cosf(-5*3.14159f/6), -r - r*sinf(-5*3.14159f/6), 0.0f, ATOM_C);
    addAtom(mol, -1.5f + r*cosf(5*3.14159f/6), -r - r*sinf(5*3.14159f/6), 0.0f, ATOM_C);

    // Second benzene ring (pendant on atom 2)
    addAtom(mol, 1.5f + r*cosf(3.14159f/2), -r - r*sinf(3.14159f/2), 0.0f, ATOM_C);
    addAtom(mol, 1.5f + r*cosf(3.14159f/6), -r - r*sinf(3.14159f/6), 0.0f, ATOM_C);
    addAtom(mol, 1.5f + r*cosf(-3.14159f/6), -r - r*sinf(-3.14159f/6), 0.0f, ATOM_C);
    addAtom(mol, 1.5f + r*cosf(-3.14159f/2), -r - r*sinf(-3.14159f/2), 0.0f, ATOM_C);
    addAtom(mol, 1.5f + r*cosf(-5*3.14159f/6), -r - r*sinf(-5*3.14159f/6), 0.0f, ATOM_C);
    addAtom(mol, 1.5f + r*cosf(5*3.14159f/6), -r - r*sinf(5*3.14159f/6), 0.0f, ATOM_C);

    // Backbone bonds
    addBond(mol, 0, 1, 1);
    addBond(mol, 1, 2, 1);
    addBond(mol, 2, 3, 1);

    // First ring bonds
    addBond(mol, 0, 4, 1);  // Connect to backbone
    addBond(mol, 4, 5, 2);
    addBond(mol, 5, 6, 1);
    addBond(mol, 6, 7, 2);
    addBond(mol, 7, 8, 1);
    addBond(mol, 8, 9, 2);
    addBond(mol, 9, 4, 1);

    // Second ring bonds
    addBond(mol, 2, 10, 1);  // Connect to backbone
    addBond(mol, 10, 11, 2);
    addBond(mol, 11, 12, 1);
    addBond(mol, 12, 13, 2);
    addBond(mol, 13, 14, 1);
    addBond(mol, 14, 15, 2);
    addBond(mol, 15, 10, 1);

    centerMolecule(mol);
}

// Build PTFE Trimer - Teflon chain segment
void buildPTFETrimer(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "PTFE Trimer/Teflon (-CF2-)6");

    // Six CF2 units
    for (int i = 0; i < 6; i++) {
        float x = i * 1.3f;
        addAtom(mol, x, 0.0f, 0.0f, ATOM_C);
        addAtom(mol, x, 0.9f, 0.5f, ATOM_F);
        addAtom(mol, x, -0.9f, -0.5f, ATOM_F);
    }

    // C-C bonds
    for (int i = 0; i < 5; i++) {
        addBond(mol, i*3, (i+1)*3, 1);
    }
    // C-F bonds
    for (int i = 0; i < 6; i++) {
        addBond(mol, i*3, i*3+1, 1);
        addBond(mol, i*3, i*3+2, 1);
    }

    centerMolecule(mol);
}

// ============== TRANS FATS ==============

// Build Elaidic Acid (C18H34O2) - main industrial trans fat
void buildElaidicAcid(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Elaidic Acid/trans-C18:1");

    // 18-carbon chain with trans double bond at C9
    // COOH at one end
    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_C);     // 0 COOH
    addAtom(mol, 0.7f, 1.0f, 0.0f, ATOM_O);     // 1 =O
    addAtom(mol, 0.7f, -1.0f, 0.0f, ATOM_O);    // 2 OH
    addAtom(mol, 1.5f, -1.0f, 0.0f, ATOM_H);    // 3

    // Carbon chain C2-C8
    float x = 1.3f;
    for (int i = 0; i < 7; i++) {
        addAtom(mol, x + i*1.2f, (i%2)*0.4f, 0.0f, ATOM_C);  // 4-10
    }

    // Trans double bond C9=C10
    addAtom(mol, x + 7*1.2f, 0.0f, 0.0f, ATOM_C);      // 11 C9
    addAtom(mol, x + 7*1.2f + 1.3f, 0.4f, 0.0f, ATOM_C); // 12 C10 (trans = same side)

    // Carbon chain C11-C18
    for (int i = 0; i < 8; i++) {
        addAtom(mol, x + 8.5f*1.2f + i*1.2f, (i%2)*0.4f + 0.4f, 0.0f, ATOM_C);  // 13-20
    }

    // Bonds
    addBond(mol, 0, 1, 2);   // C=O
    addBond(mol, 0, 2, 1);   // C-OH
    addBond(mol, 2, 3, 1);   // O-H
    addBond(mol, 0, 4, 1);   // COOH to chain
    for (int i = 4; i < 11; i++) {
        addBond(mol, i, i+1, 1);
    }
    addBond(mol, 11, 12, 2);  // Trans C=C
    for (int i = 12; i < 20; i++) {
        addBond(mol, i, i+1, 1);
    }

    centerMolecule(mol);
}

// Build Oleic Acid (C18H34O2) - cis fat (for comparison)
void buildOleicAcid(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Oleic Acid/cis-C18:1");

    // Similar to elaidic but cis configuration causes bend
    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_C);     // COOH
    addAtom(mol, 0.7f, 1.0f, 0.0f, ATOM_O);     // =O
    addAtom(mol, 0.7f, -1.0f, 0.0f, ATOM_O);    // OH
    addAtom(mol, 1.5f, -1.0f, 0.0f, ATOM_H);

    // Chain up to double bond
    float x = 1.3f;
    for (int i = 0; i < 7; i++) {
        addAtom(mol, x + i*1.2f, (i%2)*0.3f, 0.0f, ATOM_C);
    }

    // Cis double bond - creates 30 degree bend
    addAtom(mol, x + 7*1.2f, 0.0f, 0.0f, ATOM_C);
    addAtom(mol, x + 7*1.2f + 1.0f, -0.8f, 0.0f, ATOM_C);  // Cis = opposite side

    // Chain continues at angle
    for (int i = 0; i < 8; i++) {
        addAtom(mol, x + 8*1.2f + i*1.0f, -1.2f - (i%2)*0.3f, 0.0f, ATOM_C);
    }

    // Bonds
    addBond(mol, 0, 1, 2);
    addBond(mol, 0, 2, 1);
    addBond(mol, 2, 3, 1);
    addBond(mol, 0, 4, 1);
    for (int i = 4; i < 11; i++) {
        addBond(mol, i, i+1, 1);
    }
    addBond(mol, 11, 12, 2);  // Cis C=C
    for (int i = 12; i < 20; i++) {
        addBond(mol, i, i+1, 1);
    }

    centerMolecule(mol);
}

// Build Vaccenic Acid (C18H34O2) - natural trans fat (ruminant)
void buildVaccenicAcid(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Vaccenic Acid/trans-11");

    // trans-11-octadecenoic acid
    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_C);
    addAtom(mol, 0.7f, 1.0f, 0.0f, ATOM_O);
    addAtom(mol, 0.7f, -1.0f, 0.0f, ATOM_O);
    addAtom(mol, 1.5f, -1.0f, 0.0f, ATOM_H);

    // Chain C2-C10
    float x = 1.3f;
    for (int i = 0; i < 9; i++) {
        addAtom(mol, x + i*1.2f, (i%2)*0.4f, 0.0f, ATOM_C);
    }

    // Trans double bond C11=C12
    addAtom(mol, x + 9*1.2f, 0.0f, 0.0f, ATOM_C);
    addAtom(mol, x + 9*1.2f + 1.3f, 0.4f, 0.0f, ATOM_C);

    // Chain C13-C18
    for (int i = 0; i < 6; i++) {
        addAtom(mol, x + 10.5f*1.2f + i*1.2f, (i%2)*0.4f + 0.4f, 0.0f, ATOM_C);
    }

    // Bonds
    addBond(mol, 0, 1, 2);
    addBond(mol, 0, 2, 1);
    addBond(mol, 2, 3, 1);
    addBond(mol, 0, 4, 1);
    for (int i = 4; i < 13; i++) {
        addBond(mol, i, i+1, 1);
    }
    addBond(mol, 13, 14, 2);  // Trans C=C
    for (int i = 14; i < 19; i++) {
        addBond(mol, i, i+1, 1);
    }

    centerMolecule(mol);
}

// Build trans-Palmitoleic Acid (C16H30O2)
void buildTransPalmitoleicAcid(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "trans-Palmitoleic (C16:1)");

    // trans-9-hexadecenoic acid
    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_C);
    addAtom(mol, 0.7f, 1.0f, 0.0f, ATOM_O);
    addAtom(mol, 0.7f, -1.0f, 0.0f, ATOM_O);
    addAtom(mol, 1.5f, -1.0f, 0.0f, ATOM_H);

    // Chain C2-C8
    float x = 1.3f;
    for (int i = 0; i < 7; i++) {
        addAtom(mol, x + i*1.2f, (i%2)*0.4f, 0.0f, ATOM_C);
    }

    // Trans double bond C9=C10
    addAtom(mol, x + 7*1.2f, 0.0f, 0.0f, ATOM_C);
    addAtom(mol, x + 7*1.2f + 1.3f, 0.4f, 0.0f, ATOM_C);

    // Chain C11-C16
    for (int i = 0; i < 6; i++) {
        addAtom(mol, x + 8.5f*1.2f + i*1.2f, (i%2)*0.4f + 0.4f, 0.0f, ATOM_C);
    }

    // Bonds
    addBond(mol, 0, 1, 2);
    addBond(mol, 0, 2, 1);
    addBond(mol, 2, 3, 1);
    addBond(mol, 0, 4, 1);
    for (int i = 4; i < 11; i++) {
        addBond(mol, i, i+1, 1);
    }
    addBond(mol, 11, 12, 2);
    for (int i = 12; i < 17; i++) {
        addBond(mol, i, i+1, 1);
    }

    centerMolecule(mol);
}

// Build Linoelaidic Acid (C18H32O2) - trans,trans-linoleic
void buildLinoelaidicAcid(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Linoelaidic/trans,trans-C18:2");

    // trans,trans-9,12-octadecadienoic acid
    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_C);
    addAtom(mol, 0.7f, 1.0f, 0.0f, ATOM_O);
    addAtom(mol, 0.7f, -1.0f, 0.0f, ATOM_O);
    addAtom(mol, 1.5f, -1.0f, 0.0f, ATOM_H);

    // Chain to first double bond
    float x = 1.3f;
    for (int i = 0; i < 7; i++) {
        addAtom(mol, x + i*1.2f, (i%2)*0.3f, 0.0f, ATOM_C);
    }

    // First trans double bond C9=C10
    addAtom(mol, x + 7*1.2f, 0.0f, 0.0f, ATOM_C);
    addAtom(mol, x + 8*1.2f, 0.3f, 0.0f, ATOM_C);

    // CH2 at C11
    addAtom(mol, x + 9*1.2f, 0.0f, 0.0f, ATOM_C);

    // Second trans double bond C12=C13
    addAtom(mol, x + 10*1.2f, 0.3f, 0.0f, ATOM_C);
    addAtom(mol, x + 11*1.2f, 0.0f, 0.0f, ATOM_C);

    // Rest of chain C14-C18
    for (int i = 0; i < 5; i++) {
        addAtom(mol, x + (12+i)*1.2f, (i%2)*0.3f, 0.0f, ATOM_C);
    }

    // Bonds
    addBond(mol, 0, 1, 2);
    addBond(mol, 0, 2, 1);
    addBond(mol, 2, 3, 1);
    addBond(mol, 0, 4, 1);
    for (int i = 4; i < 11; i++) {
        addBond(mol, i, i+1, 1);
    }
    addBond(mol, 11, 12, 2);  // First trans
    addBond(mol, 12, 13, 1);
    addBond(mol, 13, 14, 1);
    addBond(mol, 14, 15, 2);  // Second trans
    for (int i = 15; i < 20; i++) {
        addBond(mol, i, i+1, 1);
    }

    centerMolecule(mol);
}

// Build Rumenic Acid (C18H32O2) - CLA, cis-9,trans-11
void buildRumenicAcid(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Rumenic Acid/CLA c9,t11");

    // Conjugated linoleic acid: cis-9, trans-11
    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_C);
    addAtom(mol, 0.7f, 1.0f, 0.0f, ATOM_O);
    addAtom(mol, 0.7f, -1.0f, 0.0f, ATOM_O);
    addAtom(mol, 1.5f, -1.0f, 0.0f, ATOM_H);

    // Chain C2-C8
    float x = 1.3f;
    for (int i = 0; i < 7; i++) {
        addAtom(mol, x + i*1.2f, (i%2)*0.3f, 0.0f, ATOM_C);
    }

    // Cis double bond C9=C10 (conjugated)
    addAtom(mol, x + 7*1.2f, 0.0f, 0.0f, ATOM_C);
    addAtom(mol, x + 8*1.2f, -0.6f, 0.0f, ATOM_C);  // Cis - bends down

    // Trans double bond C11=C12 (conjugated with C9=C10)
    addAtom(mol, x + 9*1.2f, -0.3f, 0.0f, ATOM_C);
    addAtom(mol, x + 10*1.2f, 0.3f, 0.0f, ATOM_C);  // Trans - continues up

    // Rest of chain C13-C18
    for (int i = 0; i < 6; i++) {
        addAtom(mol, x + (11+i)*1.2f, (i%2)*0.3f + 0.3f, 0.0f, ATOM_C);
    }

    // Bonds
    addBond(mol, 0, 1, 2);
    addBond(mol, 0, 2, 1);
    addBond(mol, 2, 3, 1);
    addBond(mol, 0, 4, 1);
    for (int i = 4; i < 11; i++) {
        addBond(mol, i, i+1, 1);
    }
    addBond(mol, 11, 12, 2);  // Cis C9=C10
    addBond(mol, 12, 13, 1);  // Single bond in conjugated system
    addBond(mol, 13, 14, 2);  // Trans C11=C12
    for (int i = 14; i < 20; i++) {
        addBond(mol, i, i+1, 1);
    }

    centerMolecule(mol);
}

// Build Stearic Acid (C18H36O2) - saturated fat for comparison
void buildStearicAcid(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Stearic Acid/Sat C18:0");

    // Fully saturated 18-carbon fatty acid
    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_C);
    addAtom(mol, 0.7f, 1.0f, 0.0f, ATOM_O);
    addAtom(mol, 0.7f, -1.0f, 0.0f, ATOM_O);
    addAtom(mol, 1.5f, -1.0f, 0.0f, ATOM_H);

    // All 17 CH2/CH3 carbons
    float x = 1.3f;
    for (int i = 0; i < 17; i++) {
        addAtom(mol, x + i*1.2f, (i%2)*0.3f, 0.0f, ATOM_C);
    }

    // All single bonds
    addBond(mol, 0, 1, 2);
    addBond(mol, 0, 2, 1);
    addBond(mol, 2, 3, 1);
    addBond(mol, 0, 4, 1);
    for (int i = 4; i < 20; i++) {
        addBond(mol, i, i+1, 1);
    }

    centerMolecule(mol);
}

// Build Palmitic Acid (C16H32O2) - common saturated fat
void buildPalmiticAcid(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Palmitic Acid/Sat C16:0");

    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_C);
    addAtom(mol, 0.7f, 1.0f, 0.0f, ATOM_O);
    addAtom(mol, 0.7f, -1.0f, 0.0f, ATOM_O);
    addAtom(mol, 1.5f, -1.0f, 0.0f, ATOM_H);

    float x = 1.3f;
    for (int i = 0; i < 15; i++) {
        addAtom(mol, x + i*1.2f, (i%2)*0.3f, 0.0f, ATOM_C);
    }

    addBond(mol, 0, 1, 2);
    addBond(mol, 0, 2, 1);
    addBond(mol, 2, 3, 1);
    addBond(mol, 0, 4, 1);
    for (int i = 4; i < 18; i++) {
        addBond(mol, i, i+1, 1);
    }

    centerMolecule(mol);
}

// Build Linoleic Acid (C18H32O2) - essential omega-6
void buildLinoleicAcid(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Linoleic Acid/cis,cis-C18:2");

    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_C);
    addAtom(mol, 0.7f, 1.0f, 0.0f, ATOM_O);
    addAtom(mol, 0.7f, -1.0f, 0.0f, ATOM_O);
    addAtom(mol, 1.5f, -1.0f, 0.0f, ATOM_H);

    float x = 1.3f;
    // Chain to C9
    for (int i = 0; i < 7; i++) {
        addAtom(mol, x + i*1.2f, (i%2)*0.3f, 0.0f, ATOM_C);
    }
    // Cis C9=C10
    addAtom(mol, x + 7*1.2f, 0.0f, 0.0f, ATOM_C);
    addAtom(mol, x + 8*1.2f, -0.8f, 0.0f, ATOM_C);
    // C11
    addAtom(mol, x + 9*1.2f, -0.5f, 0.0f, ATOM_C);
    // Cis C12=C13
    addAtom(mol, x + 10*1.2f, -1.0f, 0.0f, ATOM_C);
    addAtom(mol, x + 11*1.2f, -1.8f, 0.0f, ATOM_C);
    // C14-C18
    for (int i = 0; i < 5; i++) {
        addAtom(mol, x + (12+i)*1.2f, -1.5f - (i%2)*0.3f, 0.0f, ATOM_C);
    }

    addBond(mol, 0, 1, 2);
    addBond(mol, 0, 2, 1);
    addBond(mol, 2, 3, 1);
    addBond(mol, 0, 4, 1);
    for (int i = 4; i < 11; i++) {
        addBond(mol, i, i+1, 1);
    }
    addBond(mol, 11, 12, 2);  // First cis
    addBond(mol, 12, 13, 1);
    addBond(mol, 13, 14, 1);
    addBond(mol, 14, 15, 2);  // Second cis
    for (int i = 15; i < 20; i++) {
        addBond(mol, i, i+1, 1);
    }

    centerMolecule(mol);
}

// Build Alpha-Linolenic Acid (C18H30O2) - essential omega-3
void buildAlphaLinolenicAcid(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "ALA/Omega-3 C18:3");

    // cis,cis,cis-9,12,15-octadecatrienoic acid
    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_C);
    addAtom(mol, 0.7f, 1.0f, 0.0f, ATOM_O);
    addAtom(mol, 0.7f, -1.0f, 0.0f, ATOM_O);
    addAtom(mol, 1.5f, -1.0f, 0.0f, ATOM_H);

    float x = 1.3f;
    // Chain to C9
    for (int i = 0; i < 7; i++) {
        addAtom(mol, x + i*1.0f, (i%2)*0.3f, 0.0f, ATOM_C);
    }
    // C9=C10 cis
    addAtom(mol, x + 7*1.0f, 0.0f, 0.0f, ATOM_C);
    addAtom(mol, x + 8*1.0f, -0.7f, 0.0f, ATOM_C);
    // C11
    addAtom(mol, x + 9*1.0f, -0.4f, 0.0f, ATOM_C);
    // C12=C13 cis
    addAtom(mol, x + 10*1.0f, -0.9f, 0.0f, ATOM_C);
    addAtom(mol, x + 11*1.0f, -1.6f, 0.0f, ATOM_C);
    // C14
    addAtom(mol, x + 12*1.0f, -1.3f, 0.0f, ATOM_C);
    // C15=C16 cis
    addAtom(mol, x + 13*1.0f, -1.8f, 0.0f, ATOM_C);
    addAtom(mol, x + 14*1.0f, -2.5f, 0.0f, ATOM_C);
    // C17-C18
    addAtom(mol, x + 15*1.0f, -2.2f, 0.0f, ATOM_C);
    addAtom(mol, x + 16*1.0f, -2.5f, 0.0f, ATOM_C);

    addBond(mol, 0, 1, 2);
    addBond(mol, 0, 2, 1);
    addBond(mol, 2, 3, 1);
    addBond(mol, 0, 4, 1);
    for (int i = 4; i < 11; i++) {
        addBond(mol, i, i+1, 1);
    }
    addBond(mol, 11, 12, 2);  // First cis
    addBond(mol, 12, 13, 1);
    addBond(mol, 13, 14, 1);
    addBond(mol, 14, 15, 2);  // Second cis
    addBond(mol, 15, 16, 1);
    addBond(mol, 16, 17, 1);
    addBond(mol, 17, 18, 2);  // Third cis
    addBond(mol, 18, 19, 1);
    addBond(mol, 19, 20, 1);

    centerMolecule(mol);
}

// ============== METAL COMPOUNDS ==============

// Build Rust/Hematite (Fe2O3) - iron(III) oxide
void buildRust(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Rust/Hematite (Fe2O3)");

    // Trigonal structure
    addAtom(mol, -0.8f, 0.0f, 0.0f, ATOM_FE);   // Fe
    addAtom(mol, 0.8f, 0.0f, 0.0f, ATOM_FE);    // Fe
    addAtom(mol, 0.0f, 1.2f, 0.0f, ATOM_O);     // O (bridging)
    addAtom(mol, -1.5f, -1.0f, 0.0f, ATOM_O);   // O
    addAtom(mol, 1.5f, -1.0f, 0.0f, ATOM_O);    // O

    addBond(mol, 0, 2, 1);  // Fe-O
    addBond(mol, 1, 2, 1);  // Fe-O
    addBond(mol, 0, 3, 2);  // Fe=O
    addBond(mol, 1, 4, 2);  // Fe=O

    centerMolecule(mol);
}

// Build Magnetite (Fe3O4) - magnetic iron oxide
void buildMagnetite(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Magnetite (Fe3O4)");

    // Inverse spinel structure (simplified)
    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_FE);    // Fe center
    addAtom(mol, -1.5f, 0.0f, 0.0f, ATOM_FE);   // Fe
    addAtom(mol, 1.5f, 0.0f, 0.0f, ATOM_FE);    // Fe
    addAtom(mol, 0.0f, 1.3f, 0.0f, ATOM_O);     // O
    addAtom(mol, 0.0f, -1.3f, 0.0f, ATOM_O);    // O
    addAtom(mol, -0.75f, 0.0f, 1.1f, ATOM_O);   // O
    addAtom(mol, 0.75f, 0.0f, -1.1f, ATOM_O);   // O

    addBond(mol, 0, 3, 1);
    addBond(mol, 0, 4, 1);
    addBond(mol, 0, 5, 1);
    addBond(mol, 0, 6, 1);
    addBond(mol, 1, 3, 1);
    addBond(mol, 1, 5, 1);
    addBond(mol, 2, 4, 1);
    addBond(mol, 2, 6, 1);

    centerMolecule(mol);
}

// Build Iron Pentacarbonyl (Fe(CO)5) - organometallic
void buildIronPentacarbonyl(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Iron Pentacarbonyl (Fe(CO)5)");

    // Trigonal bipyramidal structure
    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_FE);    // Fe center

    // Equatorial CO (3 at 120 degrees)
    float r = 1.8f;
    for (int i = 0; i < 3; i++) {
        float angle = i * 2.0f * 3.14159f / 3.0f;
        addAtom(mol, r*cosf(angle), r*sinf(angle), 0.0f, ATOM_C);
        addAtom(mol, 1.5f*r*cosf(angle), 1.5f*r*sinf(angle), 0.0f, ATOM_O);
    }

    // Axial CO (top and bottom)
    addAtom(mol, 0.0f, 0.0f, 1.8f, ATOM_C);
    addAtom(mol, 0.0f, 0.0f, 3.0f, ATOM_O);
    addAtom(mol, 0.0f, 0.0f, -1.8f, ATOM_C);
    addAtom(mol, 0.0f, 0.0f, -3.0f, ATOM_O);

    // Fe-C bonds
    addBond(mol, 0, 1, 1);
    addBond(mol, 0, 3, 1);
    addBond(mol, 0, 5, 1);
    addBond(mol, 0, 7, 1);
    addBond(mol, 0, 9, 1);

    // C≡O bonds
    addBond(mol, 1, 2, 3);
    addBond(mol, 3, 4, 3);
    addBond(mol, 5, 6, 3);
    addBond(mol, 7, 8, 3);
    addBond(mol, 9, 10, 3);

    centerMolecule(mol);
}

// Build Copper Sulfate (CuSO4) - blue vitriol
void buildCopperSulfate(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Copper Sulfate (CuSO4)");

    // Cu2+ with SO4 2-
    addAtom(mol, -2.0f, 0.0f, 0.0f, ATOM_CU);   // Cu
    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_S);     // S
    addAtom(mol, 1.2f, 0.0f, 0.0f, ATOM_O);     // O
    addAtom(mol, -0.6f, 1.0f, 0.0f, ATOM_O);    // O (bound to Cu)
    addAtom(mol, -0.6f, -1.0f, 0.0f, ATOM_O);   // O
    addAtom(mol, 0.0f, 0.0f, 1.2f, ATOM_O);     // O

    addBond(mol, 1, 2, 2);  // S=O
    addBond(mol, 1, 3, 1);  // S-O
    addBond(mol, 1, 4, 2);  // S=O
    addBond(mol, 1, 5, 1);  // S-O
    addBond(mol, 0, 3, 1);  // Cu-O coordination

    centerMolecule(mol);
}

// Build Verdigris (Cu2(OH)2CO3) - basic copper carbonate (patina)
void buildVerdigris(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Verdigris (Cu2(OH)2CO3)");

    // Two copper atoms with hydroxide and carbonate
    addAtom(mol, -1.5f, 0.0f, 0.0f, ATOM_CU);   // Cu1
    addAtom(mol, 1.5f, 0.0f, 0.0f, ATOM_CU);    // Cu2

    // Carbonate CO3
    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_C);     // C
    addAtom(mol, 0.0f, 1.2f, 0.0f, ATOM_O);     // O (bridging both Cu)
    addAtom(mol, -0.6f, -0.6f, 0.0f, ATOM_O);   // O
    addAtom(mol, 0.6f, -0.6f, 0.0f, ATOM_O);    // O

    // Hydroxide OH groups
    addAtom(mol, -2.5f, 1.0f, 0.0f, ATOM_O);    // OH
    addAtom(mol, -3.0f, 1.6f, 0.0f, ATOM_H);
    addAtom(mol, 2.5f, 1.0f, 0.0f, ATOM_O);     // OH
    addAtom(mol, 3.0f, 1.6f, 0.0f, ATOM_H);

    // Bonds
    addBond(mol, 2, 3, 2);  // C=O
    addBond(mol, 2, 4, 1);  // C-O
    addBond(mol, 2, 5, 1);  // C-O
    addBond(mol, 0, 4, 1);  // Cu-O
    addBond(mol, 1, 5, 1);  // Cu-O
    addBond(mol, 0, 6, 1);  // Cu-OH
    addBond(mol, 6, 7, 1);  // O-H
    addBond(mol, 1, 8, 1);  // Cu-OH
    addBond(mol, 8, 9, 1);  // O-H

    centerMolecule(mol);
}

// Build Ferrocene (Fe(C5H5)2) - sandwich compound
void buildFerrocene(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Ferrocene (Fe(C5H5)2)");

    // Central iron
    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_FE);    // Fe center

    // Top cyclopentadienyl ring (Cp)
    float r = 1.2f;
    float z_top = 1.7f;
    for (int i = 0; i < 5; i++) {
        float angle = i * 2.0f * 3.14159f / 5.0f;
        addAtom(mol, r*cosf(angle), r*sinf(angle), z_top, ATOM_C);
    }
    // H on top ring
    for (int i = 0; i < 5; i++) {
        float angle = i * 2.0f * 3.14159f / 5.0f;
        addAtom(mol, 1.9f*r*cosf(angle), 1.9f*r*sinf(angle), z_top+0.3f, ATOM_H);
    }

    // Bottom cyclopentadienyl ring (Cp) - staggered
    float z_bot = -1.7f;
    for (int i = 0; i < 5; i++) {
        float angle = (i + 0.5f) * 2.0f * 3.14159f / 5.0f;  // 36 degree offset
        addAtom(mol, r*cosf(angle), r*sinf(angle), z_bot, ATOM_C);
    }
    // H on bottom ring
    for (int i = 0; i < 5; i++) {
        float angle = (i + 0.5f) * 2.0f * 3.14159f / 5.0f;
        addAtom(mol, 1.9f*r*cosf(angle), 1.9f*r*sinf(angle), z_bot-0.3f, ATOM_H);
    }

    // Top ring C-C bonds (aromatic)
    for (int i = 0; i < 5; i++) {
        addBond(mol, 1 + i, 1 + ((i+1)%5), 1);  // C-C in ring
        addBond(mol, 1 + i, 6 + i, 1);          // C-H
    }

    // Bottom ring C-C bonds
    for (int i = 0; i < 5; i++) {
        addBond(mol, 11 + i, 11 + ((i+1)%5), 1);  // C-C in ring
        addBond(mol, 11 + i, 16 + i, 1);          // C-H
    }

    centerMolecule(mol);
}

// Build Cisplatin (Pt(NH3)2Cl2) - anticancer drug
void buildCisplatin(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Cisplatin (Pt(NH3)2Cl2)");

    // Square planar geometry
    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_PT);    // Pt center

    // Cis configuration: NH3 groups adjacent, Cl atoms adjacent
    addAtom(mol, 2.0f, 0.0f, 0.0f, ATOM_N);     // NH3
    addAtom(mol, 0.0f, 2.0f, 0.0f, ATOM_N);     // NH3
    addAtom(mol, -2.0f, 0.0f, 0.0f, ATOM_CL);   // Cl
    addAtom(mol, 0.0f, -2.0f, 0.0f, ATOM_CL);   // Cl

    // H atoms on NH3 groups
    addAtom(mol, 2.5f, 0.8f, 0.0f, ATOM_H);
    addAtom(mol, 2.5f, -0.4f, 0.7f, ATOM_H);
    addAtom(mol, 2.5f, -0.4f, -0.7f, ATOM_H);
    addAtom(mol, 0.8f, 2.5f, 0.0f, ATOM_H);
    addAtom(mol, -0.4f, 2.5f, 0.7f, ATOM_H);
    addAtom(mol, -0.4f, 2.5f, -0.7f, ATOM_H);

    // Pt-ligand bonds
    addBond(mol, 0, 1, 1);  // Pt-N
    addBond(mol, 0, 2, 1);  // Pt-N
    addBond(mol, 0, 3, 1);  // Pt-Cl
    addBond(mol, 0, 4, 1);  // Pt-Cl

    // N-H bonds
    addBond(mol, 1, 5, 1);
    addBond(mol, 1, 6, 1);
    addBond(mol, 1, 7, 1);
    addBond(mol, 2, 8, 1);
    addBond(mol, 2, 9, 1);
    addBond(mol, 2, 10, 1);

    centerMolecule(mol);
}

// Build Aluminum Oxide (Al2O3) - corundum/alumina
void buildAluminumOxide(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Aluminum Oxide (Al2O3)");

    // Simplified corundum structure
    addAtom(mol, -0.8f, 0.0f, 0.0f, ATOM_AL);   // Al
    addAtom(mol, 0.8f, 0.0f, 0.0f, ATOM_AL);    // Al
    addAtom(mol, 0.0f, 1.0f, 0.0f, ATOM_O);     // O bridging
    addAtom(mol, -1.3f, -0.8f, 0.0f, ATOM_O);   // O
    addAtom(mol, 1.3f, -0.8f, 0.0f, ATOM_O);    // O

    addBond(mol, 0, 2, 1);
    addBond(mol, 1, 2, 1);
    addBond(mol, 0, 3, 2);
    addBond(mol, 1, 4, 2);

    centerMolecule(mol);
}

// Build Titanium Dioxide (TiO2) - titania/white pigment
void buildTitaniumDioxide(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Titanium Dioxide (TiO2)");

    // Linear/bent structure
    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_TI);    // Ti center
    addAtom(mol, 1.6f, 0.0f, 0.0f, ATOM_O);     // O
    addAtom(mol, -1.6f, 0.0f, 0.0f, ATOM_O);    // O

    addBond(mol, 0, 1, 2);  // Ti=O
    addBond(mol, 0, 2, 2);  // Ti=O

    centerMolecule(mol);
}

// Build Carbon Dioxide (CO2)
void buildCarbonDioxide(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Carbon Dioxide (CO2)");

    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_C);
    addAtom(mol, -1.16f, 0.0f, 0.0f, ATOM_O);
    addAtom(mol, 1.16f, 0.0f, 0.0f, ATOM_O);

    addBond(mol, 0, 1, 2);
    addBond(mol, 0, 2, 2);

    centerMolecule(mol);
}

// Build Formaldehyde (CH2O)
void buildFormaldehyde(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Formaldehyde (CH2O)");

    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_C);
    addAtom(mol, 0.0f, 1.2f, 0.0f, ATOM_O);
    addAtom(mol, 0.93f, -0.54f, 0.0f, ATOM_H);
    addAtom(mol, -0.93f, -0.54f, 0.0f, ATOM_H);

    addBond(mol, 0, 1, 2);
    addBond(mol, 0, 2, 1);
    addBond(mol, 0, 3, 1);

    centerMolecule(mol);
}

// Build Acetone (C3H6O)
void buildAcetone(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Acetone (C3H6O)");

    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_C);      // Central C
    addAtom(mol, 0.0f, 1.2f, 0.0f, ATOM_O);      // =O
    addAtom(mol, -1.4f, -0.5f, 0.0f, ATOM_C);    // CH3
    addAtom(mol, 1.4f, -0.5f, 0.0f, ATOM_C);     // CH3

    // Hydrogens on CH3 groups
    addAtom(mol, -1.4f, -1.6f, 0.0f, ATOM_H);
    addAtom(mol, -2.1f, -0.1f, 0.8f, ATOM_H);
    addAtom(mol, -2.1f, -0.1f, -0.8f, ATOM_H);
    addAtom(mol, 1.4f, -1.6f, 0.0f, ATOM_H);
    addAtom(mol, 2.1f, -0.1f, 0.8f, ATOM_H);
    addAtom(mol, 2.1f, -0.1f, -0.8f, ATOM_H);

    addBond(mol, 0, 1, 2);
    addBond(mol, 0, 2, 1);
    addBond(mol, 0, 3, 1);
    addBond(mol, 2, 4, 1);
    addBond(mol, 2, 5, 1);
    addBond(mol, 2, 6, 1);
    addBond(mol, 3, 7, 1);
    addBond(mol, 3, 8, 1);
    addBond(mol, 3, 9, 1);

    centerMolecule(mol);
}

// Build Acetic Acid (CH3COOH)
void buildAceticAcid(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Acetic Acid (CH3COOH)");

    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_C);      // CH3
    addAtom(mol, 1.5f, 0.0f, 0.0f, ATOM_C);      // C=O
    addAtom(mol, 2.1f, 1.1f, 0.0f, ATOM_O);      // =O
    addAtom(mol, 2.1f, -1.1f, 0.0f, ATOM_O);     // OH
    addAtom(mol, 3.0f, -1.1f, 0.0f, ATOM_H);     // H on OH
    addAtom(mol, -0.5f, 1.0f, 0.0f, ATOM_H);
    addAtom(mol, -0.5f, -0.5f, 0.87f, ATOM_H);
    addAtom(mol, -0.5f, -0.5f, -0.87f, ATOM_H);

    addBond(mol, 0, 1, 1);
    addBond(mol, 1, 2, 2);
    addBond(mol, 1, 3, 1);
    addBond(mol, 3, 4, 1);
    addBond(mol, 0, 5, 1);
    addBond(mol, 0, 6, 1);
    addBond(mol, 0, 7, 1);

    centerMolecule(mol);
}

// Build Propane (C3H8)
void buildPropane(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Propane (C3H8)");

    addAtom(mol, -1.5f, 0.0f, 0.0f, ATOM_C);
    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_C);
    addAtom(mol, 1.5f, 0.0f, 0.0f, ATOM_C);

    // H on C1
    addAtom(mol, -2.0f, 1.0f, 0.0f, ATOM_H);
    addAtom(mol, -2.0f, -0.5f, 0.87f, ATOM_H);
    addAtom(mol, -2.0f, -0.5f, -0.87f, ATOM_H);
    // H on C2
    addAtom(mol, 0.0f, 0.6f, 0.9f, ATOM_H);
    addAtom(mol, 0.0f, 0.6f, -0.9f, ATOM_H);
    // H on C3
    addAtom(mol, 2.0f, 1.0f, 0.0f, ATOM_H);
    addAtom(mol, 2.0f, -0.5f, 0.87f, ATOM_H);
    addAtom(mol, 2.0f, -0.5f, -0.87f, ATOM_H);

    addBond(mol, 0, 1, 1);
    addBond(mol, 1, 2, 1);
    addBond(mol, 0, 3, 1);
    addBond(mol, 0, 4, 1);
    addBond(mol, 0, 5, 1);
    addBond(mol, 1, 6, 1);
    addBond(mol, 1, 7, 1);
    addBond(mol, 2, 8, 1);
    addBond(mol, 2, 9, 1);
    addBond(mol, 2, 10, 1);

    centerMolecule(mol);
}

// Build Butane (C4H10)
void buildButane(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Butane (C4H10)");

    addAtom(mol, -2.25f, 0.0f, 0.0f, ATOM_C);
    addAtom(mol, -0.75f, 0.0f, 0.0f, ATOM_C);
    addAtom(mol, 0.75f, 0.0f, 0.0f, ATOM_C);
    addAtom(mol, 2.25f, 0.0f, 0.0f, ATOM_C);

    // Hydrogens
    addAtom(mol, -2.75f, 1.0f, 0.0f, ATOM_H);
    addAtom(mol, -2.75f, -0.5f, 0.87f, ATOM_H);
    addAtom(mol, -2.75f, -0.5f, -0.87f, ATOM_H);
    addAtom(mol, -0.75f, 0.6f, 0.9f, ATOM_H);
    addAtom(mol, -0.75f, 0.6f, -0.9f, ATOM_H);
    addAtom(mol, 0.75f, 0.6f, 0.9f, ATOM_H);
    addAtom(mol, 0.75f, 0.6f, -0.9f, ATOM_H);
    addAtom(mol, 2.75f, 1.0f, 0.0f, ATOM_H);
    addAtom(mol, 2.75f, -0.5f, 0.87f, ATOM_H);
    addAtom(mol, 2.75f, -0.5f, -0.87f, ATOM_H);

    addBond(mol, 0, 1, 1);
    addBond(mol, 1, 2, 1);
    addBond(mol, 2, 3, 1);
    addBond(mol, 0, 4, 1);
    addBond(mol, 0, 5, 1);
    addBond(mol, 0, 6, 1);
    addBond(mol, 1, 7, 1);
    addBond(mol, 1, 8, 1);
    addBond(mol, 2, 9, 1);
    addBond(mol, 2, 10, 1);
    addBond(mol, 3, 11, 1);
    addBond(mol, 3, 12, 1);
    addBond(mol, 3, 13, 1);

    centerMolecule(mol);
}

// Build Cyclohexane (C6H12)
void buildCyclohexane(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Cyclohexane (C6H12)");

    float r = 1.4f;
    // Chair conformation
    for (int i = 0; i < 6; i++) {
        float angle = i * PI / 3.0f;
        float z = (i % 2 == 0) ? 0.3f : -0.3f;
        addAtom(mol, r * cosf(angle), r * sinf(angle), z, ATOM_C);
    }

    // Axial and equatorial hydrogens
    for (int i = 0; i < 6; i++) {
        float angle = i * PI / 3.0f;
        float z = (i % 2 == 0) ? 0.3f : -0.3f;
        // Axial H
        addAtom(mol, r * cosf(angle) * 0.6f, r * sinf(angle) * 0.6f, z + ((i % 2 == 0) ? 1.0f : -1.0f), ATOM_H);
        // Equatorial H
        addAtom(mol, (r + 0.9f) * cosf(angle), (r + 0.9f) * sinf(angle), z, ATOM_H);
    }

    // C-C bonds
    for (int i = 0; i < 6; i++) {
        addBond(mol, i, (i + 1) % 6, 1);
    }
    // C-H bonds
    for (int i = 0; i < 6; i++) {
        addBond(mol, i, 6 + i * 2, 1);
        addBond(mol, i, 6 + i * 2 + 1, 1);
    }

    centerMolecule(mol);
}

// Build Naphthalene (C10H8)
void buildNaphthalene(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Naphthalene (C10H8)");

    // Two fused benzene rings
    float dx = 1.2f;
    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_C);      // C1
    addAtom(mol, dx, 0.7f, 0.0f, ATOM_C);        // C2
    addAtom(mol, 2*dx, 0.0f, 0.0f, ATOM_C);      // C3
    addAtom(mol, 2*dx, -1.4f, 0.0f, ATOM_C);     // C4
    addAtom(mol, dx, -2.1f, 0.0f, ATOM_C);       // C5
    addAtom(mol, 0.0f, -1.4f, 0.0f, ATOM_C);     // C6
    // Second ring
    addAtom(mol, 3*dx, 0.7f, 0.0f, ATOM_C);      // C7
    addAtom(mol, 4*dx, 0.0f, 0.0f, ATOM_C);      // C8
    addAtom(mol, 4*dx, -1.4f, 0.0f, ATOM_C);     // C9
    addAtom(mol, 3*dx, -2.1f, 0.0f, ATOM_C);     // C10

    // Hydrogens
    addAtom(mol, -0.9f, 0.5f, 0.0f, ATOM_H);
    addAtom(mol, dx, 1.7f, 0.0f, ATOM_H);
    addAtom(mol, dx, -3.1f, 0.0f, ATOM_H);
    addAtom(mol, -0.9f, -1.9f, 0.0f, ATOM_H);
    addAtom(mol, 3*dx, 1.7f, 0.0f, ATOM_H);
    addAtom(mol, 4*dx + 0.9f, 0.5f, 0.0f, ATOM_H);
    addAtom(mol, 4*dx + 0.9f, -1.9f, 0.0f, ATOM_H);
    addAtom(mol, 3*dx, -3.1f, 0.0f, ATOM_H);

    // Bonds - first ring
    addBond(mol, 0, 1, 2);
    addBond(mol, 1, 2, 1);
    addBond(mol, 2, 3, 2);
    addBond(mol, 3, 4, 1);
    addBond(mol, 4, 5, 2);
    addBond(mol, 5, 0, 1);
    // Second ring
    addBond(mol, 2, 6, 1);
    addBond(mol, 6, 7, 2);
    addBond(mol, 7, 8, 1);
    addBond(mol, 8, 9, 2);
    addBond(mol, 9, 3, 1);
    // C-H bonds
    addBond(mol, 0, 10, 1);
    addBond(mol, 1, 11, 1);
    addBond(mol, 4, 12, 1);
    addBond(mol, 5, 13, 1);
    addBond(mol, 6, 14, 1);
    addBond(mol, 7, 15, 1);
    addBond(mol, 8, 16, 1);
    addBond(mol, 9, 17, 1);

    centerMolecule(mol);
}

// Build Urea (CH4N2O)
void buildUrea(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Urea (CH4N2O)");

    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_C);
    addAtom(mol, 0.0f, 1.25f, 0.0f, ATOM_O);
    addAtom(mol, -1.2f, -0.6f, 0.0f, ATOM_N);
    addAtom(mol, 1.2f, -0.6f, 0.0f, ATOM_N);
    addAtom(mol, -1.3f, -1.6f, 0.0f, ATOM_H);
    addAtom(mol, -2.0f, -0.1f, 0.0f, ATOM_H);
    addAtom(mol, 1.3f, -1.6f, 0.0f, ATOM_H);
    addAtom(mol, 2.0f, -0.1f, 0.0f, ATOM_H);

    addBond(mol, 0, 1, 2);
    addBond(mol, 0, 2, 1);
    addBond(mol, 0, 3, 1);
    addBond(mol, 2, 4, 1);
    addBond(mol, 2, 5, 1);
    addBond(mol, 3, 6, 1);
    addBond(mol, 3, 7, 1);

    centerMolecule(mol);
}

// Build Glycine (C2H5NO2) - simplest amino acid
void buildGlycine(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Glycine (C2H5NO2)");

    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_C);      // Alpha carbon
    addAtom(mol, -1.3f, 0.5f, 0.0f, ATOM_N);     // Amino
    addAtom(mol, 1.3f, 0.7f, 0.0f, ATOM_C);      // Carboxyl C
    addAtom(mol, 1.3f, 2.0f, 0.0f, ATOM_O);      // =O
    addAtom(mol, 2.4f, 0.0f, 0.0f, ATOM_O);      // OH
    addAtom(mol, 0.0f, -0.6f, 0.9f, ATOM_H);
    addAtom(mol, 0.0f, -0.6f, -0.9f, ATOM_H);
    addAtom(mol, -1.4f, 1.5f, 0.0f, ATOM_H);
    addAtom(mol, -2.1f, 0.0f, 0.0f, ATOM_H);
    addAtom(mol, 3.2f, 0.5f, 0.0f, ATOM_H);

    addBond(mol, 0, 1, 1);
    addBond(mol, 0, 2, 1);
    addBond(mol, 2, 3, 2);
    addBond(mol, 2, 4, 1);
    addBond(mol, 0, 5, 1);
    addBond(mol, 0, 6, 1);
    addBond(mol, 1, 7, 1);
    addBond(mol, 1, 8, 1);
    addBond(mol, 4, 9, 1);

    centerMolecule(mol);
}

// Build Alanine (C3H7NO2)
void buildAlanine(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Alanine (C3H7NO2)");

    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_C);      // Alpha carbon
    addAtom(mol, -1.3f, 0.5f, 0.0f, ATOM_N);     // Amino
    addAtom(mol, 1.3f, 0.7f, 0.0f, ATOM_C);      // Carboxyl C
    addAtom(mol, 1.3f, 2.0f, 0.0f, ATOM_O);      // =O
    addAtom(mol, 2.4f, 0.0f, 0.0f, ATOM_O);      // OH
    addAtom(mol, 0.0f, -1.5f, 0.0f, ATOM_C);     // CH3
    addAtom(mol, 0.0f, 0.5f, 0.9f, ATOM_H);
    addAtom(mol, -1.4f, 1.5f, 0.0f, ATOM_H);
    addAtom(mol, -2.1f, 0.0f, 0.0f, ATOM_H);
    addAtom(mol, 3.2f, 0.5f, 0.0f, ATOM_H);
    addAtom(mol, 0.0f, -2.1f, 0.9f, ATOM_H);
    addAtom(mol, 0.87f, -2.1f, -0.45f, ATOM_H);
    addAtom(mol, -0.87f, -2.1f, -0.45f, ATOM_H);

    addBond(mol, 0, 1, 1);
    addBond(mol, 0, 2, 1);
    addBond(mol, 2, 3, 2);
    addBond(mol, 2, 4, 1);
    addBond(mol, 0, 5, 1);
    addBond(mol, 0, 6, 1);
    addBond(mol, 1, 7, 1);
    addBond(mol, 1, 8, 1);
    addBond(mol, 4, 9, 1);
    addBond(mol, 5, 10, 1);
    addBond(mol, 5, 11, 1);
    addBond(mol, 5, 12, 1);

    centerMolecule(mol);
}

// Build Thymine (DNA base)
void buildThymine(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Thymine (DNA base)");

    // Pyrimidine ring
    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_N);      // N1
    addAtom(mol, 1.2f, 0.7f, 0.0f, ATOM_C);      // C2
    addAtom(mol, 1.2f, 2.1f, 0.0f, ATOM_N);      // N3
    addAtom(mol, 0.0f, 2.8f, 0.0f, ATOM_C);      // C4
    addAtom(mol, -1.2f, 2.1f, 0.0f, ATOM_C);     // C5
    addAtom(mol, -1.2f, 0.7f, 0.0f, ATOM_C);     // C6

    // Carbonyls
    addAtom(mol, 2.3f, 0.1f, 0.0f, ATOM_O);      // O on C2
    addAtom(mol, 0.0f, 4.1f, 0.0f, ATOM_O);      // O on C4

    // Methyl on C5
    addAtom(mol, -2.4f, 2.8f, 0.0f, ATOM_C);     // CH3

    // Hydrogens
    addAtom(mol, 0.0f, -1.0f, 0.0f, ATOM_H);     // H on N1
    addAtom(mol, 2.1f, 2.5f, 0.0f, ATOM_H);      // H on N3
    addAtom(mol, -2.1f, 0.1f, 0.0f, ATOM_H);     // H on C6
    addAtom(mol, -2.4f, 3.9f, 0.0f, ATOM_H);
    addAtom(mol, -3.1f, 2.4f, 0.8f, ATOM_H);
    addAtom(mol, -3.1f, 2.4f, -0.8f, ATOM_H);

    addBond(mol, 0, 1, 1);
    addBond(mol, 1, 2, 1);
    addBond(mol, 2, 3, 1);
    addBond(mol, 3, 4, 1);
    addBond(mol, 4, 5, 2);
    addBond(mol, 5, 0, 1);
    addBond(mol, 1, 6, 2);
    addBond(mol, 3, 7, 2);
    addBond(mol, 4, 8, 1);
    addBond(mol, 0, 9, 1);
    addBond(mol, 2, 10, 1);
    addBond(mol, 5, 11, 1);
    addBond(mol, 8, 12, 1);
    addBond(mol, 8, 13, 1);
    addBond(mol, 8, 14, 1);

    centerMolecule(mol);
}

// Build Cytosine (DNA base)
void buildCytosine(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Cytosine (DNA base)");

    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_N);      // N1
    addAtom(mol, 1.2f, 0.7f, 0.0f, ATOM_C);      // C2
    addAtom(mol, 1.2f, 2.1f, 0.0f, ATOM_N);      // N3
    addAtom(mol, 0.0f, 2.8f, 0.0f, ATOM_C);      // C4
    addAtom(mol, -1.2f, 2.1f, 0.0f, ATOM_C);     // C5
    addAtom(mol, -1.2f, 0.7f, 0.0f, ATOM_C);     // C6

    addAtom(mol, 2.3f, 0.1f, 0.0f, ATOM_O);      // O on C2
    addAtom(mol, 0.0f, 4.1f, 0.0f, ATOM_N);      // NH2 on C4

    addAtom(mol, 0.0f, -1.0f, 0.0f, ATOM_H);
    addAtom(mol, -2.1f, 2.5f, 0.0f, ATOM_H);
    addAtom(mol, -2.1f, 0.1f, 0.0f, ATOM_H);
    addAtom(mol, -0.8f, 4.6f, 0.0f, ATOM_H);
    addAtom(mol, 0.8f, 4.6f, 0.0f, ATOM_H);

    addBond(mol, 0, 1, 1);
    addBond(mol, 1, 2, 2);
    addBond(mol, 2, 3, 1);
    addBond(mol, 3, 4, 2);
    addBond(mol, 4, 5, 1);
    addBond(mol, 5, 0, 1);
    addBond(mol, 1, 6, 2);
    addBond(mol, 3, 7, 1);
    addBond(mol, 0, 8, 1);
    addBond(mol, 4, 9, 1);
    addBond(mol, 5, 10, 1);
    addBond(mol, 7, 11, 1);
    addBond(mol, 7, 12, 1);

    centerMolecule(mol);
}

// Build Guanine (DNA base)
void buildGuanine(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Guanine (DNA base)");

    // Purine ring system
    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_N);      // N1
    addAtom(mol, 1.3f, 0.3f, 0.0f, ATOM_C);      // C2
    addAtom(mol, 1.8f, 1.6f, 0.0f, ATOM_N);      // N3
    addAtom(mol, 0.9f, 2.6f, 0.0f, ATOM_C);      // C4
    addAtom(mol, -0.4f, 2.3f, 0.0f, ATOM_C);     // C5
    addAtom(mol, -0.8f, 1.0f, 0.0f, ATOM_C);     // C6

    addAtom(mol, 1.2f, 3.9f, 0.0f, ATOM_N);      // N7
    addAtom(mol, 0.0f, 4.5f, 0.0f, ATOM_C);      // C8
    addAtom(mol, -1.0f, 3.5f, 0.0f, ATOM_N);     // N9

    addAtom(mol, 2.0f, -0.6f, 0.0f, ATOM_N);     // NH2 on C2
    addAtom(mol, -2.0f, 0.6f, 0.0f, ATOM_O);     // O on C6

    addAtom(mol, -0.3f, -0.9f, 0.0f, ATOM_H);
    addAtom(mol, -0.2f, 5.5f, 0.0f, ATOM_H);
    addAtom(mol, -1.9f, 3.8f, 0.0f, ATOM_H);
    addAtom(mol, 1.6f, -1.5f, 0.0f, ATOM_H);
    addAtom(mol, 2.9f, -0.4f, 0.0f, ATOM_H);

    addBond(mol, 0, 1, 1);
    addBond(mol, 1, 2, 2);
    addBond(mol, 2, 3, 1);
    addBond(mol, 3, 4, 1);
    addBond(mol, 4, 5, 1);
    addBond(mol, 5, 0, 1);
    addBond(mol, 3, 6, 1);
    addBond(mol, 6, 7, 2);
    addBond(mol, 7, 8, 1);
    addBond(mol, 8, 4, 1);
    addBond(mol, 1, 9, 1);
    addBond(mol, 5, 10, 2);
    addBond(mol, 0, 11, 1);
    addBond(mol, 7, 12, 1);
    addBond(mol, 8, 13, 1);
    addBond(mol, 9, 14, 1);
    addBond(mol, 9, 15, 1);

    centerMolecule(mol);
}

// Build Dopamine (C8H11NO2)
void buildDopamine(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Dopamine (C8H11NO2)");

    float r = 1.4f;
    // Catechol ring (benzene with 2 OH)
    for (int i = 0; i < 6; i++) {
        float angle = i * PI / 3.0f;
        addAtom(mol, r * cosf(angle), r * sinf(angle), 0.0f, ATOM_C);
    }

    // OH groups on C3 and C4
    addAtom(mol, r * cosf(PI) - 0.8f, r * sinf(PI) + 0.6f, 0.0f, ATOM_O);
    addAtom(mol, r * cosf(4*PI/3) - 0.6f, r * sinf(4*PI/3) - 0.8f, 0.0f, ATOM_O);

    // Ethylamine chain from C1
    addAtom(mol, 2.5f, 0.5f, 0.0f, ATOM_C);      // CH2
    addAtom(mol, 3.8f, -0.2f, 0.0f, ATOM_C);     // CH2
    addAtom(mol, 5.0f, 0.5f, 0.0f, ATOM_N);      // NH2

    // Hydrogens
    float rH = 2.4f;
    addAtom(mol, rH * cosf(PI/3), rH * sinf(PI/3), 0.0f, ATOM_H);   // H on C2
    addAtom(mol, rH * cosf(5*PI/3), rH * sinf(5*PI/3), 0.0f, ATOM_H);// H on C6
    addAtom(mol, r * cosf(PI) - 1.6f, r * sinf(PI) + 0.3f, 0.0f, ATOM_H);  // H on OH
    addAtom(mol, r * cosf(4*PI/3) - 0.3f, r * sinf(4*PI/3) - 1.6f, 0.0f, ATOM_H);  // H on OH
    addAtom(mol, 2.5f, 1.1f, 0.9f, ATOM_H);
    addAtom(mol, 2.5f, 1.1f, -0.9f, ATOM_H);
    addAtom(mol, 3.8f, -0.8f, 0.9f, ATOM_H);
    addAtom(mol, 3.8f, -0.8f, -0.9f, ATOM_H);
    addAtom(mol, 5.0f, 1.1f, 0.8f, ATOM_H);
    addAtom(mol, 5.8f, 0.0f, 0.0f, ATOM_H);

    // Ring bonds
    for (int i = 0; i < 6; i++) {
        addBond(mol, i, (i + 1) % 6, (i % 2 == 0) ? 2 : 1);
    }
    // OH bonds
    addBond(mol, 3, 6, 1);
    addBond(mol, 4, 7, 1);
    // Chain bonds
    addBond(mol, 0, 8, 1);
    addBond(mol, 8, 9, 1);
    addBond(mol, 9, 10, 1);
    // H bonds
    addBond(mol, 1, 11, 1);
    addBond(mol, 5, 12, 1);
    addBond(mol, 6, 13, 1);
    addBond(mol, 7, 14, 1);
    addBond(mol, 8, 15, 1);
    addBond(mol, 8, 16, 1);
    addBond(mol, 9, 17, 1);
    addBond(mol, 9, 18, 1);
    addBond(mol, 10, 19, 1);
    addBond(mol, 10, 20, 1);

    centerMolecule(mol);
}

// Build Serotonin (C10H12N2O)
void buildSerotonin(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Serotonin (C10H12N2O)");

    // Indole ring system
    float r = 1.3f;
    // Benzene ring
    for (int i = 0; i < 6; i++) {
        float angle = i * PI / 3.0f + PI/6;
        addAtom(mol, r * cosf(angle), r * sinf(angle), 0.0f, ATOM_C);
    }
    // Pyrrole ring (fused)
    addAtom(mol, 2.3f, 0.7f, 0.0f, ATOM_C);      // C7
    addAtom(mol, 2.8f, -0.5f, 0.0f, ATOM_C);     // C8
    addAtom(mol, 1.9f, -1.3f, 0.0f, ATOM_N);     // N

    // OH on C5 position
    addAtom(mol, r * cosf(5*PI/3 + PI/6) + 0.9f, r * sinf(5*PI/3 + PI/6) - 0.5f, 0.0f, ATOM_O);

    // Ethylamine chain
    addAtom(mol, 4.1f, -0.8f, 0.0f, ATOM_C);     // CH2
    addAtom(mol, 5.0f, 0.3f, 0.0f, ATOM_C);      // CH2
    addAtom(mol, 6.3f, 0.0f, 0.0f, ATOM_N);      // NH2

    // Hydrogens (simplified)
    addAtom(mol, r * cosf(PI/3 + PI/6) * 1.7f, r * sinf(PI/3 + PI/6) * 1.7f, 0.0f, ATOM_H);
    addAtom(mol, r * cosf(2*PI/3 + PI/6) * 1.7f, r * sinf(2*PI/3 + PI/6) * 1.7f, 0.0f, ATOM_H);
    addAtom(mol, r * cosf(PI + PI/6) * 1.7f, r * sinf(PI + PI/6) * 1.7f, 0.0f, ATOM_H);
    addAtom(mol, 2.5f, 1.7f, 0.0f, ATOM_H);
    addAtom(mol, 1.9f, -2.3f, 0.0f, ATOM_H);
    addAtom(mol, r * cosf(5*PI/3 + PI/6) + 1.7f, r * sinf(5*PI/3 + PI/6) - 0.2f, 0.0f, ATOM_H);
    addAtom(mol, 4.3f, -1.4f, 0.9f, ATOM_H);
    addAtom(mol, 4.3f, -1.4f, -0.9f, ATOM_H);
    addAtom(mol, 4.8f, 0.9f, 0.9f, ATOM_H);
    addAtom(mol, 4.8f, 0.9f, -0.9f, ATOM_H);
    addAtom(mol, 6.5f, -0.5f, 0.8f, ATOM_H);
    addAtom(mol, 7.0f, 0.7f, 0.0f, ATOM_H);

    // Ring bonds
    for (int i = 0; i < 6; i++) {
        addBond(mol, i, (i + 1) % 6, (i % 2 == 0) ? 2 : 1);
    }
    addBond(mol, 0, 6, 1);
    addBond(mol, 6, 7, 2);
    addBond(mol, 7, 8, 1);
    addBond(mol, 8, 5, 1);
    addBond(mol, 4, 9, 1);
    addBond(mol, 7, 10, 1);
    addBond(mol, 10, 11, 1);
    addBond(mol, 11, 12, 1);
    // H bonds
    addBond(mol, 1, 13, 1);
    addBond(mol, 2, 14, 1);
    addBond(mol, 3, 15, 1);
    addBond(mol, 6, 16, 1);
    addBond(mol, 8, 17, 1);
    addBond(mol, 9, 18, 1);
    addBond(mol, 10, 19, 1);
    addBond(mol, 10, 20, 1);
    addBond(mol, 11, 21, 1);
    addBond(mol, 11, 22, 1);
    addBond(mol, 12, 23, 1);
    addBond(mol, 12, 24, 1);

    centerMolecule(mol);
}

// Build Nitric Oxide (NO)
void buildNitricOxide(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Nitric Oxide (NO)");

    addAtom(mol, -0.58f, 0.0f, 0.0f, ATOM_N);
    addAtom(mol, 0.58f, 0.0f, 0.0f, ATOM_O);

    addBond(mol, 0, 1, 2);

    centerMolecule(mol);
}

// Build Hydrogen Peroxide (H2O2)
void buildHydrogenPeroxide(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Hydrogen Peroxide (H2O2)");

    addAtom(mol, -0.7f, 0.0f, 0.0f, ATOM_O);
    addAtom(mol, 0.7f, 0.0f, 0.0f, ATOM_O);
    addAtom(mol, -1.2f, 0.8f, 0.3f, ATOM_H);
    addAtom(mol, 1.2f, -0.8f, -0.3f, ATOM_H);

    addBond(mol, 0, 1, 1);
    addBond(mol, 0, 2, 1);
    addBond(mol, 1, 3, 1);

    centerMolecule(mol);
}

// Build Sulfuric Acid (H2SO4)
void buildSulfuricAcid(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Sulfuric Acid (H2SO4)");

    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_S);
    addAtom(mol, 0.0f, 1.4f, 0.0f, ATOM_O);      // =O
    addAtom(mol, 0.0f, -1.4f, 0.0f, ATOM_O);     // =O
    addAtom(mol, 1.4f, 0.0f, 0.0f, ATOM_O);      // OH
    addAtom(mol, -1.4f, 0.0f, 0.0f, ATOM_O);     // OH
    addAtom(mol, 2.1f, 0.7f, 0.0f, ATOM_H);
    addAtom(mol, -2.1f, -0.7f, 0.0f, ATOM_H);

    addBond(mol, 0, 1, 2);
    addBond(mol, 0, 2, 2);
    addBond(mol, 0, 3, 1);
    addBond(mol, 0, 4, 1);
    addBond(mol, 3, 5, 1);
    addBond(mol, 4, 6, 1);

    centerMolecule(mol);
}

// Build Phosphoric Acid (H3PO4)
void buildPhosphoricAcid(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Phosphoric Acid (H3PO4)");

    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_P);
    addAtom(mol, 0.0f, 1.5f, 0.0f, ATOM_O);      // =O
    addAtom(mol, 1.3f, -0.75f, 0.0f, ATOM_O);    // OH
    addAtom(mol, -1.3f, -0.75f, 0.0f, ATOM_O);   // OH
    addAtom(mol, 0.0f, -0.5f, 1.4f, ATOM_O);     // OH
    addAtom(mol, 2.0f, -0.3f, 0.0f, ATOM_H);
    addAtom(mol, -2.0f, -0.3f, 0.0f, ATOM_H);
    addAtom(mol, 0.0f, 0.0f, 2.1f, ATOM_H);

    addBond(mol, 0, 1, 2);
    addBond(mol, 0, 2, 1);
    addBond(mol, 0, 3, 1);
    addBond(mol, 0, 4, 1);
    addBond(mol, 2, 5, 1);
    addBond(mol, 3, 6, 1);
    addBond(mol, 4, 7, 1);

    centerMolecule(mol);
}

// Build Toluene (C7H8) - methylbenzene
void buildToluene(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Toluene (C7H8)");

    float r = 1.4f;
    // Benzene ring
    for (int i = 0; i < 6; i++) {
        float angle = i * PI / 3.0f;
        addAtom(mol, r * cosf(angle), r * sinf(angle), 0.0f, ATOM_C);
    }

    // Methyl group on C1
    addAtom(mol, 2.6f, 0.0f, 0.0f, ATOM_C);

    // H on benzene (5 positions)
    float rH = 2.4f;
    for (int i = 1; i < 6; i++) {
        float angle = i * PI / 3.0f;
        addAtom(mol, rH * cosf(angle), rH * sinf(angle), 0.0f, ATOM_H);
    }

    // H on methyl
    addAtom(mol, 3.2f, 0.9f, 0.0f, ATOM_H);
    addAtom(mol, 3.2f, -0.45f, 0.8f, ATOM_H);
    addAtom(mol, 3.2f, -0.45f, -0.8f, ATOM_H);

    // Ring bonds
    for (int i = 0; i < 6; i++) {
        addBond(mol, i, (i + 1) % 6, (i % 2 == 0) ? 2 : 1);
    }
    addBond(mol, 0, 6, 1);
    // C-H on ring
    for (int i = 1; i < 6; i++) {
        addBond(mol, i, 6 + i, 1);
    }
    // Methyl H
    addBond(mol, 6, 12, 1);
    addBond(mol, 6, 13, 1);
    addBond(mol, 6, 14, 1);

    centerMolecule(mol);
}

// Build Phenol (C6H5OH)
void buildPhenol(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Phenol (C6H5OH)");

    float r = 1.4f;
    // Benzene ring
    for (int i = 0; i < 6; i++) {
        float angle = i * PI / 3.0f;
        addAtom(mol, r * cosf(angle), r * sinf(angle), 0.0f, ATOM_C);
    }

    // OH group on C1
    addAtom(mol, 2.4f, 0.0f, 0.0f, ATOM_O);
    addAtom(mol, 3.1f, 0.7f, 0.0f, ATOM_H);

    // H on benzene (5 positions)
    float rH = 2.4f;
    for (int i = 1; i < 6; i++) {
        float angle = i * PI / 3.0f;
        addAtom(mol, rH * cosf(angle), rH * sinf(angle), 0.0f, ATOM_H);
    }

    // Ring bonds
    for (int i = 0; i < 6; i++) {
        addBond(mol, i, (i + 1) % 6, (i % 2 == 0) ? 2 : 1);
    }
    addBond(mol, 0, 6, 1);
    addBond(mol, 6, 7, 1);
    for (int i = 1; i < 6; i++) {
        addBond(mol, i, 7 + i, 1);
    }

    centerMolecule(mol);
}

// Build Acetylene (C2H2)
void buildAcetylene(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Acetylene (C2H2)");

    addAtom(mol, -0.6f, 0.0f, 0.0f, ATOM_C);
    addAtom(mol, 0.6f, 0.0f, 0.0f, ATOM_C);
    addAtom(mol, -1.65f, 0.0f, 0.0f, ATOM_H);
    addAtom(mol, 1.65f, 0.0f, 0.0f, ATOM_H);

    addBond(mol, 0, 1, 3);
    addBond(mol, 0, 2, 1);
    addBond(mol, 1, 3, 1);

    centerMolecule(mol);
}

// Build Ethylene (C2H4)
void buildEthylene(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Ethylene (C2H4)");

    addAtom(mol, -0.67f, 0.0f, 0.0f, ATOM_C);
    addAtom(mol, 0.67f, 0.0f, 0.0f, ATOM_C);
    addAtom(mol, -1.23f, 0.92f, 0.0f, ATOM_H);
    addAtom(mol, -1.23f, -0.92f, 0.0f, ATOM_H);
    addAtom(mol, 1.23f, 0.92f, 0.0f, ATOM_H);
    addAtom(mol, 1.23f, -0.92f, 0.0f, ATOM_H);

    addBond(mol, 0, 1, 2);
    addBond(mol, 0, 2, 1);
    addBond(mol, 0, 3, 1);
    addBond(mol, 1, 4, 1);
    addBond(mol, 1, 5, 1);

    centerMolecule(mol);
}

// Build Hydrogen Cyanide (HCN)
void buildHydrogenCyanide(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Hydrogen Cyanide (HCN)");

    addAtom(mol, -1.15f, 0.0f, 0.0f, ATOM_H);
    addAtom(mol, -0.08f, 0.0f, 0.0f, ATOM_C);
    addAtom(mol, 1.08f, 0.0f, 0.0f, ATOM_N);

    addBond(mol, 0, 1, 1);
    addBond(mol, 1, 2, 3);

    centerMolecule(mol);
}

// Build Hydrogen Sulfide (H2S)
void buildHydrogenSulfide(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Hydrogen Sulfide (H2S)");

    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_S);
    addAtom(mol, 0.96f, 0.62f, 0.0f, ATOM_H);
    addAtom(mol, -0.96f, 0.62f, 0.0f, ATOM_H);

    addBond(mol, 0, 1, 1);
    addBond(mol, 0, 2, 1);

    centerMolecule(mol);
}

// Build Chloroform (CHCl3)
void buildChloroform(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Chloroform (CHCl3)");

    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_C);
    addAtom(mol, 0.0f, 1.07f, 0.0f, ATOM_H);
    addAtom(mol, 1.45f, -0.53f, 0.0f, ATOM_CL);
    addAtom(mol, -0.73f, -0.53f, 1.26f, ATOM_CL);
    addAtom(mol, -0.73f, -0.53f, -1.26f, ATOM_CL);

    addBond(mol, 0, 1, 1);
    addBond(mol, 0, 2, 1);
    addBond(mol, 0, 3, 1);
    addBond(mol, 0, 4, 1);

    centerMolecule(mol);
}

// Build Iodine (I2)
void buildIodine(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Iodine (I2)");

    addAtom(mol, -1.33f, 0.0f, 0.0f, ATOM_I);
    addAtom(mol, 1.33f, 0.0f, 0.0f, ATOM_I);

    addBond(mol, 0, 1, 1);

    centerMolecule(mol);
}

// Build Ozone (O3)
void buildOzone(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Ozone (O3)");

    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_O);
    addAtom(mol, 1.1f, 0.6f, 0.0f, ATOM_O);
    addAtom(mol, -1.1f, 0.6f, 0.0f, ATOM_O);

    addBond(mol, 0, 1, 2);
    addBond(mol, 0, 2, 1);

    centerMolecule(mol);
}

// Build Oxygen (O2)
void buildOxygen(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Oxygen (O2)");

    addAtom(mol, -0.6f, 0.0f, 0.0f, ATOM_O);
    addAtom(mol, 0.6f, 0.0f, 0.0f, ATOM_O);

    addBond(mol, 0, 1, 2);

    centerMolecule(mol);
}

// Build Nitrogen (N2)
void buildNitrogen(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Nitrogen (N2)");

    addAtom(mol, -0.55f, 0.0f, 0.0f, ATOM_N);
    addAtom(mol, 0.55f, 0.0f, 0.0f, ATOM_N);

    addBond(mol, 0, 1, 3);

    centerMolecule(mol);
}

// Build Hydrogen (H2)
void buildHydrogen(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Hydrogen (H2)");

    addAtom(mol, -0.37f, 0.0f, 0.0f, ATOM_H);
    addAtom(mol, 0.37f, 0.0f, 0.0f, ATOM_H);

    addBond(mol, 0, 1, 1);

    centerMolecule(mol);
}

// Build Carbon Monoxide (CO)
void buildCarbonMonoxide(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Carbon Monoxide (CO)");

    addAtom(mol, -0.56f, 0.0f, 0.0f, ATOM_C);
    addAtom(mol, 0.56f, 0.0f, 0.0f, ATOM_O);

    addBond(mol, 0, 1, 3);

    centerMolecule(mol);
}

// Build Nitrous Oxide (N2O)
void buildNitrousOxide(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Nitrous Oxide (N2O)");

    addAtom(mol, -1.13f, 0.0f, 0.0f, ATOM_N);
    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_N);
    addAtom(mol, 1.19f, 0.0f, 0.0f, ATOM_O);

    addBond(mol, 0, 1, 2);
    addBond(mol, 1, 2, 2);

    centerMolecule(mol);
}

// Build Sulfur Dioxide (SO2)
void buildSulfurDioxide(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Sulfur Dioxide (SO2)");

    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_S);
    addAtom(mol, 1.2f, 0.7f, 0.0f, ATOM_O);
    addAtom(mol, -1.2f, 0.7f, 0.0f, ATOM_O);

    addBond(mol, 0, 1, 2);
    addBond(mol, 0, 2, 2);

    centerMolecule(mol);
}

// Build Hydrogen Chloride (HCl)
void buildHydrogenChloride(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Hydrogen Chloride (HCl)");

    addAtom(mol, -0.64f, 0.0f, 0.0f, ATOM_H);
    addAtom(mol, 0.64f, 0.0f, 0.0f, ATOM_CL);

    addBond(mol, 0, 1, 1);

    centerMolecule(mol);
}

// Build Nitric Acid (HNO3)
void buildNitricAcid(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Nitric Acid (HNO3)");

    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_N);
    addAtom(mol, 1.2f, 0.5f, 0.0f, ATOM_O);      // =O
    addAtom(mol, -1.0f, 0.8f, 0.0f, ATOM_O);     // =O
    addAtom(mol, 0.3f, -1.2f, 0.0f, ATOM_O);     // OH
    addAtom(mol, 1.1f, -1.5f, 0.0f, ATOM_H);

    addBond(mol, 0, 1, 2);
    addBond(mol, 0, 2, 2);
    addBond(mol, 0, 3, 1);
    addBond(mol, 3, 4, 1);

    centerMolecule(mol);
}

// Build Methanol (CH3OH)
void buildMethanol(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Methanol (CH3OH)");

    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_C);
    addAtom(mol, 1.4f, 0.0f, 0.0f, ATOM_O);
    addAtom(mol, 2.1f, 0.8f, 0.0f, ATOM_H);
    addAtom(mol, -0.5f, 1.0f, 0.0f, ATOM_H);
    addAtom(mol, -0.5f, -0.5f, 0.87f, ATOM_H);
    addAtom(mol, -0.5f, -0.5f, -0.87f, ATOM_H);

    addBond(mol, 0, 1, 1);
    addBond(mol, 1, 2, 1);
    addBond(mol, 0, 3, 1);
    addBond(mol, 0, 4, 1);
    addBond(mol, 0, 5, 1);

    centerMolecule(mol);
}

// Build Ethane (C2H6)
void buildEthane(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Ethane (C2H6)");

    addAtom(mol, -0.77f, 0.0f, 0.0f, ATOM_C);
    addAtom(mol, 0.77f, 0.0f, 0.0f, ATOM_C);
    // H on C1
    addAtom(mol, -1.15f, 1.0f, 0.0f, ATOM_H);
    addAtom(mol, -1.15f, -0.5f, 0.87f, ATOM_H);
    addAtom(mol, -1.15f, -0.5f, -0.87f, ATOM_H);
    // H on C2
    addAtom(mol, 1.15f, 1.0f, 0.0f, ATOM_H);
    addAtom(mol, 1.15f, -0.5f, 0.87f, ATOM_H);
    addAtom(mol, 1.15f, -0.5f, -0.87f, ATOM_H);

    addBond(mol, 0, 1, 1);
    addBond(mol, 0, 2, 1);
    addBond(mol, 0, 3, 1);
    addBond(mol, 0, 4, 1);
    addBond(mol, 1, 5, 1);
    addBond(mol, 1, 6, 1);
    addBond(mol, 1, 7, 1);

    centerMolecule(mol);
}

// Build Propene (C3H6)
void buildPropene(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Propene (C3H6)");

    addAtom(mol, -1.3f, 0.0f, 0.0f, ATOM_C);     // CH3
    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_C);      // CH=
    addAtom(mol, 1.3f, 0.0f, 0.0f, ATOM_C);      // =CH2
    // H on CH3
    addAtom(mol, -1.7f, 1.0f, 0.0f, ATOM_H);
    addAtom(mol, -1.7f, -0.5f, 0.87f, ATOM_H);
    addAtom(mol, -1.7f, -0.5f, -0.87f, ATOM_H);
    // H on middle C
    addAtom(mol, 0.0f, 1.1f, 0.0f, ATOM_H);
    // H on =CH2
    addAtom(mol, 1.9f, 0.9f, 0.0f, ATOM_H);
    addAtom(mol, 1.9f, -0.9f, 0.0f, ATOM_H);

    addBond(mol, 0, 1, 1);
    addBond(mol, 1, 2, 2);
    addBond(mol, 0, 3, 1);
    addBond(mol, 0, 4, 1);
    addBond(mol, 0, 5, 1);
    addBond(mol, 1, 6, 1);
    addBond(mol, 2, 7, 1);
    addBond(mol, 2, 8, 1);

    centerMolecule(mol);
}

// Build Isopropanol / 2-Propanol (C3H8O)
void buildIsopropanol(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Isopropanol (C3H8O)");

    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_C);      // Central C with OH
    addAtom(mol, -1.5f, 0.0f, 0.0f, ATOM_C);     // CH3
    addAtom(mol, 1.5f, 0.0f, 0.0f, ATOM_C);      // CH3
    addAtom(mol, 0.0f, 1.4f, 0.0f, ATOM_O);      // OH
    addAtom(mol, 0.0f, 2.2f, 0.0f, ATOM_H);      // H on OH
    addAtom(mol, 0.0f, -1.1f, 0.0f, ATOM_H);     // H on central C
    // H on CH3 (left)
    addAtom(mol, -2.0f, 1.0f, 0.0f, ATOM_H);
    addAtom(mol, -2.0f, -0.5f, 0.87f, ATOM_H);
    addAtom(mol, -2.0f, -0.5f, -0.87f, ATOM_H);
    // H on CH3 (right)
    addAtom(mol, 2.0f, 1.0f, 0.0f, ATOM_H);
    addAtom(mol, 2.0f, -0.5f, 0.87f, ATOM_H);
    addAtom(mol, 2.0f, -0.5f, -0.87f, ATOM_H);

    addBond(mol, 0, 1, 1);
    addBond(mol, 0, 2, 1);
    addBond(mol, 0, 3, 1);
    addBond(mol, 3, 4, 1);
    addBond(mol, 0, 5, 1);
    addBond(mol, 1, 6, 1);
    addBond(mol, 1, 7, 1);
    addBond(mol, 1, 8, 1);
    addBond(mol, 2, 9, 1);
    addBond(mol, 2, 10, 1);
    addBond(mol, 2, 11, 1);

    centerMolecule(mol);
}

// Build Ethylene Glycol (C2H6O2)
void buildEthyleneGlycol(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Ethylene Glycol (C2H6O2)");

    addAtom(mol, -0.75f, 0.0f, 0.0f, ATOM_C);
    addAtom(mol, 0.75f, 0.0f, 0.0f, ATOM_C);
    addAtom(mol, -1.5f, 1.1f, 0.0f, ATOM_O);     // OH
    addAtom(mol, 1.5f, 1.1f, 0.0f, ATOM_O);      // OH
    addAtom(mol, -2.3f, 1.0f, 0.0f, ATOM_H);     // H on OH
    addAtom(mol, 2.3f, 1.0f, 0.0f, ATOM_H);      // H on OH
    // H on carbons
    addAtom(mol, -1.1f, -0.5f, 0.9f, ATOM_H);
    addAtom(mol, -1.1f, -0.5f, -0.9f, ATOM_H);
    addAtom(mol, 1.1f, -0.5f, 0.9f, ATOM_H);
    addAtom(mol, 1.1f, -0.5f, -0.9f, ATOM_H);

    addBond(mol, 0, 1, 1);
    addBond(mol, 0, 2, 1);
    addBond(mol, 1, 3, 1);
    addBond(mol, 2, 4, 1);
    addBond(mol, 3, 5, 1);
    addBond(mol, 0, 6, 1);
    addBond(mol, 0, 7, 1);
    addBond(mol, 1, 8, 1);
    addBond(mol, 1, 9, 1);

    centerMolecule(mol);
}

// Build Glycerol (C3H8O3)
void buildGlycerol(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Glycerol (C3H8O3)");

    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_C);      // Middle C
    addAtom(mol, -1.5f, 0.0f, 0.0f, ATOM_C);     // Left C
    addAtom(mol, 1.5f, 0.0f, 0.0f, ATOM_C);      // Right C
    addAtom(mol, 0.0f, 1.4f, 0.0f, ATOM_O);      // Middle OH
    addAtom(mol, -2.3f, 1.0f, 0.0f, ATOM_O);     // Left OH
    addAtom(mol, 2.3f, 1.0f, 0.0f, ATOM_O);      // Right OH
    addAtom(mol, 0.0f, 2.2f, 0.0f, ATOM_H);
    addAtom(mol, -3.1f, 0.8f, 0.0f, ATOM_H);
    addAtom(mol, 3.1f, 0.8f, 0.0f, ATOM_H);
    addAtom(mol, 0.0f, -1.1f, 0.0f, ATOM_H);
    addAtom(mol, -1.9f, -0.5f, 0.9f, ATOM_H);
    addAtom(mol, -1.9f, -0.5f, -0.9f, ATOM_H);
    addAtom(mol, 1.9f, -0.5f, 0.9f, ATOM_H);
    addAtom(mol, 1.9f, -0.5f, -0.9f, ATOM_H);

    addBond(mol, 0, 1, 1);
    addBond(mol, 0, 2, 1);
    addBond(mol, 0, 3, 1);
    addBond(mol, 1, 4, 1);
    addBond(mol, 2, 5, 1);
    addBond(mol, 3, 6, 1);
    addBond(mol, 4, 7, 1);
    addBond(mol, 5, 8, 1);
    addBond(mol, 0, 9, 1);
    addBond(mol, 1, 10, 1);
    addBond(mol, 1, 11, 1);
    addBond(mol, 2, 12, 1);
    addBond(mol, 2, 13, 1);

    centerMolecule(mol);
}

// Build Acetaldehyde / Ethanal (C2H4O)
void buildAcetaldehyde(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Acetaldehyde (C2H4O)");

    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_C);      // CHO
    addAtom(mol, 1.5f, 0.0f, 0.0f, ATOM_C);      // CH3
    addAtom(mol, -0.6f, 1.1f, 0.0f, ATOM_O);     // =O
    addAtom(mol, -0.6f, -1.0f, 0.0f, ATOM_H);    // H on CHO
    addAtom(mol, 2.0f, 1.0f, 0.0f, ATOM_H);
    addAtom(mol, 2.0f, -0.5f, 0.87f, ATOM_H);
    addAtom(mol, 2.0f, -0.5f, -0.87f, ATOM_H);

    addBond(mol, 0, 1, 1);
    addBond(mol, 0, 2, 2);
    addBond(mol, 0, 3, 1);
    addBond(mol, 1, 4, 1);
    addBond(mol, 1, 5, 1);
    addBond(mol, 1, 6, 1);

    centerMolecule(mol);
}

// Build Formic Acid (HCOOH)
void buildFormicAcid(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Formic Acid (HCOOH)");

    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_C);
    addAtom(mol, -0.6f, 1.1f, 0.0f, ATOM_O);     // =O
    addAtom(mol, 0.8f, -1.0f, 0.0f, ATOM_O);     // OH
    addAtom(mol, -0.9f, -0.6f, 0.0f, ATOM_H);    // H on C
    addAtom(mol, 1.6f, -0.6f, 0.0f, ATOM_H);     // H on OH

    addBond(mol, 0, 1, 2);
    addBond(mol, 0, 2, 1);
    addBond(mol, 0, 3, 1);
    addBond(mol, 2, 4, 1);

    centerMolecule(mol);
}

// Build Lactic Acid (C3H6O3)
void buildLacticAcid(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Lactic Acid (C3H6O3)");

    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_C);      // Central C with OH
    addAtom(mol, -1.5f, 0.0f, 0.0f, ATOM_C);     // CH3
    addAtom(mol, 1.5f, 0.0f, 0.0f, ATOM_C);      // COOH
    addAtom(mol, 0.0f, 1.4f, 0.0f, ATOM_O);      // OH
    addAtom(mol, 2.1f, 1.1f, 0.0f, ATOM_O);      // =O
    addAtom(mol, 2.1f, -1.1f, 0.0f, ATOM_O);     // OH
    addAtom(mol, 0.0f, 2.2f, 0.0f, ATOM_H);
    addAtom(mol, 2.9f, -1.0f, 0.0f, ATOM_H);
    addAtom(mol, 0.0f, -1.1f, 0.0f, ATOM_H);
    addAtom(mol, -2.0f, 1.0f, 0.0f, ATOM_H);
    addAtom(mol, -2.0f, -0.5f, 0.87f, ATOM_H);
    addAtom(mol, -2.0f, -0.5f, -0.87f, ATOM_H);

    addBond(mol, 0, 1, 1);
    addBond(mol, 0, 2, 1);
    addBond(mol, 0, 3, 1);
    addBond(mol, 2, 4, 2);
    addBond(mol, 2, 5, 1);
    addBond(mol, 3, 6, 1);
    addBond(mol, 5, 7, 1);
    addBond(mol, 0, 8, 1);
    addBond(mol, 1, 9, 1);
    addBond(mol, 1, 10, 1);
    addBond(mol, 1, 11, 1);

    centerMolecule(mol);
}

// Build Ethyl Acetate (C4H8O2)
void buildEthylAcetate(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Ethyl Acetate (C4H8O2)");

    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_C);      // C=O
    addAtom(mol, -1.5f, 0.0f, 0.0f, ATOM_C);     // CH3 (acetyl)
    addAtom(mol, 0.6f, 1.1f, 0.0f, ATOM_O);      // =O
    addAtom(mol, 0.8f, -1.2f, 0.0f, ATOM_O);     // O-
    addAtom(mol, 2.2f, -1.2f, 0.0f, ATOM_C);     // CH2
    addAtom(mol, 3.0f, -2.4f, 0.0f, ATOM_C);     // CH3 (ethyl)
    // H on acetyl CH3
    addAtom(mol, -2.0f, 1.0f, 0.0f, ATOM_H);
    addAtom(mol, -2.0f, -0.5f, 0.87f, ATOM_H);
    addAtom(mol, -2.0f, -0.5f, -0.87f, ATOM_H);
    // H on CH2
    addAtom(mol, 2.6f, -0.7f, 0.9f, ATOM_H);
    addAtom(mol, 2.6f, -0.7f, -0.9f, ATOM_H);
    // H on ethyl CH3
    addAtom(mol, 4.1f, -2.3f, 0.0f, ATOM_H);
    addAtom(mol, 2.7f, -2.9f, 0.9f, ATOM_H);
    addAtom(mol, 2.7f, -2.9f, -0.9f, ATOM_H);

    addBond(mol, 0, 1, 1);
    addBond(mol, 0, 2, 2);
    addBond(mol, 0, 3, 1);
    addBond(mol, 3, 4, 1);
    addBond(mol, 4, 5, 1);
    addBond(mol, 1, 6, 1);
    addBond(mol, 1, 7, 1);
    addBond(mol, 1, 8, 1);
    addBond(mol, 4, 9, 1);
    addBond(mol, 4, 10, 1);
    addBond(mol, 5, 11, 1);
    addBond(mol, 5, 12, 1);
    addBond(mol, 5, 13, 1);

    centerMolecule(mol);
}

// Build Acetonitrile (C2H3N)
void buildAcetonitrile(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Acetonitrile (C2H3N)");

    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_C);      // CH3
    addAtom(mol, 1.46f, 0.0f, 0.0f, ATOM_C);     // C
    addAtom(mol, 2.62f, 0.0f, 0.0f, ATOM_N);     // N
    addAtom(mol, -0.5f, 1.0f, 0.0f, ATOM_H);
    addAtom(mol, -0.5f, -0.5f, 0.87f, ATOM_H);
    addAtom(mol, -0.5f, -0.5f, -0.87f, ATOM_H);

    addBond(mol, 0, 1, 1);
    addBond(mol, 1, 2, 3);
    addBond(mol, 0, 3, 1);
    addBond(mol, 0, 4, 1);
    addBond(mol, 0, 5, 1);

    centerMolecule(mol);
}

// Build DMSO / Dimethyl Sulfoxide (C2H6OS)
void buildDMSO(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "DMSO (C2H6OS)");

    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_S);
    addAtom(mol, 0.0f, 1.5f, 0.0f, ATOM_O);
    addAtom(mol, -1.5f, -0.5f, 0.0f, ATOM_C);    // CH3
    addAtom(mol, 1.5f, -0.5f, 0.0f, ATOM_C);     // CH3
    // H on CH3 (left)
    addAtom(mol, -2.0f, 0.5f, 0.0f, ATOM_H);
    addAtom(mol, -2.0f, -1.0f, 0.87f, ATOM_H);
    addAtom(mol, -2.0f, -1.0f, -0.87f, ATOM_H);
    // H on CH3 (right)
    addAtom(mol, 2.0f, 0.5f, 0.0f, ATOM_H);
    addAtom(mol, 2.0f, -1.0f, 0.87f, ATOM_H);
    addAtom(mol, 2.0f, -1.0f, -0.87f, ATOM_H);

    addBond(mol, 0, 1, 2);
    addBond(mol, 0, 2, 1);
    addBond(mol, 0, 3, 1);
    addBond(mol, 2, 4, 1);
    addBond(mol, 2, 5, 1);
    addBond(mol, 2, 6, 1);
    addBond(mol, 3, 7, 1);
    addBond(mol, 3, 8, 1);
    addBond(mol, 3, 9, 1);

    centerMolecule(mol);
}

// Build Dichloromethane (CH2Cl2)
void buildDichloromethane(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Dichloromethane (CH2Cl2)");

    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_C);
    addAtom(mol, 1.3f, 0.7f, 0.0f, ATOM_CL);
    addAtom(mol, -1.3f, 0.7f, 0.0f, ATOM_CL);
    addAtom(mol, 0.0f, -0.6f, 0.9f, ATOM_H);
    addAtom(mol, 0.0f, -0.6f, -0.9f, ATOM_H);

    addBond(mol, 0, 1, 1);
    addBond(mol, 0, 2, 1);
    addBond(mol, 0, 3, 1);
    addBond(mol, 0, 4, 1);

    centerMolecule(mol);
}

// Build Chlorobenzene (C6H5Cl)
void buildChlorobenzene(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Chlorobenzene (C6H5Cl)");

    float r = 1.4f;
    // Benzene ring
    for (int i = 0; i < 6; i++) {
        float angle = i * PI / 3.0f;
        addAtom(mol, r * cosf(angle), r * sinf(angle), 0.0f, ATOM_C);
    }

    // Chlorine on C1
    addAtom(mol, 2.6f, 0.0f, 0.0f, ATOM_CL);

    // H on benzene (5 positions)
    float rH = 2.4f;
    for (int i = 1; i < 6; i++) {
        float angle = i * PI / 3.0f;
        addAtom(mol, rH * cosf(angle), rH * sinf(angle), 0.0f, ATOM_H);
    }

    // Ring bonds
    for (int i = 0; i < 6; i++) {
        addBond(mol, i, (i + 1) % 6, (i % 2 == 0) ? 2 : 1);
    }
    addBond(mol, 0, 6, 1);
    for (int i = 1; i < 6; i++) {
        addBond(mol, i, 6 + i, 1);
    }

    centerMolecule(mol);
}

// Build Nitrobenzene (C6H5NO2)
void buildNitrobenzene(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Nitrobenzene (C6H5NO2)");

    float r = 1.4f;
    // Benzene ring
    for (int i = 0; i < 6; i++) {
        float angle = i * PI / 3.0f;
        addAtom(mol, r * cosf(angle), r * sinf(angle), 0.0f, ATOM_C);
    }

    // NO2 group on C1
    addAtom(mol, 2.6f, 0.0f, 0.0f, ATOM_N);
    addAtom(mol, 3.3f, 1.0f, 0.0f, ATOM_O);
    addAtom(mol, 3.3f, -1.0f, 0.0f, ATOM_O);

    // H on benzene (5 positions)
    float rH = 2.4f;
    for (int i = 1; i < 6; i++) {
        float angle = i * PI / 3.0f;
        addAtom(mol, rH * cosf(angle), rH * sinf(angle), 0.0f, ATOM_H);
    }

    // Ring bonds
    for (int i = 0; i < 6; i++) {
        addBond(mol, i, (i + 1) % 6, (i % 2 == 0) ? 2 : 1);
    }
    addBond(mol, 0, 6, 1);
    addBond(mol, 6, 7, 2);
    addBond(mol, 6, 8, 2);
    for (int i = 1; i < 6; i++) {
        addBond(mol, i, 8 + i, 1);
    }

    centerMolecule(mol);
}

// Build Aniline (C6H5NH2)
void buildAniline(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Aniline (C6H5NH2)");

    float r = 1.4f;
    // Benzene ring
    for (int i = 0; i < 6; i++) {
        float angle = i * PI / 3.0f;
        addAtom(mol, r * cosf(angle), r * sinf(angle), 0.0f, ATOM_C);
    }

    // NH2 group on C1
    addAtom(mol, 2.6f, 0.0f, 0.0f, ATOM_N);
    addAtom(mol, 3.1f, 0.9f, 0.0f, ATOM_H);
    addAtom(mol, 3.1f, -0.9f, 0.0f, ATOM_H);

    // H on benzene (5 positions)
    float rH = 2.4f;
    for (int i = 1; i < 6; i++) {
        float angle = i * PI / 3.0f;
        addAtom(mol, rH * cosf(angle), rH * sinf(angle), 0.0f, ATOM_H);
    }

    // Ring bonds
    for (int i = 0; i < 6; i++) {
        addBond(mol, i, (i + 1) % 6, (i % 2 == 0) ? 2 : 1);
    }
    addBond(mol, 0, 6, 1);
    addBond(mol, 6, 7, 1);
    addBond(mol, 6, 8, 1);
    for (int i = 1; i < 6; i++) {
        addBond(mol, i, 8 + i, 1);
    }

    centerMolecule(mol);
}

// Build Styrene (C8H8)
void buildStyrene(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Styrene (C8H8)");

    float r = 1.4f;
    // Benzene ring
    for (int i = 0; i < 6; i++) {
        float angle = i * PI / 3.0f;
        addAtom(mol, r * cosf(angle), r * sinf(angle), 0.0f, ATOM_C);
    }

    // Vinyl group (-CH=CH2)
    addAtom(mol, 2.6f, 0.0f, 0.0f, ATOM_C);      // =CH-
    addAtom(mol, 3.9f, 0.0f, 0.0f, ATOM_C);      // =CH2

    // H on vinyl
    addAtom(mol, 2.6f, 1.1f, 0.0f, ATOM_H);
    addAtom(mol, 4.5f, 0.9f, 0.0f, ATOM_H);
    addAtom(mol, 4.5f, -0.9f, 0.0f, ATOM_H);

    // H on benzene (5 positions)
    float rH = 2.4f;
    for (int i = 1; i < 6; i++) {
        float angle = i * PI / 3.0f;
        addAtom(mol, rH * cosf(angle), rH * sinf(angle), 0.0f, ATOM_H);
    }

    // Ring bonds
    for (int i = 0; i < 6; i++) {
        addBond(mol, i, (i + 1) % 6, (i % 2 == 0) ? 2 : 1);
    }
    addBond(mol, 0, 6, 1);
    addBond(mol, 6, 7, 2);
    addBond(mol, 6, 8, 1);
    addBond(mol, 7, 9, 1);
    addBond(mol, 7, 10, 1);
    for (int i = 1; i < 6; i++) {
        addBond(mol, i, 10 + i, 1);
    }

    centerMolecule(mol);
}

// Build Benzoic Acid (C7H6O2)
void buildBenzoicAcid(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Benzoic Acid (C7H6O2)");

    float r = 1.4f;
    // Benzene ring
    for (int i = 0; i < 6; i++) {
        float angle = i * PI / 3.0f;
        addAtom(mol, r * cosf(angle), r * sinf(angle), 0.0f, ATOM_C);
    }

    // COOH group
    addAtom(mol, 2.6f, 0.0f, 0.0f, ATOM_C);
    addAtom(mol, 3.2f, 1.1f, 0.0f, ATOM_O);      // =O
    addAtom(mol, 3.2f, -1.1f, 0.0f, ATOM_O);     // OH
    addAtom(mol, 4.0f, -1.0f, 0.0f, ATOM_H);     // H on OH

    // H on benzene (5 positions)
    float rH = 2.4f;
    for (int i = 1; i < 6; i++) {
        float angle = i * PI / 3.0f;
        addAtom(mol, rH * cosf(angle), rH * sinf(angle), 0.0f, ATOM_H);
    }

    // Ring bonds
    for (int i = 0; i < 6; i++) {
        addBond(mol, i, (i + 1) % 6, (i % 2 == 0) ? 2 : 1);
    }
    addBond(mol, 0, 6, 1);
    addBond(mol, 6, 7, 2);
    addBond(mol, 6, 8, 1);
    addBond(mol, 8, 9, 1);
    for (int i = 1; i < 6; i++) {
        addBond(mol, i, 9 + i, 1);
    }

    centerMolecule(mol);
}

// Build Valine (C5H11NO2)
void buildValine(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Valine (C5H11NO2)");

    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_C);      // Alpha C
    addAtom(mol, -1.3f, 0.5f, 0.0f, ATOM_N);     // NH2
    addAtom(mol, 1.3f, 0.7f, 0.0f, ATOM_C);      // COOH C
    addAtom(mol, 1.3f, 2.0f, 0.0f, ATOM_O);      // =O
    addAtom(mol, 2.4f, 0.0f, 0.0f, ATOM_O);      // OH
    addAtom(mol, 0.0f, -1.5f, 0.0f, ATOM_C);     // Side chain CH
    addAtom(mol, -1.3f, -2.2f, 0.0f, ATOM_C);    // CH3
    addAtom(mol, 1.3f, -2.2f, 0.0f, ATOM_C);     // CH3
    // H atoms
    addAtom(mol, 0.0f, 0.5f, 0.9f, ATOM_H);      // H on alpha C
    addAtom(mol, -1.4f, 1.5f, 0.0f, ATOM_H);     // H on NH2
    addAtom(mol, -2.1f, 0.0f, 0.0f, ATOM_H);     // H on NH2
    addAtom(mol, 3.2f, 0.5f, 0.0f, ATOM_H);      // H on OH
    addAtom(mol, 0.0f, -2.1f, 0.9f, ATOM_H);     // H on side chain CH
    // H on CH3 (left)
    addAtom(mol, -1.3f, -3.3f, 0.0f, ATOM_H);
    addAtom(mol, -2.0f, -1.8f, 0.8f, ATOM_H);
    addAtom(mol, -2.0f, -1.8f, -0.8f, ATOM_H);
    // H on CH3 (right)
    addAtom(mol, 1.3f, -3.3f, 0.0f, ATOM_H);
    addAtom(mol, 2.0f, -1.8f, 0.8f, ATOM_H);
    addAtom(mol, 2.0f, -1.8f, -0.8f, ATOM_H);

    addBond(mol, 0, 1, 1);
    addBond(mol, 0, 2, 1);
    addBond(mol, 2, 3, 2);
    addBond(mol, 2, 4, 1);
    addBond(mol, 0, 5, 1);
    addBond(mol, 5, 6, 1);
    addBond(mol, 5, 7, 1);
    addBond(mol, 0, 8, 1);
    addBond(mol, 1, 9, 1);
    addBond(mol, 1, 10, 1);
    addBond(mol, 4, 11, 1);
    addBond(mol, 5, 12, 1);
    addBond(mol, 6, 13, 1);
    addBond(mol, 6, 14, 1);
    addBond(mol, 6, 15, 1);
    addBond(mol, 7, 16, 1);
    addBond(mol, 7, 17, 1);
    addBond(mol, 7, 18, 1);

    centerMolecule(mol);
}

// Build Leucine (C6H13NO2)
void buildLeucine(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Leucine (C6H13NO2)");

    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_C);      // Alpha C
    addAtom(mol, -1.3f, 0.5f, 0.0f, ATOM_N);     // NH2
    addAtom(mol, 1.3f, 0.7f, 0.0f, ATOM_C);      // COOH C
    addAtom(mol, 1.3f, 2.0f, 0.0f, ATOM_O);      // =O
    addAtom(mol, 2.4f, 0.0f, 0.0f, ATOM_O);      // OH
    addAtom(mol, 0.0f, -1.5f, 0.0f, ATOM_C);     // CH2
    addAtom(mol, 0.0f, -3.0f, 0.0f, ATOM_C);     // CH
    addAtom(mol, -1.3f, -3.7f, 0.0f, ATOM_C);    // CH3
    addAtom(mol, 1.3f, -3.7f, 0.0f, ATOM_C);     // CH3
    // H atoms
    addAtom(mol, 0.0f, 0.5f, 0.9f, ATOM_H);      // H on alpha C
    addAtom(mol, -1.4f, 1.5f, 0.0f, ATOM_H);     // H on NH2
    addAtom(mol, -2.1f, 0.0f, 0.0f, ATOM_H);     // H on NH2
    addAtom(mol, 3.2f, 0.5f, 0.0f, ATOM_H);      // H on OH
    addAtom(mol, -0.9f, -1.8f, 0.0f, ATOM_H);    // H on CH2
    addAtom(mol, 0.9f, -1.8f, 0.0f, ATOM_H);     // H on CH2
    addAtom(mol, 0.0f, -3.6f, 0.9f, ATOM_H);     // H on CH
    // H on CH3 (left)
    addAtom(mol, -1.3f, -4.8f, 0.0f, ATOM_H);
    addAtom(mol, -2.0f, -3.3f, 0.8f, ATOM_H);
    addAtom(mol, -2.0f, -3.3f, -0.8f, ATOM_H);
    // H on CH3 (right)
    addAtom(mol, 1.3f, -4.8f, 0.0f, ATOM_H);
    addAtom(mol, 2.0f, -3.3f, 0.8f, ATOM_H);
    addAtom(mol, 2.0f, -3.3f, -0.8f, ATOM_H);

    addBond(mol, 0, 1, 1);
    addBond(mol, 0, 2, 1);
    addBond(mol, 2, 3, 2);
    addBond(mol, 2, 4, 1);
    addBond(mol, 0, 5, 1);
    addBond(mol, 5, 6, 1);
    addBond(mol, 6, 7, 1);
    addBond(mol, 6, 8, 1);
    addBond(mol, 0, 9, 1);
    addBond(mol, 1, 10, 1);
    addBond(mol, 1, 11, 1);
    addBond(mol, 4, 12, 1);
    addBond(mol, 5, 13, 1);
    addBond(mol, 5, 14, 1);
    addBond(mol, 6, 15, 1);
    addBond(mol, 7, 16, 1);
    addBond(mol, 7, 17, 1);
    addBond(mol, 7, 18, 1);
    addBond(mol, 8, 19, 1);
    addBond(mol, 8, 20, 1);
    addBond(mol, 8, 21, 1);

    centerMolecule(mol);
}

// Build tert-Butanol (C4H10O)
void buildTertButanol(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "tert-Butanol (C4H10O)");

    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_C);      // Central C
    addAtom(mol, 0.0f, 1.5f, 0.0f, ATOM_O);      // OH
    addAtom(mol, 0.0f, 2.3f, 0.0f, ATOM_H);
    addAtom(mol, 1.4f, -0.5f, 0.0f, ATOM_C);     // CH3
    addAtom(mol, -0.7f, -0.5f, 1.2f, ATOM_C);    // CH3
    addAtom(mol, -0.7f, -0.5f, -1.2f, ATOM_C);   // CH3
    // H on CH3 groups
    addAtom(mol, 2.0f, 0.4f, 0.0f, ATOM_H);
    addAtom(mol, 1.9f, -1.0f, 0.87f, ATOM_H);
    addAtom(mol, 1.9f, -1.0f, -0.87f, ATOM_H);
    addAtom(mol, -0.1f, 0.0f, 2.0f, ATOM_H);
    addAtom(mol, -0.2f, -1.5f, 1.2f, ATOM_H);
    addAtom(mol, -1.7f, -0.3f, 1.4f, ATOM_H);
    addAtom(mol, -0.1f, 0.0f, -2.0f, ATOM_H);
    addAtom(mol, -0.2f, -1.5f, -1.2f, ATOM_H);
    addAtom(mol, -1.7f, -0.3f, -1.4f, ATOM_H);

    addBond(mol, 0, 1, 1);
    addBond(mol, 1, 2, 1);
    addBond(mol, 0, 3, 1);
    addBond(mol, 0, 4, 1);
    addBond(mol, 0, 5, 1);
    addBond(mol, 3, 6, 1);
    addBond(mol, 3, 7, 1);
    addBond(mol, 3, 8, 1);
    addBond(mol, 4, 9, 1);
    addBond(mol, 4, 10, 1);
    addBond(mol, 4, 11, 1);
    addBond(mol, 5, 12, 1);
    addBond(mol, 5, 13, 1);
    addBond(mol, 5, 14, 1);

    centerMolecule(mol);
}

// Build 1-Butanol (C4H10O)
void buildButanol(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "1-Butanol (C4H10O)");

    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_C);
    addAtom(mol, 1.5f, 0.0f, 0.0f, ATOM_C);
    addAtom(mol, 3.0f, 0.0f, 0.0f, ATOM_C);
    addAtom(mol, 4.5f, 0.0f, 0.0f, ATOM_C);
    addAtom(mol, 5.9f, 0.0f, 0.0f, ATOM_O);
    addAtom(mol, 6.5f, 0.8f, 0.0f, ATOM_H);
    // H atoms
    addAtom(mol, -0.5f, 1.0f, 0.0f, ATOM_H);
    addAtom(mol, -0.5f, -0.5f, 0.87f, ATOM_H);
    addAtom(mol, -0.5f, -0.5f, -0.87f, ATOM_H);
    addAtom(mol, 1.5f, 0.6f, 0.9f, ATOM_H);
    addAtom(mol, 1.5f, 0.6f, -0.9f, ATOM_H);
    addAtom(mol, 3.0f, 0.6f, 0.9f, ATOM_H);
    addAtom(mol, 3.0f, 0.6f, -0.9f, ATOM_H);
    addAtom(mol, 4.5f, 0.6f, 0.9f, ATOM_H);
    addAtom(mol, 4.5f, 0.6f, -0.9f, ATOM_H);

    addBond(mol, 0, 1, 1);
    addBond(mol, 1, 2, 1);
    addBond(mol, 2, 3, 1);
    addBond(mol, 3, 4, 1);
    addBond(mol, 4, 5, 1);
    addBond(mol, 0, 6, 1);
    addBond(mol, 0, 7, 1);
    addBond(mol, 0, 8, 1);
    addBond(mol, 1, 9, 1);
    addBond(mol, 1, 10, 1);
    addBond(mol, 2, 11, 1);
    addBond(mol, 2, 12, 1);
    addBond(mol, 3, 13, 1);
    addBond(mol, 3, 14, 1);

    centerMolecule(mol);
}

// Build Diethyl Ether (C4H10O)
void buildDiethylEther(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Diethyl Ether (C4H10O)");

    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_O);
    addAtom(mol, -1.4f, 0.0f, 0.0f, ATOM_C);     // CH2
    addAtom(mol, 1.4f, 0.0f, 0.0f, ATOM_C);      // CH2
    addAtom(mol, -2.9f, 0.0f, 0.0f, ATOM_C);     // CH3
    addAtom(mol, 2.9f, 0.0f, 0.0f, ATOM_C);      // CH3
    // H on CH2 groups
    addAtom(mol, -1.4f, 0.6f, 0.9f, ATOM_H);
    addAtom(mol, -1.4f, 0.6f, -0.9f, ATOM_H);
    addAtom(mol, 1.4f, 0.6f, 0.9f, ATOM_H);
    addAtom(mol, 1.4f, 0.6f, -0.9f, ATOM_H);
    // H on CH3 groups
    addAtom(mol, -3.4f, 1.0f, 0.0f, ATOM_H);
    addAtom(mol, -3.4f, -0.5f, 0.87f, ATOM_H);
    addAtom(mol, -3.4f, -0.5f, -0.87f, ATOM_H);
    addAtom(mol, 3.4f, 1.0f, 0.0f, ATOM_H);
    addAtom(mol, 3.4f, -0.5f, 0.87f, ATOM_H);
    addAtom(mol, 3.4f, -0.5f, -0.87f, ATOM_H);

    addBond(mol, 0, 1, 1);
    addBond(mol, 0, 2, 1);
    addBond(mol, 1, 3, 1);
    addBond(mol, 2, 4, 1);
    addBond(mol, 1, 5, 1);
    addBond(mol, 1, 6, 1);
    addBond(mol, 2, 7, 1);
    addBond(mol, 2, 8, 1);
    addBond(mol, 3, 9, 1);
    addBond(mol, 3, 10, 1);
    addBond(mol, 3, 11, 1);
    addBond(mol, 4, 12, 1);
    addBond(mol, 4, 13, 1);
    addBond(mol, 4, 14, 1);

    centerMolecule(mol);
}

// Build MTBE - Methyl tert-butyl ether (C5H12O)
void buildMTBE(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "MTBE (C5H12O)");

    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_C);      // Central C
    addAtom(mol, 0.0f, 1.4f, 0.0f, ATOM_O);
    addAtom(mol, 0.0f, 2.8f, 0.0f, ATOM_C);      // OCH3
    addAtom(mol, 1.4f, -0.5f, 0.0f, ATOM_C);     // CH3
    addAtom(mol, -0.7f, -0.5f, 1.2f, ATOM_C);    // CH3
    addAtom(mol, -0.7f, -0.5f, -1.2f, ATOM_C);   // CH3
    // H on OCH3
    addAtom(mol, 0.5f, 3.3f, 0.87f, ATOM_H);
    addAtom(mol, 0.5f, 3.3f, -0.87f, ATOM_H);
    addAtom(mol, -1.0f, 3.1f, 0.0f, ATOM_H);
    // H on tert-butyl CH3 groups
    addAtom(mol, 2.0f, 0.4f, 0.0f, ATOM_H);
    addAtom(mol, 1.9f, -1.0f, 0.87f, ATOM_H);
    addAtom(mol, 1.9f, -1.0f, -0.87f, ATOM_H);
    addAtom(mol, -0.1f, 0.0f, 2.0f, ATOM_H);
    addAtom(mol, -0.2f, -1.5f, 1.2f, ATOM_H);
    addAtom(mol, -1.7f, -0.3f, 1.4f, ATOM_H);
    addAtom(mol, -0.1f, 0.0f, -2.0f, ATOM_H);
    addAtom(mol, -0.2f, -1.5f, -1.2f, ATOM_H);
    addAtom(mol, -1.7f, -0.3f, -1.4f, ATOM_H);

    addBond(mol, 0, 1, 1);
    addBond(mol, 1, 2, 1);
    addBond(mol, 0, 3, 1);
    addBond(mol, 0, 4, 1);
    addBond(mol, 0, 5, 1);
    addBond(mol, 2, 6, 1);
    addBond(mol, 2, 7, 1);
    addBond(mol, 2, 8, 1);
    addBond(mol, 3, 9, 1);
    addBond(mol, 3, 10, 1);
    addBond(mol, 3, 11, 1);
    addBond(mol, 4, 12, 1);
    addBond(mol, 4, 13, 1);
    addBond(mol, 4, 14, 1);
    addBond(mol, 5, 15, 1);
    addBond(mol, 5, 16, 1);
    addBond(mol, 5, 17, 1);

    centerMolecule(mol);
}

// Build THF - Tetrahydrofuran (C4H8O)
void buildTHF(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "THF (C4H8O)");

    float r = 1.2f;
    // 5-membered ring
    addAtom(mol, r * cosf(0), r * sinf(0), 0.0f, ATOM_O);
    addAtom(mol, r * cosf(2*PI/5), r * sinf(2*PI/5), 0.2f, ATOM_C);
    addAtom(mol, r * cosf(4*PI/5), r * sinf(4*PI/5), -0.2f, ATOM_C);
    addAtom(mol, r * cosf(6*PI/5), r * sinf(6*PI/5), 0.2f, ATOM_C);
    addAtom(mol, r * cosf(8*PI/5), r * sinf(8*PI/5), -0.2f, ATOM_C);
    // H atoms
    addAtom(mol, r * cosf(2*PI/5) + 0.5f, r * sinf(2*PI/5) + 0.8f, 0.5f, ATOM_H);
    addAtom(mol, r * cosf(2*PI/5) + 0.5f, r * sinf(2*PI/5) + 0.8f, -0.5f, ATOM_H);
    addAtom(mol, r * cosf(4*PI/5) - 0.8f, r * sinf(4*PI/5) + 0.5f, 0.5f, ATOM_H);
    addAtom(mol, r * cosf(4*PI/5) - 0.8f, r * sinf(4*PI/5) + 0.5f, -0.5f, ATOM_H);
    addAtom(mol, r * cosf(6*PI/5) - 0.8f, r * sinf(6*PI/5) - 0.5f, 0.5f, ATOM_H);
    addAtom(mol, r * cosf(6*PI/5) - 0.8f, r * sinf(6*PI/5) - 0.5f, -0.5f, ATOM_H);
    addAtom(mol, r * cosf(8*PI/5) + 0.5f, r * sinf(8*PI/5) - 0.8f, 0.5f, ATOM_H);
    addAtom(mol, r * cosf(8*PI/5) + 0.5f, r * sinf(8*PI/5) - 0.8f, -0.5f, ATOM_H);

    addBond(mol, 0, 1, 1);
    addBond(mol, 1, 2, 1);
    addBond(mol, 2, 3, 1);
    addBond(mol, 3, 4, 1);
    addBond(mol, 4, 0, 1);
    addBond(mol, 1, 5, 1);
    addBond(mol, 1, 6, 1);
    addBond(mol, 2, 7, 1);
    addBond(mol, 2, 8, 1);
    addBond(mol, 3, 9, 1);
    addBond(mol, 3, 10, 1);
    addBond(mol, 4, 11, 1);
    addBond(mol, 4, 12, 1);

    centerMolecule(mol);
}

// Build 1,4-Dioxane (C4H8O2)
void buildDioxane(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "1,4-Dioxane (C4H8O2)");

    float r = 1.4f;
    // 6-membered ring with O at 1 and 4 positions
    addAtom(mol, r * cosf(0), r * sinf(0), 0.3f, ATOM_O);
    addAtom(mol, r * cosf(PI/3), r * sinf(PI/3), -0.3f, ATOM_C);
    addAtom(mol, r * cosf(2*PI/3), r * sinf(2*PI/3), 0.3f, ATOM_C);
    addAtom(mol, r * cosf(PI), r * sinf(PI), -0.3f, ATOM_O);
    addAtom(mol, r * cosf(4*PI/3), r * sinf(4*PI/3), 0.3f, ATOM_C);
    addAtom(mol, r * cosf(5*PI/3), r * sinf(5*PI/3), -0.3f, ATOM_C);
    // H atoms
    addAtom(mol, r * cosf(PI/3) * 1.6f, r * sinf(PI/3) * 1.6f, 0.5f, ATOM_H);
    addAtom(mol, r * cosf(PI/3) * 1.6f, r * sinf(PI/3) * 1.6f, -1.1f, ATOM_H);
    addAtom(mol, r * cosf(2*PI/3) * 1.6f, r * sinf(2*PI/3) * 1.6f, 1.1f, ATOM_H);
    addAtom(mol, r * cosf(2*PI/3) * 1.6f, r * sinf(2*PI/3) * 1.6f, -0.5f, ATOM_H);
    addAtom(mol, r * cosf(4*PI/3) * 1.6f, r * sinf(4*PI/3) * 1.6f, 1.1f, ATOM_H);
    addAtom(mol, r * cosf(4*PI/3) * 1.6f, r * sinf(4*PI/3) * 1.6f, -0.5f, ATOM_H);
    addAtom(mol, r * cosf(5*PI/3) * 1.6f, r * sinf(5*PI/3) * 1.6f, 0.5f, ATOM_H);
    addAtom(mol, r * cosf(5*PI/3) * 1.6f, r * sinf(5*PI/3) * 1.6f, -1.1f, ATOM_H);

    addBond(mol, 0, 1, 1);
    addBond(mol, 1, 2, 1);
    addBond(mol, 2, 3, 1);
    addBond(mol, 3, 4, 1);
    addBond(mol, 4, 5, 1);
    addBond(mol, 5, 0, 1);
    addBond(mol, 1, 6, 1);
    addBond(mol, 1, 7, 1);
    addBond(mol, 2, 8, 1);
    addBond(mol, 2, 9, 1);
    addBond(mol, 4, 10, 1);
    addBond(mol, 4, 11, 1);
    addBond(mol, 5, 12, 1);
    addBond(mol, 5, 13, 1);

    centerMolecule(mol);
}

// Build DMF - Dimethylformamide (C3H7NO)
void buildDMF(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "DMF (C3H7NO)");

    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_C);      // C=O
    addAtom(mol, 0.0f, 1.2f, 0.0f, ATOM_O);      // =O
    addAtom(mol, 1.2f, -0.6f, 0.0f, ATOM_N);     // N
    addAtom(mol, -1.2f, -0.6f, 0.0f, ATOM_H);    // H on C
    addAtom(mol, 1.2f, -2.1f, 0.0f, ATOM_C);     // CH3
    addAtom(mol, 2.4f, 0.1f, 0.0f, ATOM_C);      // CH3
    // H on CH3 groups
    addAtom(mol, 0.7f, -2.6f, 0.87f, ATOM_H);
    addAtom(mol, 0.7f, -2.6f, -0.87f, ATOM_H);
    addAtom(mol, 2.2f, -2.4f, 0.0f, ATOM_H);
    addAtom(mol, 2.9f, -0.4f, 0.87f, ATOM_H);
    addAtom(mol, 2.9f, -0.4f, -0.87f, ATOM_H);
    addAtom(mol, 2.9f, 1.0f, 0.0f, ATOM_H);

    addBond(mol, 0, 1, 2);
    addBond(mol, 0, 2, 1);
    addBond(mol, 0, 3, 1);
    addBond(mol, 2, 4, 1);
    addBond(mol, 2, 5, 1);
    addBond(mol, 4, 6, 1);
    addBond(mol, 4, 7, 1);
    addBond(mol, 4, 8, 1);
    addBond(mol, 5, 9, 1);
    addBond(mol, 5, 10, 1);
    addBond(mol, 5, 11, 1);

    centerMolecule(mol);
}

// Build Carbon Tetrachloride (CCl4)
void buildCarbonTetrachloride(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Carbon Tetrachloride (CCl4)");

    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_C);
    addAtom(mol, 1.2f, 0.85f, 0.0f, ATOM_CL);
    addAtom(mol, -1.2f, 0.85f, 0.0f, ATOM_CL);
    addAtom(mol, 0.0f, -0.6f, 1.2f, ATOM_CL);
    addAtom(mol, 0.0f, -0.6f, -1.2f, ATOM_CL);

    addBond(mol, 0, 1, 1);
    addBond(mol, 0, 2, 1);
    addBond(mol, 0, 3, 1);
    addBond(mol, 0, 4, 1);

    centerMolecule(mol);
}

// Build Methyl Acetate (C3H6O2)
void buildMethylAcetate(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Methyl Acetate (C3H6O2)");

    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_C);      // C=O
    addAtom(mol, -1.5f, 0.0f, 0.0f, ATOM_C);     // CH3 (acetyl)
    addAtom(mol, 0.6f, 1.1f, 0.0f, ATOM_O);      // =O
    addAtom(mol, 0.8f, -1.2f, 0.0f, ATOM_O);     // O-
    addAtom(mol, 2.2f, -1.2f, 0.0f, ATOM_C);     // CH3 (methyl)
    // H on acetyl CH3
    addAtom(mol, -2.0f, 1.0f, 0.0f, ATOM_H);
    addAtom(mol, -2.0f, -0.5f, 0.87f, ATOM_H);
    addAtom(mol, -2.0f, -0.5f, -0.87f, ATOM_H);
    // H on methyl CH3
    addAtom(mol, 2.7f, -0.2f, 0.0f, ATOM_H);
    addAtom(mol, 2.7f, -1.7f, 0.87f, ATOM_H);
    addAtom(mol, 2.7f, -1.7f, -0.87f, ATOM_H);

    addBond(mol, 0, 1, 1);
    addBond(mol, 0, 2, 2);
    addBond(mol, 0, 3, 1);
    addBond(mol, 3, 4, 1);
    addBond(mol, 1, 5, 1);
    addBond(mol, 1, 6, 1);
    addBond(mol, 1, 7, 1);
    addBond(mol, 4, 8, 1);
    addBond(mol, 4, 9, 1);
    addBond(mol, 4, 10, 1);

    centerMolecule(mol);
}

// Build Acetic Anhydride (C4H6O3)
void buildAceticAnhydride(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Acetic Anhydride (C4H6O3)");

    addAtom(mol, -1.5f, 0.0f, 0.0f, ATOM_C);     // C=O (left)
    addAtom(mol, -2.1f, 1.1f, 0.0f, ATOM_O);     // =O
    addAtom(mol, -3.0f, -0.5f, 0.0f, ATOM_C);    // CH3
    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_O);      // Central O
    addAtom(mol, 1.5f, 0.0f, 0.0f, ATOM_C);      // C=O (right)
    addAtom(mol, 2.1f, 1.1f, 0.0f, ATOM_O);      // =O
    addAtom(mol, 3.0f, -0.5f, 0.0f, ATOM_C);     // CH3
    // H atoms
    addAtom(mol, -3.5f, 0.5f, 0.0f, ATOM_H);
    addAtom(mol, -3.5f, -1.0f, 0.87f, ATOM_H);
    addAtom(mol, -3.5f, -1.0f, -0.87f, ATOM_H);
    addAtom(mol, 3.5f, 0.5f, 0.0f, ATOM_H);
    addAtom(mol, 3.5f, -1.0f, 0.87f, ATOM_H);
    addAtom(mol, 3.5f, -1.0f, -0.87f, ATOM_H);

    addBond(mol, 0, 1, 2);
    addBond(mol, 0, 2, 1);
    addBond(mol, 0, 3, 1);
    addBond(mol, 3, 4, 1);
    addBond(mol, 4, 5, 2);
    addBond(mol, 4, 6, 1);
    addBond(mol, 2, 7, 1);
    addBond(mol, 2, 8, 1);
    addBond(mol, 2, 9, 1);
    addBond(mol, 6, 10, 1);
    addBond(mol, 6, 11, 1);
    addBond(mol, 6, 12, 1);

    centerMolecule(mol);
}

// Build Propionic Acid (C3H6O2)
void buildPropionicAcid(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Propionic Acid (C3H6O2)");

    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_C);      // CH3
    addAtom(mol, 1.5f, 0.0f, 0.0f, ATOM_C);      // CH2
    addAtom(mol, 3.0f, 0.0f, 0.0f, ATOM_C);      // C=O
    addAtom(mol, 3.6f, 1.1f, 0.0f, ATOM_O);      // =O
    addAtom(mol, 3.6f, -1.1f, 0.0f, ATOM_O);     // OH
    addAtom(mol, 4.5f, -1.0f, 0.0f, ATOM_H);
    addAtom(mol, -0.5f, 1.0f, 0.0f, ATOM_H);
    addAtom(mol, -0.5f, -0.5f, 0.87f, ATOM_H);
    addAtom(mol, -0.5f, -0.5f, -0.87f, ATOM_H);
    addAtom(mol, 1.5f, 0.6f, 0.9f, ATOM_H);
    addAtom(mol, 1.5f, 0.6f, -0.9f, ATOM_H);

    addBond(mol, 0, 1, 1);
    addBond(mol, 1, 2, 1);
    addBond(mol, 2, 3, 2);
    addBond(mol, 2, 4, 1);
    addBond(mol, 4, 5, 1);
    addBond(mol, 0, 6, 1);
    addBond(mol, 0, 7, 1);
    addBond(mol, 0, 8, 1);
    addBond(mol, 1, 9, 1);
    addBond(mol, 1, 10, 1);

    centerMolecule(mol);
}

// Build Butyric Acid (C4H8O2)
void buildButyricAcid(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Butyric Acid (C4H8O2)");

    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_C);      // CH3
    addAtom(mol, 1.5f, 0.0f, 0.0f, ATOM_C);      // CH2
    addAtom(mol, 3.0f, 0.0f, 0.0f, ATOM_C);      // CH2
    addAtom(mol, 4.5f, 0.0f, 0.0f, ATOM_C);      // C=O
    addAtom(mol, 5.1f, 1.1f, 0.0f, ATOM_O);      // =O
    addAtom(mol, 5.1f, -1.1f, 0.0f, ATOM_O);     // OH
    addAtom(mol, 6.0f, -1.0f, 0.0f, ATOM_H);
    addAtom(mol, -0.5f, 1.0f, 0.0f, ATOM_H);
    addAtom(mol, -0.5f, -0.5f, 0.87f, ATOM_H);
    addAtom(mol, -0.5f, -0.5f, -0.87f, ATOM_H);
    addAtom(mol, 1.5f, 0.6f, 0.9f, ATOM_H);
    addAtom(mol, 1.5f, 0.6f, -0.9f, ATOM_H);
    addAtom(mol, 3.0f, 0.6f, 0.9f, ATOM_H);
    addAtom(mol, 3.0f, 0.6f, -0.9f, ATOM_H);

    addBond(mol, 0, 1, 1);
    addBond(mol, 1, 2, 1);
    addBond(mol, 2, 3, 1);
    addBond(mol, 3, 4, 2);
    addBond(mol, 3, 5, 1);
    addBond(mol, 5, 6, 1);
    addBond(mol, 0, 7, 1);
    addBond(mol, 0, 8, 1);
    addBond(mol, 0, 9, 1);
    addBond(mol, 1, 10, 1);
    addBond(mol, 1, 11, 1);
    addBond(mol, 2, 12, 1);
    addBond(mol, 2, 13, 1);

    centerMolecule(mol);
}

// Build Succinic Acid (C4H6O4)
void buildSuccinicAcid(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Succinic Acid (C4H6O4)");

    addAtom(mol, -1.5f, 0.0f, 0.0f, ATOM_C);     // COOH
    addAtom(mol, -2.1f, 1.1f, 0.0f, ATOM_O);
    addAtom(mol, -2.1f, -1.1f, 0.0f, ATOM_O);
    addAtom(mol, -3.0f, -1.0f, 0.0f, ATOM_H);
    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_C);      // CH2
    addAtom(mol, 1.5f, 0.0f, 0.0f, ATOM_C);      // CH2
    addAtom(mol, 3.0f, 0.0f, 0.0f, ATOM_C);      // COOH
    addAtom(mol, 3.6f, 1.1f, 0.0f, ATOM_O);
    addAtom(mol, 3.6f, -1.1f, 0.0f, ATOM_O);
    addAtom(mol, 4.5f, -1.0f, 0.0f, ATOM_H);
    addAtom(mol, 0.0f, 0.6f, 0.9f, ATOM_H);
    addAtom(mol, 0.0f, 0.6f, -0.9f, ATOM_H);
    addAtom(mol, 1.5f, 0.6f, 0.9f, ATOM_H);
    addAtom(mol, 1.5f, 0.6f, -0.9f, ATOM_H);

    addBond(mol, 0, 1, 2);
    addBond(mol, 0, 2, 1);
    addBond(mol, 2, 3, 1);
    addBond(mol, 0, 4, 1);
    addBond(mol, 4, 5, 1);
    addBond(mol, 5, 6, 1);
    addBond(mol, 6, 7, 2);
    addBond(mol, 6, 8, 1);
    addBond(mol, 8, 9, 1);
    addBond(mol, 4, 10, 1);
    addBond(mol, 4, 11, 1);
    addBond(mol, 5, 12, 1);
    addBond(mol, 5, 13, 1);

    centerMolecule(mol);
}

// Build Benzaldehyde (C7H6O)
void buildBenzaldehyde(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Benzaldehyde (C7H6O)");

    float r = 1.4f;
    for (int i = 0; i < 6; i++) {
        float angle = i * PI / 3.0f;
        addAtom(mol, r * cosf(angle), r * sinf(angle), 0.0f, ATOM_C);
    }
    // CHO group
    addAtom(mol, 2.6f, 0.0f, 0.0f, ATOM_C);
    addAtom(mol, 3.2f, 1.1f, 0.0f, ATOM_O);
    addAtom(mol, 3.2f, -1.0f, 0.0f, ATOM_H);
    // H on benzene
    float rH = 2.4f;
    for (int i = 1; i < 6; i++) {
        float angle = i * PI / 3.0f;
        addAtom(mol, rH * cosf(angle), rH * sinf(angle), 0.0f, ATOM_H);
    }

    for (int i = 0; i < 6; i++) {
        addBond(mol, i, (i + 1) % 6, (i % 2 == 0) ? 2 : 1);
    }
    addBond(mol, 0, 6, 1);
    addBond(mol, 6, 7, 2);
    addBond(mol, 6, 8, 1);
    for (int i = 1; i < 6; i++) {
        addBond(mol, i, 8 + i, 1);
    }

    centerMolecule(mol);
}

// Build Bromobenzene (C6H5Br)
void buildBromobenzene(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Bromobenzene (C6H5Br)");

    float r = 1.4f;
    for (int i = 0; i < 6; i++) {
        float angle = i * PI / 3.0f;
        addAtom(mol, r * cosf(angle), r * sinf(angle), 0.0f, ATOM_C);
    }
    addAtom(mol, 2.8f, 0.0f, 0.0f, ATOM_BR);
    float rH = 2.4f;
    for (int i = 1; i < 6; i++) {
        float angle = i * PI / 3.0f;
        addAtom(mol, rH * cosf(angle), rH * sinf(angle), 0.0f, ATOM_H);
    }

    for (int i = 0; i < 6; i++) {
        addBond(mol, i, (i + 1) % 6, (i % 2 == 0) ? 2 : 1);
    }
    addBond(mol, 0, 6, 1);
    for (int i = 1; i < 6; i++) {
        addBond(mol, i, 6 + i, 1);
    }

    centerMolecule(mol);
}

// Build p-Xylene (C8H10)
void buildPXylene(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "p-Xylene (C8H10)");

    float r = 1.4f;
    for (int i = 0; i < 6; i++) {
        float angle = i * PI / 3.0f;
        addAtom(mol, r * cosf(angle), r * sinf(angle), 0.0f, ATOM_C);
    }
    // CH3 at positions 1 and 4
    addAtom(mol, 2.6f, 0.0f, 0.0f, ATOM_C);
    addAtom(mol, -2.6f, 0.0f, 0.0f, ATOM_C);
    // H on CH3 groups
    addAtom(mol, 3.2f, 0.9f, 0.0f, ATOM_H);
    addAtom(mol, 3.2f, -0.45f, 0.8f, ATOM_H);
    addAtom(mol, 3.2f, -0.45f, -0.8f, ATOM_H);
    addAtom(mol, -3.2f, 0.9f, 0.0f, ATOM_H);
    addAtom(mol, -3.2f, -0.45f, 0.8f, ATOM_H);
    addAtom(mol, -3.2f, -0.45f, -0.8f, ATOM_H);
    // H on ring
    float rH = 2.4f;
    addAtom(mol, rH * cosf(PI/3), rH * sinf(PI/3), 0.0f, ATOM_H);
    addAtom(mol, rH * cosf(2*PI/3), rH * sinf(2*PI/3), 0.0f, ATOM_H);
    addAtom(mol, rH * cosf(4*PI/3), rH * sinf(4*PI/3), 0.0f, ATOM_H);
    addAtom(mol, rH * cosf(5*PI/3), rH * sinf(5*PI/3), 0.0f, ATOM_H);

    for (int i = 0; i < 6; i++) {
        addBond(mol, i, (i + 1) % 6, (i % 2 == 0) ? 2 : 1);
    }
    addBond(mol, 0, 6, 1);
    addBond(mol, 3, 7, 1);
    addBond(mol, 6, 8, 1);
    addBond(mol, 6, 9, 1);
    addBond(mol, 6, 10, 1);
    addBond(mol, 7, 11, 1);
    addBond(mol, 7, 12, 1);
    addBond(mol, 7, 13, 1);
    addBond(mol, 1, 14, 1);
    addBond(mol, 2, 15, 1);
    addBond(mol, 4, 16, 1);
    addBond(mol, 5, 17, 1);

    centerMolecule(mol);
}

// Build Anisole (C7H8O)
void buildAnisole(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Anisole (C7H8O)");

    float r = 1.4f;
    for (int i = 0; i < 6; i++) {
        float angle = i * PI / 3.0f;
        addAtom(mol, r * cosf(angle), r * sinf(angle), 0.0f, ATOM_C);
    }
    addAtom(mol, 2.4f, 0.0f, 0.0f, ATOM_O);
    addAtom(mol, 3.6f, 0.0f, 0.0f, ATOM_C);      // CH3
    addAtom(mol, 4.1f, 1.0f, 0.0f, ATOM_H);
    addAtom(mol, 4.1f, -0.5f, 0.87f, ATOM_H);
    addAtom(mol, 4.1f, -0.5f, -0.87f, ATOM_H);
    float rH = 2.4f;
    for (int i = 1; i < 6; i++) {
        float angle = i * PI / 3.0f;
        addAtom(mol, rH * cosf(angle), rH * sinf(angle), 0.0f, ATOM_H);
    }

    for (int i = 0; i < 6; i++) {
        addBond(mol, i, (i + 1) % 6, (i % 2 == 0) ? 2 : 1);
    }
    addBond(mol, 0, 6, 1);
    addBond(mol, 6, 7, 1);
    addBond(mol, 7, 8, 1);
    addBond(mol, 7, 9, 1);
    addBond(mol, 7, 10, 1);
    for (int i = 1; i < 6; i++) {
        addBond(mol, i, 10 + i, 1);
    }

    centerMolecule(mol);
}

// Build Phenylacetylene (C8H6)
void buildPhenylacetylene(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Phenylacetylene (C8H6)");

    float r = 1.4f;
    for (int i = 0; i < 6; i++) {
        float angle = i * PI / 3.0f;
        addAtom(mol, r * cosf(angle), r * sinf(angle), 0.0f, ATOM_C);
    }
    addAtom(mol, 2.6f, 0.0f, 0.0f, ATOM_C);      // C triple bond
    addAtom(mol, 3.8f, 0.0f, 0.0f, ATOM_C);
    addAtom(mol, 4.85f, 0.0f, 0.0f, ATOM_H);
    float rH = 2.4f;
    for (int i = 1; i < 6; i++) {
        float angle = i * PI / 3.0f;
        addAtom(mol, rH * cosf(angle), rH * sinf(angle), 0.0f, ATOM_H);
    }

    for (int i = 0; i < 6; i++) {
        addBond(mol, i, (i + 1) % 6, (i % 2 == 0) ? 2 : 1);
    }
    addBond(mol, 0, 6, 1);
    addBond(mol, 6, 7, 3);
    addBond(mol, 7, 8, 1);
    for (int i = 1; i < 6; i++) {
        addBond(mol, i, 8 + i, 1);
    }

    centerMolecule(mol);
}

// Build Fructose (C6H12O6)
void buildFructose(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Fructose (C6H12O6)");

    // Furanose ring form
    float r = 1.3f;
    addAtom(mol, r, 0.0f, 0.3f, ATOM_C);          // C2
    addAtom(mol, r * 0.31f, r * 0.95f, -0.2f, ATOM_C);   // C3
    addAtom(mol, -r * 0.81f, r * 0.59f, 0.3f, ATOM_C);   // C4
    addAtom(mol, -r * 0.81f, -r * 0.59f, -0.2f, ATOM_C); // C5
    addAtom(mol, r * 0.31f, -r * 0.95f, 0.1f, ATOM_O);   // Ring O
    // C1 and C6
    addAtom(mol, r + 1.3f, 0.3f, 0.0f, ATOM_C);   // C1 (CH2OH)
    addAtom(mol, -r * 0.81f - 1.2f, -r * 0.59f, 0.0f, ATOM_C);  // C6 (CH2OH)
    // OH groups
    addAtom(mol, r * 1.5f, 0.0f, 1.0f, ATOM_O);
    addAtom(mol, r * 0.5f, r * 1.5f, -0.8f, ATOM_O);
    addAtom(mol, -r * 1.3f, r * 1.0f, 1.0f, ATOM_O);
    addAtom(mol, r + 1.8f, 1.2f, 0.0f, ATOM_O);
    addAtom(mol, -r * 0.81f - 1.8f, -r * 0.59f - 1.0f, 0.0f, ATOM_O);
    // Simplified H atoms
    addAtom(mol, r * 1.5f + 0.8f, 0.0f, 1.0f, ATOM_H);
    addAtom(mol, r * 0.5f + 0.8f, r * 1.5f, -0.8f, ATOM_H);
    addAtom(mol, -r * 1.3f - 0.8f, r * 1.0f, 1.0f, ATOM_H);
    addAtom(mol, r + 1.8f + 0.8f, 1.2f, 0.0f, ATOM_H);
    addAtom(mol, -r * 0.81f - 1.8f - 0.8f, -r * 0.59f - 1.0f, 0.0f, ATOM_H);

    addBond(mol, 0, 1, 1);
    addBond(mol, 1, 2, 1);
    addBond(mol, 2, 3, 1);
    addBond(mol, 3, 4, 1);
    addBond(mol, 4, 0, 1);
    addBond(mol, 0, 5, 1);
    addBond(mol, 3, 6, 1);
    addBond(mol, 0, 7, 1);
    addBond(mol, 1, 8, 1);
    addBond(mol, 2, 9, 1);
    addBond(mol, 5, 10, 1);
    addBond(mol, 6, 11, 1);
    addBond(mol, 7, 12, 1);
    addBond(mol, 8, 13, 1);
    addBond(mol, 9, 14, 1);
    addBond(mol, 10, 15, 1);
    addBond(mol, 11, 16, 1);

    centerMolecule(mol);
}

// Build Ribose (C5H10O5)
void buildRibose(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Ribose (C5H10O5)");

    float r = 1.2f;
    // Furanose ring
    addAtom(mol, r, 0.0f, 0.2f, ATOM_C);          // C1
    addAtom(mol, r * 0.31f, r * 0.95f, -0.2f, ATOM_C);   // C2
    addAtom(mol, -r * 0.81f, r * 0.59f, 0.2f, ATOM_C);   // C3
    addAtom(mol, -r * 0.81f, -r * 0.59f, -0.2f, ATOM_C); // C4
    addAtom(mol, r * 0.31f, -r * 0.95f, 0.1f, ATOM_O);   // Ring O
    addAtom(mol, -r * 0.81f - 1.2f, -r * 0.59f, 0.0f, ATOM_C);  // C5
    // OH groups
    addAtom(mol, r + 0.9f, 0.7f, 0.0f, ATOM_O);
    addAtom(mol, r * 0.5f, r * 1.5f, -0.5f, ATOM_O);
    addAtom(mol, -r * 1.3f, r * 1.0f, 0.5f, ATOM_O);
    addAtom(mol, -r * 0.81f - 1.8f, -r * 0.59f - 0.8f, 0.0f, ATOM_O);
    // H on OH
    addAtom(mol, r + 1.7f, 0.5f, 0.0f, ATOM_H);
    addAtom(mol, r * 0.5f + 0.8f, r * 1.5f, -0.5f, ATOM_H);
    addAtom(mol, -r * 1.3f - 0.8f, r * 1.0f, 0.5f, ATOM_H);
    addAtom(mol, -r * 0.81f - 2.6f, -r * 0.59f - 0.6f, 0.0f, ATOM_H);

    addBond(mol, 0, 1, 1);
    addBond(mol, 1, 2, 1);
    addBond(mol, 2, 3, 1);
    addBond(mol, 3, 4, 1);
    addBond(mol, 4, 0, 1);
    addBond(mol, 3, 5, 1);
    addBond(mol, 0, 6, 1);
    addBond(mol, 1, 7, 1);
    addBond(mol, 2, 8, 1);
    addBond(mol, 5, 9, 1);
    addBond(mol, 6, 10, 1);
    addBond(mol, 7, 11, 1);
    addBond(mol, 8, 12, 1);
    addBond(mol, 9, 13, 1);

    centerMolecule(mol);
}

// Build Deoxyribose (C5H10O4)
void buildDeoxyribose(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Deoxyribose (C5H10O4)");

    float r = 1.2f;
    addAtom(mol, r, 0.0f, 0.2f, ATOM_C);          // C1
    addAtom(mol, r * 0.31f, r * 0.95f, -0.2f, ATOM_C);   // C2
    addAtom(mol, -r * 0.81f, r * 0.59f, 0.2f, ATOM_C);   // C3
    addAtom(mol, -r * 0.81f, -r * 0.59f, -0.2f, ATOM_C); // C4
    addAtom(mol, r * 0.31f, -r * 0.95f, 0.1f, ATOM_O);   // Ring O
    addAtom(mol, -r * 0.81f - 1.2f, -r * 0.59f, 0.0f, ATOM_C);  // C5
    // OH groups (one less than ribose)
    addAtom(mol, r + 0.9f, 0.7f, 0.0f, ATOM_O);
    addAtom(mol, -r * 1.3f, r * 1.0f, 0.5f, ATOM_O);
    addAtom(mol, -r * 0.81f - 1.8f, -r * 0.59f - 0.8f, 0.0f, ATOM_O);
    // H atoms
    addAtom(mol, r + 1.7f, 0.5f, 0.0f, ATOM_H);
    addAtom(mol, -r * 1.3f - 0.8f, r * 1.0f, 0.5f, ATOM_H);
    addAtom(mol, -r * 0.81f - 2.6f, -r * 0.59f - 0.6f, 0.0f, ATOM_H);
    addAtom(mol, r * 0.5f, r * 1.5f, 0.5f, ATOM_H);
    addAtom(mol, r * 0.5f, r * 1.5f, -0.9f, ATOM_H);

    addBond(mol, 0, 1, 1);
    addBond(mol, 1, 2, 1);
    addBond(mol, 2, 3, 1);
    addBond(mol, 3, 4, 1);
    addBond(mol, 4, 0, 1);
    addBond(mol, 3, 5, 1);
    addBond(mol, 0, 6, 1);
    addBond(mol, 2, 7, 1);
    addBond(mol, 5, 8, 1);
    addBond(mol, 6, 9, 1);
    addBond(mol, 7, 10, 1);
    addBond(mol, 8, 11, 1);
    addBond(mol, 1, 12, 1);
    addBond(mol, 1, 13, 1);

    centerMolecule(mol);
}

// Build Isoleucine (C6H13NO2)
void buildIsoleucine(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Isoleucine (C6H13NO2)");

    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_C);      // Alpha C
    addAtom(mol, -1.3f, 0.5f, 0.0f, ATOM_N);     // NH2
    addAtom(mol, 1.3f, 0.7f, 0.0f, ATOM_C);      // COOH C
    addAtom(mol, 1.3f, 2.0f, 0.0f, ATOM_O);      // =O
    addAtom(mol, 2.4f, 0.0f, 0.0f, ATOM_O);      // OH
    addAtom(mol, 0.0f, -1.5f, 0.0f, ATOM_C);     // Beta C (CH)
    addAtom(mol, -1.3f, -2.2f, 0.0f, ATOM_C);    // CH3
    addAtom(mol, 1.3f, -2.2f, 0.0f, ATOM_C);     // CH2
    addAtom(mol, 1.3f, -3.7f, 0.0f, ATOM_C);     // CH3
    // H atoms
    addAtom(mol, 0.0f, 0.5f, 0.9f, ATOM_H);
    addAtom(mol, -1.4f, 1.5f, 0.0f, ATOM_H);
    addAtom(mol, -2.1f, 0.0f, 0.0f, ATOM_H);
    addAtom(mol, 3.2f, 0.5f, 0.0f, ATOM_H);
    addAtom(mol, 0.0f, -2.1f, 0.9f, ATOM_H);
    addAtom(mol, -1.3f, -3.3f, 0.0f, ATOM_H);
    addAtom(mol, -2.0f, -1.8f, 0.8f, ATOM_H);
    addAtom(mol, -2.0f, -1.8f, -0.8f, ATOM_H);
    addAtom(mol, 2.0f, -1.8f, 0.8f, ATOM_H);
    addAtom(mol, 2.0f, -1.8f, -0.8f, ATOM_H);
    addAtom(mol, 1.3f, -4.3f, 0.9f, ATOM_H);
    addAtom(mol, 2.0f, -4.0f, -0.5f, ATOM_H);
    addAtom(mol, 0.4f, -4.0f, -0.5f, ATOM_H);

    addBond(mol, 0, 1, 1);
    addBond(mol, 0, 2, 1);
    addBond(mol, 2, 3, 2);
    addBond(mol, 2, 4, 1);
    addBond(mol, 0, 5, 1);
    addBond(mol, 5, 6, 1);
    addBond(mol, 5, 7, 1);
    addBond(mol, 7, 8, 1);
    addBond(mol, 0, 9, 1);
    addBond(mol, 1, 10, 1);
    addBond(mol, 1, 11, 1);
    addBond(mol, 4, 12, 1);
    addBond(mol, 5, 13, 1);
    addBond(mol, 6, 14, 1);
    addBond(mol, 6, 15, 1);
    addBond(mol, 6, 16, 1);
    addBond(mol, 7, 17, 1);
    addBond(mol, 7, 18, 1);
    addBond(mol, 8, 19, 1);
    addBond(mol, 8, 20, 1);
    addBond(mol, 8, 21, 1);

    centerMolecule(mol);
}

// Build Serine (C3H7NO3)
void buildSerine(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Serine (C3H7NO3)");

    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_C);      // Alpha C
    addAtom(mol, -1.3f, 0.5f, 0.0f, ATOM_N);     // NH2
    addAtom(mol, 1.3f, 0.7f, 0.0f, ATOM_C);      // COOH C
    addAtom(mol, 1.3f, 2.0f, 0.0f, ATOM_O);      // =O
    addAtom(mol, 2.4f, 0.0f, 0.0f, ATOM_O);      // OH
    addAtom(mol, 0.0f, -1.5f, 0.0f, ATOM_C);     // CH2
    addAtom(mol, 0.0f, -2.9f, 0.0f, ATOM_O);     // OH
    addAtom(mol, 0.0f, 0.5f, 0.9f, ATOM_H);
    addAtom(mol, -1.4f, 1.5f, 0.0f, ATOM_H);
    addAtom(mol, -2.1f, 0.0f, 0.0f, ATOM_H);
    addAtom(mol, 3.2f, 0.5f, 0.0f, ATOM_H);
    addAtom(mol, -0.9f, -1.8f, 0.0f, ATOM_H);
    addAtom(mol, 0.9f, -1.8f, 0.0f, ATOM_H);
    addAtom(mol, 0.8f, -3.3f, 0.0f, ATOM_H);

    addBond(mol, 0, 1, 1);
    addBond(mol, 0, 2, 1);
    addBond(mol, 2, 3, 2);
    addBond(mol, 2, 4, 1);
    addBond(mol, 0, 5, 1);
    addBond(mol, 5, 6, 1);
    addBond(mol, 0, 7, 1);
    addBond(mol, 1, 8, 1);
    addBond(mol, 1, 9, 1);
    addBond(mol, 4, 10, 1);
    addBond(mol, 5, 11, 1);
    addBond(mol, 5, 12, 1);
    addBond(mol, 6, 13, 1);

    centerMolecule(mol);
}

// Build Threonine (C4H9NO3)
void buildThreonine(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Threonine (C4H9NO3)");

    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_C);      // Alpha C
    addAtom(mol, -1.3f, 0.5f, 0.0f, ATOM_N);     // NH2
    addAtom(mol, 1.3f, 0.7f, 0.0f, ATOM_C);      // COOH C
    addAtom(mol, 1.3f, 2.0f, 0.0f, ATOM_O);      // =O
    addAtom(mol, 2.4f, 0.0f, 0.0f, ATOM_O);      // OH
    addAtom(mol, 0.0f, -1.5f, 0.0f, ATOM_C);     // CH (with OH)
    addAtom(mol, 0.0f, -2.9f, 0.0f, ATOM_C);     // CH3
    addAtom(mol, 1.3f, -1.8f, 0.0f, ATOM_O);     // OH
    addAtom(mol, 0.0f, 0.5f, 0.9f, ATOM_H);
    addAtom(mol, -1.4f, 1.5f, 0.0f, ATOM_H);
    addAtom(mol, -2.1f, 0.0f, 0.0f, ATOM_H);
    addAtom(mol, 3.2f, 0.5f, 0.0f, ATOM_H);
    addAtom(mol, -0.9f, -1.8f, 0.0f, ATOM_H);
    addAtom(mol, 2.0f, -1.4f, 0.0f, ATOM_H);
    addAtom(mol, 0.0f, -3.5f, 0.9f, ATOM_H);
    addAtom(mol, 0.9f, -3.2f, -0.5f, ATOM_H);
    addAtom(mol, -0.9f, -3.2f, -0.5f, ATOM_H);

    addBond(mol, 0, 1, 1);
    addBond(mol, 0, 2, 1);
    addBond(mol, 2, 3, 2);
    addBond(mol, 2, 4, 1);
    addBond(mol, 0, 5, 1);
    addBond(mol, 5, 6, 1);
    addBond(mol, 5, 7, 1);
    addBond(mol, 0, 8, 1);
    addBond(mol, 1, 9, 1);
    addBond(mol, 1, 10, 1);
    addBond(mol, 4, 11, 1);
    addBond(mol, 5, 12, 1);
    addBond(mol, 7, 13, 1);
    addBond(mol, 6, 14, 1);
    addBond(mol, 6, 15, 1);
    addBond(mol, 6, 16, 1);

    centerMolecule(mol);
}

// Build Aspartic Acid (C4H7NO4)
void buildAsparticAcid(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Aspartic Acid (C4H7NO4)");

    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_C);      // Alpha C
    addAtom(mol, -1.3f, 0.5f, 0.0f, ATOM_N);     // NH2
    addAtom(mol, 1.3f, 0.7f, 0.0f, ATOM_C);      // COOH C
    addAtom(mol, 1.3f, 2.0f, 0.0f, ATOM_O);
    addAtom(mol, 2.4f, 0.0f, 0.0f, ATOM_O);
    addAtom(mol, 0.0f, -1.5f, 0.0f, ATOM_C);     // CH2
    addAtom(mol, 0.0f, -3.0f, 0.0f, ATOM_C);     // Side COOH
    addAtom(mol, 0.0f, -4.2f, 0.0f, ATOM_O);
    addAtom(mol, 1.1f, -3.3f, 0.0f, ATOM_O);
    addAtom(mol, 0.0f, 0.5f, 0.9f, ATOM_H);
    addAtom(mol, -1.4f, 1.5f, 0.0f, ATOM_H);
    addAtom(mol, -2.1f, 0.0f, 0.0f, ATOM_H);
    addAtom(mol, 3.2f, 0.5f, 0.0f, ATOM_H);
    addAtom(mol, -0.9f, -1.8f, 0.0f, ATOM_H);
    addAtom(mol, 0.9f, -1.8f, 0.0f, ATOM_H);
    addAtom(mol, 1.9f, -2.9f, 0.0f, ATOM_H);

    addBond(mol, 0, 1, 1);
    addBond(mol, 0, 2, 1);
    addBond(mol, 2, 3, 2);
    addBond(mol, 2, 4, 1);
    addBond(mol, 0, 5, 1);
    addBond(mol, 5, 6, 1);
    addBond(mol, 6, 7, 2);
    addBond(mol, 6, 8, 1);
    addBond(mol, 0, 9, 1);
    addBond(mol, 1, 10, 1);
    addBond(mol, 1, 11, 1);
    addBond(mol, 4, 12, 1);
    addBond(mol, 5, 13, 1);
    addBond(mol, 5, 14, 1);
    addBond(mol, 8, 15, 1);

    centerMolecule(mol);
}

// Build Glutamic Acid (C5H9NO4)
void buildGlutamicAcid(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Glutamic Acid (C5H9NO4)");

    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_C);      // Alpha C
    addAtom(mol, -1.3f, 0.5f, 0.0f, ATOM_N);     // NH2
    addAtom(mol, 1.3f, 0.7f, 0.0f, ATOM_C);      // COOH C
    addAtom(mol, 1.3f, 2.0f, 0.0f, ATOM_O);
    addAtom(mol, 2.4f, 0.0f, 0.0f, ATOM_O);
    addAtom(mol, 0.0f, -1.5f, 0.0f, ATOM_C);     // CH2
    addAtom(mol, 0.0f, -3.0f, 0.0f, ATOM_C);     // CH2
    addAtom(mol, 0.0f, -4.5f, 0.0f, ATOM_C);     // Side COOH
    addAtom(mol, 0.0f, -5.7f, 0.0f, ATOM_O);
    addAtom(mol, 1.1f, -4.8f, 0.0f, ATOM_O);
    addAtom(mol, 0.0f, 0.5f, 0.9f, ATOM_H);
    addAtom(mol, -1.4f, 1.5f, 0.0f, ATOM_H);
    addAtom(mol, -2.1f, 0.0f, 0.0f, ATOM_H);
    addAtom(mol, 3.2f, 0.5f, 0.0f, ATOM_H);
    addAtom(mol, -0.9f, -1.8f, 0.0f, ATOM_H);
    addAtom(mol, 0.9f, -1.8f, 0.0f, ATOM_H);
    addAtom(mol, -0.9f, -3.3f, 0.0f, ATOM_H);
    addAtom(mol, 0.9f, -3.3f, 0.0f, ATOM_H);
    addAtom(mol, 1.9f, -4.4f, 0.0f, ATOM_H);

    addBond(mol, 0, 1, 1);
    addBond(mol, 0, 2, 1);
    addBond(mol, 2, 3, 2);
    addBond(mol, 2, 4, 1);
    addBond(mol, 0, 5, 1);
    addBond(mol, 5, 6, 1);
    addBond(mol, 6, 7, 1);
    addBond(mol, 7, 8, 2);
    addBond(mol, 7, 9, 1);
    addBond(mol, 0, 10, 1);
    addBond(mol, 1, 11, 1);
    addBond(mol, 1, 12, 1);
    addBond(mol, 4, 13, 1);
    addBond(mol, 5, 14, 1);
    addBond(mol, 5, 15, 1);
    addBond(mol, 6, 16, 1);
    addBond(mol, 6, 17, 1);
    addBond(mol, 9, 18, 1);

    centerMolecule(mol);
}

// Build Lysine (C6H14N2O2)
void buildLysine(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Lysine (C6H14N2O2)");

    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_C);      // Alpha C
    addAtom(mol, -1.3f, 0.5f, 0.0f, ATOM_N);     // NH2
    addAtom(mol, 1.3f, 0.7f, 0.0f, ATOM_C);      // COOH C
    addAtom(mol, 1.3f, 2.0f, 0.0f, ATOM_O);
    addAtom(mol, 2.4f, 0.0f, 0.0f, ATOM_O);
    addAtom(mol, 0.0f, -1.5f, 0.0f, ATOM_C);     // CH2
    addAtom(mol, 0.0f, -3.0f, 0.0f, ATOM_C);     // CH2
    addAtom(mol, 0.0f, -4.5f, 0.0f, ATOM_C);     // CH2
    addAtom(mol, 0.0f, -6.0f, 0.0f, ATOM_C);     // CH2
    addAtom(mol, 0.0f, -7.5f, 0.0f, ATOM_N);     // NH2
    addAtom(mol, 0.0f, 0.5f, 0.9f, ATOM_H);
    addAtom(mol, -1.4f, 1.5f, 0.0f, ATOM_H);
    addAtom(mol, -2.1f, 0.0f, 0.0f, ATOM_H);
    addAtom(mol, 3.2f, 0.5f, 0.0f, ATOM_H);
    addAtom(mol, 0.8f, -8.0f, 0.0f, ATOM_H);
    addAtom(mol, -0.8f, -8.0f, 0.0f, ATOM_H);

    addBond(mol, 0, 1, 1);
    addBond(mol, 0, 2, 1);
    addBond(mol, 2, 3, 2);
    addBond(mol, 2, 4, 1);
    addBond(mol, 0, 5, 1);
    addBond(mol, 5, 6, 1);
    addBond(mol, 6, 7, 1);
    addBond(mol, 7, 8, 1);
    addBond(mol, 8, 9, 1);
    addBond(mol, 0, 10, 1);
    addBond(mol, 1, 11, 1);
    addBond(mol, 1, 12, 1);
    addBond(mol, 4, 13, 1);
    addBond(mol, 9, 14, 1);
    addBond(mol, 9, 15, 1);

    centerMolecule(mol);
}

// Build Histidine (C6H9N3O2)
void buildHistidine(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Histidine (C6H9N3O2)");

    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_C);      // Alpha C
    addAtom(mol, -1.3f, 0.5f, 0.0f, ATOM_N);     // NH2
    addAtom(mol, 1.3f, 0.7f, 0.0f, ATOM_C);      // COOH C
    addAtom(mol, 1.3f, 2.0f, 0.0f, ATOM_O);
    addAtom(mol, 2.4f, 0.0f, 0.0f, ATOM_O);
    addAtom(mol, 0.0f, -1.5f, 0.0f, ATOM_C);     // CH2
    // Imidazole ring
    addAtom(mol, 0.0f, -3.0f, 0.0f, ATOM_C);     // C
    addAtom(mol, -1.0f, -3.7f, 0.0f, ATOM_N);    // N
    addAtom(mol, -0.6f, -5.0f, 0.0f, ATOM_C);    // C
    addAtom(mol, 0.8f, -5.0f, 0.0f, ATOM_N);     // NH
    addAtom(mol, 1.1f, -3.6f, 0.0f, ATOM_C);     // C
    // H atoms
    addAtom(mol, 0.0f, 0.5f, 0.9f, ATOM_H);
    addAtom(mol, -1.4f, 1.5f, 0.0f, ATOM_H);
    addAtom(mol, -2.1f, 0.0f, 0.0f, ATOM_H);
    addAtom(mol, 3.2f, 0.5f, 0.0f, ATOM_H);
    addAtom(mol, -1.1f, -5.9f, 0.0f, ATOM_H);
    addAtom(mol, 1.3f, -5.8f, 0.0f, ATOM_H);
    addAtom(mol, 2.0f, -3.2f, 0.0f, ATOM_H);

    addBond(mol, 0, 1, 1);
    addBond(mol, 0, 2, 1);
    addBond(mol, 2, 3, 2);
    addBond(mol, 2, 4, 1);
    addBond(mol, 0, 5, 1);
    addBond(mol, 5, 6, 1);
    addBond(mol, 6, 7, 1);
    addBond(mol, 7, 8, 2);
    addBond(mol, 8, 9, 1);
    addBond(mol, 9, 10, 1);
    addBond(mol, 10, 6, 2);
    addBond(mol, 0, 11, 1);
    addBond(mol, 1, 12, 1);
    addBond(mol, 1, 13, 1);
    addBond(mol, 4, 14, 1);
    addBond(mol, 8, 15, 1);
    addBond(mol, 9, 16, 1);
    addBond(mol, 10, 17, 1);

    centerMolecule(mol);
}

// Build Phenylalanine (C9H11NO2)
void buildPhenylalanine(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Phenylalanine (C9H11NO2)");

    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_C);      // Alpha C
    addAtom(mol, -1.3f, 0.5f, 0.0f, ATOM_N);     // NH2
    addAtom(mol, 1.3f, 0.7f, 0.0f, ATOM_C);      // COOH C
    addAtom(mol, 1.3f, 2.0f, 0.0f, ATOM_O);
    addAtom(mol, 2.4f, 0.0f, 0.0f, ATOM_O);
    addAtom(mol, 0.0f, -1.5f, 0.0f, ATOM_C);     // CH2
    // Benzene ring
    float r = 1.4f;
    float cx = 0.0f, cy = -3.5f;
    for (int i = 0; i < 6; i++) {
        float angle = i * PI / 3.0f;
        addAtom(mol, cx + r * cosf(angle), cy + r * sinf(angle), 0.0f, ATOM_C);
    }
    // H atoms
    addAtom(mol, 0.0f, 0.5f, 0.9f, ATOM_H);
    addAtom(mol, -1.4f, 1.5f, 0.0f, ATOM_H);
    addAtom(mol, -2.1f, 0.0f, 0.0f, ATOM_H);
    addAtom(mol, 3.2f, 0.5f, 0.0f, ATOM_H);
    float rH = 2.4f;
    for (int i = 0; i < 6; i++) {
        float angle = i * PI / 3.0f;
        addAtom(mol, cx + rH * cosf(angle), cy + rH * sinf(angle), 0.0f, ATOM_H);
    }

    addBond(mol, 0, 1, 1);
    addBond(mol, 0, 2, 1);
    addBond(mol, 2, 3, 2);
    addBond(mol, 2, 4, 1);
    addBond(mol, 0, 5, 1);
    addBond(mol, 5, 6, 1);
    // Ring bonds
    for (int i = 0; i < 6; i++) {
        addBond(mol, 6 + i, 6 + (i + 1) % 6, (i % 2 == 0) ? 2 : 1);
    }
    addBond(mol, 0, 12, 1);
    addBond(mol, 1, 13, 1);
    addBond(mol, 1, 14, 1);
    addBond(mol, 4, 15, 1);
    for (int i = 0; i < 6; i++) {
        addBond(mol, 6 + i, 16 + i, 1);
    }

    centerMolecule(mol);
}

// Build Tyrosine (C9H11NO3)
void buildTyrosine(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Tyrosine (C9H11NO3)");

    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_C);      // Alpha C
    addAtom(mol, -1.3f, 0.5f, 0.0f, ATOM_N);     // NH2
    addAtom(mol, 1.3f, 0.7f, 0.0f, ATOM_C);      // COOH C
    addAtom(mol, 1.3f, 2.0f, 0.0f, ATOM_O);
    addAtom(mol, 2.4f, 0.0f, 0.0f, ATOM_O);
    addAtom(mol, 0.0f, -1.5f, 0.0f, ATOM_C);     // CH2
    // Benzene ring
    float r = 1.4f;
    float cx = 0.0f, cy = -3.5f;
    for (int i = 0; i < 6; i++) {
        float angle = i * PI / 3.0f;
        addAtom(mol, cx + r * cosf(angle), cy + r * sinf(angle), 0.0f, ATOM_C);
    }
    // OH on para position
    addAtom(mol, cx - 2.6f, cy, 0.0f, ATOM_O);
    addAtom(mol, cx - 3.4f, cy, 0.0f, ATOM_H);
    // H atoms
    addAtom(mol, 0.0f, 0.5f, 0.9f, ATOM_H);
    addAtom(mol, -1.4f, 1.5f, 0.0f, ATOM_H);
    addAtom(mol, -2.1f, 0.0f, 0.0f, ATOM_H);
    addAtom(mol, 3.2f, 0.5f, 0.0f, ATOM_H);
    float rH = 2.4f;
    for (int i = 0; i < 6; i++) {
        if (i == 3) continue;  // Skip para position (has OH)
        float angle = i * PI / 3.0f;
        addAtom(mol, cx + rH * cosf(angle), cy + rH * sinf(angle), 0.0f, ATOM_H);
    }

    addBond(mol, 0, 1, 1);
    addBond(mol, 0, 2, 1);
    addBond(mol, 2, 3, 2);
    addBond(mol, 2, 4, 1);
    addBond(mol, 0, 5, 1);
    addBond(mol, 5, 6, 1);
    for (int i = 0; i < 6; i++) {
        addBond(mol, 6 + i, 6 + (i + 1) % 6, (i % 2 == 0) ? 2 : 1);
    }
    addBond(mol, 9, 12, 1);  // C-OH
    addBond(mol, 12, 13, 1);
    addBond(mol, 0, 14, 1);
    addBond(mol, 1, 15, 1);
    addBond(mol, 1, 16, 1);
    addBond(mol, 4, 17, 1);
    int hIdx = 18;
    for (int i = 0; i < 6; i++) {
        if (i == 3) continue;
        addBond(mol, 6 + i, hIdx++, 1);
    }

    centerMolecule(mol);
}

// Build Tryptophan (C11H12N2O2)
void buildTryptophan(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Tryptophan (C11H12N2O2)");

    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_C);      // Alpha C
    addAtom(mol, -1.3f, 0.5f, 0.0f, ATOM_N);     // NH2
    addAtom(mol, 1.3f, 0.7f, 0.0f, ATOM_C);      // COOH C
    addAtom(mol, 1.3f, 2.0f, 0.0f, ATOM_O);
    addAtom(mol, 2.4f, 0.0f, 0.0f, ATOM_O);
    addAtom(mol, 0.0f, -1.5f, 0.0f, ATOM_C);     // CH2
    // Indole ring system
    addAtom(mol, 0.0f, -3.0f, 0.0f, ATOM_C);     // C3 pyrrole
    addAtom(mol, -1.0f, -3.7f, 0.0f, ATOM_C);    // C3a
    addAtom(mol, 1.1f, -3.7f, 0.0f, ATOM_C);     // C2
    addAtom(mol, 1.0f, -5.0f, 0.0f, ATOM_N);     // N1
    // Benzene fused ring
    addAtom(mol, -1.0f, -5.1f, 0.0f, ATOM_C);    // C7a
    addAtom(mol, -2.2f, -3.2f, 0.0f, ATOM_C);    // C4
    addAtom(mol, -3.3f, -4.0f, 0.0f, ATOM_C);    // C5
    addAtom(mol, -3.2f, -5.4f, 0.0f, ATOM_C);    // C6
    addAtom(mol, -2.0f, -6.0f, 0.0f, ATOM_C);    // C7
    // H atoms
    addAtom(mol, 0.0f, 0.5f, 0.9f, ATOM_H);
    addAtom(mol, -1.4f, 1.5f, 0.0f, ATOM_H);
    addAtom(mol, -2.1f, 0.0f, 0.0f, ATOM_H);
    addAtom(mol, 3.2f, 0.5f, 0.0f, ATOM_H);
    addAtom(mol, 2.0f, -3.3f, 0.0f, ATOM_H);
    addAtom(mol, 1.7f, -5.6f, 0.0f, ATOM_H);

    addBond(mol, 0, 1, 1);
    addBond(mol, 0, 2, 1);
    addBond(mol, 2, 3, 2);
    addBond(mol, 2, 4, 1);
    addBond(mol, 0, 5, 1);
    addBond(mol, 5, 6, 1);
    addBond(mol, 6, 7, 1);
    addBond(mol, 6, 8, 2);
    addBond(mol, 8, 9, 1);
    addBond(mol, 9, 10, 1);
    addBond(mol, 10, 7, 1);
    addBond(mol, 7, 11, 2);
    addBond(mol, 11, 12, 1);
    addBond(mol, 12, 13, 2);
    addBond(mol, 13, 14, 1);
    addBond(mol, 14, 10, 2);
    addBond(mol, 0, 15, 1);
    addBond(mol, 1, 16, 1);
    addBond(mol, 1, 17, 1);
    addBond(mol, 4, 18, 1);
    addBond(mol, 8, 19, 1);
    addBond(mol, 9, 20, 1);

    centerMolecule(mol);
}

// Build Proline (C5H9NO2)
void buildProline(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Proline (C5H9NO2)");

    // 5-membered ring with N
    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_C);      // Alpha C
    addAtom(mol, -1.2f, 0.0f, 0.0f, ATOM_N);     // N (in ring)
    addAtom(mol, 1.3f, 0.7f, 0.0f, ATOM_C);      // COOH C
    addAtom(mol, 1.3f, 2.0f, 0.0f, ATOM_O);
    addAtom(mol, 2.4f, 0.0f, 0.0f, ATOM_O);
    addAtom(mol, 0.0f, -1.5f, 0.0f, ATOM_C);     // Ring C
    addAtom(mol, -1.4f, -1.5f, 0.0f, ATOM_C);    // Ring C
    addAtom(mol, -1.8f, -0.05f, 0.0f, ATOM_H);   // H on N
    // H atoms
    addAtom(mol, 0.0f, 0.6f, 0.9f, ATOM_H);
    addAtom(mol, 3.2f, 0.5f, 0.0f, ATOM_H);
    addAtom(mol, 0.5f, -2.0f, 0.9f, ATOM_H);
    addAtom(mol, 0.5f, -2.0f, -0.9f, ATOM_H);
    addAtom(mol, -1.9f, -2.0f, 0.9f, ATOM_H);
    addAtom(mol, -1.9f, -2.0f, -0.9f, ATOM_H);

    addBond(mol, 0, 1, 1);
    addBond(mol, 0, 2, 1);
    addBond(mol, 2, 3, 2);
    addBond(mol, 2, 4, 1);
    addBond(mol, 0, 5, 1);
    addBond(mol, 5, 6, 1);
    addBond(mol, 6, 1, 1);
    addBond(mol, 1, 7, 1);
    addBond(mol, 0, 8, 1);
    addBond(mol, 4, 9, 1);
    addBond(mol, 5, 10, 1);
    addBond(mol, 5, 11, 1);
    addBond(mol, 6, 12, 1);
    addBond(mol, 6, 13, 1);

    centerMolecule(mol);
}

// Build Cysteine (C3H7NO2S)
void buildCysteine(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Cysteine (C3H7NO2S)");

    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_C);      // Alpha C
    addAtom(mol, -1.3f, 0.5f, 0.0f, ATOM_N);     // NH2
    addAtom(mol, 1.3f, 0.7f, 0.0f, ATOM_C);      // COOH C
    addAtom(mol, 1.3f, 2.0f, 0.0f, ATOM_O);
    addAtom(mol, 2.4f, 0.0f, 0.0f, ATOM_O);
    addAtom(mol, 0.0f, -1.5f, 0.0f, ATOM_C);     // CH2
    addAtom(mol, 0.0f, -3.2f, 0.0f, ATOM_S);     // SH
    addAtom(mol, 0.0f, 0.5f, 0.9f, ATOM_H);
    addAtom(mol, -1.4f, 1.5f, 0.0f, ATOM_H);
    addAtom(mol, -2.1f, 0.0f, 0.0f, ATOM_H);
    addAtom(mol, 3.2f, 0.5f, 0.0f, ATOM_H);
    addAtom(mol, -0.9f, -1.8f, 0.0f, ATOM_H);
    addAtom(mol, 0.9f, -1.8f, 0.0f, ATOM_H);
    addAtom(mol, 1.0f, -3.8f, 0.0f, ATOM_H);

    addBond(mol, 0, 1, 1);
    addBond(mol, 0, 2, 1);
    addBond(mol, 2, 3, 2);
    addBond(mol, 2, 4, 1);
    addBond(mol, 0, 5, 1);
    addBond(mol, 5, 6, 1);
    addBond(mol, 0, 7, 1);
    addBond(mol, 1, 8, 1);
    addBond(mol, 1, 9, 1);
    addBond(mol, 4, 10, 1);
    addBond(mol, 5, 11, 1);
    addBond(mol, 5, 12, 1);
    addBond(mol, 6, 13, 1);

    centerMolecule(mol);
}

// Build Methionine (C5H11NO2S)
void buildMethionine(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Methionine (C5H11NO2S)");

    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_C);      // Alpha C
    addAtom(mol, -1.3f, 0.5f, 0.0f, ATOM_N);     // NH2
    addAtom(mol, 1.3f, 0.7f, 0.0f, ATOM_C);      // COOH C
    addAtom(mol, 1.3f, 2.0f, 0.0f, ATOM_O);
    addAtom(mol, 2.4f, 0.0f, 0.0f, ATOM_O);
    addAtom(mol, 0.0f, -1.5f, 0.0f, ATOM_C);     // CH2
    addAtom(mol, 0.0f, -3.0f, 0.0f, ATOM_C);     // CH2
    addAtom(mol, 0.0f, -4.7f, 0.0f, ATOM_S);     // S
    addAtom(mol, 0.0f, -6.3f, 0.0f, ATOM_C);     // CH3
    addAtom(mol, 0.0f, 0.5f, 0.9f, ATOM_H);
    addAtom(mol, -1.4f, 1.5f, 0.0f, ATOM_H);
    addAtom(mol, -2.1f, 0.0f, 0.0f, ATOM_H);
    addAtom(mol, 3.2f, 0.5f, 0.0f, ATOM_H);
    addAtom(mol, 0.5f, -6.8f, 0.87f, ATOM_H);
    addAtom(mol, 0.5f, -6.8f, -0.87f, ATOM_H);
    addAtom(mol, -1.0f, -6.6f, 0.0f, ATOM_H);

    addBond(mol, 0, 1, 1);
    addBond(mol, 0, 2, 1);
    addBond(mol, 2, 3, 2);
    addBond(mol, 2, 4, 1);
    addBond(mol, 0, 5, 1);
    addBond(mol, 5, 6, 1);
    addBond(mol, 6, 7, 1);
    addBond(mol, 7, 8, 1);
    addBond(mol, 0, 9, 1);
    addBond(mol, 1, 10, 1);
    addBond(mol, 1, 11, 1);
    addBond(mol, 4, 12, 1);
    addBond(mol, 8, 13, 1);
    addBond(mol, 8, 14, 1);
    addBond(mol, 8, 15, 1);

    centerMolecule(mol);
}

// Build Pyruvate / Pyruvic Acid (C3H4O3)
void buildPyruvate(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Pyruvate (C3H4O3)");

    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_C);      // CH3
    addAtom(mol, 1.5f, 0.0f, 0.0f, ATOM_C);      // C=O (ketone)
    addAtom(mol, 1.5f, 1.2f, 0.0f, ATOM_O);      // =O
    addAtom(mol, 3.0f, 0.0f, 0.0f, ATOM_C);      // COOH
    addAtom(mol, 3.6f, 1.1f, 0.0f, ATOM_O);      // =O
    addAtom(mol, 3.6f, -1.1f, 0.0f, ATOM_O);     // OH
    addAtom(mol, 4.5f, -1.0f, 0.0f, ATOM_H);
    addAtom(mol, -0.5f, 1.0f, 0.0f, ATOM_H);
    addAtom(mol, -0.5f, -0.5f, 0.87f, ATOM_H);
    addAtom(mol, -0.5f, -0.5f, -0.87f, ATOM_H);

    addBond(mol, 0, 1, 1);
    addBond(mol, 1, 2, 2);
    addBond(mol, 1, 3, 1);
    addBond(mol, 3, 4, 2);
    addBond(mol, 3, 5, 1);
    addBond(mol, 5, 6, 1);
    addBond(mol, 0, 7, 1);
    addBond(mol, 0, 8, 1);
    addBond(mol, 0, 9, 1);

    centerMolecule(mol);
}

// Build Arginine (C6H14N4O2)
void buildArginine(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Arginine (C6H14N4O2)");

    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_C);      // Alpha C
    addAtom(mol, -1.3f, 0.5f, 0.0f, ATOM_N);     // NH2
    addAtom(mol, 1.3f, 0.7f, 0.0f, ATOM_C);      // COOH C
    addAtom(mol, 1.3f, 2.0f, 0.0f, ATOM_O);
    addAtom(mol, 2.4f, 0.0f, 0.0f, ATOM_O);
    addAtom(mol, 0.0f, -1.5f, 0.0f, ATOM_C);     // CH2
    addAtom(mol, 0.0f, -3.0f, 0.0f, ATOM_C);     // CH2
    addAtom(mol, 0.0f, -4.5f, 0.0f, ATOM_C);     // CH2
    addAtom(mol, 0.0f, -6.0f, 0.0f, ATOM_N);     // NH
    addAtom(mol, 0.0f, -7.3f, 0.0f, ATOM_C);     // Guanidinium C
    addAtom(mol, -1.1f, -8.0f, 0.0f, ATOM_N);    // NH2
    addAtom(mol, 1.1f, -8.0f, 0.0f, ATOM_N);     // NH2
    // H atoms
    addAtom(mol, 0.0f, 0.5f, 0.9f, ATOM_H);
    addAtom(mol, -1.4f, 1.5f, 0.0f, ATOM_H);
    addAtom(mol, -2.1f, 0.0f, 0.0f, ATOM_H);
    addAtom(mol, 3.2f, 0.5f, 0.0f, ATOM_H);
    addAtom(mol, 0.8f, -6.3f, 0.0f, ATOM_H);
    addAtom(mol, -1.1f, -9.0f, 0.0f, ATOM_H);
    addAtom(mol, -1.9f, -7.5f, 0.0f, ATOM_H);
    addAtom(mol, 1.1f, -9.0f, 0.0f, ATOM_H);
    addAtom(mol, 1.9f, -7.5f, 0.0f, ATOM_H);

    addBond(mol, 0, 1, 1);
    addBond(mol, 0, 2, 1);
    addBond(mol, 2, 3, 2);
    addBond(mol, 2, 4, 1);
    addBond(mol, 0, 5, 1);
    addBond(mol, 5, 6, 1);
    addBond(mol, 6, 7, 1);
    addBond(mol, 7, 8, 1);
    addBond(mol, 8, 9, 1);
    addBond(mol, 9, 10, 1);
    addBond(mol, 9, 11, 2);
    addBond(mol, 0, 12, 1);
    addBond(mol, 1, 13, 1);
    addBond(mol, 1, 14, 1);
    addBond(mol, 4, 15, 1);
    addBond(mol, 8, 16, 1);
    addBond(mol, 10, 17, 1);
    addBond(mol, 10, 18, 1);
    addBond(mol, 11, 19, 1);
    addBond(mol, 11, 20, 1);

    centerMolecule(mol);
}

// Build Asparagine (C4H8N2O3) - 17 atoms
void buildAsparagine(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Asparagine (C4H8N2O3)");

    // Backbone
    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_C);      // C0: Alpha carbon
    addAtom(mol, -1.3f, 0.5f, 0.0f, ATOM_N);     // N1: Amino group
    addAtom(mol, 1.3f, 0.7f, 0.0f, ATOM_C);      // C2: Carboxyl carbon
    addAtom(mol, 1.3f, 2.0f, 0.0f, ATOM_O);      // O3: Carboxyl =O
    addAtom(mol, 2.4f, 0.0f, 0.0f, ATOM_O);      // O4: Carboxyl -OH

    // Side chain: -CH2-C(=O)-NH2
    addAtom(mol, 0.0f, -1.5f, 0.0f, ATOM_C);     // C5: Beta carbon (CH2)
    addAtom(mol, 0.0f, -3.0f, 0.0f, ATOM_C);     // C6: Amide carbon
    addAtom(mol, 1.1f, -3.7f, 0.0f, ATOM_O);     // O7: Amide =O
    addAtom(mol, -1.1f, -3.7f, 0.0f, ATOM_N);    // N8: Amide -NH2

    // Hydrogens
    addAtom(mol, 0.0f, 0.5f, 0.9f, ATOM_H);      // H9: on alpha C
    addAtom(mol, -1.4f, 1.5f, 0.0f, ATOM_H);     // H10: on NH2
    addAtom(mol, -2.1f, 0.0f, 0.0f, ATOM_H);     // H11: on NH2
    addAtom(mol, 3.2f, 0.5f, 0.0f, ATOM_H);      // H12: on COOH
    addAtom(mol, 0.9f, -1.5f, 0.5f, ATOM_H);     // H13: on CH2
    addAtom(mol, -0.9f, -1.5f, 0.5f, ATOM_H);    // H14: on CH2
    addAtom(mol, -1.1f, -4.7f, 0.0f, ATOM_H);    // H15: on amide NH2
    addAtom(mol, -1.9f, -3.2f, 0.0f, ATOM_H);    // H16: on amide NH2

    // Bonds
    addBond(mol, 0, 1, 1);   // Alpha C - NH2
    addBond(mol, 0, 2, 1);   // Alpha C - COOH
    addBond(mol, 2, 3, 2);   // C=O (carboxyl)
    addBond(mol, 2, 4, 1);   // C-OH (carboxyl)
    addBond(mol, 0, 5, 1);   // Alpha C - CH2
    addBond(mol, 5, 6, 1);   // CH2 - amide C
    addBond(mol, 6, 7, 2);   // C=O (amide)
    addBond(mol, 6, 8, 1);   // C-NH2 (amide)
    addBond(mol, 0, 9, 1);   // Alpha C - H
    addBond(mol, 1, 10, 1);  // NH2 - H
    addBond(mol, 1, 11, 1);  // NH2 - H
    addBond(mol, 4, 12, 1);  // OH - H
    addBond(mol, 5, 13, 1);  // CH2 - H
    addBond(mol, 5, 14, 1);  // CH2 - H
    addBond(mol, 8, 15, 1);  // Amide NH2 - H
    addBond(mol, 8, 16, 1);  // Amide NH2 - H

    centerMolecule(mol);
}

// Build Glutamine (C5H10N2O3) - 20 atoms
void buildGlutamine(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Glutamine (C5H10N2O3)");

    // Backbone
    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_C);      // C0: Alpha carbon
    addAtom(mol, -1.3f, 0.5f, 0.0f, ATOM_N);     // N1: Amino group
    addAtom(mol, 1.3f, 0.7f, 0.0f, ATOM_C);      // C2: Carboxyl carbon
    addAtom(mol, 1.3f, 2.0f, 0.0f, ATOM_O);      // O3: Carboxyl =O
    addAtom(mol, 2.4f, 0.0f, 0.0f, ATOM_O);      // O4: Carboxyl -OH

    // Side chain: -CH2-CH2-C(=O)-NH2
    addAtom(mol, 0.0f, -1.5f, 0.0f, ATOM_C);     // C5: Beta carbon (CH2)
    addAtom(mol, 0.0f, -3.0f, 0.0f, ATOM_C);     // C6: Gamma carbon (CH2)
    addAtom(mol, 0.0f, -4.5f, 0.0f, ATOM_C);     // C7: Amide carbon
    addAtom(mol, 1.1f, -5.2f, 0.0f, ATOM_O);     // O8: Amide =O
    addAtom(mol, -1.1f, -5.2f, 0.0f, ATOM_N);    // N9: Amide -NH2

    // Hydrogens
    addAtom(mol, 0.0f, 0.5f, 0.9f, ATOM_H);      // H10: on alpha C
    addAtom(mol, -1.4f, 1.5f, 0.0f, ATOM_H);     // H11: on NH2
    addAtom(mol, -2.1f, 0.0f, 0.0f, ATOM_H);     // H12: on NH2
    addAtom(mol, 3.2f, 0.5f, 0.0f, ATOM_H);      // H13: on COOH
    addAtom(mol, 0.9f, -1.5f, 0.5f, ATOM_H);     // H14: on beta CH2
    addAtom(mol, -0.9f, -1.5f, 0.5f, ATOM_H);    // H15: on beta CH2
    addAtom(mol, 0.9f, -3.0f, 0.5f, ATOM_H);     // H16: on gamma CH2
    addAtom(mol, -0.9f, -3.0f, 0.5f, ATOM_H);    // H17: on gamma CH2
    addAtom(mol, -1.1f, -6.2f, 0.0f, ATOM_H);    // H18: on amide NH2
    addAtom(mol, -1.9f, -4.7f, 0.0f, ATOM_H);    // H19: on amide NH2

    // Bonds
    addBond(mol, 0, 1, 1);   // Alpha C - NH2
    addBond(mol, 0, 2, 1);   // Alpha C - COOH
    addBond(mol, 2, 3, 2);   // C=O (carboxyl)
    addBond(mol, 2, 4, 1);   // C-OH (carboxyl)
    addBond(mol, 0, 5, 1);   // Alpha C - beta CH2
    addBond(mol, 5, 6, 1);   // Beta CH2 - gamma CH2
    addBond(mol, 6, 7, 1);   // Gamma CH2 - amide C
    addBond(mol, 7, 8, 2);   // C=O (amide)
    addBond(mol, 7, 9, 1);   // C-NH2 (amide)
    addBond(mol, 0, 10, 1);  // Alpha C - H
    addBond(mol, 1, 11, 1);  // NH2 - H
    addBond(mol, 1, 12, 1);  // NH2 - H
    addBond(mol, 4, 13, 1);  // OH - H
    addBond(mol, 5, 14, 1);  // Beta CH2 - H
    addBond(mol, 5, 15, 1);  // Beta CH2 - H
    addBond(mol, 6, 16, 1);  // Gamma CH2 - H
    addBond(mol, 6, 17, 1);  // Gamma CH2 - H
    addBond(mol, 9, 18, 1);  // Amide NH2 - H
    addBond(mol, 9, 19, 1);  // Amide NH2 - H

    centerMolecule(mol);
}

// ============== VITAMINS ==============

// Build Vitamin C - L-Ascorbic Acid (C6H8O6)
void buildAscorbicAcid(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Vitamin C (C6H8O6)");

    // Lactone ring (5-membered)
    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_C);      // C1
    addAtom(mol, 1.2f, 0.7f, 0.0f, ATOM_C);      // C2
    addAtom(mol, 2.0f, -0.3f, 0.0f, ATOM_C);     // C3
    addAtom(mol, 1.3f, -1.4f, 0.0f, ATOM_O);     // O (ring)
    addAtom(mol, 0.0f, -1.2f, 0.0f, ATOM_C);     // C4 (lactone carbonyl)
    addAtom(mol, -1.0f, -1.8f, 0.0f, ATOM_O);    // =O (lactone)
    // Side chain
    addAtom(mol, 3.4f, -0.3f, 0.0f, ATOM_C);     // C5
    addAtom(mol, 4.2f, 0.9f, 0.0f, ATOM_C);      // C6
    // Hydroxyl groups
    addAtom(mol, -0.8f, 0.8f, 0.0f, ATOM_O);     // OH on C1
    addAtom(mol, 1.2f, 2.0f, 0.0f, ATOM_O);      // OH on C2
    addAtom(mol, 3.8f, -1.5f, 0.0f, ATOM_O);     // OH on C5
    addAtom(mol, 5.5f, 0.7f, 0.0f, ATOM_O);      // OH on C6
    // H atoms on OH
    addAtom(mol, -1.6f, 0.4f, 0.0f, ATOM_H);
    addAtom(mol, 1.2f, 2.7f, 0.0f, ATOM_H);
    addAtom(mol, 4.7f, -1.6f, 0.0f, ATOM_H);
    addAtom(mol, 5.8f, 1.5f, 0.0f, ATOM_H);
    // Other H atoms
    addAtom(mol, 3.6f, 0.2f, 0.9f, ATOM_H);
    addAtom(mol, 3.8f, 1.8f, 0.0f, ATOM_H);
    addAtom(mol, 4.0f, 0.9f, -0.9f, ATOM_H);

    addBond(mol, 0, 1, 2);
    addBond(mol, 1, 2, 1);
    addBond(mol, 2, 3, 1);
    addBond(mol, 3, 4, 1);
    addBond(mol, 4, 0, 1);
    addBond(mol, 4, 5, 2);
    addBond(mol, 2, 6, 1);
    addBond(mol, 6, 7, 1);
    addBond(mol, 0, 8, 1);
    addBond(mol, 1, 9, 1);
    addBond(mol, 6, 10, 1);
    addBond(mol, 7, 11, 1);
    addBond(mol, 8, 12, 1);
    addBond(mol, 9, 13, 1);
    addBond(mol, 10, 14, 1);
    addBond(mol, 11, 15, 1);
    addBond(mol, 6, 16, 1);
    addBond(mol, 7, 17, 1);
    addBond(mol, 7, 18, 1);

    centerMolecule(mol);
}

// Build Thiamine - Vitamin B1 (C12H17N4OS)
void buildThiamine(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Vitamin B1 (C12H17N4OS)");

    // Pyrimidine ring
    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_N);      // N1
    addAtom(mol, 1.2f, 0.7f, 0.0f, ATOM_C);      // C2
    addAtom(mol, 2.4f, 0.0f, 0.0f, ATOM_N);      // N3
    addAtom(mol, 2.4f, -1.4f, 0.0f, ATOM_C);     // C4
    addAtom(mol, 1.2f, -2.1f, 0.0f, ATOM_C);     // C5
    addAtom(mol, 0.0f, -1.4f, 0.0f, ATOM_C);     // C6
    // Amino group on C4
    addAtom(mol, 3.6f, -2.0f, 0.0f, ATOM_N);     // NH2
    // Methyl on C2
    addAtom(mol, 1.2f, 2.1f, 0.0f, ATOM_C);      // CH3
    // Methylene bridge
    addAtom(mol, 1.2f, -3.6f, 0.0f, ATOM_C);     // CH2
    // Thiazole ring
    addAtom(mol, 1.2f, -5.0f, 0.0f, ATOM_C);     // C
    addAtom(mol, 0.0f, -5.7f, 0.0f, ATOM_N);     // N (positive)
    addAtom(mol, 0.0f, -7.1f, 0.0f, ATOM_C);     // C
    addAtom(mol, 1.2f, -7.8f, 0.0f, ATOM_C);     // C
    addAtom(mol, 2.2f, -6.8f, 0.0f, ATOM_S);     // S
    // Hydroxyethyl side chain
    addAtom(mol, 1.2f, -9.3f, 0.0f, ATOM_C);     // CH2
    addAtom(mol, 1.2f, -10.7f, 0.0f, ATOM_C);    // CH2
    addAtom(mol, 1.2f, -12.0f, 0.0f, ATOM_O);    // OH
    // Methyl on thiazole
    addAtom(mol, 2.4f, -4.3f, 0.0f, ATOM_C);     // CH3
    // H atoms
    addAtom(mol, -0.9f, -1.8f, 0.0f, ATOM_H);
    addAtom(mol, 3.6f, -3.0f, 0.0f, ATOM_H);
    addAtom(mol, 4.4f, -1.5f, 0.0f, ATOM_H);
    addAtom(mol, -0.9f, -7.5f, 0.0f, ATOM_H);
    addAtom(mol, 1.2f, -12.8f, 0.0f, ATOM_H);

    addBond(mol, 0, 1, 2);
    addBond(mol, 1, 2, 1);
    addBond(mol, 2, 3, 2);
    addBond(mol, 3, 4, 1);
    addBond(mol, 4, 5, 2);
    addBond(mol, 5, 0, 1);
    addBond(mol, 3, 6, 1);
    addBond(mol, 1, 7, 1);
    addBond(mol, 4, 8, 1);
    addBond(mol, 8, 9, 1);
    addBond(mol, 9, 10, 2);
    addBond(mol, 10, 11, 1);
    addBond(mol, 11, 12, 2);
    addBond(mol, 12, 13, 1);
    addBond(mol, 13, 9, 1);
    addBond(mol, 12, 14, 1);
    addBond(mol, 14, 15, 1);
    addBond(mol, 15, 16, 1);
    addBond(mol, 9, 17, 1);
    addBond(mol, 5, 18, 1);
    addBond(mol, 6, 19, 1);
    addBond(mol, 6, 20, 1);
    addBond(mol, 11, 21, 1);
    addBond(mol, 16, 22, 1);

    centerMolecule(mol);
}

// Build Riboflavin - Vitamin B2 (C17H20N4O6)
void buildRiboflavin(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Vitamin B2 (C17H20N4O6)");

    // Isoalloxazine ring system (flavin)
    // Benzene ring (fused)
    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_C);      // C5a
    addAtom(mol, 1.2f, 0.7f, 0.0f, ATOM_C);      // C6
    addAtom(mol, 2.4f, 0.0f, 0.0f, ATOM_C);      // C7
    addAtom(mol, 2.4f, -1.4f, 0.0f, ATOM_C);     // C8
    addAtom(mol, 1.2f, -2.1f, 0.0f, ATOM_C);     // C9
    addAtom(mol, 0.0f, -1.4f, 0.0f, ATOM_C);     // C9a
    // Pyrazine ring
    addAtom(mol, -1.2f, 0.7f, 0.0f, ATOM_N);     // N5
    addAtom(mol, -1.2f, -2.1f, 0.0f, ATOM_N);    // N10
    // Pyrimidine ring
    addAtom(mol, -2.4f, 0.0f, 0.0f, ATOM_C);     // C4a
    addAtom(mol, -2.4f, -1.4f, 0.0f, ATOM_C);    // C10a
    addAtom(mol, -3.6f, 0.7f, 0.0f, ATOM_C);     // C4 (=O)
    addAtom(mol, -3.6f, 2.0f, 0.0f, ATOM_O);     // =O
    addAtom(mol, -4.8f, 0.0f, 0.0f, ATOM_N);     // N3
    addAtom(mol, -4.8f, -1.4f, 0.0f, ATOM_C);    // C2 (=O)
    addAtom(mol, -4.8f, -2.6f, 0.0f, ATOM_O);    // =O
    addAtom(mol, -3.6f, -2.1f, 0.0f, ATOM_N);    // N1
    // Methyl groups
    addAtom(mol, 1.2f, 2.1f, 0.0f, ATOM_C);      // CH3 on C6
    addAtom(mol, 2.4f, -2.8f, 0.0f, ATOM_C);     // CH3 on C8
    // Ribityl chain start
    addAtom(mol, -1.2f, -3.5f, 0.0f, ATOM_C);    // CH2
    addAtom(mol, -1.2f, -4.9f, 0.0f, ATOM_C);    // CHOH
    addAtom(mol, -1.2f, -6.3f, 0.0f, ATOM_C);    // CHOH
    addAtom(mol, -1.2f, -7.7f, 0.0f, ATOM_C);    // CHOH
    addAtom(mol, -1.2f, -9.1f, 0.0f, ATOM_C);    // CH2OH
    // Ribityl OH groups
    addAtom(mol, 0.0f, -5.2f, 0.0f, ATOM_O);
    addAtom(mol, -2.4f, -6.6f, 0.0f, ATOM_O);
    addAtom(mol, 0.0f, -8.0f, 0.0f, ATOM_O);
    addAtom(mol, -1.2f, -10.4f, 0.0f, ATOM_O);
    // Some H atoms
    addAtom(mol, 3.3f, 0.5f, 0.0f, ATOM_H);
    addAtom(mol, -5.6f, 0.5f, 0.0f, ATOM_H);

    addBond(mol, 0, 1, 2);
    addBond(mol, 1, 2, 1);
    addBond(mol, 2, 3, 2);
    addBond(mol, 3, 4, 1);
    addBond(mol, 4, 5, 2);
    addBond(mol, 5, 0, 1);
    addBond(mol, 0, 6, 1);
    addBond(mol, 5, 7, 1);
    addBond(mol, 6, 8, 2);
    addBond(mol, 7, 9, 1);
    addBond(mol, 8, 9, 1);
    addBond(mol, 8, 10, 1);
    addBond(mol, 10, 11, 2);
    addBond(mol, 10, 12, 1);
    addBond(mol, 12, 13, 1);
    addBond(mol, 13, 14, 2);
    addBond(mol, 13, 15, 1);
    addBond(mol, 15, 9, 2);
    addBond(mol, 1, 16, 1);
    addBond(mol, 3, 17, 1);
    addBond(mol, 7, 18, 1);
    addBond(mol, 18, 19, 1);
    addBond(mol, 19, 20, 1);
    addBond(mol, 20, 21, 1);
    addBond(mol, 21, 22, 1);
    addBond(mol, 19, 23, 1);
    addBond(mol, 20, 24, 1);
    addBond(mol, 21, 25, 1);
    addBond(mol, 22, 26, 1);
    addBond(mol, 2, 27, 1);
    addBond(mol, 12, 28, 1);

    centerMolecule(mol);
}

// Build Niacin - Vitamin B3 / Nicotinic Acid (C6H5NO2)
void buildNiacin(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Vitamin B3 (C6H5NO2)");

    // Pyridine ring
    float r = 1.4f;
    addAtom(mol, r * cosf(0), r * sinf(0), 0.0f, ATOM_N);        // N1
    addAtom(mol, r * cosf(PI/3), r * sinf(PI/3), 0.0f, ATOM_C);  // C2
    addAtom(mol, r * cosf(2*PI/3), r * sinf(2*PI/3), 0.0f, ATOM_C); // C3
    addAtom(mol, r * cosf(PI), r * sinf(PI), 0.0f, ATOM_C);      // C4
    addAtom(mol, r * cosf(4*PI/3), r * sinf(4*PI/3), 0.0f, ATOM_C); // C5
    addAtom(mol, r * cosf(5*PI/3), r * sinf(5*PI/3), 0.0f, ATOM_C); // C6
    // Carboxylic acid on C3
    addAtom(mol, r * cosf(2*PI/3) * 1.8f, r * sinf(2*PI/3) * 1.8f + 0.6f, 0.0f, ATOM_C); // COOH C
    addAtom(mol, r * cosf(2*PI/3) * 2.2f, r * sinf(2*PI/3) * 2.2f + 1.5f, 0.0f, ATOM_O); // =O
    addAtom(mol, r * cosf(2*PI/3) * 2.4f, r * sinf(2*PI/3) * 1.2f, 0.0f, ATOM_O);        // OH
    addAtom(mol, r * cosf(2*PI/3) * 3.2f, r * sinf(2*PI/3) * 1.0f, 0.0f, ATOM_H);        // H on OH
    // H atoms on ring
    float rH = 2.4f;
    addAtom(mol, rH * cosf(PI/3), rH * sinf(PI/3), 0.0f, ATOM_H);
    addAtom(mol, rH * cosf(PI), rH * sinf(PI), 0.0f, ATOM_H);
    addAtom(mol, rH * cosf(4*PI/3), rH * sinf(4*PI/3), 0.0f, ATOM_H);
    addAtom(mol, rH * cosf(5*PI/3), rH * sinf(5*PI/3), 0.0f, ATOM_H);

    addBond(mol, 0, 1, 2);
    addBond(mol, 1, 2, 1);
    addBond(mol, 2, 3, 2);
    addBond(mol, 3, 4, 1);
    addBond(mol, 4, 5, 2);
    addBond(mol, 5, 0, 1);
    addBond(mol, 2, 6, 1);
    addBond(mol, 6, 7, 2);
    addBond(mol, 6, 8, 1);
    addBond(mol, 8, 9, 1);
    addBond(mol, 1, 10, 1);
    addBond(mol, 3, 11, 1);
    addBond(mol, 4, 12, 1);
    addBond(mol, 5, 13, 1);

    centerMolecule(mol);
}

// Build Pantothenic Acid - Vitamin B5 (C9H17NO5)
void buildPanthothenicAcid(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Vitamin B5 (C9H17NO5)");

    // Beta-alanine part
    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_C);      // COOH C
    addAtom(mol, 0.0f, 1.2f, 0.0f, ATOM_O);      // =O
    addAtom(mol, 1.1f, -0.6f, 0.0f, ATOM_O);     // OH
    addAtom(mol, -1.3f, -0.7f, 0.0f, ATOM_C);    // CH2
    addAtom(mol, -2.6f, 0.0f, 0.0f, ATOM_C);     // CH2
    addAtom(mol, -3.9f, -0.7f, 0.0f, ATOM_N);    // NH
    // Amide linkage
    addAtom(mol, -5.2f, 0.0f, 0.0f, ATOM_C);     // C=O
    addAtom(mol, -5.2f, 1.2f, 0.0f, ATOM_O);     // =O
    // Pantoic acid part
    addAtom(mol, -6.5f, -0.7f, 0.0f, ATOM_C);    // C (with OH)
    addAtom(mol, -6.5f, -2.1f, 0.0f, ATOM_O);    // OH
    addAtom(mol, -7.8f, 0.0f, 0.0f, ATOM_C);     // C (quaternary)
    addAtom(mol, -7.8f, 1.4f, 0.0f, ATOM_C);     // CH3
    addAtom(mol, -7.8f, -1.4f, 0.0f, ATOM_C);    // CH3
    addAtom(mol, -9.1f, 0.0f, 0.0f, ATOM_C);     // CH2
    addAtom(mol, -10.4f, 0.0f, 0.0f, ATOM_O);    // OH
    // H atoms
    addAtom(mol, 1.8f, -0.1f, 0.0f, ATOM_H);
    addAtom(mol, -3.9f, -1.6f, 0.0f, ATOM_H);
    addAtom(mol, -6.5f, -2.8f, 0.0f, ATOM_H);
    addAtom(mol, -11.1f, 0.0f, 0.0f, ATOM_H);

    addBond(mol, 0, 1, 2);
    addBond(mol, 0, 2, 1);
    addBond(mol, 0, 3, 1);
    addBond(mol, 3, 4, 1);
    addBond(mol, 4, 5, 1);
    addBond(mol, 5, 6, 1);
    addBond(mol, 6, 7, 2);
    addBond(mol, 6, 8, 1);
    addBond(mol, 8, 9, 1);
    addBond(mol, 8, 10, 1);
    addBond(mol, 10, 11, 1);
    addBond(mol, 10, 12, 1);
    addBond(mol, 10, 13, 1);
    addBond(mol, 13, 14, 1);
    addBond(mol, 2, 15, 1);
    addBond(mol, 5, 16, 1);
    addBond(mol, 9, 17, 1);
    addBond(mol, 14, 18, 1);

    centerMolecule(mol);
}

// Build Pyridoxine - Vitamin B6 (C8H11NO3)
void buildPyridoxine(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Vitamin B6 (C8H11NO3)");

    // Pyridine ring
    float r = 1.4f;
    addAtom(mol, r * cosf(0), r * sinf(0), 0.0f, ATOM_N);        // N1
    addAtom(mol, r * cosf(PI/3), r * sinf(PI/3), 0.0f, ATOM_C);  // C2
    addAtom(mol, r * cosf(2*PI/3), r * sinf(2*PI/3), 0.0f, ATOM_C); // C3
    addAtom(mol, r * cosf(PI), r * sinf(PI), 0.0f, ATOM_C);      // C4
    addAtom(mol, r * cosf(4*PI/3), r * sinf(4*PI/3), 0.0f, ATOM_C); // C5
    addAtom(mol, r * cosf(5*PI/3), r * sinf(5*PI/3), 0.0f, ATOM_C); // C6
    // CH3 on C2
    addAtom(mol, r * cosf(PI/3) * 1.8f, r * sinf(PI/3) * 1.8f, 0.0f, ATOM_C);
    // CH2OH on C4
    addAtom(mol, r * cosf(PI) * 1.8f, r * sinf(PI) * 1.8f, 0.0f, ATOM_C);
    addAtom(mol, r * cosf(PI) * 2.8f, r * sinf(PI) * 2.0f, 0.0f, ATOM_O);
    // OH on C3
    addAtom(mol, r * cosf(2*PI/3) * 1.8f, r * sinf(2*PI/3) * 1.8f, 0.0f, ATOM_O);
    // CH2OH on C5
    addAtom(mol, r * cosf(4*PI/3) * 1.8f, r * sinf(4*PI/3) * 1.8f, 0.0f, ATOM_C);
    addAtom(mol, r * cosf(4*PI/3) * 2.8f, r * sinf(4*PI/3) * 2.0f, 0.0f, ATOM_O);
    // H atoms
    addAtom(mol, r * cosf(5*PI/3) * 1.7f, r * sinf(5*PI/3) * 1.7f, 0.0f, ATOM_H);
    addAtom(mol, r * cosf(2*PI/3) * 2.5f, r * sinf(2*PI/3) * 2.5f, 0.0f, ATOM_H);
    addAtom(mol, r * cosf(PI) * 3.5f, r * sinf(PI) * 2.3f, 0.0f, ATOM_H);
    addAtom(mol, r * cosf(4*PI/3) * 3.5f, r * sinf(4*PI/3) * 2.3f, 0.0f, ATOM_H);

    addBond(mol, 0, 1, 2);
    addBond(mol, 1, 2, 1);
    addBond(mol, 2, 3, 2);
    addBond(mol, 3, 4, 1);
    addBond(mol, 4, 5, 2);
    addBond(mol, 5, 0, 1);
    addBond(mol, 1, 6, 1);
    addBond(mol, 3, 7, 1);
    addBond(mol, 7, 8, 1);
    addBond(mol, 2, 9, 1);
    addBond(mol, 4, 10, 1);
    addBond(mol, 10, 11, 1);
    addBond(mol, 5, 12, 1);
    addBond(mol, 9, 13, 1);
    addBond(mol, 8, 14, 1);
    addBond(mol, 11, 15, 1);

    centerMolecule(mol);
}

// Build Biotin - Vitamin B7 (C10H16N2O3S)
void buildBiotin(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Vitamin B7 (C10H16N2O3S)");

    // Ureido ring (imidazolidinone)
    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_C);      // C=O (urea)
    addAtom(mol, 0.0f, 1.2f, 0.0f, ATOM_O);      // =O
    addAtom(mol, -1.1f, -0.7f, 0.0f, ATOM_N);    // NH
    addAtom(mol, 1.1f, -0.7f, 0.0f, ATOM_N);     // NH
    addAtom(mol, -1.1f, -2.1f, 0.0f, ATOM_C);    // CH
    addAtom(mol, 1.1f, -2.1f, 0.0f, ATOM_C);     // CH
    // Tetrahydrothiophene ring
    addAtom(mol, 0.0f, -2.8f, 0.0f, ATOM_S);     // S
    addAtom(mol, -1.1f, -3.5f, 0.0f, ATOM_C);    // CH
    addAtom(mol, 1.1f, -3.5f, 0.0f, ATOM_C);     // CH2
    // Valeric acid side chain
    addAtom(mol, -1.1f, -5.0f, 0.0f, ATOM_C);    // CH2
    addAtom(mol, -1.1f, -6.5f, 0.0f, ATOM_C);    // CH2
    addAtom(mol, -1.1f, -8.0f, 0.0f, ATOM_C);    // CH2
    addAtom(mol, -1.1f, -9.5f, 0.0f, ATOM_C);    // CH2
    addAtom(mol, -1.1f, -11.0f, 0.0f, ATOM_C);   // COOH
    addAtom(mol, -1.1f, -12.2f, 0.0f, ATOM_O);   // =O
    addAtom(mol, 0.0f, -11.5f, 0.0f, ATOM_O);    // OH
    // H atoms
    addAtom(mol, -1.8f, -0.2f, 0.0f, ATOM_H);
    addAtom(mol, 1.8f, -0.2f, 0.0f, ATOM_H);
    addAtom(mol, 0.7f, -11.9f, 0.0f, ATOM_H);

    addBond(mol, 0, 1, 2);
    addBond(mol, 0, 2, 1);
    addBond(mol, 0, 3, 1);
    addBond(mol, 2, 4, 1);
    addBond(mol, 3, 5, 1);
    addBond(mol, 4, 6, 1);
    addBond(mol, 5, 6, 1);
    addBond(mol, 6, 7, 1);
    addBond(mol, 6, 8, 1);
    addBond(mol, 4, 7, 1);
    addBond(mol, 5, 8, 1);
    addBond(mol, 7, 9, 1);
    addBond(mol, 9, 10, 1);
    addBond(mol, 10, 11, 1);
    addBond(mol, 11, 12, 1);
    addBond(mol, 12, 13, 1);
    addBond(mol, 13, 14, 2);
    addBond(mol, 13, 15, 1);
    addBond(mol, 2, 16, 1);
    addBond(mol, 3, 17, 1);
    addBond(mol, 15, 18, 1);

    centerMolecule(mol);
}

// Build Folic Acid - Vitamin B9 (C19H19N7O6) - Simplified representation
void buildFolicAcid(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Vitamin B9 (C19H19N7O6)");

    // Pteridine ring system
    // Pyrimidine
    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_N);      // N1
    addAtom(mol, 1.2f, 0.7f, 0.0f, ATOM_C);      // C2
    addAtom(mol, 2.4f, 0.0f, 0.0f, ATOM_N);      // N3
    addAtom(mol, 2.4f, -1.4f, 0.0f, ATOM_C);     // C4
    addAtom(mol, 1.2f, -2.1f, 0.0f, ATOM_C);     // C4a
    addAtom(mol, 0.0f, -1.4f, 0.0f, ATOM_C);     // C8a
    // Pyrazine fused
    addAtom(mol, 1.2f, -3.5f, 0.0f, ATOM_N);     // N5
    addAtom(mol, 2.4f, -4.2f, 0.0f, ATOM_C);     // C6
    addAtom(mol, 3.6f, -3.5f, 0.0f, ATOM_C);     // C7
    addAtom(mol, 3.6f, -2.1f, 0.0f, ATOM_N);     // N8
    // Amino on C2
    addAtom(mol, 1.2f, 2.1f, 0.0f, ATOM_N);      // NH2
    // =O on C4
    addAtom(mol, 3.6f, -0.7f, 0.0f, ATOM_O);     // =O
    // p-Aminobenzoic acid part (simplified)
    addAtom(mol, 2.4f, -5.6f, 0.0f, ATOM_C);     // CH2
    addAtom(mol, 2.4f, -7.0f, 0.0f, ATOM_N);     // NH
    // Benzene ring
    addAtom(mol, 2.4f, -8.4f, 0.0f, ATOM_C);
    addAtom(mol, 1.2f, -9.1f, 0.0f, ATOM_C);
    addAtom(mol, 1.2f, -10.5f, 0.0f, ATOM_C);
    addAtom(mol, 2.4f, -11.2f, 0.0f, ATOM_C);
    addAtom(mol, 3.6f, -10.5f, 0.0f, ATOM_C);
    addAtom(mol, 3.6f, -9.1f, 0.0f, ATOM_C);
    // COOH on benzene
    addAtom(mol, 2.4f, -12.6f, 0.0f, ATOM_C);
    addAtom(mol, 1.2f, -13.3f, 0.0f, ATOM_O);
    addAtom(mol, 3.6f, -13.3f, 0.0f, ATOM_N);    // NH to glutamate
    // Glutamate (simplified)
    addAtom(mol, 4.8f, -12.6f, 0.0f, ATOM_C);    // CH
    addAtom(mol, 6.0f, -13.3f, 0.0f, ATOM_C);    // COOH
    addAtom(mol, 6.0f, -14.5f, 0.0f, ATOM_O);
    addAtom(mol, 7.1f, -12.6f, 0.0f, ATOM_O);

    addBond(mol, 0, 1, 2);
    addBond(mol, 1, 2, 1);
    addBond(mol, 2, 3, 2);
    addBond(mol, 3, 4, 1);
    addBond(mol, 4, 5, 2);
    addBond(mol, 5, 0, 1);
    addBond(mol, 4, 6, 1);
    addBond(mol, 6, 7, 2);
    addBond(mol, 7, 8, 1);
    addBond(mol, 8, 9, 2);
    addBond(mol, 9, 3, 1);
    addBond(mol, 1, 10, 1);
    addBond(mol, 3, 11, 1);
    addBond(mol, 7, 12, 1);
    addBond(mol, 12, 13, 1);
    addBond(mol, 13, 14, 1);
    addBond(mol, 14, 15, 2);
    addBond(mol, 15, 16, 1);
    addBond(mol, 16, 17, 2);
    addBond(mol, 17, 18, 1);
    addBond(mol, 18, 19, 2);
    addBond(mol, 19, 14, 1);
    addBond(mol, 17, 20, 1);
    addBond(mol, 20, 21, 2);
    addBond(mol, 20, 22, 1);
    addBond(mol, 22, 23, 1);
    addBond(mol, 23, 24, 1);
    addBond(mol, 24, 25, 2);
    addBond(mol, 24, 26, 1);

    centerMolecule(mol);
}

// Build Retinol - Vitamin A (C20H30O)
void buildRetinol(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Vitamin A (C20H30O)");

    // Cyclohexene ring (beta-ionone)
    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_C);      // C1
    addAtom(mol, 1.3f, 0.7f, 0.0f, ATOM_C);      // C2
    addAtom(mol, 2.6f, 0.0f, 0.0f, ATOM_C);      // C3
    addAtom(mol, 2.6f, -1.4f, 0.0f, ATOM_C);     // C4
    addAtom(mol, 1.3f, -2.1f, 0.0f, ATOM_C);     // C5 (C=C)
    addAtom(mol, 0.0f, -1.4f, 0.0f, ATOM_C);     // C6 (C=C)
    // Methyls on C1
    addAtom(mol, -1.1f, 0.7f, 0.0f, ATOM_C);     // CH3
    addAtom(mol, 0.0f, 1.0f, 1.0f, ATOM_C);      // CH3
    // Methyl on C5
    addAtom(mol, 1.3f, -3.5f, 0.0f, ATOM_C);     // CH3
    // Polyene chain
    addAtom(mol, -1.3f, -2.1f, 0.0f, ATOM_C);    // C7
    addAtom(mol, -2.6f, -1.4f, 0.0f, ATOM_C);    // C8
    addAtom(mol, -3.9f, -2.1f, 0.0f, ATOM_C);    // C9
    addAtom(mol, -3.9f, -3.5f, 0.0f, ATOM_C);    // CH3 on C9
    addAtom(mol, -5.2f, -1.4f, 0.0f, ATOM_C);    // C10
    addAtom(mol, -6.5f, -2.1f, 0.0f, ATOM_C);    // C11
    addAtom(mol, -7.8f, -1.4f, 0.0f, ATOM_C);    // C12
    addAtom(mol, -9.1f, -2.1f, 0.0f, ATOM_C);    // C13
    addAtom(mol, -9.1f, -3.5f, 0.0f, ATOM_C);    // CH3 on C13
    addAtom(mol, -10.4f, -1.4f, 0.0f, ATOM_C);   // C14
    addAtom(mol, -11.7f, -2.1f, 0.0f, ATOM_C);   // C15 (CH2OH)
    addAtom(mol, -13.0f, -1.4f, 0.0f, ATOM_O);   // OH
    addAtom(mol, -13.7f, -2.0f, 0.0f, ATOM_H);

    addBond(mol, 0, 1, 1);
    addBond(mol, 1, 2, 1);
    addBond(mol, 2, 3, 1);
    addBond(mol, 3, 4, 1);
    addBond(mol, 4, 5, 2);
    addBond(mol, 5, 0, 1);
    addBond(mol, 0, 6, 1);
    addBond(mol, 0, 7, 1);
    addBond(mol, 4, 8, 1);
    addBond(mol, 5, 9, 1);
    addBond(mol, 9, 10, 2);
    addBond(mol, 10, 11, 1);
    addBond(mol, 11, 12, 1);
    addBond(mol, 11, 13, 2);
    addBond(mol, 13, 14, 1);
    addBond(mol, 14, 15, 2);
    addBond(mol, 15, 16, 1);
    addBond(mol, 16, 17, 2);
    addBond(mol, 16, 18, 1);
    addBond(mol, 17, 19, 1);
    addBond(mol, 19, 20, 1);
    addBond(mol, 20, 21, 1);

    centerMolecule(mol);
}

// Build Beta-Carotene (C40H56) - Provitamin A - Simplified
void buildBetaCarotene(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Beta-Carotene (C40H56)");

    // Simplified: Two beta-ionone rings connected by polyene chain
    // Left ring
    addAtom(mol, -8.0f, 0.0f, 0.0f, ATOM_C);
    addAtom(mol, -7.0f, 1.2f, 0.0f, ATOM_C);
    addAtom(mol, -5.6f, 1.0f, 0.0f, ATOM_C);
    addAtom(mol, -5.0f, -0.3f, 0.0f, ATOM_C);
    addAtom(mol, -5.8f, -1.5f, 0.0f, ATOM_C);
    addAtom(mol, -7.2f, -1.3f, 0.0f, ATOM_C);
    // Central polyene (simplified)
    addAtom(mol, -3.6f, -0.5f, 0.0f, ATOM_C);
    addAtom(mol, -2.4f, 0.2f, 0.0f, ATOM_C);
    addAtom(mol, -1.0f, 0.0f, 0.0f, ATOM_C);
    addAtom(mol, 0.2f, 0.7f, 0.0f, ATOM_C);
    addAtom(mol, 1.6f, 0.5f, 0.0f, ATOM_C);
    addAtom(mol, 2.8f, 1.2f, 0.0f, ATOM_C);
    addAtom(mol, 4.2f, 1.0f, 0.0f, ATOM_C);
    // Right ring
    addAtom(mol, 5.4f, 0.3f, 0.0f, ATOM_C);
    addAtom(mol, 5.8f, -1.1f, 0.0f, ATOM_C);
    addAtom(mol, 7.2f, -1.3f, 0.0f, ATOM_C);
    addAtom(mol, 8.0f, 0.0f, 0.0f, ATOM_C);
    addAtom(mol, 7.0f, 1.2f, 0.0f, ATOM_C);
    addAtom(mol, 5.6f, 1.4f, 0.0f, ATOM_C);
    // Methyls
    addAtom(mol, -8.8f, 0.8f, 0.0f, ATOM_C);
    addAtom(mol, 8.8f, -0.8f, 0.0f, ATOM_C);
    addAtom(mol, -1.0f, -1.4f, 0.0f, ATOM_C);
    addAtom(mol, 1.6f, -0.9f, 0.0f, ATOM_C);

    addBond(mol, 0, 1, 1);
    addBond(mol, 1, 2, 1);
    addBond(mol, 2, 3, 1);
    addBond(mol, 3, 4, 2);
    addBond(mol, 4, 5, 1);
    addBond(mol, 5, 0, 1);
    addBond(mol, 3, 6, 1);
    addBond(mol, 6, 7, 2);
    addBond(mol, 7, 8, 1);
    addBond(mol, 8, 9, 2);
    addBond(mol, 9, 10, 1);
    addBond(mol, 10, 11, 2);
    addBond(mol, 11, 12, 1);
    addBond(mol, 12, 13, 2);
    addBond(mol, 13, 14, 1);
    addBond(mol, 14, 15, 1);
    addBond(mol, 15, 16, 1);
    addBond(mol, 16, 17, 1);
    addBond(mol, 17, 18, 1);
    addBond(mol, 18, 13, 1);
    addBond(mol, 0, 19, 1);
    addBond(mol, 16, 20, 1);
    addBond(mol, 8, 21, 1);
    addBond(mol, 10, 22, 1);

    centerMolecule(mol);
}

// Build Cholecalciferol - Vitamin D3 (C27H44O)
void buildCholecalciferol(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Vitamin D3 (C27H44O)");

    // Secosteroid structure - simplified
    // Ring A (opened)
    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_C);      // C1
    addAtom(mol, 1.3f, 0.7f, 0.0f, ATOM_C);      // C2
    addAtom(mol, 2.6f, 0.0f, 0.0f, ATOM_C);      // C3 (with OH)
    addAtom(mol, 2.6f, -1.4f, 0.0f, ATOM_O);     // OH
    addAtom(mol, 2.6f, -2.2f, 0.0f, ATOM_H);
    addAtom(mol, 3.9f, 0.7f, 0.0f, ATOM_C);      // C4
    addAtom(mol, 3.9f, 2.1f, 0.0f, ATOM_C);      // C5
    addAtom(mol, 2.6f, 2.8f, 0.0f, ATOM_C);      // C6
    addAtom(mol, 1.3f, 2.1f, 0.0f, ATOM_C);      // C7
    addAtom(mol, 0.0f, 2.8f, 0.0f, ATOM_C);      // C8
    // Ring C
    addAtom(mol, -1.3f, 2.1f, 0.0f, ATOM_C);     // C9
    addAtom(mol, -1.3f, 0.7f, 0.0f, ATOM_C);     // C10
    addAtom(mol, -2.6f, 2.8f, 0.0f, ATOM_C);     // C11
    addAtom(mol, -3.9f, 2.1f, 0.0f, ATOM_C);     // C12
    addAtom(mol, -3.9f, 0.7f, 0.0f, ATOM_C);     // C13
    addAtom(mol, -2.6f, 0.0f, 0.0f, ATOM_C);     // C14
    // Ring D
    addAtom(mol, -5.2f, 0.0f, 0.0f, ATOM_C);     // C15
    addAtom(mol, -6.5f, 0.7f, 0.0f, ATOM_C);     // C16
    addAtom(mol, -6.5f, 2.1f, 0.0f, ATOM_C);     // C17
    addAtom(mol, -5.2f, 2.8f, 0.0f, ATOM_C);     // C18 (CH3)
    // Side chain
    addAtom(mol, -7.8f, 2.8f, 0.0f, ATOM_C);     // C20
    addAtom(mol, -7.8f, 4.2f, 0.0f, ATOM_C);     // CH3
    addAtom(mol, -9.1f, 2.1f, 0.0f, ATOM_C);     // C22
    addAtom(mol, -10.4f, 2.8f, 0.0f, ATOM_C);    // C23
    addAtom(mol, -11.7f, 2.1f, 0.0f, ATOM_C);    // C24
    addAtom(mol, -13.0f, 2.8f, 0.0f, ATOM_C);    // C25
    addAtom(mol, -14.3f, 2.1f, 0.0f, ATOM_C);    // C26 (CH3)
    addAtom(mol, -13.0f, 4.2f, 0.0f, ATOM_C);    // C27 (CH3)

    addBond(mol, 0, 1, 1);
    addBond(mol, 1, 2, 1);
    addBond(mol, 2, 3, 1);
    addBond(mol, 3, 4, 1);
    addBond(mol, 2, 5, 1);
    addBond(mol, 5, 6, 2);
    addBond(mol, 6, 7, 1);
    addBond(mol, 7, 8, 2);
    addBond(mol, 8, 9, 1);
    addBond(mol, 9, 10, 1);
    addBond(mol, 10, 11, 1);
    addBond(mol, 11, 0, 1);
    addBond(mol, 9, 12, 1);
    addBond(mol, 12, 13, 1);
    addBond(mol, 13, 14, 1);
    addBond(mol, 14, 15, 1);
    addBond(mol, 15, 10, 1);
    addBond(mol, 14, 16, 1);
    addBond(mol, 16, 17, 1);
    addBond(mol, 17, 18, 1);
    addBond(mol, 18, 13, 1);
    addBond(mol, 13, 19, 1);
    addBond(mol, 18, 20, 1);
    addBond(mol, 20, 21, 1);
    addBond(mol, 20, 22, 1);
    addBond(mol, 22, 23, 1);
    addBond(mol, 23, 24, 1);
    addBond(mol, 24, 25, 1);
    addBond(mol, 25, 26, 1);
    addBond(mol, 25, 27, 1);

    centerMolecule(mol);
}

// Build Alpha-Tocopherol - Vitamin E (C29H50O2) - Simplified
void buildAlphaTocopherol(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Vitamin E (C29H50O2)");

    // Chromanol ring (benzene + pyran)
    float r = 1.4f;
    // Benzene part
    addAtom(mol, r, 0.0f, 0.0f, ATOM_C);
    addAtom(mol, r * 0.5f, r * 0.866f, 0.0f, ATOM_C);
    addAtom(mol, -r * 0.5f, r * 0.866f, 0.0f, ATOM_C);
    addAtom(mol, -r, 0.0f, 0.0f, ATOM_C);
    addAtom(mol, -r * 0.5f, -r * 0.866f, 0.0f, ATOM_C);
    addAtom(mol, r * 0.5f, -r * 0.866f, 0.0f, ATOM_C);
    // Pyran oxygen
    addAtom(mol, r * 1.8f, 0.0f, 0.0f, ATOM_O);
    // Pyran carbons
    addAtom(mol, r * 2.2f, r * 0.7f, 0.0f, ATOM_C);
    addAtom(mol, r * 1.5f, r * 1.5f, 0.0f, ATOM_C);
    // OH on benzene
    addAtom(mol, -r * 1.8f, 0.0f, 0.0f, ATOM_O);
    addAtom(mol, -r * 2.5f, 0.0f, 0.0f, ATOM_H);
    // Methyls on benzene
    addAtom(mol, r * 0.7f, r * 1.5f, 0.0f, ATOM_C);
    addAtom(mol, -r * 0.7f, r * 1.5f, 0.0f, ATOM_C);
    addAtom(mol, -r * 0.7f, -r * 1.5f, 0.0f, ATOM_C);
    // Phytyl tail (simplified)
    addAtom(mol, r * 3.5f, r * 0.7f, 0.0f, ATOM_C);
    addAtom(mol, r * 4.8f, r * 0.0f, 0.0f, ATOM_C);
    addAtom(mol, r * 6.1f, r * 0.7f, 0.0f, ATOM_C);
    addAtom(mol, r * 7.4f, r * 0.0f, 0.0f, ATOM_C);
    addAtom(mol, r * 8.7f, r * 0.7f, 0.0f, ATOM_C);
    addAtom(mol, r * 10.0f, r * 0.0f, 0.0f, ATOM_C);
    addAtom(mol, r * 3.5f, r * 2.0f, 0.0f, ATOM_C);  // CH3 branches
    addAtom(mol, r * 6.1f, r * 2.0f, 0.0f, ATOM_C);
    addAtom(mol, r * 8.7f, r * 2.0f, 0.0f, ATOM_C);

    addBond(mol, 0, 1, 2);
    addBond(mol, 1, 2, 1);
    addBond(mol, 2, 3, 2);
    addBond(mol, 3, 4, 1);
    addBond(mol, 4, 5, 2);
    addBond(mol, 5, 0, 1);
    addBond(mol, 0, 6, 1);
    addBond(mol, 6, 7, 1);
    addBond(mol, 7, 8, 1);
    addBond(mol, 8, 1, 1);
    addBond(mol, 3, 9, 1);
    addBond(mol, 9, 10, 1);
    addBond(mol, 1, 11, 1);
    addBond(mol, 2, 12, 1);
    addBond(mol, 4, 13, 1);
    addBond(mol, 7, 14, 1);
    addBond(mol, 14, 15, 1);
    addBond(mol, 15, 16, 1);
    addBond(mol, 16, 17, 1);
    addBond(mol, 17, 18, 1);
    addBond(mol, 18, 19, 1);
    addBond(mol, 14, 20, 1);
    addBond(mol, 16, 21, 1);
    addBond(mol, 18, 22, 1);

    centerMolecule(mol);
}

// Build Phylloquinone - Vitamin K1 (C31H46O2) - Simplified
void buildPhylloquinone(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Vitamin K1 (C31H46O2)");

    // Naphthoquinone core
    // First benzene ring
    float r = 1.4f;
    addAtom(mol, r, 0.0f, 0.0f, ATOM_C);
    addAtom(mol, r * 0.5f, r * 0.866f, 0.0f, ATOM_C);
    addAtom(mol, -r * 0.5f, r * 0.866f, 0.0f, ATOM_C);
    addAtom(mol, -r, 0.0f, 0.0f, ATOM_C);
    addAtom(mol, -r * 0.5f, -r * 0.866f, 0.0f, ATOM_C);
    addAtom(mol, r * 0.5f, -r * 0.866f, 0.0f, ATOM_C);
    // Quinone ring (fused)
    addAtom(mol, r * 1.5f, r * 0.866f, 0.0f, ATOM_C);   // C=O
    addAtom(mol, r * 2.0f, r * 1.5f, 0.0f, ATOM_O);
    addAtom(mol, r * 2.0f, r * 0.0f, 0.0f, ATOM_C);
    addAtom(mol, r * 2.0f, -r * 1.0f, 0.0f, ATOM_C);
    addAtom(mol, r * 1.5f, -r * 0.866f, 0.0f, ATOM_C);  // C=O
    addAtom(mol, r * 2.0f, -r * 1.8f, 0.0f, ATOM_O);
    // Methyl on quinone
    addAtom(mol, r * 3.3f, r * 0.0f, 0.0f, ATOM_C);
    // Phytyl tail (simplified)
    addAtom(mol, r * 3.3f, -r * 1.2f, 0.0f, ATOM_C);
    addAtom(mol, r * 4.6f, -r * 1.2f, 0.0f, ATOM_C);
    addAtom(mol, r * 5.9f, -r * 1.2f, 0.0f, ATOM_C);
    addAtom(mol, r * 7.2f, -r * 1.2f, 0.0f, ATOM_C);
    addAtom(mol, r * 8.5f, -r * 1.2f, 0.0f, ATOM_C);
    // H atoms
    addAtom(mol, -r * 1.7f, 0.0f, 0.0f, ATOM_H);
    addAtom(mol, -r * 0.8f, r * 1.6f, 0.0f, ATOM_H);
    addAtom(mol, -r * 0.8f, -r * 1.6f, 0.0f, ATOM_H);

    addBond(mol, 0, 1, 2);
    addBond(mol, 1, 2, 1);
    addBond(mol, 2, 3, 2);
    addBond(mol, 3, 4, 1);
    addBond(mol, 4, 5, 2);
    addBond(mol, 5, 0, 1);
    addBond(mol, 0, 6, 1);
    addBond(mol, 6, 7, 2);
    addBond(mol, 6, 8, 1);
    addBond(mol, 8, 9, 2);
    addBond(mol, 9, 10, 1);
    addBond(mol, 10, 11, 2);
    addBond(mol, 10, 5, 1);
    addBond(mol, 8, 12, 1);
    addBond(mol, 9, 13, 1);
    addBond(mol, 13, 14, 1);
    addBond(mol, 14, 15, 2);
    addBond(mol, 15, 16, 1);
    addBond(mol, 16, 17, 1);
    addBond(mol, 3, 18, 1);
    addBond(mol, 2, 19, 1);
    addBond(mol, 4, 20, 1);

    centerMolecule(mol);
}

// Build Nicotinamide (Niacinamide) - B3 form (C6H6N2O)
void buildNicotinamide(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Nicotinamide (C6H6N2O)");

    // Pyridine ring
    float r = 1.4f;
    addAtom(mol, r * cosf(0), r * sinf(0), 0.0f, ATOM_N);
    addAtom(mol, r * cosf(PI/3), r * sinf(PI/3), 0.0f, ATOM_C);
    addAtom(mol, r * cosf(2*PI/3), r * sinf(2*PI/3), 0.0f, ATOM_C);
    addAtom(mol, r * cosf(PI), r * sinf(PI), 0.0f, ATOM_C);
    addAtom(mol, r * cosf(4*PI/3), r * sinf(4*PI/3), 0.0f, ATOM_C);
    addAtom(mol, r * cosf(5*PI/3), r * sinf(5*PI/3), 0.0f, ATOM_C);
    // Amide on C3
    addAtom(mol, r * cosf(2*PI/3) * 1.8f, r * sinf(2*PI/3) * 1.8f + 0.6f, 0.0f, ATOM_C);
    addAtom(mol, r * cosf(2*PI/3) * 2.2f, r * sinf(2*PI/3) * 2.2f + 1.5f, 0.0f, ATOM_O);
    addAtom(mol, r * cosf(2*PI/3) * 2.4f, r * sinf(2*PI/3) * 1.2f, 0.0f, ATOM_N);
    // H atoms
    addAtom(mol, r * cosf(2*PI/3) * 3.0f, r * sinf(2*PI/3) * 0.8f, 0.0f, ATOM_H);
    addAtom(mol, r * cosf(2*PI/3) * 2.8f, r * sinf(2*PI/3) * 1.8f, 0.0f, ATOM_H);
    float rH = 2.4f;
    addAtom(mol, rH * cosf(PI/3), rH * sinf(PI/3), 0.0f, ATOM_H);
    addAtom(mol, rH * cosf(PI), rH * sinf(PI), 0.0f, ATOM_H);
    addAtom(mol, rH * cosf(4*PI/3), rH * sinf(4*PI/3), 0.0f, ATOM_H);
    addAtom(mol, rH * cosf(5*PI/3), rH * sinf(5*PI/3), 0.0f, ATOM_H);

    addBond(mol, 0, 1, 2);
    addBond(mol, 1, 2, 1);
    addBond(mol, 2, 3, 2);
    addBond(mol, 3, 4, 1);
    addBond(mol, 4, 5, 2);
    addBond(mol, 5, 0, 1);
    addBond(mol, 2, 6, 1);
    addBond(mol, 6, 7, 2);
    addBond(mol, 6, 8, 1);
    addBond(mol, 8, 9, 1);
    addBond(mol, 8, 10, 1);
    addBond(mol, 1, 11, 1);
    addBond(mol, 3, 12, 1);
    addBond(mol, 4, 13, 1);
    addBond(mol, 5, 14, 1);

    centerMolecule(mol);
}

// ============== ADDITIONAL COMPOUNDS ==============

// Build Cocaine (C17H21NO4) - Tropane alkaloid
void buildCocaine(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Cocaine (C17H21NO4)");

    // Tropane ring system (bicyclic)
    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_N);      // N (bridgehead)
    addAtom(mol, 1.3f, 0.7f, 0.0f, ATOM_C);      // C2
    addAtom(mol, 2.6f, 0.0f, 0.0f, ATOM_C);      // C3 (with ester)
    addAtom(mol, 2.6f, -1.4f, 0.0f, ATOM_C);     // C4
    addAtom(mol, 1.3f, -2.1f, 0.0f, ATOM_C);     // C5
    addAtom(mol, 0.0f, -1.4f, 0.0f, ATOM_C);     // C6
    addAtom(mol, -1.3f, -0.7f, 0.0f, ATOM_C);    // C7
    // Bridge carbon
    addAtom(mol, 1.3f, -0.7f, 1.2f, ATOM_C);     // C1 (bridge)
    // N-methyl
    addAtom(mol, -0.5f, 1.0f, 0.0f, ATOM_C);     // CH3 on N
    // Benzoyl ester at C3
    addAtom(mol, 3.9f, 0.7f, 0.0f, ATOM_C);      // C=O ester
    addAtom(mol, 3.9f, 1.9f, 0.0f, ATOM_O);      // =O
    addAtom(mol, 5.2f, 0.0f, 0.0f, ATOM_O);      // O-benzene
    // Benzene ring
    addAtom(mol, 6.5f, 0.7f, 0.0f, ATOM_C);
    addAtom(mol, 7.7f, 0.0f, 0.0f, ATOM_C);
    addAtom(mol, 8.9f, 0.7f, 0.0f, ATOM_C);
    addAtom(mol, 8.9f, 2.1f, 0.0f, ATOM_C);
    addAtom(mol, 7.7f, 2.8f, 0.0f, ATOM_C);
    addAtom(mol, 6.5f, 2.1f, 0.0f, ATOM_C);
    // Methyl ester at C2
    addAtom(mol, 1.3f, 2.1f, 0.0f, ATOM_C);      // C=O
    addAtom(mol, 0.2f, 2.7f, 0.0f, ATOM_O);      // =O
    addAtom(mol, 2.4f, 2.7f, 0.0f, ATOM_O);      // O-CH3
    addAtom(mol, 2.4f, 4.0f, 0.0f, ATOM_C);      // CH3

    addBond(mol, 0, 1, 1);
    addBond(mol, 1, 2, 1);
    addBond(mol, 2, 3, 1);
    addBond(mol, 3, 4, 1);
    addBond(mol, 4, 5, 1);
    addBond(mol, 5, 6, 1);
    addBond(mol, 6, 0, 1);
    addBond(mol, 1, 7, 1);
    addBond(mol, 4, 7, 1);
    addBond(mol, 0, 8, 1);
    addBond(mol, 2, 9, 1);
    addBond(mol, 9, 10, 2);
    addBond(mol, 9, 11, 1);
    addBond(mol, 11, 12, 1);
    addBond(mol, 12, 13, 2);
    addBond(mol, 13, 14, 1);
    addBond(mol, 14, 15, 2);
    addBond(mol, 15, 16, 1);
    addBond(mol, 16, 17, 2);
    addBond(mol, 17, 12, 1);
    addBond(mol, 1, 18, 1);
    addBond(mol, 18, 19, 2);
    addBond(mol, 18, 20, 1);
    addBond(mol, 20, 21, 1);

    centerMolecule(mol);
}

// Build Heroin / Diacetylmorphine (C21H23NO5)
void buildHeroin(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Heroin (C21H23NO5)");

    // Morphine core - pentacyclic structure (simplified)
    // Ring A (benzene)
    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_C);
    addAtom(mol, 1.2f, 0.7f, 0.0f, ATOM_C);
    addAtom(mol, 2.4f, 0.0f, 0.0f, ATOM_C);
    addAtom(mol, 2.4f, -1.4f, 0.0f, ATOM_C);
    addAtom(mol, 1.2f, -2.1f, 0.0f, ATOM_C);
    addAtom(mol, 0.0f, -1.4f, 0.0f, ATOM_C);
    // Ring B (cyclohexene fused)
    addAtom(mol, -1.3f, 0.7f, 0.0f, ATOM_C);
    addAtom(mol, -1.3f, 2.1f, 0.0f, ATOM_C);
    addAtom(mol, 0.0f, 2.8f, 0.0f, ATOM_C);
    addAtom(mol, 1.2f, 2.1f, 0.0f, ATOM_C);
    // Ring C with oxygen bridge
    addAtom(mol, -2.6f, 0.0f, 0.0f, ATOM_C);
    addAtom(mol, -2.6f, -1.4f, 0.0f, ATOM_O);    // O bridge
    addAtom(mol, -1.3f, -2.1f, 0.0f, ATOM_C);
    // Ring D (piperidine with N)
    addAtom(mol, -2.6f, 2.8f, 0.0f, ATOM_C);
    addAtom(mol, -3.9f, 2.1f, 0.0f, ATOM_C);
    addAtom(mol, -3.9f, 0.7f, 0.0f, ATOM_N);     // N
    addAtom(mol, -4.5f, 0.0f, 0.0f, ATOM_C);     // N-CH3
    // Acetyl groups (diacetyl = heroin)
    addAtom(mol, 3.6f, 0.7f, 0.0f, ATOM_O);      // O-acetyl 1
    addAtom(mol, 4.8f, 0.0f, 0.0f, ATOM_C);      // C=O
    addAtom(mol, 4.8f, -1.2f, 0.0f, ATOM_O);     // =O
    addAtom(mol, 6.0f, 0.7f, 0.0f, ATOM_C);      // CH3
    addAtom(mol, 3.6f, -2.1f, 0.0f, ATOM_O);     // O-acetyl 2
    addAtom(mol, 4.8f, -2.8f, 0.0f, ATOM_C);     // C=O
    addAtom(mol, 4.8f, -4.0f, 0.0f, ATOM_O);     // =O
    addAtom(mol, 6.0f, -2.1f, 0.0f, ATOM_C);     // CH3

    addBond(mol, 0, 1, 2);
    addBond(mol, 1, 2, 1);
    addBond(mol, 2, 3, 2);
    addBond(mol, 3, 4, 1);
    addBond(mol, 4, 5, 2);
    addBond(mol, 5, 0, 1);
    addBond(mol, 0, 6, 1);
    addBond(mol, 6, 7, 1);
    addBond(mol, 7, 8, 2);
    addBond(mol, 8, 9, 1);
    addBond(mol, 9, 1, 1);
    addBond(mol, 6, 10, 1);
    addBond(mol, 10, 11, 1);
    addBond(mol, 11, 12, 1);
    addBond(mol, 12, 5, 1);
    addBond(mol, 7, 13, 1);
    addBond(mol, 13, 14, 1);
    addBond(mol, 14, 15, 1);
    addBond(mol, 15, 10, 1);
    addBond(mol, 15, 16, 1);
    addBond(mol, 2, 17, 1);
    addBond(mol, 17, 18, 1);
    addBond(mol, 18, 19, 2);
    addBond(mol, 18, 20, 1);
    addBond(mol, 3, 21, 1);
    addBond(mol, 21, 22, 1);
    addBond(mol, 22, 23, 2);
    addBond(mol, 22, 24, 1);

    centerMolecule(mol);
}

// Build Fentanyl (C22H28N2O) - Synthetic opioid
void buildFentanyl(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Fentanyl/Sublimaze (C22H28N2O)");

    // Piperidine ring
    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_N);      // N1
    addAtom(mol, 1.3f, 0.7f, 0.0f, ATOM_C);      // C2
    addAtom(mol, 2.6f, 0.0f, 0.0f, ATOM_C);      // C3
    addAtom(mol, 2.6f, -1.4f, 0.0f, ATOM_C);     // C4 (with phenethyl)
    addAtom(mol, 1.3f, -2.1f, 0.0f, ATOM_C);     // C5
    addAtom(mol, 0.0f, -1.4f, 0.0f, ATOM_C);     // C6
    // Phenethyl group on C4
    addAtom(mol, 3.9f, -2.1f, 0.0f, ATOM_C);     // CH2
    addAtom(mol, 5.2f, -1.4f, 0.0f, ATOM_C);     // CH2
    // Benzene ring
    addAtom(mol, 6.5f, -2.1f, 0.0f, ATOM_C);
    addAtom(mol, 7.7f, -1.4f, 0.0f, ATOM_C);
    addAtom(mol, 8.9f, -2.1f, 0.0f, ATOM_C);
    addAtom(mol, 8.9f, -3.5f, 0.0f, ATOM_C);
    addAtom(mol, 7.7f, -4.2f, 0.0f, ATOM_C);
    addAtom(mol, 6.5f, -3.5f, 0.0f, ATOM_C);
    // N-phenyl-propanamide on N1
    addAtom(mol, -1.3f, 0.7f, 0.0f, ATOM_C);     // C=O
    addAtom(mol, -1.3f, 2.0f, 0.0f, ATOM_O);     // =O
    addAtom(mol, -2.6f, 0.0f, 0.0f, ATOM_C);     // CH2CH3
    addAtom(mol, -3.9f, 0.7f, 0.0f, ATOM_C);     // CH3
    // Aniline (N-phenyl)
    addAtom(mol, -1.3f, -2.1f, 0.0f, ATOM_N);    // N-aniline
    addAtom(mol, -2.6f, -2.8f, 0.0f, ATOM_C);    // benzene
    addAtom(mol, -2.6f, -4.2f, 0.0f, ATOM_C);
    addAtom(mol, -3.8f, -4.9f, 0.0f, ATOM_C);
    addAtom(mol, -5.0f, -4.2f, 0.0f, ATOM_C);
    addAtom(mol, -5.0f, -2.8f, 0.0f, ATOM_C);
    addAtom(mol, -3.8f, -2.1f, 0.0f, ATOM_C);

    addBond(mol, 0, 1, 1);
    addBond(mol, 1, 2, 1);
    addBond(mol, 2, 3, 1);
    addBond(mol, 3, 4, 1);
    addBond(mol, 4, 5, 1);
    addBond(mol, 5, 0, 1);
    addBond(mol, 3, 6, 1);
    addBond(mol, 6, 7, 1);
    addBond(mol, 7, 8, 1);
    addBond(mol, 8, 9, 2);
    addBond(mol, 9, 10, 1);
    addBond(mol, 10, 11, 2);
    addBond(mol, 11, 12, 1);
    addBond(mol, 12, 13, 2);
    addBond(mol, 13, 8, 1);
    addBond(mol, 0, 14, 1);
    addBond(mol, 14, 15, 2);
    addBond(mol, 14, 16, 1);
    addBond(mol, 16, 17, 1);
    addBond(mol, 5, 18, 1);
    addBond(mol, 18, 19, 1);
    addBond(mol, 19, 20, 2);
    addBond(mol, 20, 21, 1);
    addBond(mol, 21, 22, 2);
    addBond(mol, 22, 23, 1);
    addBond(mol, 23, 24, 2);
    addBond(mol, 24, 19, 1);

    centerMolecule(mol);
}

// Build Propofol (C12H18O) - General anesthetic
void buildPropofol(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Propofol/Diprivan (C12H18O)");

    // Benzene ring with OH and two isopropyl groups
    float r = 1.4f;
    addAtom(mol, r * cosf(0), r * sinf(0), 0.0f, ATOM_C);           // C1 (OH)
    addAtom(mol, r * cosf(PI/3), r * sinf(PI/3), 0.0f, ATOM_C);     // C2 (iPr)
    addAtom(mol, r * cosf(2*PI/3), r * sinf(2*PI/3), 0.0f, ATOM_C); // C3
    addAtom(mol, r * cosf(PI), r * sinf(PI), 0.0f, ATOM_C);         // C4
    addAtom(mol, r * cosf(4*PI/3), r * sinf(4*PI/3), 0.0f, ATOM_C); // C5
    addAtom(mol, r * cosf(5*PI/3), r * sinf(5*PI/3), 0.0f, ATOM_C); // C6 (iPr)
    // OH on C1
    addAtom(mol, r * 1.8f * cosf(0), r * 1.8f * sinf(0), 0.0f, ATOM_O);
    addAtom(mol, r * 2.5f * cosf(0), r * 2.5f * sinf(0), 0.0f, ATOM_H);
    // Isopropyl on C2
    addAtom(mol, r * 1.8f * cosf(PI/3), r * 1.8f * sinf(PI/3), 0.0f, ATOM_C);  // CH
    addAtom(mol, r * 2.2f * cosf(PI/4), r * 2.2f * sinf(PI/4), 0.8f, ATOM_C);  // CH3
    addAtom(mol, r * 2.2f * cosf(PI/4), r * 2.2f * sinf(PI/4), -0.8f, ATOM_C); // CH3
    // Isopropyl on C6
    addAtom(mol, r * 1.8f * cosf(5*PI/3), r * 1.8f * sinf(5*PI/3), 0.0f, ATOM_C);  // CH
    addAtom(mol, r * 2.2f * cosf(11*PI/6), r * 2.2f * sinf(11*PI/6), 0.8f, ATOM_C);  // CH3
    addAtom(mol, r * 2.2f * cosf(11*PI/6), r * 2.2f * sinf(11*PI/6), -0.8f, ATOM_C); // CH3
    // H atoms on ring
    addAtom(mol, r * 1.7f * cosf(2*PI/3), r * 1.7f * sinf(2*PI/3), 0.0f, ATOM_H);
    addAtom(mol, r * 1.7f * cosf(PI), r * 1.7f * sinf(PI), 0.0f, ATOM_H);
    addAtom(mol, r * 1.7f * cosf(4*PI/3), r * 1.7f * sinf(4*PI/3), 0.0f, ATOM_H);

    addBond(mol, 0, 1, 2);
    addBond(mol, 1, 2, 1);
    addBond(mol, 2, 3, 2);
    addBond(mol, 3, 4, 1);
    addBond(mol, 4, 5, 2);
    addBond(mol, 5, 0, 1);
    addBond(mol, 0, 6, 1);
    addBond(mol, 6, 7, 1);
    addBond(mol, 1, 8, 1);
    addBond(mol, 8, 9, 1);
    addBond(mol, 8, 10, 1);
    addBond(mol, 5, 11, 1);
    addBond(mol, 11, 12, 1);
    addBond(mol, 11, 13, 1);
    addBond(mol, 2, 14, 1);
    addBond(mol, 3, 15, 1);
    addBond(mol, 4, 16, 1);

    centerMolecule(mol);
}

// Build THC - Tetrahydrocannabinol (C21H30O2)
void buildTHC(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "THC (C21H30O2)");

    // Dibenzopyran core
    // Benzene ring A
    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_C);
    addAtom(mol, 1.2f, 0.7f, 0.0f, ATOM_C);
    addAtom(mol, 2.4f, 0.0f, 0.0f, ATOM_C);
    addAtom(mol, 2.4f, -1.4f, 0.0f, ATOM_C);
    addAtom(mol, 1.2f, -2.1f, 0.0f, ATOM_C);
    addAtom(mol, 0.0f, -1.4f, 0.0f, ATOM_C);
    // Pyran ring (with O)
    addAtom(mol, -1.2f, 0.7f, 0.0f, ATOM_C);
    addAtom(mol, -2.4f, 0.0f, 0.0f, ATOM_O);     // O in pyran
    addAtom(mol, -2.4f, -1.4f, 0.0f, ATOM_C);
    addAtom(mol, -1.2f, -2.1f, 0.0f, ATOM_C);
    // Cyclohexene ring C
    addAtom(mol, -1.2f, -3.5f, 0.0f, ATOM_C);
    addAtom(mol, 0.0f, -4.2f, 0.0f, ATOM_C);
    addAtom(mol, 1.2f, -3.5f, 0.0f, ATOM_C);
    // OH on benzene
    addAtom(mol, 3.6f, 0.7f, 0.0f, ATOM_O);
    addAtom(mol, 4.3f, 0.2f, 0.0f, ATOM_H);
    // Pentyl chain
    addAtom(mol, 3.6f, -2.1f, 0.0f, ATOM_C);
    addAtom(mol, 4.9f, -1.4f, 0.0f, ATOM_C);
    addAtom(mol, 6.2f, -2.1f, 0.0f, ATOM_C);
    addAtom(mol, 7.5f, -1.4f, 0.0f, ATOM_C);
    addAtom(mol, 8.8f, -2.1f, 0.0f, ATOM_C);
    // Methyls on cyclohexene
    addAtom(mol, -2.4f, -4.2f, 0.0f, ATOM_C);    // gem-dimethyl
    addAtom(mol, -2.4f, -3.0f, 1.0f, ATOM_C);
    addAtom(mol, -1.2f, 2.0f, 0.0f, ATOM_C);     // CH3 on C

    addBond(mol, 0, 1, 2);
    addBond(mol, 1, 2, 1);
    addBond(mol, 2, 3, 2);
    addBond(mol, 3, 4, 1);
    addBond(mol, 4, 5, 2);
    addBond(mol, 5, 0, 1);
    addBond(mol, 0, 6, 1);
    addBond(mol, 6, 7, 1);
    addBond(mol, 7, 8, 1);
    addBond(mol, 8, 9, 1);
    addBond(mol, 9, 5, 1);
    addBond(mol, 9, 10, 1);
    addBond(mol, 10, 11, 2);
    addBond(mol, 11, 12, 1);
    addBond(mol, 12, 4, 1);
    addBond(mol, 2, 13, 1);
    addBond(mol, 13, 14, 1);
    addBond(mol, 3, 15, 1);
    addBond(mol, 15, 16, 1);
    addBond(mol, 16, 17, 1);
    addBond(mol, 17, 18, 1);
    addBond(mol, 18, 19, 1);
    addBond(mol, 10, 20, 1);
    addBond(mol, 10, 21, 1);
    addBond(mol, 6, 22, 1);

    centerMolecule(mol);
}

// Build Creatine (C4H9N3O2) - Sports supplement
void buildCreatine(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Creatine (C4H9N3O2)");

    // Guanidino group
    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_C);      // Central C of guanidino
    addAtom(mol, 0.0f, 1.2f, 0.0f, ATOM_N);      // =NH
    addAtom(mol, -1.1f, -0.7f, 0.0f, ATOM_N);    // NH2
    addAtom(mol, 1.1f, -0.7f, 0.0f, ATOM_N);     // N-CH3
    addAtom(mol, 2.4f, 0.0f, 0.0f, ATOM_C);      // CH3 on N
    // Acetic acid part
    addAtom(mol, 1.1f, -2.1f, 0.0f, ATOM_C);     // CH2
    addAtom(mol, 2.4f, -2.8f, 0.0f, ATOM_C);     // COOH
    addAtom(mol, 2.4f, -4.0f, 0.0f, ATOM_O);     // =O
    addAtom(mol, 3.5f, -2.1f, 0.0f, ATOM_O);     // OH
    // H atoms
    addAtom(mol, 0.0f, 2.0f, 0.0f, ATOM_H);      // H on =NH
    addAtom(mol, -1.8f, -0.2f, 0.0f, ATOM_H);    // H on NH2
    addAtom(mol, -1.4f, -1.5f, 0.0f, ATOM_H);    // H on NH2
    addAtom(mol, 4.3f, -2.5f, 0.0f, ATOM_H);     // H on OH

    addBond(mol, 0, 1, 2);
    addBond(mol, 0, 2, 1);
    addBond(mol, 0, 3, 1);
    addBond(mol, 3, 4, 1);
    addBond(mol, 3, 5, 1);
    addBond(mol, 5, 6, 1);
    addBond(mol, 6, 7, 2);
    addBond(mol, 6, 8, 1);
    addBond(mol, 1, 9, 1);
    addBond(mol, 2, 10, 1);
    addBond(mol, 2, 11, 1);
    addBond(mol, 8, 12, 1);

    centerMolecule(mol);
}

// Build Octane (C8H18) - Representative gasoline component
void buildOctane(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Octane/Gasoline (C8H18)");

    // Linear 8-carbon chain
    for (int i = 0; i < 8; i++) {
        float x = i * 1.5f;
        float y = (i % 2 == 0) ? 0.0f : 0.3f;  // Slight zigzag
        addAtom(mol, x, y, 0.0f, ATOM_C);
    }
    // H atoms on each carbon (simplified - 3H on terminals, 2H on others)
    // Terminal CH3 groups
    addAtom(mol, -0.8f, 0.5f, 0.0f, ATOM_H);
    addAtom(mol, -0.8f, -0.5f, 0.0f, ATOM_H);
    addAtom(mol, -0.5f, 0.0f, 0.8f, ATOM_H);
    addAtom(mol, 11.3f, 0.8f, 0.0f, ATOM_H);
    addAtom(mol, 11.3f, -0.2f, 0.0f, ATOM_H);
    addAtom(mol, 11.0f, 0.3f, 0.8f, ATOM_H);
    // CH2 groups (2H each for carbons 1-6)
    for (int i = 1; i < 7; i++) {
        float x = i * 1.5f;
        float y = (i % 2 == 0) ? 0.0f : 0.3f;
        addAtom(mol, x, y + 0.8f, 0.0f, ATOM_H);
        addAtom(mol, x, y - 0.8f, 0.0f, ATOM_H);
    }

    // C-C backbone bonds
    for (int i = 0; i < 7; i++) {
        addBond(mol, i, i + 1, 1);
    }
    // C-H bonds for terminal methyls
    addBond(mol, 0, 8, 1);
    addBond(mol, 0, 9, 1);
    addBond(mol, 0, 10, 1);
    addBond(mol, 7, 11, 1);
    addBond(mol, 7, 12, 1);
    addBond(mol, 7, 13, 1);
    // C-H bonds for CH2 groups
    for (int i = 0; i < 6; i++) {
        addBond(mol, i + 1, 14 + i * 2, 1);
        addBond(mol, i + 1, 15 + i * 2, 1);
    }

    centerMolecule(mol);
}

// ============== STATINS & NSAIDs ==============

// Build Simvastatin (C25H38O5) - Statin cholesterol drug
void buildSimvastatin(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Simvastatin/Zocor (C25H38O5)");

    // Decalin-like fused ring system (simplified)
    // Ring A (cyclohexene)
    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_C);      // C1
    addAtom(mol, 1.3f, 0.7f, 0.0f, ATOM_C);      // C2
    addAtom(mol, 2.6f, 0.0f, 0.0f, ATOM_C);      // C3
    addAtom(mol, 2.6f, -1.4f, 0.0f, ATOM_C);     // C4
    addAtom(mol, 1.3f, -2.1f, 0.0f, ATOM_C);     // C5
    addAtom(mol, 0.0f, -1.4f, 0.0f, ATOM_C);     // C6
    // Ring B (fused cyclohexane)
    addAtom(mol, -1.3f, 0.7f, 0.0f, ATOM_C);     // C7
    addAtom(mol, -2.6f, 0.0f, 0.0f, ATOM_C);     // C8
    addAtom(mol, -2.6f, -1.4f, 0.0f, ATOM_C);    // C9
    addAtom(mol, -1.3f, -2.1f, 0.0f, ATOM_C);    // C10
    // Lactone ring
    addAtom(mol, 3.9f, 0.7f, 0.0f, ATOM_C);      // C11
    addAtom(mol, 5.2f, 0.0f, 0.0f, ATOM_O);      // O (lactone)
    addAtom(mol, 5.2f, -1.4f, 0.0f, ATOM_C);     // C12 (C=O)
    addAtom(mol, 5.2f, -2.6f, 0.0f, ATOM_O);     // =O
    addAtom(mol, 3.9f, -0.7f, 0.0f, ATOM_C);     // C13
    // Hydroxy groups
    addAtom(mol, 3.9f, 2.0f, 0.0f, ATOM_O);      // OH
    addAtom(mol, 3.9f, -2.0f, 0.0f, ATOM_O);     // OH
    // Side chain ester
    addAtom(mol, -3.9f, 0.7f, 0.0f, ATOM_O);     // O-ester
    addAtom(mol, -5.2f, 0.0f, 0.0f, ATOM_C);     // C=O
    addAtom(mol, -5.2f, -1.2f, 0.0f, ATOM_O);    // =O
    addAtom(mol, -6.5f, 0.7f, 0.0f, ATOM_C);     // CH(CH3)2
    addAtom(mol, -7.8f, 0.0f, 0.0f, ATOM_C);     // CH3
    addAtom(mol, -6.5f, 2.0f, 0.0f, ATOM_C);     // CH3
    // Methyls
    addAtom(mol, 1.3f, 2.0f, 0.0f, ATOM_C);      // CH3
    addAtom(mol, -1.3f, -3.5f, 0.0f, ATOM_C);    // CH3

    addBond(mol, 0, 1, 1);
    addBond(mol, 1, 2, 2);
    addBond(mol, 2, 3, 1);
    addBond(mol, 3, 4, 1);
    addBond(mol, 4, 5, 1);
    addBond(mol, 5, 0, 1);
    addBond(mol, 0, 6, 1);
    addBond(mol, 6, 7, 1);
    addBond(mol, 7, 8, 1);
    addBond(mol, 8, 9, 1);
    addBond(mol, 9, 5, 1);
    addBond(mol, 2, 10, 1);
    addBond(mol, 10, 11, 1);
    addBond(mol, 11, 12, 1);
    addBond(mol, 12, 13, 2);
    addBond(mol, 12, 14, 1);
    addBond(mol, 14, 3, 1);
    addBond(mol, 10, 15, 1);
    addBond(mol, 14, 16, 1);
    addBond(mol, 7, 17, 1);
    addBond(mol, 17, 18, 1);
    addBond(mol, 18, 19, 2);
    addBond(mol, 18, 20, 1);
    addBond(mol, 20, 21, 1);
    addBond(mol, 20, 22, 1);
    addBond(mol, 1, 23, 1);
    addBond(mol, 9, 24, 1);

    centerMolecule(mol);
}

// Build Ibuprofen (C13H18O2) - NSAID
void buildIbuprofen(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Ibuprofen/Advil (C13H18O2)");

    // Benzene ring
    float r = 1.4f;
    addAtom(mol, r * cosf(0), r * sinf(0), 0.0f, ATOM_C);           // C1
    addAtom(mol, r * cosf(PI/3), r * sinf(PI/3), 0.0f, ATOM_C);     // C2
    addAtom(mol, r * cosf(2*PI/3), r * sinf(2*PI/3), 0.0f, ATOM_C); // C3
    addAtom(mol, r * cosf(PI), r * sinf(PI), 0.0f, ATOM_C);         // C4 (isobutyl)
    addAtom(mol, r * cosf(4*PI/3), r * sinf(4*PI/3), 0.0f, ATOM_C); // C5
    addAtom(mol, r * cosf(5*PI/3), r * sinf(5*PI/3), 0.0f, ATOM_C); // C6
    // Propionic acid on C1
    addAtom(mol, r * 1.8f, 0.0f, 0.0f, ATOM_C);   // CH(CH3)
    addAtom(mol, r * 2.4f, 1.0f, 0.0f, ATOM_C);   // CH3
    addAtom(mol, r * 3.0f, -0.7f, 0.0f, ATOM_C);  // COOH
    addAtom(mol, r * 3.0f, -1.9f, 0.0f, ATOM_O);  // =O
    addAtom(mol, r * 4.1f, 0.0f, 0.0f, ATOM_O);   // OH
    addAtom(mol, r * 4.8f, -0.5f, 0.0f, ATOM_H);
    // Isobutyl on C4
    addAtom(mol, r * cosf(PI) * 1.8f, r * sinf(PI) * 1.8f, 0.0f, ATOM_C);  // CH2
    addAtom(mol, r * cosf(PI) * 2.8f, r * sinf(PI) * 1.2f, 0.0f, ATOM_C);  // CH
    addAtom(mol, r * cosf(PI) * 3.5f, r * sinf(PI) * 0.5f, 0.0f, ATOM_C);  // CH3
    addAtom(mol, r * cosf(PI) * 3.5f, r * sinf(PI) * 1.9f, 0.0f, ATOM_C);  // CH3
    // H atoms on ring
    addAtom(mol, r * 1.7f * cosf(PI/3), r * 1.7f * sinf(PI/3), 0.0f, ATOM_H);
    addAtom(mol, r * 1.7f * cosf(2*PI/3), r * 1.7f * sinf(2*PI/3), 0.0f, ATOM_H);
    addAtom(mol, r * 1.7f * cosf(4*PI/3), r * 1.7f * sinf(4*PI/3), 0.0f, ATOM_H);
    addAtom(mol, r * 1.7f * cosf(5*PI/3), r * 1.7f * sinf(5*PI/3), 0.0f, ATOM_H);

    addBond(mol, 0, 1, 2);
    addBond(mol, 1, 2, 1);
    addBond(mol, 2, 3, 2);
    addBond(mol, 3, 4, 1);
    addBond(mol, 4, 5, 2);
    addBond(mol, 5, 0, 1);
    addBond(mol, 0, 6, 1);
    addBond(mol, 6, 7, 1);
    addBond(mol, 6, 8, 1);
    addBond(mol, 8, 9, 2);
    addBond(mol, 8, 10, 1);
    addBond(mol, 10, 11, 1);
    addBond(mol, 3, 12, 1);
    addBond(mol, 12, 13, 1);
    addBond(mol, 13, 14, 1);
    addBond(mol, 13, 15, 1);
    addBond(mol, 1, 16, 1);
    addBond(mol, 2, 17, 1);
    addBond(mol, 4, 18, 1);
    addBond(mol, 5, 19, 1);

    centerMolecule(mol);
}

// Build Naproxen (C14H14O3) - NSAID
void buildNaproxen(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Naproxen/Aleve (C14H14O3)");

    // Naphthalene ring system
    // Ring 1
    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_C);
    addAtom(mol, 1.2f, 0.7f, 0.0f, ATOM_C);
    addAtom(mol, 2.4f, 0.0f, 0.0f, ATOM_C);
    addAtom(mol, 2.4f, -1.4f, 0.0f, ATOM_C);
    addAtom(mol, 1.2f, -2.1f, 0.0f, ATOM_C);
    addAtom(mol, 0.0f, -1.4f, 0.0f, ATOM_C);
    // Ring 2 (fused)
    addAtom(mol, 3.6f, 0.7f, 0.0f, ATOM_C);
    addAtom(mol, 4.8f, 0.0f, 0.0f, ATOM_C);
    addAtom(mol, 4.8f, -1.4f, 0.0f, ATOM_C);
    addAtom(mol, 3.6f, -2.1f, 0.0f, ATOM_C);
    // Methoxy on C6 position
    addAtom(mol, 3.6f, 2.0f, 0.0f, ATOM_O);      // O
    addAtom(mol, 3.6f, 3.3f, 0.0f, ATOM_C);      // CH3
    // Propionic acid on C2 position
    addAtom(mol, -1.3f, 0.7f, 0.0f, ATOM_C);     // CH(CH3)
    addAtom(mol, -1.3f, 2.0f, 0.0f, ATOM_C);     // CH3
    addAtom(mol, -2.6f, 0.0f, 0.0f, ATOM_C);     // COOH
    addAtom(mol, -2.6f, -1.2f, 0.0f, ATOM_O);    // =O
    addAtom(mol, -3.7f, 0.7f, 0.0f, ATOM_O);     // OH
    addAtom(mol, -4.4f, 0.2f, 0.0f, ATOM_H);
    // H atoms
    addAtom(mol, 1.2f, 1.7f, 0.0f, ATOM_H);
    addAtom(mol, 1.2f, -3.1f, 0.0f, ATOM_H);
    addAtom(mol, 5.7f, 0.5f, 0.0f, ATOM_H);
    addAtom(mol, 5.7f, -1.9f, 0.0f, ATOM_H);
    addAtom(mol, 3.6f, -3.1f, 0.0f, ATOM_H);

    addBond(mol, 0, 1, 2);
    addBond(mol, 1, 2, 1);
    addBond(mol, 2, 3, 2);
    addBond(mol, 3, 4, 1);
    addBond(mol, 4, 5, 2);
    addBond(mol, 5, 0, 1);
    addBond(mol, 2, 6, 1);
    addBond(mol, 6, 7, 2);
    addBond(mol, 7, 8, 1);
    addBond(mol, 8, 9, 2);
    addBond(mol, 9, 3, 1);
    addBond(mol, 6, 10, 1);
    addBond(mol, 10, 11, 1);
    addBond(mol, 0, 12, 1);
    addBond(mol, 12, 13, 1);
    addBond(mol, 12, 14, 1);
    addBond(mol, 14, 15, 2);
    addBond(mol, 14, 16, 1);
    addBond(mol, 16, 17, 1);
    addBond(mol, 1, 18, 1);
    addBond(mol, 4, 19, 1);
    addBond(mol, 7, 20, 1);
    addBond(mol, 8, 21, 1);
    addBond(mol, 9, 22, 1);

    centerMolecule(mol);
}

// Build Diclofenac (C14H11Cl2NO2) - NSAID
void buildDiclofenac(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Diclofenac/Voltaren (C14H11Cl2NO2)");

    // Phenylacetic acid ring
    float r = 1.4f;
    addAtom(mol, r * cosf(0), r * sinf(0), 0.0f, ATOM_C);
    addAtom(mol, r * cosf(PI/3), r * sinf(PI/3), 0.0f, ATOM_C);
    addAtom(mol, r * cosf(2*PI/3), r * sinf(2*PI/3), 0.0f, ATOM_C);
    addAtom(mol, r * cosf(PI), r * sinf(PI), 0.0f, ATOM_C);
    addAtom(mol, r * cosf(4*PI/3), r * sinf(4*PI/3), 0.0f, ATOM_C);
    addAtom(mol, r * cosf(5*PI/3), r * sinf(5*PI/3), 0.0f, ATOM_C);
    // Acetic acid on C1
    addAtom(mol, r * 1.8f, 0.0f, 0.0f, ATOM_C);   // CH2
    addAtom(mol, r * 3.0f, 0.0f, 0.0f, ATOM_C);   // COOH
    addAtom(mol, r * 3.5f, 1.0f, 0.0f, ATOM_O);   // =O
    addAtom(mol, r * 3.5f, -1.0f, 0.0f, ATOM_O);  // OH
    addAtom(mol, r * 4.3f, -1.3f, 0.0f, ATOM_H);
    // NH bridge on C2
    addAtom(mol, r * cosf(PI/3) * 1.8f, r * sinf(PI/3) * 1.8f, 0.0f, ATOM_N);
    // Dichlorophenyl ring
    addAtom(mol, r * cosf(PI/3) * 2.8f, r * sinf(PI/3) * 2.8f, 0.0f, ATOM_C);
    addAtom(mol, r * cosf(PI/3) * 3.3f + 0.6f, r * sinf(PI/3) * 3.3f + 1.0f, 0.0f, ATOM_C);
    addAtom(mol, r * cosf(PI/3) * 3.8f + 1.2f, r * sinf(PI/3) * 3.3f + 0.3f, 0.0f, ATOM_C);  // Cl
    addAtom(mol, r * cosf(PI/3) * 3.8f + 1.2f, r * sinf(PI/3) * 2.3f - 0.4f, 0.0f, ATOM_C);
    addAtom(mol, r * cosf(PI/3) * 3.3f + 0.6f, r * sinf(PI/3) * 1.8f - 0.4f, 0.0f, ATOM_C);
    addAtom(mol, r * cosf(PI/3) * 2.8f, r * sinf(PI/3) * 2.3f, 0.0f, ATOM_C);   // Cl
    // Chlorines
    addAtom(mol, r * cosf(PI/3) * 4.5f + 1.8f, r * sinf(PI/3) * 3.5f + 0.5f, 0.0f, ATOM_CL);
    addAtom(mol, r * cosf(PI/3) * 2.3f - 0.6f, r * sinf(PI/3) * 2.3f - 0.6f, 0.0f, ATOM_CL);
    // H on NH
    addAtom(mol, r * cosf(PI/3) * 1.8f + 0.5f, r * sinf(PI/3) * 1.8f + 0.7f, 0.0f, ATOM_H);

    addBond(mol, 0, 1, 2);
    addBond(mol, 1, 2, 1);
    addBond(mol, 2, 3, 2);
    addBond(mol, 3, 4, 1);
    addBond(mol, 4, 5, 2);
    addBond(mol, 5, 0, 1);
    addBond(mol, 0, 6, 1);
    addBond(mol, 6, 7, 1);
    addBond(mol, 7, 8, 2);
    addBond(mol, 7, 9, 1);
    addBond(mol, 9, 10, 1);
    addBond(mol, 1, 11, 1);
    addBond(mol, 11, 12, 1);
    addBond(mol, 12, 13, 2);
    addBond(mol, 13, 14, 1);
    addBond(mol, 14, 15, 2);
    addBond(mol, 15, 16, 1);
    addBond(mol, 16, 17, 2);
    addBond(mol, 17, 12, 1);
    addBond(mol, 14, 18, 1);
    addBond(mol, 17, 19, 1);
    addBond(mol, 11, 20, 1);

    centerMolecule(mol);
}

// Build Indomethacin (C19H16ClNO4) - NSAID
void buildIndomethacin(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Indomethacin/Indocin (C19H16ClNO4)");

    // Indole ring system
    // Benzene
    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_C);
    addAtom(mol, 1.2f, 0.7f, 0.0f, ATOM_C);
    addAtom(mol, 2.4f, 0.0f, 0.0f, ATOM_C);
    addAtom(mol, 2.4f, -1.4f, 0.0f, ATOM_C);
    addAtom(mol, 1.2f, -2.1f, 0.0f, ATOM_C);
    addAtom(mol, 0.0f, -1.4f, 0.0f, ATOM_C);
    // Pyrrole fused
    addAtom(mol, -1.0f, 0.7f, 0.0f, ATOM_C);
    addAtom(mol, -1.0f, -0.7f, 0.0f, ATOM_C);
    addAtom(mol, -2.2f, 0.0f, 0.0f, ATOM_N);
    // Methoxy on benzene
    addAtom(mol, 3.6f, 0.7f, 0.0f, ATOM_O);
    addAtom(mol, 4.8f, 0.0f, 0.0f, ATOM_C);
    // N-p-chlorobenzoyl
    addAtom(mol, -3.4f, 0.7f, 0.0f, ATOM_C);     // C=O
    addAtom(mol, -3.4f, 2.0f, 0.0f, ATOM_O);     // =O
    // p-Chlorobenzene
    addAtom(mol, -4.6f, 0.0f, 0.0f, ATOM_C);
    addAtom(mol, -5.8f, 0.7f, 0.0f, ATOM_C);
    addAtom(mol, -7.0f, 0.0f, 0.0f, ATOM_C);
    addAtom(mol, -7.0f, -1.4f, 0.0f, ATOM_C);    // para-Cl
    addAtom(mol, -5.8f, -2.1f, 0.0f, ATOM_C);
    addAtom(mol, -4.6f, -1.4f, 0.0f, ATOM_C);
    addAtom(mol, -8.2f, -2.1f, 0.0f, ATOM_CL);   // Cl
    // Acetic acid on pyrrole C3
    addAtom(mol, -1.0f, -2.1f, 0.0f, ATOM_C);    // CH2
    addAtom(mol, -1.0f, -3.5f, 0.0f, ATOM_C);    // COOH
    addAtom(mol, 0.1f, -4.1f, 0.0f, ATOM_O);     // =O
    addAtom(mol, -2.1f, -4.1f, 0.0f, ATOM_O);    // OH
    // Methyl on pyrrole C2
    addAtom(mol, -1.0f, 2.0f, 0.0f, ATOM_C);

    addBond(mol, 0, 1, 2);
    addBond(mol, 1, 2, 1);
    addBond(mol, 2, 3, 2);
    addBond(mol, 3, 4, 1);
    addBond(mol, 4, 5, 2);
    addBond(mol, 5, 0, 1);
    addBond(mol, 0, 6, 1);
    addBond(mol, 5, 7, 1);
    addBond(mol, 6, 8, 1);
    addBond(mol, 8, 7, 1);
    addBond(mol, 2, 9, 1);
    addBond(mol, 9, 10, 1);
    addBond(mol, 8, 11, 1);
    addBond(mol, 11, 12, 2);
    addBond(mol, 11, 13, 1);
    addBond(mol, 13, 14, 2);
    addBond(mol, 14, 15, 1);
    addBond(mol, 15, 16, 2);
    addBond(mol, 16, 17, 1);
    addBond(mol, 17, 18, 2);
    addBond(mol, 18, 13, 1);
    addBond(mol, 16, 19, 1);
    addBond(mol, 7, 20, 1);
    addBond(mol, 20, 21, 1);
    addBond(mol, 21, 22, 2);
    addBond(mol, 21, 23, 1);
    addBond(mol, 6, 24, 1);

    centerMolecule(mol);
}

// Build Celecoxib (C17H14F3N3O2S) - COX-2 selective NSAID
void buildCelecoxib(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Celecoxib/Celebrex (C17H14F3N3O2S)");

    // Pyrazole ring
    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_N);      // N1
    addAtom(mol, 1.0f, 0.8f, 0.0f, ATOM_N);      // N2
    addAtom(mol, 2.2f, 0.2f, 0.0f, ATOM_C);      // C3
    addAtom(mol, 2.0f, -1.2f, 0.0f, ATOM_C);     // C4
    addAtom(mol, 0.6f, -1.2f, 0.0f, ATOM_C);     // C5
    // Tolyl (4-methylphenyl) on C3
    addAtom(mol, 3.5f, 0.9f, 0.0f, ATOM_C);      // benzene
    addAtom(mol, 4.7f, 0.2f, 0.0f, ATOM_C);
    addAtom(mol, 5.9f, 0.9f, 0.0f, ATOM_C);
    addAtom(mol, 5.9f, 2.3f, 0.0f, ATOM_C);      // para-CH3
    addAtom(mol, 4.7f, 3.0f, 0.0f, ATOM_C);
    addAtom(mol, 3.5f, 2.3f, 0.0f, ATOM_C);
    addAtom(mol, 7.1f, 3.0f, 0.0f, ATOM_C);      // CH3
    // Trifluoromethylphenyl on C5
    addAtom(mol, 0.0f, -2.5f, 0.0f, ATOM_C);     // benzene
    addAtom(mol, -1.2f, -3.2f, 0.0f, ATOM_C);
    addAtom(mol, -1.2f, -4.6f, 0.0f, ATOM_C);
    addAtom(mol, 0.0f, -5.3f, 0.0f, ATOM_C);     // para-CF3
    addAtom(mol, 1.2f, -4.6f, 0.0f, ATOM_C);
    addAtom(mol, 1.2f, -3.2f, 0.0f, ATOM_C);
    // CF3 group
    addAtom(mol, 0.0f, -6.7f, 0.0f, ATOM_C);
    addAtom(mol, -1.1f, -7.3f, 0.0f, ATOM_F);
    addAtom(mol, 0.0f, -7.9f, 0.0f, ATOM_F);
    addAtom(mol, 1.1f, -7.3f, 0.0f, ATOM_F);
    // Sulfonamide on N1
    addAtom(mol, -1.3f, 0.7f, 0.0f, ATOM_S);     // S
    addAtom(mol, -1.3f, 2.0f, 0.0f, ATOM_O);     // =O
    addAtom(mol, -2.5f, 0.0f, 0.0f, ATOM_O);     // =O
    addAtom(mol, -1.3f, 0.7f, 1.5f, ATOM_N);     // NH2

    addBond(mol, 0, 1, 1);
    addBond(mol, 1, 2, 2);
    addBond(mol, 2, 3, 1);
    addBond(mol, 3, 4, 2);
    addBond(mol, 4, 0, 1);
    addBond(mol, 2, 5, 1);
    addBond(mol, 5, 6, 2);
    addBond(mol, 6, 7, 1);
    addBond(mol, 7, 8, 2);
    addBond(mol, 8, 9, 1);
    addBond(mol, 9, 10, 2);
    addBond(mol, 10, 5, 1);
    addBond(mol, 8, 11, 1);
    addBond(mol, 4, 12, 1);
    addBond(mol, 12, 13, 2);
    addBond(mol, 13, 14, 1);
    addBond(mol, 14, 15, 2);
    addBond(mol, 15, 16, 1);
    addBond(mol, 16, 17, 2);
    addBond(mol, 17, 12, 1);
    addBond(mol, 15, 18, 1);
    addBond(mol, 18, 19, 1);
    addBond(mol, 18, 20, 1);
    addBond(mol, 18, 21, 1);
    addBond(mol, 0, 22, 1);
    addBond(mol, 22, 23, 2);
    addBond(mol, 22, 24, 2);
    addBond(mol, 22, 25, 1);

    centerMolecule(mol);
}

// Build Meloxicam (C14H13N3O4S2) - NSAID
void buildMeloxicam(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Meloxicam/Mobic (C14H13N3O4S2)");

    // Benzothiazine core
    // Benzene ring
    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_C);
    addAtom(mol, 1.2f, 0.7f, 0.0f, ATOM_C);
    addAtom(mol, 2.4f, 0.0f, 0.0f, ATOM_C);
    addAtom(mol, 2.4f, -1.4f, 0.0f, ATOM_C);
    addAtom(mol, 1.2f, -2.1f, 0.0f, ATOM_C);
    addAtom(mol, 0.0f, -1.4f, 0.0f, ATOM_C);
    // Thiazine fused
    addAtom(mol, -1.2f, 0.7f, 0.0f, ATOM_S);     // S
    addAtom(mol, -2.4f, 0.0f, 0.0f, ATOM_C);     // C
    addAtom(mol, -2.4f, -1.4f, 0.0f, ATOM_C);    // C (enol)
    addAtom(mol, -1.2f, -2.1f, 0.0f, ATOM_N);    // N
    // Sulfonamide SO2 on S
    addAtom(mol, -1.2f, 2.0f, 0.0f, ATOM_O);     // =O
    addAtom(mol, -0.5f, 0.7f, 1.2f, ATOM_O);     // =O
    // Amide C=O
    addAtom(mol, -3.6f, 0.7f, 0.0f, ATOM_O);     // =O on C8
    // Thiazole ring on enol C
    addAtom(mol, -3.6f, -2.1f, 0.0f, ATOM_N);
    addAtom(mol, -4.8f, -1.4f, 0.0f, ATOM_C);
    addAtom(mol, -6.0f, -2.1f, 0.0f, ATOM_S);
    addAtom(mol, -5.5f, -3.5f, 0.0f, ATOM_C);
    addAtom(mol, -4.0f, -3.5f, 0.0f, ATOM_C);
    // Methyl on thiazole
    addAtom(mol, -4.8f, 0.0f, 0.0f, ATOM_C);     // CH3
    // OH (enol)
    addAtom(mol, -2.4f, -2.6f, 0.0f, ATOM_O);
    // N-methyl
    addAtom(mol, -1.2f, -3.5f, 0.0f, ATOM_C);

    addBond(mol, 0, 1, 2);
    addBond(mol, 1, 2, 1);
    addBond(mol, 2, 3, 2);
    addBond(mol, 3, 4, 1);
    addBond(mol, 4, 5, 2);
    addBond(mol, 5, 0, 1);
    addBond(mol, 0, 6, 1);
    addBond(mol, 6, 7, 1);
    addBond(mol, 7, 8, 2);
    addBond(mol, 8, 9, 1);
    addBond(mol, 9, 5, 1);
    addBond(mol, 6, 10, 2);
    addBond(mol, 6, 11, 2);
    addBond(mol, 7, 12, 2);
    addBond(mol, 8, 13, 1);
    addBond(mol, 13, 14, 2);
    addBond(mol, 14, 15, 1);
    addBond(mol, 15, 16, 1);
    addBond(mol, 16, 17, 2);
    addBond(mol, 17, 13, 1);
    addBond(mol, 14, 18, 1);
    addBond(mol, 8, 19, 1);
    addBond(mol, 9, 20, 1);

    centerMolecule(mol);
}

// Build Acetaminophen/Paracetamol (C8H9NO2) - Tylenol
void buildAcetaminophen(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Acetaminophen/Tylenol (C8H9NO2)");

    // Benzene ring (para-substituted)
    float r = 1.4f;
    addAtom(mol, r * cosf(0), r * sinf(0), 0.0f, ATOM_C);           // C1 (OH)
    addAtom(mol, r * cosf(PI/3), r * sinf(PI/3), 0.0f, ATOM_C);     // C2
    addAtom(mol, r * cosf(2*PI/3), r * sinf(2*PI/3), 0.0f, ATOM_C); // C3
    addAtom(mol, r * cosf(PI), r * sinf(PI), 0.0f, ATOM_C);         // C4 (NHCOCH3)
    addAtom(mol, r * cosf(4*PI/3), r * sinf(4*PI/3), 0.0f, ATOM_C); // C5
    addAtom(mol, r * cosf(5*PI/3), r * sinf(5*PI/3), 0.0f, ATOM_C); // C6
    // Hydroxyl on C1 (para to amide)
    addAtom(mol, r * 1.8f * cosf(0), r * 1.8f * sinf(0), 0.0f, ATOM_O);  // OH
    addAtom(mol, r * 2.5f * cosf(0), r * 2.5f * sinf(0), 0.0f, ATOM_H);  // H
    // Acetamide on C4 (para position)
    addAtom(mol, r * 1.8f * cosf(PI), r * 1.8f * sinf(PI), 0.0f, ATOM_N);   // NH
    addAtom(mol, r * 1.8f * cosf(PI) - 0.5f, r * 1.8f * sinf(PI) + 0.8f, 0.0f, ATOM_H);  // H on N
    addAtom(mol, r * 2.8f * cosf(PI), r * 2.8f * sinf(PI), 0.0f, ATOM_C);   // C=O
    addAtom(mol, r * 2.8f * cosf(PI), r * 2.8f * sinf(PI) + 1.2f, 0.0f, ATOM_O);  // =O
    addAtom(mol, r * 4.0f * cosf(PI), r * 4.0f * sinf(PI), 0.0f, ATOM_C);   // CH3
    // H atoms on ring
    addAtom(mol, r * 1.7f * cosf(PI/3), r * 1.7f * sinf(PI/3), 0.0f, ATOM_H);
    addAtom(mol, r * 1.7f * cosf(2*PI/3), r * 1.7f * sinf(2*PI/3), 0.0f, ATOM_H);
    addAtom(mol, r * 1.7f * cosf(4*PI/3), r * 1.7f * sinf(4*PI/3), 0.0f, ATOM_H);
    addAtom(mol, r * 1.7f * cosf(5*PI/3), r * 1.7f * sinf(5*PI/3), 0.0f, ATOM_H);

    addBond(mol, 0, 1, 2);
    addBond(mol, 1, 2, 1);
    addBond(mol, 2, 3, 2);
    addBond(mol, 3, 4, 1);
    addBond(mol, 4, 5, 2);
    addBond(mol, 5, 0, 1);
    addBond(mol, 0, 6, 1);
    addBond(mol, 6, 7, 1);
    addBond(mol, 3, 8, 1);
    addBond(mol, 8, 9, 1);
    addBond(mol, 8, 10, 1);
    addBond(mol, 10, 11, 2);
    addBond(mol, 10, 12, 1);
    addBond(mol, 1, 13, 1);
    addBond(mol, 2, 14, 1);
    addBond(mol, 4, 15, 1);
    addBond(mol, 5, 16, 1);

    centerMolecule(mol);
}

// ============== STEROID HORMONES ==============

// Helper: Build the steroid core (gonane skeleton) - 4 fused rings (A, B, C, D)
// Returns indices: 0-5 = ring A, 6-9 shared with A = ring B, 10-12 shared = ring C, 13-16 = ring D
void buildSteroidCore(Molecule* mol) {
    // Ring A (cyclohexane) - carbons 0-5
    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_C);      // C0 (C1)
    addAtom(mol, 1.4f, 0.0f, 0.5f, ATOM_C);      // C1 (C2)
    addAtom(mol, 2.5f, 0.0f, -0.4f, ATOM_C);     // C2 (C3)
    addAtom(mol, 2.5f, 1.3f, -1.0f, ATOM_C);     // C3 (C4)
    addAtom(mol, 1.4f, 2.0f, -0.5f, ATOM_C);     // C4 (C5)
    addAtom(mol, 0.0f, 1.5f, 0.0f, ATOM_C);      // C5 (C10)

    // Ring B (cyclohexane) shares C4, C5 - carbons 6-9
    addAtom(mol, -1.2f, 2.2f, 0.5f, ATOM_C);     // C6 (C6)
    addAtom(mol, -2.4f, 1.5f, 0.0f, ATOM_C);     // C7 (C7)
    addAtom(mol, -2.4f, 0.0f, 0.5f, ATOM_C);     // C8 (C8)
    addAtom(mol, -1.2f, -0.5f, 0.0f, ATOM_C);    // C9 (C9)

    // Ring C (cyclohexane) shares C7, C8 - carbons 10-13
    addAtom(mol, -3.6f, -0.5f, 0.0f, ATOM_C);    // C10 (C11)
    addAtom(mol, -4.8f, 0.2f, 0.5f, ATOM_C);     // C11 (C12)
    addAtom(mol, -4.8f, 1.5f, 0.0f, ATOM_C);     // C12 (C13)
    addAtom(mol, -3.6f, 2.2f, 0.5f, ATOM_C);     // C13 (C14)

    // Ring D (cyclopentane) shares C12, C13 - carbons 14-16
    addAtom(mol, -5.2f, 2.8f, -0.5f, ATOM_C);    // C14 (C15)
    addAtom(mol, -6.0f, 2.0f, -1.2f, ATOM_C);    // C15 (C16)
    addAtom(mol, -5.8f, 0.6f, -0.7f, ATOM_C);    // C16 (C17)

    // Angular methyl at C10 (C18)
    addAtom(mol, -1.2f, 1.5f, 1.5f, ATOM_C);     // C17 (C19 - angular methyl)

    // Angular methyl at C13 (C19)
    addAtom(mol, -3.6f, 2.0f, 2.0f, ATOM_C);     // C18 (C18 - angular methyl)

    // Ring A bonds
    addBond(mol, 0, 1, 1); addBond(mol, 1, 2, 1); addBond(mol, 2, 3, 1);
    addBond(mol, 3, 4, 1); addBond(mol, 4, 5, 1); addBond(mol, 5, 0, 1);

    // Ring B bonds (shares 5-9 edge conceptually, connects to ring A)
    addBond(mol, 5, 6, 1); addBond(mol, 6, 7, 1); addBond(mol, 7, 8, 1);
    addBond(mol, 8, 9, 1); addBond(mol, 9, 0, 1);

    // Ring C bonds (shares 8-13 with ring B)
    addBond(mol, 8, 10, 1); addBond(mol, 10, 11, 1); addBond(mol, 11, 12, 1);
    addBond(mol, 12, 13, 1); addBond(mol, 13, 7, 1);

    // Ring D bonds (5-membered, shares 12-13 with ring C)
    addBond(mol, 12, 14, 1); addBond(mol, 14, 15, 1); addBond(mol, 15, 16, 1);
    addBond(mol, 16, 11, 1);

    // Angular methyl bonds
    addBond(mol, 9, 17, 1);  // C19 methyl at junction
    addBond(mol, 13, 18, 1); // C18 methyl at junction
}

void buildTestosterone(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Testosterone (C19H28O2)");

    buildSteroidCore(mol);  // 19 carbons (0-18)

    // C3 ketone (=O on C2, index 2)
    addAtom(mol, 2.5f, -1.0f, -1.0f, ATOM_O);   // O at C3 (index 19)
    addBond(mol, 2, 19, 2);  // C=O double bond

    // C4-C5 double bond (indices 3-4)
    // Change bond 3-4 to double
    mol->bonds[3].order = 2;

    // C17 hydroxyl (on C16, index 16)
    addAtom(mol, -6.8f, 0.0f, -1.2f, ATOM_O);   // OH at C17 (index 20)
    addBond(mol, 16, 20, 1);

    // Add hydrogens to reach C19H28O2 (need 28 H)
    // Simplified - add key hydrogens
    int hStart = mol->numAtoms;
    addAtom(mol, -7.5f, 0.3f, -0.8f, ATOM_H);   // H on OH
    addBond(mol, 20, hStart, 1);

    // Add more H atoms at various positions
    float hPositions[][3] = {
        {0.0f, -0.8f, 0.7f}, {1.4f, -0.8f, 1.2f}, {3.3f, -0.5f, 0.0f},
        {1.4f, 3.0f, -0.8f}, {-1.2f, 3.2f, 0.8f}, {-0.8f, 2.2f, 1.3f},
        {-3.0f, -0.5f, -0.8f}, {-4.0f, -1.5f, 0.3f}, {-5.5f, -0.3f, 1.2f},
        {-5.5f, 2.0f, 0.7f}, {-4.8f, 3.5f, -1.0f}, {-6.8f, 2.5f, -1.8f},
        {-1.0f, 0.8f, 2.2f}, {-1.8f, 2.2f, 1.8f}, {-0.5f, 1.8f, 1.8f},
        {-3.2f, 1.3f, 2.7f}, {-4.3f, 2.7f, 2.3f}, {-3.0f, 2.8f, 2.3f},
        {-1.8f, -1.3f, -0.5f}, {0.3f, 2.0f, -0.8f}, {-2.8f, 1.8f, -0.8f},
        {-5.0f, 3.3f, 0.2f}, {-4.2f, -0.8f, 0.8f}, {-6.5f, 0.3f, 0.0f},
        {3.2f, 1.8f, -1.5f}, {2.0f, 1.0f, 1.0f}, {-0.7f, -0.3f, -0.7f}
    };
    for (int i = 0; i < 27 && mol->numAtoms < MAX_ATOMS; i++) {
        int idx = mol->numAtoms;
        addAtom(mol, hPositions[i][0], hPositions[i][1], hPositions[i][2], ATOM_H);
        addBond(mol, i % 19, idx, 1);
    }

    centerMolecule(mol);
}

void buildDHT(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "DHT (C19H30O2)");

    buildSteroidCore(mol);  // 19 carbons

    // C3 ketone
    addAtom(mol, 2.5f, -1.0f, -1.0f, ATOM_O);
    addBond(mol, 2, 19, 2);

    // C17 hydroxyl (no C4-C5 double bond in DHT, unlike testosterone)
    addAtom(mol, -6.8f, 0.0f, -1.2f, ATOM_O);
    addBond(mol, 16, 20, 1);

    // Add hydroxyl H
    addAtom(mol, -7.5f, 0.3f, -0.8f, ATOM_H);
    addBond(mol, 20, 21, 1);

    // Add remaining hydrogens (30 total - 1 on OH = 29 more needed, simplified)
    float hPos[][3] = {
        {0.0f, -0.8f, 0.7f}, {1.4f, -0.8f, 1.2f}, {3.3f, -0.5f, 0.0f},
        {3.2f, 1.8f, -1.5f}, {1.4f, 3.0f, -0.8f}, {-1.2f, 3.2f, 0.8f},
        {-3.0f, -0.5f, -0.8f}, {-4.0f, -1.5f, 0.3f}, {-5.5f, 2.0f, 0.7f},
        {-4.8f, 3.5f, -1.0f}, {-6.8f, 2.5f, -1.8f}, {-1.0f, 0.8f, 2.2f},
        {-3.2f, 1.3f, 2.7f}, {-4.3f, 2.7f, 2.3f}, {-1.8f, -1.3f, -0.5f},
        {0.3f, 2.0f, -0.8f}, {-2.8f, 1.8f, -0.8f}, {-5.0f, 3.3f, 0.2f},
        {1.8f, 1.5f, 0.2f}, {2.8f, 0.8f, -1.2f}, {-0.7f, 0.8f, -0.7f},
        {-1.8f, 2.8f, 0.0f}, {-2.0f, 0.5f, 0.8f}, {-4.2f, 0.5f, -0.7f},
        {-5.5f, 1.0f, 1.2f}, {-6.2f, 1.5f, -0.2f}, {-5.2f, 2.2f, -1.2f},
        {0.8f, 0.5f, 0.5f}, {-0.5f, -0.5f, 0.5f}
    };
    for (int i = 0; i < 29 && mol->numAtoms < MAX_ATOMS; i++) {
        int idx = mol->numAtoms;
        addAtom(mol, hPos[i][0], hPos[i][1], hPos[i][2], ATOM_H);
        addBond(mol, i % 19, idx, 1);
    }

    centerMolecule(mol);
}

void buildAndrostenedione(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Androstenedione (C19H26O2)");

    buildSteroidCore(mol);

    // C3 ketone
    addAtom(mol, 2.5f, -1.0f, -1.0f, ATOM_O);
    addBond(mol, 2, 19, 2);

    // C4-C5 double bond
    mol->bonds[3].order = 2;

    // C17 ketone (instead of hydroxyl)
    addAtom(mol, -6.8f, 0.0f, -1.2f, ATOM_O);
    addBond(mol, 16, 20, 2);  // Double bond for ketone

    // Add hydrogens (26 total)
    float hPos[][3] = {
        {0.0f, -0.8f, 0.7f}, {1.4f, -0.8f, 1.2f}, {3.3f, -0.5f, 0.0f},
        {1.4f, 3.0f, -0.8f}, {-1.2f, 3.2f, 0.8f}, {-3.0f, -0.5f, -0.8f},
        {-4.0f, -1.5f, 0.3f}, {-5.5f, 2.0f, 0.7f}, {-4.8f, 3.5f, -1.0f},
        {-6.8f, 2.5f, -1.8f}, {-1.0f, 0.8f, 2.2f}, {-3.2f, 1.3f, 2.7f},
        {-4.3f, 2.7f, 2.3f}, {-1.8f, -1.3f, -0.5f}, {0.3f, 2.0f, -0.8f},
        {-2.8f, 1.8f, -0.8f}, {-5.0f, 3.3f, 0.2f}, {-0.7f, 0.8f, -0.7f},
        {-1.8f, 2.8f, 0.0f}, {-2.0f, 0.5f, 0.8f}, {-4.2f, 0.5f, -0.7f},
        {-5.5f, 1.0f, 1.2f}, {-5.2f, 2.2f, -1.2f}, {0.8f, 0.5f, 0.5f},
        {3.2f, 1.8f, -1.5f}, {-0.5f, -0.5f, 0.5f}
    };
    for (int i = 0; i < 26 && mol->numAtoms < MAX_ATOMS; i++) {
        int idx = mol->numAtoms;
        addAtom(mol, hPos[i][0], hPos[i][1], hPos[i][2], ATOM_H);
        addBond(mol, i % 19, idx, 1);
    }

    centerMolecule(mol);
}

void buildEstradiol(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Estradiol/E2 (C18H24O2)");

    // Estrogens have aromatic A ring and no C19 methyl
    // Ring A (aromatic benzene-like) - carbons 0-5
    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_C);
    addAtom(mol, 1.4f, 0.0f, 0.0f, ATOM_C);
    addAtom(mol, 2.1f, 1.2f, 0.0f, ATOM_C);
    addAtom(mol, 1.4f, 2.4f, 0.0f, ATOM_C);
    addAtom(mol, 0.0f, 2.4f, 0.0f, ATOM_C);
    addAtom(mol, -0.7f, 1.2f, 0.0f, ATOM_C);

    // Ring B - carbons 6-9
    addAtom(mol, -2.1f, 1.2f, 0.5f, ATOM_C);
    addAtom(mol, -2.8f, 0.0f, 0.0f, ATOM_C);
    addAtom(mol, -2.1f, -1.2f, 0.5f, ATOM_C);
    addAtom(mol, -0.7f, -1.2f, 0.0f, ATOM_C);

    // Ring C - carbons 10-13
    addAtom(mol, -2.8f, -2.4f, 0.0f, ATOM_C);
    addAtom(mol, -4.2f, -2.4f, 0.5f, ATOM_C);
    addAtom(mol, -4.9f, -1.2f, 0.0f, ATOM_C);
    addAtom(mol, -4.2f, 0.0f, 0.5f, ATOM_C);

    // Ring D (5-membered) - carbons 14-16
    addAtom(mol, -5.6f, 0.0f, -0.5f, ATOM_C);
    addAtom(mol, -6.3f, -1.2f, 0.0f, ATOM_C);
    addAtom(mol, -5.6f, -2.4f, -0.5f, ATOM_C);

    // C18 angular methyl (only one in estrogens)
    addAtom(mol, -4.2f, 0.0f, 2.0f, ATOM_C);

    // Ring A bonds (aromatic)
    addBond(mol, 0, 1, 2); addBond(mol, 1, 2, 1); addBond(mol, 2, 3, 2);
    addBond(mol, 3, 4, 1); addBond(mol, 4, 5, 2); addBond(mol, 5, 0, 1);

    // Ring B bonds
    addBond(mol, 5, 6, 1); addBond(mol, 6, 7, 1); addBond(mol, 7, 8, 1);
    addBond(mol, 8, 9, 1); addBond(mol, 9, 0, 1);

    // Ring C bonds
    addBond(mol, 8, 10, 1); addBond(mol, 10, 11, 1); addBond(mol, 11, 12, 1);
    addBond(mol, 12, 13, 1); addBond(mol, 13, 7, 1);

    // Ring D bonds
    addBond(mol, 12, 14, 1); addBond(mol, 14, 15, 1); addBond(mol, 15, 16, 1);
    addBond(mol, 16, 11, 1);

    // Angular methyl bond
    addBond(mol, 13, 17, 1);

    // C3 hydroxyl (phenolic OH on aromatic ring)
    addAtom(mol, 1.4f, 3.6f, 0.0f, ATOM_O);
    addBond(mol, 3, 18, 1);

    // C17 hydroxyl
    addAtom(mol, -5.6f, 1.2f, -1.0f, ATOM_O);
    addBond(mol, 14, 19, 1);

    // Hydroxyl hydrogens
    addAtom(mol, 2.2f, 3.9f, 0.0f, ATOM_H);
    addBond(mol, 18, 20, 1);
    addAtom(mol, -6.3f, 1.5f, -0.5f, ATOM_H);
    addBond(mol, 19, 21, 1);

    // Add remaining hydrogens
    float hPos[][3] = {
        {2.0f, -0.9f, 0.0f}, {3.2f, 1.2f, 0.0f}, {-0.5f, 3.3f, 0.0f},
        {-2.1f, 1.2f, 1.6f}, {-2.6f, 2.1f, 0.2f}, {-2.1f, -1.2f, 1.6f},
        {-0.7f, -1.2f, -1.1f}, {-0.2f, -2.1f, 0.3f}, {-2.3f, -3.3f, 0.3f},
        {-4.2f, -2.4f, 1.6f}, {-4.7f, -3.3f, 0.2f}, {-6.0f, -0.5f, -1.3f},
        {-7.0f, -0.9f, 0.7f}, {-7.0f, -1.5f, -0.7f}, {-5.2f, -3.0f, 0.2f},
        {-6.2f, -2.9f, -1.0f}, {-3.5f, -0.5f, 2.5f}, {-5.0f, 0.5f, 2.3f},
        {-4.0f, 0.9f, 2.3f}, {-3.5f, -0.9f, -0.5f}, {-2.8f, -0.5f, -1.0f}
    };
    for (int i = 0; i < 21 && mol->numAtoms < MAX_ATOMS; i++) {
        int idx = mol->numAtoms;
        addAtom(mol, hPos[i][0], hPos[i][1], hPos[i][2], ATOM_H);
        addBond(mol, i % 18, idx, 1);
    }

    centerMolecule(mol);
}

void buildEstrone(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Estrone/E1 (C18H22O2)");

    // Similar to estradiol but C17 is ketone
    // Ring A (aromatic)
    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_C);
    addAtom(mol, 1.4f, 0.0f, 0.0f, ATOM_C);
    addAtom(mol, 2.1f, 1.2f, 0.0f, ATOM_C);
    addAtom(mol, 1.4f, 2.4f, 0.0f, ATOM_C);
    addAtom(mol, 0.0f, 2.4f, 0.0f, ATOM_C);
    addAtom(mol, -0.7f, 1.2f, 0.0f, ATOM_C);

    // Rings B, C, D
    addAtom(mol, -2.1f, 1.2f, 0.5f, ATOM_C);
    addAtom(mol, -2.8f, 0.0f, 0.0f, ATOM_C);
    addAtom(mol, -2.1f, -1.2f, 0.5f, ATOM_C);
    addAtom(mol, -0.7f, -1.2f, 0.0f, ATOM_C);
    addAtom(mol, -2.8f, -2.4f, 0.0f, ATOM_C);
    addAtom(mol, -4.2f, -2.4f, 0.5f, ATOM_C);
    addAtom(mol, -4.9f, -1.2f, 0.0f, ATOM_C);
    addAtom(mol, -4.2f, 0.0f, 0.5f, ATOM_C);
    addAtom(mol, -5.6f, 0.0f, -0.5f, ATOM_C);
    addAtom(mol, -6.3f, -1.2f, 0.0f, ATOM_C);
    addAtom(mol, -5.6f, -2.4f, -0.5f, ATOM_C);
    addAtom(mol, -4.2f, 0.0f, 2.0f, ATOM_C);  // Angular methyl

    // Aromatic ring A bonds
    addBond(mol, 0, 1, 2); addBond(mol, 1, 2, 1); addBond(mol, 2, 3, 2);
    addBond(mol, 3, 4, 1); addBond(mol, 4, 5, 2); addBond(mol, 5, 0, 1);

    // Other ring bonds
    addBond(mol, 5, 6, 1); addBond(mol, 6, 7, 1); addBond(mol, 7, 8, 1);
    addBond(mol, 8, 9, 1); addBond(mol, 9, 0, 1);
    addBond(mol, 8, 10, 1); addBond(mol, 10, 11, 1); addBond(mol, 11, 12, 1);
    addBond(mol, 12, 13, 1); addBond(mol, 13, 7, 1);
    addBond(mol, 12, 14, 1); addBond(mol, 14, 15, 1); addBond(mol, 15, 16, 1);
    addBond(mol, 16, 11, 1);
    addBond(mol, 13, 17, 1);

    // C3 hydroxyl
    addAtom(mol, 1.4f, 3.6f, 0.0f, ATOM_O);
    addBond(mol, 3, 18, 1);

    // C17 ketone (double bond O)
    addAtom(mol, -5.6f, 1.2f, -1.0f, ATOM_O);
    addBond(mol, 14, 19, 2);

    // Hydroxyl hydrogen
    addAtom(mol, 2.2f, 3.9f, 0.0f, ATOM_H);
    addBond(mol, 18, 20, 1);

    // Add remaining hydrogens (22 - 1 on phenol = 21)
    float hPos[][3] = {
        {2.0f, -0.9f, 0.0f}, {3.2f, 1.2f, 0.0f}, {-0.5f, 3.3f, 0.0f},
        {-2.1f, 1.2f, 1.6f}, {-2.6f, 2.1f, 0.2f}, {-2.1f, -1.2f, 1.6f},
        {-0.7f, -1.2f, -1.1f}, {-0.2f, -2.1f, 0.3f}, {-2.3f, -3.3f, 0.3f},
        {-4.2f, -2.4f, 1.6f}, {-4.7f, -3.3f, 0.2f}, {-7.0f, -0.9f, 0.7f},
        {-7.0f, -1.5f, -0.7f}, {-5.2f, -3.0f, 0.2f}, {-6.2f, -2.9f, -1.0f},
        {-3.5f, -0.5f, 2.5f}, {-5.0f, 0.5f, 2.3f}, {-4.0f, 0.9f, 2.3f},
        {-3.5f, -0.9f, -0.5f}, {-2.8f, -0.5f, -1.0f}, {-6.0f, -0.5f, -1.3f}
    };
    for (int i = 0; i < 21 && mol->numAtoms < MAX_ATOMS; i++) {
        int idx = mol->numAtoms;
        addAtom(mol, hPos[i][0], hPos[i][1], hPos[i][2], ATOM_H);
        addBond(mol, i % 18, idx, 1);
    }

    centerMolecule(mol);
}

void buildEstriol(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Estriol/E3 (C18H24O3)");

    // Same base as estradiol but with extra OH at C16
    // Ring A (aromatic)
    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_C);
    addAtom(mol, 1.4f, 0.0f, 0.0f, ATOM_C);
    addAtom(mol, 2.1f, 1.2f, 0.0f, ATOM_C);
    addAtom(mol, 1.4f, 2.4f, 0.0f, ATOM_C);
    addAtom(mol, 0.0f, 2.4f, 0.0f, ATOM_C);
    addAtom(mol, -0.7f, 1.2f, 0.0f, ATOM_C);

    // Rings B, C, D
    addAtom(mol, -2.1f, 1.2f, 0.5f, ATOM_C);
    addAtom(mol, -2.8f, 0.0f, 0.0f, ATOM_C);
    addAtom(mol, -2.1f, -1.2f, 0.5f, ATOM_C);
    addAtom(mol, -0.7f, -1.2f, 0.0f, ATOM_C);
    addAtom(mol, -2.8f, -2.4f, 0.0f, ATOM_C);
    addAtom(mol, -4.2f, -2.4f, 0.5f, ATOM_C);
    addAtom(mol, -4.9f, -1.2f, 0.0f, ATOM_C);
    addAtom(mol, -4.2f, 0.0f, 0.5f, ATOM_C);
    addAtom(mol, -5.6f, 0.0f, -0.5f, ATOM_C);
    addAtom(mol, -6.3f, -1.2f, 0.0f, ATOM_C);
    addAtom(mol, -5.6f, -2.4f, -0.5f, ATOM_C);
    addAtom(mol, -4.2f, 0.0f, 2.0f, ATOM_C);  // Angular methyl

    // Ring bonds
    addBond(mol, 0, 1, 2); addBond(mol, 1, 2, 1); addBond(mol, 2, 3, 2);
    addBond(mol, 3, 4, 1); addBond(mol, 4, 5, 2); addBond(mol, 5, 0, 1);
    addBond(mol, 5, 6, 1); addBond(mol, 6, 7, 1); addBond(mol, 7, 8, 1);
    addBond(mol, 8, 9, 1); addBond(mol, 9, 0, 1);
    addBond(mol, 8, 10, 1); addBond(mol, 10, 11, 1); addBond(mol, 11, 12, 1);
    addBond(mol, 12, 13, 1); addBond(mol, 13, 7, 1);
    addBond(mol, 12, 14, 1); addBond(mol, 14, 15, 1); addBond(mol, 15, 16, 1);
    addBond(mol, 16, 11, 1);
    addBond(mol, 13, 17, 1);

    // C3 hydroxyl (phenolic)
    addAtom(mol, 1.4f, 3.6f, 0.0f, ATOM_O);
    addBond(mol, 3, 18, 1);

    // C16 hydroxyl
    addAtom(mol, -7.0f, -1.2f, 1.0f, ATOM_O);
    addBond(mol, 15, 19, 1);

    // C17 hydroxyl
    addAtom(mol, -5.6f, 1.2f, -1.0f, ATOM_O);
    addBond(mol, 14, 20, 1);

    // Hydroxyl hydrogens
    addAtom(mol, 2.2f, 3.9f, 0.0f, ATOM_H);
    addBond(mol, 18, 21, 1);
    addAtom(mol, -7.7f, -0.8f, 0.6f, ATOM_H);
    addBond(mol, 19, 22, 1);
    addAtom(mol, -6.3f, 1.5f, -0.5f, ATOM_H);
    addBond(mol, 20, 23, 1);

    // Add remaining hydrogens
    float hPos[][3] = {
        {2.0f, -0.9f, 0.0f}, {3.2f, 1.2f, 0.0f}, {-0.5f, 3.3f, 0.0f},
        {-2.1f, 1.2f, 1.6f}, {-2.6f, 2.1f, 0.2f}, {-2.1f, -1.2f, 1.6f},
        {-0.7f, -1.2f, -1.1f}, {-0.2f, -2.1f, 0.3f}, {-2.3f, -3.3f, 0.3f},
        {-4.2f, -2.4f, 1.6f}, {-4.7f, -3.3f, 0.2f}, {-6.0f, -0.5f, -1.3f},
        {-5.2f, -3.0f, 0.2f}, {-6.2f, -2.9f, -1.0f}, {-3.5f, -0.5f, 2.5f},
        {-5.0f, 0.5f, 2.3f}, {-4.0f, 0.9f, 2.3f}, {-3.5f, -0.9f, -0.5f},
        {-2.8f, -0.5f, -1.0f}, {-6.5f, -1.8f, -0.7f}, {-5.8f, -2.8f, 0.3f}
    };
    for (int i = 0; i < 21 && mol->numAtoms < MAX_ATOMS; i++) {
        int idx = mol->numAtoms;
        addAtom(mol, hPos[i][0], hPos[i][1], hPos[i][2], ATOM_H);
        addBond(mol, i % 18, idx, 1);
    }

    centerMolecule(mol);
}

void buildProgesterone(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Progesterone (C21H30O2)");

    buildSteroidCore(mol);  // 19 carbons (0-18)

    // C3 ketone
    addAtom(mol, 2.5f, -1.0f, -1.0f, ATOM_O);
    addBond(mol, 2, 19, 2);

    // C4-C5 double bond
    mol->bonds[3].order = 2;

    // C17 acetyl group (-COCH3)
    addAtom(mol, -6.8f, 0.0f, -1.2f, ATOM_C);   // Carbonyl carbon (C20)
    addBond(mol, 16, 20, 1);
    addAtom(mol, -7.5f, 0.5f, -2.0f, ATOM_O);   // Carbonyl oxygen
    addBond(mol, 20, 21, 2);
    addAtom(mol, -7.5f, -1.0f, -0.5f, ATOM_C);  // Methyl (C21)
    addBond(mol, 20, 22, 1);

    // Methyl hydrogens on C21
    addAtom(mol, -8.3f, -0.5f, 0.0f, ATOM_H);
    addAtom(mol, -7.9f, -1.6f, -1.2f, ATOM_H);
    addAtom(mol, -6.9f, -1.6f, 0.1f, ATOM_H);
    addBond(mol, 22, 23, 1);
    addBond(mol, 22, 24, 1);
    addBond(mol, 22, 25, 1);

    // Add remaining hydrogens
    float hPos[][3] = {
        {0.0f, -0.8f, 0.7f}, {1.4f, -0.8f, 1.2f}, {3.3f, -0.5f, 0.0f},
        {1.4f, 3.0f, -0.8f}, {-1.2f, 3.2f, 0.8f}, {-3.0f, -0.5f, -0.8f},
        {-4.0f, -1.5f, 0.3f}, {-5.5f, 2.0f, 0.7f}, {-4.8f, 3.5f, -1.0f},
        {-6.8f, 2.5f, -1.8f}, {-1.0f, 0.8f, 2.2f}, {-3.2f, 1.3f, 2.7f},
        {-4.3f, 2.7f, 2.3f}, {-1.8f, -1.3f, -0.5f}, {0.3f, 2.0f, -0.8f},
        {-2.8f, 1.8f, -0.8f}, {-5.0f, 3.3f, 0.2f}, {-0.7f, 0.8f, -0.7f},
        {-1.8f, 2.8f, 0.0f}, {-2.0f, 0.5f, 0.8f}, {-4.2f, 0.5f, -0.7f},
        {-5.5f, 1.0f, 1.2f}, {0.8f, 0.5f, 0.5f}, {3.2f, 1.8f, -1.5f}
    };
    for (int i = 0; i < 24 && mol->numAtoms < MAX_ATOMS; i++) {
        int idx = mol->numAtoms;
        addAtom(mol, hPos[i][0], hPos[i][1], hPos[i][2], ATOM_H);
        addBond(mol, i % 19, idx, 1);
    }

    centerMolecule(mol);
}

void buildCortisol(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Cortisol (C21H30O5)");

    buildSteroidCore(mol);  // 19 carbons (0-18)

    // C3 ketone
    addAtom(mol, 2.5f, -1.0f, -1.0f, ATOM_O);
    addBond(mol, 2, 19, 2);

    // C4-C5 double bond
    mol->bonds[3].order = 2;

    // C11 hydroxyl
    addAtom(mol, -4.8f, -0.5f, 1.5f, ATOM_O);
    addBond(mol, 10, 20, 1);

    // C17 hydroxyl
    addAtom(mol, -5.8f, 1.0f, -1.5f, ATOM_O);
    addBond(mol, 16, 21, 1);

    // C17 side chain: -COCH2OH
    addAtom(mol, -6.8f, 0.0f, -0.5f, ATOM_C);   // Carbonyl carbon (C20)
    addBond(mol, 16, 22, 1);
    addAtom(mol, -7.5f, 0.5f, 0.3f, ATOM_O);    // Carbonyl oxygen
    addBond(mol, 22, 23, 2);
    addAtom(mol, -7.2f, -1.2f, -1.0f, ATOM_C);  // CH2 (C21)
    addBond(mol, 22, 24, 1);
    addAtom(mol, -8.3f, -1.5f, -0.3f, ATOM_O);  // Primary OH
    addBond(mol, 24, 25, 1);

    // Hydroxyl hydrogens
    addAtom(mol, -5.5f, 0.0f, 2.0f, ATOM_H);
    addBond(mol, 20, 26, 1);
    addAtom(mol, -5.3f, 1.5f, -2.0f, ATOM_H);
    addBond(mol, 21, 27, 1);
    addAtom(mol, -8.8f, -0.8f, 0.0f, ATOM_H);
    addBond(mol, 25, 28, 1);

    // CH2 hydrogens
    addAtom(mol, -6.5f, -1.8f, -1.5f, ATOM_H);
    addAtom(mol, -7.5f, -1.0f, -1.9f, ATOM_H);
    addBond(mol, 24, 29, 1);
    addBond(mol, 24, 30, 1);

    // Add remaining hydrogens
    float hPos[][3] = {
        {0.0f, -0.8f, 0.7f}, {1.4f, -0.8f, 1.2f}, {3.3f, -0.5f, 0.0f},
        {1.4f, 3.0f, -0.8f}, {-1.2f, 3.2f, 0.8f}, {-3.0f, -0.5f, -0.8f},
        {-5.5f, 2.0f, 0.7f}, {-4.8f, 3.5f, -1.0f}, {-6.8f, 2.5f, -1.8f},
        {-1.0f, 0.8f, 2.2f}, {-3.2f, 1.3f, 2.7f}, {-4.3f, 2.7f, 2.3f},
        {-1.8f, -1.3f, -0.5f}, {0.3f, 2.0f, -0.8f}, {-2.8f, 1.8f, -0.8f},
        {-5.0f, 3.3f, 0.2f}, {0.8f, 0.5f, 0.5f}, {3.2f, 1.8f, -1.5f},
        {-4.0f, -1.5f, 0.3f}
    };
    for (int i = 0; i < 19 && mol->numAtoms < MAX_ATOMS; i++) {
        int idx = mol->numAtoms;
        addAtom(mol, hPos[i][0], hPos[i][1], hPos[i][2], ATOM_H);
        addBond(mol, i % 19, idx, 1);
    }

    centerMolecule(mol);
}

void buildCortisone(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Cortisone (C21H28O5)");

    buildSteroidCore(mol);

    // C3 ketone
    addAtom(mol, 2.5f, -1.0f, -1.0f, ATOM_O);
    addBond(mol, 2, 19, 2);

    // C4-C5 double bond
    mol->bonds[3].order = 2;

    // C11 ketone (not hydroxyl like cortisol)
    addAtom(mol, -4.8f, -0.5f, 1.5f, ATOM_O);
    addBond(mol, 10, 20, 2);

    // C17 hydroxyl
    addAtom(mol, -5.8f, 1.0f, -1.5f, ATOM_O);
    addBond(mol, 16, 21, 1);

    // C17 side chain
    addAtom(mol, -6.8f, 0.0f, -0.5f, ATOM_C);
    addBond(mol, 16, 22, 1);
    addAtom(mol, -7.5f, 0.5f, 0.3f, ATOM_O);
    addBond(mol, 22, 23, 2);
    addAtom(mol, -7.2f, -1.2f, -1.0f, ATOM_C);
    addBond(mol, 22, 24, 1);
    addAtom(mol, -8.3f, -1.5f, -0.3f, ATOM_O);
    addBond(mol, 24, 25, 1);

    // Hydroxyl hydrogens
    addAtom(mol, -5.3f, 1.5f, -2.0f, ATOM_H);
    addBond(mol, 21, 26, 1);
    addAtom(mol, -8.8f, -0.8f, 0.0f, ATOM_H);
    addBond(mol, 25, 27, 1);

    // CH2 hydrogens
    addAtom(mol, -6.5f, -1.8f, -1.5f, ATOM_H);
    addAtom(mol, -7.5f, -1.0f, -1.9f, ATOM_H);
    addBond(mol, 24, 28, 1);
    addBond(mol, 24, 29, 1);

    // Add remaining hydrogens (28 - 2 OH H - 2 CH2 H = 24 more)
    float hPos[][3] = {
        {0.0f, -0.8f, 0.7f}, {1.4f, -0.8f, 1.2f}, {3.3f, -0.5f, 0.0f},
        {1.4f, 3.0f, -0.8f}, {-1.2f, 3.2f, 0.8f}, {-3.0f, -0.5f, -0.8f},
        {-5.5f, 2.0f, 0.7f}, {-4.8f, 3.5f, -1.0f}, {-6.8f, 2.5f, -1.8f},
        {-1.0f, 0.8f, 2.2f}, {-3.2f, 1.3f, 2.7f}, {-4.3f, 2.7f, 2.3f},
        {-1.8f, -1.3f, -0.5f}, {0.3f, 2.0f, -0.8f}, {-2.8f, 1.8f, -0.8f},
        {-5.0f, 3.3f, 0.2f}, {0.8f, 0.5f, 0.5f}, {3.2f, 1.8f, -1.5f},
        {-4.0f, -1.5f, 0.3f}, {-2.5f, -0.8f, 0.5f}, {-0.5f, 0.5f, -0.8f},
        {-3.8f, 0.5f, 0.0f}, {-5.2f, -0.3f, -0.5f}, {-6.2f, 1.8f, 0.0f}
    };
    for (int i = 0; i < 24 && mol->numAtoms < MAX_ATOMS; i++) {
        int idx = mol->numAtoms;
        addAtom(mol, hPos[i][0], hPos[i][1], hPos[i][2], ATOM_H);
        addBond(mol, i % 19, idx, 1);
    }

    centerMolecule(mol);
}

void buildAldosterone(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Aldosterone (C21H28O5)");

    buildSteroidCore(mol);

    // C3 ketone
    addAtom(mol, 2.5f, -1.0f, -1.0f, ATOM_O);
    addBond(mol, 2, 19, 2);

    // C4-C5 double bond
    mol->bonds[3].order = 2;

    // C11 hydroxyl
    addAtom(mol, -4.8f, -0.5f, 1.5f, ATOM_O);
    addBond(mol, 10, 20, 1);

    // C18 aldehyde (CHO replaces angular methyl at C13)
    // Remove the existing methyl hydrogen assumption, add aldehyde
    addAtom(mol, -3.6f, 3.0f, 2.5f, ATOM_O);  // Aldehyde oxygen
    addBond(mol, 18, 21, 2);

    // C17 side chain
    addAtom(mol, -6.8f, 0.0f, -0.5f, ATOM_C);
    addBond(mol, 16, 22, 1);
    addAtom(mol, -7.5f, 0.5f, 0.3f, ATOM_O);
    addBond(mol, 22, 23, 2);
    addAtom(mol, -7.2f, -1.2f, -1.0f, ATOM_C);
    addBond(mol, 22, 24, 1);
    addAtom(mol, -8.3f, -1.5f, -0.3f, ATOM_O);
    addBond(mol, 24, 25, 1);

    // Hydroxyl hydrogens
    addAtom(mol, -5.5f, 0.0f, 2.0f, ATOM_H);
    addBond(mol, 20, 26, 1);
    addAtom(mol, -8.8f, -0.8f, 0.0f, ATOM_H);
    addBond(mol, 25, 27, 1);

    // Aldehyde H
    addAtom(mol, -3.0f, 3.5f, 1.8f, ATOM_H);
    addBond(mol, 18, 28, 1);

    // CH2 hydrogens
    addAtom(mol, -6.5f, -1.8f, -1.5f, ATOM_H);
    addAtom(mol, -7.5f, -1.0f, -1.9f, ATOM_H);
    addBond(mol, 24, 29, 1);
    addBond(mol, 24, 30, 1);

    // Add remaining hydrogens
    float hPos[][3] = {
        {0.0f, -0.8f, 0.7f}, {1.4f, -0.8f, 1.2f}, {3.3f, -0.5f, 0.0f},
        {1.4f, 3.0f, -0.8f}, {-1.2f, 3.2f, 0.8f}, {-3.0f, -0.5f, -0.8f},
        {-5.5f, 2.0f, 0.7f}, {-4.8f, 3.5f, -1.0f}, {-6.8f, 2.5f, -1.8f},
        {-1.0f, 0.8f, 2.2f}, {-1.8f, -1.3f, -0.5f}, {0.3f, 2.0f, -0.8f},
        {-2.8f, 1.8f, -0.8f}, {-5.0f, 3.3f, 0.2f}, {0.8f, 0.5f, 0.5f},
        {3.2f, 1.8f, -1.5f}, {-4.0f, -1.5f, 0.3f}, {-2.5f, -0.8f, 0.5f},
        {-0.5f, 0.5f, -0.8f}, {-3.8f, 0.5f, 0.0f}, {-5.2f, -0.3f, -0.5f}
    };
    for (int i = 0; i < 21 && mol->numAtoms < MAX_ATOMS; i++) {
        int idx = mol->numAtoms;
        addAtom(mol, hPos[i][0], hPos[i][1], hPos[i][2], ATOM_H);
        addBond(mol, i % 19, idx, 1);
    }

    centerMolecule(mol);
}

// Build Acetylcholine (C7H16NO2+) - Memory & muscle control neurotransmitter
void buildAcetylcholine(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Acetylcholine (C7H16NO2)");

    // Choline moiety: (CH3)3N+-CH2-CH2-O-
    // Quaternary ammonium nitrogen
    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_N);         // 0: N+

    // Three methyl groups on nitrogen
    addAtom(mol, 1.0f, 0.8f, 0.5f, ATOM_C);         // 1: CH3
    addAtom(mol, 0.8f, -1.0f, -0.3f, ATOM_C);       // 2: CH3
    addAtom(mol, -0.5f, 0.5f, -1.2f, ATOM_C);       // 3: CH3

    // Ethyl bridge to ester
    addAtom(mol, -1.3f, -0.5f, 0.7f, ATOM_C);       // 4: CH2
    addAtom(mol, -2.5f, 0.2f, 1.3f, ATOM_C);        // 5: CH2

    // Ester oxygen
    addAtom(mol, -3.6f, -0.5f, 0.8f, ATOM_O);       // 6: O (ester)

    // Acetyl group: CH3-C(=O)-
    addAtom(mol, -4.8f, 0.2f, 1.2f, ATOM_C);        // 7: C=O
    addAtom(mol, -5.0f, 1.3f, 1.8f, ATOM_O);        // 8: =O
    addAtom(mol, -6.0f, -0.5f, 0.8f, ATOM_C);       // 9: CH3

    // Hydrogens on methyl groups attached to N
    addAtom(mol, 1.8f, 0.3f, 0.9f, ATOM_H);         // 10
    addAtom(mol, 0.7f, 1.6f, 1.0f, ATOM_H);         // 11
    addAtom(mol, 1.4f, 1.2f, -0.3f, ATOM_H);        // 12

    addAtom(mol, 1.6f, -1.0f, 0.3f, ATOM_H);        // 13
    addAtom(mol, 1.0f, -1.0f, -1.3f, ATOM_H);       // 14
    addAtom(mol, 0.3f, -1.9f, -0.1f, ATOM_H);       // 15

    addAtom(mol, -0.2f, 1.4f, -1.7f, ATOM_H);       // 16
    addAtom(mol, -1.4f, 0.7f, -0.9f, ATOM_H);       // 17
    addAtom(mol, -0.7f, -0.2f, -1.9f, ATOM_H);      // 18

    // Hydrogens on CH2 groups
    addAtom(mol, -1.5f, -1.4f, 0.2f, ATOM_H);       // 19
    addAtom(mol, -1.0f, -0.9f, 1.6f, ATOM_H);       // 20
    addAtom(mol, -2.3f, 1.1f, 1.8f, ATOM_H);        // 21
    addAtom(mol, -2.8f, 0.6f, 0.4f, ATOM_H);        // 22

    // Hydrogens on terminal CH3
    addAtom(mol, -6.8f, 0.1f, 1.2f, ATOM_H);        // 23
    addAtom(mol, -6.2f, -0.6f, -0.2f, ATOM_H);      // 24
    addAtom(mol, -5.9f, -1.4f, 1.3f, ATOM_H);       // 25

    // Bonds
    addBond(mol, 0, 1, 1);  // N-CH3
    addBond(mol, 0, 2, 1);  // N-CH3
    addBond(mol, 0, 3, 1);  // N-CH3
    addBond(mol, 0, 4, 1);  // N-CH2
    addBond(mol, 4, 5, 1);  // CH2-CH2
    addBond(mol, 5, 6, 1);  // CH2-O
    addBond(mol, 6, 7, 1);  // O-C=O
    addBond(mol, 7, 8, 2);  // C=O
    addBond(mol, 7, 9, 1);  // C-CH3

    // H bonds for methyl groups on N
    addBond(mol, 1, 10, 1);
    addBond(mol, 1, 11, 1);
    addBond(mol, 1, 12, 1);
    addBond(mol, 2, 13, 1);
    addBond(mol, 2, 14, 1);
    addBond(mol, 2, 15, 1);
    addBond(mol, 3, 16, 1);
    addBond(mol, 3, 17, 1);
    addBond(mol, 3, 18, 1);

    // H bonds for CH2 groups
    addBond(mol, 4, 19, 1);
    addBond(mol, 4, 20, 1);
    addBond(mol, 5, 21, 1);
    addBond(mol, 5, 22, 1);

    // H bonds for terminal CH3
    addBond(mol, 9, 23, 1);
    addBond(mol, 9, 24, 1);
    addBond(mol, 9, 25, 1);

    centerMolecule(mol);
}

// Build GABA (C4H9NO2) - Main inhibitory neurotransmitter
void buildGABA(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "GABA (C4H9NO2)");

    // Gamma-aminobutyric acid: H2N-CH2-CH2-CH2-COOH
    // Amino group
    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_N);         // 0: NH2

    // Carbon chain
    addAtom(mol, 1.4f, 0.3f, 0.5f, ATOM_C);         // 1: CH2 (gamma)
    addAtom(mol, 2.5f, -0.5f, -0.2f, ATOM_C);       // 2: CH2 (beta)
    addAtom(mol, 3.9f, -0.2f, 0.3f, ATOM_C);        // 3: CH2 (alpha)

    // Carboxyl group
    addAtom(mol, 5.0f, -0.9f, -0.4f, ATOM_C);       // 4: COOH
    addAtom(mol, 5.0f, -2.0f, -1.0f, ATOM_O);       // 5: =O
    addAtom(mol, 6.1f, -0.2f, -0.2f, ATOM_O);       // 6: OH

    // Hydrogens on amino group
    addAtom(mol, -0.5f, 0.8f, 0.4f, ATOM_H);        // 7
    addAtom(mol, -0.5f, -0.8f, 0.3f, ATOM_H);       // 8

    // Hydrogens on CH2 groups
    addAtom(mol, 1.5f, 1.3f, 0.3f, ATOM_H);         // 9
    addAtom(mol, 1.5f, 0.2f, 1.5f, ATOM_H);         // 10
    addAtom(mol, 2.4f, -1.5f, 0.0f, ATOM_H);        // 11
    addAtom(mol, 2.4f, -0.4f, -1.2f, ATOM_H);       // 12
    addAtom(mol, 4.0f, 0.8f, 0.1f, ATOM_H);         // 13
    addAtom(mol, 4.0f, -0.3f, 1.3f, ATOM_H);        // 14

    // Hydrogen on COOH
    addAtom(mol, 6.8f, -0.7f, -0.6f, ATOM_H);       // 15

    // Bonds
    addBond(mol, 0, 1, 1);  // N-CH2
    addBond(mol, 1, 2, 1);  // CH2-CH2
    addBond(mol, 2, 3, 1);  // CH2-CH2
    addBond(mol, 3, 4, 1);  // CH2-COOH
    addBond(mol, 4, 5, 2);  // C=O
    addBond(mol, 4, 6, 1);  // C-OH

    // H bonds
    addBond(mol, 0, 7, 1);
    addBond(mol, 0, 8, 1);
    addBond(mol, 1, 9, 1);
    addBond(mol, 1, 10, 1);
    addBond(mol, 2, 11, 1);
    addBond(mol, 2, 12, 1);
    addBond(mol, 3, 13, 1);
    addBond(mol, 3, 14, 1);
    addBond(mol, 6, 15, 1);

    centerMolecule(mol);
}

// Build Glutamate (C5H9NO4) - Main excitatory neurotransmitter
void buildGlutamate(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Glutamate (C5H9NO4)");

    // Glutamic acid: HOOC-CH2-CH2-CH(NH2)-COOH
    // Alpha carbon with amino group
    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_C);         // 0: CH (alpha)
    addAtom(mol, 0.5f, 1.3f, 0.3f, ATOM_N);         // 1: NH2

    // Alpha carboxyl group
    addAtom(mol, 1.3f, -0.7f, -0.5f, ATOM_C);       // 2: COOH
    addAtom(mol, 1.5f, -1.8f, -0.1f, ATOM_O);       // 3: =O
    addAtom(mol, 2.2f, -0.1f, -1.3f, ATOM_O);       // 4: OH

    // Side chain CH2-CH2-COOH
    addAtom(mol, -1.2f, -0.3f, 0.8f, ATOM_C);       // 5: CH2 (beta)
    addAtom(mol, -2.4f, 0.4f, 0.3f, ATOM_C);        // 6: CH2 (gamma)

    // Gamma carboxyl group
    addAtom(mol, -3.6f, 0.0f, 1.0f, ATOM_C);        // 7: COOH
    addAtom(mol, -3.8f, -1.0f, 1.7f, ATOM_O);       // 8: =O
    addAtom(mol, -4.6f, 0.8f, 0.8f, ATOM_O);        // 9: OH

    // Hydrogens
    addAtom(mol, -0.3f, 0.3f, -0.9f, ATOM_H);       // 10: H on alpha C
    addAtom(mol, 0.0f, 2.0f, -0.1f, ATOM_H);        // 11: H on NH2
    addAtom(mol, 1.4f, 1.5f, 0.0f, ATOM_H);         // 12: H on NH2
    addAtom(mol, 2.9f, -0.6f, -1.6f, ATOM_H);       // 13: H on alpha COOH
    addAtom(mol, -1.0f, 0.0f, 1.8f, ATOM_H);        // 14: H on beta CH2
    addAtom(mol, -1.4f, -1.3f, 0.9f, ATOM_H);       // 15: H on beta CH2
    addAtom(mol, -2.2f, 0.2f, -0.7f, ATOM_H);       // 16: H on gamma CH2
    addAtom(mol, -2.6f, 1.4f, 0.5f, ATOM_H);        // 17: H on gamma CH2
    addAtom(mol, -5.3f, 0.4f, 1.3f, ATOM_H);        // 18: H on gamma COOH

    // Bonds
    addBond(mol, 0, 1, 1);  // CH-NH2
    addBond(mol, 0, 2, 1);  // CH-COOH
    addBond(mol, 2, 3, 2);  // C=O
    addBond(mol, 2, 4, 1);  // C-OH
    addBond(mol, 0, 5, 1);  // CH-CH2
    addBond(mol, 5, 6, 1);  // CH2-CH2
    addBond(mol, 6, 7, 1);  // CH2-COOH
    addBond(mol, 7, 8, 2);  // C=O
    addBond(mol, 7, 9, 1);  // C-OH

    // H bonds
    addBond(mol, 0, 10, 1);
    addBond(mol, 1, 11, 1);
    addBond(mol, 1, 12, 1);
    addBond(mol, 4, 13, 1);
    addBond(mol, 5, 14, 1);
    addBond(mol, 5, 15, 1);
    addBond(mol, 6, 16, 1);
    addBond(mol, 6, 17, 1);
    addBond(mol, 9, 18, 1);

    centerMolecule(mol);
}

// Build Norepinephrine/Epinephrine (C8H11NO3) - Fight or flight hormone/neurotransmitter
void buildNorepinephrine(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Norepinephrine (C8H11NO3)");

    // Catechol ring (benzene with 2 OH groups)
    float r = 1.4f;
    for (int i = 0; i < 6; i++) {
        float angle = i * PI / 3.0f;
        addAtom(mol, r * cosf(angle), r * sinf(angle), 0.0f, ATOM_C);
    }

    // OH groups on C3 and C4 positions (meta and para to side chain)
    addAtom(mol, r * cosf(PI) - 0.9f, r * sinf(PI) + 0.5f, 0.0f, ATOM_O);      // 6: OH
    addAtom(mol, r * cosf(4*PI/3) - 0.5f, r * sinf(4*PI/3) - 0.9f, 0.0f, ATOM_O); // 7: OH

    // Side chain: -CH(OH)-CH2-NH2 from C1
    addAtom(mol, 2.5f, 0.3f, 0.0f, ATOM_C);         // 8: CH(OH)
    addAtom(mol, 3.0f, 0.8f, 1.2f, ATOM_O);         // 9: OH on side chain
    addAtom(mol, 3.5f, -0.5f, -0.8f, ATOM_C);       // 10: CH2
    addAtom(mol, 4.8f, 0.2f, -0.5f, ATOM_N);        // 11: NH2

    // Hydrogens on ring (C2, C5, C6 positions)
    float rH = 2.4f;
    addAtom(mol, rH * cosf(PI/3), rH * sinf(PI/3), 0.0f, ATOM_H);     // 12: H on C2
    addAtom(mol, rH * cosf(5*PI/3), rH * sinf(5*PI/3), 0.0f, ATOM_H); // 13: H on C6
    addAtom(mol, rH * cosf(0), rH * sinf(0) - 0.5f, 0.0f, ATOM_H);    // 14: H on C5 (near chain)

    // Hydrogens on OH groups
    addAtom(mol, r * cosf(PI) - 1.7f, r * sinf(PI) + 0.2f, 0.0f, ATOM_H);       // 15: H on OH
    addAtom(mol, r * cosf(4*PI/3) - 0.2f, r * sinf(4*PI/3) - 1.7f, 0.0f, ATOM_H);// 16: H on OH

    // Hydrogen on side chain CH(OH)
    addAtom(mol, 2.6f, 1.0f, -0.7f, ATOM_H);        // 17: H on CH
    addAtom(mol, 3.7f, 0.5f, 1.6f, ATOM_H);         // 18: H on side chain OH

    // Hydrogens on CH2
    addAtom(mol, 3.6f, -1.5f, -0.5f, ATOM_H);       // 19
    addAtom(mol, 3.2f, -0.5f, -1.8f, ATOM_H);       // 20

    // Hydrogens on NH2
    addAtom(mol, 5.1f, 0.8f, -1.2f, ATOM_H);        // 21
    addAtom(mol, 5.5f, -0.4f, -0.2f, ATOM_H);       // 22

    // Ring bonds (alternating single/double for aromatic)
    for (int i = 0; i < 6; i++) {
        addBond(mol, i, (i + 1) % 6, (i % 2 == 0) ? 2 : 1);
    }

    // OH bonds to ring
    addBond(mol, 3, 6, 1);   // C4-OH
    addBond(mol, 4, 7, 1);   // C3-OH

    // Side chain bonds
    addBond(mol, 0, 8, 1);   // C1-CH(OH)
    addBond(mol, 8, 9, 1);   // CH-OH
    addBond(mol, 8, 10, 1);  // CH-CH2
    addBond(mol, 10, 11, 1); // CH2-NH2

    // H bonds on ring
    addBond(mol, 1, 12, 1);  // C2-H
    addBond(mol, 5, 13, 1);  // C6-H
    addBond(mol, 2, 14, 1);  // C5-H (corrected position)

    // H bonds on OH groups
    addBond(mol, 6, 15, 1);
    addBond(mol, 7, 16, 1);

    // H bonds on side chain
    addBond(mol, 8, 17, 1);  // CH-H
    addBond(mol, 9, 18, 1);  // OH-H
    addBond(mol, 10, 19, 1); // CH2-H
    addBond(mol, 10, 20, 1); // CH2-H
    addBond(mol, 11, 21, 1); // NH2-H
    addBond(mol, 11, 22, 1); // NH2-H

    centerMolecule(mol);
}

// Build ATP (C10H16N5O13P3) - Universal energy currency of cells
void buildATP(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "ATP (C10H16N5O13P3)");

    // Adenine base (purine ring system)
    // 6-membered pyrimidine ring fused with 5-membered imidazole ring
    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_N);         // 0: N1
    addAtom(mol, 1.2f, 0.5f, 0.0f, ATOM_C);         // 1: C2
    addAtom(mol, 2.2f, -0.3f, 0.0f, ATOM_N);        // 2: N3
    addAtom(mol, 2.0f, -1.6f, 0.0f, ATOM_C);        // 3: C4
    addAtom(mol, 0.7f, -2.0f, 0.0f, ATOM_C);        // 4: C5
    addAtom(mol, -0.2f, -1.0f, 0.0f, ATOM_C);       // 5: C6
    addAtom(mol, -1.5f, -1.2f, 0.0f, ATOM_N);       // 6: NH2 (amino)
    addAtom(mol, 3.0f, -2.5f, 0.0f, ATOM_N);        // 7: N7
    addAtom(mol, 2.5f, -3.7f, 0.0f, ATOM_C);        // 8: C8
    addAtom(mol, 1.2f, -3.5f, 0.0f, ATOM_N);        // 9: N9

    // Ribose sugar attached at N9
    addAtom(mol, 0.5f, -4.7f, 0.3f, ATOM_C);        // 10: C1' (anomeric)
    addAtom(mol, 1.0f, -5.8f, -0.5f, ATOM_C);       // 11: C2'
    addAtom(mol, 0.2f, -5.5f, -1.8f, ATOM_C);       // 12: C3'
    addAtom(mol, -1.0f, -4.7f, -1.3f, ATOM_C);      // 13: C4'
    addAtom(mol, -0.9f, -4.5f, 0.2f, ATOM_O);       // 14: O (ring)
    addAtom(mol, 0.8f, -7.1f, 0.0f, ATOM_O);        // 15: 2'-OH
    addAtom(mol, 0.5f, -6.2f, -2.8f, ATOM_O);       // 16: 3'-OH
    addAtom(mol, -2.2f, -5.3f, -1.8f, ATOM_C);      // 17: C5'
    addAtom(mol, -3.3f, -4.5f, -1.3f, ATOM_O);      // 18: 5'-O to phosphate

    // Triphosphate chain
    addAtom(mol, -4.6f, -5.0f, -1.5f, ATOM_P);      // 19: P-alpha
    addAtom(mol, -4.5f, -6.3f, -2.2f, ATOM_O);      // 20: =O
    addAtom(mol, -5.3f, -5.2f, -0.2f, ATOM_O);      // 21: -O-
    addAtom(mol, -5.5f, -4.0f, -2.3f, ATOM_O);      // 22: O to P-beta

    addAtom(mol, -6.8f, -4.3f, -3.0f, ATOM_P);      // 23: P-beta
    addAtom(mol, -6.6f, -5.5f, -3.8f, ATOM_O);      // 24: =O
    addAtom(mol, -7.5f, -4.5f, -1.7f, ATOM_O);      // 25: -O-
    addAtom(mol, -7.6f, -3.2f, -3.8f, ATOM_O);      // 26: O to P-gamma

    addAtom(mol, -8.9f, -3.5f, -4.5f, ATOM_P);      // 27: P-gamma
    addAtom(mol, -8.7f, -4.7f, -5.3f, ATOM_O);      // 28: =O
    addAtom(mol, -9.6f, -3.7f, -3.2f, ATOM_O);      // 29: -O-
    addAtom(mol, -9.6f, -2.4f, -5.3f, ATOM_O);      // 30: -O-

    // Key hydrogens
    addAtom(mol, 1.3f, 1.5f, 0.0f, ATOM_H);         // 31: H on C2
    addAtom(mol, 3.2f, -4.4f, 0.0f, ATOM_H);        // 32: H on C8
    addAtom(mol, -1.9f, -0.4f, 0.3f, ATOM_H);       // 33: H on NH2
    addAtom(mol, -2.0f, -2.0f, -0.3f, ATOM_H);      // 34: H on NH2
    addAtom(mol, 0.8f, -4.8f, 1.3f, ATOM_H);        // 35: H on C1'
    addAtom(mol, 2.0f, -5.7f, -0.8f, ATOM_H);       // 36: H on C2'
    addAtom(mol, 0.0f, -5.0f, -2.5f, ATOM_H);       // 37: H on C3'
    addAtom(mol, -1.0f, -3.7f, -1.7f, ATOM_H);      // 38: H on C4'
    addAtom(mol, -2.1f, -6.3f, -1.5f, ATOM_H);      // 39: H on C5'
    addAtom(mol, -2.4f, -5.3f, -2.8f, ATOM_H);      // 40: H on C5'
    addAtom(mol, 1.5f, -7.5f, -0.3f, ATOM_H);       // 41: H on 2'-OH
    addAtom(mol, 0.0f, -6.9f, -3.0f, ATOM_H);       // 42: H on 3'-OH

    // Adenine ring bonds
    addBond(mol, 0, 1, 2);   // N1=C2
    addBond(mol, 1, 2, 1);   // C2-N3
    addBond(mol, 2, 3, 2);   // N3=C4
    addBond(mol, 3, 4, 1);   // C4-C5
    addBond(mol, 4, 5, 2);   // C5=C6
    addBond(mol, 5, 0, 1);   // C6-N1
    addBond(mol, 5, 6, 1);   // C6-NH2
    addBond(mol, 3, 7, 1);   // C4-N7
    addBond(mol, 7, 8, 2);   // N7=C8
    addBond(mol, 8, 9, 1);   // C8-N9
    addBond(mol, 9, 4, 1);   // N9-C5

    // Ribose bonds
    addBond(mol, 9, 10, 1);  // N9-C1'
    addBond(mol, 10, 11, 1); // C1'-C2'
    addBond(mol, 11, 12, 1); // C2'-C3'
    addBond(mol, 12, 13, 1); // C3'-C4'
    addBond(mol, 13, 14, 1); // C4'-O
    addBond(mol, 14, 10, 1); // O-C1'
    addBond(mol, 11, 15, 1); // C2'-OH
    addBond(mol, 12, 16, 1); // C3'-OH
    addBond(mol, 13, 17, 1); // C4'-C5'
    addBond(mol, 17, 18, 1); // C5'-O

    // Triphosphate bonds
    addBond(mol, 18, 19, 1); // O-P-alpha
    addBond(mol, 19, 20, 2); // P=O
    addBond(mol, 19, 21, 1); // P-O-
    addBond(mol, 19, 22, 1); // P-O-P
    addBond(mol, 22, 23, 1); // O-P-beta
    addBond(mol, 23, 24, 2); // P=O
    addBond(mol, 23, 25, 1); // P-O-
    addBond(mol, 23, 26, 1); // P-O-P
    addBond(mol, 26, 27, 1); // O-P-gamma
    addBond(mol, 27, 28, 2); // P=O
    addBond(mol, 27, 29, 1); // P-O-
    addBond(mol, 27, 30, 1); // P-O-

    // Hydrogen bonds
    addBond(mol, 1, 31, 1);
    addBond(mol, 8, 32, 1);
    addBond(mol, 6, 33, 1);
    addBond(mol, 6, 34, 1);
    addBond(mol, 10, 35, 1);
    addBond(mol, 11, 36, 1);
    addBond(mol, 12, 37, 1);
    addBond(mol, 13, 38, 1);
    addBond(mol, 17, 39, 1);
    addBond(mol, 17, 40, 1);
    addBond(mol, 15, 41, 1);
    addBond(mol, 16, 42, 1);

    centerMolecule(mol);
}

// Build ADP (C10H15N5O10P2) - ATP's partner in energy transfer
void buildADP(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "ADP (C10H15N5O10P2)");

    // Adenine base (same as ATP)
    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_N);         // 0: N1
    addAtom(mol, 1.2f, 0.5f, 0.0f, ATOM_C);         // 1: C2
    addAtom(mol, 2.2f, -0.3f, 0.0f, ATOM_N);        // 2: N3
    addAtom(mol, 2.0f, -1.6f, 0.0f, ATOM_C);        // 3: C4
    addAtom(mol, 0.7f, -2.0f, 0.0f, ATOM_C);        // 4: C5
    addAtom(mol, -0.2f, -1.0f, 0.0f, ATOM_C);       // 5: C6
    addAtom(mol, -1.5f, -1.2f, 0.0f, ATOM_N);       // 6: NH2
    addAtom(mol, 3.0f, -2.5f, 0.0f, ATOM_N);        // 7: N7
    addAtom(mol, 2.5f, -3.7f, 0.0f, ATOM_C);        // 8: C8
    addAtom(mol, 1.2f, -3.5f, 0.0f, ATOM_N);        // 9: N9

    // Ribose sugar
    addAtom(mol, 0.5f, -4.7f, 0.3f, ATOM_C);        // 10: C1'
    addAtom(mol, 1.0f, -5.8f, -0.5f, ATOM_C);       // 11: C2'
    addAtom(mol, 0.2f, -5.5f, -1.8f, ATOM_C);       // 12: C3'
    addAtom(mol, -1.0f, -4.7f, -1.3f, ATOM_C);      // 13: C4'
    addAtom(mol, -0.9f, -4.5f, 0.2f, ATOM_O);       // 14: O ring
    addAtom(mol, 0.8f, -7.1f, 0.0f, ATOM_O);        // 15: 2'-OH
    addAtom(mol, 0.5f, -6.2f, -2.8f, ATOM_O);       // 16: 3'-OH
    addAtom(mol, -2.2f, -5.3f, -1.8f, ATOM_C);      // 17: C5'
    addAtom(mol, -3.3f, -4.5f, -1.3f, ATOM_O);      // 18: 5'-O

    // Diphosphate chain (only 2 phosphates)
    addAtom(mol, -4.6f, -5.0f, -1.5f, ATOM_P);      // 19: P-alpha
    addAtom(mol, -4.5f, -6.3f, -2.2f, ATOM_O);      // 20: =O
    addAtom(mol, -5.3f, -5.2f, -0.2f, ATOM_O);      // 21: -O-
    addAtom(mol, -5.5f, -4.0f, -2.3f, ATOM_O);      // 22: O to P-beta

    addAtom(mol, -6.8f, -4.3f, -3.0f, ATOM_P);      // 23: P-beta
    addAtom(mol, -6.6f, -5.5f, -3.8f, ATOM_O);      // 24: =O
    addAtom(mol, -7.5f, -4.5f, -1.7f, ATOM_O);      // 25: -O-
    addAtom(mol, -7.6f, -3.2f, -3.8f, ATOM_O);      // 26: -O- (terminal)

    // Key hydrogens
    addAtom(mol, 1.3f, 1.5f, 0.0f, ATOM_H);         // 27: H on C2
    addAtom(mol, 3.2f, -4.4f, 0.0f, ATOM_H);        // 28: H on C8
    addAtom(mol, -1.9f, -0.4f, 0.3f, ATOM_H);       // 29: H on NH2
    addAtom(mol, -2.0f, -2.0f, -0.3f, ATOM_H);      // 30: H on NH2
    addAtom(mol, 0.8f, -4.8f, 1.3f, ATOM_H);        // 31: H on C1'
    addAtom(mol, 2.0f, -5.7f, -0.8f, ATOM_H);       // 32: H on C2'
    addAtom(mol, 0.0f, -5.0f, -2.5f, ATOM_H);       // 33: H on C3'
    addAtom(mol, -1.0f, -3.7f, -1.7f, ATOM_H);      // 34: H on C4'
    addAtom(mol, 1.5f, -7.5f, -0.3f, ATOM_H);       // 35: H on 2'-OH
    addAtom(mol, 0.0f, -6.9f, -3.0f, ATOM_H);       // 36: H on 3'-OH

    // Adenine ring bonds
    addBond(mol, 0, 1, 2);
    addBond(mol, 1, 2, 1);
    addBond(mol, 2, 3, 2);
    addBond(mol, 3, 4, 1);
    addBond(mol, 4, 5, 2);
    addBond(mol, 5, 0, 1);
    addBond(mol, 5, 6, 1);
    addBond(mol, 3, 7, 1);
    addBond(mol, 7, 8, 2);
    addBond(mol, 8, 9, 1);
    addBond(mol, 9, 4, 1);

    // Ribose bonds
    addBond(mol, 9, 10, 1);
    addBond(mol, 10, 11, 1);
    addBond(mol, 11, 12, 1);
    addBond(mol, 12, 13, 1);
    addBond(mol, 13, 14, 1);
    addBond(mol, 14, 10, 1);
    addBond(mol, 11, 15, 1);
    addBond(mol, 12, 16, 1);
    addBond(mol, 13, 17, 1);
    addBond(mol, 17, 18, 1);

    // Diphosphate bonds
    addBond(mol, 18, 19, 1);
    addBond(mol, 19, 20, 2);
    addBond(mol, 19, 21, 1);
    addBond(mol, 19, 22, 1);
    addBond(mol, 22, 23, 1);
    addBond(mol, 23, 24, 2);
    addBond(mol, 23, 25, 1);
    addBond(mol, 23, 26, 1);

    // Hydrogen bonds
    addBond(mol, 1, 27, 1);
    addBond(mol, 8, 28, 1);
    addBond(mol, 6, 29, 1);
    addBond(mol, 6, 30, 1);
    addBond(mol, 10, 31, 1);
    addBond(mol, 11, 32, 1);
    addBond(mol, 12, 33, 1);
    addBond(mol, 13, 34, 1);
    addBond(mol, 15, 35, 1);
    addBond(mol, 16, 36, 1);

    centerMolecule(mol);
}

// Build NAD+/NADH (C21H27N7O14P2) - Cellular respiration coenzyme
void buildNADH(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "NADH (C21H29N7O14P2)");

    // Nicotinamide ring (reduced form)
    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_C);         // 0: C2
    addAtom(mol, 1.2f, 0.5f, 0.0f, ATOM_C);         // 1: C3
    addAtom(mol, 2.3f, -0.3f, 0.0f, ATOM_C);        // 2: C4
    addAtom(mol, 2.2f, -1.6f, 0.0f, ATOM_C);        // 3: C5
    addAtom(mol, 1.0f, -2.1f, 0.0f, ATOM_C);        // 4: C6
    addAtom(mol, -0.1f, -1.3f, 0.0f, ATOM_N);       // 5: N1

    // Carboxamide on C3
    addAtom(mol, 1.3f, 1.9f, 0.0f, ATOM_C);         // 6: C=O
    addAtom(mol, 0.3f, 2.6f, 0.3f, ATOM_O);         // 7: =O
    addAtom(mol, 2.5f, 2.5f, -0.3f, ATOM_N);        // 8: NH2

    // Nicotinamide ribose
    addAtom(mol, -1.4f, -1.8f, 0.3f, ATOM_C);       // 9: C1'
    addAtom(mol, -1.8f, -3.0f, -0.5f, ATOM_C);      // 10: C2'
    addAtom(mol, -2.8f, -3.7f, 0.3f, ATOM_C);       // 11: C3'
    addAtom(mol, -3.8f, -2.6f, 0.6f, ATOM_C);       // 12: C4'
    addAtom(mol, -3.0f, -1.5f, 0.0f, ATOM_O);       // 13: O ring
    addAtom(mol, -0.7f, -3.8f, -0.8f, ATOM_O);      // 14: 2'-OH
    addAtom(mol, -3.3f, -4.8f, -0.4f, ATOM_O);      // 15: 3'-OH
    addAtom(mol, -5.0f, -2.8f, -0.2f, ATOM_C);      // 16: C5'
    addAtom(mol, -6.0f, -1.8f, 0.2f, ATOM_O);       // 17: 5'-O

    // Diphosphate bridge
    addAtom(mol, -7.3f, -2.2f, -0.3f, ATOM_P);      // 18: P1
    addAtom(mol, -7.2f, -3.5f, -1.0f, ATOM_O);      // 19: =O
    addAtom(mol, -8.0f, -2.3f, 1.0f, ATOM_O);       // 20: -O-
    addAtom(mol, -8.2f, -1.2f, -1.1f, ATOM_O);      // 21: O-P-O bridge

    addAtom(mol, -9.5f, -1.5f, -1.8f, ATOM_P);      // 22: P2
    addAtom(mol, -9.3f, -2.8f, -2.5f, ATOM_O);      // 23: =O
    addAtom(mol, -10.2f, -1.7f, -0.5f, ATOM_O);     // 24: -O-
    addAtom(mol, -10.3f, -0.5f, -2.6f, ATOM_O);     // 25: 5'-O to adenosine

    // Adenosine ribose
    addAtom(mol, -11.4f, 0.2f, -2.2f, ATOM_C);      // 26: C5'
    addAtom(mol, -12.5f, -0.5f, -1.5f, ATOM_C);     // 27: C4'
    addAtom(mol, -12.3f, -0.3f, 0.0f, ATOM_O);      // 28: O ring
    addAtom(mol, -13.8f, 0.0f, -1.9f, ATOM_C);      // 29: C3'
    addAtom(mol, -13.8f, 0.2f, -3.3f, ATOM_O);      // 30: 3'-OH
    addAtom(mol, -14.5f, 1.1f, -1.2f, ATOM_C);      // 31: C2'
    addAtom(mol, -15.7f, 1.5f, -1.8f, ATOM_O);      // 32: 2'-OH
    addAtom(mol, -13.4f, 1.0f, 0.1f, ATOM_C);       // 33: C1'

    // Adenine base
    addAtom(mol, -13.2f, 1.2f, 1.5f, ATOM_N);       // 34: N9
    addAtom(mol, -14.2f, 1.8f, 2.3f, ATOM_C);       // 35: C4
    addAtom(mol, -13.8f, 2.0f, 3.6f, ATOM_N);       // 36: N3
    addAtom(mol, -14.8f, 2.6f, 4.3f, ATOM_C);       // 37: C2
    addAtom(mol, -16.0f, 2.9f, 3.7f, ATOM_N);       // 38: N1
    addAtom(mol, -16.3f, 2.7f, 2.4f, ATOM_C);       // 39: C6
    addAtom(mol, -17.5f, 3.0f, 1.8f, ATOM_N);       // 40: NH2
    addAtom(mol, -15.4f, 2.1f, 1.7f, ATOM_C);       // 41: C5
    addAtom(mol, -15.4f, 1.8f, 0.3f, ATOM_N);       // 42: N7
    addAtom(mol, -14.2f, 1.3f, 0.1f, ATOM_C);       // 43: C8

    // Key hydrogens
    addAtom(mol, -0.8f, 0.6f, 0.0f, ATOM_H);        // 44
    addAtom(mol, 3.2f, 0.1f, 0.0f, ATOM_H);         // 45
    addAtom(mol, 3.1f, -2.2f, 0.0f, ATOM_H);        // 46
    addAtom(mol, 0.9f, -3.1f, 0.0f, ATOM_H);        // 47
    addAtom(mol, 2.6f, 3.5f, -0.3f, ATOM_H);        // 48
    addAtom(mol, 3.3f, 2.0f, -0.5f, ATOM_H);        // 49
    addAtom(mol, -14.6f, 2.8f, 5.3f, ATOM_H);       // 50
    addAtom(mol, -17.7f, 2.8f, 0.9f, ATOM_H);       // 51
    addAtom(mol, -18.2f, 3.4f, 2.3f, ATOM_H);       // 52

    // Nicotinamide ring bonds
    addBond(mol, 0, 1, 1);
    addBond(mol, 1, 2, 2);
    addBond(mol, 2, 3, 1);
    addBond(mol, 3, 4, 2);
    addBond(mol, 4, 5, 1);
    addBond(mol, 5, 0, 1);
    addBond(mol, 1, 6, 1);
    addBond(mol, 6, 7, 2);
    addBond(mol, 6, 8, 1);

    // Nicotinamide ribose bonds
    addBond(mol, 5, 9, 1);
    addBond(mol, 9, 10, 1);
    addBond(mol, 10, 11, 1);
    addBond(mol, 11, 12, 1);
    addBond(mol, 12, 13, 1);
    addBond(mol, 13, 9, 1);
    addBond(mol, 10, 14, 1);
    addBond(mol, 11, 15, 1);
    addBond(mol, 12, 16, 1);
    addBond(mol, 16, 17, 1);

    // Diphosphate bonds
    addBond(mol, 17, 18, 1);
    addBond(mol, 18, 19, 2);
    addBond(mol, 18, 20, 1);
    addBond(mol, 18, 21, 1);
    addBond(mol, 21, 22, 1);
    addBond(mol, 22, 23, 2);
    addBond(mol, 22, 24, 1);
    addBond(mol, 22, 25, 1);

    // Adenosine ribose bonds
    addBond(mol, 25, 26, 1);
    addBond(mol, 26, 27, 1);
    addBond(mol, 27, 28, 1);
    addBond(mol, 27, 29, 1);
    addBond(mol, 29, 30, 1);
    addBond(mol, 29, 31, 1);
    addBond(mol, 31, 32, 1);
    addBond(mol, 31, 33, 1);
    addBond(mol, 33, 28, 1);

    // Adenine bonds
    addBond(mol, 33, 34, 1);
    addBond(mol, 34, 35, 1);
    addBond(mol, 35, 36, 2);
    addBond(mol, 36, 37, 1);
    addBond(mol, 37, 38, 2);
    addBond(mol, 38, 39, 1);
    addBond(mol, 39, 40, 1);
    addBond(mol, 39, 41, 2);
    addBond(mol, 41, 35, 1);
    addBond(mol, 41, 42, 1);
    addBond(mol, 42, 43, 2);
    addBond(mol, 43, 34, 1);

    // Hydrogen bonds
    addBond(mol, 0, 44, 1);
    addBond(mol, 2, 45, 1);
    addBond(mol, 3, 46, 1);
    addBond(mol, 4, 47, 1);
    addBond(mol, 8, 48, 1);
    addBond(mol, 8, 49, 1);
    addBond(mol, 37, 50, 1);
    addBond(mol, 40, 51, 1);
    addBond(mol, 40, 52, 1);

    centerMolecule(mol);
}

// Build Penicillin G (C16H18N2O4S) - The first antibiotic
void buildPenicillinG(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Penicillin G (C16H18N2O4S)");

    // Beta-lactam ring (4-membered)
    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_C);         // 0: C (carbonyl)
    addAtom(mol, 0.0f, 1.2f, 0.5f, ATOM_N);         // 1: N
    addAtom(mol, 1.2f, 1.5f, 0.0f, ATOM_C);         // 2: C (chiral)
    addAtom(mol, 1.5f, 0.3f, -0.5f, ATOM_C);        // 3: C (chiral)
    addAtom(mol, -1.0f, -0.5f, 0.3f, ATOM_O);       // 4: =O (lactam)

    // Thiazolidine ring (5-membered) fused to beta-lactam
    addAtom(mol, 2.5f, 2.3f, 0.3f, ATOM_C);         // 5: C (gem-dimethyl)
    addAtom(mol, 2.8f, 1.0f, -0.3f, ATOM_S);        // 6: S

    // Carboxylic acid on thiazolidine
    addAtom(mol, 1.8f, -0.8f, -1.0f, ATOM_C);       // 7: COOH
    addAtom(mol, 1.0f, -1.8f, -1.3f, ATOM_O);       // 8: =O
    addAtom(mol, 3.0f, -1.0f, -1.5f, ATOM_O);       // 9: OH

    // Acylamino side chain -NH-C(=O)-CH2-phenyl
    addAtom(mol, -0.5f, 2.2f, 1.2f, ATOM_C);        // 10: C=O (amide)
    addAtom(mol, -1.5f, 2.0f, 2.0f, ATOM_O);        // 11: =O
    addAtom(mol, 0.0f, 3.5f, 1.5f, ATOM_C);         // 12: CH2

    // Benzyl group (phenyl ring)
    addAtom(mol, -0.5f, 4.6f, 0.7f, ATOM_C);        // 13: C1
    addAtom(mol, -1.8f, 4.5f, 0.2f, ATOM_C);        // 14: C2
    addAtom(mol, -2.3f, 5.5f, -0.5f, ATOM_C);       // 15: C3
    addAtom(mol, -1.5f, 6.6f, -0.7f, ATOM_C);       // 16: C4
    addAtom(mol, -0.2f, 6.7f, -0.2f, ATOM_C);       // 17: C5
    addAtom(mol, 0.3f, 5.7f, 0.5f, ATOM_C);         // 18: C6

    // Gem-dimethyl groups on C5
    addAtom(mol, 2.3f, 3.3f, -0.7f, ATOM_C);        // 19: CH3
    addAtom(mol, 3.5f, 2.8f, 1.2f, ATOM_C);         // 20: CH3

    // Key hydrogens
    addAtom(mol, 1.0f, 2.2f, -0.7f, ATOM_H);        // 21: H on C2
    addAtom(mol, 1.5f, 0.0f, 0.5f, ATOM_H);         // 22: H on C3
    addAtom(mol, 3.5f, -1.7f, -1.2f, ATOM_H);       // 23: H on COOH
    addAtom(mol, 0.8f, 3.8f, 2.2f, ATOM_H);         // 24: H on CH2
    addAtom(mol, 0.5f, 3.3f, 2.3f, ATOM_H);         // 25: H on CH2
    addAtom(mol, -2.4f, 3.6f, 0.3f, ATOM_H);        // 26: H on phenyl
    addAtom(mol, -3.3f, 5.4f, -0.9f, ATOM_H);       // 27: H on phenyl
    addAtom(mol, -1.9f, 7.4f, -1.3f, ATOM_H);       // 28: H on phenyl
    addAtom(mol, 0.4f, 7.6f, -0.3f, ATOM_H);        // 29: H on phenyl
    addAtom(mol, 1.3f, 5.8f, 0.9f, ATOM_H);         // 30: H on phenyl

    // Beta-lactam ring bonds
    addBond(mol, 0, 1, 1);   // C-N
    addBond(mol, 1, 2, 1);   // N-C
    addBond(mol, 2, 3, 1);   // C-C
    addBond(mol, 3, 0, 1);   // C-C (closing ring)
    addBond(mol, 0, 4, 2);   // C=O

    // Thiazolidine ring bonds
    addBond(mol, 2, 5, 1);   // C-C
    addBond(mol, 5, 6, 1);   // C-S
    addBond(mol, 6, 3, 1);   // S-C

    // Carboxylic acid
    addBond(mol, 3, 7, 1);   // C-COOH
    addBond(mol, 7, 8, 2);   // C=O
    addBond(mol, 7, 9, 1);   // C-OH

    // Side chain
    addBond(mol, 1, 10, 1);  // N-C=O
    addBond(mol, 10, 11, 2); // C=O
    addBond(mol, 10, 12, 1); // C-CH2
    addBond(mol, 12, 13, 1); // CH2-phenyl

    // Benzyl ring bonds
    addBond(mol, 13, 14, 2);
    addBond(mol, 14, 15, 1);
    addBond(mol, 15, 16, 2);
    addBond(mol, 16, 17, 1);
    addBond(mol, 17, 18, 2);
    addBond(mol, 18, 13, 1);

    // Gem-dimethyl bonds
    addBond(mol, 5, 19, 1);
    addBond(mol, 5, 20, 1);

    // Hydrogen bonds
    addBond(mol, 2, 21, 1);
    addBond(mol, 3, 22, 1);
    addBond(mol, 9, 23, 1);
    addBond(mol, 12, 24, 1);
    addBond(mol, 12, 25, 1);
    addBond(mol, 14, 26, 1);
    addBond(mol, 15, 27, 1);
    addBond(mol, 16, 28, 1);
    addBond(mol, 17, 29, 1);
    addBond(mol, 18, 30, 1);

    centerMolecule(mol);
}

// Build Amoxicillin (C16H19N3O5S) - Most prescribed antibiotic
void buildAmoxicillin(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Amoxicillin (C16H19N3O5S)");

    // Beta-lactam ring (4-membered)
    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_C);         // 0: C (carbonyl)
    addAtom(mol, 0.0f, 1.2f, 0.5f, ATOM_N);         // 1: N
    addAtom(mol, 1.2f, 1.5f, 0.0f, ATOM_C);         // 2: C (chiral)
    addAtom(mol, 1.5f, 0.3f, -0.5f, ATOM_C);        // 3: C (chiral)
    addAtom(mol, -1.0f, -0.5f, 0.3f, ATOM_O);       // 4: =O (lactam)

    // Thiazolidine ring (5-membered) fused to beta-lactam
    addAtom(mol, 2.5f, 2.3f, 0.3f, ATOM_C);         // 5: C (gem-dimethyl)
    addAtom(mol, 2.8f, 1.0f, -0.3f, ATOM_S);        // 6: S

    // Carboxylic acid on thiazolidine
    addAtom(mol, 1.8f, -0.8f, -1.0f, ATOM_C);       // 7: COOH
    addAtom(mol, 1.0f, -1.8f, -1.3f, ATOM_O);       // 8: =O
    addAtom(mol, 3.0f, -1.0f, -1.5f, ATOM_O);       // 9: OH

    // Acylamino side chain with alpha-amino group
    addAtom(mol, -0.5f, 2.2f, 1.2f, ATOM_C);        // 10: C=O (amide)
    addAtom(mol, -1.5f, 2.0f, 2.0f, ATOM_O);        // 11: =O
    addAtom(mol, 0.0f, 3.5f, 1.5f, ATOM_C);         // 12: CH (alpha carbon)
    addAtom(mol, 1.2f, 3.8f, 2.3f, ATOM_N);         // 13: NH2 (alpha amino)

    // 4-Hydroxyphenyl group
    addAtom(mol, -0.5f, 4.6f, 0.7f, ATOM_C);        // 14: C1
    addAtom(mol, -1.8f, 4.5f, 0.2f, ATOM_C);        // 15: C2
    addAtom(mol, -2.3f, 5.5f, -0.5f, ATOM_C);       // 16: C3
    addAtom(mol, -1.5f, 6.6f, -0.7f, ATOM_C);       // 17: C4 (with OH)
    addAtom(mol, -0.2f, 6.7f, -0.2f, ATOM_C);       // 18: C5
    addAtom(mol, 0.3f, 5.7f, 0.5f, ATOM_C);         // 19: C6
    addAtom(mol, -2.0f, 7.7f, -1.4f, ATOM_O);       // 20: OH (para position)

    // Gem-dimethyl groups on C5
    addAtom(mol, 2.3f, 3.3f, -0.7f, ATOM_C);        // 21: CH3
    addAtom(mol, 3.5f, 2.8f, 1.2f, ATOM_C);         // 22: CH3

    // Key hydrogens
    addAtom(mol, 1.0f, 2.2f, -0.7f, ATOM_H);        // 23: H on C2
    addAtom(mol, 1.5f, 0.0f, 0.5f, ATOM_H);         // 24: H on C3
    addAtom(mol, 3.5f, -1.7f, -1.2f, ATOM_H);       // 25: H on COOH
    addAtom(mol, 0.3f, 3.3f, 2.4f, ATOM_H);         // 26: H on CH
    addAtom(mol, 1.0f, 4.5f, 2.9f, ATOM_H);         // 27: H on NH2
    addAtom(mol, 1.8f, 3.2f, 2.6f, ATOM_H);         // 28: H on NH2
    addAtom(mol, -2.4f, 3.6f, 0.3f, ATOM_H);        // 29: H on phenyl
    addAtom(mol, -3.3f, 5.4f, -0.9f, ATOM_H);       // 30: H on phenyl
    addAtom(mol, 0.4f, 7.6f, -0.3f, ATOM_H);        // 31: H on phenyl
    addAtom(mol, 1.3f, 5.8f, 0.9f, ATOM_H);         // 32: H on phenyl
    addAtom(mol, -1.4f, 8.4f, -1.5f, ATOM_H);       // 33: H on phenol OH

    // Beta-lactam ring bonds
    addBond(mol, 0, 1, 1);
    addBond(mol, 1, 2, 1);
    addBond(mol, 2, 3, 1);
    addBond(mol, 3, 0, 1);
    addBond(mol, 0, 4, 2);

    // Thiazolidine ring bonds
    addBond(mol, 2, 5, 1);
    addBond(mol, 5, 6, 1);
    addBond(mol, 6, 3, 1);

    // Carboxylic acid
    addBond(mol, 3, 7, 1);
    addBond(mol, 7, 8, 2);
    addBond(mol, 7, 9, 1);

    // Side chain
    addBond(mol, 1, 10, 1);
    addBond(mol, 10, 11, 2);
    addBond(mol, 10, 12, 1);
    addBond(mol, 12, 13, 1);  // Alpha amino
    addBond(mol, 12, 14, 1);  // To phenyl

    // Hydroxyphenyl ring bonds
    addBond(mol, 14, 15, 2);
    addBond(mol, 15, 16, 1);
    addBond(mol, 16, 17, 2);
    addBond(mol, 17, 18, 1);
    addBond(mol, 18, 19, 2);
    addBond(mol, 19, 14, 1);
    addBond(mol, 17, 20, 1);  // Phenol OH

    // Gem-dimethyl bonds
    addBond(mol, 5, 21, 1);
    addBond(mol, 5, 22, 1);

    // Hydrogen bonds
    addBond(mol, 2, 23, 1);
    addBond(mol, 3, 24, 1);
    addBond(mol, 9, 25, 1);
    addBond(mol, 12, 26, 1);
    addBond(mol, 13, 27, 1);
    addBond(mol, 13, 28, 1);
    addBond(mol, 15, 29, 1);
    addBond(mol, 16, 30, 1);
    addBond(mol, 18, 31, 1);
    addBond(mol, 19, 32, 1);
    addBond(mol, 20, 33, 1);

    centerMolecule(mol);
}

// Build Vancomycin (simplified core) - Last resort antibiotic (glycopeptide)
void buildVancomycin(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Vancomycin Core (C66H75Cl2N9O24)");

    // Vancomycin is very large - this represents the core heptapeptide backbone
    // with key structural features (chlorinated aromatics, sugars simplified)

    // Central aromatic ring 1 (chlorinated phenyl)
    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_C);         // 0
    addAtom(mol, 1.2f, 0.7f, 0.0f, ATOM_C);         // 1
    addAtom(mol, 1.2f, 2.1f, 0.0f, ATOM_C);         // 2
    addAtom(mol, 0.0f, 2.8f, 0.0f, ATOM_C);         // 3
    addAtom(mol, -1.2f, 2.1f, 0.0f, ATOM_C);        // 4
    addAtom(mol, -1.2f, 0.7f, 0.0f, ATOM_C);        // 5
    addAtom(mol, 2.4f, 2.8f, 0.0f, ATOM_CL);        // 6: Cl
    addAtom(mol, 0.0f, 4.2f, 0.0f, ATOM_O);         // 7: O ether link

    // Central aromatic ring 2 (second chlorinated phenyl)
    addAtom(mol, 0.0f, 5.5f, 0.0f, ATOM_C);         // 8
    addAtom(mol, 1.2f, 6.2f, 0.0f, ATOM_C);         // 9
    addAtom(mol, 1.2f, 7.6f, 0.0f, ATOM_C);         // 10
    addAtom(mol, 0.0f, 8.3f, 0.0f, ATOM_C);         // 11
    addAtom(mol, -1.2f, 7.6f, 0.0f, ATOM_C);        // 12
    addAtom(mol, -1.2f, 6.2f, 0.0f, ATOM_C);        // 13
    addAtom(mol, 2.4f, 8.3f, 0.0f, ATOM_CL);        // 14: Cl

    // Peptide backbone segment 1
    addAtom(mol, -2.4f, 0.0f, 0.0f, ATOM_C);        // 15: CH
    addAtom(mol, -2.4f, -1.2f, 0.8f, ATOM_N);       // 16: NH
    addAtom(mol, -3.6f, -1.8f, 0.5f, ATOM_C);       // 17: C=O
    addAtom(mol, -3.8f, -2.8f, 1.2f, ATOM_O);       // 18: =O

    // Peptide backbone segment 2
    addAtom(mol, -4.6f, -1.3f, -0.5f, ATOM_C);      // 19: CH
    addAtom(mol, -5.8f, -0.5f, -0.3f, ATOM_N);      // 20: NH
    addAtom(mol, -6.5f, 0.3f, -1.2f, ATOM_C);       // 21: C=O
    addAtom(mol, -6.3f, 0.3f, -2.4f, ATOM_O);       // 22: =O

    // Hydroxyl groups (glycosylation sites - simplified)
    addAtom(mol, -2.4f, 7.6f, 0.0f, ATOM_O);        // 23: OH
    addAtom(mol, 2.4f, 0.0f, 0.0f, ATOM_O);         // 24: OH

    // Sugar unit 1 (simplified as pyranose ring)
    addAtom(mol, -3.6f, 8.3f, 0.5f, ATOM_C);        // 25
    addAtom(mol, -4.8f, 7.6f, 0.0f, ATOM_C);        // 26
    addAtom(mol, -4.8f, 6.2f, 0.5f, ATOM_C);        // 27
    addAtom(mol, -3.6f, 5.5f, 0.0f, ATOM_C);        // 28
    addAtom(mol, -2.4f, 6.2f, 0.5f, ATOM_O);        // 29: ring O
    addAtom(mol, -3.6f, 4.1f, 0.5f, ATOM_O);        // 30: OH

    // Amino group (characteristic of vancosamine)
    addAtom(mol, -6.0f, 8.3f, 0.5f, ATOM_N);        // 31: NH2 (methylated)
    addAtom(mol, -6.5f, 9.2f, 1.3f, ATOM_C);        // 32: CH3
    addAtom(mol, -6.5f, 7.5f, -0.5f, ATOM_C);       // 33: CH3

    // Key hydrogens
    addAtom(mol, 1.2f, -0.3f, 0.0f, ATOM_H);        // 34
    addAtom(mol, -2.2f, 2.6f, 0.0f, ATOM_H);        // 35
    addAtom(mol, 1.2f, 5.6f, 0.0f, ATOM_H);         // 36
    addAtom(mol, 0.0f, 9.3f, 0.0f, ATOM_H);         // 37
    addAtom(mol, -1.6f, -1.6f, 1.4f, ATOM_H);       // 38: NH
    addAtom(mol, -5.8f, -0.3f, 0.6f, ATOM_H);       // 39: NH

    // Ring 1 bonds
    addBond(mol, 0, 1, 2);
    addBond(mol, 1, 2, 1);
    addBond(mol, 2, 3, 2);
    addBond(mol, 3, 4, 1);
    addBond(mol, 4, 5, 2);
    addBond(mol, 5, 0, 1);
    addBond(mol, 2, 6, 1);   // Cl
    addBond(mol, 3, 7, 1);   // O ether

    // Ring 2 bonds
    addBond(mol, 7, 8, 1);   // ether link
    addBond(mol, 8, 9, 2);
    addBond(mol, 9, 10, 1);
    addBond(mol, 10, 11, 2);
    addBond(mol, 11, 12, 1);
    addBond(mol, 12, 13, 2);
    addBond(mol, 13, 8, 1);
    addBond(mol, 10, 14, 1); // Cl

    // Peptide backbone
    addBond(mol, 5, 15, 1);
    addBond(mol, 15, 16, 1);
    addBond(mol, 16, 17, 1);
    addBond(mol, 17, 18, 2);
    addBond(mol, 17, 19, 1);
    addBond(mol, 19, 20, 1);
    addBond(mol, 20, 21, 1);
    addBond(mol, 21, 22, 2);

    // Hydroxyl groups
    addBond(mol, 12, 23, 1);
    addBond(mol, 1, 24, 1);

    // Sugar unit
    addBond(mol, 23, 25, 1);
    addBond(mol, 25, 26, 1);
    addBond(mol, 26, 27, 1);
    addBond(mol, 27, 28, 1);
    addBond(mol, 28, 29, 1);
    addBond(mol, 29, 25, 1);
    addBond(mol, 28, 30, 1);

    // Amino group
    addBond(mol, 26, 31, 1);
    addBond(mol, 31, 32, 1);
    addBond(mol, 31, 33, 1);

    // Hydrogen bonds
    addBond(mol, 0, 34, 1);
    addBond(mol, 4, 35, 1);
    addBond(mol, 9, 36, 1);
    addBond(mol, 11, 37, 1);
    addBond(mol, 16, 38, 1);
    addBond(mol, 20, 39, 1);

    centerMolecule(mol);
}

// Build TNT (C7H5N3O6) - 2,4,6-Trinitrotoluene - Classic explosive
void buildTNT(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "TNT (C7H5N3O6)");

    // Benzene ring
    float r = 1.4f;
    for (int i = 0; i < 6; i++) {
        float angle = i * PI / 3.0f;
        addAtom(mol, r * cosf(angle), r * sinf(angle), 0.0f, ATOM_C);
    }

    // Methyl group on C1 (position 0)
    addAtom(mol, 2.5f, 0.0f, 0.0f, ATOM_C);         // 6: CH3

    // Nitro group on C2 (position 1) - ortho to methyl
    addAtom(mol, r * cosf(PI/3) + 0.8f, r * sinf(PI/3) + 0.8f, 0.0f, ATOM_N);   // 7: N
    addAtom(mol, r * cosf(PI/3) + 1.5f, r * sinf(PI/3) + 0.3f, 0.3f, ATOM_O);   // 8: O
    addAtom(mol, r * cosf(PI/3) + 0.5f, r * sinf(PI/3) + 1.8f, -0.3f, ATOM_O);  // 9: O

    // Nitro group on C4 (position 3) - para to methyl
    addAtom(mol, r * cosf(PI) - 0.8f, r * sinf(PI), 0.0f, ATOM_N);              // 10: N
    addAtom(mol, r * cosf(PI) - 1.5f, r * sinf(PI) + 0.8f, 0.3f, ATOM_O);       // 11: O
    addAtom(mol, r * cosf(PI) - 1.5f, r * sinf(PI) - 0.8f, -0.3f, ATOM_O);      // 12: O

    // Nitro group on C6 (position 5) - ortho to methyl
    addAtom(mol, r * cosf(5*PI/3) + 0.8f, r * sinf(5*PI/3) - 0.8f, 0.0f, ATOM_N); // 13: N
    addAtom(mol, r * cosf(5*PI/3) + 1.5f, r * sinf(5*PI/3) - 0.3f, 0.3f, ATOM_O); // 14: O
    addAtom(mol, r * cosf(5*PI/3) + 0.5f, r * sinf(5*PI/3) - 1.8f, -0.3f, ATOM_O);// 15: O

    // Hydrogens on methyl group
    addAtom(mol, 3.1f, 0.0f, 0.9f, ATOM_H);         // 16
    addAtom(mol, 3.1f, 0.8f, -0.4f, ATOM_H);        // 17
    addAtom(mol, 3.1f, -0.8f, -0.4f, ATOM_H);       // 18

    // Hydrogens on ring (C3 and C5 - positions 2 and 4)
    float rH = 2.4f;
    addAtom(mol, rH * cosf(2*PI/3), rH * sinf(2*PI/3), 0.0f, ATOM_H);   // 19: H on C3
    addAtom(mol, rH * cosf(4*PI/3), rH * sinf(4*PI/3), 0.0f, ATOM_H);   // 20: H on C5

    // Ring bonds (alternating for aromatic)
    for (int i = 0; i < 6; i++) {
        addBond(mol, i, (i + 1) % 6, (i % 2 == 0) ? 2 : 1);
    }

    // Methyl bond
    addBond(mol, 0, 6, 1);

    // Nitro group bonds
    addBond(mol, 1, 7, 1);   // C-N
    addBond(mol, 7, 8, 2);   // N=O
    addBond(mol, 7, 9, 1);   // N-O

    addBond(mol, 3, 10, 1);  // C-N
    addBond(mol, 10, 11, 2); // N=O
    addBond(mol, 10, 12, 1); // N-O

    addBond(mol, 5, 13, 1);  // C-N
    addBond(mol, 13, 14, 2); // N=O
    addBond(mol, 13, 15, 1); // N-O

    // Hydrogen bonds
    addBond(mol, 6, 16, 1);
    addBond(mol, 6, 17, 1);
    addBond(mol, 6, 18, 1);
    addBond(mol, 2, 19, 1);
    addBond(mol, 4, 20, 1);

    centerMolecule(mol);
}

// Build Nitroglycerin (C3H5N3O9) - Dynamite ingredient, also heart medicine!
void buildNitroglycerin(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Nitroglycerin (C3H5N3O9)");

    // Glycerol backbone: CH2-CH-CH2
    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_C);         // 0: C1
    addAtom(mol, 1.5f, 0.0f, 0.3f, ATOM_C);         // 1: C2 (central)
    addAtom(mol, 3.0f, 0.0f, 0.0f, ATOM_C);         // 2: C3

    // Nitrate ester on C1: -O-NO2
    addAtom(mol, -0.8f, -1.0f, 0.5f, ATOM_O);       // 3: O ester
    addAtom(mol, -1.8f, -0.5f, 1.2f, ATOM_N);       // 4: N
    addAtom(mol, -2.8f, -1.2f, 1.5f, ATOM_O);       // 5: O
    addAtom(mol, -1.8f, 0.7f, 1.5f, ATOM_O);        // 6: O

    // Nitrate ester on C2: -O-NO2
    addAtom(mol, 1.5f, 1.3f, 0.8f, ATOM_O);         // 7: O ester
    addAtom(mol, 1.5f, 2.5f, 0.2f, ATOM_N);         // 8: N
    addAtom(mol, 0.5f, 3.2f, 0.5f, ATOM_O);         // 9: O
    addAtom(mol, 2.5f, 3.0f, -0.4f, ATOM_O);        // 10: O

    // Nitrate ester on C3: -O-NO2
    addAtom(mol, 3.8f, -1.0f, 0.5f, ATOM_O);        // 11: O ester
    addAtom(mol, 4.8f, -0.5f, 1.2f, ATOM_N);        // 12: N
    addAtom(mol, 5.8f, -1.2f, 1.5f, ATOM_O);        // 13: O
    addAtom(mol, 4.8f, 0.7f, 1.5f, ATOM_O);         // 14: O

    // Hydrogens on carbons
    addAtom(mol, 0.0f, 0.0f, -1.0f, ATOM_H);        // 15: H on C1
    addAtom(mol, -0.3f, 0.9f, 0.4f, ATOM_H);        // 16: H on C1
    addAtom(mol, 1.5f, -0.5f, 1.3f, ATOM_H);        // 17: H on C2
    addAtom(mol, 3.0f, 0.0f, -1.0f, ATOM_H);        // 18: H on C3
    addAtom(mol, 3.3f, 0.9f, 0.4f, ATOM_H);         // 19: H on C3

    // Backbone bonds
    addBond(mol, 0, 1, 1);   // C1-C2
    addBond(mol, 1, 2, 1);   // C2-C3

    // Nitrate ester on C1
    addBond(mol, 0, 3, 1);   // C-O
    addBond(mol, 3, 4, 1);   // O-N
    addBond(mol, 4, 5, 2);   // N=O
    addBond(mol, 4, 6, 1);   // N-O

    // Nitrate ester on C2
    addBond(mol, 1, 7, 1);   // C-O
    addBond(mol, 7, 8, 1);   // O-N
    addBond(mol, 8, 9, 2);   // N=O
    addBond(mol, 8, 10, 1);  // N-O

    // Nitrate ester on C3
    addBond(mol, 2, 11, 1);  // C-O
    addBond(mol, 11, 12, 1); // O-N
    addBond(mol, 12, 13, 2); // N=O
    addBond(mol, 12, 14, 1); // N-O

    // Hydrogen bonds
    addBond(mol, 0, 15, 1);
    addBond(mol, 0, 16, 1);
    addBond(mol, 1, 17, 1);
    addBond(mol, 2, 18, 1);
    addBond(mol, 2, 19, 1);

    centerMolecule(mol);
}

// Build RDX (C3H6N6O6) - Cyclotrimethylenetrinitramine - Military explosive
void buildRDX(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "RDX (C3H6N6O6)");

    // 1,3,5-triazine ring with N-NO2 groups
    // Six-membered ring alternating C and N
    float r = 1.3f;

    // Ring atoms: C-N-C-N-C-N (positions 0,1,2,3,4,5)
    addAtom(mol, r * cosf(0), r * sinf(0), 0.0f, ATOM_C);           // 0: C
    addAtom(mol, r * cosf(PI/3), r * sinf(PI/3), 0.0f, ATOM_N);     // 1: N (with NO2)
    addAtom(mol, r * cosf(2*PI/3), r * sinf(2*PI/3), 0.0f, ATOM_C); // 2: C
    addAtom(mol, r * cosf(PI), r * sinf(PI), 0.0f, ATOM_N);         // 3: N (with NO2)
    addAtom(mol, r * cosf(4*PI/3), r * sinf(4*PI/3), 0.0f, ATOM_C); // 4: C
    addAtom(mol, r * cosf(5*PI/3), r * sinf(5*PI/3), 0.0f, ATOM_N); // 5: N (with NO2)

    // Nitro group on N1 (position 1)
    addAtom(mol, r * cosf(PI/3) + 0.5f, r * sinf(PI/3) + 1.2f, 0.3f, ATOM_N);   // 6: N
    addAtom(mol, r * cosf(PI/3) + 1.3f, r * sinf(PI/3) + 1.8f, 0.0f, ATOM_O);   // 7: O
    addAtom(mol, r * cosf(PI/3) - 0.3f, r * sinf(PI/3) + 2.0f, 0.8f, ATOM_O);   // 8: O

    // Nitro group on N3 (position 3)
    addAtom(mol, r * cosf(PI) - 1.2f, r * sinf(PI), 0.5f, ATOM_N);              // 9: N
    addAtom(mol, r * cosf(PI) - 1.8f, r * sinf(PI) + 0.9f, 0.3f, ATOM_O);       // 10: O
    addAtom(mol, r * cosf(PI) - 1.8f, r * sinf(PI) - 0.9f, 0.8f, ATOM_O);       // 11: O

    // Nitro group on N5 (position 5)
    addAtom(mol, r * cosf(5*PI/3) + 0.5f, r * sinf(5*PI/3) - 1.2f, 0.3f, ATOM_N); // 12: N
    addAtom(mol, r * cosf(5*PI/3) + 1.3f, r * sinf(5*PI/3) - 1.8f, 0.0f, ATOM_O); // 13: O
    addAtom(mol, r * cosf(5*PI/3) - 0.3f, r * sinf(5*PI/3) - 2.0f, 0.8f, ATOM_O); // 14: O

    // Hydrogens on carbons (2 each)
    addAtom(mol, r * cosf(0) + 0.8f, r * sinf(0) + 0.5f, 0.8f, ATOM_H);         // 15
    addAtom(mol, r * cosf(0) + 0.8f, r * sinf(0) - 0.5f, -0.8f, ATOM_H);        // 16
    addAtom(mol, r * cosf(2*PI/3) - 0.3f, r * sinf(2*PI/3) + 0.9f, 0.8f, ATOM_H);   // 17
    addAtom(mol, r * cosf(2*PI/3) - 0.3f, r * sinf(2*PI/3) + 0.3f, -0.9f, ATOM_H);  // 18
    addAtom(mol, r * cosf(4*PI/3) - 0.3f, r * sinf(4*PI/3) - 0.9f, 0.8f, ATOM_H);   // 19
    addAtom(mol, r * cosf(4*PI/3) - 0.3f, r * sinf(4*PI/3) - 0.3f, -0.9f, ATOM_H);  // 20

    // Ring bonds
    addBond(mol, 0, 1, 1);
    addBond(mol, 1, 2, 1);
    addBond(mol, 2, 3, 1);
    addBond(mol, 3, 4, 1);
    addBond(mol, 4, 5, 1);
    addBond(mol, 5, 0, 1);

    // Nitro group bonds
    addBond(mol, 1, 6, 1);   // N-N
    addBond(mol, 6, 7, 2);   // N=O
    addBond(mol, 6, 8, 1);   // N-O

    addBond(mol, 3, 9, 1);   // N-N
    addBond(mol, 9, 10, 2);  // N=O
    addBond(mol, 9, 11, 1);  // N-O

    addBond(mol, 5, 12, 1);  // N-N
    addBond(mol, 12, 13, 2); // N=O
    addBond(mol, 12, 14, 1); // N-O

    // Hydrogen bonds
    addBond(mol, 0, 15, 1);
    addBond(mol, 0, 16, 1);
    addBond(mol, 2, 17, 1);
    addBond(mol, 2, 18, 1);
    addBond(mol, 4, 19, 1);
    addBond(mol, 4, 20, 1);

    centerMolecule(mol);
}

// Build Sucrose (C12H22O11) - Table sugar (glucose + fructose disaccharide)
void buildSucrose(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Sucrose (C12H22O11)");

    // Glucose ring (6-membered pyranose)
    // C1-C5 + ring oxygen
    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_C);         // 0: C1 (anomeric)
    addAtom(mol, 1.2f, 0.7f, 0.3f, ATOM_C);         // 1: C2
    addAtom(mol, 2.4f, 0.0f, 0.0f, ATOM_C);         // 2: C3
    addAtom(mol, 2.4f, -1.4f, 0.5f, ATOM_C);        // 3: C4
    addAtom(mol, 1.2f, -2.1f, 0.2f, ATOM_C);        // 4: C5
    addAtom(mol, 0.0f, -1.4f, 0.5f, ATOM_O);        // 5: ring O

    // Glucose hydroxyl groups and CH2OH
    addAtom(mol, -1.0f, 0.5f, 0.5f, ATOM_O);        // 6: OH on C1
    addAtom(mol, 1.2f, 2.1f, 0.0f, ATOM_O);         // 7: OH on C2
    addAtom(mol, 3.5f, 0.7f, 0.3f, ATOM_O);         // 8: OH on C3
    addAtom(mol, 3.5f, -2.1f, 0.2f, ATOM_O);        // 9: OH on C4
    addAtom(mol, 1.2f, -3.5f, 0.5f, ATOM_C);        // 10: C6 (CH2OH)
    addAtom(mol, 1.2f, -4.2f, 1.5f, ATOM_O);        // 11: OH on C6

    // Fructose ring (5-membered furanose)
    // Connected via glycosidic bond from glucose C1 to fructose C2
    addAtom(mol, -2.2f, 0.0f, -0.3f, ATOM_C);       // 12: C1 (CH2OH)
    addAtom(mol, -3.0f, 1.2f, 0.0f, ATOM_C);        // 13: C2 (anomeric, ketose)
    addAtom(mol, -4.4f, 0.8f, -0.3f, ATOM_C);       // 14: C3
    addAtom(mol, -4.7f, -0.5f, 0.4f, ATOM_C);       // 15: C4
    addAtom(mol, -3.4f, -1.2f, 0.2f, ATOM_C);       // 16: C5
    addAtom(mol, -2.8f, -0.5f, -0.9f, ATOM_O);      // 17: ring O

    // Fructose hydroxyl groups and CH2OH
    addAtom(mol, -1.8f, 0.5f, -1.3f, ATOM_O);       // 18: OH on C1
    addAtom(mol, -5.3f, 1.7f, 0.0f, ATOM_O);        // 19: OH on C3
    addAtom(mol, -5.8f, -1.2f, 0.2f, ATOM_O);       // 20: OH on C4
    addAtom(mol, -3.4f, -2.6f, 0.5f, ATOM_C);       // 21: C6 (CH2OH)
    addAtom(mol, -3.4f, -3.3f, -0.5f, ATOM_O);      // 22: OH on C6

    // Key hydrogens (selected for clarity)
    addAtom(mol, 0.0f, 0.0f, 1.0f, ATOM_H);         // 23: H on C1
    addAtom(mol, 1.2f, 0.7f, 1.3f, ATOM_H);         // 24: H on C2
    addAtom(mol, 2.4f, 0.0f, -1.0f, ATOM_H);        // 25: H on C3
    addAtom(mol, 2.4f, -1.4f, 1.5f, ATOM_H);        // 26: H on C4
    addAtom(mol, 1.2f, -2.1f, -0.8f, ATOM_H);       // 27: H on C5

    // Glucose ring bonds
    addBond(mol, 0, 1, 1);
    addBond(mol, 1, 2, 1);
    addBond(mol, 2, 3, 1);
    addBond(mol, 3, 4, 1);
    addBond(mol, 4, 5, 1);
    addBond(mol, 5, 0, 1);  // ring closure

    // Glucose substituents
    addBond(mol, 0, 6, 1);  // C1-OH (anomeric)
    addBond(mol, 1, 7, 1);  // C2-OH
    addBond(mol, 2, 8, 1);  // C3-OH
    addBond(mol, 3, 9, 1);  // C4-OH
    addBond(mol, 4, 10, 1); // C5-C6
    addBond(mol, 10, 11, 1);// C6-OH

    // Glycosidic bond (glucose C1 O to fructose C2)
    addBond(mol, 6, 13, 1);

    // Fructose ring bonds
    addBond(mol, 12, 13, 1);
    addBond(mol, 13, 14, 1);
    addBond(mol, 14, 15, 1);
    addBond(mol, 15, 16, 1);
    addBond(mol, 16, 17, 1);
    addBond(mol, 17, 12, 1); // ring closure

    // Fructose substituents
    addBond(mol, 12, 18, 1); // C1-OH
    addBond(mol, 14, 19, 1); // C3-OH
    addBond(mol, 15, 20, 1); // C4-OH
    addBond(mol, 16, 21, 1); // C5-C6
    addBond(mol, 21, 22, 1); // C6-OH

    // Selected hydrogen bonds
    addBond(mol, 0, 23, 1);
    addBond(mol, 1, 24, 1);
    addBond(mol, 2, 25, 1);
    addBond(mol, 3, 26, 1);
    addBond(mol, 4, 27, 1);

    centerMolecule(mol);
}

// Build Aspartame (C14H18N2O5) - Diet soda sweetener (aspartyl-phenylalanine methyl ester)
void buildAspartame(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Aspartame (C14H18N2O5)");

    // Aspartic acid part (left side)
    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_C);         // 0: C alpha
    addAtom(mol, -0.5f, 1.3f, 0.3f, ATOM_N);        // 1: N amino
    addAtom(mol, 0.8f, -0.7f, 1.2f, ATOM_C);        // 2: C beta
    addAtom(mol, 0.3f, -0.5f, 2.5f, ATOM_C);        // 3: C gamma (carboxyl)
    addAtom(mol, -0.7f, 0.2f, 2.8f, ATOM_O);        // 4: O (C=O)
    addAtom(mol, 0.9f, -1.2f, 3.5f, ATOM_O);        // 5: OH

    // Peptide bond to phenylalanine
    addAtom(mol, 1.2f, -0.3f, -0.8f, ATOM_C);       // 6: C carbonyl
    addAtom(mol, 1.0f, -0.5f, -2.0f, ATOM_O);       // 7: O (C=O)
    addAtom(mol, 2.3f, 0.4f, -0.4f, ATOM_N);        // 8: N amide

    // Phenylalanine part
    addAtom(mol, 3.5f, 0.5f, -1.2f, ATOM_C);        // 9: C alpha
    addAtom(mol, 4.3f, -0.8f, -1.5f, ATOM_C);       // 10: C beta

    // Phenyl ring
    float r = 1.4f;
    float baseX = 5.5f;
    float baseY = -0.8f;
    for (int i = 0; i < 6; i++) {
        float angle = i * PI / 3.0f;
        addAtom(mol, baseX + r * cosf(angle), baseY + r * sinf(angle), 0.0f, ATOM_C);  // 11-16
    }

    // Methyl ester on C-terminus
    addAtom(mol, 4.3f, 1.5f, -0.5f, ATOM_C);        // 17: C carbonyl
    addAtom(mol, 4.2f, 2.7f, -0.8f, ATOM_O);        // 18: O (C=O)
    addAtom(mol, 5.3f, 1.0f, 0.3f, ATOM_O);         // 19: O ester
    addAtom(mol, 6.3f, 1.8f, 0.6f, ATOM_C);         // 20: CH3

    // Key hydrogens
    addAtom(mol, -1.4f, 1.2f, 0.7f, ATOM_H);        // 21: NH2
    addAtom(mol, -0.3f, 2.1f, -0.3f, ATOM_H);       // 22: NH2
    addAtom(mol, 2.2f, 1.3f, 0.0f, ATOM_H);         // 23: NH
    addAtom(mol, 0.0f, -0.3f, -0.9f, ATOM_H);       // 24: H on C alpha
    addAtom(mol, 3.3f, 1.0f, -2.1f, ATOM_H);        // 25: H on C alpha

    // Aspartic acid bonds
    addBond(mol, 0, 1, 1);   // C-N
    addBond(mol, 0, 2, 1);   // C-C beta
    addBond(mol, 2, 3, 1);   // C beta - C gamma
    addBond(mol, 3, 4, 2);   // C=O
    addBond(mol, 3, 5, 1);   // C-OH

    // Peptide bond
    addBond(mol, 0, 6, 1);   // C alpha - C carbonyl
    addBond(mol, 6, 7, 2);   // C=O
    addBond(mol, 6, 8, 1);   // C-N amide

    // Phenylalanine bonds
    addBond(mol, 8, 9, 1);   // N-C alpha
    addBond(mol, 9, 10, 1);  // C alpha - C beta
    addBond(mol, 10, 11, 1); // C beta - phenyl ring

    // Phenyl ring
    for (int i = 0; i < 6; i++) {
        addBond(mol, 11 + i, 11 + ((i + 1) % 6), (i % 2 == 0) ? 2 : 1);
    }

    // Ester bonds
    addBond(mol, 9, 17, 1);  // C alpha - C carbonyl
    addBond(mol, 17, 18, 2); // C=O
    addBond(mol, 17, 19, 1); // C-O ester
    addBond(mol, 19, 20, 1); // O-CH3

    // Hydrogen bonds
    addBond(mol, 1, 21, 1);
    addBond(mol, 1, 22, 1);
    addBond(mol, 8, 23, 1);
    addBond(mol, 0, 24, 1);
    addBond(mol, 9, 25, 1);
}

// Build DDT (C14H9Cl5) - Dichlorodiphenyltrichloroethane - Historic insecticide
void buildDDT(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "DDT (C14H9Cl5)");

    // Central carbon with CCl3 group
    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_C);         // 0: Central CH
    addAtom(mol, 0.0f, 1.5f, 0.0f, ATOM_C);         // 1: CCl3 carbon
    addAtom(mol, -1.0f, 2.2f, 0.8f, ATOM_CL);       // 2: Cl
    addAtom(mol, 1.0f, 2.2f, 0.8f, ATOM_CL);        // 3: Cl
    addAtom(mol, 0.0f, 2.3f, -1.0f, ATOM_CL);       // 4: Cl

    // First para-chlorophenyl ring
    float r = 1.4f;
    float offsetX = -2.0f;
    for (int i = 0; i < 6; i++) {
        float angle = i * PI / 3.0f;
        addAtom(mol, offsetX + r * cosf(angle), r * sinf(angle) - 1.0f, 0.0f, ATOM_C);
    }
    // Para-chlorine on ring 1
    addAtom(mol, offsetX - 2.4f, -1.0f, 0.0f, ATOM_CL);  // 11: Cl para

    // Second para-chlorophenyl ring
    offsetX = 2.0f;
    for (int i = 0; i < 6; i++) {
        float angle = i * PI / 3.0f;
        addAtom(mol, offsetX + r * cosf(angle), r * sinf(angle) - 1.0f, 0.0f, ATOM_C);
    }
    // Para-chlorine on ring 2
    addAtom(mol, offsetX + 2.4f, -1.0f, 0.0f, ATOM_CL);  // 18: Cl para

    // Hydrogen on central carbon
    addAtom(mol, 0.0f, 0.0f, -1.0f, ATOM_H);        // 19: H on central C

    // Hydrogens on ring 1 (4 H's - ortho and meta positions)
    float rH = 2.4f;
    addAtom(mol, -2.0f + rH * cosf(PI/3), rH * sinf(PI/3) - 1.0f, 0.0f, ATOM_H);   // 20
    addAtom(mol, -2.0f + rH * cosf(2*PI/3), rH * sinf(2*PI/3) - 1.0f, 0.0f, ATOM_H); // 21
    addAtom(mol, -2.0f + rH * cosf(4*PI/3), rH * sinf(4*PI/3) - 1.0f, 0.0f, ATOM_H); // 22
    addAtom(mol, -2.0f + rH * cosf(5*PI/3), rH * sinf(5*PI/3) - 1.0f, 0.0f, ATOM_H); // 23

    // Hydrogens on ring 2 (4 H's)
    addAtom(mol, 2.0f + rH * cosf(PI/3), rH * sinf(PI/3) - 1.0f, 0.0f, ATOM_H);    // 24
    addAtom(mol, 2.0f + rH * cosf(2*PI/3), rH * sinf(2*PI/3) - 1.0f, 0.0f, ATOM_H);  // 25
    addAtom(mol, 2.0f + rH * cosf(4*PI/3), rH * sinf(4*PI/3) - 1.0f, 0.0f, ATOM_H);  // 26
    addAtom(mol, 2.0f + rH * cosf(5*PI/3), rH * sinf(5*PI/3) - 1.0f, 0.0f, ATOM_H);  // 27

    // Central bonds
    addBond(mol, 0, 1, 1);   // CH-CCl3
    addBond(mol, 1, 2, 1);   // C-Cl
    addBond(mol, 1, 3, 1);   // C-Cl
    addBond(mol, 1, 4, 1);   // C-Cl

    // Connect central C to rings
    addBond(mol, 0, 5, 1);   // To ring 1 C1
    addBond(mol, 0, 12, 1);  // To ring 2 C1

    // Ring 1 bonds
    for (int i = 0; i < 6; i++) {
        addBond(mol, 5 + i, 5 + (i + 1) % 6, (i % 2 == 0) ? 2 : 1);
    }
    addBond(mol, 8, 11, 1);  // Para-Cl

    // Ring 2 bonds
    for (int i = 0; i < 6; i++) {
        addBond(mol, 12 + i, 12 + (i + 1) % 6, (i % 2 == 0) ? 2 : 1);
    }
    addBond(mol, 15, 18, 1); // Para-Cl

    // Hydrogen bonds
    addBond(mol, 0, 19, 1);
    addBond(mol, 6, 20, 1);
    addBond(mol, 7, 21, 1);
    addBond(mol, 9, 22, 1);
    addBond(mol, 10, 23, 1);
    addBond(mol, 13, 24, 1);
    addBond(mol, 14, 25, 1);
    addBond(mol, 16, 26, 1);
    addBond(mol, 17, 27, 1);

    centerMolecule(mol);
}

// Build Saccharin (C7H5NO3S) - Original artificial sweetener (benzoic sulfimide)
void buildSaccharin(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Saccharin (C7H5NO3S)");

    // Benzene ring
    float r = 1.4f;
    for (int i = 0; i < 6; i++) {
        float angle = i * PI / 3.0f;
        addAtom(mol, r * cosf(angle), r * sinf(angle), 0.0f, ATOM_C);  // 0-5
    }

    // Fused heterocycle (isothiazole S-dioxide ring)
    // Connected to C1-C2 of benzene
    addAtom(mol, r * cosf(PI/3) + 1.0f, r * sinf(PI/3) + 1.0f, 0.5f, ATOM_S);   // 6: S
    addAtom(mol, r * cosf(PI/3) + 0.3f, r * sinf(PI/3) + 2.2f, 0.3f, ATOM_O);   // 7: O on S (sulfonyl)
    addAtom(mol, r * cosf(PI/3) + 2.2f, r * sinf(PI/3) + 1.0f, 0.3f, ATOM_O);   // 8: O on S (sulfonyl)
    addAtom(mol, r * cosf(0) + 0.8f, r * sinf(0) + 1.3f, 0.5f, ATOM_N);         // 9: N
    addAtom(mol, r * cosf(0) + 1.5f, r * sinf(0) + 0.3f, 0.3f, ATOM_C);         // 10: C=O
    addAtom(mol, r * cosf(0) + 2.7f, r * sinf(0) + 0.5f, 0.5f, ATOM_O);         // 11: O (C=O)

    // Hydrogens on benzene
    float rH = 2.4f;
    addAtom(mol, rH * cosf(2*PI/3), rH * sinf(2*PI/3), 0.0f, ATOM_H);   // 12: H
    addAtom(mol, rH * cosf(PI), rH * sinf(PI), 0.0f, ATOM_H);           // 13: H
    addAtom(mol, rH * cosf(4*PI/3), rH * sinf(4*PI/3), 0.0f, ATOM_H);   // 14: H
    addAtom(mol, rH * cosf(5*PI/3), rH * sinf(5*PI/3), 0.0f, ATOM_H);   // 15: H

    // Hydrogen on nitrogen (acidic)
    addAtom(mol, r * cosf(0) + 0.5f, r * sinf(0) + 2.2f, 0.8f, ATOM_H); // 16: NH

    // Benzene ring bonds
    for (int i = 0; i < 6; i++) {
        addBond(mol, i, (i + 1) % 6, (i % 2 == 0) ? 2 : 1);
    }

    // Fused ring bonds
    addBond(mol, 1, 6, 1);   // C-S
    addBond(mol, 6, 7, 2);   // S=O
    addBond(mol, 6, 8, 2);   // S=O
    addBond(mol, 6, 9, 1);   // S-N
    addBond(mol, 9, 10, 1);  // N-C
    addBond(mol, 10, 11, 2); // C=O
    addBond(mol, 10, 0, 1);  // C-C (ring fusion)

    // Hydrogen bonds
    addBond(mol, 2, 12, 1);
    addBond(mol, 3, 13, 1);
    addBond(mol, 4, 14, 1);
    addBond(mol, 5, 15, 1);
    addBond(mol, 9, 16, 1);
}

// Build Glyphosate (C3H8NO5P) - World's most used herbicide (Roundup)
void buildGlyphosate(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Glyphosate (C3H8NO5P)");

    // Structure: HOOC-CH2-NH-CH2-PO(OH)2
    // Carboxylic acid group
    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_C);         // 0: COOH carbon
    addAtom(mol, -1.0f, 0.5f, 0.5f, ATOM_O);        // 1: =O
    addAtom(mol, 0.5f, -1.0f, 0.5f, ATOM_O);        // 2: OH

    // CH2 (alpha)
    addAtom(mol, 1.3f, 0.5f, -0.3f, ATOM_C);        // 3: CH2

    // NH (secondary amine)
    addAtom(mol, 2.5f, -0.2f, 0.2f, ATOM_N);        // 4: NH

    // CH2 (connected to phosphonate)
    addAtom(mol, 3.7f, 0.5f, -0.3f, ATOM_C);        // 5: CH2

    // Phosphonate group PO(OH)2
    addAtom(mol, 5.0f, -0.2f, 0.2f, ATOM_P);        // 6: P
    addAtom(mol, 5.5f, -1.3f, -0.5f, ATOM_O);       // 7: =O
    addAtom(mol, 6.0f, 0.8f, -0.5f, ATOM_O);        // 8: OH
    addAtom(mol, 5.0f, 0.3f, 1.6f, ATOM_O);         // 9: OH

    // Hydrogens
    addAtom(mol, 0.0f, -1.5f, 1.0f, ATOM_H);        // 10: H on COOH
    addAtom(mol, 1.2f, 1.5f, 0.0f, ATOM_H);         // 11: H on CH2
    addAtom(mol, 1.4f, 0.5f, -1.3f, ATOM_H);        // 12: H on CH2
    addAtom(mol, 2.5f, -1.1f, 0.0f, ATOM_H);        // 13: H on NH
    addAtom(mol, 3.6f, 1.5f, 0.0f, ATOM_H);         // 14: H on CH2
    addAtom(mol, 3.8f, 0.5f, -1.3f, ATOM_H);        // 15: H on CH2
    addAtom(mol, 6.8f, 0.5f, -0.8f, ATOM_H);        // 16: H on P-OH
    addAtom(mol, 5.5f, 0.0f, 2.2f, ATOM_H);         // 17: H on P-OH

    // Bonds
    addBond(mol, 0, 1, 2);   // C=O
    addBond(mol, 0, 2, 1);   // C-OH
    addBond(mol, 0, 3, 1);   // C-CH2
    addBond(mol, 3, 4, 1);   // CH2-NH
    addBond(mol, 4, 5, 1);   // NH-CH2
    addBond(mol, 5, 6, 1);   // CH2-P
    addBond(mol, 6, 7, 2);   // P=O
    addBond(mol, 6, 8, 1);   // P-OH
    addBond(mol, 6, 9, 1);   // P-OH

    // Hydrogen bonds
    addBond(mol, 2, 10, 1);
    addBond(mol, 3, 11, 1);
    addBond(mol, 3, 12, 1);
    addBond(mol, 4, 13, 1);
    addBond(mol, 5, 14, 1);
    addBond(mol, 5, 15, 1);
    addBond(mol, 8, 16, 1);
    addBond(mol, 9, 17, 1);

    centerMolecule(mol);
}

// Build Sucralose (C12H19Cl3O8) - Splenda (chlorinated sucrose)
void buildSucralose(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Sucralose (C12H19Cl3O8)");

    // Similar to sucrose but with 3 Cl substitutions (positions 4-OH, 1'-OH, 6'-OH)
    // Glucose ring (6-membered pyranose)
    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_C);         // 0: C1 (anomeric)
    addAtom(mol, 1.2f, 0.7f, 0.3f, ATOM_C);         // 1: C2
    addAtom(mol, 2.4f, 0.0f, 0.0f, ATOM_C);         // 2: C3
    addAtom(mol, 2.4f, -1.4f, 0.5f, ATOM_C);        // 3: C4
    addAtom(mol, 1.2f, -2.1f, 0.2f, ATOM_C);        // 4: C5
    addAtom(mol, 0.0f, -1.4f, 0.5f, ATOM_O);        // 5: ring O

    // Glucose substituents (one Cl on C4)
    addAtom(mol, -1.0f, 0.5f, 0.5f, ATOM_O);        // 6: OH on C1
    addAtom(mol, 1.2f, 2.1f, 0.0f, ATOM_O);         // 7: OH on C2
    addAtom(mol, 3.5f, 0.7f, 0.3f, ATOM_O);         // 8: OH on C3
    addAtom(mol, 3.5f, -2.1f, 0.2f, ATOM_CL);       // 9: Cl on C4 (substitution!)
    addAtom(mol, 1.2f, -3.5f, 0.5f, ATOM_C);        // 10: C6 (CH2OH)
    addAtom(mol, 1.2f, -4.2f, 1.5f, ATOM_O);        // 11: OH on C6

    // Fructose ring (5-membered furanose)
    addAtom(mol, -2.2f, 0.0f, -0.3f, ATOM_C);       // 12: C1 (CH2Cl - substitution!)
    addAtom(mol, -3.0f, 1.2f, 0.0f, ATOM_C);        // 13: C2 (anomeric)
    addAtom(mol, -4.4f, 0.8f, -0.3f, ATOM_C);       // 14: C3
    addAtom(mol, -4.7f, -0.5f, 0.4f, ATOM_C);       // 15: C4
    addAtom(mol, -3.4f, -1.2f, 0.2f, ATOM_C);       // 16: C5
    addAtom(mol, -2.8f, -0.5f, -0.9f, ATOM_O);      // 17: ring O

    // Fructose substituents (Cl on C1' and C6')
    addAtom(mol, -1.8f, 0.5f, -1.3f, ATOM_CL);      // 18: Cl on C1' (substitution!)
    addAtom(mol, -5.3f, 1.7f, 0.0f, ATOM_O);        // 19: OH on C3'
    addAtom(mol, -5.8f, -1.2f, 0.2f, ATOM_O);       // 20: OH on C4'
    addAtom(mol, -3.4f, -2.6f, 0.5f, ATOM_C);       // 21: C6' (CH2Cl)
    addAtom(mol, -3.4f, -3.3f, -0.5f, ATOM_CL);     // 22: Cl on C6' (substitution!)

    // Selected hydrogens
    addAtom(mol, 0.0f, 0.0f, 1.0f, ATOM_H);         // 23: H on C1
    addAtom(mol, 1.2f, 0.7f, 1.3f, ATOM_H);         // 24: H on C2
    addAtom(mol, 2.4f, 0.0f, -1.0f, ATOM_H);        // 25: H on C3
    addAtom(mol, 2.4f, -1.4f, 1.5f, ATOM_H);        // 26: H on C4
    addAtom(mol, 1.2f, -2.1f, -0.8f, ATOM_H);       // 27: H on C5

    // Glucose ring bonds
    addBond(mol, 0, 1, 1);
    addBond(mol, 1, 2, 1);
    addBond(mol, 2, 3, 1);
    addBond(mol, 3, 4, 1);
    addBond(mol, 4, 5, 1);
    addBond(mol, 5, 0, 1);

    // Glucose substituents
    addBond(mol, 0, 6, 1);
    addBond(mol, 1, 7, 1);
    addBond(mol, 2, 8, 1);
    addBond(mol, 3, 9, 1);  // C4-Cl
    addBond(mol, 4, 10, 1);
    addBond(mol, 10, 11, 1);

    // Glycosidic bond
    addBond(mol, 6, 13, 1);

    // Fructose ring bonds
    addBond(mol, 12, 13, 1);
    addBond(mol, 13, 14, 1);
    addBond(mol, 14, 15, 1);
    addBond(mol, 15, 16, 1);
    addBond(mol, 16, 17, 1);
    addBond(mol, 17, 12, 1);

    // Fructose substituents
    addBond(mol, 12, 18, 1); // C1'-Cl
    addBond(mol, 14, 19, 1);
    addBond(mol, 15, 20, 1);
    addBond(mol, 16, 21, 1);
    addBond(mol, 21, 22, 1); // C6'-Cl

    // Hydrogen bonds
    addBond(mol, 0, 23, 1);
    addBond(mol, 1, 24, 1);
    addBond(mol, 2, 25, 1);
    addBond(mol, 3, 26, 1);
    addBond(mol, 4, 27, 1);
}

// Build Malathion (C10H19O6PS2) - Common organophosphate insecticide
void buildMalathion(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Malathion (C10H19O6PS2)");

    // Thiophosphate core
    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_P);         // 0: P
    addAtom(mol, 0.0f, 1.5f, 0.0f, ATOM_S);         // 1: =S (thiophosphoryl)
    addAtom(mol, -1.3f, -0.5f, 0.8f, ATOM_O);       // 2: O-CH3
    addAtom(mol, 1.3f, -0.5f, 0.8f, ATOM_O);        // 3: O-CH3
    addAtom(mol, 0.0f, -0.8f, -1.3f, ATOM_S);       // 4: S-CH (thioether)

    // Two methoxy groups
    addAtom(mol, -2.5f, 0.0f, 0.5f, ATOM_C);        // 5: CH3
    addAtom(mol, 2.5f, 0.0f, 0.5f, ATOM_C);         // 6: CH3

    // Succinate backbone: -S-CH(COOCH3)-CH2-COOCH3
    addAtom(mol, 0.0f, -2.2f, -1.8f, ATOM_C);       // 7: CH (chiral center)
    addAtom(mol, -1.3f, -2.8f, -1.3f, ATOM_C);      // 8: C=O (ester 1)
    addAtom(mol, -1.5f, -4.0f, -1.0f, ATOM_O);      // 9: =O
    addAtom(mol, -2.3f, -2.0f, -1.5f, ATOM_O);      // 10: O-CH3
    addAtom(mol, -3.5f, -2.5f, -1.2f, ATOM_C);      // 11: CH3

    addAtom(mol, 1.2f, -3.0f, -1.5f, ATOM_C);       // 12: CH2
    addAtom(mol, 2.4f, -2.5f, -2.0f, ATOM_C);       // 13: C=O (ester 2)
    addAtom(mol, 2.6f, -1.3f, -2.3f, ATOM_O);       // 14: =O
    addAtom(mol, 3.4f, -3.3f, -2.0f, ATOM_O);       // 15: O-CH3
    addAtom(mol, 4.6f, -2.8f, -2.5f, ATOM_C);       // 16: CH3

    // Hydrogens on methoxy groups (simplified - 3 each)
    addAtom(mol, -3.0f, 0.8f, 0.0f, ATOM_H);        // 17
    addAtom(mol, -2.8f, -0.8f, 0.0f, ATOM_H);       // 18
    addAtom(mol, -2.8f, 0.0f, 1.5f, ATOM_H);        // 19
    addAtom(mol, 3.0f, 0.8f, 0.0f, ATOM_H);         // 20
    addAtom(mol, 2.8f, -0.8f, 0.0f, ATOM_H);        // 21
    addAtom(mol, 2.8f, 0.0f, 1.5f, ATOM_H);         // 22

    // H on CH
    addAtom(mol, 0.0f, -2.2f, -2.8f, ATOM_H);       // 23

    // H's on CH2
    addAtom(mol, 1.0f, -3.8f, -0.8f, ATOM_H);       // 24
    addAtom(mol, 1.4f, -3.5f, -2.3f, ATOM_H);       // 25

    // Bonds
    addBond(mol, 0, 1, 2);   // P=S
    addBond(mol, 0, 2, 1);   // P-O
    addBond(mol, 0, 3, 1);   // P-O
    addBond(mol, 0, 4, 1);   // P-S
    addBond(mol, 2, 5, 1);   // O-CH3
    addBond(mol, 3, 6, 1);   // O-CH3
    addBond(mol, 4, 7, 1);   // S-CH
    addBond(mol, 7, 8, 1);   // CH-C=O
    addBond(mol, 8, 9, 2);   // C=O
    addBond(mol, 8, 10, 1);  // C-O
    addBond(mol, 10, 11, 1); // O-CH3
    addBond(mol, 7, 12, 1);  // CH-CH2
    addBond(mol, 12, 13, 1); // CH2-C=O
    addBond(mol, 13, 14, 2); // C=O
    addBond(mol, 13, 15, 1); // C-O
    addBond(mol, 15, 16, 1); // O-CH3

    // Hydrogen bonds
    addBond(mol, 5, 17, 1);
    addBond(mol, 5, 18, 1);
    addBond(mol, 5, 19, 1);
    addBond(mol, 6, 20, 1);
    addBond(mol, 6, 21, 1);
    addBond(mol, 6, 22, 1);
    addBond(mol, 7, 23, 1);
    addBond(mol, 12, 24, 1);
    addBond(mol, 12, 25, 1);

    centerMolecule(mol);
}

// ============== FRAGRANCES/FLAVORS ==============

// Build Vanillin (C8H8O3) - Vanilla flavor/fragrance
void buildVanillin(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Vanillin (C8H8O3)");

    // Benzene ring (phenyl)
    float r = 1.4f;
    for (int i = 0; i < 6; i++) {
        float angle = i * PI / 3.0f;
        addAtom(mol, r * cosf(angle), r * sinf(angle), 0.0f, ATOM_C);  // 0-5
    }

    // Substituents on benzene:
    // Position 1 (C0): -CHO (aldehyde)
    // Position 3 (C2): -OCH3 (methoxy)
    // Position 4 (C3): -OH (hydroxyl)

    // Aldehyde group at C0
    addAtom(mol, r + 1.2f, 0.0f, 0.0f, ATOM_C);         // 6: CHO carbon
    addAtom(mol, r + 1.8f, 1.0f, 0.0f, ATOM_O);         // 7: =O
    addAtom(mol, r + 1.8f, -0.9f, 0.0f, ATOM_H);        // 8: H on CHO

    // Methoxy group at C2 (meta to aldehyde)
    addAtom(mol, r * cosf(2*PI/3) - 1.2f, r * sinf(2*PI/3), 0.0f, ATOM_O);  // 9: O
    addAtom(mol, r * cosf(2*PI/3) - 2.4f, r * sinf(2*PI/3), 0.0f, ATOM_C);  // 10: CH3

    // Hydroxyl at C3 (para to aldehyde)
    addAtom(mol, r * cosf(PI) - 1.2f, 0.0f, 0.0f, ATOM_O);    // 11: OH oxygen
    addAtom(mol, r * cosf(PI) - 1.8f, 0.8f, 0.0f, ATOM_H);    // 12: H on OH

    // Hydrogens on benzene (positions 2, 5, 6 / C1, C4, C5)
    float rH = 2.4f;
    addAtom(mol, rH * cosf(PI/3), rH * sinf(PI/3), 0.0f, ATOM_H);       // 13: H on C1
    addAtom(mol, rH * cosf(4*PI/3), rH * sinf(4*PI/3), 0.0f, ATOM_H);   // 14: H on C4
    addAtom(mol, rH * cosf(5*PI/3), rH * sinf(5*PI/3), 0.0f, ATOM_H);   // 15: H on C5

    // Methyl hydrogens
    addAtom(mol, r * cosf(2*PI/3) - 2.8f, r * sinf(2*PI/3) + 0.9f, 0.0f, ATOM_H);   // 16
    addAtom(mol, r * cosf(2*PI/3) - 2.8f, r * sinf(2*PI/3) - 0.5f, 0.8f, ATOM_H);   // 17
    addAtom(mol, r * cosf(2*PI/3) - 2.8f, r * sinf(2*PI/3) - 0.5f, -0.8f, ATOM_H);  // 18

    // Benzene ring bonds (alternating single/double)
    for (int i = 0; i < 6; i++) {
        addBond(mol, i, (i + 1) % 6, (i % 2 == 0) ? 2 : 1);
    }

    // Aldehyde bonds
    addBond(mol, 0, 6, 1);   // C-CHO
    addBond(mol, 6, 7, 2);   // C=O
    addBond(mol, 6, 8, 1);   // C-H

    // Methoxy bonds
    addBond(mol, 2, 9, 1);   // C-O
    addBond(mol, 9, 10, 1);  // O-CH3

    // Hydroxyl bonds
    addBond(mol, 3, 11, 1);  // C-O
    addBond(mol, 11, 12, 1); // O-H

    // Benzene hydrogen bonds
    addBond(mol, 1, 13, 1);
    addBond(mol, 4, 14, 1);
    addBond(mol, 5, 15, 1);

    // Methyl hydrogen bonds
    addBond(mol, 10, 16, 1);
    addBond(mol, 10, 17, 1);
    addBond(mol, 10, 18, 1);

    centerMolecule(mol);
}

// Build Limonene (C10H16) - Citrus scent (D-limonene shown)
void buildLimonene(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Limonene (C10H16)");

    // Cyclohexene ring with isopropenyl substituent
    // Ring carbons (chair-like conformation)
    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_C);         // 0: C1 (double bond)
    addAtom(mol, 1.3f, 0.5f, 0.3f, ATOM_C);         // 1: C2 (double bond)
    addAtom(mol, 2.3f, -0.5f, 0.0f, ATOM_C);        // 2: C3
    addAtom(mol, 2.0f, -1.8f, -0.5f, ATOM_C);       // 3: C4 (chiral, with substituent)
    addAtom(mol, 0.6f, -2.0f, -0.3f, ATOM_C);       // 4: C5
    addAtom(mol, -0.4f, -1.0f, 0.2f, ATOM_C);       // 5: C6 (methyl attached)

    // Methyl group on C6
    addAtom(mol, -1.8f, -1.3f, 0.5f, ATOM_C);       // 6: CH3 on C6

    // Isopropenyl group on C4: =C(CH3)2 -> actually =CH-CH3 (prop-1-en-2-yl)
    addAtom(mol, 2.8f, -2.8f, -0.8f, ATOM_C);       // 7: =CH2 carbon
    addAtom(mol, 2.5f, -2.5f, -2.2f, ATOM_C);       // 8: CH3 on isopropenyl

    // Hydrogens on ring
    addAtom(mol, -0.5f, 0.8f, -0.5f, ATOM_H);       // 9: H on C1
    addAtom(mol, 1.5f, 1.5f, 0.6f, ATOM_H);         // 10: H on C2
    addAtom(mol, 3.2f, -0.2f, 0.5f, ATOM_H);        // 11: H on C3
    addAtom(mol, 2.8f, -0.6f, -0.9f, ATOM_H);       // 12: H on C3
    addAtom(mol, 2.3f, -1.8f, 0.5f, ATOM_H);        // 13: H on C4
    addAtom(mol, 0.3f, -2.2f, -1.3f, ATOM_H);       // 14: H on C5
    addAtom(mol, 0.4f, -2.9f, 0.3f, ATOM_H);        // 15: H on C5
    addAtom(mol, -0.1f, -0.8f, 1.2f, ATOM_H);       // 16: H on C6

    // Methyl hydrogens on C6-CH3
    addAtom(mol, -2.0f, -2.3f, 0.8f, ATOM_H);       // 17
    addAtom(mol, -2.4f, -0.9f, -0.3f, ATOM_H);      // 18
    addAtom(mol, -2.1f, -0.7f, 1.3f, ATOM_H);       // 19

    // Isopropenyl =CH2 hydrogens
    addAtom(mol, 3.8f, -2.6f, -0.5f, ATOM_H);       // 20
    addAtom(mol, 2.6f, -3.8f, -0.5f, ATOM_H);       // 21

    // Isopropenyl CH3 hydrogens
    addAtom(mol, 1.5f, -2.7f, -2.5f, ATOM_H);       // 22
    addAtom(mol, 3.2f, -3.2f, -2.6f, ATOM_H);       // 23
    addAtom(mol, 2.7f, -1.5f, -2.5f, ATOM_H);       // 24

    // Ring bonds
    addBond(mol, 0, 1, 2);   // C1=C2 double bond
    addBond(mol, 1, 2, 1);   // C2-C3
    addBond(mol, 2, 3, 1);   // C3-C4
    addBond(mol, 3, 4, 1);   // C4-C5
    addBond(mol, 4, 5, 1);   // C5-C6
    addBond(mol, 5, 0, 1);   // C6-C1

    // Substituent bonds
    addBond(mol, 5, 6, 1);   // C6-CH3
    addBond(mol, 3, 7, 1);   // C4-isopropenyl
    addBond(mol, 7, 8, 2);   // C=CH2 (actually shown as single for =C)

    // Ring hydrogen bonds
    addBond(mol, 0, 9, 1);
    addBond(mol, 1, 10, 1);
    addBond(mol, 2, 11, 1);
    addBond(mol, 2, 12, 1);
    addBond(mol, 3, 13, 1);
    addBond(mol, 4, 14, 1);
    addBond(mol, 4, 15, 1);
    addBond(mol, 5, 16, 1);

    // Methyl hydrogen bonds
    addBond(mol, 6, 17, 1);
    addBond(mol, 6, 18, 1);
    addBond(mol, 6, 19, 1);

    // Isopropenyl hydrogen bonds
    addBond(mol, 7, 20, 1);
    addBond(mol, 7, 21, 1);
    addBond(mol, 8, 22, 1);
    addBond(mol, 8, 23, 1);
    addBond(mol, 8, 24, 1);

    centerMolecule(mol);
}

// Build Menthol (C10H20O) - Mint/cooling sensation
void buildMenthol(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Menthol (C10H20O)");

    // Cyclohexane ring with three substituents (chair conformation)
    // C1-OH, C2-isopropyl, C5-methyl
    addAtom(mol, 0.0f, 0.0f, 0.5f, ATOM_C);         // 0: C1 (OH attached)
    addAtom(mol, 1.4f, 0.3f, 0.0f, ATOM_C);         // 1: C2 (isopropyl attached)
    addAtom(mol, 2.2f, -0.9f, 0.4f, ATOM_C);        // 2: C3
    addAtom(mol, 1.5f, -2.2f, 0.0f, ATOM_C);        // 3: C4
    addAtom(mol, 0.1f, -2.0f, 0.5f, ATOM_C);        // 4: C5 (methyl attached)
    addAtom(mol, -0.7f, -0.8f, 0.1f, ATOM_C);       // 5: C6

    // Hydroxyl on C1
    addAtom(mol, -0.5f, 1.2f, 0.0f, ATOM_O);        // 6: OH oxygen
    addAtom(mol, -0.1f, 1.9f, 0.5f, ATOM_H);        // 7: H on OH

    // Isopropyl on C2
    addAtom(mol, 2.0f, 1.6f, 0.4f, ATOM_C);         // 8: CH (isopropyl)
    addAtom(mol, 1.3f, 2.7f, -0.2f, ATOM_C);        // 9: CH3
    addAtom(mol, 3.4f, 1.8f, 0.0f, ATOM_C);         // 10: CH3

    // Methyl on C5
    addAtom(mol, -0.6f, -3.2f, 0.0f, ATOM_C);       // 11: CH3

    // Ring hydrogens (axial/equatorial)
    addAtom(mol, 0.0f, 0.0f, 1.6f, ATOM_H);         // 12: H on C1
    addAtom(mol, 1.4f, 0.3f, -1.1f, ATOM_H);        // 13: H on C2
    addAtom(mol, 3.2f, -0.8f, 0.0f, ATOM_H);        // 14: H on C3
    addAtom(mol, 2.3f, -1.0f, 1.5f, ATOM_H);        // 15: H on C3
    addAtom(mol, 1.5f, -2.3f, -1.1f, ATOM_H);       // 16: H on C4
    addAtom(mol, 2.0f, -3.1f, 0.4f, ATOM_H);        // 17: H on C4
    addAtom(mol, 0.1f, -2.0f, 1.6f, ATOM_H);        // 18: H on C5
    addAtom(mol, -1.7f, -0.7f, 0.5f, ATOM_H);       // 19: H on C6
    addAtom(mol, -0.8f, -0.8f, -1.0f, ATOM_H);      // 20: H on C6

    // Isopropyl CH hydrogen
    addAtom(mol, 1.9f, 1.7f, 1.5f, ATOM_H);         // 21: H on CH

    // Isopropyl CH3 hydrogens (first methyl)
    addAtom(mol, 0.3f, 2.5f, 0.0f, ATOM_H);         // 22
    addAtom(mol, 1.4f, 3.6f, 0.3f, ATOM_H);         // 23
    addAtom(mol, 1.5f, 2.8f, -1.3f, ATOM_H);        // 24

    // Isopropyl CH3 hydrogens (second methyl)
    addAtom(mol, 3.5f, 1.9f, -1.1f, ATOM_H);        // 25
    addAtom(mol, 3.9f, 2.6f, 0.5f, ATOM_H);         // 26
    addAtom(mol, 3.9f, 0.9f, 0.3f, ATOM_H);         // 27

    // C5 methyl hydrogens
    addAtom(mol, -1.6f, -3.1f, 0.3f, ATOM_H);       // 28
    addAtom(mol, -0.5f, -3.3f, -1.1f, ATOM_H);      // 29
    addAtom(mol, -0.3f, -4.1f, 0.5f, ATOM_H);       // 30

    // Ring bonds
    addBond(mol, 0, 1, 1);
    addBond(mol, 1, 2, 1);
    addBond(mol, 2, 3, 1);
    addBond(mol, 3, 4, 1);
    addBond(mol, 4, 5, 1);
    addBond(mol, 5, 0, 1);

    // Hydroxyl bonds
    addBond(mol, 0, 6, 1);   // C-O
    addBond(mol, 6, 7, 1);   // O-H

    // Isopropyl bonds
    addBond(mol, 1, 8, 1);   // C-CH
    addBond(mol, 8, 9, 1);   // CH-CH3
    addBond(mol, 8, 10, 1);  // CH-CH3

    // Methyl bond
    addBond(mol, 4, 11, 1);  // C-CH3

    // Ring hydrogen bonds
    addBond(mol, 0, 12, 1);
    addBond(mol, 1, 13, 1);
    addBond(mol, 2, 14, 1);
    addBond(mol, 2, 15, 1);
    addBond(mol, 3, 16, 1);
    addBond(mol, 3, 17, 1);
    addBond(mol, 4, 18, 1);
    addBond(mol, 5, 19, 1);
    addBond(mol, 5, 20, 1);

    // Isopropyl hydrogen bonds
    addBond(mol, 8, 21, 1);
    addBond(mol, 9, 22, 1);
    addBond(mol, 9, 23, 1);
    addBond(mol, 9, 24, 1);
    addBond(mol, 10, 25, 1);
    addBond(mol, 10, 26, 1);
    addBond(mol, 10, 27, 1);

    // C5 methyl hydrogen bonds
    addBond(mol, 11, 28, 1);
    addBond(mol, 11, 29, 1);
    addBond(mol, 11, 30, 1);

    centerMolecule(mol);
}

// Build Cinnamaldehyde (C9H8O) - Cinnamon flavor/fragrance
void buildCinnamaldehyde(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Cinnamaldehyde (C9H8O)");

    // Benzene ring
    float r = 1.4f;
    for (int i = 0; i < 6; i++) {
        float angle = i * PI / 3.0f;
        addAtom(mol, r * cosf(angle), r * sinf(angle), 0.0f, ATOM_C);  // 0-5
    }

    // Propenal chain: -CH=CH-CHO (trans configuration)
    // Attached at C0 (position 0)
    addAtom(mol, r + 1.3f, 0.0f, 0.0f, ATOM_C);         // 6: =CH (alpha)
    addAtom(mol, r + 2.5f, 0.7f, 0.0f, ATOM_C);         // 7: =CH (beta)
    addAtom(mol, r + 3.8f, 0.2f, 0.0f, ATOM_C);         // 8: CHO
    addAtom(mol, r + 4.8f, 0.9f, 0.0f, ATOM_O);         // 9: =O

    // Hydrogens on benzene ring
    float rH = 2.4f;
    addAtom(mol, rH * cosf(PI/3), rH * sinf(PI/3), 0.0f, ATOM_H);       // 10: H on C1
    addAtom(mol, rH * cosf(2*PI/3), rH * sinf(2*PI/3), 0.0f, ATOM_H);   // 11: H on C2
    addAtom(mol, rH * cosf(PI), rH * sinf(PI), 0.0f, ATOM_H);           // 12: H on C3
    addAtom(mol, rH * cosf(4*PI/3), rH * sinf(4*PI/3), 0.0f, ATOM_H);   // 13: H on C4
    addAtom(mol, rH * cosf(5*PI/3), rH * sinf(5*PI/3), 0.0f, ATOM_H);   // 14: H on C5

    // Hydrogens on propenal chain
    addAtom(mol, r + 1.3f, -1.0f, 0.0f, ATOM_H);        // 15: H on alpha carbon
    addAtom(mol, r + 2.5f, 1.7f, 0.0f, ATOM_H);         // 16: H on beta carbon
    addAtom(mol, r + 3.9f, -0.8f, 0.0f, ATOM_H);        // 17: H on CHO

    // Benzene ring bonds (alternating single/double for resonance)
    for (int i = 0; i < 6; i++) {
        addBond(mol, i, (i + 1) % 6, (i % 2 == 0) ? 2 : 1);
    }

    // Propenal chain bonds
    addBond(mol, 0, 6, 1);   // C(ring)-CH
    addBond(mol, 6, 7, 2);   // CH=CH (double bond)
    addBond(mol, 7, 8, 1);   // CH-CHO
    addBond(mol, 8, 9, 2);   // C=O (aldehyde)

    // Benzene hydrogen bonds
    addBond(mol, 1, 10, 1);
    addBond(mol, 2, 11, 1);
    addBond(mol, 3, 12, 1);
    addBond(mol, 4, 13, 1);
    addBond(mol, 5, 14, 1);

    // Chain hydrogen bonds
    addBond(mol, 6, 15, 1);
    addBond(mol, 7, 16, 1);
    addBond(mol, 8, 17, 1);

    centerMolecule(mol);
}

// Build syn-Propanethial-S-oxide (C3H6OS) - Onion lachrymatory factor (makes you cry)
void buildPropanethialSoxide(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Propanethial-S-oxide (C3H6OS)");

    // Structure: CH3-CH2-CH=S=O (syn configuration)
    // The molecule has a thial S-oxide group (-CH=S(O)-)

    // Carbon chain
    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_C);         // 0: CH3 (methyl)
    addAtom(mol, 1.5f, 0.0f, 0.0f, ATOM_C);         // 1: CH2 (methylene)
    addAtom(mol, 2.5f, 1.2f, 0.0f, ATOM_C);         // 2: CH (thial carbon)

    // Sulfoxide group (S=O)
    addAtom(mol, 3.8f, 1.0f, 0.5f, ATOM_S);         // 3: S (sulfur)
    addAtom(mol, 4.5f, 2.2f, 0.8f, ATOM_O);         // 4: O (oxide)

    // Hydrogens on CH3
    addAtom(mol, -0.5f, 0.9f, 0.3f, ATOM_H);        // 5
    addAtom(mol, -0.5f, -0.7f, 0.6f, ATOM_H);       // 6
    addAtom(mol, -0.3f, -0.2f, -1.0f, ATOM_H);      // 7

    // Hydrogens on CH2
    addAtom(mol, 1.7f, -0.5f, -0.9f, ATOM_H);       // 8
    addAtom(mol, 1.7f, -0.7f, 0.8f, ATOM_H);        // 9

    // Hydrogen on thial CH
    addAtom(mol, 2.2f, 2.2f, -0.3f, ATOM_H);        // 10

    // Bonds
    addBond(mol, 0, 1, 1);   // CH3-CH2
    addBond(mol, 1, 2, 1);   // CH2-CH
    addBond(mol, 2, 3, 2);   // CH=S (double bond)
    addBond(mol, 3, 4, 2);   // S=O (double bond)

    // CH3 hydrogen bonds
    addBond(mol, 0, 5, 1);
    addBond(mol, 0, 6, 1);
    addBond(mol, 0, 7, 1);

    // CH2 hydrogen bonds
    addBond(mol, 1, 8, 1);
    addBond(mol, 1, 9, 1);

    // CH hydrogen bond
    addBond(mol, 2, 10, 1);

    centerMolecule(mol);
}

// ============== ENVIRONMENTAL/CLIMATE MOLECULES ==============

// Build CFC-12 / Freon-12 (CCl2F2) - Ozone-depleting refrigerant
void buildCFC12(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "CFC-12/Freon-12 (CCl2F2)");

    // Tetrahedral carbon center with 2 Cl and 2 F
    // Central carbon
    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_C);         // 0: C

    // Two chlorine atoms (larger, green)
    addAtom(mol, 1.2f, 1.0f, 0.5f, ATOM_CL);        // 1: Cl
    addAtom(mol, -1.2f, 1.0f, 0.5f, ATOM_CL);       // 2: Cl

    // Two fluorine atoms (smaller, light green)
    addAtom(mol, 0.6f, -1.0f, -0.8f, ATOM_F);       // 3: F
    addAtom(mol, -0.6f, -0.5f, 1.0f, ATOM_F);       // 4: F

    // Bonds (all single bonds to central carbon)
    addBond(mol, 0, 1, 1);   // C-Cl
    addBond(mol, 0, 2, 1);   // C-Cl
    addBond(mol, 0, 3, 1);   // C-F
    addBond(mol, 0, 4, 1);   // C-F

    centerMolecule(mol);
}

// Build SF6 - Sulfur Hexafluoride (most potent greenhouse gas)
void buildSF6(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Sulfur Hexafluoride (SF6)");

    // Octahedral geometry - sulfur at center, 6 fluorines at vertices
    // Central sulfur
    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_S);         // 0: S

    // Six fluorine atoms in octahedral arrangement
    float d = 1.6f;  // S-F bond length
    addAtom(mol, d, 0.0f, 0.0f, ATOM_F);            // 1: F (+x)
    addAtom(mol, -d, 0.0f, 0.0f, ATOM_F);           // 2: F (-x)
    addAtom(mol, 0.0f, d, 0.0f, ATOM_F);            // 3: F (+y)
    addAtom(mol, 0.0f, -d, 0.0f, ATOM_F);           // 4: F (-y)
    addAtom(mol, 0.0f, 0.0f, d, ATOM_F);            // 5: F (+z)
    addAtom(mol, 0.0f, 0.0f, -d, ATOM_F);           // 6: F (-z)

    // Bonds (all single bonds from sulfur to fluorines)
    addBond(mol, 0, 1, 1);   // S-F
    addBond(mol, 0, 2, 1);   // S-F
    addBond(mol, 0, 3, 1);   // S-F
    addBond(mol, 0, 4, 1);   // S-F
    addBond(mol, 0, 5, 1);   // S-F
    addBond(mol, 0, 6, 1);   // S-F

    centerMolecule(mol);
}

// ============== ANESTHETICS ==============

// Build Lidocaine (C14H22N2O) - Local anesthetic (Xylocaine)
void buildLidocaine(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Lidocaine (C14H22N2O)");

    // Structure: 2,6-dimethylaniline linked via amide to diethylaminoethyl
    // Benzene ring (2,6-dimethylaniline part)
    float r = 1.4f;
    for (int i = 0; i < 6; i++) {
        float angle = i * PI / 3.0f;
        addAtom(mol, r * cosf(angle), r * sinf(angle), 0.0f, ATOM_C);  // 0-5
    }

    // Methyl groups at positions 2 and 6 (ortho to NH)
    addAtom(mol, r * cosf(PI/3) + 1.0f, r * sinf(PI/3) + 0.5f, 0.0f, ATOM_C);   // 6: CH3
    addAtom(mol, r * cosf(5*PI/3) + 1.0f, r * sinf(5*PI/3) - 0.5f, 0.0f, ATOM_C); // 7: CH3

    // Amide linkage: -NH-C(=O)-CH2-
    addAtom(mol, r + 1.3f, 0.0f, 0.0f, ATOM_N);         // 8: NH (amide)
    addAtom(mol, r + 2.5f, 0.0f, 0.0f, ATOM_C);         // 9: C=O
    addAtom(mol, r + 3.0f, 1.0f, 0.0f, ATOM_O);         // 10: =O
    addAtom(mol, r + 3.5f, -0.8f, 0.0f, ATOM_C);        // 11: CH2

    // Diethylamino group: -N(CH2CH3)2
    addAtom(mol, r + 4.8f, -0.5f, 0.0f, ATOM_N);        // 12: N tertiary
    addAtom(mol, r + 5.5f, 0.7f, 0.5f, ATOM_C);         // 13: CH2 (ethyl 1)
    addAtom(mol, r + 6.8f, 0.5f, 0.8f, ATOM_C);         // 14: CH3 (ethyl 1)
    addAtom(mol, r + 5.5f, -1.5f, -0.5f, ATOM_C);       // 15: CH2 (ethyl 2)
    addAtom(mol, r + 6.8f, -1.8f, -0.8f, ATOM_C);       // 16: CH3 (ethyl 2)

    // Key hydrogens (simplified)
    addAtom(mol, r + 1.5f, 0.0f, 0.9f, ATOM_H);         // 17: H on NH
    float rH = 2.4f;
    addAtom(mol, rH * cosf(2*PI/3), rH * sinf(2*PI/3), 0.0f, ATOM_H);   // 18: H on C2
    addAtom(mol, rH * cosf(PI), rH * sinf(PI), 0.0f, ATOM_H);           // 19: H on C3
    addAtom(mol, rH * cosf(4*PI/3), rH * sinf(4*PI/3), 0.0f, ATOM_H);   // 20: H on C4

    // Benzene ring bonds
    for (int i = 0; i < 6; i++) {
        addBond(mol, i, (i + 1) % 6, (i % 2 == 0) ? 2 : 1);
    }

    // Methyl bonds
    addBond(mol, 1, 6, 1);   // C-CH3
    addBond(mol, 5, 7, 1);   // C-CH3

    // Amide linkage bonds
    addBond(mol, 0, 8, 1);   // C(ring)-N
    addBond(mol, 8, 9, 1);   // N-C=O
    addBond(mol, 9, 10, 2);  // C=O
    addBond(mol, 9, 11, 1);  // C-CH2

    // Diethylamino bonds
    addBond(mol, 11, 12, 1); // CH2-N
    addBond(mol, 12, 13, 1); // N-CH2
    addBond(mol, 13, 14, 1); // CH2-CH3
    addBond(mol, 12, 15, 1); // N-CH2
    addBond(mol, 15, 16, 1); // CH2-CH3

    // Hydrogen bonds
    addBond(mol, 8, 17, 1);
    addBond(mol, 2, 18, 1);
    addBond(mol, 3, 19, 1);
    addBond(mol, 4, 20, 1);

    centerMolecule(mol);
}

// Build Ketamine (C13H16ClNO) - Dissociative anesthetic
void buildKetamine(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Ketamine (C13H16ClNO)");

    // Structure: 2-(2-chlorophenyl)-2-(methylamino)cyclohexanone
    // Cyclohexanone ring
    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_C);         // 0: C2 (chiral, with substituents)
    addAtom(mol, 1.3f, 0.5f, 0.5f, ATOM_C);         // 1: C3
    addAtom(mol, 2.3f, -0.5f, 0.3f, ATOM_C);        // 2: C4
    addAtom(mol, 2.0f, -1.8f, -0.2f, ATOM_C);       // 3: C5
    addAtom(mol, 0.6f, -2.0f, -0.5f, ATOM_C);       // 4: C6
    addAtom(mol, -0.4f, -1.0f, -0.3f, ATOM_C);      // 5: C1 (ketone)
    addAtom(mol, -1.6f, -1.2f, -0.5f, ATOM_O);      // 6: =O (ketone)

    // 2-chlorophenyl ring attached to C2
    float r = 1.4f;
    float baseX = -0.5f, baseY = 1.5f;
    for (int i = 0; i < 6; i++) {
        float angle = i * PI / 3.0f + PI/2;
        addAtom(mol, baseX + r * cosf(angle), baseY + r * sinf(angle), 0.0f, ATOM_C);  // 7-12
    }

    // Chlorine at ortho position (C2 of phenyl)
    addAtom(mol, baseX + r * cosf(PI/2 + PI/3) - 0.8f, baseY + r * sinf(PI/2 + PI/3) + 0.8f, 0.0f, ATOM_CL); // 13: Cl

    // Methylamino group on C2 of cyclohexanone
    addAtom(mol, 0.5f, 0.5f, -1.3f, ATOM_N);        // 14: NH
    addAtom(mol, 0.3f, 1.5f, -2.0f, ATOM_C);        // 15: CH3

    // Key hydrogens
    addAtom(mol, 1.0f, 0.0f, -1.8f, ATOM_H);        // 16: H on N
    float rH = 2.4f;
    addAtom(mol, baseX + rH * cosf(PI/2 + 2*PI/3), baseY + rH * sinf(PI/2 + 2*PI/3), 0.0f, ATOM_H);  // 17
    addAtom(mol, baseX + rH * cosf(PI/2 + PI), baseY + rH * sinf(PI/2 + PI), 0.0f, ATOM_H);          // 18
    addAtom(mol, baseX + rH * cosf(PI/2 + 4*PI/3), baseY + rH * sinf(PI/2 + 4*PI/3), 0.0f, ATOM_H);  // 19
    addAtom(mol, baseX + rH * cosf(PI/2 + 5*PI/3), baseY + rH * sinf(PI/2 + 5*PI/3), 0.0f, ATOM_H);  // 20

    // Cyclohexanone ring bonds
    addBond(mol, 0, 1, 1);
    addBond(mol, 1, 2, 1);
    addBond(mol, 2, 3, 1);
    addBond(mol, 3, 4, 1);
    addBond(mol, 4, 5, 1);
    addBond(mol, 5, 0, 1);
    addBond(mol, 5, 6, 2);   // C=O

    // Phenyl ring bonds
    for (int i = 0; i < 6; i++) {
        addBond(mol, 7 + i, 7 + ((i + 1) % 6), (i % 2 == 0) ? 2 : 1);
    }

    // Connections
    addBond(mol, 0, 7, 1);   // cyclohexanone-phenyl
    addBond(mol, 8, 13, 1);  // C-Cl
    addBond(mol, 0, 14, 1);  // C-N
    addBond(mol, 14, 15, 1); // N-CH3

    // Hydrogen bonds
    addBond(mol, 14, 16, 1);
    addBond(mol, 9, 17, 1);
    addBond(mol, 10, 18, 1);
    addBond(mol, 11, 19, 1);
    addBond(mol, 12, 20, 1);

    centerMolecule(mol);
}

// Build Sevoflurane (C4H3F7O) - Inhaled general anesthetic
void buildSevoflurane(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Sevoflurane (C4H3F7O)");

    // Structure: (CF3)2CH-O-CH2F (fluoromethyl 2,2,2-trifluoro-1-[trifluoromethyl]ethyl ether)
    // Central CH connected to two CF3 groups and ether oxygen
    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_C);         // 0: CH central

    // First CF3 group
    addAtom(mol, -1.3f, 0.8f, 0.0f, ATOM_C);        // 1: CF3
    addAtom(mol, -1.8f, 1.3f, 1.0f, ATOM_F);        // 2: F
    addAtom(mol, -2.2f, 0.3f, -0.5f, ATOM_F);       // 3: F
    addAtom(mol, -1.0f, 1.8f, -0.7f, ATOM_F);       // 4: F

    // Second CF3 group
    addAtom(mol, -0.3f, -1.5f, 0.3f, ATOM_C);       // 5: CF3
    addAtom(mol, -1.3f, -2.0f, -0.3f, ATOM_F);      // 6: F
    addAtom(mol, 0.5f, -2.3f, -0.2f, ATOM_F);       // 7: F
    addAtom(mol, -0.5f, -1.7f, 1.5f, ATOM_F);       // 8: F

    // Ether linkage: -O-CH2F
    addAtom(mol, 1.3f, 0.5f, -0.3f, ATOM_O);        // 9: O (ether)
    addAtom(mol, 2.5f, 0.2f, 0.2f, ATOM_C);         // 10: CH2F
    addAtom(mol, 3.5f, 0.8f, -0.3f, ATOM_F);        // 11: F on CH2F

    // Hydrogens
    addAtom(mol, 0.3f, 0.3f, 1.0f, ATOM_H);         // 12: H on central CH
    addAtom(mol, 2.6f, -0.8f, 0.0f, ATOM_H);        // 13: H on CH2F
    addAtom(mol, 2.6f, 0.5f, 1.2f, ATOM_H);         // 14: H on CH2F

    // Bonds
    addBond(mol, 0, 1, 1);   // CH-CF3
    addBond(mol, 1, 2, 1);   // C-F
    addBond(mol, 1, 3, 1);   // C-F
    addBond(mol, 1, 4, 1);   // C-F

    addBond(mol, 0, 5, 1);   // CH-CF3
    addBond(mol, 5, 6, 1);   // C-F
    addBond(mol, 5, 7, 1);   // C-F
    addBond(mol, 5, 8, 1);   // C-F

    addBond(mol, 0, 9, 1);   // CH-O
    addBond(mol, 9, 10, 1);  // O-CH2F
    addBond(mol, 10, 11, 1); // C-F

    // Hydrogen bonds
    addBond(mol, 0, 12, 1);
    addBond(mol, 10, 13, 1);
    addBond(mol, 10, 14, 1);

    centerMolecule(mol);
}

// ============== PSYCHEDELICS ==============

// Build LSD (C20H25N3O) - Lysergic acid diethylamide
void buildLSD(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "LSD (C20H25N3O)");

    // Ergoline core - tetracyclic structure (simplified representation)
    // Ring A: Benzene ring
    float r = 1.4f;
    for (int i = 0; i < 6; i++) {
        float angle = i * PI / 3.0f;
        addAtom(mol, r * cosf(angle), r * sinf(angle), 0.0f, ATOM_C);  // 0-5
    }

    // Ring B: Pyrrole (5-membered with N)
    addAtom(mol, r * cosf(0) + 1.2f, r * sinf(0) + 0.7f, 0.3f, ATOM_C);     // 6
    addAtom(mol, r * cosf(0) + 2.0f, r * sinf(0) - 0.3f, 0.5f, ATOM_C);     // 7
    addAtom(mol, r * cosf(5*PI/3) + 1.0f, r * sinf(5*PI/3) - 0.3f, 0.3f, ATOM_N); // 8: indole N

    // Ring C: Cyclohexene
    addAtom(mol, r * cosf(0) + 2.8f, r * sinf(0) + 0.5f, 0.2f, ATOM_C);     // 9
    addAtom(mol, r * cosf(0) + 3.5f, r * sinf(0) - 0.5f, -0.3f, ATOM_C);    // 10
    addAtom(mol, r * cosf(0) + 3.0f, r * sinf(0) - 1.5f, -0.5f, ATOM_C);    // 11

    // Ring D: Piperidine with N
    addAtom(mol, r * cosf(0) + 2.0f, r * sinf(0) - 1.8f, -0.2f, ATOM_N);    // 12: piperidine N
    addAtom(mol, r * cosf(0) + 1.5f, r * sinf(0) - 2.8f, 0.0f, ATOM_C);     // 13: N-CH3 methyl

    // Diethylamide group: -C(=O)-N(C2H5)2
    addAtom(mol, r * cosf(0) + 4.5f, r * sinf(0) + 1.0f, 0.0f, ATOM_C);     // 14: C=O
    addAtom(mol, r * cosf(0) + 4.8f, r * sinf(0) + 2.2f, 0.3f, ATOM_O);     // 15: =O
    addAtom(mol, r * cosf(0) + 5.5f, r * sinf(0) + 0.2f, -0.3f, ATOM_N);    // 16: N (amide)
    addAtom(mol, r * cosf(0) + 6.5f, r * sinf(0) + 0.8f, 0.3f, ATOM_C);     // 17: CH2
    addAtom(mol, r * cosf(0) + 7.5f, r * sinf(0) + 0.2f, 0.0f, ATOM_C);     // 18: CH3
    addAtom(mol, r * cosf(0) + 5.8f, r * sinf(0) - 1.0f, -0.8f, ATOM_C);    // 19: CH2
    addAtom(mol, r * cosf(0) + 6.8f, r * sinf(0) - 1.5f, -0.5f, ATOM_C);    // 20: CH3

    // Key hydrogens
    addAtom(mol, r * cosf(5*PI/3) + 1.3f, r * sinf(5*PI/3) - 1.0f, 0.5f, ATOM_H); // 21: H on indole N
    float rH = 2.4f;
    addAtom(mol, rH * cosf(2*PI/3), rH * sinf(2*PI/3), 0.0f, ATOM_H);       // 22
    addAtom(mol, rH * cosf(PI), rH * sinf(PI), 0.0f, ATOM_H);               // 23
    addAtom(mol, rH * cosf(4*PI/3), rH * sinf(4*PI/3), 0.0f, ATOM_H);       // 24

    // Ring A bonds (benzene)
    for (int i = 0; i < 6; i++) {
        addBond(mol, i, (i + 1) % 6, (i % 2 == 0) ? 2 : 1);
    }

    // Ring B bonds (pyrrole fusion)
    addBond(mol, 0, 6, 1);
    addBond(mol, 6, 7, 2);
    addBond(mol, 7, 8, 1);
    addBond(mol, 8, 5, 1);

    // Ring C bonds
    addBond(mol, 6, 9, 1);
    addBond(mol, 9, 10, 1);
    addBond(mol, 10, 11, 1);
    addBond(mol, 11, 7, 1);

    // Ring D bonds
    addBond(mol, 11, 12, 1);
    addBond(mol, 12, 8, 1);
    addBond(mol, 12, 13, 1);  // N-CH3

    // Amide bonds
    addBond(mol, 9, 14, 1);
    addBond(mol, 14, 15, 2);  // C=O
    addBond(mol, 14, 16, 1);  // C-N
    addBond(mol, 16, 17, 1);
    addBond(mol, 17, 18, 1);
    addBond(mol, 16, 19, 1);
    addBond(mol, 19, 20, 1);

    // Hydrogen bonds
    addBond(mol, 8, 21, 1);
    addBond(mol, 2, 22, 1);
    addBond(mol, 3, 23, 1);
    addBond(mol, 4, 24, 1);

    centerMolecule(mol);
}

// Build Psilocybin (C12H17N2O4P) - Magic mushroom compound
void buildPsilocybin(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Psilocybin (C12H17N2O4P)");

    // Indole core (benzene fused with pyrrole)
    // Benzene ring
    float r = 1.4f;
    for (int i = 0; i < 6; i++) {
        float angle = i * PI / 3.0f;
        addAtom(mol, r * cosf(angle), r * sinf(angle), 0.0f, ATOM_C);  // 0-5
    }

    // Pyrrole ring (5-membered, fused at C0-C5)
    addAtom(mol, r * cosf(0) + 1.0f, r * sinf(0) + 0.8f, 0.0f, ATOM_C);     // 6: C3
    addAtom(mol, r * cosf(0) + 1.8f, r * sinf(0) - 0.2f, 0.0f, ATOM_C);     // 7: C2
    addAtom(mol, r * cosf(5*PI/3) + 0.8f, r * sinf(5*PI/3) - 0.3f, 0.0f, ATOM_N); // 8: N1 (indole)

    // Phosphate ester at position 4 (on C6)
    addAtom(mol, r * cosf(0) + 0.8f, r * sinf(0) + 2.0f, 0.0f, ATOM_O);     // 9: O (ether to P)
    addAtom(mol, r * cosf(0) + 1.5f, r * sinf(0) + 3.2f, 0.0f, ATOM_P);     // 10: P
    addAtom(mol, r * cosf(0) + 0.5f, r * sinf(0) + 4.2f, 0.0f, ATOM_O);     // 11: =O
    addAtom(mol, r * cosf(0) + 2.5f, r * sinf(0) + 3.8f, 0.8f, ATOM_O);     // 12: OH
    addAtom(mol, r * cosf(0) + 2.2f, r * sinf(0) + 3.0f, -1.0f, ATOM_O);    // 13: OH

    // Ethylamine side chain with dimethylamino: -CH2-CH2-N(CH3)2
    addAtom(mol, r * cosf(0) + 2.8f, r * sinf(0) + 0.3f, 0.0f, ATOM_C);     // 14: CH2
    addAtom(mol, r * cosf(0) + 4.0f, r * sinf(0) - 0.3f, 0.0f, ATOM_C);     // 15: CH2
    addAtom(mol, r * cosf(0) + 5.2f, r * sinf(0) + 0.3f, 0.0f, ATOM_N);     // 16: N (dimethyl)
    addAtom(mol, r * cosf(0) + 5.8f, r * sinf(0) + 1.5f, 0.0f, ATOM_C);     // 17: CH3
    addAtom(mol, r * cosf(0) + 6.2f, r * sinf(0) - 0.5f, 0.0f, ATOM_C);     // 18: CH3

    // Key hydrogens
    addAtom(mol, r * cosf(5*PI/3) + 1.0f, r * sinf(5*PI/3) - 1.2f, 0.0f, ATOM_H); // 19: H on indole N
    addAtom(mol, r * cosf(0) + 3.0f, r * sinf(0) + 4.5f, 1.0f, ATOM_H);     // 20: H on P-OH
    addAtom(mol, r * cosf(0) + 2.8f, r * sinf(0) + 2.5f, -1.5f, ATOM_H);    // 21: H on P-OH
    float rH = 2.4f;
    addAtom(mol, rH * cosf(PI/3), rH * sinf(PI/3), 0.0f, ATOM_H);           // 22
    addAtom(mol, rH * cosf(2*PI/3), rH * sinf(2*PI/3), 0.0f, ATOM_H);       // 23
    addAtom(mol, rH * cosf(PI), rH * sinf(PI), 0.0f, ATOM_H);               // 24
    addAtom(mol, rH * cosf(4*PI/3), rH * sinf(4*PI/3), 0.0f, ATOM_H);       // 25

    // Benzene ring bonds
    for (int i = 0; i < 6; i++) {
        addBond(mol, i, (i + 1) % 6, (i % 2 == 0) ? 2 : 1);
    }

    // Pyrrole fusion bonds
    addBond(mol, 0, 6, 1);
    addBond(mol, 6, 7, 2);
    addBond(mol, 7, 8, 1);
    addBond(mol, 8, 5, 1);

    // Phosphate ester bonds
    addBond(mol, 6, 9, 1);   // C-O
    addBond(mol, 9, 10, 1);  // O-P
    addBond(mol, 10, 11, 2); // P=O
    addBond(mol, 10, 12, 1); // P-OH
    addBond(mol, 10, 13, 1); // P-OH

    // Side chain bonds
    addBond(mol, 7, 14, 1);  // C-CH2
    addBond(mol, 14, 15, 1); // CH2-CH2
    addBond(mol, 15, 16, 1); // CH2-N
    addBond(mol, 16, 17, 1); // N-CH3
    addBond(mol, 16, 18, 1); // N-CH3

    // Hydrogen bonds
    addBond(mol, 8, 19, 1);
    addBond(mol, 12, 20, 1);
    addBond(mol, 13, 21, 1);
    addBond(mol, 1, 22, 1);
    addBond(mol, 2, 23, 1);
    addBond(mol, 3, 24, 1);
    addBond(mol, 4, 25, 1);

    centerMolecule(mol);
}

// Build Mescaline (C11H17NO3) - Peyote cactus alkaloid
void buildMescaline(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Mescaline (C11H17NO3)");

    // 3,4,5-trimethoxyphenethylamine
    // Benzene ring
    float r = 1.4f;
    for (int i = 0; i < 6; i++) {
        float angle = i * PI / 3.0f;
        addAtom(mol, r * cosf(angle), r * sinf(angle), 0.0f, ATOM_C);  // 0-5
    }

    // Three methoxy groups at positions 3, 4, 5 (C2, C3, C4)
    // Methoxy at C2 (position 3)
    addAtom(mol, r * cosf(2*PI/3) - 1.0f, r * sinf(2*PI/3) + 0.5f, 0.0f, ATOM_O);   // 6: O
    addAtom(mol, r * cosf(2*PI/3) - 2.2f, r * sinf(2*PI/3) + 0.3f, 0.0f, ATOM_C);   // 7: CH3

    // Methoxy at C3 (position 4)
    addAtom(mol, r * cosf(PI) - 1.2f, 0.0f, 0.0f, ATOM_O);                          // 8: O
    addAtom(mol, r * cosf(PI) - 2.4f, 0.0f, 0.0f, ATOM_C);                          // 9: CH3

    // Methoxy at C4 (position 5)
    addAtom(mol, r * cosf(4*PI/3) - 1.0f, r * sinf(4*PI/3) - 0.5f, 0.0f, ATOM_O);   // 10: O
    addAtom(mol, r * cosf(4*PI/3) - 2.2f, r * sinf(4*PI/3) - 0.3f, 0.0f, ATOM_C);   // 11: CH3

    // Ethylamine side chain at C0: -CH2-CH2-NH2
    addAtom(mol, r + 1.3f, 0.0f, 0.0f, ATOM_C);         // 12: CH2
    addAtom(mol, r + 2.6f, 0.5f, 0.0f, ATOM_C);         // 13: CH2
    addAtom(mol, r + 3.8f, 0.0f, 0.0f, ATOM_N);         // 14: NH2

    // Hydrogens on benzene (positions 2 and 6)
    float rH = 2.4f;
    addAtom(mol, rH * cosf(PI/3), rH * sinf(PI/3), 0.0f, ATOM_H);           // 15: H on C1
    addAtom(mol, rH * cosf(5*PI/3), rH * sinf(5*PI/3), 0.0f, ATOM_H);       // 16: H on C5

    // Hydrogens on amine
    addAtom(mol, r + 4.3f, 0.5f, 0.7f, ATOM_H);         // 17: H on NH2
    addAtom(mol, r + 4.3f, 0.3f, -0.8f, ATOM_H);        // 18: H on NH2

    // Methyl hydrogens (simplified - one per CH3)
    addAtom(mol, r * cosf(2*PI/3) - 2.6f, r * sinf(2*PI/3) + 1.2f, 0.0f, ATOM_H);   // 19
    addAtom(mol, r * cosf(PI) - 2.8f, 0.9f, 0.0f, ATOM_H);                          // 20
    addAtom(mol, r * cosf(4*PI/3) - 2.6f, r * sinf(4*PI/3) - 1.2f, 0.0f, ATOM_H);   // 21

    // Benzene ring bonds
    for (int i = 0; i < 6; i++) {
        addBond(mol, i, (i + 1) % 6, (i % 2 == 0) ? 2 : 1);
    }

    // Methoxy bonds
    addBond(mol, 2, 6, 1);   // C-O
    addBond(mol, 6, 7, 1);   // O-CH3
    addBond(mol, 3, 8, 1);   // C-O
    addBond(mol, 8, 9, 1);   // O-CH3
    addBond(mol, 4, 10, 1);  // C-O
    addBond(mol, 10, 11, 1); // O-CH3

    // Ethylamine bonds
    addBond(mol, 0, 12, 1);  // C(ring)-CH2
    addBond(mol, 12, 13, 1); // CH2-CH2
    addBond(mol, 13, 14, 1); // CH2-NH2

    // Hydrogen bonds
    addBond(mol, 1, 15, 1);
    addBond(mol, 5, 16, 1);
    addBond(mol, 14, 17, 1);
    addBond(mol, 14, 18, 1);
    addBond(mol, 7, 19, 1);
    addBond(mol, 9, 20, 1);
    addBond(mol, 11, 21, 1);

    centerMolecule(mol);
}

// Random molecule generator
float randf() { return (float)rand() / RAND_MAX; }

void buildRandomMolecule(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Random Molecule");

    // Start with a carbon backbone
    int backboneLength = 4 + rand() % 6;

    // Add backbone carbons
    for (int i = 0; i < backboneLength; i++) {
        float angle = randf() * 0.5f - 0.25f;
        float x = i * 1.5f + randf() * 0.3f;
        float y = sinf(i * 0.8f) * 0.8f + randf() * 0.3f;
        float z = cosf(i * 0.5f) * 0.5f + randf() * 0.3f;
        addAtom(mol, x, y, z, ATOM_C);
    }

    // Connect backbone
    for (int i = 0; i < backboneLength - 1; i++) {
        int order = (rand() % 4 == 0) ? 2 : 1;
        addBond(mol, i, i + 1, order);
    }

    // Maybe add a ring
    if (backboneLength >= 5 && rand() % 2 == 0) {
        addBond(mol, 0, backboneLength - 1, 1);
    }

    // Add functional groups
    for (int i = 0; i < backboneLength; i++) {
        int numH = 2 + rand() % 2;

        // Random chance of heteroatom
        if (rand() % 4 == 0) {
            int hetero = (rand() % 3 == 0) ? ATOM_N : ATOM_O;
            float angle = randf() * TWO_PI;
            float x = mol->atoms[i].x + cosf(angle) * 1.3f;
            float y = mol->atoms[i].y + sinf(angle) * 1.3f;
            float z = mol->atoms[i].z + (randf() - 0.5f) * 0.8f;
            int heteroIdx = mol->numAtoms;
            addAtom(mol, x, y, z, hetero);
            addBond(mol, i, heteroIdx, (rand() % 2 == 0) ? 2 : 1);

            // Add H to O or N
            if (rand() % 2 == 0) {
                float hx = x + (randf() - 0.5f) * 1.0f;
                float hy = y + 0.8f;
                float hz = z + (randf() - 0.5f) * 0.5f;
                int hIdx = mol->numAtoms;
                addAtom(mol, hx, hy, hz, ATOM_H);
                addBond(mol, heteroIdx, hIdx, 1);
            }
            numH--;
        }

        // Add hydrogens
        for (int h = 0; h < numH && mol->numAtoms < MAX_ATOMS - 1; h++) {
            float angle = h * PI + randf() * 0.5f;
            float hx = mol->atoms[i].x + cosf(angle) * 0.9f;
            float hy = mol->atoms[i].y + (randf() - 0.5f) * 0.8f;
            float hz = mol->atoms[i].z + sinf(angle) * 0.9f;
            int hIdx = mol->numAtoms;
            addAtom(mol, hx, hy, hz, ATOM_H);
            addBond(mol, i, hIdx, 1);
        }
    }

    centerMolecule(mol);
}

// Build Cubane (C8H8) - cube-shaped hydrocarbon
void buildCubane(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Cubane (C8H8)");

    // Cubane has 8 carbons at the vertices of a cube
    // Cube edge length ~1.55 Angstroms for C-C bond
    float s = 0.9f;  // Half edge length

    // 8 Carbon atoms at cube vertices
    addAtom(mol, -s, -s, -s, ATOM_C);  // 0
    addAtom(mol,  s, -s, -s, ATOM_C);  // 1
    addAtom(mol,  s,  s, -s, ATOM_C);  // 2
    addAtom(mol, -s,  s, -s, ATOM_C);  // 3
    addAtom(mol, -s, -s,  s, ATOM_C);  // 4
    addAtom(mol,  s, -s,  s, ATOM_C);  // 5
    addAtom(mol,  s,  s,  s, ATOM_C);  // 6
    addAtom(mol, -s,  s,  s, ATOM_C);  // 7

    // 8 Hydrogen atoms pointing outward from each carbon
    float h = 1.5f;
    addAtom(mol, -h, -h, -h, ATOM_H);  // 8
    addAtom(mol,  h, -h, -h, ATOM_H);  // 9
    addAtom(mol,  h,  h, -h, ATOM_H);  // 10
    addAtom(mol, -h,  h, -h, ATOM_H);  // 11
    addAtom(mol, -h, -h,  h, ATOM_H);  // 12
    addAtom(mol,  h, -h,  h, ATOM_H);  // 13
    addAtom(mol,  h,  h,  h, ATOM_H);  // 14
    addAtom(mol, -h,  h,  h, ATOM_H);  // 15

    // C-C bonds (12 edges of the cube)
    // Bottom face
    addBond(mol, 0, 1, 1);
    addBond(mol, 1, 2, 1);
    addBond(mol, 2, 3, 1);
    addBond(mol, 3, 0, 1);
    // Top face
    addBond(mol, 4, 5, 1);
    addBond(mol, 5, 6, 1);
    addBond(mol, 6, 7, 1);
    addBond(mol, 7, 4, 1);
    // Vertical edges
    addBond(mol, 0, 4, 1);
    addBond(mol, 1, 5, 1);
    addBond(mol, 2, 6, 1);
    addBond(mol, 3, 7, 1);

    // C-H bonds
    for (int i = 0; i < 8; i++) {
        addBond(mol, i, i + 8, 1);
    }

    centerMolecule(mol);
}

// Build Buckminsterfullerene (C60) - spherical carbon cage
void buildBuckminsterfullerene(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Buckminsterfullerene C60");

    // C60 has 60 carbon atoms arranged in a truncated icosahedron
    // (like a soccer ball - 12 pentagons and 20 hexagons)
    // Each carbon bonds to exactly 3 others, giving 90 bonds total

    float r = 2.5f;  // Radius scaling factor for visualization

    // Golden ratio
    float phi = (1.0f + sqrtf(5.0f)) / 2.0f;  // ~1.618

    // Truncated icosahedron vertices (normalized)
    // These come from three types of coordinate permutations:
    // Type A: (0, ±1, ±3φ) and cyclic permutations - 12 vertices
    // Type B: (±2, ±(1+2φ), ±φ) and cyclic permutations - 24 vertices
    // Type C: (±1, ±(2+φ), ±2φ) and cyclic permutations - 24 vertices

    float coords[60][3];
    int idx = 0;

    // Normalization factor (all vertices lie on same sphere)
    float norm = sqrtf(9.0f * phi * phi + 1.0f);

    // Type A: (0, ±1, ±3φ) - 4 vertices, then cyclic permutations
    float a_coords[4][3] = {
        {0, 1, 3*phi}, {0, 1, -3*phi}, {0, -1, 3*phi}, {0, -1, -3*phi}
    };
    for (int i = 0; i < 4; i++) {
        coords[idx][0] = a_coords[i][0] / norm;
        coords[idx][1] = a_coords[i][1] / norm;
        coords[idx][2] = a_coords[i][2] / norm;
        idx++;
    }
    // Cyclic permutation (y, z, x)
    for (int i = 0; i < 4; i++) {
        coords[idx][0] = a_coords[i][1] / norm;
        coords[idx][1] = a_coords[i][2] / norm;
        coords[idx][2] = a_coords[i][0] / norm;
        idx++;
    }
    // Cyclic permutation (z, x, y)
    for (int i = 0; i < 4; i++) {
        coords[idx][0] = a_coords[i][2] / norm;
        coords[idx][1] = a_coords[i][0] / norm;
        coords[idx][2] = a_coords[i][1] / norm;
        idx++;
    }

    // Type B: (±2, ±(1+2φ), ±φ) - 8 vertices, then cyclic permutations
    float b1 = 2.0f, b2 = 1.0f + 2.0f * phi, b3 = phi;
    float b_signs[8][3] = {
        {1,1,1}, {1,1,-1}, {1,-1,1}, {1,-1,-1},
        {-1,1,1}, {-1,1,-1}, {-1,-1,1}, {-1,-1,-1}
    };
    for (int i = 0; i < 8; i++) {
        coords[idx][0] = b_signs[i][0] * b1 / norm;
        coords[idx][1] = b_signs[i][1] * b2 / norm;
        coords[idx][2] = b_signs[i][2] * b3 / norm;
        idx++;
    }
    for (int i = 0; i < 8; i++) {
        coords[idx][0] = b_signs[i][1] * b2 / norm;
        coords[idx][1] = b_signs[i][2] * b3 / norm;
        coords[idx][2] = b_signs[i][0] * b1 / norm;
        idx++;
    }
    for (int i = 0; i < 8; i++) {
        coords[idx][0] = b_signs[i][2] * b3 / norm;
        coords[idx][1] = b_signs[i][0] * b1 / norm;
        coords[idx][2] = b_signs[i][1] * b2 / norm;
        idx++;
    }

    // Type C: (±1, ±(2+φ), ±2φ) - 8 vertices, then cyclic permutations
    float c1 = 1.0f, c2 = 2.0f + phi, c3 = 2.0f * phi;
    for (int i = 0; i < 8; i++) {
        coords[idx][0] = b_signs[i][0] * c1 / norm;
        coords[idx][1] = b_signs[i][1] * c2 / norm;
        coords[idx][2] = b_signs[i][2] * c3 / norm;
        idx++;
    }
    for (int i = 0; i < 8; i++) {
        coords[idx][0] = b_signs[i][1] * c2 / norm;
        coords[idx][1] = b_signs[i][2] * c3 / norm;
        coords[idx][2] = b_signs[i][0] * c1 / norm;
        idx++;
    }
    for (int i = 0; i < 8; i++) {
        coords[idx][0] = b_signs[i][2] * c3 / norm;
        coords[idx][1] = b_signs[i][0] * c1 / norm;
        coords[idx][2] = b_signs[i][1] * c2 / norm;
        idx++;
    }

    // Add all 60 carbon atoms
    for (int i = 0; i < 60; i++) {
        addAtom(mol, coords[i][0] * r, coords[i][1] * r, coords[i][2] * r, ATOM_C);
    }

    // C60 has 90 bonds - each carbon bonds to exactly 3 others
    // Bond length in normalized coords: pentagon edges and hexagon edges
    // Find the bond threshold by looking at nearest neighbor distance
    float bondThreshold = 0.42f;  // Adjusted for normalized coordinates

    int bondCount = 0;
    for (int i = 0; i < 60; i++) {
        for (int j = i + 1; j < 60; j++) {
            float dx = coords[i][0] - coords[j][0];
            float dy = coords[i][1] - coords[j][1];
            float dz = coords[i][2] - coords[j][2];
            float dist = sqrtf(dx*dx + dy*dy + dz*dz);
            if (dist < bondThreshold) {
                addBond(mol, i, j, 1);
                bondCount++;
            }
        }
    }

    centerMolecule(mol);
}

// ============== MOLECULE REGISTRY (233 molecules) ==============

typedef void (*MoleculeBuilder)(Molecule*);

struct MoleculeInfo {
    MoleculeBuilder builder;
    const char* name;
    int category;
    const char* description;
};

static MoleculeInfo molecules[] = {
    // === SIMPLE MOLECULES (0-31 original + gases) ===
    { buildWater, "Water", CAT_SIMPLE, "H2O - Universal solvent" },
    { buildMethane, "Methane", CAT_SIMPLE, "CH4 - Natural gas" },
    { buildBenzene, "Benzene", CAT_ORGANIC, "C6H6 - Aromatic ring" },
    { buildEthanol, "Ethanol", CAT_ORGANIC, "C2H5OH - Alcohol" },
    { buildCaffeine, "Caffeine", CAT_PHARMA, "C8H10N4O2 - Stimulant" },
    { buildAdenine, "Adenine", CAT_ENERGY, "DNA nucleobase" },
    { buildGlucose, "Glucose", CAT_SUGARS, "C6H12O6 - Blood sugar" },
    { buildAspirin, "Aspirin/Bayer", CAT_PHARMA, "Pain reliever" },
    { buildAmmonia, "Ammonia", CAT_SIMPLE, "NH3 - Cleaning agent" },
    { buildCarbonDioxide, "CO2", CAT_SIMPLE, "Carbon dioxide" },
    { buildFormaldehyde, "Formaldehyde", CAT_ORGANIC, "CH2O - Preservative" },
    { buildAcetone, "Acetone", CAT_ORGANIC, "C3H6O - Solvent" },
    { buildAceticAcid, "Acetic Acid", CAT_ACIDS, "CH3COOH - Vinegar" },
    { buildPropane, "Propane", CAT_ORGANIC, "C3H8 - Fuel" },
    { buildButane, "Butane", CAT_ORGANIC, "C4H10 - Lighter fuel" },
    { buildCyclohexane, "Cyclohexane", CAT_ORGANIC, "C6H12 - Cyclic alkane" },
    { buildNaphthalene, "Naphthalene", CAT_ORGANIC, "C10H8 - Mothballs" },
    { buildUrea, "Urea", CAT_ORGANIC, "CH4N2O - Metabolic waste" },
    { buildGlycine, "Glycine", CAT_AMINO_ACIDS, "Simplest amino acid" },
    { buildAlanine, "Alanine", CAT_AMINO_ACIDS, "Nonpolar amino acid" },
    { buildThymine, "Thymine", CAT_ENERGY, "DNA nucleobase (T)" },
    { buildCytosine, "Cytosine", CAT_ENERGY, "DNA/RNA nucleobase (C)" },
    { buildGuanine, "Guanine", CAT_ENERGY, "DNA/RNA nucleobase (G)" },
    { buildDopamine, "Dopamine", CAT_NEUROTRANS, "Reward neurotransmitter" },
    { buildSerotonin, "Serotonin", CAT_NEUROTRANS, "Mood neurotransmitter" },
    { buildNitricOxide, "NO", CAT_SIMPLE, "Nitric oxide" },
    { buildHydrogenPeroxide, "H2O2", CAT_HOUSEHOLD, "Hydrogen peroxide" },
    { buildSulfuricAcid, "H2SO4", CAT_ACIDS, "Sulfuric acid" },
    { buildPhosphoricAcid, "H3PO4", CAT_ACIDS, "Phosphoric acid" },
    { buildToluene, "Toluene", CAT_ORGANIC, "C7H8 - Solvent" },
    { buildPhenol, "Phenol", CAT_ORGANIC, "C6H5OH - Carbolic acid" },
    { buildAcetylene, "Acetylene", CAT_ORGANIC, "C2H2 - Welding gas" },
    // === BATCH 1: Gases and solvents (32-63) ===
    { buildOxygen, "O2", CAT_SIMPLE, "Molecular oxygen" },
    { buildNitrogen, "N2", CAT_SIMPLE, "Molecular nitrogen" },
    { buildHydrogen, "H2", CAT_SIMPLE, "Molecular hydrogen" },
    { buildOzone, "O3", CAT_SIMPLE, "Ozone" },
    { buildCarbonMonoxide, "CO", CAT_SIMPLE, "Carbon monoxide" },
    { buildNitrousOxide, "N2O", CAT_SIMPLE, "Laughing gas" },
    { buildSulfurDioxide, "SO2", CAT_SIMPLE, "Sulfur dioxide" },
    { buildHydrogenChloride, "HCl", CAT_ACIDS, "Hydrochloric acid" },
    { buildNitricAcid, "HNO3", CAT_ACIDS, "Nitric acid" },
    { buildMethanol, "Methanol", CAT_ORGANIC, "Wood alcohol" },
    { buildEthane, "Ethane", CAT_ORGANIC, "C2H6" },
    { buildPropene, "Propene", CAT_ORGANIC, "Propylene" },
    { buildIsopropanol, "Isopropanol", CAT_ORGANIC, "Rubbing alcohol" },
    { buildEthyleneGlycol, "Ethylene Glycol", CAT_ORGANIC, "Antifreeze" },
    { buildGlycerol, "Glycerol", CAT_ORGANIC, "Glycerin" },
    { buildAcetaldehyde, "Acetaldehyde", CAT_ORGANIC, "Ethanal" },
    { buildFormicAcid, "Formic Acid", CAT_ACIDS, "Ant venom" },
    { buildLacticAcid, "Lactic Acid", CAT_ACIDS, "Muscle fatigue" },
    { buildEthylAcetate, "Ethyl Acetate", CAT_ORGANIC, "Nail polish solvent" },
    { buildAcetonitrile, "Acetonitrile", CAT_ORGANIC, "Polar solvent" },
    { buildDMSO, "DMSO", CAT_ORGANIC, "Dimethyl sulfoxide" },
    { buildDichloromethane, "CH2Cl2", CAT_ORGANIC, "Dichloromethane" },
    { buildChlorobenzene, "Chlorobenzene", CAT_ORGANIC, "Aromatic halide" },
    { buildNitrobenzene, "Nitrobenzene", CAT_ORGANIC, "Almond odor" },
    { buildAniline, "Aniline", CAT_ORGANIC, "Dye precursor" },
    { buildStyrene, "Styrene", CAT_PLASTICS, "Polystyrene monomer" },
    { buildBenzoicAcid, "Benzoic Acid", CAT_ACIDS, "Preservative" },
    { buildValine, "Valine", CAT_AMINO_ACIDS, "Branched-chain AA" },
    { buildLeucine, "Leucine", CAT_AMINO_ACIDS, "Essential AA" },
    { buildEthylene, "Ethylene", CAT_ORGANIC, "Plant hormone" },
    { buildHydrogenSulfide, "H2S", CAT_SIMPLE, "Rotten egg smell" },
    { buildChloroform, "CHCl3", CAT_ORGANIC, "Chloroform" },
    // === BATCH 2: More organics (64-101) ===
    { buildTertButanol, "tert-Butanol", CAT_ORGANIC, "Tertiary alcohol" },
    { buildButanol, "1-Butanol", CAT_ORGANIC, "Butyl alcohol" },
    { buildDiethylEther, "Diethyl Ether", CAT_ORGANIC, "Classic anesthetic" },
    { buildMTBE, "MTBE", CAT_ORGANIC, "Fuel additive" },
    { buildTHF, "THF", CAT_ORGANIC, "Tetrahydrofuran" },
    { buildDioxane, "1,4-Dioxane", CAT_ORGANIC, "Cyclic ether" },
    { buildDMF, "DMF", CAT_ORGANIC, "Dimethylformamide" },
    { buildCarbonTetrachloride, "CCl4", CAT_ORGANIC, "Carbon tet" },
    { buildMethylAcetate, "Methyl Acetate", CAT_ORGANIC, "Ester solvent" },
    { buildAceticAnhydride, "Acetic Anhydride", CAT_ORGANIC, "Acetylating agent" },
    { buildPropionicAcid, "Propionic Acid", CAT_ACIDS, "Preservative" },
    { buildButyricAcid, "Butyric Acid", CAT_ACIDS, "Rancid butter" },
    { buildSuccinicAcid, "Succinic Acid", CAT_ACIDS, "Krebs cycle" },
    { buildBenzaldehyde, "Benzaldehyde", CAT_FLAVORS, "Almond flavor" },
    { buildBromobenzene, "Bromobenzene", CAT_ORGANIC, "Aromatic halide" },
    { buildPXylene, "p-Xylene", CAT_ORGANIC, "Para-xylene" },
    { buildAnisole, "Anisole", CAT_ORGANIC, "Methoxybenzene" },
    { buildPhenylacetylene, "Phenylacetylene", CAT_ORGANIC, "Aromatic alkyne" },
    { buildFructose, "Fructose", CAT_SUGARS, "Fruit sugar" },
    { buildRibose, "Ribose", CAT_SUGARS, "RNA sugar" },
    { buildDeoxyribose, "Deoxyribose", CAT_SUGARS, "DNA sugar" },
    { buildIsoleucine, "Isoleucine", CAT_AMINO_ACIDS, "Branched-chain AA" },
    { buildSerine, "Serine", CAT_AMINO_ACIDS, "Polar amino acid" },
    { buildThreonine, "Threonine", CAT_AMINO_ACIDS, "Essential AA" },
    { buildAsparticAcid, "Aspartic Acid", CAT_AMINO_ACIDS, "Acidic AA" },
    { buildGlutamicAcid, "Glutamic Acid", CAT_AMINO_ACIDS, "MSG precursor" },
    { buildLysine, "Lysine", CAT_AMINO_ACIDS, "Basic AA" },
    { buildHistidine, "Histidine", CAT_AMINO_ACIDS, "Aromatic AA" },
    { buildPhenylalanine, "Phenylalanine", CAT_AMINO_ACIDS, "Aromatic AA" },
    { buildTyrosine, "Tyrosine", CAT_AMINO_ACIDS, "Aromatic AA" },
    { buildTryptophan, "Tryptophan", CAT_AMINO_ACIDS, "Aromatic AA" },
    { buildProline, "Proline", CAT_AMINO_ACIDS, "Cyclic AA" },
    { buildCysteine, "Cysteine", CAT_AMINO_ACIDS, "Sulfur AA" },
    { buildMethionine, "Methionine", CAT_AMINO_ACIDS, "Sulfur AA" },
    { buildPyruvate, "Pyruvate", CAT_ENERGY, "Glycolysis product" },
    { buildArginine, "Arginine", CAT_AMINO_ACIDS, "Basic AA" },
    { buildAsparagine, "Asparagine", CAT_AMINO_ACIDS, "Polar AA" },
    { buildGlutamine, "Glutamine", CAT_AMINO_ACIDS, "Polar AA" },
    // === VITAMINS (102-115) ===
    { buildAscorbicAcid, "Vitamin C", CAT_VITAMINS, "Ascorbic acid" },
    { buildThiamine, "Vitamin B1", CAT_VITAMINS, "Thiamine" },
    { buildRiboflavin, "Vitamin B2", CAT_VITAMINS, "Riboflavin" },
    { buildNiacin, "Vitamin B3", CAT_VITAMINS, "Niacin" },
    { buildPanthothenicAcid, "Vitamin B5", CAT_VITAMINS, "Pantothenic acid" },
    { buildPyridoxine, "Vitamin B6", CAT_VITAMINS, "Pyridoxine" },
    { buildBiotin, "Vitamin B7", CAT_VITAMINS, "Biotin" },
    { buildFolicAcid, "Vitamin B9", CAT_VITAMINS, "Folic acid" },
    { buildRetinol, "Vitamin A", CAT_VITAMINS, "Retinol" },
    { buildBetaCarotene, "Beta-Carotene", CAT_VITAMINS, "Provitamin A" },
    { buildCholecalciferol, "Vitamin D3", CAT_VITAMINS, "Cholecalciferol" },
    { buildAlphaTocopherol, "Vitamin E", CAT_VITAMINS, "Alpha-tocopherol" },
    { buildPhylloquinone, "Vitamin K1", CAT_VITAMINS, "Phylloquinone" },
    { buildNicotinamide, "Nicotinamide", CAT_VITAMINS, "B3 amide form" },
    // === CONTROLLED SUBSTANCES (116-122) ===
    { buildCocaine, "Cocaine", CAT_PHARMA, "Stimulant alkaloid" },
    { buildHeroin, "Heroin", CAT_PHARMA, "Opioid" },
    { buildFentanyl, "Fentanyl/Sublimaze", CAT_PHARMA, "Synthetic opioid" },
    { buildPropofol, "Propofol/Diprivan", CAT_PHARMA, "IV anesthetic" },
    { buildTHC, "THC", CAT_PHARMA, "Cannabis active" },
    { buildCreatine, "Creatine", CAT_ENERGY, "Muscle energy" },
    { buildOctane, "Octane", CAT_ORGANIC, "Gasoline component" },
    // === NSAIDS & STATINS (123-130) ===
    { buildSimvastatin, "Simvastatin/Zocor", CAT_PHARMA, "Cholesterol drug" },
    { buildIbuprofen, "Ibuprofen/Advil", CAT_PHARMA, "NSAID" },
    { buildNaproxen, "Naproxen/Aleve", CAT_PHARMA, "NSAID" },
    { buildDiclofenac, "Diclofenac/Voltaren", CAT_PHARMA, "NSAID" },
    { buildIndomethacin, "Indomethacin/Indocin", CAT_PHARMA, "NSAID" },
    { buildCelecoxib, "Celecoxib/Celebrex", CAT_PHARMA, "COX-2 inhibitor" },
    { buildMeloxicam, "Meloxicam/Mobic", CAT_PHARMA, "NSAID" },
    { buildAcetaminophen, "Tylenol", CAT_PHARMA, "Paracetamol" },
    // === STEROID HORMONES (131-140) ===
    { buildTestosterone, "Testosterone", CAT_HORMONES, "Male sex hormone" },
    { buildDHT, "DHT", CAT_HORMONES, "Dihydrotestosterone" },
    { buildAndrostenedione, "Androstenedione", CAT_HORMONES, "Androgen precursor" },
    { buildEstradiol, "Estradiol/E2", CAT_HORMONES, "Primary estrogen" },
    { buildEstrone, "Estrone/E1", CAT_HORMONES, "Postmenopausal estrogen" },
    { buildEstriol, "Estriol/E3", CAT_HORMONES, "Pregnancy estrogen" },
    { buildProgesterone, "Progesterone", CAT_HORMONES, "Pregnancy hormone" },
    { buildCortisol, "Cortisol", CAT_HORMONES, "Stress hormone" },
    { buildCortisone, "Cortisone", CAT_HORMONES, "Cortisol precursor" },
    { buildAldosterone, "Aldosterone", CAT_HORMONES, "Salt balance" },
    // === HOUSEHOLD CHEMICALS (141-150) ===
    { buildBleach, "Bleach/NaOCl", CAT_HOUSEHOLD, "Disinfectant" },
    { buildLye, "Lye/NaOH", CAT_HOUSEHOLD, "Caustic soda" },
    { buildTableSalt, "Table Salt/NaCl", CAT_HOUSEHOLD, "Sodium chloride" },
    { buildBakingSoda, "Baking Soda/NaHCO3", CAT_HOUSEHOLD, "Sodium bicarbonate" },
    { buildWashingSoda, "Washing Soda/Na2CO3", CAT_HOUSEHOLD, "Sodium carbonate" },
    { buildSodiumFluoride, "Sodium Fluoride/NaF", CAT_HOUSEHOLD, "Toothpaste" },
    { buildSodiumNitrate, "Sodium Nitrate/NaNO3", CAT_HOUSEHOLD, "Preservative" },
    { buildSodiumNitrite, "Sodium Nitrite/NaNO2", CAT_HOUSEHOLD, "Meat curing" },
    { buildSodiumSulfate, "Sodium Sulfate/Na2SO4", CAT_HOUSEHOLD, "Detergent filler" },
    { buildMSG, "MSG", CAT_HOUSEHOLD, "Flavor enhancer" },
    // === ACIDS (151-161) ===
    { buildCitricAcid, "Citric Acid", CAT_ACIDS, "Citrus fruits" },
    { buildCarbonicAcid, "Carbonic Acid/H2CO3", CAT_ACIDS, "CO2 in water" },
    { buildBoricAcid, "Boric Acid/H3BO3", CAT_ACIDS, "Antiseptic" },
    { buildOxalicAcid, "Oxalic Acid", CAT_ACIDS, "Rust remover" },
    { buildTartaricAcid, "Tartaric Acid", CAT_ACIDS, "Wine acid" },
    { buildMalicAcid, "Malic Acid", CAT_ACIDS, "Apple acid" },
    { buildHydrofluoricAcid, "Hydrofluoric Acid/HF", CAT_ACIDS, "Glass etching" },
    { buildHydrobromicAcid, "Hydrobromic Acid/HBr", CAT_ACIDS, "Strong acid" },
    { buildHydroiodicAcid, "Hydroiodic Acid/HI", CAT_ACIDS, "Strong acid" },
    { buildPerchloricAcid, "Perchloric Acid/HClO4", CAT_ACIDS, "Strongest acid" },
    { buildAcrylicAcid, "Acrylic Acid", CAT_ACIDS, "Polymer precursor" },
    // === PLASTICS (162-179) ===
    { buildVinylChloride, "Vinyl Chloride/PVC", CAT_PLASTICS, "PVC monomer" },
    { buildMethylMethacrylate, "MMA/Plexiglas", CAT_PLASTICS, "Acrylic monomer" },
    { buildTetrafluoroethylene, "TFE/Teflon", CAT_PLASTICS, "PTFE monomer" },
    { buildCaprolactam, "Caprolactam/Nylon-6", CAT_PLASTICS, "Nylon-6 monomer" },
    { buildAcrylonitrile, "Acrylonitrile/ABS", CAT_PLASTICS, "ABS component" },
    { buildButadiene, "1,3-Butadiene/ABS", CAT_PLASTICS, "Rubber/ABS" },
    { buildBisphenolA, "Bisphenol A/BPA", CAT_PLASTICS, "Polycarbonate" },
    { buildTerephthalicAcid, "Terephthalic Acid/PET", CAT_PLASTICS, "PET monomer" },
    { buildVinylAcetate, "Vinyl Acetate/EVA", CAT_PLASTICS, "EVA monomer" },
    { buildDimethylsiloxane, "PDMS/Silicone", CAT_PLASTICS, "Silicone monomer" },
    { buildAdipicAcid, "Adipic Acid/Nylon-6,6", CAT_PLASTICS, "Nylon monomer" },
    { buildHexamethylenediamine, "HMDA/Nylon-6,6", CAT_PLASTICS, "Nylon monomer" },
    { buildIsoprene, "Isoprene/Rubber", CAT_PLASTICS, "Natural rubber" },
    { buildPEDimer, "PE Dimer/HDPE", CAT_PLASTICS, "Polyethylene" },
    { buildPPDimer, "PP Dimer", CAT_PLASTICS, "Polypropylene" },
    { buildPVCTrimer, "PVC Trimer", CAT_PLASTICS, "PVC oligomer" },
    { buildPSDimer, "PS Dimer", CAT_PLASTICS, "Polystyrene" },
    { buildPTFETrimer, "PTFE Trimer/Teflon", CAT_PLASTICS, "Teflon oligomer" },
    // === FATTY ACIDS (180-189) ===
    { buildElaidicAcid, "Elaidic Acid/trans-9", CAT_FATS, "Industrial trans fat" },
    { buildOleicAcid, "Oleic Acid/cis-9", CAT_FATS, "Olive oil" },
    { buildVaccenicAcid, "Vaccenic Acid/trans-11", CAT_FATS, "Natural trans fat" },
    { buildTransPalmitoleicAcid, "trans-Palmitoleic", CAT_FATS, "Dairy fat" },
    { buildLinoelaidicAcid, "Linoelaidic/trans,trans", CAT_FATS, "Hydrogenation" },
    { buildRumenicAcid, "Rumenic Acid/CLA", CAT_FATS, "Conjugated linoleic" },
    { buildStearicAcid, "Stearic Acid/C18:0", CAT_FATS, "Saturated fat" },
    { buildPalmiticAcid, "Palmitic Acid/C16:0", CAT_FATS, "Palm oil" },
    { buildLinoleicAcid, "Linoleic Acid/Omega-6", CAT_FATS, "Essential FA" },
    { buildAlphaLinolenicAcid, "ALA/Omega-3", CAT_FATS, "Plant omega-3" },
    // === METAL COMPOUNDS (190-198) ===
    { buildRust, "Rust/Fe2O3", CAT_METALS, "Iron oxide" },
    { buildMagnetite, "Magnetite/Fe3O4", CAT_METALS, "Magnetic iron" },
    { buildIronPentacarbonyl, "Iron Pentacarbonyl", CAT_METALS, "Fe(CO)5" },
    { buildCopperSulfate, "Copper Sulfate/CuSO4", CAT_METALS, "Blue vitriol" },
    { buildVerdigris, "Verdigris/Patina", CAT_METALS, "Copper carbonate" },
    { buildFerrocene, "Ferrocene", CAT_METALS, "Sandwich compound" },
    { buildCisplatin, "Cisplatin", CAT_METALS, "Cancer drug" },
    { buildAluminumOxide, "Alumina/Al2O3", CAT_METALS, "Aluminum oxide" },
    { buildTitaniumDioxide, "Titania/TiO2", CAT_METALS, "White pigment" },
    // === NEUROTRANSMITTERS (199-202) ===
    { buildAcetylcholine, "Acetylcholine", CAT_NEUROTRANS, "Muscle/memory" },
    { buildGABA, "GABA", CAT_NEUROTRANS, "Inhibitory NT" },
    { buildGlutamate, "Glutamate", CAT_NEUROTRANS, "Excitatory NT" },
    { buildNorepinephrine, "Norepinephrine", CAT_NEUROTRANS, "Fight-or-flight" },
    // === ENERGY MOLECULES (203-205) ===
    { buildATP, "ATP", CAT_ENERGY, "Energy currency" },
    { buildADP, "ADP", CAT_ENERGY, "Adenosine diphosphate" },
    { buildNADH, "NADH", CAT_ENERGY, "Electron carrier" },
    // === ANTIBIOTICS (206-208) ===
    { buildPenicillinG, "Penicillin G", CAT_PHARMA, "Beta-lactam" },
    { buildAmoxicillin, "Amoxicillin", CAT_PHARMA, "Broad-spectrum" },
    { buildVancomycin, "Vancomycin", CAT_PHARMA, "Glycopeptide" },
    // === EXPLOSIVES (209-211) ===
    { buildTNT, "TNT", CAT_OTHER, "Trinitrotoluene" },
    { buildNitroglycerin, "Nitroglycerin", CAT_PHARMA, "Explosive/angina" },
    { buildRDX, "RDX", CAT_OTHER, "Military explosive" },
    // === SWEETENERS (212-215) ===
    { buildSucrose, "Sucrose/Table Sugar", CAT_SUGARS, "Disaccharide" },
    { buildAspartame, "Aspartame", CAT_OTHER, "NutraSweet" },
    { buildSaccharin, "Saccharin", CAT_OTHER, "Sweet'N Low" },
    { buildSucralose, "Sucralose/Splenda", CAT_OTHER, "Chlorinated sugar" },
    // === PESTICIDES (216-218) ===
    { buildDDT, "DDT", CAT_OTHER, "Pesticide" },
    { buildGlyphosate, "Glyphosate/Roundup", CAT_OTHER, "Herbicide" },
    { buildMalathion, "Malathion", CAT_OTHER, "Insecticide" },
    // === FRAGRANCES (219-222) ===
    { buildVanillin, "Vanillin/Vanilla", CAT_FLAVORS, "Vanilla flavor" },
    { buildLimonene, "Limonene/Citrus", CAT_FLAVORS, "Citrus scent" },
    { buildMenthol, "Menthol/Mint", CAT_FLAVORS, "Cooling sensation" },
    { buildCinnamaldehyde, "Cinnamaldehyde/Cinnamon", CAT_FLAVORS, "Cinnamon" },
    // === LACHRYMATORY (223) ===
    { buildPropanethialSoxide, "Onion Factor/Tears", CAT_FLAVORS, "Makes you cry" },
    // === ENVIRONMENTAL (224-225) ===
    { buildCFC12, "CFC-12/Freon", CAT_OTHER, "Ozone depleter" },
    { buildSF6, "SF6/Greenhouse", CAT_OTHER, "Potent GHG" },
    // === ANESTHETICS (226-228) ===
    { buildLidocaine, "Lidocaine/Xylocaine", CAT_PHARMA, "Local anesthetic" },
    { buildKetamine, "Ketamine", CAT_PHARMA, "Dissociative" },
    { buildSevoflurane, "Sevoflurane", CAT_PHARMA, "Inhalation anesthetic" },
    // === PSYCHEDELICS (229-231) ===
    { buildLSD, "LSD", CAT_PHARMA, "Psychedelic" },
    { buildPsilocybin, "Psilocybin/Shrooms", CAT_PHARMA, "Magic mushroom" },
    { buildMescaline, "Mescaline/Peyote", CAT_PHARMA, "Peyote cactus" },
    // === RANDOM (232) ===
    { buildRandomMolecule, "Random", CAT_OTHER, "Random structure" },
    // === EXOTIC STRUCTURES (233-234) ===
    { buildCubane, "Cubane", CAT_ORGANIC, "C8H8 - Cube-shaped" },
    { buildBuckminsterfullerene, "Buckyball/C60", CAT_ORGANIC, "Carbon soccer ball" },
};

static const int NUM_MOLECULES = sizeof(molecules) / sizeof(molecules[0]);

// ============== PUBLIC API ==============

int molecule_get_count() {
    return NUM_MOLECULES;
}

const char* molecule_get_name(int index) {
    if (index < 0 || index >= NUM_MOLECULES) return "Unknown";
    return molecules[index].name;
}

int molecule_get_category(int index) {
    if (index < 0 || index >= NUM_MOLECULES) return CAT_OTHER;
    return molecules[index].category;
}

const char* molecule_get_description(int index) {
    if (index < 0 || index >= NUM_MOLECULES) return "";
    return molecules[index].description;
}

void molecule_build(int index, Molecule* mol) {
    if (index < 0 || index >= NUM_MOLECULES) index = 0;
    molecules[index].builder(mol);
}

void molecule_build_random(Molecule* mol) {
    buildRandomMolecule(mol);
}
