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
    "Controlled Substances", // CAT_CONTROLLED
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
    float radii[ATOM_TYPE_COUNT] = {
        0.25f, 0.40f, 0.38f, 0.35f, 0.45f, 0.45f, 0.45f, 0.50f, 0.35f, 0.55f,
        0.55f, 0.48f, 0.42f, 0.55f, 0.50f, 0.50f, 0.52f, 0.55f,
        0.62f // Re (approx)
    };
    a->radius = radii[(type >= 0 && type < ATOM_TYPE_COUNT) ? type : 0];
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
    // Missing H on C3 to satisfy C6H8O6
    addAtom(mol, 2.6f, -0.9f, 0.0f, ATOM_H);

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
    addBond(mol, 2, 19, 1);

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

// Build Lactose (C12H22O11) - Milk sugar (glucose + galactose disaccharide)
void buildLactose(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Lactose (C12H22O11)");

    // Glucose ring (6-membered pyranose)
    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_C);         // 0: C1 (anomeric)
    addAtom(mol, 1.2f, 0.7f, 0.3f, ATOM_C);         // 1: C2
    addAtom(mol, 2.4f, 0.0f, 0.0f, ATOM_C);         // 2: C3
    addAtom(mol, 2.4f, -1.4f, 0.5f, ATOM_C);        // 3: C4
    addAtom(mol, 1.2f, -2.1f, 0.2f, ATOM_C);        // 4: C5
    addAtom(mol, 0.0f, -1.4f, 0.5f, ATOM_O);        // 5: ring O

    // Glucose hydroxyl groups and CH2OH
    addAtom(mol, -1.0f, 0.5f, 0.5f, ATOM_O);        // 6: OH on C1 (anomeric)
    addAtom(mol, 1.2f, 2.1f, 0.0f, ATOM_O);         // 7: OH on C2
    addAtom(mol, 3.5f, 0.7f, 0.3f, ATOM_O);         // 8: OH on C3
    addAtom(mol, 3.5f, -2.1f, 0.2f, ATOM_O);        // 9: OH on C4
    addAtom(mol, 1.2f, -3.5f, 0.5f, ATOM_C);        // 10: C6 (CH2OH)
    addAtom(mol, 1.2f, -4.2f, 1.5f, ATOM_O);        // 11: OH on C6

    // Galactose ring (6-membered pyranose, glucose C4 epimer)
    // Connected via glycosidic bond from glucose C1 to galactose C4
    addAtom(mol, -2.2f, 0.5f, -0.3f, ATOM_C);       // 12: C1 (CH2OH)
    addAtom(mol, -3.4f, -0.2f, 0.0f, ATOM_C);       // 13: C2
    addAtom(mol, -4.6f, 0.5f, 0.7f, ATOM_C);        // 14: C3 (galactose - flipped stereochemistry)
    addAtom(mol, -4.6f, 1.9f, -0.6f, ATOM_C);       // 15: C4 (galactose epimer)
    addAtom(mol, -3.4f, 2.6f, -0.1f, ATOM_C);       // 16: C5
    addAtom(mol, -2.2f, 1.9f, 0.2f, ATOM_O);        // 17: ring O

    // Galactose hydroxyl groups and CH2OH
    addAtom(mol, -1.0f, -0.3f, 0.0f, ATOM_O);       // 18: OH on C1
    addAtom(mol, -3.4f, -1.6f, -0.5f, ATOM_O);      // 19: OH on C2
    addAtom(mol, -5.8f, 1.2f, 1.0f, ATOM_O);        // 20: OH on C3 (galactose position)
    addAtom(mol, -5.8f, 1.2f, -0.9f, ATOM_O);       // 21: OH on C4 (galactose position)
    addAtom(mol, -3.4f, 4.0f, 0.2f, ATOM_C);        // 22: C6 (CH2OH)
    addAtom(mol, -3.4f, 4.7f, -0.8f, ATOM_O);       // 23: OH on C6

    // Selected hydrogens for clarity
    addAtom(mol, 0.0f, 0.0f, 1.0f, ATOM_H);         // 24: H on glucose C1
    addAtom(mol, 1.2f, 0.7f, 1.3f, ATOM_H);         // 25: H on glucose C2
    addAtom(mol, -4.6f, 1.9f, 1.2f, ATOM_H);        // 26: H on galactose C4

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

    // Glycosidic bond (glucose C1 O to galactose C4)
    addBond(mol, 6, 15, 1);

    // Galactose ring bonds
    addBond(mol, 12, 13, 1);
    addBond(mol, 13, 14, 1);
    addBond(mol, 14, 15, 1);
    addBond(mol, 15, 16, 1);
    addBond(mol, 16, 17, 1);
    addBond(mol, 17, 12, 1); // ring closure

    // Galactose substituents
    addBond(mol, 12, 18, 1); // C1-OH
    addBond(mol, 13, 19, 1); // C2-OH
    addBond(mol, 14, 20, 1); // C3-OH
    addBond(mol, 15, 21, 1); // C4-OH
    addBond(mol, 16, 22, 1); // C5-C6
    addBond(mol, 22, 23, 1); // C6-OH

    // Selected hydrogen bonds
    addBond(mol, 0, 24, 1);
    addBond(mol, 1, 25, 1);
    addBond(mol, 15, 26, 1);

    centerMolecule(mol);
}

// Build Maltose (C12H22O11) - Malt sugar (two glucose units)
void buildMaltose(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Maltose (C12H22O11)");

    // Glucose ring 1 (6-membered pyranose)
    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_C);         // 0: C1 (anomeric)
    addAtom(mol, 1.2f, 0.7f, 0.3f, ATOM_C);         // 1: C2
    addAtom(mol, 2.4f, 0.0f, 0.0f, ATOM_C);         // 2: C3
    addAtom(mol, 2.4f, -1.4f, 0.5f, ATOM_C);        // 3: C4
    addAtom(mol, 1.2f, -2.1f, 0.2f, ATOM_C);        // 4: C5
    addAtom(mol, 0.0f, -1.4f, 0.5f, ATOM_O);        // 5: ring O

    // Glucose 1 hydroxyl groups and CH2OH
    addAtom(mol, -1.0f, 0.5f, 0.5f, ATOM_O);        // 6: OH on C1 (anomeric)
    addAtom(mol, 1.2f, 2.1f, 0.0f, ATOM_O);         // 7: OH on C2
    addAtom(mol, 3.5f, 0.7f, 0.3f, ATOM_O);         // 8: OH on C3
    addAtom(mol, 3.5f, -2.1f, 0.2f, ATOM_O);        // 9: OH on C4
    addAtom(mol, 1.2f, -3.5f, 0.5f, ATOM_C);        // 10: C6 (CH2OH)
    addAtom(mol, 1.2f, -4.2f, 1.5f, ATOM_O);        // 11: OH on C6

    // Glucose ring 2 (6-membered pyranose)
    // Connected via glycosidic bond from glucose1 C1 to glucose2 C4
    addAtom(mol, -2.2f, 0.5f, -0.3f, ATOM_C);       // 12: C1 (CH2OH)
    addAtom(mol, -3.4f, -0.2f, 0.0f, ATOM_C);       // 13: C2
    addAtom(mol, -4.6f, 0.5f, -0.3f, ATOM_C);       // 14: C3
    addAtom(mol, -4.6f, 1.9f, 0.2f, ATOM_C);        // 15: C4 (glycosidic linkage)
    addAtom(mol, -3.4f, 2.6f, -0.1f, ATOM_C);       // 16: C5
    addAtom(mol, -2.2f, 1.9f, 0.2f, ATOM_O);        // 17: ring O

    // Glucose 2 hydroxyl groups and CH2OH
    addAtom(mol, -1.0f, -0.3f, 0.0f, ATOM_O);       // 18: OH on C1
    addAtom(mol, -3.4f, -1.6f, -0.5f, ATOM_O);      // 19: OH on C2
    addAtom(mol, -5.8f, -0.2f, -0.6f, ATOM_O);      // 20: OH on C3
    addAtom(mol, -5.8f, 2.6f, 0.5f, ATOM_O);        // 21: OH on C4
    addAtom(mol, -3.4f, 4.0f, 0.2f, ATOM_C);        // 22: C6 (CH2OH)
    addAtom(mol, -3.4f, 4.7f, -0.8f, ATOM_O);       // 23: OH on C6

    // Selected hydrogens
    addAtom(mol, 0.0f, 0.0f, 1.0f, ATOM_H);         // 24: H on glucose1 C1
    addAtom(mol, 1.2f, 0.7f, 1.3f, ATOM_H);         // 25: H on glucose1 C2
    addAtom(mol, -4.6f, 1.9f, 1.2f, ATOM_H);        // 26: H on glucose2 C4

    // Glucose 1 ring bonds
    addBond(mol, 0, 1, 1);
    addBond(mol, 1, 2, 1);
    addBond(mol, 2, 3, 1);
    addBond(mol, 3, 4, 1);
    addBond(mol, 4, 5, 1);
    addBond(mol, 5, 0, 1);  // ring closure

    // Glucose 1 substituents
    addBond(mol, 0, 6, 1);  // C1-OH (anomeric)
    addBond(mol, 1, 7, 1);  // C2-OH
    addBond(mol, 2, 8, 1);  // C3-OH
    addBond(mol, 3, 9, 1);  // C4-OH
    addBond(mol, 4, 10, 1); // C5-C6
    addBond(mol, 10, 11, 1);// C6-OH

    // Glycosidic bond (glucose1 C1 O to glucose2 C4)
    addBond(mol, 6, 15, 1);

    // Glucose 2 ring bonds
    addBond(mol, 12, 13, 1);
    addBond(mol, 13, 14, 1);
    addBond(mol, 14, 15, 1);
    addBond(mol, 15, 16, 1);
    addBond(mol, 16, 17, 1);
    addBond(mol, 17, 12, 1); // ring closure

    // Glucose 2 substituents
    addBond(mol, 12, 18, 1); // C1-OH
    addBond(mol, 13, 19, 1); // C2-OH
    addBond(mol, 14, 20, 1); // C3-OH
    addBond(mol, 15, 21, 1); // C4-OH
    addBond(mol, 16, 22, 1); // C5-C6
    addBond(mol, 22, 23, 1); // C6-OH

    // Selected hydrogen bonds
    addBond(mol, 0, 24, 1);
    addBond(mol, 1, 25, 1);
    addBond(mol, 15, 26, 1);

    centerMolecule(mol);
}

// Build Galactose (C6H12O6) - Milk sugar monosaccharide (glucose epimer)
void buildGalactose(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Galactose (C6H12O6)");

    // Galactose ring (6-membered pyranose)
    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_C);         // 0: C1 (anomeric)
    addAtom(mol, 1.2f, 0.7f, 0.3f, ATOM_C);         // 1: C2
    addAtom(mol, 2.4f, 0.0f, 0.0f, ATOM_C);         // 2: C3
    addAtom(mol, 2.4f, -1.4f, 0.5f, ATOM_C);        // 3: C4 (differs from glucose)
    addAtom(mol, 1.2f, -2.1f, 0.2f, ATOM_C);        // 4: C5
    addAtom(mol, 0.0f, -1.4f, 0.5f, ATOM_O);        // 5: ring O

    // Hydroxyl groups and CH2OH
    addAtom(mol, -1.0f, 0.5f, 0.5f, ATOM_O);        // 6: OH on C1 (anomeric)
    addAtom(mol, 1.2f, 2.1f, 0.0f, ATOM_O);         // 7: OH on C2
    addAtom(mol, 3.5f, 0.7f, 0.3f, ATOM_O);         // 8: OH on C3
    addAtom(mol, 3.5f, -2.1f, 0.2f, ATOM_O);        // 9: OH on C4 (axial, differs from glucose)
    addAtom(mol, 1.2f, -3.5f, 0.5f, ATOM_C);        // 10: C6 (CH2OH)
    addAtom(mol, 1.2f, -4.2f, 1.5f, ATOM_O);        // 11: OH on C6

    // Hydrogens
    addAtom(mol, 0.0f, 0.0f, 1.0f, ATOM_H);         // 12: H on C1
    addAtom(mol, 1.2f, 0.7f, 1.3f, ATOM_H);         // 13: H on C2
    addAtom(mol, 2.4f, 0.0f, -1.0f, ATOM_H);        // 14: H on C3
    addAtom(mol, 2.4f, -1.4f, 1.5f, ATOM_H);        // 15: H on C4
    addAtom(mol, 1.2f, -2.1f, -0.8f, ATOM_H);       // 16: H on C5

    // Ring bonds
    addBond(mol, 0, 1, 1);
    addBond(mol, 1, 2, 1);
    addBond(mol, 2, 3, 1);
    addBond(mol, 3, 4, 1);
    addBond(mol, 4, 5, 1);
    addBond(mol, 5, 0, 1);  // ring closure

    // Hydroxyl and CH2OH bonds
    addBond(mol, 0, 6, 1);  // C1-OH (anomeric)
    addBond(mol, 1, 7, 1);  // C2-OH
    addBond(mol, 2, 8, 1);  // C3-OH
    addBond(mol, 3, 9, 1);  // C4-OH
    addBond(mol, 4, 10, 1); // C5-C6
    addBond(mol, 10, 11, 1);// C6-OH

    // Hydrogen bonds
    addBond(mol, 0, 12, 1);
    addBond(mol, 1, 13, 1);
    addBond(mol, 2, 14, 1);
    addBond(mol, 3, 15, 1);
    addBond(mol, 4, 16, 1);

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

// Build PCP - Phencyclidine (C17H25N)
void buildPCP(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "PCP (C17H25N)");

    // 1-phenylcyclohexyl piperidine
    // Cyclohexyl ring
    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_C);         // 0: C1 (with phenyl)
    addAtom(mol, 1.4f, 0.5f, 0.0f, ATOM_C);         // 1: C2
    addAtom(mol, 2.0f, 1.8f, 0.5f, ATOM_C);         // 2: C3
    addAtom(mol, 1.2f, 2.8f, 0.8f, ATOM_C);         // 3: C4
    addAtom(mol, -0.2f, 2.3f, 0.8f, ATOM_C);        // 4: C5
    addAtom(mol, -0.8f, 1.0f, 0.3f, ATOM_C);        // 5: C6

    // Phenyl ring attached to C1
    float r = 1.4f;
    float baseX = -1.2f, baseY = -1.0f;
    for (int i = 0; i < 6; i++) {
        float angle = i * PI / 3.0f;
        addAtom(mol, baseX + r * cosf(angle), baseY + r * sinf(angle), 0.0f, ATOM_C);  // 6-11
    }

    // Piperidine ring attached to C1
    addAtom(mol, 0.2f, -0.5f, 1.4f, ATOM_N);        // 12: N
    addAtom(mol, -0.8f, -1.3f, 1.6f, ATOM_C);       // 13: CH2
    addAtom(mol, -1.6f, -2.2f, 0.8f, ATOM_C);       // 14: CH2
    addAtom(mol, -0.8f, -3.0f, -0.2f, ATOM_C);      // 15: CH2
    addAtom(mol, 0.6f, -2.8f, -0.4f, ATOM_C);       // 16: CH2
    addAtom(mol, 1.4f, -1.9f, 0.4f, ATOM_C);        // 17: CH2

    // Hydrogens on cyclohexyl ring (C1-C5)
    addAtom(mol, 2.0f, -0.2f, -0.5f, ATOM_H);       // 18: H on C1
    addAtom(mol, 1.6f, 0.8f, 0.7f, ATOM_H);         // 19: H on C1
    addAtom(mol, 2.9f, 2.1f, 0.0f, ATOM_H);         // 20: H on C2
    addAtom(mol, 2.3f, 2.2f, 1.4f, ATOM_H);         // 21: H on C2
    addAtom(mol, 2.0f, 3.3f, -0.1f, ATOM_H);        // 22: H on C3
    addAtom(mol, 2.7f, 3.5f, 1.0f, ATOM_H);         // 23: H on C3
    addAtom(mol, 1.5f, 3.8f, 1.2f, ATOM_H);         // 24: H on C4
    addAtom(mol, 1.3f, 2.5f, 1.6f, ATOM_H);         // 25: H on C4
    addAtom(mol, -0.9f, 3.0f, 1.4f, ATOM_H);        // 26: H on C5
    addAtom(mol, -0.4f, 2.8f, 0.0f, ATOM_H);        // 27: H on C5

    // Hydrogens on phenyl ring
    addAtom(mol, -2.0f, -1.3f, 0.0f, ATOM_H);       // 28: H on C7
    addAtom(mol, -2.2f, -0.2f, -0.6f, ATOM_H);      // 29: H on C8
    addAtom(mol, -1.8f, 0.8f, -0.8f, ATOM_H);       // 30: H on C9
    addAtom(mol, -0.8f, 1.3f, -0.5f, ATOM_H);       // 31: H on C10
    addAtom(mol, 0.2f, 0.7f, 0.2f, ATOM_H);         // 32: H on C11

    // Hydrogens on piperidine ring (C13-C17)
    addAtom(mol, -1.3f, -0.8f, 2.3f, ATOM_H);       // 33: H on C13
    addAtom(mol, -0.5f, -1.8f, 2.0f, ATOM_H);       // 34: H on C13
    addAtom(mol, -2.5f, -1.8f, 0.5f, ATOM_H);       // 35: H on C14
    addAtom(mol, -1.9f, -2.9f, 1.3f, ATOM_H);       // 36: H on C14
    addAtom(mol, -1.2f, -3.9f, -0.1f, ATOM_H);      // 37: H on C15
    addAtom(mol, -0.8f, -2.8f, -1.1f, ATOM_H);      // 38: H on C15
    addAtom(mol, 1.0f, -3.6f, -0.8f, ATOM_H);       // 39: H on C16
    addAtom(mol, 0.8f, -2.3f, -1.2f, ATOM_H);       // 40: H on C16
    addAtom(mol, 2.3f, -1.8f, 0.2f, ATOM_H);        // 41: H on C17
    addAtom(mol, 1.6f, -1.0f, 1.1f, ATOM_H);        // 42: H on C17

    // Cyclohexyl ring bonds
    addBond(mol, 0, 1, 1);
    addBond(mol, 1, 2, 1);
    addBond(mol, 2, 3, 1);
    addBond(mol, 3, 4, 1);
    addBond(mol, 4, 5, 1);
    addBond(mol, 5, 0, 1);

    // Phenyl ring bonds
    for (int i = 0; i < 6; i++) {
        addBond(mol, 6 + i, 6 + ((i + 1) % 6), (i % 2 == 0) ? 2 : 1);
    }

    // Phenyl attachment
    addBond(mol, 0, 6, 1);

    // Piperidine ring bonds
    addBond(mol, 0, 12, 1);  // C1-N
    addBond(mol, 12, 13, 1); // N-CH2
    addBond(mol, 13, 14, 1); // CH2-CH2
    addBond(mol, 14, 15, 1); // CH2-CH2
    addBond(mol, 15, 16, 1); // CH2-CH2
    addBond(mol, 16, 17, 1); // CH2-CH2
    addBond(mol, 17, 12, 1); // CH2-N (closing ring)

    // Hydrogen bonds on cyclohexyl (C1-C5)
    addBond(mol, 1, 18, 1);  // C1-H
    addBond(mol, 1, 19, 1);  // C1-H
    addBond(mol, 2, 20, 1);  // C2-H
    addBond(mol, 2, 21, 1);  // C2-H
    addBond(mol, 3, 22, 1);  // C3-H
    addBond(mol, 3, 23, 1);  // C3-H
    addBond(mol, 4, 24, 1);  // C4-H
    addBond(mol, 4, 25, 1);  // C4-H
    addBond(mol, 5, 26, 1);  // C5-H
    addBond(mol, 5, 27, 1);  // C5-H

    // Hydrogen bonds on phenyl ring
    addBond(mol, 7, 28, 1);  // C7-H
    addBond(mol, 8, 29, 1);  // C8-H
    addBond(mol, 9, 30, 1);  // C9-H
    addBond(mol, 10, 31, 1); // C10-H
    addBond(mol, 11, 32, 1); // C11-H

    // Hydrogen bonds on piperidine (C13-C17)
    addBond(mol, 13, 33, 1); // C13-H
    addBond(mol, 13, 34, 1); // C13-H
    addBond(mol, 14, 35, 1); // C14-H
    addBond(mol, 14, 36, 1); // C14-H
    addBond(mol, 15, 37, 1); // C15-H
    addBond(mol, 15, 38, 1); // C15-H
    addBond(mol, 16, 39, 1); // C16-H
    addBond(mol, 16, 40, 1); // C16-H
    addBond(mol, 17, 41, 1); // C17-H
    addBond(mol, 17, 42, 1); // C17-H

    centerMolecule(mol);
}

// Build MDMA - 3,4-Methylenedioxymethamphetamine (C11H15NO2)
void buildMDMA(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "MDMA/Ecstasy (C11H15NO2)");

    // Methylenedioxy benzene ring
    float r = 1.4f;
    for (int i = 0; i < 6; i++) {
        float angle = i * PI / 3.0f;
        addAtom(mol, r * cosf(angle), r * sinf(angle), 0.0f, ATOM_C);  // 0-5
    }

    // Methylenedioxy bridge (OCH2O) at positions 2,3 (C3,C4)
    addAtom(mol, r * cosf(2*PI/3) - 0.8f, r * sinf(2*PI/3) + 0.8f, 0.5f, ATOM_O);     // 6: O
    addAtom(mol, r * cosf(2*PI/3) - 0.8f, r * sinf(2*PI/3) + 0.8f, -0.5f, ATOM_C);    // 7: CH2
    addAtom(mol, r * cosf(PI) - 0.8f, r * sinf(PI) + 0.8f, -0.5f, ATOM_O);            // 8: O

    // Methamphetamine side chain: -CH2-CH(CH3)-NH(CH3) at C0
    addAtom(mol, r + 1.3f, 0.0f, 0.0f, ATOM_C);         // 9: CH2
    addAtom(mol, r + 2.6f, 0.5f, 0.0f, ATOM_C);         // 10: CH (chiral center)
    addAtom(mol, r + 2.8f, 1.9f, 0.0f, ATOM_C);         // 11: CH3 (methyl on chiral)
    addAtom(mol, r + 3.8f, -0.3f, 0.0f, ATOM_N);        // 12: NH
    addAtom(mol, r + 5.1f, 0.3f, 0.0f, ATOM_C);         // 13: CH3 (N-methyl)

    // Hydrogens on benzene
    float rH = 2.4f;
    addAtom(mol, rH * cosf(PI/3), rH * sinf(PI/3), 0.0f, ATOM_H);           // 14: H on C1
    addAtom(mol, rH * cosf(5*PI/3), rH * sinf(5*PI/3), 0.0f, ATOM_H);       // 15: H on C5

    // Hydrogen on amine
    addAtom(mol, r + 4.0f, -1.2f, 0.0f, ATOM_H);        // 16: H on NH

    // Benzene ring bonds
    for (int i = 0; i < 6; i++) {
        addBond(mol, i, (i + 1) % 6, (i % 2 == 0) ? 2 : 1);
    }

    // Methylenedioxy bonds
    addBond(mol, 2, 6, 1);   // C-O
    addBond(mol, 6, 7, 1);   // O-CH2
    addBond(mol, 7, 8, 1);   // CH2-O
    addBond(mol, 8, 3, 1);   // O-C

    // Side chain bonds
    addBond(mol, 0, 9, 1);   // C(ring)-CH2
    addBond(mol, 9, 10, 1);  // CH2-CH
    addBond(mol, 10, 11, 1); // CH-CH3
    addBond(mol, 10, 12, 1); // CH-NH
    addBond(mol, 12, 13, 1); // NH-CH3

    // Hydrogen bonds
    addBond(mol, 1, 14, 1);
    addBond(mol, 5, 15, 1);
    addBond(mol, 12, 16, 1);

    centerMolecule(mol);
}

// Build Methylone - Bath Salts (C11H13NO2)
void buildMethylone(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Methylone/Bath Salts (C11H13NO2)");

    // Methylenedioxy benzene ring
    float r = 1.4f;
    for (int i = 0; i < 6; i++) {
        float angle = i * PI / 3.0f;
        addAtom(mol, r * cosf(angle), r * sinf(angle), 0.0f, ATOM_C);  // 0-5
    }

    // Methylenedioxy bridge (OCH2O) at positions 2,3 (C3,C4)
    addAtom(mol, r * cosf(2*PI/3) - 0.8f, r * sinf(2*PI/3) + 0.8f, 0.5f, ATOM_O);     // 6: O
    addAtom(mol, r * cosf(2*PI/3) - 0.8f, r * sinf(2*PI/3) + 0.8f, -0.5f, ATOM_C);    // 7: CH2
    addAtom(mol, r * cosf(PI) - 0.8f, r * sinf(PI) + 0.8f, -0.5f, ATOM_O);            // 8: O

    // Cathinone side chain: -CH2-CO-CH(CH3)-NH2 at C0
    addAtom(mol, r + 1.3f, 0.0f, 0.0f, ATOM_C);         // 9: CH2
    addAtom(mol, r + 2.6f, 0.5f, 0.0f, ATOM_C);         // 10: C=O (carbonyl)
    addAtom(mol, r + 2.6f, 0.5f, 1.3f, ATOM_O);         // 11: O (ketone)
    addAtom(mol, r + 3.8f, -0.3f, -0.5f, ATOM_C);       // 12: CH
    addAtom(mol, r + 4.0f, 0.6f, -1.8f, ATOM_C);        // 13: CH3
    addAtom(mol, r + 5.2f, -0.9f, -0.2f, ATOM_N);       // 14: NH2

    // Hydrogens on benzene
    float rH = 2.4f;
    addAtom(mol, rH * cosf(PI/3), rH * sinf(PI/3), 0.0f, ATOM_H);           // 15: H on C1
    addAtom(mol, rH * cosf(5*PI/3), rH * sinf(5*PI/3), 0.0f, ATOM_H);       // 16: H on C5

    // Hydrogens on amine
    addAtom(mol, r + 5.7f, -1.6f, -0.7f, ATOM_H);       // 17: H on NH2
    addAtom(mol, r + 5.7f, -0.8f, 0.5f, ATOM_H);        // 18: H on NH2

    // Hydrogen on chiral center
    addAtom(mol, r + 3.8f, -1.2f, -0.5f, ATOM_H);       // 19: H on CH

    // Benzene ring bonds
    for (int i = 0; i < 6; i++) {
        addBond(mol, i, (i + 1) % 6, (i % 2 == 0) ? 2 : 1);
    }

    // Methylenedioxy bonds
    addBond(mol, 2, 6, 1);   // C-O
    addBond(mol, 6, 7, 1);   // O-CH2
    addBond(mol, 7, 8, 1);   // CH2-O
    addBond(mol, 8, 3, 1);   // O-C

    // Side chain bonds
    addBond(mol, 0, 9, 1);   // C(ring)-CH2
    addBond(mol, 9, 10, 1);  // CH2-C(=O)
    addBond(mol, 10, 11, 2); // C=O (ketone)
    addBond(mol, 10, 12, 1); // C-CH
    addBond(mol, 12, 13, 1); // CH-CH3
    addBond(mol, 12, 14, 1); // CH-NH2

    // Hydrogen bonds
    addBond(mol, 1, 15, 1);
    addBond(mol, 5, 16, 1);
    addBond(mol, 14, 17, 1);
    addBond(mol, 14, 18, 1);
    addBond(mol, 12, 19, 1);

    centerMolecule(mol);
}

// Build Valium - Diazepam (C16H13ClN2O)
void buildValium(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Valium/Diazepam (C16H13ClN2O)");

    // Fused ring system: benzene-7-membered ring-benzene
    // First benzene ring (top)
    float r = 1.4f;
    for (int i = 0; i < 6; i++) {
        float angle = i * PI / 3.0f + PI/2;
        addAtom(mol, r * cosf(angle), r * sinf(angle) + 1.2f, 0.0f, ATOM_C);  // 0-5
    }

    // Seven-membered middle ring
    addAtom(mol, 0.0f, 1.2f, 0.0f, ATOM_C);           // 6: shared with ring 1
    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_C);           // 7
    addAtom(mol, -1.4f, -0.7f, 0.0f, ATOM_C);         // 8
    addAtom(mol, -1.4f, -2.0f, 0.0f, ATOM_N);         // 9: N
    addAtom(mol, 0.0f, -2.7f, 0.0f, ATOM_C);          // 10
    addAtom(mol, 1.4f, -2.0f, 0.0f, ATOM_N);          // 11: N
    addAtom(mol, 1.4f, -0.7f, 0.0f, ATOM_C);          // 12: shared with ring 2

    // Second benzene ring (bottom)
    for (int i = 0; i < 6; i++) {
        float angle = i * PI / 3.0f + PI/2;
        addAtom(mol, r * cosf(angle), r * sinf(angle) - 1.2f, 0.0f, ATOM_C);  // 13-18
    }

    // Chlorine at position on first ring
    addAtom(mol, r * cosf(PI/3 + PI/2) + 0.8f, r * sinf(PI/3 + PI/2) + 1.2f, 0.0f, ATOM_CL);  // 19

    // Methyl group on N
    addAtom(mol, -2.6f, -2.8f, 0.0f, ATOM_C);         // 20: CH3

    // Carbonyl oxygen on middle ring
    addAtom(mol, 0.0f, -4.0f, 0.0f, ATOM_O);          // 21: C=O

    // Hydrogens (simplified)
    float rH = 2.4f;
    addAtom(mol, rH * cosf(PI/6 + PI/2), rH * sinf(PI/6 + PI/2) + 1.2f, 0.0f, ATOM_H);  // 22
    addAtom(mol, rH * cosf(5*PI/6 + PI/2), rH * sinf(5*PI/6 + PI/2) + 1.2f, 0.0f, ATOM_H);  // 23

    // First ring bonds
    for (int i = 0; i < 6; i++) {
        addBond(mol, i, (i + 1) % 6, (i % 2 == 0) ? 2 : 1);
    }

    // Second ring bonds
    for (int i = 0; i < 6; i++) {
        addBond(mol, 13 + i, 13 + ((i + 1) % 6), (i % 2 == 0) ? 2 : 1);
    }

    // Middle ring bonds
    addBond(mol, 6, 7, 1);
    addBond(mol, 7, 8, 1);
    addBond(mol, 8, 9, 1);
    addBond(mol, 9, 10, 1);
    addBond(mol, 10, 11, 1);
    addBond(mol, 11, 12, 1);
    addBond(mol, 12, 6, 1);

    // Ring connections
    addBond(mol, 0, 6, 1);   // First ring to middle
    addBond(mol, 12, 13, 1); // Middle to second ring

    // Chlorine bond
    addBond(mol, 1, 19, 1);

    // Methyl bond
    addBond(mol, 9, 20, 1);

    // Carbonyl bond
    addBond(mol, 10, 21, 2);

    // Hydrogen bonds
    addBond(mol, 5, 22, 1);
    addBond(mol, 4, 23, 1);

    centerMolecule(mol);
}

// Build Methamphetamine - Crystal Meth (C10H15N)
void buildMethamphetamine(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Methamphetamine/Crystal Meth (C10H15N)");

    // Benzene ring
    float r = 1.4f;
    for (int i = 0; i < 6; i++) {
        float angle = i * PI / 3.0f;
        addAtom(mol, r * cosf(angle), r * sinf(angle), 0.0f, ATOM_C);  // 0-5
    }

    // Methamphetamine side chain: -CH2-CH(CH3)-NH(CH3)
    addAtom(mol, r + 1.3f, 0.0f, 0.0f, ATOM_C);       // 6: CH2
    addAtom(mol, r + 2.6f, 0.5f, 0.0f, ATOM_C);       // 7: CH (chiral center)
    addAtom(mol, r + 2.8f, 1.9f, 0.0f, ATOM_C);       // 8: CH3 (on chiral)
    addAtom(mol, r + 3.8f, -0.3f, 0.0f, ATOM_N);      // 9: NH
    addAtom(mol, r + 5.1f, 0.3f, 0.0f, ATOM_C);       // 10: CH3 (N-methyl)

    // Hydrogens on benzene ring (all 6 carbons)
    float rH = 2.4f;
    for (int i = 0; i < 6; i++) {
        float angle = i * PI / 3.0f;
        addAtom(mol, rH * cosf(angle), rH * sinf(angle), 0.0f, ATOM_H);  // 11-16: H on benzene
    }

    // Hydrogens on CH2 group
    addAtom(mol, r + 1.5f, -0.8f, 0.0f, ATOM_H);      // 17: H on CH2
    addAtom(mol, r + 1.1f, 0.8f, 0.0f, ATOM_H);       // 18: H on CH2

    // Hydrogen on CH chiral center
    addAtom(mol, r + 2.6f, 0.5f, -1.0f, ATOM_H);      // 19: H on CH

    // Hydrogens on CH3 (on chiral center)
    addAtom(mol, r + 3.4f, 2.4f, 0.0f, ATOM_H);       // 20: H on CH3
    addAtom(mol, r + 3.0f, 1.9f, -0.8f, ATOM_H);      // 21: H on CH3
    addAtom(mol, r + 2.2f, 2.3f, 0.5f, ATOM_H);       // 22: H on CH3

    // Hydrogen on amine
    addAtom(mol, r + 4.0f, -1.2f, 0.0f, ATOM_H);      // 23: H on NH

    // Hydrogens on N-CH3 methyl group
    addAtom(mol, r + 5.5f, 1.2f, 0.0f, ATOM_H);       // 24: H on N-CH3
    addAtom(mol, r + 5.3f, -0.4f, 0.0f, ATOM_H);      // 25: H on N-CH3
    addAtom(mol, r + 5.8f, 0.4f, -0.8f, ATOM_H);      // 26: H on N-CH3

    // Benzene ring bonds
    for (int i = 0; i < 6; i++) {
        addBond(mol, i, (i + 1) % 6, (i % 2 == 0) ? 2 : 1);
    }

    // Side chain bonds
    addBond(mol, 0, 6, 1);   // C(ring)-CH2
    addBond(mol, 6, 7, 1);   // CH2-CH
    addBond(mol, 7, 8, 1);   // CH-CH3
    addBond(mol, 7, 9, 1);   // CH-NH
    addBond(mol, 9, 10, 1);  // NH-CH3

    // Hydrogen bonds on benzene ring
    for (int i = 0; i < 6; i++) {
        addBond(mol, i, 11 + i, 1);  // H on each benzene carbon
    }

    // Hydrogen bonds on side chain
    addBond(mol, 6, 17, 1);   // H on CH2
    addBond(mol, 6, 18, 1);   // H on CH2
    addBond(mol, 7, 19, 1);   // H on CH
    addBond(mol, 8, 20, 1);   // H on CH3
    addBond(mol, 8, 21, 1);   // H on CH3
    addBond(mol, 8, 22, 1);   // H on CH3
    addBond(mol, 9, 23, 1);   // H on NH
    addBond(mol, 10, 24, 1);  // H on N-CH3
    addBond(mol, 10, 25, 1);  // H on N-CH3
    addBond(mol, 10, 26, 1);  // H on N-CH3

    centerMolecule(mol);
}

// Build Codeine (C18H21NO3)
void buildCodeine(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Codeine (C18H21NO3)");

    // Phenanthrene-like core (fused ring system)
    // Top benzene-like ring
    float r = 1.4f;
    for (int i = 0; i < 6; i++) {
        float angle = i * PI / 3.0f + PI/2;
        addAtom(mol, r * cosf(angle) - 1.0f, r * sinf(angle) + 1.5f, 0.0f, ATOM_C);  // 0-5
    }

    // Middle fused ring
    addAtom(mol, -0.5f, 1.5f, 0.0f, ATOM_C);          // 6
    addAtom(mol, -0.5f, 0.2f, 0.0f, ATOM_C);          // 7
    addAtom(mol, -1.9f, -0.5f, 0.0f, ATOM_C);         // 8
    addAtom(mol, -1.9f, -1.8f, 0.0f, ATOM_N);         // 9: N
    addAtom(mol, -0.5f, -2.5f, 0.0f, ATOM_C);         // 10
    addAtom(mol, 0.9f, -1.8f, 0.0f, ATOM_C);          // 11

    // Bottom benzene ring
    for (int i = 0; i < 6; i++) {
        float angle = i * PI / 3.0f + PI/2;
        addAtom(mol, r * cosf(angle) + 1.0f, r * sinf(angle) - 1.5f, 0.0f, ATOM_C);  // 12-17
    }

    // Methoxy group at position 3
    addAtom(mol, -2.4f, 2.4f, 0.0f, ATOM_O);          // 18: O (methoxy)
    addAtom(mol, -3.6f, 2.4f, 0.0f, ATOM_C);          // 19: CH3

    // Hydroxyl group at position 6
    addAtom(mol, 2.4f, -2.4f, 0.0f, ATOM_O);          // 20: OH
    addAtom(mol, 3.2f, -3.2f, 0.0f, ATOM_H);          // 21: H

    // N-methyl group
    addAtom(mol, -2.8f, -2.6f, 0.0f, ATOM_C);         // 22: CH3 on N

    // Top ring bonds
    for (int i = 0; i < 6; i++) {
        addBond(mol, i, (i + 1) % 6, (i % 2 == 0) ? 2 : 1);
    }

    // Bottom ring bonds
    for (int i = 0; i < 6; i++) {
        addBond(mol, 12 + i, 12 + ((i + 1) % 6), (i % 2 == 0) ? 2 : 1);
    }

    // Middle ring bonds
    addBond(mol, 6, 7, 1);
    addBond(mol, 7, 8, 1);
    addBond(mol, 8, 9, 1);
    addBond(mol, 9, 10, 1);
    addBond(mol, 10, 11, 1);
    addBond(mol, 11, 6, 1);

    // Ring connections
    addBond(mol, 0, 6, 1);
    addBond(mol, 11, 12, 1);

    // Methoxy bonds
    addBond(mol, 5, 18, 1);
    addBond(mol, 18, 19, 1);

    // Hydroxyl bond
    addBond(mol, 17, 20, 1);
    addBond(mol, 20, 21, 1);

    // N-methyl bond
    addBond(mol, 9, 22, 1);

    centerMolecule(mol);
}

// Build Methadone (C21H27NO)
void buildMethadone(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Methadone (C21H27NO)");

    // Two phenyl rings connected by a chain
    // First phenyl ring (left)
    float r = 1.4f;
    float baseX1 = -5.0f;
    for (int i = 0; i < 6; i++) {
        float angle = i * PI / 3.0f;
        addAtom(mol, baseX1 + r * cosf(angle), r * sinf(angle), 0.0f, ATOM_C);  // 0-5
    }

    // Second phenyl ring (right)
    float baseX2 = 5.0f;
    for (int i = 0; i < 6; i++) {
        float angle = i * PI / 3.0f;
        addAtom(mol, baseX2 + r * cosf(angle), r * sinf(angle), 0.0f, ATOM_C);  // 6-11
    }

    // Connecting chain: (C2H5)C(CH3)(C2H5)-CO-CH2-CH2-N(CH3)2
    addAtom(mol, -2.5f, 0.0f, 0.0f, ATOM_C);          // 12: Quaternary C
    addAtom(mol, -1.2f, -1.2f, 0.0f, ATOM_C);         // 13: CH2
    addAtom(mol, 0.2f, -1.2f, 0.0f, ATOM_C);          // 14: C=O (carbonyl)
    addAtom(mol, 0.2f, -1.2f, 1.3f, ATOM_O);          // 15: O (ketone)
    addAtom(mol, 1.5f, -1.2f, -0.5f, ATOM_C);         // 16: CH2
    addAtom(mol, 2.8f, -1.2f, -0.5f, ATOM_C);         // 17: CH2
    addAtom(mol, 4.1f, -1.2f, -0.5f, ATOM_N);         // 18: N

    // Methyl groups and ethyl substituents on quaternary C
    addAtom(mol, -2.5f, 0.0f, 1.4f, ATOM_C);          // 19: CH3
    addAtom(mol, -2.7f, -1.4f, -0.5f, ATOM_C);        // 20: CH2 (ethyl 1)
    addAtom(mol, -3.5f, -2.2f, 0.2f, ATOM_C);         // 21: CH3 (ethyl 1)
    addAtom(mol, -2.7f, 1.4f, -0.5f, ATOM_C);         // 22: CH2 (ethyl 2)
    addAtom(mol, -3.5f, 2.2f, 0.2f, ATOM_C);          // 23: CH3 (ethyl 2)

    // N-methyl groups
    addAtom(mol, 4.3f, -2.6f, -0.5f, ATOM_C);         // 24: CH3 (N-methyl 1)
    addAtom(mol, 4.3f, 0.2f, -0.5f, ATOM_C);          // 25: CH3 (N-methyl 2)

    // Hydrogens on phenyl rings
    float rH = 2.4f;
    addAtom(mol, baseX1 + rH * cosf(PI/3), rH * sinf(PI/3), 0.0f, ATOM_H);     // 26
    addAtom(mol, baseX2 + rH * cosf(PI/3), rH * sinf(PI/3), 0.0f, ATOM_H);     // 27

    // First phenyl ring bonds
    for (int i = 0; i < 6; i++) {
        addBond(mol, i, (i + 1) % 6, (i % 2 == 0) ? 2 : 1);
    }

    // Second phenyl ring bonds
    for (int i = 0; i < 6; i++) {
        addBond(mol, 6 + i, 6 + ((i + 1) % 6), (i % 2 == 0) ? 2 : 1);
    }

    // Chain bonds
    addBond(mol, 0, 12, 1);   // Phenyl1-quaternary C
    addBond(mol, 12, 13, 1);  // Quaternary C-CH2
    addBond(mol, 13, 14, 1);  // CH2-C=O
    addBond(mol, 14, 15, 2);  // C=O
    addBond(mol, 14, 16, 1);  // C-CH2
    addBond(mol, 16, 17, 1);  // CH2-CH2
    addBond(mol, 17, 18, 1);  // CH2-N
    addBond(mol, 6, 18, 1);   // Phenyl2-N (connection through chain)

    // Substituent bonds
    addBond(mol, 12, 19, 1);  // Quaternary C-CH3
    addBond(mol, 12, 20, 1);  // Quaternary C-CH2
    addBond(mol, 20, 21, 1);  // CH2-CH3
    addBond(mol, 12, 22, 1);  // Quaternary C-CH2
    addBond(mol, 22, 23, 1);  // CH2-CH3

    // N-methyl bonds
    addBond(mol, 18, 24, 1);  // N-CH3
    addBond(mol, 18, 25, 1);  // N-CH3

    // Hydrogen bonds
    addBond(mol, 5, 26, 1);
    addBond(mol, 11, 27, 1);

    centerMolecule(mol);
}

// Build Diphenhydramine/Benadryl (C17H21NO) - 1st generation antihistamine
void buildDiphenhydramine(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Diphenhydramine/Benadryl (C17H21NO)");

    // Two phenyl rings connected by ethoxy chain
    // First phenyl ring (left)
    float r = 1.4f;
    float baseX1 = -4.0f;
    for (int i = 0; i < 6; i++) {
        float angle = i * PI / 3.0f;
        addAtom(mol, baseX1 + r * cosf(angle), r * sinf(angle), 0.0f, ATOM_C);  // 0-5
    }

    // Second phenyl ring (right)
    float baseX2 = 4.0f;
    for (int i = 0; i < 6; i++) {
        float angle = i * PI / 3.0f;
        addAtom(mol, baseX2 + r * cosf(angle), r * sinf(angle), 0.0f, ATOM_C);  // 6-11
    }

    // Connecting chain: Ph-CH(OCH2CH2N(CH3)2)-Ph
    addAtom(mol, 0.0f, 1.0f, 0.0f, ATOM_C);            // 12: CH (chiral)
    addAtom(mol, 0.0f, 2.3f, 0.0f, ATOM_O);            // 13: O (ether)
    addAtom(mol, -1.3f, 3.0f, 0.0f, ATOM_C);           // 14: CH2
    addAtom(mol, -1.3f, 4.3f, 0.0f, ATOM_C);           // 15: CH2
    addAtom(mol, 0.0f, 5.0f, 0.0f, ATOM_N);            // 16: N
    addAtom(mol, -0.5f, 6.2f, 0.0f, ATOM_C);           // 17: CH3 (N-methyl)
    addAtom(mol, 1.3f, 5.2f, 0.0f, ATOM_C);            // 18: CH3 (N-methyl)

    // Hydrogens on phenyl rings
    float rH = 2.4f;
    addAtom(mol, baseX1 + rH * cosf(PI/3), rH * sinf(PI/3), 0.0f, ATOM_H);     // 19
    addAtom(mol, baseX2 + rH * cosf(PI/3), rH * sinf(PI/3), 0.0f, ATOM_H);     // 20

    // First phenyl ring bonds
    for (int i = 0; i < 6; i++) {
        addBond(mol, i, (i + 1) % 6, (i % 2 == 0) ? 2 : 1);
    }

    // Second phenyl ring bonds
    for (int i = 0; i < 6; i++) {
        addBond(mol, 6 + i, 6 + ((i + 1) % 6), (i % 2 == 0) ? 2 : 1);
    }

    // Chain bonds
    addBond(mol, 0, 12, 1);   // Phenyl1-CH
    addBond(mol, 12, 6, 1);   // CH-Phenyl2
    addBond(mol, 12, 13, 1);  // CH-O
    addBond(mol, 13, 14, 1);  // O-CH2
    addBond(mol, 14, 15, 1);  // CH2-CH2
    addBond(mol, 15, 16, 1);  // CH2-N
    addBond(mol, 16, 17, 1);  // N-CH3
    addBond(mol, 16, 18, 1);  // N-CH3

    // Hydrogen bonds
    addBond(mol, 5, 19, 1);
    addBond(mol, 11, 20, 1);

    centerMolecule(mol);
}

// Build Chlorpheniramine/Chlor-Trimeton (C16H19ClN2) - 1st generation antihistamine
void buildChlorpheniramine(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Chlorpheniramine/Chlor-Trimeton (C16H19ClN2)");

    // Dibenzothiazepine-like core with piperidine
    // Benzene ring 1
    float r = 1.4f;
    for (int i = 0; i < 6; i++) {
        float angle = i * PI / 3.0f + PI/2;
        addAtom(mol, r * cosf(angle) - 1.0f, r * sinf(angle) + 1.5f, 0.0f, ATOM_C);  // 0-5
    }

    // Middle fused ring with N
    addAtom(mol, -0.5f, 1.5f, 0.0f, ATOM_C);           // 6
    addAtom(mol, -0.5f, 0.2f, 0.0f, ATOM_C);           // 7
    addAtom(mol, -1.9f, -0.5f, 0.0f, ATOM_C);          // 8
    addAtom(mol, -1.9f, -1.8f, 0.0f, ATOM_N);          // 9: N (ring N)
    addAtom(mol, -0.5f, -2.5f, 0.0f, ATOM_C);          // 10
    addAtom(mol, 0.9f, -1.8f, 0.0f, ATOM_C);           // 11

    // Benzene ring 2
    for (int i = 0; i < 6; i++) {
        float angle = i * PI / 3.0f + PI/2;
        addAtom(mol, r * cosf(angle) + 1.0f, r * sinf(angle) - 1.5f, 0.0f, ATOM_C);  // 12-17
    }

    // Chlorine substituent
    addAtom(mol, -2.4f, 2.4f, 0.0f, ATOM_CL);          // 18: Cl

    // Piperidine ring attached to N
    addAtom(mol, -3.2f, -2.6f, 0.0f, ATOM_C);          // 19: CH2
    addAtom(mol, -4.5f, -2.0f, 0.0f, ATOM_C);          // 20: CH2
    addAtom(mol, -4.7f, -0.6f, 0.0f, ATOM_C);          // 21: CH2
    addAtom(mol, -3.5f, 0.2f, 0.0f, ATOM_N);           // 22: N (piperidine)
    addAtom(mol, -2.2f, -0.4f, 0.0f, ATOM_C);          // 23: CH2
    addAtom(mol, -2.0f, -1.8f, 0.0f, ATOM_C);          // 24: CH2

    // Methyl on piperidine N
    addAtom(mol, -3.5f, 1.6f, 0.0f, ATOM_C);           // 25: CH3

    // Hydrogens on aromatic rings
    float rH = 2.4f;
    addAtom(mol, rH * cosf(PI/6 + PI/2) - 1.0f, rH * sinf(PI/6 + PI/2) + 1.5f, 0.0f, ATOM_H);  // 26
    addAtom(mol, rH * cosf(5*PI/6 + PI/2) + 1.0f, rH * sinf(5*PI/6 + PI/2) - 1.5f, 0.0f, ATOM_H);  // 27

    // Benzene ring 1 bonds
    for (int i = 0; i < 6; i++) {
        addBond(mol, i, (i + 1) % 6, (i % 2 == 0) ? 2 : 1);
    }

    // Benzene ring 2 bonds
    for (int i = 0; i < 6; i++) {
        addBond(mol, 12 + i, 12 + ((i + 1) % 6), (i % 2 == 0) ? 2 : 1);
    }

    // Middle ring bonds
    addBond(mol, 6, 7, 1);
    addBond(mol, 7, 8, 1);
    addBond(mol, 8, 9, 1);
    addBond(mol, 9, 10, 1);
    addBond(mol, 10, 11, 1);
    addBond(mol, 11, 6, 1);

    // Ring connections
    addBond(mol, 0, 6, 1);
    addBond(mol, 11, 12, 1);

    // Chlorine bond
    addBond(mol, 5, 18, 1);

    // Piperidine ring bonds
    addBond(mol, 9, 19, 1);
    addBond(mol, 19, 20, 1);
    addBond(mol, 20, 21, 1);
    addBond(mol, 21, 22, 1);
    addBond(mol, 22, 23, 1);
    addBond(mol, 23, 24, 1);
    addBond(mol, 24, 9, 1);

    // Methyl on piperidine
    addBond(mol, 22, 25, 1);

    // Hydrogen bonds
    addBond(mol, 4, 26, 1);
    addBond(mol, 16, 27, 1);

    centerMolecule(mol);
}

// Build Loratadine/Claritin (C22H23ClN2O) - 2nd generation antihistamine
void buildLoratadine(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Loratadine/Claritin (C22H23ClN2O)");

    // Tricyclic core: benzimidazole fused to quinoline
    // Main aromatic system
    float r = 1.4f;

    // Benzene ring
    for (int i = 0; i < 6; i++) {
        float angle = i * PI / 3.0f;
        addAtom(mol, r * cosf(angle), r * sinf(angle) + 1.0f, 0.0f, ATOM_C);  // 0-5
    }

    // Imidazole ring
    addAtom(mol, 1.2f, 0.0f, 0.0f, ATOM_C);            // 6
    addAtom(mol, 0.6f, -1.0f, 0.0f, ATOM_N);           // 7: N
    addAtom(mol, -0.6f, -1.0f, 0.0f, ATOM_C);          // 8
    addAtom(mol, -1.2f, 0.0f, 0.0f, ATOM_N);           // 9: N

    // Chlorophenyl substituent
    float baseX = -3.0f;
    for (int i = 0; i < 6; i++) {
        float angle = i * PI / 3.0f;
        addAtom(mol, baseX + r * cosf(angle), r * sinf(angle) - 2.0f, 0.0f, ATOM_C);  // 10-15
    }

    // Chlorine on substituent
    addAtom(mol, baseX - 1.0f, r * sinf(-PI/3) - 2.0f - 0.8f, 0.0f, ATOM_CL);  // 16

    // Piperidine ring
    addAtom(mol, 2.8f, 0.5f, 0.0f, ATOM_C);            // 17: CH2
    addAtom(mol, 4.0f, 1.3f, 0.0f, ATOM_C);            // 18: CH2
    addAtom(mol, 3.5f, 2.7f, 0.0f, ATOM_N);            // 19: N
    addAtom(mol, 2.3f, 2.5f, 0.0f, ATOM_C);            // 20: CH2
    addAtom(mol, 1.1f, 3.3f, 0.0f, ATOM_C);            // 21: CH2

    // Carbonyl and side chain
    addAtom(mol, 4.7f, 3.5f, 0.0f, ATOM_C);            // 22: C=O
    addAtom(mol, 4.7f, 3.5f, 1.3f, ATOM_O);            // 23: O (carbonyl)

    // Hydrogens
    float rH = 2.4f;
    addAtom(mol, rH * cosf(PI/3), rH * sinf(PI/3) + 1.0f, 0.0f, ATOM_H);  // 24

    // Benzene ring bonds
    for (int i = 0; i < 6; i++) {
        addBond(mol, i, (i + 1) % 6, (i % 2 == 0) ? 2 : 1);
    }

    // Chlorophenyl ring bonds
    for (int i = 0; i < 6; i++) {
        addBond(mol, 10 + i, 10 + ((i + 1) % 6), (i % 2 == 0) ? 2 : 1);
    }

    // Imidazole bonds
    addBond(mol, 1, 6, 1);
    addBond(mol, 6, 7, 1);
    addBond(mol, 7, 8, 1);
    addBond(mol, 8, 9, 1);
    addBond(mol, 9, 4, 1);

    // Connection to chlorophenyl
    addBond(mol, 8, 10, 1);

    // Chlorine bond
    addBond(mol, 12, 16, 1);

    // Piperidine bonds
    addBond(mol, 6, 17, 1);
    addBond(mol, 17, 18, 1);
    addBond(mol, 18, 19, 1);
    addBond(mol, 19, 20, 1);
    addBond(mol, 20, 21, 1);
    addBond(mol, 21, 6, 1);

    // Carbonyl bonds
    addBond(mol, 19, 22, 1);
    addBond(mol, 22, 23, 2);

    // Hydrogen bonds
    addBond(mol, 5, 24, 1);

    centerMolecule(mol);
}

// Build Cetirizine/Zyrtec (C21H25ClN2O) - 2nd generation antihistamine
void buildCetirizine(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Cetirizine/Zyrtec (C21H25ClN2O)");

    // Similar to loratadine but with piperazine side chain and carboxylic acid
    // Benzimidazole core
    float r = 1.4f;

    // Benzene ring
    for (int i = 0; i < 6; i++) {
        float angle = i * PI / 3.0f;
        addAtom(mol, r * cosf(angle), r * sinf(angle) + 1.0f, 0.0f, ATOM_C);  // 0-5
    }

    // Imidazole ring
    addAtom(mol, 1.2f, 0.0f, 0.0f, ATOM_C);            // 6
    addAtom(mol, 0.6f, -1.0f, 0.0f, ATOM_N);           // 7: N
    addAtom(mol, -0.6f, -1.0f, 0.0f, ATOM_C);          // 8
    addAtom(mol, -1.2f, 0.0f, 0.0f, ATOM_N);           // 9: N

    // Chlorophenyl substituent
    float baseX = -3.0f;
    for (int i = 0; i < 6; i++) {
        float angle = i * PI / 3.0f;
        addAtom(mol, baseX + r * cosf(angle), r * sinf(angle) - 2.0f, 0.0f, ATOM_C);  // 10-15
    }

    // Chlorine on substituent
    addAtom(mol, baseX - 1.0f, r * sinf(-PI/3) - 2.0f - 0.8f, 0.0f, ATOM_CL);  // 16

    // Piperazine ring
    addAtom(mol, 2.8f, 0.5f, 0.0f, ATOM_C);            // 17: CH2
    addAtom(mol, 4.0f, 1.3f, 0.0f, ATOM_C);            // 18: CH2
    addAtom(mol, 4.0f, 2.7f, 0.0f, ATOM_N);            // 19: N
    addAtom(mol, 2.8f, 3.5f, 0.0f, ATOM_C);            // 20: CH2
    addAtom(mol, 1.6f, 2.7f, 0.0f, ATOM_C);            // 21: CH2

    // Carboxylic acid side chain
    addAtom(mol, 4.0f, 4.1f, 0.0f, ATOM_C);            // 22: CH2
    addAtom(mol, 5.3f, 4.8f, 0.0f, ATOM_C);            // 23: COOH
    addAtom(mol, 5.3f, 4.8f, 1.3f, ATOM_O);            // 24: O (carbonyl)
    addAtom(mol, 6.5f, 5.6f, 0.0f, ATOM_O);            // 25: OH
    addAtom(mol, 7.3f, 5.6f, 0.0f, ATOM_H);            // 26: H (on OH)

    // Hydrogens
    float rH = 2.4f;
    addAtom(mol, rH * cosf(PI/3), rH * sinf(PI/3) + 1.0f, 0.0f, ATOM_H);  // 27

    // Benzene ring bonds
    for (int i = 0; i < 6; i++) {
        addBond(mol, i, (i + 1) % 6, (i % 2 == 0) ? 2 : 1);
    }

    // Chlorophenyl ring bonds
    for (int i = 0; i < 6; i++) {
        addBond(mol, 10 + i, 10 + ((i + 1) % 6), (i % 2 == 0) ? 2 : 1);
    }

    // Imidazole bonds
    addBond(mol, 1, 6, 1);
    addBond(mol, 6, 7, 1);
    addBond(mol, 7, 8, 1);
    addBond(mol, 8, 9, 1);
    addBond(mol, 9, 4, 1);

    // Connection to chlorophenyl
    addBond(mol, 8, 10, 1);

    // Chlorine bond
    addBond(mol, 12, 16, 1);

    // Piperazine bonds
    addBond(mol, 6, 17, 1);
    addBond(mol, 17, 18, 1);
    addBond(mol, 18, 19, 1);
    addBond(mol, 19, 20, 1);
    addBond(mol, 20, 21, 1);
    addBond(mol, 21, 6, 1);

    // Side chain bonds
    addBond(mol, 19, 22, 1);
    addBond(mol, 22, 23, 1);
    addBond(mol, 23, 24, 2);
    addBond(mol, 23, 25, 1);
    addBond(mol, 25, 26, 1);

    // Hydrogen bonds
    addBond(mol, 5, 27, 1);

    centerMolecule(mol);
}

// Build Fexofenadine/Allegra (C32H39ClN2O4) - 2nd generation antihistamine
void buildFexofenadine(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Fexofenadine/Allegra (C32H39ClN2O4)");

    // Larger non-sedating antihistamine with two aromatic rings
    // First phenyl ring
    float r = 1.4f;
    float baseX1 = -5.0f;
    for (int i = 0; i < 6; i++) {
        float angle = i * PI / 3.0f;
        addAtom(mol, baseX1 + r * cosf(angle), r * sinf(angle), 0.0f, ATOM_C);  // 0-5
    }

    // Second phenyl ring
    float baseX2 = 5.0f;
    for (int i = 0; i < 6; i++) {
        float angle = i * PI / 3.0f;
        addAtom(mol, baseX2 + r * cosf(angle), r * sinf(angle), 0.0f, ATOM_C);  // 6-11
    }

    // Central piperidine ring
    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_C);            // 12
    addAtom(mol, -1.3f, 0.7f, 0.0f, ATOM_C);           // 13
    addAtom(mol, -1.3f, 2.1f, 0.0f, ATOM_N);           // 14: N
    addAtom(mol, 0.0f, 2.8f, 0.0f, ATOM_C);            // 15
    addAtom(mol, 1.3f, 2.1f, 0.0f, ATOM_C);            // 16
    addAtom(mol, 1.3f, 0.7f, 0.0f, ATOM_C);            // 17

    // Chlorine on first phenyl
    addAtom(mol, baseX1 - 1.0f, r * sinf(-PI/3) - 0.8f, 0.0f, ATOM_CL);  // 18

    // Carboxylic acid side chain on second phenyl
    addAtom(mol, baseX2 + 2.0f, r * sinf(PI/3) + 1.0f, 0.0f, ATOM_C);     // 19: CH2
    addAtom(mol, baseX2 + 3.3f, r * sinf(PI/3) + 1.7f, 0.0f, ATOM_C);     // 20: COOH
    addAtom(mol, baseX2 + 3.3f, r * sinf(PI/3) + 1.7f, 1.3f, ATOM_O);     // 21: O
    addAtom(mol, baseX2 + 4.5f, r * sinf(PI/3) + 2.5f, 0.0f, ATOM_O);     // 22: OH
    addAtom(mol, baseX2 + 5.3f, r * sinf(PI/3) + 2.5f, 0.0f, ATOM_H);     // 23: H

    // t-Butyl group on piperidine
    addAtom(mol, -2.6f, 2.8f, 0.0f, ATOM_C);           // 24: C(quaternary)
    addAtom(mol, -3.8f, 2.0f, 0.0f, ATOM_C);           // 25: CH3
    addAtom(mol, -2.8f, 4.2f, 0.0f, ATOM_C);           // 26: CH3
    addAtom(mol, -3.2f, 3.0f, 1.4f, ATOM_C);           // 27: CH3

    // Hydrogens on aromatic rings
    float rH = 2.4f;
    addAtom(mol, baseX1 + rH * cosf(PI/3), rH * sinf(PI/3), 0.0f, ATOM_H);     // 28
    addAtom(mol, baseX2 + rH * cosf(PI/3), rH * sinf(PI/3), 0.0f, ATOM_H);     // 29

    // First phenyl ring bonds
    for (int i = 0; i < 6; i++) {
        addBond(mol, i, (i + 1) % 6, (i % 2 == 0) ? 2 : 1);
    }

    // Second phenyl ring bonds
    for (int i = 0; i < 6; i++) {
        addBond(mol, 6 + i, 6 + ((i + 1) % 6), (i % 2 == 0) ? 2 : 1);
    }

    // Piperidine ring bonds
    addBond(mol, 12, 13, 1);
    addBond(mol, 13, 14, 1);
    addBond(mol, 14, 15, 1);
    addBond(mol, 15, 16, 1);
    addBond(mol, 16, 17, 1);
    addBond(mol, 17, 12, 1);

    // Connections to phenyl rings
    addBond(mol, 0, 12, 1);
    addBond(mol, 6, 17, 1);

    // Chlorine bond
    addBond(mol, 2, 18, 1);

    // Carboxylic acid bonds
    addBond(mol, 8, 19, 1);
    addBond(mol, 19, 20, 1);
    addBond(mol, 20, 21, 2);
    addBond(mol, 20, 22, 1);
    addBond(mol, 22, 23, 1);

    // t-Butyl bonds
    addBond(mol, 14, 24, 1);
    addBond(mol, 24, 25, 1);
    addBond(mol, 24, 26, 1);
    addBond(mol, 24, 27, 1);

    // Hydrogen bonds
    addBond(mol, 5, 28, 1);
    addBond(mol, 11, 29, 1);

    centerMolecule(mol);
}

// Build Pseudoephedrine (C9H13NO) - Decongestant
void buildPseudoephedrine(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Pseudoephedrine (C9H13NO)");

    // Benzene ring
    float r = 1.4f;
    for (int i = 0; i < 6; i++) {
        float angle = i * PI / 3.0f;
        addAtom(mol, r * cosf(angle), r * sinf(angle), 0.0f, ATOM_C);  // 0-5
    }

    // Side chain: -CH(OH)-CH(NHCH3)-CH3
    addAtom(mol, r + 1.3f, 0.0f, 0.0f, ATOM_C);        // 6: CH(OH)
    addAtom(mol, r + 1.5f, -0.5f, 1.3f, ATOM_O);       // 7: OH
    addAtom(mol, r + 2.6f, 0.7f, 0.0f, ATOM_C);        // 8: CH(NHCH3)
    addAtom(mol, r + 2.8f, 2.1f, 0.0f, ATOM_C);        // 9: CH3
    addAtom(mol, r + 3.8f, 0.0f, 0.0f, ATOM_N);        // 10: NH
    addAtom(mol, r + 5.1f, 0.5f, 0.0f, ATOM_C);        // 11: CH3 (N-methyl)

    // Hydrogens on benzene
    float rH = 2.4f;
    addAtom(mol, rH * cosf(PI/3), rH * sinf(PI/3), 0.0f, ATOM_H);     // 12
    addAtom(mol, rH * cosf(5*PI/3), rH * sinf(5*PI/3), 0.0f, ATOM_H); // 13

    // Hydrogens on OH and NH
    addAtom(mol, r + 2.0f, -0.9f, 1.8f, ATOM_H);       // 14: H on OH
    addAtom(mol, r + 4.0f, -1.0f, 0.0f, ATOM_H);       // 15: H on NH

    // Benzene ring bonds
    for (int i = 0; i < 6; i++) {
        addBond(mol, i, (i + 1) % 6, (i % 2 == 0) ? 2 : 1);
    }

    // Side chain bonds
    addBond(mol, 0, 6, 1);   // Benzene-CH(OH)
    addBond(mol, 6, 7, 1);   // CH-OH
    addBond(mol, 6, 8, 1);   // CH-CH
    addBond(mol, 8, 9, 1);   // CH-CH3
    addBond(mol, 8, 10, 1);  // CH-NH
    addBond(mol, 10, 11, 1); // NH-CH3

    // Hydrogen bonds
    addBond(mol, 1, 12, 1);
    addBond(mol, 5, 13, 1);
    addBond(mol, 7, 14, 1);
    addBond(mol, 10, 15, 1);

    centerMolecule(mol);
}

// Build Lisinopril (C21H31N3O5) - ACE Inhibitor
void buildLisinopril(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Lisinopril (C21H31N3O5)");

    // Biphenyl-like core with carboxylic acid and amine groups
    float r = 1.4f;
    for (int i = 0; i < 6; i++) {
        float angle = i * PI / 3.0f;
        addAtom(mol, r * cosf(angle) - 2.0f, r * sinf(angle), 0.0f, ATOM_C);  // 0-5
    }

    for (int i = 0; i < 6; i++) {
        float angle = i * PI / 3.0f;
        addAtom(mol, r * cosf(angle) + 2.0f, r * sinf(angle), 0.0f, ATOM_C);  // 6-11
    }

    // Carboxylic acid group
    addAtom(mol, -3.2f, 1.4f, 0.0f, ATOM_C);       // 12: C=O
    addAtom(mol, -3.2f, 1.4f, 1.3f, ATOM_O);       // 13: O
    addAtom(mol, -4.4f, 2.2f, 0.0f, ATOM_O);       // 14: OH
    addAtom(mol, -5.2f, 2.2f, 0.0f, ATOM_H);       // 15: H

    // Amine side chain
    addAtom(mol, 3.2f, 1.4f, 0.0f, ATOM_N);        // 16: NH2
    addAtom(mol, 4.5f, 0.8f, 0.0f, ATOM_C);        // 17: CH3
    addAtom(mol, 4.5f, -0.6f, 0.0f, ATOM_H);       // 18: H on N

    // Ring bonds
    for (int i = 0; i < 6; i++) {
        addBond(mol, i, (i + 1) % 6, (i % 2 == 0) ? 2 : 1);
        addBond(mol, 6 + i, 6 + ((i + 1) % 6), (i % 2 == 0) ? 2 : 1);
    }

    // Ring connection
    addBond(mol, 0, 6, 1);

    // Carboxylic acid bonds
    addBond(mol, 5, 12, 1);
    addBond(mol, 12, 13, 2);
    addBond(mol, 12, 14, 1);
    addBond(mol, 14, 15, 1);

    // Amine bonds
    addBond(mol, 11, 16, 1);
    addBond(mol, 16, 17, 1);
    addBond(mol, 16, 18, 1);

    centerMolecule(mol);
}

// Build Enalapril (C20H32N2O5) - ACE Inhibitor
void buildEnalapril(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Enalapril (C20H32N2O5)");

    // Phenyl ring core
    float r = 1.4f;
    for (int i = 0; i < 6; i++) {
        float angle = i * PI / 3.0f;
        addAtom(mol, r * cosf(angle), r * sinf(angle), 0.0f, ATOM_C);  // 0-5
    }

    // Aliphatic chain with ester and amine
    addAtom(mol, r + 1.3f, 0.0f, 0.0f, ATOM_C);    // 6: CH2
    addAtom(mol, r + 2.6f, 0.7f, 0.0f, ATOM_C);    // 7: C=O (ester)
    addAtom(mol, r + 2.6f, 0.7f, 1.3f, ATOM_O);    // 8: O
    addAtom(mol, r + 3.8f, -0.1f, 0.0f, ATOM_O);   // 9: O (ester oxygen)
    addAtom(mol, r + 5.1f, 0.5f, 0.0f, ATOM_C);    // 10: CH3 (ester methyl)

    // Amine side chain
    addAtom(mol, r + 2.8f, -1.5f, 0.0f, ATOM_N);   // 11: N
    addAtom(mol, r + 4.1f, -2.2f, 0.0f, ATOM_C);   // 12: CH3
    addAtom(mol, r + 2.8f, -2.4f, 0.0f, ATOM_H);   // 13: H on N

    // Hydrogens on ring
    float rH = 2.4f;
    addAtom(mol, rH * cosf(PI/3), rH * sinf(PI/3), 0.0f, ATOM_H);     // 14

    // Ring bonds
    for (int i = 0; i < 6; i++) {
        addBond(mol, i, (i + 1) % 6, (i % 2 == 0) ? 2 : 1);
    }

    // Chain bonds
    addBond(mol, 0, 6, 1);
    addBond(mol, 6, 7, 1);
    addBond(mol, 7, 8, 2);
    addBond(mol, 7, 9, 1);
    addBond(mol, 9, 10, 1);
    addBond(mol, 7, 11, 1);
    addBond(mol, 11, 12, 1);
    addBond(mol, 11, 13, 1);

    // Ring hydrogen
    addBond(mol, 5, 14, 1);

    centerMolecule(mol);
}

// Build Ramipril (C23H32N2O5) - ACE Inhibitor
void buildRamipril(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Ramipril (C23H32N2O5)");

    // Similar to Enalapril but with additional substituent
    float r = 1.4f;
    for (int i = 0; i < 6; i++) {
        float angle = i * PI / 3.0f;
        addAtom(mol, r * cosf(angle), r * sinf(angle), 0.0f, ATOM_C);  // 0-5
    }

    // Aliphatic chain
    addAtom(mol, r + 1.3f, 0.0f, 0.0f, ATOM_C);    // 6: CH2
    addAtom(mol, r + 2.6f, 0.7f, 0.0f, ATOM_C);    // 7: C=O
    addAtom(mol, r + 2.6f, 0.7f, 1.3f, ATOM_O);    // 8: O
    addAtom(mol, r + 3.9f, 0.0f, 0.0f, ATOM_C);    // 9: CH2
    addAtom(mol, r + 5.2f, 0.7f, 0.0f, ATOM_C);    // 10: C=O (carboxylic)
    addAtom(mol, r + 5.2f, 0.7f, 1.3f, ATOM_O);    // 11: O
    addAtom(mol, r + 6.4f, 1.5f, 0.0f, ATOM_O);    // 12: OH
    addAtom(mol, r + 7.2f, 1.5f, 0.0f, ATOM_H);    // 13: H

    // Amine
    addAtom(mol, r + 2.8f, -1.5f, 0.0f, ATOM_N);   // 14: N
    addAtom(mol, r + 4.1f, -2.2f, 0.0f, ATOM_C);   // 15: CH3

    // Ring bonds
    for (int i = 0; i < 6; i++) {
        addBond(mol, i, (i + 1) % 6, (i % 2 == 0) ? 2 : 1);
    }

    // Chain bonds
    addBond(mol, 0, 6, 1);
    addBond(mol, 6, 7, 1);
    addBond(mol, 7, 8, 2);
    addBond(mol, 7, 9, 1);
    addBond(mol, 9, 10, 1);
    addBond(mol, 10, 11, 2);
    addBond(mol, 10, 12, 1);
    addBond(mol, 12, 13, 1);
    addBond(mol, 7, 14, 1);
    addBond(mol, 14, 15, 1);

    centerMolecule(mol);
}

// Build Losartan (C22H23ClN6O) - ARB
void buildLosartan(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Losartan (C22H23ClN6O)");

    // Imidazole-tetrazole core
    float r = 1.4f;

    // Main phenyl ring
    for (int i = 0; i < 6; i++) {
        float angle = i * PI / 3.0f;
        addAtom(mol, r * cosf(angle) - 1.5f, r * sinf(angle), 0.0f, ATOM_C);  // 0-5
    }

    // Imidazole ring attached to phenyl
    addAtom(mol, -0.5f, 1.2f, 0.0f, ATOM_C);       // 6
    addAtom(mol, 0.5f, 0.5f, 0.0f, ATOM_N);        // 7: N
    addAtom(mol, 0.5f, -0.8f, 0.0f, ATOM_C);       // 8
    addAtom(mol, -0.5f, -1.0f, 0.0f, ATOM_N);      // 9: N

    // Tetrazole ring
    addAtom(mol, 1.8f, -1.0f, 0.0f, ATOM_C);       // 10: C
    addAtom(mol, 3.0f, -0.3f, 0.0f, ATOM_N);       // 11: N
    addAtom(mol, 3.8f, -1.2f, 0.0f, ATOM_N);       // 12: N
    addAtom(mol, 3.0f, -2.3f, 0.0f, ATOM_N);       // 13: N
    addAtom(mol, 1.8f, -2.0f, 0.0f, ATOM_N);       // 14: N

    // Butyl side chain
    addAtom(mol, -2.8f, 1.2f, 0.0f, ATOM_C);       // 15: CH2
    addAtom(mol, -4.1f, 0.5f, 0.0f, ATOM_C);       // 16: CH2
    addAtom(mol, -5.3f, 1.3f, 0.0f, ATOM_C);       // 17: CH2
    addAtom(mol, -6.5f, 0.6f, 0.0f, ATOM_C);       // 18: CH3

    // Chlorine
    addAtom(mol, -2.4f, -1.4f, 0.0f, ATOM_CL);     // 19: Cl

    // Carboxylic acid
    addAtom(mol, 4.4f, 0.8f, 0.0f, ATOM_C);        // 20: C=O
    addAtom(mol, 4.4f, 0.8f, 1.3f, ATOM_O);        // 21: O
    addAtom(mol, 5.6f, 1.6f, 0.0f, ATOM_O);        // 22: OH
    addAtom(mol, 6.4f, 1.6f, 0.0f, ATOM_H);        // 23: H

    // Ring bonds
    for (int i = 0; i < 6; i++) {
        addBond(mol, i, (i + 1) % 6, (i % 2 == 0) ? 2 : 1);
    }

    // Imidazole bonds
    addBond(mol, 0, 6, 1);
    addBond(mol, 6, 7, 1);
    addBond(mol, 7, 8, 1);
    addBond(mol, 8, 9, 1);
    addBond(mol, 9, 6, 1);

    // Tetrazole bonds
    addBond(mol, 8, 10, 1);
    addBond(mol, 10, 11, 2);
    addBond(mol, 11, 12, 1);
    addBond(mol, 12, 13, 2);
    addBond(mol, 13, 14, 1);
    addBond(mol, 14, 10, 1);

    // Side chain
    addBond(mol, 6, 15, 1);
    addBond(mol, 15, 16, 1);
    addBond(mol, 16, 17, 1);
    addBond(mol, 17, 18, 1);

    // Chlorine
    addBond(mol, 4, 19, 1);

    // Carboxylic acid
    addBond(mol, 11, 20, 1);
    addBond(mol, 20, 21, 2);
    addBond(mol, 20, 22, 1);
    addBond(mol, 22, 23, 1);

    centerMolecule(mol);
}

// Build Valsartan (C24H29N5O3) - ARB
void buildValsartan(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Valsartan (C24H29N5O3)");

    // Similar ARB structure to Losartan
    float r = 1.4f;

    // Main phenyl ring
    for (int i = 0; i < 6; i++) {
        float angle = i * PI / 3.0f;
        addAtom(mol, r * cosf(angle) - 1.5f, r * sinf(angle), 0.0f, ATOM_C);  // 0-5
    }

    // Imidazole core
    addAtom(mol, -0.5f, 1.2f, 0.0f, ATOM_C);       // 6
    addAtom(mol, 0.5f, 0.5f, 0.0f, ATOM_N);        // 7: N
    addAtom(mol, 0.5f, -0.8f, 0.0f, ATOM_C);       // 8
    addAtom(mol, -0.5f, -1.0f, 0.0f, ATOM_N);      // 9: N

    // Carboxylic acid attachment
    addAtom(mol, 1.8f, -1.0f, 0.0f, ATOM_C);       // 10
    addAtom(mol, 3.0f, -0.3f, 0.0f, ATOM_C);       // 11: C=O
    addAtom(mol, 3.0f, -0.3f, 1.3f, ATOM_O);       // 12: O
    addAtom(mol, 4.2f, 0.5f, 0.0f, ATOM_O);        // 13: OH
    addAtom(mol, 5.0f, 0.5f, 0.0f, ATOM_H);        // 14: H

    // Isopropyl side chain
    addAtom(mol, -2.8f, 1.2f, 0.0f, ATOM_C);       // 15: CH
    addAtom(mol, -3.5f, 2.5f, 0.0f, ATOM_C);       // 16: CH3
    addAtom(mol, -4.0f, 0.3f, 0.0f, ATOM_C);       // 17: CH3

    // Methyltetrazole
    addAtom(mol, 1.8f, -2.3f, 0.0f, ATOM_N);       // 18: N
    addAtom(mol, 3.0f, -3.0f, 0.0f, ATOM_N);       // 19: N
    addAtom(mol, 4.2f, -2.3f, 0.0f, ATOM_N);       // 20: N
    addAtom(mol, 3.8f, -1.0f, 0.0f, ATOM_N);       // 21: N
    addAtom(mol, 0.8f, -3.0f, 0.0f, ATOM_C);       // 22: CH3

    // Ring bonds
    for (int i = 0; i < 6; i++) {
        addBond(mol, i, (i + 1) % 6, (i % 2 == 0) ? 2 : 1);
    }

    // Imidazole bonds
    addBond(mol, 0, 6, 1);
    addBond(mol, 6, 7, 1);
    addBond(mol, 7, 8, 1);
    addBond(mol, 8, 9, 1);
    addBond(mol, 9, 6, 1);

    // Carboxylic acid bonds
    addBond(mol, 8, 10, 1);
    addBond(mol, 10, 11, 1);
    addBond(mol, 11, 12, 2);
    addBond(mol, 11, 13, 1);
    addBond(mol, 13, 14, 1);

    // Isopropyl bonds
    addBond(mol, 6, 15, 1);
    addBond(mol, 15, 16, 1);
    addBond(mol, 15, 17, 1);

    // Tetrazole bonds
    addBond(mol, 10, 18, 1);
    addBond(mol, 18, 19, 1);
    addBond(mol, 19, 20, 2);
    addBond(mol, 20, 21, 1);
    addBond(mol, 21, 18, 1);
    addBond(mol, 18, 22, 1);

    centerMolecule(mol);
}

// Build Telmisartan (C33H30N4O2) - ARB
void buildTelmisartan(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Telmisartan (C33H30N4O2)");

    // Biphenyl with imidazole and benzimidazole
    float r = 1.4f;

    // First phenyl ring
    for (int i = 0; i < 6; i++) {
        float angle = i * PI / 3.0f;
        addAtom(mol, r * cosf(angle) - 2.0f, r * sinf(angle), 0.0f, ATOM_C);  // 0-5
    }

    // Second phenyl ring
    for (int i = 0; i < 6; i++) {
        float angle = i * PI / 3.0f;
        addAtom(mol, r * cosf(angle) + 2.0f, r * sinf(angle), 0.0f, ATOM_C);  // 6-11
    }

    // Imidazole core
    addAtom(mol, 0.0f, 1.4f, 0.0f, ATOM_C);        // 12
    addAtom(mol, 1.2f, 0.7f, 0.0f, ATOM_N);        // 13: N
    addAtom(mol, 1.2f, -0.6f, 0.0f, ATOM_C);       // 14
    addAtom(mol, 0.0f, -1.3f, 0.0f, ATOM_N);       // 15: N

    // Carboxylic acid
    addAtom(mol, 2.5f, -1.3f, 0.0f, ATOM_C);       // 16: C=O
    addAtom(mol, 2.5f, -1.3f, 1.3f, ATOM_O);       // 17: O
    addAtom(mol, 3.7f, -2.1f, 0.0f, ATOM_O);       // 18: OH
    addAtom(mol, 4.5f, -2.1f, 0.0f, ATOM_H);       // 19: H

    // Propyl side chain
    addAtom(mol, -1.2f, 2.3f, 0.0f, ATOM_C);       // 20: CH2
    addAtom(mol, -2.5f, 3.0f, 0.0f, ATOM_C);       // 21: CH2
    addAtom(mol, -3.7f, 2.3f, 0.0f, ATOM_C);       // 22: CH3

    // Ring bonds
    for (int i = 0; i < 6; i++) {
        addBond(mol, i, (i + 1) % 6, (i % 2 == 0) ? 2 : 1);
        addBond(mol, 6 + i, 6 + ((i + 1) % 6), (i % 2 == 0) ? 2 : 1);
    }

    // Ring connection
    addBond(mol, 0, 6, 1);

    // Imidazole bonds
    addBond(mol, 0, 12, 1);
    addBond(mol, 12, 13, 1);
    addBond(mol, 13, 14, 1);
    addBond(mol, 14, 15, 1);
    addBond(mol, 15, 12, 1);

    // Carboxylic acid bonds
    addBond(mol, 14, 16, 1);
    addBond(mol, 16, 17, 2);
    addBond(mol, 16, 18, 1);
    addBond(mol, 18, 19, 1);

    // Side chain
    addBond(mol, 12, 20, 1);
    addBond(mol, 20, 21, 1);
    addBond(mol, 21, 22, 1);

    centerMolecule(mol);
}

// Build Metoprolol (C15H25NO3) - Beta Blocker
void buildMetoprolol(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Metoprolol (C15H25NO3)");

    // Catechol-like core (two hydroxy groups on benzene)
    float r = 1.4f;
    for (int i = 0; i < 6; i++) {
        float angle = i * PI / 3.0f;
        addAtom(mol, r * cosf(angle), r * sinf(angle), 0.0f, ATOM_C);  // 0-5
    }

    // Hydroxyl groups (beta-blocker characteristic)
    addAtom(mol, r * cosf(PI/3) + 0.8f, r * sinf(PI/3) + 0.8f, 0.0f, ATOM_O);     // 6: OH
    addAtom(mol, r * cosf(PI/3) + 1.6f, r * sinf(PI/3) + 1.6f, 0.0f, ATOM_H);     // 7: H
    addAtom(mol, r * cosf(2*PI/3) - 0.8f, r * sinf(2*PI/3) + 0.8f, 0.0f, ATOM_O); // 8: OH
    addAtom(mol, r * cosf(2*PI/3) - 1.6f, r * sinf(2*PI/3) + 1.6f, 0.0f, ATOM_H); // 9: H

    // Ethyl ether side chain
    addAtom(mol, -r - 1.3f, 0.0f, 0.0f, ATOM_O);   // 10: O (ether)
    addAtom(mol, -r - 2.6f, 0.0f, 0.0f, ATOM_C);   // 11: CH2
    addAtom(mol, -r - 3.9f, 0.0f, 0.0f, ATOM_C);   // 12: CH3

    // Isopropyl amine side chain
    addAtom(mol, r + 1.3f, 0.0f, 0.0f, ATOM_C);    // 13: CH
    addAtom(mol, r + 1.5f, 1.4f, 0.0f, ATOM_C);    // 14: CH3
    addAtom(mol, r + 2.6f, -0.7f, 0.0f, ATOM_C);   // 15: CH3
    addAtom(mol, r + 2.6f, 0.7f, 0.0f, ATOM_N);    // 16: N
    addAtom(mol, r + 3.9f, 1.4f, 0.0f, ATOM_C);    // 17: CH3
    addAtom(mol, r + 2.6f, 1.8f, 0.0f, ATOM_H);    // 18: H on N

    // Ring bonds
    for (int i = 0; i < 6; i++) {
        addBond(mol, i, (i + 1) % 6, (i % 2 == 0) ? 2 : 1);
    }

    // Hydroxyl bonds
    addBond(mol, 1, 6, 1);
    addBond(mol, 6, 7, 1);
    addBond(mol, 2, 8, 1);
    addBond(mol, 8, 9, 1);

    // Ether side chain
    addBond(mol, 5, 10, 1);
    addBond(mol, 10, 11, 1);
    addBond(mol, 11, 12, 1);

    // Amine side chain
    addBond(mol, 0, 13, 1);
    addBond(mol, 13, 14, 1);
    addBond(mol, 13, 15, 1);
    addBond(mol, 13, 16, 1);
    addBond(mol, 16, 17, 1);
    addBond(mol, 16, 18, 1);

    centerMolecule(mol);
}

// Build Atenolol (C14H22N2O3) - Beta Blocker
void buildAtenolol(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Atenolol (C14H22N2O3)");

    // Para-substituted benzene with hydroxyl and amine side chains
    float r = 1.4f;
    for (int i = 0; i < 6; i++) {
        float angle = i * PI / 3.0f;
        addAtom(mol, r * cosf(angle), r * sinf(angle), 0.0f, ATOM_C);  // 0-5
    }

    // Hydroxyl group at para position
    addAtom(mol, r + 0.8f, 0.0f, 0.0f, ATOM_O);    // 6: OH
    addAtom(mol, r + 1.6f, 0.0f, 0.0f, ATOM_H);    // 7: H

    // Propyl amine side chain
    addAtom(mol, -r - 1.3f, 0.0f, 0.0f, ATOM_C);   // 8: CH2
    addAtom(mol, -r - 2.6f, 0.0f, 0.0f, ATOM_C);   // 9: CH2
    addAtom(mol, -r - 3.9f, 0.0f, 0.0f, ATOM_C);   // 10: CH2
    addAtom(mol, -r - 5.2f, 0.0f, 0.0f, ATOM_N);   // 11: NH2
    addAtom(mol, -r - 6.1f, 0.7f, 0.0f, ATOM_H);   // 12: H on N
    addAtom(mol, -r - 6.1f, -0.7f, 0.0f, ATOM_H);  // 13: H on N

    // Carbamate group
    addAtom(mol, 0.0f, -r - 1.3f, 0.0f, ATOM_C);   // 14: C=O
    addAtom(mol, 0.0f, -r - 1.3f, 1.3f, ATOM_O);   // 15: O
    addAtom(mol, 1.2f, -r - 2.0f, 0.0f, ATOM_N);   // 16: N
    addAtom(mol, 2.0f, -r - 2.0f, 0.0f, ATOM_H);   // 17: H

    // Ring bonds
    for (int i = 0; i < 6; i++) {
        addBond(mol, i, (i + 1) % 6, (i % 2 == 0) ? 2 : 1);
    }

    // Hydroxyl bond
    addBond(mol, 0, 6, 1);
    addBond(mol, 6, 7, 1);

    // Propyl amine chain
    addBond(mol, 5, 8, 1);
    addBond(mol, 8, 9, 1);
    addBond(mol, 9, 10, 1);
    addBond(mol, 10, 11, 1);
    addBond(mol, 11, 12, 1);
    addBond(mol, 11, 13, 1);

    // Carbamate
    addBond(mol, 3, 14, 1);
    addBond(mol, 14, 15, 2);
    addBond(mol, 14, 16, 1);
    addBond(mol, 16, 17, 1);

    centerMolecule(mol);
}

// Build Carvedilol (C24H26N2O4) - Beta Blocker
void buildCarvedilol(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Carvedilol (C24H26N2O4)");

    // Carbazole-like core (tricyclic aromatic)
    float r = 1.4f;

    // First aromatic ring
    for (int i = 0; i < 6; i++) {
        float angle = i * PI / 3.0f;
        addAtom(mol, r * cosf(angle) - 1.0f, r * sinf(angle), 0.0f, ATOM_C);  // 0-5
    }

    // Second aromatic ring (fused)
    for (int i = 0; i < 6; i++) {
        float angle = i * PI / 3.0f;
        addAtom(mol, r * cosf(angle) + 1.0f, r * sinf(angle), 0.0f, ATOM_C);  // 6-11
    }

    // Propyl ether side chain
    addAtom(mol, 0.0f, -r - 1.3f, 0.0f, ATOM_O);   // 12: O (ether)
    addAtom(mol, -1.3f, -r - 2.0f, 0.0f, ATOM_C);  // 13: CH2
    addAtom(mol, -2.6f, -r - 1.3f, 0.0f, ATOM_C);  // 14: CH2
    addAtom(mol, -3.9f, -r - 2.0f, 0.0f, ATOM_C);  // 15: CH2

    // Secondary amine group
    addAtom(mol, -5.2f, -r - 1.3f, 0.0f, ATOM_N);  // 16: NH
    addAtom(mol, -6.5f, -r - 2.0f, 0.0f, ATOM_C);  // 17: CH (isopropyl)
    addAtom(mol, -7.2f, -r - 0.6f, 0.0f, ATOM_C);  // 18: CH3
    addAtom(mol, -7.2f, -r - 3.4f, 0.0f, ATOM_C);  // 19: CH3
    addAtom(mol, -6.5f, -r - 2.0f, -1.0f, ATOM_H); // 20: H on N

    // Hydroxyalkoxy side chain
    addAtom(mol, 0.0f, r + 1.3f, 0.0f, ATOM_O);    // 21: O
    addAtom(mol, 1.3f, r + 2.0f, 0.0f, ATOM_C);    // 22: CH
    addAtom(mol, 2.6f, r + 1.3f, 0.0f, ATOM_C);    // 23: CH2
    addAtom(mol, 2.6f, r + 1.3f, 1.3f, ATOM_O);    // 24: OH
    addAtom(mol, 3.4f, r + 1.3f, 1.8f, ATOM_H);    // 25: H on OH
    addAtom(mol, 1.3f, r + 3.4f, 0.0f, ATOM_C);    // 26: CH3

    // Ring bonds
    for (int i = 0; i < 6; i++) {
        addBond(mol, i, (i + 1) % 6, (i % 2 == 0) ? 2 : 1);
        addBond(mol, 6 + i, 6 + ((i + 1) % 6), (i % 2 == 0) ? 2 : 1);
    }

    // Ring connection
    addBond(mol, 0, 6, 1);

    // Ether side chain
    addBond(mol, 3, 12, 1);
    addBond(mol, 12, 13, 1);
    addBond(mol, 13, 14, 1);
    addBond(mol, 14, 15, 1);
    addBond(mol, 15, 16, 1);
    addBond(mol, 16, 17, 1);
    addBond(mol, 17, 18, 1);
    addBond(mol, 17, 19, 1);
    addBond(mol, 16, 20, 1);

    // Hydroxyalkoxy side chain
    addBond(mol, 2, 21, 1);
    addBond(mol, 21, 22, 1);
    addBond(mol, 22, 23, 1);
    addBond(mol, 23, 24, 1);
    addBond(mol, 24, 25, 1);
    addBond(mol, 22, 26, 1);

    centerMolecule(mol);
}

// Build Amlodipine (C26H33ClN2O8) - Calcium Channel Blocker
void buildAmlodipine(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Amlodipine (C26H33ClN2O8)");

    // Dihydropyridine ring with two phenyl groups
    // Central dihydropyridine ring
    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_C);        // 0
    addAtom(mol, 1.3f, 0.7f, 0.0f, ATOM_C);        // 1
    addAtom(mol, 1.3f, 2.1f, 0.0f, ATOM_C);        // 2
    addAtom(mol, 0.0f, 2.8f, 0.0f, ATOM_N);        // 3: N
    addAtom(mol, -1.3f, 2.1f, 0.0f, ATOM_C);       // 4
    addAtom(mol, -1.3f, 0.7f, 0.0f, ATOM_C);       // 5

    // Phenyl ring 1 (2-chlorophenyl)
    float r = 1.4f;
    float baseX1 = 2.8f, baseY1 = 3.5f;
    for (int i = 0; i < 6; i++) {
        float angle = i * PI / 3.0f;
        addAtom(mol, baseX1 + r * cosf(angle), baseY1 + r * sinf(angle), 0.0f, ATOM_C);  // 6-11
    }

    // Chlorine on first phenyl
    addAtom(mol, baseX1 + r * cosf(PI/3) + 0.8f, baseY1 + r * sinf(PI/3) + 0.8f, 0.0f, ATOM_CL);  // 12

    // Phenyl ring 2
    float baseX2 = -2.8f, baseY2 = 3.5f;
    for (int i = 0; i < 6; i++) {
        float angle = i * PI / 3.0f;
        addAtom(mol, baseX2 + r * cosf(angle), baseY2 + r * sinf(angle), 0.0f, ATOM_C);  // 13-18
    }

    // Ester groups
    addAtom(mol, 2.6f, -1.2f, 0.0f, ATOM_C);       // 19: C=O (ester)
    addAtom(mol, 2.6f, -1.2f, 1.3f, ATOM_O);       // 20: O
    addAtom(mol, 3.8f, -1.9f, 0.0f, ATOM_O);       // 21: O (ester)
    addAtom(mol, 5.1f, -1.3f, 0.0f, ATOM_C);       // 22: CH3

    addAtom(mol, -2.6f, -1.2f, 0.0f, ATOM_C);      // 23: C=O (ester)
    addAtom(mol, -2.6f, -1.2f, 1.3f, ATOM_O);      // 24: O
    addAtom(mol, -3.8f, -1.9f, 0.0f, ATOM_O);      // 25: O (ester)
    addAtom(mol, -5.1f, -1.3f, 0.0f, ATOM_C);      // 26: CH3

    // Central ring bonds
    addBond(mol, 0, 1, 1);
    addBond(mol, 1, 2, 1);
    addBond(mol, 2, 3, 1);
    addBond(mol, 3, 4, 1);
    addBond(mol, 4, 5, 1);
    addBond(mol, 5, 0, 1);

    // Phenyl ring 1 bonds
    for (int i = 0; i < 6; i++) {
        addBond(mol, 6 + i, 6 + ((i + 1) % 6), (i % 2 == 0) ? 2 : 1);
    }

    // Phenyl ring 2 bonds
    for (int i = 0; i < 6; i++) {
        addBond(mol, 13 + i, 13 + ((i + 1) % 6), (i % 2 == 0) ? 2 : 1);
    }

    // Phenyl attachments
    addBond(mol, 2, 6, 1);
    addBond(mol, 4, 13, 1);

    // Chlorine
    addBond(mol, 8, 12, 1);

    // Ester bonds 1
    addBond(mol, 0, 19, 1);
    addBond(mol, 19, 20, 2);
    addBond(mol, 19, 21, 1);
    addBond(mol, 21, 22, 1);

    // Ester bonds 2
    addBond(mol, 5, 23, 1);
    addBond(mol, 23, 24, 2);
    addBond(mol, 23, 25, 1);
    addBond(mol, 25, 26, 1);

    centerMolecule(mol);
}

// Build Diltiazem (C26H36N2O4S) - Calcium Channel Blocker
void buildDiltiazem(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Diltiazem (C26H36N2O4S)");

    // Benzothiazepine-like tricyclic core
    float r = 1.4f;

    // Main ring system
    for (int i = 0; i < 6; i++) {
        float angle = i * PI / 3.0f;
        addAtom(mol, r * cosf(angle), r * sinf(angle) + 1.0f, 0.0f, ATOM_C);  // 0-5
    }

    // Seven-membered middle ring with S and N
    addAtom(mol, 0.0f, 1.0f, 0.0f, ATOM_S);        // 6: S
    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_C);        // 7
    addAtom(mol, -1.3f, -0.7f, 0.0f, ATOM_N);      // 8: N
    addAtom(mol, -1.3f, -2.0f, 0.0f, ATOM_C);      // 9
    addAtom(mol, 0.0f, -2.7f, 0.0f, ATOM_C);       // 10
    addAtom(mol, 1.3f, -2.0f, 0.0f, ATOM_C);       // 11
    addAtom(mol, 1.3f, -0.7f, 0.0f, ATOM_C);       // 12

    // Second aromatic ring
    for (int i = 0; i < 6; i++) {
        float angle = i * PI / 3.0f;
        addAtom(mol, r * cosf(angle), r * sinf(angle) - 1.5f, 0.0f, ATOM_C);  // 13-18
    }

    // Acetyl amine side chain
    addAtom(mol, -2.6f, -2.8f, 0.0f, ATOM_C);      // 19: C=O (acetyl)
    addAtom(mol, -2.6f, -2.8f, 1.3f, ATOM_O);      // 20: O
    addAtom(mol, -3.8f, -3.5f, 0.0f, ATOM_C);      // 21: CH3

    // N-methyl group
    addAtom(mol, -2.0f, -0.2f, 0.0f, ATOM_C);      // 22: CH3

    // Methoxy groups
    addAtom(mol, 2.8f, 1.9f, 0.0f, ATOM_O);        // 23: OCH3
    addAtom(mol, 4.1f, 2.6f, 0.0f, ATOM_C);        // 24: CH3
    addAtom(mol, -2.8f, 1.9f, 0.0f, ATOM_O);       // 25: OCH3
    addAtom(mol, -4.1f, 2.6f, 0.0f, ATOM_C);       // 26: CH3

    // Ring bonds
    for (int i = 0; i < 6; i++) {
        addBond(mol, i, (i + 1) % 6, (i % 2 == 0) ? 2 : 1);
        addBond(mol, 13 + i, 13 + ((i + 1) % 6), (i % 2 == 0) ? 2 : 1);
    }

    // Middle ring bonds
    addBond(mol, 0, 6, 1);
    addBond(mol, 6, 7, 1);
    addBond(mol, 7, 8, 1);
    addBond(mol, 8, 9, 1);
    addBond(mol, 9, 10, 1);
    addBond(mol, 10, 11, 1);
    addBond(mol, 11, 12, 1);
    addBond(mol, 12, 6, 1);

    // Ring connections
    addBond(mol, 0, 13, 1);
    addBond(mol, 12, 13, 1);

    // Acetyl amine bonds
    addBond(mol, 9, 19, 1);
    addBond(mol, 19, 20, 2);
    addBond(mol, 19, 21, 1);

    // N-methyl bond
    addBond(mol, 8, 22, 1);

    // Methoxy bonds
    addBond(mol, 1, 23, 1);
    addBond(mol, 23, 24, 1);
    addBond(mol, 5, 25, 1);
    addBond(mol, 25, 26, 1);

    centerMolecule(mol);
}

// Build Nifedipine (C17H18N2O6) - Calcium Channel Blocker
void buildNifedipine(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Nifedipine (C17H18N2O6)");

    // Dihydropyridine ring with nitro and ester groups
    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_C);        // 0
    addAtom(mol, 1.3f, 0.7f, 0.0f, ATOM_C);        // 1
    addAtom(mol, 1.3f, 2.1f, 0.0f, ATOM_C);        // 2
    addAtom(mol, 0.0f, 2.8f, 0.0f, ATOM_N);        // 3: N
    addAtom(mol, -1.3f, 2.1f, 0.0f, ATOM_C);       // 4
    addAtom(mol, -1.3f, 0.7f, 0.0f, ATOM_C);       // 5

    // Phenyl ring with nitro group
    float r = 1.4f;
    float baseX = 2.8f, baseY = 3.5f;
    for (int i = 0; i < 6; i++) {
        float angle = i * PI / 3.0f;
        addAtom(mol, baseX + r * cosf(angle), baseY + r * sinf(angle), 0.0f, ATOM_C);  // 6-11
    }

    // Nitro group (NO2) on phenyl
    addAtom(mol, baseX + r * cosf(2*PI/3) - 0.8f, baseY + r * sinf(2*PI/3) - 0.8f, 0.0f, ATOM_N);  // 12: N
    addAtom(mol, baseX + r * cosf(2*PI/3) - 0.8f, baseY + r * sinf(2*PI/3) - 2.1f, 0.0f, ATOM_O);  // 13: O
    addAtom(mol, baseX + r * cosf(2*PI/3) + 0.5f, baseY + r * sinf(2*PI/3) - 0.2f, 0.0f, ATOM_O);  // 14: O

    // Methyl ester at position 1
    addAtom(mol, 2.6f, -1.2f, 0.0f, ATOM_C);       // 15: C=O
    addAtom(mol, 2.6f, -1.2f, 1.3f, ATOM_O);       // 16: O
    addAtom(mol, 3.8f, -1.9f, 0.0f, ATOM_O);       // 17: O (ester)
    addAtom(mol, 5.1f, -1.3f, 0.0f, ATOM_C);       // 18: CH3

    // Methyl ester at position 3
    addAtom(mol, -2.6f, -1.2f, 0.0f, ATOM_C);      // 19: C=O
    addAtom(mol, -2.6f, -1.2f, 1.3f, ATOM_O);      // 20: O
    addAtom(mol, -3.8f, -1.9f, 0.0f, ATOM_O);      // 21: O (ester)
    addAtom(mol, -5.1f, -1.3f, 0.0f, ATOM_C);      // 22: CH3

    // Central ring bonds
    addBond(mol, 0, 1, 1);
    addBond(mol, 1, 2, 1);
    addBond(mol, 2, 3, 1);
    addBond(mol, 3, 4, 1);
    addBond(mol, 4, 5, 1);
    addBond(mol, 5, 0, 1);

    // Phenyl ring bonds
    for (int i = 0; i < 6; i++) {
        addBond(mol, 6 + i, 6 + ((i + 1) % 6), (i % 2 == 0) ? 2 : 1);
    }

    // Phenyl attachment
    addBond(mol, 2, 6, 1);

    // Nitro bonds
    addBond(mol, 8, 12, 1);
    addBond(mol, 12, 13, 2);
    addBond(mol, 12, 14, 2);

    // Ester 1 bonds
    addBond(mol, 0, 15, 1);
    addBond(mol, 15, 16, 2);
    addBond(mol, 15, 17, 1);
    addBond(mol, 17, 18, 1);

    // Ester 2 bonds
    addBond(mol, 5, 19, 1);
    addBond(mol, 19, 20, 2);
    addBond(mol, 19, 21, 1);
    addBond(mol, 21, 22, 1);

    centerMolecule(mol);
}

// Build Hydrochlorothiazide/HCTZ (C7H8ClN3O4S2) - Diuretic
void buildHydrochlorothiazide(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Hydrochlorothiazide/HCTZ (C7H8ClN3O4S2)");

    // Benzothiazide core with sulfonamide
    float r = 1.4f;
    for (int i = 0; i < 6; i++) {
        float angle = i * PI / 3.0f;
        addAtom(mol, r * cosf(angle), r * sinf(angle), 0.0f, ATOM_C);  // 0-5
    }

    // Thiazide ring
    addAtom(mol, 0.0f, -1.4f, 0.0f, ATOM_S);       // 6: S
    addAtom(mol, 1.3f, -2.1f, 0.0f, ATOM_C);       // 7: C=O
    addAtom(mol, 1.3f, -2.1f, 1.3f, ATOM_O);       // 8: O
    addAtom(mol, 2.6f, -2.8f, 0.0f, ATOM_N);       // 9: N
    addAtom(mol, 3.9f, -2.1f, 0.0f, ATOM_C);       // 10: Cl attachment point

    // Chlorine substituent
    addAtom(mol, 5.2f, -2.8f, 0.0f, ATOM_CL);      // 11: Cl

    // Sulfonamide group SO2NH2
    addAtom(mol, 0.0f, -2.8f, 0.0f, ATOM_S);       // 12: S (sulfonamide)
    addAtom(mol, 0.0f, -4.1f, 0.0f, ATOM_O);       // 13: O (sulfonyl)
    addAtom(mol, 1.3f, -3.4f, 0.0f, ATOM_O);       // 14: O (sulfonyl)
    addAtom(mol, -1.3f, -3.4f, 0.0f, ATOM_N);      // 15: NH2
    addAtom(mol, -1.8f, -4.3f, 0.0f, ATOM_H);      // 16: H
    addAtom(mol, -2.1f, -2.8f, 0.0f, ATOM_H);      // 17: H

    // Benzene ring bonds
    for (int i = 0; i < 6; i++) {
        addBond(mol, i, (i + 1) % 6, (i % 2 == 0) ? 2 : 1);
    }

    // Thiazide ring bonds
    addBond(mol, 2, 6, 1);
    addBond(mol, 6, 7, 1);
    addBond(mol, 7, 8, 2);
    addBond(mol, 7, 9, 1);
    addBond(mol, 9, 10, 1);
    addBond(mol, 10, 5, 1);

    // Chlorine
    addBond(mol, 10, 11, 1);

    // Sulfonamide
    addBond(mol, 3, 12, 1);
    addBond(mol, 12, 13, 2);
    addBond(mol, 12, 14, 2);
    addBond(mol, 12, 15, 1);
    addBond(mol, 15, 16, 1);
    addBond(mol, 15, 17, 1);

    centerMolecule(mol);
}

// Build Furosemide/Lasix (C12H11ClN2O5S) - Diuretic
void buildFurosemide(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Furosemide/Lasix (C12H11ClN2O5S)");

    // Arylsulfonylurea core
    float r = 1.4f;

    // First phenyl ring (with Cl)
    for (int i = 0; i < 6; i++) {
        float angle = i * PI / 3.0f;
        addAtom(mol, r * cosf(angle) - 2.0f, r * sinf(angle), 0.0f, ATOM_C);  // 0-5
    }

    // Chlorine on first ring
    addAtom(mol, -2.0f + r * cosf(2*PI/3) - 0.8f, r * sinf(2*PI/3) - 0.8f, 0.0f, ATOM_CL);  // 6

    // Sulfone bridge
    addAtom(mol, 0.0f, 1.4f, 0.0f, ATOM_S);        // 7: S
    addAtom(mol, 0.0f, 1.4f, 1.3f, ATOM_O);        // 8: O (sulfonyl)
    addAtom(mol, 1.3f, 2.1f, 0.0f, ATOM_O);        // 9: O (sulfonyl)

    // Second phenyl ring (with amino group)
    float baseX2 = 2.5f;
    for (int i = 0; i < 6; i++) {
        float angle = i * PI / 3.0f;
        addAtom(mol, baseX2 + r * cosf(angle), r * sinf(angle), 0.0f, ATOM_C);  // 10-15
    }

    // Amino group NH2 on second ring
    addAtom(mol, baseX2 + r * cosf(2*PI/3) - 0.8f, r * sinf(2*PI/3) - 0.8f, 0.0f, ATOM_N);  // 16: N
    addAtom(mol, baseX2 + r * cosf(2*PI/3) - 1.3f, r * sinf(2*PI/3) - 1.8f, 0.0f, ATOM_H);  // 17: H
    addAtom(mol, baseX2 + r * cosf(2*PI/3) + 0.5f, r * sinf(2*PI/3) - 0.2f, 0.0f, ATOM_H);  // 18: H

    // Carboxylic acid group
    addAtom(mol, baseX2 + r + 1.3f, 0.0f, 0.0f, ATOM_C);      // 19: C=O
    addAtom(mol, baseX2 + r + 1.3f, 0.0f, 1.3f, ATOM_O);      // 20: O
    addAtom(mol, baseX2 + r + 2.6f, 0.0f, 0.0f, ATOM_O);      // 21: OH
    addAtom(mol, baseX2 + r + 3.4f, 0.0f, 0.0f, ATOM_H);      // 22: H

    // Ring 1 bonds
    for (int i = 0; i < 6; i++) {
        addBond(mol, i, (i + 1) % 6, (i % 2 == 0) ? 2 : 1);
    }

    // Ring 2 bonds
    for (int i = 0; i < 6; i++) {
        addBond(mol, 10 + i, 10 + ((i + 1) % 6), (i % 2 == 0) ? 2 : 1);
    }

    // Chlorine
    addBond(mol, 2, 6, 1);

    // Sulfone bridge
    addBond(mol, 5, 7, 1);
    addBond(mol, 7, 8, 2);
    addBond(mol, 7, 9, 2);
    addBond(mol, 7, 10, 1);

    // Amino group
    addBond(mol, 12, 16, 1);
    addBond(mol, 16, 17, 1);
    addBond(mol, 16, 18, 1);

    // Carboxylic acid
    addBond(mol, 15, 19, 1);
    addBond(mol, 19, 20, 2);
    addBond(mol, 19, 21, 1);
    addBond(mol, 21, 22, 1);

    centerMolecule(mol);
}

// Build Spironolactone (C24H32O4S) - Diuretic
void buildSpironolactone(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Spironolactone (C24H32O4S)");

    // Steroid-like core with lactone ring
    // Simplified testosterone-like backbone
    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_C);        // 0: A ring
    addAtom(mol, 1.3f, 0.7f, 0.0f, ATOM_C);        // 1
    addAtom(mol, 2.6f, 0.0f, 0.0f, ATOM_C);        // 2
    addAtom(mol, 1.3f, -1.3f, 0.0f, ATOM_C);       // 3
    addAtom(mol, 0.0f, -0.6f, 0.0f, ATOM_C);       // 4

    // B ring
    addAtom(mol, 1.3f, 2.1f, 0.0f, ATOM_C);        // 5
    addAtom(mol, 2.6f, 2.8f, 0.0f, ATOM_C);        // 6
    addAtom(mol, 3.9f, 2.1f, 0.0f, ATOM_C);        // 7
    addAtom(mol, 2.6f, 0.7f, 0.0f, ATOM_C);        // 8: shared

    // C ring
    addAtom(mol, 3.9f, 2.1f, 0.0f, ATOM_C);        // 9
    addAtom(mol, 5.2f, 2.8f, 0.0f, ATOM_C);        // 10
    addAtom(mol, 5.2f, 1.4f, 0.0f, ATOM_C);        // 11

    // D ring with ketone and hydroxyl
    addAtom(mol, 3.9f, 0.7f, 0.0f, ATOM_C);        // 12: ketone C
    addAtom(mol, 3.9f, 0.7f, 1.3f, ATOM_O);        // 13: C=O (ketone)
    addAtom(mol, 5.2f, 0.0f, 0.0f, ATOM_C);        // 14: C-OH
    addAtom(mol, 6.5f, 0.7f, 0.0f, ATOM_O);        // 15: OH
    addAtom(mol, 7.3f, 0.7f, 0.0f, ATOM_H);        // 16: H on OH

    // Lactone ring side chain
    addAtom(mol, 2.6f, -2.6f, 0.0f, ATOM_C);       // 17: CH2
    addAtom(mol, 3.9f, -3.3f, 0.0f, ATOM_C);       // 18: CH2
    addAtom(mol, 5.2f, -2.6f, 0.0f, ATOM_C);       // 19: C=O (lactone)
    addAtom(mol, 5.2f, -2.6f, 1.3f, ATOM_O);       // 20: O (carbonyl)
    addAtom(mol, 6.5f, -3.3f, 0.0f, ATOM_O);       // 21: O (ether in ring)
    addAtom(mol, 6.5f, -2.0f, 0.0f, ATOM_C);       // 22: CH2 (lactone ring)

    // Acetyl group (with thio-link)
    addAtom(mol, 7.8f, -2.7f, 0.0f, ATOM_C);       // 23: CH3
    addAtom(mol, 7.1f, -2.7f, 0.8f, ATOM_S);       // 24: S between 22 and CH3 (simplified thioether)

    // Hydrogens (approximate placement)
    // CH3 on acetyl
    addAtom(mol, 8.6f, -2.2f, 0.6f, ATOM_H);       // 25
    addAtom(mol, 8.6f, -3.4f, 0.0f, ATOM_H);       // 26
    addAtom(mol, 8.2f, -2.2f, -0.8f, ATOM_H);      // 27
    // CH2 groups in lactone chain
    addAtom(mol, 2.6f, -2.6f, 1.0f, ATOM_H);       // 28 (C17)
    addAtom(mol, 2.0f, -3.1f, -0.6f, ATOM_H);      // 29 (C17)
    addAtom(mol, 3.9f, -3.3f, 1.0f, ATOM_H);       // 30 (C18)
    addAtom(mol, 3.3f, -3.8f, -0.6f, ATOM_H);      // 31 (C18)
    addAtom(mol, 6.5f, -2.0f, 1.0f, ATOM_H);       // 32 (C22)
    addAtom(mol, 6.5f, -1.4f, -0.6f, ATOM_H);      // 33 (C22)
    // A/B rings
    addAtom(mol, -0.8f, 0.5f, 0.0f, ATOM_H);       // 34 (C0)
    addAtom(mol, -0.6f, -1.0f, 0.0f, ATOM_H);      // 35 (C0)
    addAtom(mol, 1.3f, 0.7f, 1.0f, ATOM_H);        // 36 (C1)
    addAtom(mol, 2.6f, 0.0f, 1.0f, ATOM_H);        // 37 (C2)
    addAtom(mol, 1.3f, -1.3f, 1.0f, ATOM_H);       // 38 (C3)
    addAtom(mol, 1.3f, -1.3f, -1.0f, ATOM_H);      // 39 (C3)
    addAtom(mol, -0.4f, -0.9f, 1.0f, ATOM_H);      // 40 (C4)
    addAtom(mol, 0.9f, 2.6f, 1.0f, ATOM_H);        // 41 (C5)
    addAtom(mol, 0.9f, 2.6f, -1.0f, ATOM_H);       // 42 (C5)
    addAtom(mol, 2.6f, 2.8f, 1.0f, ATOM_H);        // 43 (C6)
    addAtom(mol, 3.9f, 2.1f, 1.0f, ATOM_H);        // 44 (C7)
    addAtom(mol, 2.6f, 0.7f, 1.0f, ATOM_H);        // 45 (C8)
    // C ring
    addAtom(mol, 3.9f, 3.0f, 0.8f, ATOM_H);        // 46 (C9)
    addAtom(mol, 5.2f, 2.8f, 1.0f, ATOM_H);        // 47 (C10)
    addAtom(mol, 5.2f, 2.8f, -1.0f, ATOM_H);       // 48 (C10)
    addAtom(mol, 5.2f, 1.4f, 1.0f, ATOM_H);        // 49 (C11)
    addAtom(mol, 5.2f, 1.4f, -1.0f, ATOM_H);       // 50 (C11)
    // D ring side carbon
    addAtom(mol, 5.2f, 0.0f, 1.0f, ATOM_H);        // 51 (C14)
    addAtom(mol, 5.8f, -0.6f, -0.6f, ATOM_H);      // 52 (C14)
    // Extra hydrogens to match formula
    addAtom(mol, 2.6f, 2.0f, -1.0f, ATOM_H);       // 53 (C9 extra)
    addAtom(mol, 2.6f, 2.8f, -1.2f, ATOM_H);       // 54 (C6 extra)
    addAtom(mol, 2.6f, 0.3f, -1.0f, ATOM_H);       // 55 (C8 extra)

    // Ring bonds (simplified)
    addBond(mol, 0, 1, 1);
    addBond(mol, 1, 2, 1);
    addBond(mol, 2, 3, 1);
    addBond(mol, 3, 4, 1);
    addBond(mol, 4, 0, 1);
    addBond(mol, 1, 5, 1);
    addBond(mol, 5, 6, 1);
    addBond(mol, 6, 7, 1);
    addBond(mol, 7, 8, 1);
    addBond(mol, 8, 2, 1);

    // Complete additional ring bonds (simplified connectivity)
    addBond(mol, 7, 9, 1);
    addBond(mol, 9, 10, 1);
    addBond(mol, 10, 11, 1);
    addBond(mol, 11, 12, 1);
    addBond(mol, 12, 8, 1);

    // Ketone bond
    addBond(mol, 12, 13, 2);

    // Hydroxyl bonds
    addBond(mol, 14, 15, 1);
    addBond(mol, 15, 16, 1);
    addBond(mol, 12, 14, 1);

    // Lactone side chain
    addBond(mol, 4, 17, 1);
    addBond(mol, 17, 18, 1);
    addBond(mol, 18, 19, 1);
    addBond(mol, 19, 20, 2);
    addBond(mol, 19, 21, 1);
    addBond(mol, 21, 22, 1);
    addBond(mol, 22, 18, 1);
    addBond(mol, 22, 24, 1);   // 22-S
    addBond(mol, 24, 23, 1);   // S-CH3

    // H bonds
    addBond(mol, 23, 25, 1); addBond(mol, 23, 26, 1); addBond(mol, 23, 27, 1);
    addBond(mol, 17, 28, 1); addBond(mol, 17, 29, 1);
    addBond(mol, 18, 30, 1); addBond(mol, 18, 31, 1);
    addBond(mol, 22, 32, 1); addBond(mol, 22, 33, 1);
    addBond(mol, 0, 34, 1); addBond(mol, 0, 35, 1);
    addBond(mol, 1, 36, 1);
    addBond(mol, 2, 37, 1);
    addBond(mol, 3, 38, 1); addBond(mol, 3, 39, 1);
    addBond(mol, 4, 40, 1);
    addBond(mol, 5, 41, 1); addBond(mol, 5, 42, 1);
    addBond(mol, 6, 43, 1);
    addBond(mol, 7, 44, 1);
    addBond(mol, 8, 45, 1);
    addBond(mol, 9, 46, 1); addBond(mol, 10, 47, 1); addBond(mol, 10, 48, 1);
    addBond(mol, 11, 49, 1); addBond(mol, 11, 50, 1);
    addBond(mol, 14, 51, 1); addBond(mol, 14, 52, 1);
    addBond(mol, 9, 53, 1);
    addBond(mol, 6, 54, 1);
    addBond(mol, 8, 55, 1);

    centerMolecule(mol);
}

// Build Oseltamivir/Tamiflu (C15H31N3O8) - Antiviral
void buildOseltamivir(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Oseltamivir/Tamiflu (C15H31N3O8)");

    // Cyclohexene core with side chains
    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_C);        // 0
    addAtom(mol, 1.3f, 0.7f, 0.0f, ATOM_C);        // 1
    addAtom(mol, 2.6f, 0.0f, 0.0f, ATOM_C);        // 2: C=
    addAtom(mol, 2.6f, -1.4f, 0.0f, ATOM_C);       // 3: =C
    addAtom(mol, 1.3f, -2.1f, 0.0f, ATOM_C);       // 4
    addAtom(mol, 0.0f, -1.4f, 0.0f, ATOM_C);       // 5

    // Amino group on ring
    addAtom(mol, -1.3f, 0.7f, 0.0f, ATOM_N);       // 6: N
    addAtom(mol, -2.6f, 0.0f, 0.0f, ATOM_C);       // 7: CH3
    addAtom(mol, -1.5f, 2.1f, 0.0f, ATOM_H);       // 8: H on N

    // Acetyl ester side chain
    addAtom(mol, 3.9f, 0.7f, 0.0f, ATOM_O);        // 9: O (ether)
    addAtom(mol, 5.2f, 0.0f, 0.0f, ATOM_C);        // 10: C=O
    addAtom(mol, 5.2f, 0.0f, 1.3f, ATOM_O);        // 11: O (carbonyl)
    addAtom(mol, 6.5f, -0.7f, 0.0f, ATOM_C);       // 12: CH3

    // Carboxylic acid side chain
    addAtom(mol, 1.3f, 2.1f, 0.0f, ATOM_C);        // 13: CH2
    addAtom(mol, 0.0f, 2.8f, 0.0f, ATOM_C);        // 14: COOH
    addAtom(mol, 0.0f, 2.8f, 1.3f, ATOM_O);        // 15: O (carbonyl)
    addAtom(mol, -1.3f, 3.5f, 0.0f, ATOM_O);       // 16: OH
    addAtom(mol, -2.1f, 3.5f, 0.0f, ATOM_H);       // 17: H

    // Guanidinium group
    addAtom(mol, 2.6f, -2.8f, 0.0f, ATOM_N);       // 18: N
    addAtom(mol, 3.9f, -3.5f, 0.0f, ATOM_C);       // 19: C
    addAtom(mol, 5.2f, -2.8f, 0.0f, ATOM_N);       // 20: N
    addAtom(mol, 3.9f, -4.8f, 0.0f, ATOM_N);       // 21: N

    // Ring bonds
    addBond(mol, 0, 1, 1);
    addBond(mol, 1, 2, 1);
    addBond(mol, 2, 3, 2);
    addBond(mol, 3, 4, 1);
    addBond(mol, 4, 5, 1);
    addBond(mol, 5, 0, 1);

    // Amino side chain
    addBond(mol, 0, 6, 1);
    addBond(mol, 6, 7, 1);
    addBond(mol, 6, 8, 1);

    // Ester side chain
    addBond(mol, 2, 9, 1);
    addBond(mol, 9, 10, 1);
    addBond(mol, 10, 11, 2);
    addBond(mol, 10, 12, 1);

    // Carboxylic acid side chain
    addBond(mol, 1, 13, 1);
    addBond(mol, 13, 14, 1);
    addBond(mol, 14, 15, 2);
    addBond(mol, 14, 16, 1);
    addBond(mol, 16, 17, 1);

    // Guanidinium
    addBond(mol, 4, 18, 1);
    addBond(mol, 18, 19, 1);
    addBond(mol, 19, 20, 1);
    addBond(mol, 19, 21, 2);

    centerMolecule(mol);
}

// Build Nirmatrelvir (C23H32F2N2O2) - Paxlovid component
void buildNirmatrelvir(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Nirmatrelvir/Paxlovid (C23H32F2N2O2)");

    // Protease inhibitor with peptide-like core
    // Cyclohexyl ring
    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_C);        // 0
    addAtom(mol, 1.4f, 0.5f, 0.0f, ATOM_C);        // 1
    addAtom(mol, 2.0f, 1.8f, 0.5f, ATOM_C);        // 2
    addAtom(mol, 1.2f, 2.8f, 0.8f, ATOM_C);        // 3
    addAtom(mol, -0.2f, 2.3f, 0.8f, ATOM_C);       // 4
    addAtom(mol, -0.8f, 1.0f, 0.3f, ATOM_C);       // 5

    // Phenyl ring with fluorine
    float r = 1.4f;
    for (int i = 0; i < 6; i++) {
        float angle = i * PI / 3.0f;
        addAtom(mol, r * cosf(angle) - 2.5f, r * sinf(angle), 0.0f, ATOM_C);  // 6-11
    }

    // Fluorine atoms on phenyl
    addAtom(mol, -2.5f + r * cosf(PI/3) + 0.8f, r * sinf(PI/3) + 0.8f, 0.0f, ATOM_F);  // 12
    addAtom(mol, -2.5f + r * cosf(2*PI/3) - 0.8f, r * sinf(2*PI/3) + 0.8f, 0.0f, ATOM_F);  // 13

    // Carboxylic acid ester
    addAtom(mol, 2.6f, -1.2f, 0.0f, ATOM_C);       // 14: C=O
    addAtom(mol, 2.6f, -1.2f, 1.3f, ATOM_O);       // 15: O
    addAtom(mol, 3.8f, -1.9f, 0.0f, ATOM_O);       // 16: O (ester)
    addAtom(mol, 5.1f, -1.3f, 0.0f, ATOM_C);       // 17: CH3

    // Amide group
    addAtom(mol, -0.8f, -1.3f, 0.0f, ATOM_C);      // 18: C=O
    addAtom(mol, -0.8f, -1.3f, 1.3f, ATOM_O);      // 19: O
    addAtom(mol, -2.1f, -2.0f, 0.0f, ATOM_N);      // 20: N
    addAtom(mol, -3.4f, -1.3f, 0.0f, ATOM_C);      // 21: CH3
    addAtom(mol, -2.1f, -3.4f, 0.0f, ATOM_H);      // 22: H

    // Ring bonds
    addBond(mol, 0, 1, 1);
    addBond(mol, 1, 2, 1);
    addBond(mol, 2, 3, 1);
    addBond(mol, 3, 4, 1);
    addBond(mol, 4, 5, 1);
    addBond(mol, 5, 0, 1);

    // Phenyl ring bonds
    for (int i = 0; i < 6; i++) {
        addBond(mol, 6 + i, 6 + ((i + 1) % 6), (i % 2 == 0) ? 2 : 1);
    }

    // Phenyl attachment
    addBond(mol, 5, 6, 1);

    // Fluorine bonds
    addBond(mol, 8, 12, 1);
    addBond(mol, 9, 13, 1);

    // Ester bonds
    addBond(mol, 0, 14, 1);
    addBond(mol, 14, 15, 2);
    addBond(mol, 14, 16, 1);
    addBond(mol, 16, 17, 1);

    // Amide bonds
    addBond(mol, 3, 18, 1);
    addBond(mol, 18, 19, 2);
    addBond(mol, 18, 20, 1);
    addBond(mol, 20, 21, 1);
    addBond(mol, 20, 22, 1);

    centerMolecule(mol);
}

// Build Ritonavir (C37H48N6O5S2) - Paxlovid component/Protease inhibitor
void buildRitonavir(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Ritonavir (C37H48N6O5S2)");

    // Complex peptide-based protease inhibitor - simplified
    // Cyclohexyl ring
    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_C);        // 0
    addAtom(mol, 1.4f, 0.5f, 0.0f, ATOM_C);        // 1
    addAtom(mol, 2.0f, 1.8f, 0.5f, ATOM_C);        // 2
    addAtom(mol, 1.2f, 2.8f, 0.8f, ATOM_C);        // 3
    addAtom(mol, -0.2f, 2.3f, 0.8f, ATOM_C);       // 4
    addAtom(mol, -0.8f, 1.0f, 0.3f, ATOM_C);       // 5

    // Two thiazole rings (simplified)
    addAtom(mol, 2.6f, -1.2f, 0.0f, ATOM_C);       // 6: thiazole1
    addAtom(mol, 3.9f, -1.9f, 0.0f, ATOM_S);       // 7: S
    addAtom(mol, 5.2f, -1.2f, 0.0f, ATOM_C);       // 8
    addAtom(mol, 5.2f, 0.1f, 0.0f, ATOM_N);        // 9: N

    addAtom(mol, -1.3f, -1.3f, 0.0f, ATOM_C);      // 10: thiazole2
    addAtom(mol, -2.6f, -2.0f, 0.0f, ATOM_S);      // 11: S
    addAtom(mol, -3.9f, -1.3f, 0.0f, ATOM_C);      // 12
    addAtom(mol, -3.9f, 0.0f, 0.0f, ATOM_N);       // 13: N

    // Hydroxyl group (HIV protease interaction)
    addAtom(mol, 3.0f, 3.5f, 0.0f, ATOM_O);        // 14: OH
    addAtom(mol, 3.8f, 3.5f, 0.0f, ATOM_H);        // 15: H

    // Two isopropyl groups
    addAtom(mol, 6.5f, -1.9f, 0.0f, ATOM_C);       // 16: CH
    addAtom(mol, 7.2f, -0.5f, 0.0f, ATOM_C);       // 17: CH3
    addAtom(mol, 7.2f, -3.3f, 0.0f, ATOM_C);       // 18: CH3

    addAtom(mol, -5.2f, -2.0f, 0.0f, ATOM_C);      // 19: CH
    addAtom(mol, -5.9f, -0.6f, 0.0f, ATOM_C);      // 20: CH3
    addAtom(mol, -5.9f, -3.4f, 0.0f, ATOM_C);      // 21: CH3

    // Ring bonds
    addBond(mol, 0, 1, 1);
    addBond(mol, 1, 2, 1);
    addBond(mol, 2, 3, 1);
    addBond(mol, 3, 4, 1);
    addBond(mol, 4, 5, 1);
    addBond(mol, 5, 0, 1);

    // Thiazole 1 bonds
    addBond(mol, 0, 6, 1);
    addBond(mol, 6, 7, 1);
    addBond(mol, 7, 8, 1);
    addBond(mol, 8, 9, 1);
    addBond(mol, 9, 6, 1);

    // Thiazole 2 bonds
    addBond(mol, 5, 10, 1);
    addBond(mol, 10, 11, 1);
    addBond(mol, 11, 12, 1);
    addBond(mol, 12, 13, 1);
    addBond(mol, 13, 10, 1);

    // Hydroxyl
    addBond(mol, 3, 14, 1);
    addBond(mol, 14, 15, 1);

    // Isopropyl groups
    addBond(mol, 9, 16, 1);
    addBond(mol, 16, 17, 1);
    addBond(mol, 16, 18, 1);

    addBond(mol, 13, 19, 1);
    addBond(mol, 19, 20, 1);
    addBond(mol, 19, 21, 1);

    centerMolecule(mol);
}

// Build Molnupiravir/Lagevrio (C13H19N3O8) - Antiviral
void buildMolnupiravir(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Molnupiravir/Lagevrio (C13H19N3O8)");

    // Ribose-like sugar with modified cytosine
    // Pentose ring (furanose)
    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_C);        // 0: C1
    addAtom(mol, 1.3f, 0.7f, 0.0f, ATOM_O);        // 1: O (in ring)
    addAtom(mol, 0.6f, 2.0f, 0.0f, ATOM_C);        // 2: C2
    addAtom(mol, -0.8f, 1.8f, 0.0f, ATOM_C);       // 3: C3
    addAtom(mol, -0.8f, 0.4f, 0.0f, ATOM_C);       // 4: C4

    // Hydroxyl groups on sugar
    addAtom(mol, 2.0f, -0.5f, 0.0f, ATOM_O);       // 5: OH (C1)
    addAtom(mol, 2.8f, -0.5f, 0.0f, ATOM_H);       // 6: H
    addAtom(mol, 1.8f, 3.2f, 0.0f, ATOM_O);        // 7: OH (C2)
    addAtom(mol, 2.6f, 3.2f, 0.0f, ATOM_H);        // 8: H
    addAtom(mol, -1.6f, 2.8f, 0.0f, ATOM_O);       // 9: OH (C3)
    addAtom(mol, -2.4f, 2.8f, 0.0f, ATOM_H);       // 10: H

    // Phosphate group
    addAtom(mol, -1.5f, -0.8f, 0.0f, ATOM_P);      // 11: P
    addAtom(mol, -1.5f, -0.8f, 1.3f, ATOM_O);      // 12: O (double bond)
    addAtom(mol, -2.8f, -1.5f, 0.0f, ATOM_O);      // 13: O (single)
    addAtom(mol, -0.2f, -1.5f, 0.0f, ATOM_O);      // 14: O (single)

    // Modified cytosine base
    addAtom(mol, 0.6f, -1.4f, 0.0f, ATOM_C);       // 15: C (aromatic)
    addAtom(mol, 1.9f, -2.1f, 0.0f, ATOM_N);       // 16: N
    addAtom(mol, 3.2f, -1.4f, 0.0f, ATOM_C);       // 17: C (aromatic)
    addAtom(mol, 3.2f, -0.1f, 0.0f, ATOM_C);       // 18: C=O
    addAtom(mol, 3.2f, -0.1f, 1.3f, ATOM_O);       // 19: O
    addAtom(mol, 1.9f, 0.6f, 0.0f, ATOM_N);        // 20: N

    // Isopropylidene group (modification)
    addAtom(mol, 4.5f, -2.1f, 0.0f, ATOM_C);       // 21: C (quaternary)
    addAtom(mol, 5.2f, -0.7f, 0.0f, ATOM_C);       // 22: CH3
    addAtom(mol, 5.2f, -3.5f, 0.0f, ATOM_C);       // 23: CH3

    // Ring bonds
    addBond(mol, 0, 1, 1);
    addBond(mol, 1, 2, 1);
    addBond(mol, 2, 3, 1);
    addBond(mol, 3, 4, 1);
    addBond(mol, 4, 0, 1);

    // Hydroxyl bonds
    addBond(mol, 0, 5, 1);
    addBond(mol, 5, 6, 1);
    addBond(mol, 2, 7, 1);
    addBond(mol, 7, 8, 1);
    addBond(mol, 3, 9, 1);
    addBond(mol, 9, 10, 1);

    // Phosphate bonds
    addBond(mol, 4, 11, 1);
    addBond(mol, 11, 12, 2);
    addBond(mol, 11, 13, 1);
    addBond(mol, 11, 14, 1);

    // Cytosine base bonds
    addBond(mol, 0, 15, 1);
    addBond(mol, 15, 16, 1);
    addBond(mol, 16, 17, 1);
    addBond(mol, 17, 18, 1);
    addBond(mol, 18, 19, 2);
    addBond(mol, 18, 20, 1);
    addBond(mol, 20, 15, 1);

    // Isopropylidene bonds
    addBond(mol, 17, 21, 1);
    addBond(mol, 21, 22, 1);
    addBond(mol, 21, 23, 1);

    centerMolecule(mol);
}

// Build Remdesivir/Veklury (C27H35N6O8P) - Antiviral
void buildRemdesivir(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Remdesivir/Veklury (C27H35N6O8P)");

    // Adenosine analog with phosphoramide prodrug
    // Pentose ring
    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_C);        // 0: C1
    addAtom(mol, 1.3f, 0.7f, 0.0f, ATOM_O);        // 1: O (in ring)
    addAtom(mol, 0.6f, 2.0f, 0.0f, ATOM_C);        // 2: C2
    addAtom(mol, -0.8f, 1.8f, 0.0f, ATOM_C);       // 3: C3
    addAtom(mol, -0.8f, 0.4f, 0.0f, ATOM_C);       // 4: C4

    // Hydroxyl groups
    addAtom(mol, 2.0f, -0.5f, 0.0f, ATOM_O);       // 5: OH
    addAtom(mol, 2.8f, -0.5f, 0.0f, ATOM_H);       // 6: H
    addAtom(mol, 1.8f, 3.2f, 0.0f, ATOM_O);        // 7: OH
    addAtom(mol, 2.6f, 3.2f, 0.0f, ATOM_H);        // 8: H
    addAtom(mol, -1.6f, 2.8f, 0.0f, ATOM_O);       // 9: OH
    addAtom(mol, -2.4f, 2.8f, 0.0f, ATOM_H);       // 10: H

    // Purine base (adenine-like)
    addAtom(mol, 0.6f, -1.4f, 0.0f, ATOM_C);       // 11: C (aromatic)
    addAtom(mol, 1.9f, -2.1f, 0.0f, ATOM_N);       // 12: N
    addAtom(mol, 3.2f, -1.4f, 0.0f, ATOM_C);       // 13: C (aromatic)
    addAtom(mol, 3.2f, -0.1f, 0.0f, ATOM_N);       // 14: N
    addAtom(mol, 1.9f, 0.6f, 0.0f, ATOM_N);        // 15: N

    // Six-membered ring of purine
    addAtom(mol, 4.5f, -2.1f, 0.0f, ATOM_C);       // 16
    addAtom(mol, 4.5f, -3.4f, 0.0f, ATOM_N);       // 17: N
    addAtom(mol, 3.2f, -4.1f, 0.0f, ATOM_C);       // 18: C
    addAtom(mol, 1.9f, -3.4f, 0.0f, ATOM_N);       // 19: N

    // Phosphoramide group
    addAtom(mol, -1.5f, -0.8f, 0.0f, ATOM_P);      // 20: P
    addAtom(mol, -1.5f, -0.8f, 1.3f, ATOM_O);      // 21: O (double bond)
    addAtom(mol, -2.8f, -1.5f, 0.0f, ATOM_O);      // 22: O
    addAtom(mol, -0.2f, -1.5f, 0.0f, ATOM_N);      // 23: N

    // Isopropyl group on phosphoramide
    addAtom(mol, -3.6f, -1.5f, 0.0f, ATOM_C);      // 24: CH
    addAtom(mol, -4.3f, -0.1f, 0.0f, ATOM_C);      // 25: CH3
    addAtom(mol, -4.3f, -2.9f, 0.0f, ATOM_C);      // 26: CH3

    // Aliphatic side chain
    addAtom(mol, -0.9f, -2.8f, 0.0f, ATOM_C);      // 27: CH3

    // Ring bonds
    addBond(mol, 0, 1, 1);
    addBond(mol, 1, 2, 1);
    addBond(mol, 2, 3, 1);
    addBond(mol, 3, 4, 1);
    addBond(mol, 4, 0, 1);

    // Hydroxyl bonds
    addBond(mol, 0, 5, 1);
    addBond(mol, 5, 6, 1);
    addBond(mol, 2, 7, 1);
    addBond(mol, 7, 8, 1);
    addBond(mol, 3, 9, 1);
    addBond(mol, 9, 10, 1);

    // Base bonds
    addBond(mol, 0, 11, 1);
    addBond(mol, 11, 12, 1);
    addBond(mol, 12, 13, 1);
    addBond(mol, 13, 14, 1);
    addBond(mol, 14, 15, 1);
    addBond(mol, 15, 11, 1);
    addBond(mol, 13, 16, 1);
    addBond(mol, 16, 17, 1);
    addBond(mol, 17, 18, 1);
    addBond(mol, 18, 19, 1);
    addBond(mol, 19, 12, 1);

    // Phosphoramide bonds
    addBond(mol, 4, 20, 1);
    addBond(mol, 20, 21, 2);
    addBond(mol, 20, 22, 1);
    addBond(mol, 20, 23, 1);
    addBond(mol, 22, 24, 1);
    addBond(mol, 24, 25, 1);
    addBond(mol, 24, 26, 1);
    addBond(mol, 23, 27, 1);

    centerMolecule(mol);
}

// Build Dexamethasone (C22H29FO5) - Corticosteroid
void buildDexamethasone(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Dexamethasone (C22H29FO5)");

    // Steroid nucleus with four fused rings
    // Ring A
    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_C);        // 0
    addAtom(mol, 1.4f, 0.5f, 0.0f, ATOM_C);        // 1
    addAtom(mol, 2.1f, 1.8f, 0.0f, ATOM_C);        // 2
    addAtom(mol, 0.7f, 2.3f, 0.0f, ATOM_C);        // 3
    addAtom(mol, -0.7f, 1.8f, 0.0f, ATOM_C);       // 4
    addAtom(mol, -0.7f, 0.4f, 0.0f, ATOM_C);       // 5

    // Ring B (fused)
    addAtom(mol, 1.4f, 2.0f, 0.0f, ATOM_C);        // 6: shared C
    addAtom(mol, 2.8f, 2.5f, 0.0f, ATOM_C);        // 7
    addAtom(mol, 2.8f, 3.9f, 0.0f, ATOM_C);        // 8
    addAtom(mol, 1.4f, 4.4f, 0.0f, ATOM_C);        // 9

    // Ring C (fused)
    addAtom(mol, 0.0f, 3.7f, 0.0f, ATOM_C);        // 10: shared C
    addAtom(mol, -1.4f, 4.2f, 0.0f, ATOM_C);       // 11
    addAtom(mol, -2.1f, 2.9f, 0.0f, ATOM_C);       // 12

    // Ring D (fused)
    addAtom(mol, -2.1f, 1.5f, 0.0f, ATOM_C);       // 13: shared C
    addAtom(mol, -3.5f, 1.0f, 0.0f, ATOM_C);       // 14
    addAtom(mol, -4.2f, 2.3f, 0.0f, ATOM_C);       // 15

    // Ketone at C3
    addAtom(mol, 0.7f, 3.6f, 0.0f, ATOM_C);        // 16: C=O
    addAtom(mol, 0.7f, 3.6f, 1.3f, ATOM_O);        // 17: O

    // Hydroxyl at C17
    addAtom(mol, 0.7f, 5.7f, 0.0f, ATOM_O);        // 18: OH
    addAtom(mol, 1.5f, 5.7f, 0.0f, ATOM_H);        // 19: H

    // Fluorine at C9
    addAtom(mol, 0.0f, 5.7f, 0.0f, ATOM_F);        // 20: F

    // Hydroxyl at C11
    addAtom(mol, -2.1f, 5.6f, 0.0f, ATOM_O);       // 21: OH
    addAtom(mol, -2.9f, 5.6f, 0.0f, ATOM_H);       // 22: H

    // Methyl groups (angular methyls)
    addAtom(mol, -0.7f, -1.0f, 0.0f, ATOM_C);      // 23: C18 methyl
    addAtom(mol, 2.8f, 0.0f, 0.0f, ATOM_C);        // 24: C19 methyl

    // Ethyl side chain at C17
    addAtom(mol, 2.1f, 5.2f, 0.0f, ATOM_C);        // 25: CH2
    addAtom(mol, 3.4f, 5.9f, 0.0f, ATOM_C);        // 26: CH3

    // Ring A bonds
    addBond(mol, 0, 1, 1);
    addBond(mol, 1, 2, 1);
    addBond(mol, 2, 3, 1);
    addBond(mol, 3, 4, 1);
    addBond(mol, 4, 5, 1);
    addBond(mol, 5, 0, 1);

    // Ring B bonds
    addBond(mol, 6, 7, 1);
    addBond(mol, 7, 8, 1);
    addBond(mol, 8, 9, 1);
    addBond(mol, 9, 10, 1);

    // Ring C bonds
    addBond(mol, 10, 11, 1);
    addBond(mol, 11, 12, 1);
    addBond(mol, 12, 13, 1);

    // Ring D bonds
    addBond(mol, 13, 14, 1);
    addBond(mol, 14, 15, 1);
    addBond(mol, 15, 12, 1);

    // Ketone
    addBond(mol, 3, 16, 1);
    addBond(mol, 16, 17, 2);

    // Hydroxyl at C17
    addBond(mol, 9, 18, 1);
    addBond(mol, 18, 19, 1);

    // Fluorine at C9
    addBond(mol, 9, 20, 1);

    // Hydroxyl at C11
    addBond(mol, 11, 21, 1);
    addBond(mol, 21, 22, 1);

    // Angular methyls
    addBond(mol, 0, 23, 1);
    addBond(mol, 5, 24, 1);

    // Ethyl side chain
    addBond(mol, 9, 25, 1);
    addBond(mol, 25, 26, 1);

    centerMolecule(mol);
}

// Fluconazole (C13H12FN3O) - Triazole antifungal
void buildFluconazole(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Fluconazole (C13H12FN3O)");

    // Imidazole ring (N-containing)
    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_N);        // 0: N1
    addAtom(mol, 1.2f, 0.5f, 0.0f, ATOM_C);        // 1: C
    addAtom(mol, 0.8f, 1.7f, 0.0f, ATOM_N);        // 2: N
    addAtom(mol, -0.5f, 1.5f, 0.0f, ATOM_N);       // 3: N

    // Triazole linker to phenyl
    addAtom(mol, -1.0f, 2.8f, 0.0f, ATOM_C);       // 4: Linker C

    // Phenyl ring
    float r = 1.4f;
    addAtom(mol, -2.5f, 3.1f, 0.0f, ATOM_C);       // 5: C1
    addAtom(mol, -3.2f, 4.3f, 0.0f, ATOM_C);       // 6: C2
    addAtom(mol, -2.5f, 5.5f, 0.0f, ATOM_C);       // 7: C3 (para to linker)
    addAtom(mol, -1.1f, 5.5f, 0.0f, ATOM_C);       // 8: C4
    addAtom(mol, -0.4f, 4.3f, 0.0f, ATOM_C);       // 9: C5
    addAtom(mol, -1.1f, 3.1f, 0.0f, ATOM_C);       // 10: C6 (connected to linker)

    // Fluorine on phenyl
    addAtom(mol, -3.2f, 6.7f, 0.0f, ATOM_F);       // 11: F

    // Alcohol side chain
    addAtom(mol, 2.5f, 0.3f, 0.0f, ATOM_C);        // 12: OCH2CH2OH linker
    addAtom(mol, 3.5f, 1.2f, 0.0f, ATOM_C);        // 13: CH2
    addAtom(mol, 4.8f, 0.8f, 0.0f, ATOM_C);        // 14: CH2
    addAtom(mol, 6.0f, 1.7f, 0.0f, ATOM_O);        // 15: OH
    addAtom(mol, 6.8f, 1.0f, 0.0f, ATOM_H);        // 16: H

    // Hydrogens on imidazole and rings
    addAtom(mol, 1.9f, -0.3f, 0.0f, ATOM_H);       // 17

    // Imidazole ring bonds
    addBond(mol, 0, 1, 1);
    addBond(mol, 1, 2, 2);
    addBond(mol, 2, 3, 1);
    addBond(mol, 3, 0, 1);

    // Linker from triazole to phenyl
    addBond(mol, 3, 4, 1);
    addBond(mol, 4, 5, 1);

    // Phenyl ring bonds
    addBond(mol, 5, 6, 1);
    addBond(mol, 6, 7, 2);
    addBond(mol, 7, 8, 1);
    addBond(mol, 8, 9, 2);
    addBond(mol, 9, 10, 1);
    addBond(mol, 10, 4, 2);

    // Fluorine bond
    addBond(mol, 7, 11, 1);

    // Alcohol side chain
    addBond(mol, 1, 12, 1);
    addBond(mol, 12, 13, 1);
    addBond(mol, 13, 14, 1);
    addBond(mol, 14, 15, 1);
    addBond(mol, 15, 16, 1);

    // H on imidazole
    addBond(mol, 1, 17, 1);

    centerMolecule(mol);
}

// Methotrexate (C20H22N8O5) - Simplified representation
void buildMethotrexate(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Methotrexate (C20H22N8O5)");

    // Pteridine ring (fused pyrimidine-pyrazine)
    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_N);        // 0
    addAtom(mol, 1.2f, 0.5f, 0.0f, ATOM_C);        // 1
    addAtom(mol, 1.2f, 1.8f, 0.0f, ATOM_N);        // 2
    addAtom(mol, 0.0f, 2.3f, 0.0f, ATOM_C);        // 3
    addAtom(mol, -1.2f, 1.8f, 0.0f, ATOM_N);       // 4
    addAtom(mol, -1.2f, 0.5f, 0.0f, ATOM_C);       // 5
    addAtom(mol, 2.4f, 0.0f, 0.0f, ATOM_C);        // 6
    addAtom(mol, 2.4f, -1.3f, 0.0f, ATOM_N);       // 7
    addAtom(mol, 0.0f, -1.3f, 0.0f, ATOM_N);       // 8

    // p-Aminobenzoyl glutamate
    addAtom(mol, 0.0f, 3.5f, 0.0f, ATOM_C);        // 9: Linker C
    addAtom(mol, 1.5f, 4.2f, 0.0f, ATOM_C);        // 10: Benzene
    addAtom(mol, 2.9f, 3.8f, 0.0f, ATOM_C);        // 11
    addAtom(mol, 3.9f, 4.8f, 0.0f, ATOM_C);        // 12
    addAtom(mol, 3.4f, 6.1f, 0.0f, ATOM_C);        // 13
    addAtom(mol, 2.0f, 6.5f, 0.0f, ATOM_C);        // 14
    addAtom(mol, 1.0f, 5.5f, 0.0f, ATOM_C);        // 15

    // Amino group on benzene
    addAtom(mol, 4.8f, 4.4f, 0.0f, ATOM_N);        // 16: NH2
    addAtom(mol, 5.6f, 5.1f, 0.0f, ATOM_H);        // 17
    addAtom(mol, 5.6f, 3.6f, 0.0f, ATOM_H);        // 18

    // Glutamate tail
    addAtom(mol, 4.5f, 7.2f, 0.0f, ATOM_C);        // 19: C=O
    addAtom(mol, 5.8f, 6.8f, 0.0f, ATOM_O);        // 20: O
    addAtom(mol, -1.2f, 4.1f, 0.0f, ATOM_C);       // 21: Glu chain
    addAtom(mol, -2.4f, 3.5f, 0.0f, ATOM_C);       // 22
    addAtom(mol, -3.6f, 4.3f, 0.0f, ATOM_C);       // 23
    addAtom(mol, -4.8f, 3.6f, 0.0f, ATOM_C);       // 24
    addAtom(mol, -6.0f, 4.4f, 0.0f, ATOM_O);       // 25: COOH
    addAtom(mol, -6.8f, 3.7f, 0.0f, ATOM_O);       // 26: OH

    // Pteridine ring bonds
    addBond(mol, 0, 1, 1);
    addBond(mol, 1, 2, 2);
    addBond(mol, 2, 3, 1);
    addBond(mol, 3, 4, 2);
    addBond(mol, 4, 5, 1);
    addBond(mol, 5, 0, 1);
    addBond(mol, 1, 6, 1);
    addBond(mol, 6, 7, 2);
    addBond(mol, 7, 8, 1);
    addBond(mol, 8, 0, 1);

    // Para-aminobenzoyl linker
    addBond(mol, 3, 9, 1);
    addBond(mol, 9, 10, 1);
    addBond(mol, 10, 15, 1);
    addBond(mol, 15, 1, 1);
    addBond(mol, 12, 13, 1);
    addBond(mol, 13, 14, 1);
    addBond(mol, 14, 15, 1);
    addBond(mol, 10, 11, 1);
    addBond(mol, 11, 12, 1);

    // Amino group
    addBond(mol, 12, 16, 1);
    addBond(mol, 16, 17, 1);
    addBond(mol, 16, 18, 1);

    // Glutamate tail
    addBond(mol, 13, 19, 1);
    addBond(mol, 19, 20, 2);
    addBond(mol, 9, 21, 1);
    addBond(mol, 21, 22, 1);
    addBond(mol, 22, 23, 1);
    addBond(mol, 23, 24, 1);
    addBond(mol, 24, 25, 1);
    addBond(mol, 25, 26, 1);

    centerMolecule(mol);
}

// Hydroxychloroquine (C18H26ClN3O) - Simplified representation
void buildHydroxychloroquine(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Hydroxychloroquine (C18H26ClN3O)");

    // Quinoline ring (main therapeutic moiety)
    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_N);        // 0
    addAtom(mol, 1.2f, 0.5f, 0.0f, ATOM_C);        // 1
    addAtom(mol, 2.2f, 0.0f, 0.0f, ATOM_C);        // 2
    addAtom(mol, 2.2f, -1.3f, 0.0f, ATOM_C);       // 3
    addAtom(mol, 1.2f, -1.8f, 0.0f, ATOM_C);       // 4
    addAtom(mol, 0.2f, -1.3f, 0.0f, ATOM_C);       // 5
    addAtom(mol, 0.2f, 0.0f, 0.0f, ATOM_C);        // 6

    // Chlorine substituent on quinoline
    addAtom(mol, 3.4f, 0.7f, 0.0f, ATOM_CL);       // 7: Cl

    // Side chain ether
    addAtom(mol, 1.2f, 1.9f, 0.0f, ATOM_C);        // 8: -OCH2-
    addAtom(mol, 1.2f, 1.2f, 0.0f, ATOM_O);        // 9: O
    addAtom(mol, 2.4f, 2.4f, 0.0f, ATOM_C);        // 10: CH2

    // Aliphatic chain with diethylamino
    addAtom(mol, 3.6f, 1.9f, 0.0f, ATOM_C);        // 11: CH2
    addAtom(mol, 4.8f, 2.4f, 0.0f, ATOM_C);        // 12: CH2
    addAtom(mol, 6.0f, 1.9f, 0.0f, ATOM_N);        // 13: N
    addAtom(mol, 7.2f, 2.4f, 0.0f, ATOM_C);        // 14: CH2 (ethyl 1)
    addAtom(mol, 8.4f, 1.9f, 0.0f, ATOM_C);        // 15: CH3 (ethyl 1)
    addAtom(mol, 6.0f, 3.2f, 0.0f, ATOM_C);        // 16: CH2 (ethyl 2)
    addAtom(mol, 7.2f, 3.7f, 0.0f, ATOM_C);        // 17: CH3 (ethyl 2)

    // Hydroxyl group
    addAtom(mol, 1.2f, 3.2f, 0.0f, ATOM_O);        // 18: OH
    addAtom(mol, 0.5f, 3.8f, 0.0f, ATOM_H);        // 19: H

    // Quinoline ring bonds
    addBond(mol, 0, 1, 2);
    addBond(mol, 1, 2, 1);
    addBond(mol, 2, 3, 2);
    addBond(mol, 3, 4, 1);
    addBond(mol, 4, 5, 2);
    addBond(mol, 5, 6, 1);
    addBond(mol, 6, 0, 2);

    // Chlorine
    addBond(mol, 2, 7, 1);

    // Side chain ether to diethylamino
    addBond(mol, 1, 8, 1);
    addBond(mol, 8, 9, 1);
    addBond(mol, 9, 10, 1);
    addBond(mol, 10, 11, 1);
    addBond(mol, 11, 12, 1);
    addBond(mol, 12, 13, 1);
    addBond(mol, 13, 14, 1);
    addBond(mol, 14, 15, 1);
    addBond(mol, 13, 16, 1);
    addBond(mol, 16, 17, 1);

    // Hydroxyl
    addBond(mol, 8, 18, 1);
    addBond(mol, 18, 19, 1);

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

// Build Adamantane (C10H16) - diamond-like cage structure
void buildAdamantane(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Adamantane (C10H16)");

    // Adamantane has a diamond cubic lattice structure
    // 4 "bridgehead" carbons (CH) at tetrahedral positions
    // 6 "bridge" carbons (CH2) between them

    float s = 0.9f;  // Scale factor

    // 4 bridgehead carbons (CH) at tetrahedral positions
    addAtom(mol,  s,  s,  s, ATOM_C);  // 0
    addAtom(mol,  s, -s, -s, ATOM_C);  // 1
    addAtom(mol, -s,  s, -s, ATOM_C);  // 2
    addAtom(mol, -s, -s,  s, ATOM_C);  // 3

    // 6 bridge carbons (CH2) at edge midpoints
    addAtom(mol,  0,  0,  s*1.5f, ATOM_C);  // 4 - between 0,3
    addAtom(mol,  0,  s*1.5f,  0, ATOM_C);  // 5 - between 0,2
    addAtom(mol,  s*1.5f,  0,  0, ATOM_C);  // 6 - between 0,1
    addAtom(mol,  0,  0, -s*1.5f, ATOM_C);  // 7 - between 1,2
    addAtom(mol,  0, -s*1.5f,  0, ATOM_C);  // 8 - between 1,3
    addAtom(mol, -s*1.5f,  0,  0, ATOM_C);  // 9 - between 2,3

    // Hydrogens on bridgehead carbons (1 each, pointing outward)
    float h = 1.8f;
    addAtom(mol,  h,  h,  h, ATOM_H);  // 10 on C0
    addAtom(mol,  h, -h, -h, ATOM_H);  // 11 on C1
    addAtom(mol, -h,  h, -h, ATOM_H);  // 12 on C2
    addAtom(mol, -h, -h,  h, ATOM_H);  // 13 on C3

    // Hydrogens on bridge carbons (2 each)
    float hb = 0.9f;
    // C4 (0,0,s*1.5)
    addAtom(mol,  hb,  hb, s*1.5f + hb*0.5f, ATOM_H);  // 14
    addAtom(mol, -hb, -hb, s*1.5f + hb*0.5f, ATOM_H);  // 15
    // C5 (0,s*1.5,0)
    addAtom(mol,  hb, s*1.5f + hb*0.5f,  hb, ATOM_H);  // 16
    addAtom(mol, -hb, s*1.5f + hb*0.5f, -hb, ATOM_H);  // 17
    // C6 (s*1.5,0,0)
    addAtom(mol, s*1.5f + hb*0.5f,  hb,  hb, ATOM_H);  // 18
    addAtom(mol, s*1.5f + hb*0.5f, -hb, -hb, ATOM_H);  // 19
    // C7 (0,0,-s*1.5)
    addAtom(mol,  hb,  hb, -s*1.5f - hb*0.5f, ATOM_H);  // 20
    addAtom(mol, -hb, -hb, -s*1.5f - hb*0.5f, ATOM_H);  // 21
    // C8 (0,-s*1.5,0)
    addAtom(mol,  hb, -s*1.5f - hb*0.5f,  hb, ATOM_H);  // 22
    addAtom(mol, -hb, -s*1.5f - hb*0.5f, -hb, ATOM_H);  // 23
    // C9 (-s*1.5,0,0)
    addAtom(mol, -s*1.5f - hb*0.5f,  hb,  hb, ATOM_H);  // 24
    addAtom(mol, -s*1.5f - hb*0.5f, -hb, -hb, ATOM_H);  // 25

    // C-C bonds: bridgehead to bridge carbons (12 bonds)
    // C0 connects to C4, C5, C6
    addBond(mol, 0, 4, 1);
    addBond(mol, 0, 5, 1);
    addBond(mol, 0, 6, 1);
    // C1 connects to C6, C7, C8
    addBond(mol, 1, 6, 1);
    addBond(mol, 1, 7, 1);
    addBond(mol, 1, 8, 1);
    // C2 connects to C5, C7, C9
    addBond(mol, 2, 5, 1);
    addBond(mol, 2, 7, 1);
    addBond(mol, 2, 9, 1);
    // C3 connects to C4, C8, C9
    addBond(mol, 3, 4, 1);
    addBond(mol, 3, 8, 1);
    addBond(mol, 3, 9, 1);

    // C-H bonds
    addBond(mol, 0, 10, 1);
    addBond(mol, 1, 11, 1);
    addBond(mol, 2, 12, 1);
    addBond(mol, 3, 13, 1);
    addBond(mol, 4, 14, 1);
    addBond(mol, 4, 15, 1);
    addBond(mol, 5, 16, 1);
    addBond(mol, 5, 17, 1);
    addBond(mol, 6, 18, 1);
    addBond(mol, 6, 19, 1);
    addBond(mol, 7, 20, 1);
    addBond(mol, 7, 21, 1);
    addBond(mol, 8, 22, 1);
    addBond(mol, 8, 23, 1);
    addBond(mol, 9, 24, 1);
    addBond(mol, 9, 25, 1);

    centerMolecule(mol);
}

// Build Chlorophyll a (C55H72MgN4O5) - simplified porphyrin core
void buildChlorophyll(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Chlorophyll a (core)");

    // Chlorophyll has a porphyrin ring with Mg at center
    // We'll build the core porphyrin macrocycle with Mg
    // Full chlorophyll has a long phytol tail which we'll simplify

    // Central Magnesium atom (using Fe slot since we don't have Mg)
    // Note: Mg would be index 18+ but we'll use a placeholder
    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_FE);  // 0 - Mg center (shown as Fe color)

    // Porphyrin ring: 4 pyrrole rings connected by methine bridges
    float r1 = 2.0f;   // Inner ring radius (N atoms)
    float r2 = 2.8f;   // Outer ring radius (C atoms)
    float r3 = 3.6f;   // Methine bridge radius

    // 4 Nitrogen atoms coordinating to Mg (pyrrole nitrogens)
    for (int i = 0; i < 4; i++) {
        float angle = i * PI / 2.0f;
        addAtom(mol, r1 * cosf(angle), r1 * sinf(angle), 0.0f, ATOM_N);
    }
    // N atoms are 1, 2, 3, 4

    // Build each pyrrole ring (4 pyrroles, each with 5 atoms including N)
    // Each pyrrole: N already placed, add 4 carbons around it
    for (int p = 0; p < 4; p++) {
        float baseAngle = p * PI / 2.0f;

        // Alpha carbons (adjacent to N, inner)
        float a1 = baseAngle - 0.35f;
        float a2 = baseAngle + 0.35f;
        addAtom(mol, r2 * cosf(a1), r2 * sinf(a1), 0.0f, ATOM_C);
        addAtom(mol, r2 * cosf(a2), r2 * sinf(a2), 0.0f, ATOM_C);

        // Beta carbons (outer edge of pyrrole)
        float b1 = baseAngle - 0.25f;
        float b2 = baseAngle + 0.25f;
        float rb = 3.5f;
        addAtom(mol, rb * cosf(b1), rb * sinf(b1), 0.0f, ATOM_C);
        addAtom(mol, rb * cosf(b2), rb * sinf(b2), 0.0f, ATOM_C);
    }
    // Pyrrole carbons: 5-20 (4 per pyrrole x 4 pyrroles)

    // Methine bridges between pyrroles (4 CH bridges)
    for (int i = 0; i < 4; i++) {
        float angle = i * PI / 2.0f + PI / 4.0f;  // 45 degrees offset
        addAtom(mol, r3 * cosf(angle), r3 * sinf(angle), 0.0f, ATOM_C);
    }
    // Methine carbons: 21, 22, 23, 24

    // Add some methyl groups on the pyrroles (characteristic of chlorophyll)
    float rm = 4.5f;
    for (int i = 0; i < 4; i++) {
        float angle = i * PI / 2.0f;
        addAtom(mol, rm * cosf(angle), rm * sinf(angle), 0.3f, ATOM_C);
    }
    // Methyl carbons: 25, 26, 27, 28

    // Add hydrogens on methine bridges
    for (int i = 0; i < 4; i++) {
        float angle = i * PI / 2.0f + PI / 4.0f;
        addAtom(mol, (r3 + 0.9f) * cosf(angle), (r3 + 0.9f) * sinf(angle), 0.0f, ATOM_H);
    }
    // Hydrogens: 29, 30, 31, 32

    // Simplified phytol tail (just a few carbons to suggest the chain)
    addAtom(mol, 4.0f, -2.0f, 0.0f, ATOM_C);  // 33
    addAtom(mol, 5.2f, -2.5f, 0.0f, ATOM_C);  // 34
    addAtom(mol, 6.4f, -2.0f, 0.0f, ATOM_C);  // 35
    addAtom(mol, 7.6f, -2.5f, 0.0f, ATOM_C);  // 36
    addAtom(mol, 8.8f, -2.0f, 0.0f, ATOM_C);  // 37

    // Carbonyl oxygen (part of the isocyclic ring in real chlorophyll)
    addAtom(mol, 3.0f, -3.5f, 0.0f, ATOM_O);  // 38

    // Bonds: Mg to nitrogens
    addBond(mol, 0, 1, 1);
    addBond(mol, 0, 2, 1);
    addBond(mol, 0, 3, 1);
    addBond(mol, 0, 4, 1);

    // Pyrrole ring bonds
    for (int p = 0; p < 4; p++) {
        int n = 1 + p;           // Nitrogen index
        int c1 = 5 + p * 4;      // First alpha carbon
        int c2 = 5 + p * 4 + 1;  // Second alpha carbon
        int c3 = 5 + p * 4 + 2;  // First beta carbon
        int c4 = 5 + p * 4 + 3;  // Second beta carbon

        addBond(mol, n, c1, 1);
        addBond(mol, n, c2, 1);
        addBond(mol, c1, c3, 2);
        addBond(mol, c2, c4, 2);
        addBond(mol, c3, c4, 1);
    }

    // Methine bridge bonds (connecting pyrroles)
    for (int i = 0; i < 4; i++) {
        int methine = 21 + i;
        int pyrrole1 = i;
        int pyrrole2 = (i + 1) % 4;
        int c1 = 5 + pyrrole1 * 4 + 1;  // Alpha carbon of current pyrrole
        int c2 = 5 + pyrrole2 * 4;      // Alpha carbon of next pyrrole
        addBond(mol, methine, c1, 1);
        addBond(mol, methine, c2, 1);
    }

    // Methine to hydrogen bonds
    for (int i = 0; i < 4; i++) {
        addBond(mol, 21 + i, 29 + i, 1);
    }

    // Methyl group bonds
    for (int i = 0; i < 4; i++) {
        int beta = 5 + i * 4 + 2;  // Beta carbon
        addBond(mol, 25 + i, beta, 1);
    }

    // Phytol tail bonds
    addBond(mol, 5 + 3 * 4 + 3, 33, 1);  // Connect to last pyrrole
    addBond(mol, 33, 34, 1);
    addBond(mol, 34, 35, 1);
    addBond(mol, 35, 36, 1);
    addBond(mol, 36, 37, 1);

    // Carbonyl bond
    addBond(mol, 33, 38, 2);

    centerMolecule(mol);
}

// Build Xenon Difluoride (XeF2) - linear noble gas compound
void buildXenonDifluoride(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Xenon Difluoride (XeF2)");

    // Linear geometry: F-Xe-F (180 degrees)
    // Xe-F bond length ~2.0 Angstroms
    float bondLen = 1.5f;

    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_PT);      // 0 - Xe (using Pt color - silver)
    addAtom(mol, -bondLen, 0.0f, 0.0f, ATOM_F);   // 1 - F
    addAtom(mol, bondLen, 0.0f, 0.0f, ATOM_F);    // 2 - F

    addBond(mol, 0, 1, 1);
    addBond(mol, 0, 2, 1);

    centerMolecule(mol);
}

// Build Rhenium(III) Chloride Dimer ([Re2Cl8]2-) - metal-metal quadruple bond
void buildRheniumChlorideDimer(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "[Re2Cl8]^2- Quadruple Bond");

    // Two Re atoms with quadruple bond between them
    // Each Re has 4 Cl ligands in eclipsed square planar arrangement
    float reBond = 1.0f;   // Re-Re bond length (very short for quadruple bond)
    float reCl = 1.8f;     // Re-Cl bond length

    // Rhenium atoms
    addAtom(mol, -reBond/2, 0.0f, 0.0f, ATOM_RE);  // 0 - Re1
    addAtom(mol,  reBond/2, 0.0f, 0.0f, ATOM_RE);  // 1 - Re2

    // Chlorines on Re1 (eclipsed square planar)
    addAtom(mol, -reBond/2,  reCl,  0.0f, ATOM_CL);  // 2
    addAtom(mol, -reBond/2, -reCl,  0.0f, ATOM_CL);  // 3
    addAtom(mol, -reBond/2,  0.0f,  reCl, ATOM_CL);  // 4
    addAtom(mol, -reBond/2,  0.0f, -reCl, ATOM_CL);  // 5

    // Chlorines on Re2 (eclipsed with Re1's ligands)
    addAtom(mol,  reBond/2,  reCl,  0.0f, ATOM_CL);  // 6
    addAtom(mol,  reBond/2, -reCl,  0.0f, ATOM_CL);  // 7
    addAtom(mol,  reBond/2,  0.0f,  reCl, ATOM_CL);  // 8
    addAtom(mol,  reBond/2,  0.0f, -reCl, ATOM_CL);  // 9

    // Re-Re quadruple bond (shown as order 3 since we max at triple visually)
    addBond(mol, 0, 1, 3);

    // Re-Cl bonds
    addBond(mol, 0, 2, 1);
    addBond(mol, 0, 3, 1);
    addBond(mol, 0, 4, 1);
    addBond(mol, 0, 5, 1);
    addBond(mol, 1, 6, 1);
    addBond(mol, 1, 7, 1);
    addBond(mol, 1, 8, 1);
    addBond(mol, 1, 9, 1);

    centerMolecule(mol);
}

// Build Tungsten Hexacarbonyl (W(CO)6) - octahedral carbonyl complex
void buildTungstenHexacarbonyl(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "W(CO)6 Hexacarbonyl");

    // Perfect octahedral geometry
    // W at center, 6 CO ligands along ±x, ±y, ±z axes
    float wC = 1.5f;   // W-C bond length
    float cO = 0.9f;   // C-O bond length (triple bond)

    // Central tungsten (using Fe color - metallic)
    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_FE);  // 0 - W

    // CO ligands along +x
    addAtom(mol,  wC, 0.0f, 0.0f, ATOM_C);           // 1 - C
    addAtom(mol,  wC + cO, 0.0f, 0.0f, ATOM_O);      // 2 - O

    // CO ligands along -x
    addAtom(mol, -wC, 0.0f, 0.0f, ATOM_C);           // 3 - C
    addAtom(mol, -wC - cO, 0.0f, 0.0f, ATOM_O);      // 4 - O

    // CO ligands along +y
    addAtom(mol, 0.0f,  wC, 0.0f, ATOM_C);           // 5 - C
    addAtom(mol, 0.0f,  wC + cO, 0.0f, ATOM_O);      // 6 - O

    // CO ligands along -y
    addAtom(mol, 0.0f, -wC, 0.0f, ATOM_C);           // 7 - C
    addAtom(mol, 0.0f, -wC - cO, 0.0f, ATOM_O);      // 8 - O

    // CO ligands along +z
    addAtom(mol, 0.0f, 0.0f,  wC, ATOM_C);           // 9 - C
    addAtom(mol, 0.0f, 0.0f,  wC + cO, ATOM_O);      // 10 - O

    // CO ligands along -z
    addAtom(mol, 0.0f, 0.0f, -wC, ATOM_C);           // 11 - C
    addAtom(mol, 0.0f, 0.0f, -wC - cO, ATOM_O);      // 12 - O

    // W-C bonds (6 bonds to central W)
    addBond(mol, 0, 1, 1);
    addBond(mol, 0, 3, 1);
    addBond(mol, 0, 5, 1);
    addBond(mol, 0, 7, 1);
    addBond(mol, 0, 9, 1);
    addBond(mol, 0, 11, 1);

    // C≡O triple bonds (6 CO groups)
    addBond(mol, 1, 2, 3);
    addBond(mol, 3, 4, 3);
    addBond(mol, 5, 6, 3);
    addBond(mol, 7, 8, 3);
    addBond(mol, 9, 10, 3);
    addBond(mol, 11, 12, 3);

    centerMolecule(mol);
}

// Build o-Carborane (C2B10H12) - icosahedral cage, "sci-fi grenade"
void buildOCarborane(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "o-Carborane (C2B10H12)");

    // Icosahedron vertices using golden ratio
    // φ = (1 + √5) / 2 ≈ 1.618
    const float phi = 1.618033988749895f;
    const float scale = 1.0f;
    const float hDist = 0.8f;  // B-H and C-H bond length extension

    // 12 vertices of icosahedron: (0, ±1, ±φ), (±1, ±φ, 0), (±φ, 0, ±1)
    float verts[12][3] = {
        { 0,  1,  phi}, { 0,  1, -phi}, { 0, -1,  phi}, { 0, -1, -phi},
        { 1,  phi,  0}, {-1,  phi,  0}, { 1, -phi,  0}, {-1, -phi,  0},
        { phi,  0,  1}, {-phi,  0,  1}, { phi,  0, -1}, {-phi,  0, -1}
    };

    // Scale vertices
    for (int i = 0; i < 12; i++) {
        verts[i][0] *= scale;
        verts[i][1] *= scale;
        verts[i][2] *= scale;
    }

    // In o-carborane, carbons are at adjacent positions (vertices 0 and 4 are adjacent)
    // Atoms 0,1 = Carbon, Atoms 2-11 = Boron
    addAtom(mol, verts[0][0], verts[0][1], verts[0][2], ATOM_C);  // 0 - C
    addAtom(mol, verts[4][0], verts[4][1], verts[4][2], ATOM_C);  // 1 - C

    // Remaining 10 borons
    int boronVerts[] = {1, 2, 3, 5, 6, 7, 8, 9, 10, 11};
    for (int i = 0; i < 10; i++) {
        int v = boronVerts[i];
        addAtom(mol, verts[v][0], verts[v][1], verts[v][2], ATOM_B);  // 2-11 - B
    }

    // Add hydrogens radially outward from each cage atom
    for (int i = 0; i < 12; i++) {
        float x, y, z;
        if (i == 0) { x = verts[0][0]; y = verts[0][1]; z = verts[0][2]; }
        else if (i == 1) { x = verts[4][0]; y = verts[4][1]; z = verts[4][2]; }
        else { int v = boronVerts[i-2]; x = verts[v][0]; y = verts[v][1]; z = verts[v][2]; }

        // Normalize and extend for H position
        float len = sqrtf(x*x + y*y + z*z);
        float hx = x + (x/len) * hDist;
        float hy = y + (y/len) * hDist;
        float hz = z + (z/len) * hDist;
        addAtom(mol, hx, hy, hz, ATOM_H);  // 12-23 - H
    }

    // Icosahedron edges (30 edges connecting 12 vertices)
    // Each vertex connects to 5 neighbors
    int edges[30][2] = {
        {0,2}, {0,4}, {0,5}, {0,8}, {0,9},      // vertex 0
        {1,3}, {1,4}, {1,5}, {1,10}, {1,11},    // vertex 1
        {2,6}, {2,7}, {2,8}, {2,9},             // vertex 2
        {3,6}, {3,7}, {3,10}, {3,11},           // vertex 3
        {4,5}, {4,8}, {4,10},                   // vertex 4
        {5,9}, {5,11},                          // vertex 5
        {6,7}, {6,8}, {6,10},                   // vertex 6
        {7,9}, {7,11},                          // vertex 7
        {8,10}, {9,11}                          // vertices 8,9
    };

    // Map original vertex indices to our atom indices
    // Vertex 0 -> atom 0 (C), vertex 4 -> atom 1 (C)
    // Other vertices map to atoms 2-11 (B)
    auto vertToAtom = [&](int v) -> int {
        if (v == 0) return 0;
        if (v == 4) return 1;
        for (int i = 0; i < 10; i++) {
            if (boronVerts[i] == v) return i + 2;
        }
        return -1;
    };

    for (int i = 0; i < 30; i++) {
        int a1 = vertToAtom(edges[i][0]);
        int a2 = vertToAtom(edges[i][1]);
        if (a1 >= 0 && a2 >= 0) {
            addBond(mol, a1, a2, 1);
        }
    }

    // C-H and B-H bonds (cage atoms 0-11 to hydrogens 12-23)
    for (int i = 0; i < 12; i++) {
        addBond(mol, i, i + 12, 1);
    }

    centerMolecule(mol);
}

// Build Dodecaborate ([B12H12]2-) - perfect icosahedron, super clean
void buildDodecaborate(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "[B12H12]2- Dodecaborate");

    // Perfect icosahedron with golden ratio
    const float phi = 1.618033988749895f;
    const float scale = 1.0f;
    const float hDist = 0.85f;  // B-H bond length extension

    // 12 vertices of icosahedron
    float verts[12][3] = {
        { 0,  1,  phi}, { 0,  1, -phi}, { 0, -1,  phi}, { 0, -1, -phi},
        { 1,  phi,  0}, {-1,  phi,  0}, { 1, -phi,  0}, {-1, -phi,  0},
        { phi,  0,  1}, {-phi,  0,  1}, { phi,  0, -1}, {-phi,  0, -1}
    };

    // Scale and add boron atoms at all 12 vertices
    for (int i = 0; i < 12; i++) {
        addAtom(mol, verts[i][0] * scale, verts[i][1] * scale, verts[i][2] * scale, ATOM_B);
    }

    // Add hydrogens radially outward from each boron
    for (int i = 0; i < 12; i++) {
        float x = verts[i][0] * scale;
        float y = verts[i][1] * scale;
        float z = verts[i][2] * scale;
        float len = sqrtf(x*x + y*y + z*z);
        float hx = x + (x/len) * hDist;
        float hy = y + (y/len) * hDist;
        float hz = z + (z/len) * hDist;
        addAtom(mol, hx, hy, hz, ATOM_H);  // 12-23 - H
    }

    // All 30 icosahedral edges (B-B bonds)
    int edges[30][2] = {
        {0,2}, {0,4}, {0,5}, {0,8}, {0,9},
        {1,3}, {1,4}, {1,5}, {1,10}, {1,11},
        {2,6}, {2,7}, {2,8}, {2,9},
        {3,6}, {3,7}, {3,10}, {3,11},
        {4,5}, {4,8}, {4,10},
        {5,9}, {5,11},
        {6,7}, {6,8}, {6,10},
        {7,9}, {7,11},
        {8,10}, {9,11}
    };

    for (int i = 0; i < 30; i++) {
        addBond(mol, edges[i][0], edges[i][1], 1);
    }

    // B-H bonds (borons 0-11 to hydrogens 12-23)
    for (int i = 0; i < 12; i++) {
        addBond(mol, i, i + 12, 1);
    }

    centerMolecule(mol);
}

// Build [2]Catenane - two interlocked rings (molecular chain links)
void buildCatenane(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "[2]Catenane");

    // Two interlocked rings, each with 20 carbons
    // Ring 1: in XY plane, Ring 2: in XZ plane, offset and tilted to interlock
    const int ringSize = 20;
    const float r1 = 2.5f;  // radius of ring 1
    const float r2 = 2.5f;  // radius of ring 2
    const float kPi = 3.14159265359f;

    // Ring 1: circular in XY plane (carbons 0-19)
    for (int i = 0; i < ringSize; i++) {
        float angle = 2.0f * kPi * i / ringSize;
        float x = r1 * cosf(angle);
        float y = r1 * sinf(angle);
        float z = 0.0f;
        addAtom(mol, x, y, z, ATOM_C);
    }

    // Ring 2: circular in XZ plane, offset to pass through ring 1 (carbons 20-39)
    for (int i = 0; i < ringSize; i++) {
        float angle = 2.0f * kPi * i / ringSize;
        float x = r2 * cosf(angle);
        float y = 0.0f;
        float z = r2 * sinf(angle);
        addAtom(mol, x, y, z, ATOM_N);  // Use N for visual distinction
    }

    // Bonds for ring 1
    for (int i = 0; i < ringSize; i++) {
        addBond(mol, i, (i + 1) % ringSize, 1);
    }

    // Bonds for ring 2
    for (int i = 0; i < ringSize; i++) {
        addBond(mol, ringSize + i, ringSize + ((i + 1) % ringSize), 1);
    }

    centerMolecule(mol);
}

// Build Rotaxane - ring on an axle with bulky stoppers
void buildRotaxane(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Rotaxane");

    const float kPi = 3.14159265359f;

    // Central axle: linear chain of 12 carbons along X axis
    const int axleLength = 12;
    const float axleSpacing = 0.8f;
    float axleStart = -axleLength * axleSpacing / 2.0f;

    for (int i = 0; i < axleLength; i++) {
        float x = axleStart + i * axleSpacing;
        addAtom(mol, x, 0.0f, 0.0f, ATOM_C);  // 0-11: axle carbons
    }

    // Ring around the middle (threaded bead) - 12 nitrogens
    const int ringSize = 12;
    const float ringRadius = 1.5f;
    for (int i = 0; i < ringSize; i++) {
        float angle = 2.0f * kPi * i / ringSize;
        float y = ringRadius * cosf(angle);
        float z = ringRadius * sinf(angle);
        addAtom(mol, 0.0f, y, z, ATOM_N);  // 12-23: ring atoms
    }

    // Left stopper (bulky group) - tetrahedral arrangement of 4 oxygens
    float stopperX = axleStart - 0.5f;
    float stopperR = 1.2f;
    addAtom(mol, stopperX, stopperR, 0.0f, ATOM_O);           // 24
    addAtom(mol, stopperX, -stopperR, 0.0f, ATOM_O);          // 25
    addAtom(mol, stopperX, 0.0f, stopperR, ATOM_O);           // 26
    addAtom(mol, stopperX, 0.0f, -stopperR, ATOM_O);          // 27

    // Right stopper
    stopperX = axleStart + (axleLength - 1) * axleSpacing + 0.5f;
    addAtom(mol, stopperX, stopperR, 0.0f, ATOM_O);           // 28
    addAtom(mol, stopperX, -stopperR, 0.0f, ATOM_O);          // 29
    addAtom(mol, stopperX, 0.0f, stopperR, ATOM_O);           // 30
    addAtom(mol, stopperX, 0.0f, -stopperR, ATOM_O);          // 31

    // Axle bonds
    for (int i = 0; i < axleLength - 1; i++) {
        addBond(mol, i, i + 1, 1);
    }

    // Ring bonds
    for (int i = 0; i < ringSize; i++) {
        addBond(mol, 12 + i, 12 + ((i + 1) % ringSize), 1);
    }

    // Left stopper bonds to first axle carbon
    addBond(mol, 0, 24, 1);
    addBond(mol, 0, 25, 1);
    addBond(mol, 0, 26, 1);
    addBond(mol, 0, 27, 1);

    // Right stopper bonds to last axle carbon
    addBond(mol, axleLength - 1, 28, 1);
    addBond(mol, axleLength - 1, 29, 1);
    addBond(mol, axleLength - 1, 30, 1);
    addBond(mol, axleLength - 1, 31, 1);

    centerMolecule(mol);
}

// Build Molecular Trefoil Knot - a knotted loop
void buildTrefoilKnot(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Molecular Trefoil Knot");

    // Parametric trefoil knot:
    // x = sin(t) + 2*sin(2t)
    // y = cos(t) - 2*cos(2t)
    // z = -sin(3t)
    const int numAtoms = 60;  // Number of atoms in the knot
    const float kPi = 3.14159265359f;
    const float scale = 0.8f;

    for (int i = 0; i < numAtoms; i++) {
        float t = 2.0f * kPi * i / numAtoms;
        float x = (sinf(t) + 2.0f * sinf(2.0f * t)) * scale;
        float y = (cosf(t) - 2.0f * cosf(2.0f * t)) * scale;
        float z = -sinf(3.0f * t) * scale;

        // Alternate colors for visual effect (C, N, O pattern)
        int atomType = (i % 3 == 0) ? ATOM_C : ((i % 3 == 1) ? ATOM_N : ATOM_O);
        addAtom(mol, x, y, z, atomType);
    }

    // Connect atoms in sequence to form the knot
    for (int i = 0; i < numAtoms; i++) {
        addBond(mol, i, (i + 1) % numAtoms, 1);
    }

    centerMolecule(mol);
}

// Build [1.1.1]Propellane (C5H6) - inverted internal C-C bond under tension
void buildPropellane(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "[1.1.1]Propellane (C5H6)");

    // Structure: two bridgehead carbons connected by 3 methylene bridges
    // and a direct "inverted" bond between them
    // Bridgehead carbons at top and bottom
    float bh = 0.8f;  // bridgehead separation (half)
    float br = 1.2f;  // bridge radius

    // Bridgehead carbons (the famous inverted bond between these)
    addAtom(mol, 0.0f, 0.0f, -bh, ATOM_C);  // 0 - bottom bridgehead
    addAtom(mol, 0.0f, 0.0f,  bh, ATOM_C);  // 1 - top bridgehead

    // Three bridge carbons (CH2) arranged in triangle around z-axis
    float angle1 = 0.0f;
    float angle2 = 2.094395f;  // 120 degrees
    float angle3 = 4.188790f;  // 240 degrees

    addAtom(mol, br * cosf(angle1), br * sinf(angle1), 0.0f, ATOM_C);  // 2
    addAtom(mol, br * cosf(angle2), br * sinf(angle2), 0.0f, ATOM_C);  // 3
    addAtom(mol, br * cosf(angle3), br * sinf(angle3), 0.0f, ATOM_C);  // 4

    // Hydrogens on bridge carbons (2 per carbon, pointing outward)
    float hDist = 0.8f;
    for (int i = 0; i < 3; i++) {
        float angle = i * 2.094395f;
        float cx = br * cosf(angle);
        float cy = br * sinf(angle);
        // Two H's per bridge carbon, displaced up/down and outward
        float hx = (br + hDist) * cosf(angle);
        float hy = (br + hDist) * sinf(angle);
        addAtom(mol, hx, hy, -0.5f, ATOM_H);  // 5,7,9
        addAtom(mol, hx, hy,  0.5f, ATOM_H);  // 6,8,10
    }

    // The famous inverted bond between bridgeheads
    addBond(mol, 0, 1, 1);

    // Bridgehead to bridge bonds
    addBond(mol, 0, 2, 1); addBond(mol, 1, 2, 1);
    addBond(mol, 0, 3, 1); addBond(mol, 1, 3, 1);
    addBond(mol, 0, 4, 1); addBond(mol, 1, 4, 1);

    // C-H bonds
    addBond(mol, 2, 5, 1); addBond(mol, 2, 6, 1);
    addBond(mol, 3, 7, 1); addBond(mol, 3, 8, 1);
    addBond(mol, 4, 9, 1); addBond(mol, 4, 10, 1);

    centerMolecule(mol);
}

// Build Prismane (C6H6) - triangular prism of carbons
void buildPrismane(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Prismane (C6H6)");

    // Triangular prism: 2 triangular faces connected by 3 rectangular faces
    float r = 1.0f;   // radius of triangular face
    float h = 0.9f;   // half-height of prism

    // Bottom triangle (z = -h)
    float a1 = 0.0f, a2 = 2.094395f, a3 = 4.188790f;  // 0, 120, 240 degrees
    addAtom(mol, r * cosf(a1), r * sinf(a1), -h, ATOM_C);  // 0
    addAtom(mol, r * cosf(a2), r * sinf(a2), -h, ATOM_C);  // 1
    addAtom(mol, r * cosf(a3), r * sinf(a3), -h, ATOM_C);  // 2

    // Top triangle (z = +h)
    addAtom(mol, r * cosf(a1), r * sinf(a1),  h, ATOM_C);  // 3
    addAtom(mol, r * cosf(a2), r * sinf(a2),  h, ATOM_C);  // 4
    addAtom(mol, r * cosf(a3), r * sinf(a3),  h, ATOM_C);  // 5

    // Hydrogens pointing outward from each carbon
    float hDist = 0.8f;
    for (int i = 0; i < 3; i++) {
        float angle = i * 2.094395f;
        float hx = (r + hDist) * cosf(angle);
        float hy = (r + hDist) * sinf(angle);
        addAtom(mol, hx, hy, -h - 0.3f, ATOM_H);  // 6,7,8 - bottom H's
        addAtom(mol, hx, hy,  h + 0.3f, ATOM_H);  // 9,10,11 - top H's
    }

    // Bottom triangle bonds
    addBond(mol, 0, 1, 1); addBond(mol, 1, 2, 1); addBond(mol, 2, 0, 1);
    // Top triangle bonds
    addBond(mol, 3, 4, 1); addBond(mol, 4, 5, 1); addBond(mol, 5, 3, 1);
    // Vertical bonds (connecting triangles)
    addBond(mol, 0, 3, 1); addBond(mol, 1, 4, 1); addBond(mol, 2, 5, 1);

    // C-H bonds
    addBond(mol, 0, 6, 1); addBond(mol, 1, 7, 1); addBond(mol, 2, 8, 1);
    addBond(mol, 3, 9, 1); addBond(mol, 4, 10, 1); addBond(mol, 5, 11, 1);

    centerMolecule(mol);
}

// Build Tetrahedrane (C4H4) - tetrahedral carbon cage
void buildTetrahedrane(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Tetrahedrane (C4H4)");

    // Perfect tetrahedron vertices
    float a = 1.0f;  // scale
    // Tetrahedral coordinates:
    // (1,1,1), (1,-1,-1), (-1,1,-1), (-1,-1,1) normalized
    float v = a * 0.577350269f;  // 1/sqrt(3)

    addAtom(mol,  v,  v,  v, ATOM_C);  // 0
    addAtom(mol,  v, -v, -v, ATOM_C);  // 1
    addAtom(mol, -v,  v, -v, ATOM_C);  // 2
    addAtom(mol, -v, -v,  v, ATOM_C);  // 3

    // Hydrogens pointing outward from each carbon
    float hDist = 0.8f;
    float hv = (a + hDist) * 0.577350269f;
    addAtom(mol,  hv,  hv,  hv, ATOM_H);  // 4
    addAtom(mol,  hv, -hv, -hv, ATOM_H);  // 5
    addAtom(mol, -hv,  hv, -hv, ATOM_H);  // 6
    addAtom(mol, -hv, -hv,  hv, ATOM_H);  // 7

    // All 6 edges of tetrahedron (C-C bonds)
    addBond(mol, 0, 1, 1); addBond(mol, 0, 2, 1); addBond(mol, 0, 3, 1);
    addBond(mol, 1, 2, 1); addBond(mol, 1, 3, 1); addBond(mol, 2, 3, 1);

    // C-H bonds
    addBond(mol, 0, 4, 1); addBond(mol, 1, 5, 1);
    addBond(mol, 2, 6, 1); addBond(mol, 3, 7, 1);

    centerMolecule(mol);
}

// Build Dewar Benzene (C6H6) - bicyclic bent benzene isomer
void buildDewarBenzene(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Dewar Benzene (C6H6)");

    // Structure: bicyclo[2.2.0]hexa-2,5-diene
    // Two fused cyclobutene rings forming a "butterfly" shape
    float w = 1.0f;   // width
    float h = 0.7f;   // height of fold
    float d = 0.6f;   // depth

    // Four carbons of the central 4-membered ring (folded)
    addAtom(mol, -w/2, -d,  h, ATOM_C);  // 0 - back left (up)
    addAtom(mol,  w/2, -d,  h, ATOM_C);  // 1 - back right (up)
    addAtom(mol,  w/2,  d, -h, ATOM_C);  // 2 - front right (down)
    addAtom(mol, -w/2,  d, -h, ATOM_C);  // 3 - front left (down)

    // Two carbons bridging across (the "wings")
    addAtom(mol, -w/2, 0.0f, 0.0f, ATOM_C);  // 4 - left bridge
    addAtom(mol,  w/2, 0.0f, 0.0f, ATOM_C);  // 5 - right bridge

    // Hydrogens on each carbon
    float hDist = 0.8f;
    addAtom(mol, -w/2 - hDist, -d,  h + 0.3f, ATOM_H);  // 6
    addAtom(mol,  w/2 + hDist, -d,  h + 0.3f, ATOM_H);  // 7
    addAtom(mol,  w/2 + hDist,  d, -h - 0.3f, ATOM_H);  // 8
    addAtom(mol, -w/2 - hDist,  d, -h - 0.3f, ATOM_H);  // 9
    addAtom(mol, -w/2 - hDist, 0.0f, 0.0f, ATOM_H);     // 10
    addAtom(mol,  w/2 + hDist, 0.0f, 0.0f, ATOM_H);     // 11

    // Central 4-membered ring
    addBond(mol, 0, 1, 1); addBond(mol, 1, 2, 1);
    addBond(mol, 2, 3, 1); addBond(mol, 3, 0, 1);

    // Bridge bonds (double bond character in real molecule)
    addBond(mol, 0, 4, 2); addBond(mol, 3, 4, 1);
    addBond(mol, 1, 5, 2); addBond(mol, 2, 5, 1);

    // C-H bonds
    addBond(mol, 0, 6, 1); addBond(mol, 1, 7, 1);
    addBond(mol, 2, 8, 1); addBond(mol, 3, 9, 1);
    addBond(mol, 4, 10, 1); addBond(mol, 5, 11, 1);

    centerMolecule(mol);
}

// Build Hydroxyproline (C5H9NO3)
void buildHydroxyproline(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Hydroxyproline (C5H9NO3)");

    // 5-membered ring with N and OH on one carbon
    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_C);      // Alpha C
    addAtom(mol, -1.2f, 0.0f, 0.0f, ATOM_N);     // N (in ring)
    addAtom(mol, 1.3f, 0.7f, 0.0f, ATOM_C);      // COOH C
    addAtom(mol, 1.3f, 2.0f, 0.0f, ATOM_O);
    addAtom(mol, 2.4f, 0.0f, 0.0f, ATOM_O);
    addAtom(mol, 0.0f, -1.5f, 0.0f, ATOM_C);     // Ring C with OH
    addAtom(mol, -0.8f, -2.2f, 0.0f, ATOM_O);    // OH on ring
    addAtom(mol, -1.4f, -1.5f, 0.0f, ATOM_C);    // Ring C
    addAtom(mol, -1.8f, -0.05f, 0.0f, ATOM_H);   // H on N
    addAtom(mol, 0.0f, 0.6f, 0.9f, ATOM_H);
    addAtom(mol, 3.2f, 0.5f, 0.0f, ATOM_H);
    addAtom(mol, 0.5f, -2.0f, 0.9f, ATOM_H);
    addAtom(mol, -1.9f, -2.0f, 0.9f, ATOM_H);
    addAtom(mol, -1.9f, -2.0f, -0.9f, ATOM_H);
    addAtom(mol, -1.5f, -2.8f, 0.0f, ATOM_H);    // H on OH

    addBond(mol, 0, 1, 1);
    addBond(mol, 0, 2, 1);
    addBond(mol, 2, 3, 2);
    addBond(mol, 2, 4, 1);
    addBond(mol, 0, 5, 1);
    addBond(mol, 5, 6, 1);
    addBond(mol, 5, 7, 1);
    addBond(mol, 7, 1, 1);
    addBond(mol, 1, 8, 1);
    addBond(mol, 0, 9, 1);
    addBond(mol, 4, 10, 1);
    addBond(mol, 5, 11, 1);
    addBond(mol, 7, 12, 1);
    addBond(mol, 7, 13, 1);
    addBond(mol, 6, 14, 1);

    centerMolecule(mol);
}

// Build Niacinamide/Vitamin B3 (C6H6N2O)
void buildNiacinamide(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Niacinamide/Vitamin B3 (C6H6N2O)");

    // Pyridine ring
    float r = 1.4f;
    for (int i = 0; i < 6; i++) {
        float angle = i * PI / 3.0f;
        if (i == 0) {
            addAtom(mol, r * cosf(angle), r * sinf(angle), 0.0f, ATOM_N);  // N at position 1
        } else {
            addAtom(mol, r * cosf(angle), r * sinf(angle), 0.0f, ATOM_C);  // C at 2-6
        }
    }

    // Amide group (-CONH2) at position 3
    addAtom(mol, 2.4f, 0.8f, 0.0f, ATOM_C);      // C=O carbon
    addAtom(mol, 3.5f, 1.2f, 0.0f, ATOM_O);      // =O
    addAtom(mol, 2.4f, -0.3f, 0.0f, ATOM_N);     // NH2
    addAtom(mol, 3.3f, -0.8f, 0.0f, ATOM_H);     // H on N
    addAtom(mol, 1.5f, -0.7f, 0.0f, ATOM_H);     // H on N

    // Hydrogens on ring
    float hR = 2.4f;
    for (int i = 0; i < 6; i++) {
        if (i != 2) {  // Skip position 3 where amide is
            float angle = i * PI / 3.0f;
            addAtom(mol, hR * cosf(angle), hR * sinf(angle), 0.0f, ATOM_H);
        }
    }

    // Ring bonds
    for (int i = 0; i < 6; i++) {
        addBond(mol, i, (i + 1) % 6, (i % 2 == 0) ? 2 : 1);
    }

    // Amide bonds
    addBond(mol, 2, 6, 1);   // C-CO
    addBond(mol, 6, 7, 2);   // C=O
    addBond(mol, 6, 8, 1);   // C-N
    addBond(mol, 8, 9, 1);   // N-H
    addBond(mol, 8, 10, 1);  // N-H

    // Hydrogen bonds
    addBond(mol, 1, 11, 1);
    addBond(mol, 3, 12, 1);
    addBond(mol, 4, 13, 1);
    addBond(mol, 5, 14, 1);

    centerMolecule(mol);
}

// Build Glycerin/Glycerol (C3H8O3)
void buildGlycerin(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Glycerin/Glycerol (C3H8O3)");

    // Three carbons in a chain
    addAtom(mol, -1.5f, 0.0f, 0.0f, ATOM_C);     // C1
    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_C);      // C2 (central)
    addAtom(mol, 1.5f, 0.0f, 0.0f, ATOM_C);      // C3

    // Three OH groups
    addAtom(mol, -1.5f, 1.2f, 0.0f, ATOM_O);     // OH on C1
    addAtom(mol, 0.0f, 1.2f, 0.0f, ATOM_O);      // OH on C2
    addAtom(mol, 1.5f, 1.2f, 0.0f, ATOM_O);      // OH on C3

    // Hydrogens on carbons and oxygens
    addAtom(mol, -1.5f, -0.8f, 0.85f, ATOM_H);   // H on C1
    addAtom(mol, -1.5f, -0.8f, -0.85f, ATOM_H);  // H on C1
    addAtom(mol, 0.0f, -0.8f, 0.85f, ATOM_H);    // H on C2
    addAtom(mol, -2.3f, 1.6f, 0.0f, ATOM_H);     // H on OH
    addAtom(mol, -0.8f, 1.6f, 0.0f, ATOM_H);     // H on OH
    addAtom(mol, 2.3f, 1.6f, 0.0f, ATOM_H);      // H on OH
    addAtom(mol, 1.5f, -0.8f, 0.85f, ATOM_H);    // H on C3
    addAtom(mol, 1.5f, -0.8f, -0.85f, ATOM_H);   // H on C3

    // Bonds
    addBond(mol, 0, 1, 1);  // C-C
    addBond(mol, 1, 2, 1);  // C-C
    addBond(mol, 0, 3, 1);  // C-O
    addBond(mol, 1, 4, 1);  // C-O
    addBond(mol, 2, 5, 1);  // C-O
    addBond(mol, 0, 6, 1);  // C-H
    addBond(mol, 0, 7, 1);  // C-H
    addBond(mol, 1, 8, 1);  // C-H
    addBond(mol, 3, 9, 1);  // O-H
    addBond(mol, 4, 10, 1); // O-H
    addBond(mol, 5, 11, 1); // O-H
    addBond(mol, 2, 12, 1); // C-H
    addBond(mol, 2, 13, 1); // C-H

    centerMolecule(mol);
}

// Build Glycolic Acid (C2H4O3)
void buildGlycolicAcid(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Glycolic Acid (C2H4O3)");

    // CH2-OH and COOH groups
    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_C);      // CH2-OH carbon
    addAtom(mol, 1.5f, 0.0f, 0.0f, ATOM_C);      // Carboxylic acid carbon
    addAtom(mol, 2.2f, 1.1f, 0.0f, ATOM_O);      // =O
    addAtom(mol, 2.2f, -1.1f, 0.0f, ATOM_O);     // OH
    addAtom(mol, 0.0f, 1.2f, 0.0f, ATOM_O);      // OH on CH2
    addAtom(mol, 3.1f, -1.1f, 0.0f, ATOM_H);     // H on carboxylic acid
    addAtom(mol, 0.8f, 1.6f, 0.0f, ATOM_H);      // H on CH2 OH
    addAtom(mol, -0.9f, 0.5f, 0.0f, ATOM_H);     // H on CH2
    addAtom(mol, -0.9f, -0.5f, 0.0f, ATOM_H);    // H on CH2

    addBond(mol, 0, 1, 1);  // C-C
    addBond(mol, 1, 2, 2);  // C=O
    addBond(mol, 1, 3, 1);  // C-OH
    addBond(mol, 0, 4, 1);  // C-OH
    addBond(mol, 3, 5, 1);  // O-H
    addBond(mol, 4, 6, 1);  // O-H
    addBond(mol, 0, 7, 1);  // C-H
    addBond(mol, 0, 8, 1);  // C-H

    centerMolecule(mol);
}

// Build Salicylic Acid (C7H6O3) - BHA
void buildSalicylicAcid(Molecule* mol) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strcpy(mol->name, "Salicylic Acid (C7H6O3)");

    // Benzene ring with OH and COOH substituents
    float r = 1.4f;
    for (int i = 0; i < 6; i++) {
        float angle = i * PI / 3.0f;
        addAtom(mol, r * cosf(angle), r * sinf(angle), 0.0f, ATOM_C);  // 0-5
    }

    // COOH group at position 0
    addAtom(mol, r + 1.3f, 0.0f, 0.0f, ATOM_C);  // 6: COOH C
    addAtom(mol, r + 2.0f, 1.1f, 0.0f, ATOM_O);  // 7: =O
    addAtom(mol, r + 2.0f, -1.1f, 0.0f, ATOM_O); // 8: OH
    addAtom(mol, r + 2.9f, -1.1f, 0.0f, ATOM_H); // 9: H on COOH

    // OH group at position 1 (ortho to COOH)
    addAtom(mol, r * cosf(PI/3) + 0.5f, r * sinf(PI/3) + 0.9f, 0.0f, ATOM_O);  // 10: OH
    addAtom(mol, r * cosf(PI/3) + 1.2f, r * sinf(PI/3) + 1.4f, 0.0f, ATOM_H);  // 11: H on OH

    // Hydrogens on benzene
    float hR = 2.4f;
    addAtom(mol, hR * cosf(2*PI/3), hR * sinf(2*PI/3), 0.0f, ATOM_H);  // 12: H on C2
    addAtom(mol, hR * cosf(PI), hR * sinf(PI), 0.0f, ATOM_H);           // 13: H on C3
    addAtom(mol, hR * cosf(4*PI/3), hR * sinf(4*PI/3), 0.0f, ATOM_H);   // 14: H on C4
    addAtom(mol, hR * cosf(5*PI/3), hR * sinf(5*PI/3), 0.0f, ATOM_H);   // 15: H on C5

    // Benzene ring bonds
    for (int i = 0; i < 6; i++) {
        addBond(mol, i, (i + 1) % 6, (i % 2 == 0) ? 2 : 1);
    }

    // COOH bonds
    addBond(mol, 0, 6, 1);  // C-COOH
    addBond(mol, 6, 7, 2);  // C=O
    addBond(mol, 6, 8, 1);  // C-OH
    addBond(mol, 8, 9, 1);  // O-H

    // OH bonds
    addBond(mol, 1, 10, 1); // C-OH
    addBond(mol, 10, 11, 1);// O-H

    // Hydrogen bonds
    addBond(mol, 2, 12, 1);
    addBond(mol, 3, 13, 1);
    addBond(mol, 4, 14, 1);
    addBond(mol, 5, 15, 1);

    centerMolecule(mol);
}

// ============== MOLECULE REGISTRY (258 molecules) ==============

typedef void (*MoleculeBuilder)(Molecule*);

struct MoleculeInfo {
    MoleculeBuilder builder;
    const char* name;
    int category;
    const char* description;
    const char* longDescription;
};

static MoleculeInfo molecules[] = {
    // === SIMPLE MOLECULES (0-31 original + gases) ===
    { buildWater, "Water", CAT_SIMPLE, "H2O - Universal solvent",
      "Water is essential for all life on Earth. Made of two hydrogen atoms bonded to one oxygen atom, it covers 71% of our planet's surface. Its unique properties allow it to dissolve many substances, earning it the nickname 'universal solvent.'" },
    { buildMethane, "Methane", CAT_SIMPLE, "CH4 - Natural gas",
      "Methane is the simplest hydrocarbon and main component of natural gas. Produced by decomposing organic matter and livestock, it's both a valuable fuel and a potent greenhouse gas about 80 times more warming than CO2 over 20 years." },
    { buildBenzene, "Benzene", CAT_ORGANIC, "C6H6 - Aromatic ring",
      "Benzene is a ring-shaped molecule that forms the backbone of countless organic compounds. Discovered in 1825, its unique structure puzzled chemists until Kekule dreamed of a snake eating its tail. It's used to make plastics, rubber, and dyes." },
    { buildEthanol, "Ethanol", CAT_ORGANIC, "C2H5OH - Alcohol",
      "Ethanol is the alcohol in beer, wine, and spirits, produced by yeast fermenting sugars. Beyond beverages, it's used as fuel, solvent, and antiseptic. Your liver can process about one drink per hour, breaking ethanol into acetaldehyde." },
    { buildCaffeine, "Caffeine", CAT_PHARMA, "C8H10N4O2 - Stimulant",
      "Caffeine is the world's most popular psychoactive drug, found in coffee, tea, and chocolate. It works by blocking adenosine receptors in your brain, preventing drowsiness. About 90% of adults worldwide consume caffeine daily." },
    { buildAdenine, "Adenine", CAT_ENERGY, "DNA nucleobase",
      "Adenine is one of four nucleobases in DNA's genetic code, always pairing with thymine. Beyond DNA, adenine is part of ATP, your body's energy currency. Every cell in your body contains about one billion ATP molecules." },
    { buildGlucose, "Glucose", CAT_SUGARS, "C6H12O6 - Blood sugar",
      "Glucose is your body's primary fuel source, powering every cell from your brain to your muscles. Your blood contains about 4 grams at any time. Plants create glucose through photosynthesis, converting sunlight into stored energy." },
    { buildAspirin, "Aspirin/Bayer", CAT_PHARMA, "Pain reliever",
      "Aspirin is one of the oldest and most widely used medicines, derived from willow bark used since ancient Egypt. It reduces pain, fever, and inflammation by blocking prostaglandin production. Low doses help prevent heart attacks by thinning blood." },
    { buildAmmonia, "Ammonia", CAT_SIMPLE, "NH3 - Cleaning agent",
      "Ammonia is a pungent gas essential for life on Earth. The Haber process converts nitrogen and hydrogen into ammonia for fertilizers, feeding nearly half the world's population. Without synthetic ammonia, we couldn't grow enough food for everyone." },
    { buildCarbonDioxide, "CO2", CAT_SIMPLE, "Carbon dioxide",
      "Carbon dioxide is what you exhale with every breath and what plants breathe in. While essential for photosynthesis, rising CO2 levels from burning fossil fuels trap heat and warm our planet. Current levels are higher than any time in 800,000 years." },
    { buildFormaldehyde, "Formaldehyde", CAT_ORGANIC, "CH2O - Preservative",
      "Formaldehyde is the simplest aldehyde, used to preserve biological specimens and make resins. Your body actually produces small amounts during metabolism. It's found naturally in fruits, vegetables, and even the air we breathe." },
    { buildAcetone, "Acetone", CAT_ORGANIC, "C3H6O - Solvent",
      "Acetone is the main ingredient in nail polish remover, evaporating quickly with a distinctive sweet smell. Your body naturally produces acetone when breaking down fat for energy. Diabetics may produce excess acetone, detectable on their breath." },
    { buildAceticAcid, "Acetic Acid", CAT_ACIDS, "CH3COOH - Vinegar",
      "Acetic acid gives vinegar its sour taste and pungent smell. Created when bacteria ferment alcohol, humans have used vinegar for over 10,000 years for cooking, cleaning, and medicine. Vinegar is typically 5-8% acetic acid in water." },
    { buildPropane, "Propane", CAT_ORGANIC, "C3H8 - Fuel",
      "Propane is a clean-burning fuel used in grills, heaters, and vehicles. Stored as a liquid under pressure, it vaporizes when released. Propane produces fewer emissions than gasoline and is a byproduct of natural gas processing and oil refining." },
    { buildButane, "Butane", CAT_ORGANIC, "C4H10 - Lighter fuel",
      "Butane is the fuel in disposable lighters and portable stoves. At room temperature, it's a gas, but easily liquefies under pressure. Butane is also used as a propellant in aerosol sprays and in the production of synthetic rubber." },
    { buildCyclohexane, "Cyclohexane", CAT_ORGANIC, "C6H12 - Cyclic alkane",
      "Cyclohexane is a ring of six carbon atoms that constantly flips between 'chair' and 'boat' shapes. Used to make nylon and as a solvent, it's found in petroleum. Unlike flat benzene, cyclohexane's carbons zigzag up and down." },
    { buildNaphthalene, "Naphthalene", CAT_ORGANIC, "C10H8 - Mothballs",
      "Naphthalene is the distinctive smell of mothballs, made of two fused benzene rings. It slowly sublimes from solid to gas, killing moths and their larvae. Coal tar and petroleum contain naphthalene, which was once used in early plastics." },
    { buildUrea, "Urea", CAT_ORGANIC, "CH4N2O - Metabolic waste",
      "Urea is how your body safely excretes nitrogen from protein breakdown. Your kidneys filter about 30 grams daily into urine. First synthesized in 1828 by Wohler, urea proved organic compounds could be made without living organisms, revolutionizing chemistry." },
    { buildGlycine, "Glycine", CAT_AMINO_ACIDS, "Simplest amino acid",
      "Glycine is the smallest amino acid, found in every protein in your body. It makes up one-third of collagen, giving skin its structure. As a neurotransmitter, glycine helps regulate sleep and calm the nervous system." },
    { buildAlanine, "Alanine", CAT_AMINO_ACIDS, "Nonpolar amino acid",
      "Alanine is a simple amino acid that helps your body convert glucose to energy. During intense exercise, muscles release alanine to the liver, which converts it back to glucose. This glucose-alanine cycle helps fuel working muscles." },
    { buildThymine, "Thymine", CAT_ENERGY, "DNA nucleobase (T)",
      "Thymine is one of DNA's four letters, always pairing with adenine. Found only in DNA (not RNA), thymine helps distinguish your genetic material from viruses. The double bond in thymine can form dangerous mutations when hit by UV light." },
    { buildCytosine, "Cytosine", CAT_ENERGY, "DNA/RNA nucleobase (C)",
      "Cytosine pairs with guanine in DNA and RNA, held together by three hydrogen bonds. This makes C-G pairs stronger than A-T pairs. Cytosine can spontaneously change to uracil, causing mutations that your cells constantly repair." },
    { buildGuanine, "Guanine", CAT_ENERGY, "DNA/RNA nucleobase (G)",
      "Guanine is a DNA letter that pairs with cytosine. Its crystals create the silvery shimmer in fish scales and give some insects their iridescent colors. The name comes from guano, the bird droppings where it was first discovered." },
    { buildDopamine, "Dopamine", CAT_NEUROTRANS, "Reward neurotransmitter",
      "Dopamine is your brain's reward chemical, surging when you eat delicious food, exercise, or achieve goals. It drives motivation and pleasure, but can be hijacked by addictive substances. Parkinson's disease results from dying dopamine neurons." },
    { buildSerotonin, "Serotonin", CAT_NEUROTRANS, "Mood neurotransmitter",
      "Serotonin regulates mood, sleep, and appetite. Surprisingly, 95% of your body's serotonin is in your gut, not your brain. Many antidepressants work by increasing serotonin levels. Sunlight and exercise naturally boost serotonin production." },
    { buildNitricOxide, "NO", CAT_SIMPLE, "Nitric oxide",
      "Nitric oxide is a tiny signaling molecule your body uses to dilate blood vessels. Discovered in the 1980s, it won the Nobel Prize in 1998. Viagra works by enhancing nitric oxide's effects. Lightning also produces NO in the atmosphere." },
    { buildHydrogenPeroxide, "H2O2", CAT_HOUSEHOLD, "Hydrogen peroxide",
      "Hydrogen peroxide is water with an extra oxygen atom, making it a powerful oxidizer. Dilute solutions clean wounds by killing bacteria with oxygen bubbles. Your white blood cells produce H2O2 to destroy invading pathogens." },
    { buildSulfuricAcid, "H2SO4", CAT_ACIDS, "Sulfuric acid",
      "Sulfuric acid is the most produced chemical in the world, used in fertilizers, petroleum refining, and batteries. It's so corrosive it can char paper and dissolve metals. Acid rain forms when SO2 emissions create sulfuric acid in clouds." },
    { buildPhosphoricAcid, "H3PO4", CAT_ACIDS, "Phosphoric acid",
      "Phosphoric acid gives cola its tangy taste and prevents bacterial growth. It's also used to remove rust and make fertilizers. Your teeth can be damaged by excessive cola consumption due to this acid's erosive properties." },
    { buildToluene, "Toluene", CAT_ORGANIC, "C7H8 - Solvent",
      "Toluene is benzene with a methyl group attached, smelling like paint thinner. It's in gasoline and used to make TNT explosive. Toluene dissolves many substances, making it useful in paints, but its fumes can cause dizziness." },
    { buildPhenol, "Phenol", CAT_ORGANIC, "C6H5OH - Carbolic acid",
      "Phenol was the first surgical antiseptic, revolutionizing medicine when Joseph Lister used it in 1867. It's now used in plastics, adhesives, and aspirin production. Despite being toxic, dilute phenol is still used in throat sprays." },
    { buildAcetylene, "Acetylene", CAT_ORGANIC, "C2H2 - Welding gas",
      "Acetylene burns with the hottest flame of any fuel gas, reaching 3,500 degrees Celsius. Welders use it to cut and join metals. First made from calcium carbide and water, acetylene once powered bicycle and car headlamps." },
    // === BATCH 1: Gases and solvents (32-63) ===
    { buildOxygen, "O2", CAT_SIMPLE, "Molecular oxygen",
      "Oxygen makes up 21% of Earth's atmosphere and is essential for most life. Produced by photosynthesis 2.4 billion years ago, it transformed our planet. Your body uses oxygen to burn glucose for energy, producing water and CO2 as byproducts." },
    { buildNitrogen, "N2", CAT_SIMPLE, "Molecular nitrogen",
      "Nitrogen comprises 78% of the air you breathe but is chemically unreactive due to its strong triple bond. Plants can't use atmospheric nitrogen directly; bacteria must first 'fix' it into ammonia. Your body contains about 3% nitrogen by weight." },
    { buildHydrogen, "H2", CAT_SIMPLE, "Molecular hydrogen",
      "Hydrogen is the lightest and most abundant element in the universe, making up 75% of all matter. Stars like our Sun fuse hydrogen into helium, releasing enormous energy. Hydrogen may become a clean fuel of the future, producing only water when burned." },
    { buildOzone, "O3", CAT_SIMPLE, "Ozone",
      "Ozone has a distinctive sharp smell after thunderstorms. High in the stratosphere, it shields life from harmful UV radiation. Near the ground, it's a pollutant causing respiratory problems. The ozone hole was caused by CFCs but is now slowly healing." },
    { buildCarbonMonoxide, "CO", CAT_SIMPLE, "Carbon monoxide",
      "Carbon monoxide is the silent killer, an odorless gas from incomplete combustion. It binds to hemoglobin 200 times stronger than oxygen, preventing blood from carrying oxygen. This is why CO detectors are essential in homes with gas appliances." },
    { buildNitrousOxide, "N2O", CAT_SIMPLE, "Laughing gas",
      "Nitrous oxide is used as an anesthetic in dentistry and as a propellant in whipped cream cans. When inhaled, it produces euphoria and mild hallucinations. It's also a greenhouse gas 300 times more potent than CO2 and depletes the ozone layer." },
    { buildSulfurDioxide, "SO2", CAT_SIMPLE, "Sulfur dioxide",
      "Sulfur dioxide has a choking smell from burning coal and volcanic eruptions. It causes acid rain when it reacts with water in clouds. Winemakers use small amounts to preserve wine, and it's been used this way since Roman times." },
    { buildHydrogenChloride, "HCl", CAT_ACIDS, "Hydrochloric acid",
      "Hydrochloric acid is the acid in your stomach, strong enough to dissolve metal but prevented from damaging your stomach by a thick mucus lining. Your stomach produces about 2 liters of gastric juice daily, with HCl killing bacteria in food." },
    { buildNitricAcid, "HNO3", CAT_ACIDS, "Nitric acid",
      "Nitric acid is a powerfully corrosive acid used to make fertilizers and explosives. It turns proteins yellow through a reaction called xanthoproteic. Alchemists called it 'aqua fortis' (strong water), and mixed with hydrochloric acid, it dissolves gold." },
    { buildMethanol, "Methanol", CAT_ORGANIC, "Wood alcohol",
      "Methanol is the simplest alcohol, extremely toxic if consumed. Just 10 mL can cause blindness; 30 mL can be fatal. It's used as antifreeze and racing fuel. During Prohibition, contaminated alcohol killed thousands before people understood methanol's dangers." },
    { buildEthane, "Ethane", CAT_ORGANIC, "C2H6",
      "Ethane is the second-simplest hydrocarbon, found in natural gas. On Saturn's moon Titan, ethane exists as liquid lakes at -180 degrees Celsius. Industrially, ethane is cracked into ethylene, the most produced organic compound globally." },
    { buildPropene, "Propene", CAT_ORGANIC, "Propylene",
      "Propene is converted into polypropylene plastic, found in food containers, ropes, and car parts. It's the second most produced organic chemical after ethylene. Your body actually produces tiny amounts of propene during fat metabolism." },
    { buildIsopropanol, "Isopropanol", CAT_ORGANIC, "Rubbing alcohol",
      "Isopropanol is the rubbing alcohol in first aid kits, evaporating quickly and killing bacteria. Unlike ethanol, drinking isopropanol is dangerous and can cause blindness. It's also used to clean electronics because it evaporates without leaving residue." },
    { buildEthyleneGlycol, "Ethylene Glycol", CAT_ORGANIC, "Antifreeze",
      "Ethylene glycol is the main ingredient in antifreeze, lowering water's freezing point in car radiators. Its sweet taste makes it dangerous to pets and children. Poisoning damages kidneys as the body converts it to oxalic acid crystals." },
    { buildGlycerol, "Glycerol", CAT_ORGANIC, "Glycerin",
      "Glycerol is a sweet, syrupy liquid forming the backbone of all fats in your body. It moisturizes skin in lotions and sweetens foods. Combined with nitric acid, it becomes nitroglycerin explosive. Your body breaks down fat to release glycerol for energy." },
    { buildAcetaldehyde, "Acetaldehyde", CAT_ORGANIC, "Ethanal",
      "Acetaldehyde causes hangover symptoms when your liver breaks down alcohol. It's more toxic than ethanol itself. Some people lack the enzyme to process acetaldehyde quickly, causing facial flushing when drinking. It's also found in ripe fruits and coffee." },
    { buildFormicAcid, "Formic Acid", CAT_ACIDS, "Ant venom",
      "Formic acid gives ant stings their burning pain, named from 'formica,' Latin for ant. Fire ant stings inject formic acid under your skin. Despite being an acid, dilute formic acid is used as a food preservative and to de-ice airplanes." },
    { buildLacticAcid, "Lactic Acid", CAT_ACIDS, "Muscle fatigue",
      "Lactic acid builds up in muscles during intense exercise, contributing to that burning sensation. Your body actually uses lactate as fuel; muscles can convert it back to glucose. Yogurt and sauerkraut get their tangy taste from bacterial lactic acid production." },
    { buildEthylAcetate, "Ethyl Acetate", CAT_ORGANIC, "Nail polish solvent",
      "Ethyl acetate gives nail polish remover its fruity smell. It's safer than acetone for removing polish. Wine contains small amounts naturally, contributing to its aroma. The 'new car smell' partly comes from ethyl acetate evaporating from plastics." },
    { buildAcetonitrile, "Acetonitrile", CAT_ORGANIC, "Polar solvent",
      "Acetonitrile is a versatile solvent in chemistry labs, dissolving both polar and nonpolar substances. It's used in pharmaceutical manufacturing and DNA sequencing. Though less toxic than many solvents, it's metabolized to cyanide in the body." },
    { buildDMSO, "DMSO", CAT_ORGANIC, "Dimethyl sulfoxide",
      "DMSO penetrates skin remarkably fast, carrying dissolved substances into your bloodstream. This makes it useful for drug delivery but dangerous with toxins. You can taste garlic within minutes of putting DMSO on your skin. It's also used to preserve frozen cells." },
    { buildDichloromethane, "CH2Cl2", CAT_ORGANIC, "Dichloromethane",
      "Dichloromethane is a powerful solvent used to decaffeinate coffee and strip paint. Unlike older chlorinated solvents, it's not flammable. It evaporates rapidly, which can cause frostbite. Modern decaffeination uses supercritical CO2 instead." },
    { buildChlorobenzene, "Chlorobenzene", CAT_ORGANIC, "Aromatic halide",
      "Chlorobenzene is used as a solvent and intermediate in making other chemicals. It was historically used to make DDT pesticide. With a faint almond-like odor, it's less toxic than many chlorinated compounds but still requires careful handling." },
    { buildNitrobenzene, "Nitrobenzene", CAT_ORGANIC, "Almond odor",
      "Nitrobenzene smells like bitter almonds or old shoe polish. It's used to make aniline dyes and is highly toxic, absorbed through skin. Despite its pleasant smell, even small exposures can cause headaches, nausea, and blue discoloration of the skin." },
    { buildAniline, "Aniline", CAT_ORGANIC, "Dye precursor",
      "Aniline is the parent molecule for synthetic dyes that revolutionized the textile industry in 1856. Mauveine, the first synthetic dye, was made from aniline. It's also used to make polyurethane foam, rubber chemicals, and pharmaceuticals." },
    { buildStyrene, "Styrene", CAT_PLASTICS, "Polystyrene monomer",
      "Styrene polymerizes into polystyrene, used in foam cups, packing peanuts, and insulation. You encounter styrene daily in plastic containers and packaging. It naturally occurs in cinnamon, coffee, and strawberries, though industrial exposure can be harmful." },
    { buildBenzoicAcid, "Benzoic Acid", CAT_ACIDS, "Preservative",
      "Benzoic acid and its sodium salt prevent bacterial growth in foods and drinks. Found naturally in cranberries and cinnamon, it's one of the oldest food preservatives. Your body converts it to hippuric acid, safely eliminated in urine." },
    { buildValine, "Valine", CAT_AMINO_ACIDS, "Branched-chain AA",
      "Valine is a branched-chain amino acid essential for muscle growth and repair. Your body can't make it, so you must get it from food like eggs, meat, and dairy. Athletes use valine supplements to reduce muscle breakdown during exercise." },
    { buildLeucine, "Leucine", CAT_AMINO_ACIDS, "Essential AA",
      "Leucine is the most important amino acid for building muscle, triggering protein synthesis. It's abundant in eggs, chicken, and soybeans. Bodybuilders often supplement leucine to maximize muscle growth. Your body cannot manufacture it." },
    { buildEthylene, "Ethylene", CAT_ORGANIC, "Plant hormone",
      "Ethylene is the simplest plant hormone, triggering fruit ripening. A single rotten apple spoils the bunch because it releases ethylene. Commercially, ethylene is used to ripen bananas picked green for shipping. It's also the most produced organic chemical globally." },
    { buildHydrogenSulfide, "H2S", CAT_SIMPLE, "Rotten egg smell",
      "Hydrogen sulfide causes the smell of rotten eggs and sewers. At low concentrations, it's detectable; at high concentrations, it deadens your sense of smell and becomes lethal. Your gut bacteria produce small amounts, and it may have signaling roles in the body." },
    { buildChloroform, "CHCl3", CAT_ORGANIC, "Chloroform",
      "Chloroform was once used as an anesthetic but was abandoned due to liver toxicity. In movies, a chloroform-soaked rag knocks people out instantly, but reality requires minutes of exposure. Small amounts form naturally in swimming pools and tap water." },
    // === BATCH 2: More organics (64-101) ===
    { buildTertButanol, "tert-Butanol", CAT_ORGANIC, "Tertiary alcohol",
      "Tert-butanol is unique among alcohols because it's a solid at room temperature. Used as a solvent and octane booster in gasoline, it's also a precursor for making MTBE fuel additive. Its branched structure prevents molecules from packing tightly." },
    { buildButanol, "1-Butanol", CAT_ORGANIC, "Butyl alcohol",
      "Butanol is being developed as a biofuel that could replace gasoline more directly than ethanol. It has higher energy density and doesn't absorb water like ethanol. Bacteria can produce butanol from plant sugars through fermentation." },
    { buildDiethylEther, "Diethyl Ether", CAT_ORGANIC, "Classic anesthetic",
      "Diethyl ether was the first widely used surgical anesthetic, revolutionizing medicine in 1846. Before ether, patients endured surgery awake and restrained. It's highly flammable and forms explosive peroxides when stored, making it dangerous in modern labs." },
    { buildMTBE, "MTBE", CAT_ORGANIC, "Fuel additive",
      "MTBE was added to gasoline to boost octane and reduce air pollution. However, it contaminated groundwater across America because it dissolves easily in water and doesn't break down. Most states have now banned MTBE, replaced by ethanol." },
    { buildTHF, "THF", CAT_ORGANIC, "Tetrahydrofuran",
      "THF is a common lab solvent that dissolves almost everything. It's used to make PVC pipe cement and as a precursor for spandex fabric. THF can form explosive peroxides over time, requiring careful storage and regular testing." },
    { buildDioxane, "1,4-Dioxane", CAT_ORGANIC, "Cyclic ether",
      "Dioxane is a stabilizer in chlorinated solvents and a contaminant in some cosmetics. It's classified as a probable carcinogen but persists in groundwater because it's hard to remove. Modern manufacturing tries to minimize dioxane contamination." },
    { buildDMF, "DMF", CAT_ORGANIC, "Dimethylformamide",
      "DMF is a powerful solvent used in making synthetic leather and pharmaceuticals. It can penetrate rubber gloves and skin, carrying dissolved substances into the body. Exposure is linked to liver damage, requiring strict safety protocols." },
    { buildCarbonTetrachloride, "CCl4", CAT_ORGANIC, "Carbon tet",
      "Carbon tetrachloride was once in fire extinguishers and dry cleaning fluid. It's now banned for most uses because it destroys the ozone layer and causes liver cancer. The four chlorine atoms make it completely non-flammable but very toxic." },
    { buildMethylAcetate, "Methyl Acetate", CAT_ORGANIC, "Ester solvent",
      "Methyl acetate is a fast-evaporating solvent in glues and nail polish removers. It's less toxic than many alternatives and biodegrades quickly. The fruity smell of Juicy Fruit gum comes partly from methyl acetate and similar esters." },
    { buildAceticAnhydride, "Acetic Anhydride", CAT_ORGANIC, "Acetylating agent",
      "Acetic anhydride is used to make aspirin, acetate fibers, and heroin (which is why it's controlled). It reacts with water to form acetic acid, releasing heat. The pharmaceutical industry uses tons annually for legitimate drug manufacturing." },
    { buildPropionicAcid, "Propionic Acid", CAT_ACIDS, "Preservative",
      "Propionic acid prevents mold in bread and animal feed. It's named from Greek 'protos pion' (first fat) because it's the smallest fatty acid. Your gut bacteria produce propionic acid, which may have health benefits for metabolism." },
    { buildButyricAcid, "Butyric Acid", CAT_ACIDS, "Rancid butter",
      "Butyric acid gives rancid butter and vomit their distinctive awful smell. Paradoxically, it's important for gut health. Your intestinal bacteria produce butyric acid from fiber, feeding the cells lining your colon and reducing inflammation." },
    { buildSuccinicAcid, "Succinic Acid", CAT_ACIDS, "Krebs cycle",
      "Succinic acid is a key intermediate in the Krebs cycle, your cells' main energy-producing pathway. Found in amber (hence called 'succinum' in Latin), it's used as a food additive and in biodegradable plastics. Your mitochondria produce it constantly." },
    { buildBenzaldehyde, "Benzaldehyde", CAT_FLAVORS, "Almond flavor",
      "Benzaldehyde gives almonds and cherries their characteristic smell. It's the main component of artificial almond extract used in baking. Bitter almonds contain benzaldehyde bound to cyanide, which is why they're processed before eating." },
    { buildBromobenzene, "Bromobenzene", CAT_ORGANIC, "Aromatic halide",
      "Bromobenzene is used in organic synthesis as a starting material for many compounds. The bromine atom can be replaced with other groups through reactions. It's denser than water and has a pleasant aromatic odor despite being harmful." },
    { buildPXylene, "p-Xylene", CAT_ORGANIC, "Para-xylene",
      "Para-xylene is the key ingredient for making PET plastic bottles and polyester fabric. It's extracted from petroleum and produced in enormous quantities. The global demand for PET bottles drives continuous expansion of p-xylene production." },
    { buildAnisole, "Anisole", CAT_ORGANIC, "Methoxybenzene",
      "Anisole has a sweet, pleasant smell used in perfumes. It's the methyl ether of phenol and serves as a precursor for many flavorings and fragrances. The anise plant contains compounds similar to anisole, giving it a licorice-like aroma." },
    { buildPhenylacetylene, "Phenylacetylene", CAT_ORGANIC, "Aromatic alkyne",
      "Phenylacetylene combines a benzene ring with a reactive triple bond. It's used in organic synthesis and can undergo click chemistry reactions important in drug development. The triple bond makes it useful for building complex molecules." },
    { buildFructose, "Fructose", CAT_SUGARS, "Fruit sugar",
      "Fructose is the sweetest natural sugar, found in fruits and honey. Unlike glucose, fructose is processed mainly by your liver. High-fructose corn syrup in sodas has been linked to obesity because the liver converts excess fructose to fat." },
    { buildRibose, "Ribose", CAT_SUGARS, "RNA sugar",
      "Ribose is the sugar in RNA and ATP, your cells' energy currency. It has one more oxygen than deoxyribose in DNA. Athletes sometimes supplement ribose to boost energy, though evidence for benefits is mixed." },
    { buildDeoxyribose, "Deoxyribose", CAT_SUGARS, "DNA sugar",
      "Deoxyribose forms DNA's backbone, lacking one oxygen atom compared to ribose. This small difference makes DNA more stable than RNA, perfect for long-term genetic storage. Every cell in your body contains meters of DNA wrapped around deoxyribose sugars." },
    { buildIsoleucine, "Isoleucine", CAT_AMINO_ACIDS, "Branched-chain AA",
      "Isoleucine is a branched-chain amino acid essential for muscle metabolism and immune function. Your body can't synthesize it, requiring dietary sources like eggs, chicken, and lentils. It helps regulate blood sugar and energy levels during exercise." },
    { buildSerine, "Serine", CAT_AMINO_ACIDS, "Polar amino acid",
      "Serine is crucial for making proteins, cell membranes, and DNA. It can be modified by adding phosphate groups, which acts as an on/off switch for proteins. Silk is rich in serine, giving it unique properties." },
    { buildThreonine, "Threonine", CAT_AMINO_ACIDS, "Essential AA",
      "Threonine was the last of the 20 common amino acids to be discovered, in 1936. It's essential for creating collagen and elastin in skin. Your intestinal lining needs threonine to maintain its protective mucus layer." },
    { buildAsparticAcid, "Aspartic Acid", CAT_AMINO_ACIDS, "Acidic AA",
      "Aspartic acid is a non-essential amino acid your body can make. It plays roles in the urea cycle and neurotransmission. Combined with phenylalanine, it forms aspartame artificial sweetener. It's abundant in asparagus, hence its name." },
    { buildGlutamicAcid, "Glutamic Acid", CAT_AMINO_ACIDS, "MSG precursor",
      "Glutamic acid is the most abundant amino acid in the brain, serving as an excitatory neurotransmitter. Its sodium salt is MSG, the flavor enhancer giving umami taste. Parmesan cheese and tomatoes are naturally rich in glutamate." },
    { buildLysine, "Lysine", CAT_AMINO_ACIDS, "Basic AA",
      "Lysine is essential for producing collagen and absorbing calcium. It may help prevent cold sores by blocking arginine, which herpes virus needs. Vegetarians should ensure adequate lysine, as grains are low in this amino acid." },
    { buildHistidine, "Histidine", CAT_AMINO_ACIDS, "Aromatic AA",
      "Histidine is converted to histamine, the compound causing allergic reactions. It's essential for children's growth though adults can synthesize small amounts. The imidazole ring in histidine helps proteins bind metals and transfer protons." },
    { buildPhenylalanine, "Phenylalanine", CAT_AMINO_ACIDS, "Aromatic AA",
      "Phenylalanine is an essential amino acid your body converts to tyrosine and then to dopamine. People with phenylketonuria (PKU) can't process it safely, which is why diet sodas with aspartame carry warnings. It's abundant in meat and dairy." },
    { buildTyrosine, "Tyrosine", CAT_AMINO_ACIDS, "Aromatic AA",
      "Tyrosine is the starting material for dopamine, adrenaline, and thyroid hormones. It's also responsible for melanin pigment in skin and hair. Under stress, tyrosine supplements may help maintain mental performance by supporting neurotransmitter production." },
    { buildTryptophan, "Tryptophan", CAT_AMINO_ACIDS, "Aromatic AA",
      "Tryptophan is converted to serotonin and melatonin, regulating mood and sleep. Despite the myth, turkey doesn't make you sleepy from tryptophan; carbs and overeating are the real culprits. It's the rarest amino acid in most proteins." },
    { buildProline, "Proline", CAT_AMINO_ACIDS, "Cyclic AA",
      "Proline is unique because its side chain loops back to the backbone, creating kinks in proteins. It makes up about 10% of collagen, giving your skin and joints their structure. Gelatin desserts get their texture from proline-rich collagen." },
    { buildCysteine, "Cysteine", CAT_AMINO_ACIDS, "Sulfur AA",
      "Cysteine contains sulfur that forms disulfide bridges, holding proteins together. Your hair's strength comes from cysteine crosslinks. It's also part of glutathione, your body's main antioxidant. Permanent waves break and reform these sulfur bonds." },
    { buildMethionine, "Methionine", CAT_AMINO_ACIDS, "Sulfur AA",
      "Methionine is the start codon amino acid, beginning every protein your cells make. It's essential and found in eggs, fish, and Brazil nuts. Methionine donates methyl groups for DNA methylation, an important gene regulation mechanism." },
    { buildPyruvate, "Pyruvate", CAT_ENERGY, "Glycolysis product",
      "Pyruvate is where glycolysis ends and either fermentation or cellular respiration begins. Your muscles produce lactate from pyruvate during intense exercise. It sits at a metabolic crossroads, convertible to glucose, fat, or energy." },
    { buildArginine, "Arginine", CAT_AMINO_ACIDS, "Basic AA",
      "Arginine is converted to nitric oxide, relaxing blood vessels. It's conditionally essential, meaning you need more during growth or illness. Bodybuilders take arginine supplements hoping to boost muscle blood flow, though benefits are debated." },
    { buildAsparagine, "Asparagine", CAT_AMINO_ACIDS, "Polar AA",
      "Asparagine was the first amino acid discovered, isolated from asparagus in 1806. It's non-essential as your body can make it. When you cook high-asparagine foods at high heat with sugars, acrylamide forms, which is a potential carcinogen." },
    { buildGlutamine, "Glutamine", CAT_AMINO_ACIDS, "Polar AA",
      "Glutamine is the most abundant amino acid in your blood and muscles. Your immune cells use it as fuel, especially during illness or stress. Hospitals give glutamine to burn patients and those recovering from surgery to speed healing." },
    { buildHydroxyproline, "Hydroxyproline", CAT_AMINO_ACIDS, "Collagen AA",
      "Hydroxyproline is almost exclusively found in collagen, making up about 10% of it. Your body makes it by modifying proline after the protein is built. Vitamin C is essential for this modification, which is why scurvy victims have weak collagen." },
    // === VITAMINS (103-116) ===
    { buildAscorbicAcid, "Vitamin C", CAT_VITAMINS, "Ascorbic acid",
      "Vitamin C is essential for making collagen, healing wounds, and immune function. Humans and guinea pigs are among the few mammals that can't synthesize it. Sailors discovered citrus fruits cured scurvy, though they didn't understand why until the 1930s." },
    { buildThiamine, "Vitamin B1", CAT_VITAMINS, "Thiamine",
      "Thiamine deficiency caused beriberi, a devastating disease in populations eating polished white rice. The discovery that rice bran cured beriberi led to the vitamin concept. Thiamine is crucial for nerve function and carbohydrate metabolism." },
    { buildRiboflavin, "Vitamin B2", CAT_VITAMINS, "Riboflavin",
      "Riboflavin is the vitamin that turns your urine bright yellow after taking supplements. It's essential for energy metabolism and maintaining healthy skin. Light destroys riboflavin, which is why milk comes in opaque containers." },
    { buildNiacin, "Vitamin B3", CAT_VITAMINS, "Niacin",
      "Niacin deficiency causes pellagra, characterized by the 'four Ds': dermatitis, diarrhea, dementia, and death. High doses cause harmless but uncomfortable flushing. Niacin is part of NAD, essential for hundreds of metabolic reactions." },
    { buildPanthothenicAcid, "Vitamin B5", CAT_VITAMINS, "Pantothenic acid",
      "Pantothenic acid is found in almost all foods, hence 'pantos' meaning everywhere. It's part of coenzyme A, crucial for making and breaking down fats. Deficiency is extremely rare because it's so widespread in the diet." },
    { buildPyridoxine, "Vitamin B6", CAT_VITAMINS, "Pyridoxine",
      "Vitamin B6 is involved in over 100 enzyme reactions, mostly in protein metabolism. It helps make neurotransmitters like serotonin and dopamine. Deficiency can cause confusion and depression due to disrupted brain chemistry." },
    { buildBiotin, "Vitamin B7", CAT_VITAMINS, "Biotin",
      "Biotin is marketed for hair and nail health, though deficiency is rare. Eating raw egg whites can cause biotin deficiency because avidin protein binds biotin. Your gut bacteria produce biotin, contributing to your daily needs." },
    { buildFolicAcid, "Vitamin B9", CAT_VITAMINS, "Folic acid",
      "Folic acid prevents neural tube defects like spina bifida in developing babies. Taking it before and during early pregnancy reduces risk by 70%. Many countries now fortify flour with folic acid because defects occur before women know they're pregnant." },
    { buildRetinol, "Vitamin A", CAT_VITAMINS, "Retinol",
      "Vitamin A is essential for vision, especially night vision. Ancient Egyptians treated night blindness with liver, now known to be rich in vitamin A. Too much vitamin A is toxic; polar bear liver has killed Arctic explorers." },
    { buildBetaCarotene, "Beta-Carotene", CAT_VITAMINS, "Provitamin A",
      "Beta-carotene gives carrots and sweet potatoes their orange color. Your body converts it to vitamin A as needed, making overdose nearly impossible. Eating lots of carrots can turn your skin orange, a harmless condition called carotenemia." },
    { buildCholecalciferol, "Vitamin D3", CAT_VITAMINS, "Cholecalciferol",
      "Vitamin D3 is made in your skin when exposed to sunlight, technically making it a hormone. It's essential for calcium absorption and bone health. Many people are deficient, especially those in northern latitudes or who stay indoors." },
    { buildAlphaTocopherol, "Vitamin E", CAT_VITAMINS, "Alpha-tocopherol",
      "Vitamin E is a powerful antioxidant protecting cell membranes from damage. It's abundant in nuts, seeds, and vegetable oils. Despite supplement marketing, clinical trials haven't shown benefits for heart disease or cancer prevention." },
    { buildPhylloquinone, "Vitamin K1", CAT_VITAMINS, "Phylloquinone",
      "Vitamin K1 is essential for blood clotting, named from the German 'Koagulation.' Green leafy vegetables are the best source. People on blood thinners must maintain consistent vitamin K intake because it can interfere with their medication." },
    { buildNicotinamide, "Nicotinamide", CAT_VITAMINS, "B3 amide form",
      "Nicotinamide is a form of vitamin B3 that doesn't cause the flushing associated with niacin. It's being studied for preventing skin cancer and improving skin appearance. It's part of NAD+, a molecule important for cellular energy and longevity research." },
    { buildNiacinamide, "Niacinamide (Alt. B3)", CAT_VITAMINS, "Vitamin B3 variant",
      "Niacinamide is another name for nicotinamide, vitamin B3's amide form. Skincare products tout it for reducing pore size, redness, and supporting the skin barrier. Unlike niacin, it doesn't cause flushing, making it ideal for sensitive skin formulations." },
    // === CONTROLLED SUBSTANCES & RESEARCH COMPOUNDS (117-123) ===
    { buildCocaine, "Cocaine", CAT_CONTROLLED, "Stimulant alkaloid",
      "Cocaine comes from coca leaves, chewed in South America for thousands of years. It blocks dopamine reuptake, causing intense euphoria and addiction. Once used in Coca-Cola and as a local anesthetic, it remains medically useful for certain eye and nose surgeries." },
    { buildHeroin, "Heroin", CAT_CONTROLLED, "Opioid",
      "Heroin was originally marketed by Bayer in 1898 as a 'non-addictive' cough suppressant. It's diacetylmorphine, crossing the blood-brain barrier faster than morphine. The opioid crisis has made heroin one of the deadliest drugs, often contaminated with fentanyl." },
    { buildFentanyl, "Fentanyl/Sublimaze", CAT_CONTROLLED, "Synthetic opioid",
      "Fentanyl is 100 times more potent than morphine, used medically in patches and lollipops for severe pain. Illegally manufactured fentanyl has caused a dramatic rise in overdose deaths. Just 2 milligrams, smaller than a few grains of salt, can be lethal." },
    { buildPropofol, "Propofol/Diprivan", CAT_CONTROLLED, "IV anesthetic",
      "Propofol is the white 'milk of amnesia' used for anesthesia induction worldwide. It works within seconds and wears off quickly, ideal for short procedures. Michael Jackson's death from propofol misuse highlighted the dangers of using it outside medical settings." },
    { buildTHC, "THC", CAT_CONTROLLED, "Cannabis active",
      "THC is the main psychoactive compound in cannabis, binding to cannabinoid receptors in the brain. Your body makes similar molecules called endocannabinoids. Medical uses include treating nausea from chemotherapy and stimulating appetite in AIDS patients." },
    { buildCreatine, "Creatine", CAT_ENERGY, "Muscle energy",
      "Creatine stores high-energy phosphate in muscles, providing quick energy for short bursts of activity. It's one of the most studied and safest sports supplements, improving strength and power. Your body makes about 1 gram daily, with more from meat." },
    { buildOctane, "Octane", CAT_ORGANIC, "Gasoline component",
      "Octane is the reference standard for fuel ratings. Higher octane fuels resist premature ignition, preventing engine 'knock.' Iso-octane has a rating of 100, while n-heptane is rated 0. Most cars run fine on 87 octane; premium is often unnecessary." },
    // === NSAIDS & STATINS (123-130) ===
    { buildSimvastatin, "Simvastatin/Zocor", CAT_PHARMA, "Cholesterol drug",
      "Simvastatin lowers cholesterol by blocking an enzyme in your liver. Doctors prescribe it to reduce heart attack and stroke risk. It's one of the most prescribed medications worldwide." },
    { buildIbuprofen, "Ibuprofen/Advil", CAT_PHARMA, "NSAID",
      "Ibuprofen blocks pain-causing prostaglandins, reducing pain, fever, and inflammation. It works within 30 minutes and is safe for most people, making it the world's most popular over-the-counter painkiller." },
    { buildNaproxen, "Naproxen/Aleve", CAT_PHARMA, "NSAID",
      "Naproxen is a longer-acting pain reliever than ibuprofen; one dose lasts 12 hours. Athletes favor it for joint pain, and it's branded as Aleve. Fewer doses needed daily makes it convenient." },
    { buildDiclofenac, "Diclofenac/Voltaren", CAT_PHARMA, "NSAID",
      "Diclofenac is a potent NSAID stronger than ibuprofen for severe pain and arthritis. It reduces inflammation effectively but has higher gastrointestinal side effects. Many countries restrict it for heart safety." },
    { buildIndomethacin, "Indomethacin/Indocin", CAT_PHARMA, "NSAID",
      "Indomethacin is one of the most powerful NSAIDs, used for severe arthritis, migraines, and gout. Its strength makes it effective but requires careful dosing. It's often reserved for cases where other NSAIDs fail." },
    { buildCelecoxib, "Celecoxib/Celebrex", CAT_PHARMA, "COX-2 inhibitor",
      "Celecoxib targets COX-2 specifically, reducing pain with fewer stomach problems than traditional NSAIDs. Developed as a safer alternative, it still carries cardiovascular risks. Patients appreciate the reduced gastrointestinal side effects." },
    { buildMeloxicam, "Meloxicam/Mobic", CAT_PHARMA, "NSAID",
      "Meloxicam is a once-daily NSAID for arthritis patients. Its long half-life means steady pain relief without multiple doses. It's particularly popular for chronic pain management in seniors." },
    { buildAcetaminophen, "Tylenol", CAT_PHARMA, "Paracetamol",
      "Acetaminophen (Tylenol) reduces pain and fever but doesn't reduce inflammation like NSAIDs. It's gentler on stomachs and safe during pregnancy. Overdose damages the liver, requiring careful dosing." },
    // === STEROID HORMONES (131-140) ===
    { buildTestosterone, "Testosterone", CAT_HORMONES, "Male sex hormone",
      "Testosterone is the primary male sex hormone driving masculine traits, muscle growth, and aggression. Females produce small amounts; males produce 10 times more. It's essential for bone density and libido." },
    { buildDHT, "DHT", CAT_HORMONES, "Dihydrotestosterone",
      "Dihydrotestosterone (DHT) is testosterone's more potent form, causing male pattern baldness. Blocking DHT prevents baldness but can reduce some masculine characteristics. It's why many baldness treatments target DHT production." },
    { buildAndrostenedione, "Androstenedione", CAT_HORMONES, "Androgen precursor",
      "Androstenedione is a male hormone precursor that your body converts to testosterone. Athletes once used it as a supplement, but it was banned from sports competition. Your adrenal glands naturally produce small amounts daily." },
    { buildEstradiol, "Estradiol/E2", CAT_HORMONES, "Primary estrogen",
      "Estradiol is the most potent estrogen, responsible for female sexual traits and reproduction. Women produce about 100 micrograms daily; men produce trace amounts. Birth control pills often contain synthetic estradiol to prevent ovulation." },
    { buildEstrone, "Estrone/E1", CAT_HORMONES, "Postmenopausal estrogen",
      "Estrone is a weaker form of estrogen produced mainly after menopause. It's created in fat tissue, which is why postmenopausal weight gain increases estrone levels. Unlike estradiol, estrone remains more stable throughout the day." },
    { buildEstriol, "Estriol/E3", CAT_HORMONES, "Pregnancy estrogen",
      "Estriol is the main estrogen during pregnancy, produced in enormous quantities by the placenta. It protects the developing baby and prepares the body for delivery. Estriol levels drop dramatically after birth, causing postpartum hormone shifts." },
    { buildProgesterone, "Progesterone", CAT_HORMONES, "Pregnancy hormone",
      "Progesterone is the 'pregnancy hormone' preparing the uterus to receive a fertilized egg. It maintains pregnancy and prevents premature labor. Levels rise after ovulation, signaling the body it's ready for pregnancy." },
    { buildCortisol, "Cortisol", CAT_HORMONES, "Stress hormone",
      "Cortisol is your stress hormone, rising when you're anxious or scared. It mobilizes glucose for fight-or-flight response but chronic stress causes health problems. Cortisol peaks in early morning, helping you wake up naturally." },
    { buildCortisone, "Cortisone", CAT_HORMONES, "Cortisol precursor",
      "Cortisone is cortisol's inactive form, stored in tissue until your body needs it. Doctors prescribe cortisone shots for joint pain and inflammation. It's different from cortisol because it doesn't act until your body converts it." },
    { buildAldosterone, "Aldosterone", CAT_HORMONES, "Salt balance",
      "Aldosterone controls salt and water balance, regulating blood pressure. Your kidneys release it when sodium levels drop, making you retain water. This is why too much salt raises blood pressure over time." },
    // === HOUSEHOLD CHEMICALS (141-150) ===
    { buildBleach, "Bleach/NaOCl", CAT_HOUSEHOLD, "Disinfectant",
      "Bleach (sodium hypochlorite) disinfects surfaces and water by releasing dangerous free radicals that destroy bacteria and viruses. It's toxic if swallowed and creates poisonous gas when mixed with ammonia. Diluted properly, it's safe for sanitizing." },
    { buildLye, "Lye/NaOH", CAT_HOUSEHOLD, "Caustic soda",
      "Lye (sodium hydroxide) is one of the strongest bases, used to unclog drains and make soap. It can cause severe chemical burns, completely dissolving organic material. Never touch lye without gloves; it liquefies skin on contact." },
    { buildTableSalt, "Table Salt/NaCl", CAT_HOUSEHOLD, "Sodium chloride",
      "Table salt (sodium chloride) seasons food and is essential for nerve function and blood pressure regulation. Your body needs about 500 mg daily, but Americans consume 10 times that amount. Too much salt increases heart disease and stroke risk." },
    { buildBakingSoda, "Baking Soda/NaHCO3", CAT_HOUSEHOLD, "Sodium bicarbonate",
      "Baking soda (sodium bicarbonate) releases carbon dioxide when heated, making cakes and cookies fluffy. It's mildly alkaline, so it neutralizes acid reflux. Many people use it as a natural deodorant and cleaning paste." },
    { buildWashingSoda, "Washing Soda/Na2CO3", CAT_HOUSEHOLD, "Sodium carbonate",
      "Washing soda (sodium carbonate) is stronger than baking soda and boosts laundry detergent power. It breaks down grease and stains more effectively than regular soap. Mixed with baking soda, it makes an excellent gentle cleaner." },
    { buildSodiumFluoride, "Sodium Fluoride/NaF", CAT_HOUSEHOLD, "Toothpaste",
      "Sodium fluoride strengthens tooth enamel, preventing cavities by replacing hydroxyl with fluoride. It's added to most toothpastes and drinking water in many cities. Some people worry about fluoride despite safety evidence, calling it controversial." },
    { buildSodiumNitrate, "Sodium Nitrate/NaNO3", CAT_HOUSEHOLD, "Preservative",
      "Sodium nitrate preserves cured meats like bacon and hot dogs by killing harmful bacteria. It's been used for centuries, long before refrigeration. High consumption links to increased cancer risk, leading to health debates." },
    { buildSodiumNitrite, "Sodium Nitrite/NaNO2", CAT_HOUSEHOLD, "Meat curing",
      "Sodium nitrite is more potent than nitrate for meat preservation and faster-acting. It also prevents clostridium botulinum toxin, potentially deadly if untreated. Modern curing combines it with nitrate for safety and effectiveness." },
    { buildSodiumSulfate, "Sodium Sulfate/Na2SO4", CAT_HOUSEHOLD, "Detergent filler",
      "Sodium sulfate fills detergents to increase volume cheaply. It doesn't help clean but makes products easier to handle. It's also used in laxatives and glass manufacturing as an economical filler." },
    { buildMSG, "MSG", CAT_HOUSEHOLD, "Flavor enhancer",
      "MSG (monosodium glutamate) enhances savory flavors, making food taste more delicious. Your tongue has specific receptors detecting umami, the fifth basic taste. It's naturally abundant in tomatoes, parmesan cheese, and aged foods." },
    { buildGlycerin, "Glycerin", CAT_HOUSEHOLD, "Humectant moisturizer",
      "Glycerin (glycerol) is a sweet liquid that absorbs water from the air, making it perfect for moisturizers and cosmetics. It's a byproduct of biodiesel production and soap making. Glycerin lubricates skin, reduces wrinkles' appearance, and is used in suppositories." },
    // === ACIDS (152-167) ===
    { buildCitricAcid, "Citric Acid", CAT_ACIDS, "Citrus fruits",
      "Citric acid gives lemons and limes their sour taste. It's used in foods, drinks, and cleaning products for its acidity. Chemically produced citric acid is indistinguishable from natural sources and costs less." },
    { buildCarbonicAcid, "Carbonic Acid/H2CO3", CAT_ACIDS, "CO2 in water",
      "Carbonic acid forms when carbon dioxide dissolves in water, making fizzy drinks bubbly. It's extremely weak, so soda doesn't harm teeth directly. However, sugar in soda causes decay; the acid may increase risk." },
    { buildBoricAcid, "Boric Acid/H3BO3", CAT_ACIDS, "Antiseptic",
      "Boric acid is a weak acid used as an antiseptic in eye drops and mouthwash. It's toxic if ingested in large quantities, especially for babies. Borax, its salt form, is used in detergents and glass making." },
    { buildOxalicAcid, "Oxalic Acid", CAT_ACIDS, "Rust remover",
      "Oxalic acid dissolves rust and mineral deposits, found naturally in spinach and rhubarb. Eating too much can interfere with calcium absorption, which is why oxalate reduction is recommended for kidney stone patients. It creates insoluble salts with calcium." },
    { buildTartaricAcid, "Tartaric Acid", CAT_ACIDS, "Wine acid",
      "Tartaric acid gives wine its tart taste and occurs naturally in grapes. It precipitates as cream of tartar (potassium bitartrate) during wine aging. Bakers use cream of tartar as a leavening agent in baked goods alongside baking soda." },
    { buildMalicAcid, "Malic Acid", CAT_ACIDS, "Apple acid",
      "Malic acid makes apples and tart candies sour, providing that pucker-inducing sensation. Your body produces it during cellular respiration as an energy intermediate. Some athletes use malic acid supplements, believing it reduces fatigue, though evidence is mixed." },
    { buildGlycolicAcid, "Glycolic Acid", CAT_ACIDS, "AHA exfoliant",
      "Glycolic acid is the smallest AHA (alpha hydroxy acid), derived from sugar cane, making it highly effective at exfoliating skin. Used in skincare products to improve skin texture, reduce acne, and minimize wrinkles by removing dead skin cells. It's the mildest AHA but very effective, penetrating skin well due to its small molecular size." },
    { buildSalicylicAcid, "Salicylic Acid", CAT_ACIDS, "BHA acne fighter",
      "Salicylic acid is the primary BHA (beta hydroxy acid) used in skincare, derived from willow bark and wintergreen. Unlike AHAs, it's lipophilic (fat-soluble), allowing it to penetrate oil-clogged pores to treat acne. Found in most acne spot treatments and blackhead-fighting products worldwide." },
    { buildHydrofluoricAcid, "Hydrofluoric Acid/HF", CAT_ACIDS, "Glass etching",
      "Hydrofluoric acid (HF) dissolves glass and rock, making it unique among common acids. Despite being weak compared to hydrochloric acid, it's extremely dangerous because fluoride ions penetrate tissue deeply. A small skin exposure can cause fatal heart arrhythmias." },
    { buildHydrobromicAcid, "Hydrobromic Acid/HBr", CAT_ACIDS, "Strong acid",
      "Hydrobromic acid is a strong acid used in pharmaceuticals and organic synthesis. It's less common than hydrochloric acid but equally corrosive. Bromide salts have been used historically as sedatives, though barbiturates replaced them." },
    { buildHydroiodicAcid, "Hydroiodic Acid/HI", CAT_ACIDS, "Strong acid",
      "Hydroiodic acid (HI) is a strong, colorless acid used in laboratory and industrial synthesis. It's unstable, slowly decomposing to iodine and water. Iodide salts are important in nuclear medicine and photography." },
    { buildPerchloricAcid, "Perchloric Acid/HClO4", CAT_ACIDS, "Strongest acid",
      "Perchloric acid is the strongest of common acids, extremely corrosive and potentially explosive. It's used as a catalyst in explosives and rocket fuel production. Perchlorate salts are found in some food supplies from environmental contamination." },
    { buildAcrylicAcid, "Acrylic Acid", CAT_ACIDS, "Polymer precursor",
      "Acrylic acid polymerizes into polyacrylate, used in super glue and diapers' absorbent gels. Acrylic acid itself is corrosive and irritating. It's produced on an enormous industrial scale for various plastics." },
    // === PLASTICS (162-179) ===
    { buildVinylChloride, "Vinyl Chloride/PVC", CAT_PLASTICS, "PVC monomer",
      "Vinyl chloride polymerizes into polyvinyl chloride (PVC), the third most common plastic. Workers exposed to vinyl chloride dust develop angiosarcoma liver cancer. PVC is now made more safely, but the hazard remains with improper handling." },
    { buildMethylMethacrylate, "MMA/Plexiglas", CAT_PLASTICS, "Acrylic monomer",
      "Methyl methacrylate polymerizes into acrylic glass (Plexiglas), stronger and lighter than regular glass. It's used in aquariums, aircraft windows, and dental prosthetics. Acrylic sheets are shatter-resistant and easily molded when heated." },
    { buildTetrafluoroethylene, "TFE/Teflon", CAT_PLASTICS, "PTFE monomer",
      "TFE polymerizes into PTFE (Teflon), the non-stick coating on frying pans. It's heat-resistant to 260°C and chemically inert. Teflon was discovered accidentally when a refrigerant experiment produced this slippery powder." },
    { buildCaprolactam, "Caprolactam/Nylon-6", CAT_PLASTICS, "Nylon-6 monomer",
      "Caprolactam polymerizes into nylon-6, used in clothing fibers, carpets, and car parts. Nylon revolutionized textiles in the 1930s, replacing silk and natural fibers. Its strength and elasticity make it ideal for stockings and parachutes." },
    { buildAcrylonitrile, "Acrylonitrile/ABS", CAT_PLASTICS, "ABS component",
      "Acrylonitrile is a component of ABS plastic, found in Lego bricks and automotive dashboards. It's carcinogenic and volatile, requiring careful handling in factories. ABS combines hardness with flexibility for durable consumer products." },
    { buildButadiene, "1,3-Butadiene/ABS", CAT_PLASTICS, "Rubber/ABS",
      "Butadiene polymerizes into synthetic rubber, replacing natural rubber for car tires and hoses. World War II accelerated synthetic rubber development when natural rubber supplies were cut off. It's also used in ABS plastic and nitrile gloves." },
    { buildBisphenolA, "Bisphenol A/BPA", CAT_PLASTICS, "Polycarbonate",
      "BPA polymerizes into polycarbonate plastic, used in water bottles and sports equipment. BPA leaches into drinks and foods, acting as an endocrine disruptor. Many countries now restrict BPA in baby bottles and sippy cups." },
    { buildTerephthalicAcid, "Terephthalic Acid/PET", CAT_PLASTICS, "PET monomer",
      "Terephthalic acid polymerizes into PET plastic, used in soda bottles and polyester clothing. PET is recyclable, labeled as #1 plastic, and the most recycled plastic worldwide. Its durability and clarity make it ideal for beverage containers." },
    { buildVinylAcetate, "Vinyl Acetate/EVA", CAT_PLASTICS, "EVA monomer",
      "Vinyl acetate polymerizes into ethylene vinyl acetate (EVA), used in shoe cushioning and flexible packaging. EVA is flexible yet durable, making comfortable shoe soles. It's also used in adhesive applications and protective coatings." },
    { buildDimethylsiloxane, "PDMS/Silicone", CAT_PLASTICS, "Silicone monomer",
      "PDMS polymerizes into silicone, used in sealants, medical implants, and cooking utensils. Silicone is biocompatible and heat-resistant, ideal for breast implants and baking pans. It's water-repellent yet breathable, used in waterproofing." },
    { buildAdipicAcid, "Adipic Acid/Nylon-6,6", CAT_PLASTICS, "Nylon monomer",
      "Adipic acid is a nylon-6,6 monomer, combined with hexamethylenediamine. Nylon-6,6 is stiffer than nylon-6, used in engineering applications. It's also used in polyester-polyol production for rigid foams and coatings." },
    { buildHexamethylenediamine, "HMDA/Nylon-6,6", CAT_PLASTICS, "Nylon monomer",
      "HMDA is nylon-6,6's second monomer, reacted with adipic acid. It's harder and stiffer than nylon-6, used in tough applications like gears and bearings. HMDA production is energy-intensive but essential for engineering plastics." },
    { buildIsoprene, "Isoprene/Rubber", CAT_PLASTICS, "Natural rubber",
      "Isoprene polymerizes into natural rubber, tapped from rubber trees. Synthetic isoprene matches natural rubber's properties, making tires and elastic products. Your body produces isoprene naturally, which is why it's on your breath." },
    { buildPEDimer, "PE Dimer/HDPE", CAT_PLASTICS, "Polyethylene",
      "Polyethylene dimers represent the basic unit of polyethylene (HDPE/LDPE), the most common plastic. PE is made from ethylene monomers; the dimer shows how they link. Everything from bags to milk jugs is made from PE." },
    { buildPPDimer, "PP Dimer", CAT_PLASTICS, "Polypropylene",
      "Polypropylene (PP) dimers represent PP plastic units, used in food containers and car parts. PP is stronger than PE and withstands higher temperatures. It's the plastic your microwave containers are labeled with." },
    { buildPVCTrimer, "PVC Trimer", CAT_PLASTICS, "PVC oligomer",
      "PVC trimers show how vinyl chloride units combine in PVC plastic. PVC is rigid, making pipes and window frames, or flexible for vinyl records and shower curtains. It requires additives (plasticizers) to be flexible." },
    { buildPSDimer, "PS Dimer", CAT_PLASTICS, "Polystyrene",
      "Polystyrene dimers represent PS plastic, used in foam cups and hard cases. PS is rigid and brittle when solid, lightweight when foamed. It's made by dissolving carbon dioxide in melted polystyrene." },
    { buildPTFETrimer, "PTFE Trimer/Teflon", CAT_PLASTICS, "Teflon oligomer",
      "PTFE trimers show how TFE units bond, representing Teflon's basic structure. PTFE is the slipperiest known material, used beyond cookware in medical devices and machinery. Its perfluoro bonds make it extremely stable." },
    // === FATTY ACIDS (180-189) ===
    { buildElaidicAcid, "Elaidic Acid/trans-9", CAT_FATS, "Industrial trans fat",
      "Elaidic acid is a trans fatty acid from industrial hydrogenation. Trans fats raise bad cholesterol and lower good cholesterol, increasing heart disease risk. Many countries now ban or limit artificial trans fats in foods." },
    { buildOleicAcid, "Oleic Acid/cis-9", CAT_FATS, "Olive oil",
      "Oleic acid is the monounsaturated fat making olive oil liquid at room temperature. It's considered healthy, abundant in Mediterranean diets. Your body can synthesize oleic acid from other fats." },
    { buildVaccenicAcid, "Vaccenic Acid/trans-11", CAT_FATS, "Natural trans fat",
      "Vaccenic acid is a natural trans fat from ruminant animals like cows and sheep. Unlike industrial trans fats, natural trans fats in small amounts don't raise disease risk as much. Grass-fed dairy has higher vaccenic acid levels." },
    { buildTransPalmitoleicAcid, "trans-Palmitoleic", CAT_FATS, "Dairy fat",
      "Trans-palmitoleic is a trans fat from dairy products, naturally occurring in small amounts. Unlike industrial trans fats, dairy trans fats don't raise heart disease risk as much. Grass-fed butter has higher trans-palmitoleic levels." },
    { buildLinoelaidicAcid, "Linoelaidic/trans,trans", CAT_FATS, "Hydrogenation",
      "Linoelaidic acid has two trans double bonds, formed during hydrogenation of vegetable oils. It's more common in processed foods, raising cholesterol problems. This is why nutritionists recommend avoiding hydrogenated oils in margarine." },
    { buildRumenicAcid, "Rumenic Acid/CLA", CAT_FATS, "Conjugated linoleic",
      "Rumenic acid is a conjugated linoleic acid (CLA) from grass-fed meat and dairy. Some studies suggest CLA aids fat loss and muscle building, though effects are modest. Grass-fed animals produce more CLA than grain-fed ones." },
    { buildStearicAcid, "Stearic Acid/C18:0", CAT_FATS, "Saturated fat",
      "Stearic acid is a saturated fat found in animal fats and cocoa butter. Unlike other saturated fats, stearic acid doesn't raise LDL cholesterol. Your body converts it to oleic acid, the healthy monounsaturated fat." },
    { buildPalmiticAcid, "Palmitic Acid/C16:0", CAT_FATS, "Palm oil",
      "Palmitic acid is the main saturated fat in palm oil, abundant and cheap. It raises cholesterol and increases heart disease risk. Palm oil production drives rainforest deforestation, making it environmentally controversial." },
    { buildLinoleicAcid, "Linoleic Acid/Omega-6", CAT_FATS, "Essential FA",
      "Linoleic acid is an essential omega-6 polyunsaturated fat you must obtain from food. Found in vegetable oils and nuts, it's vital for skin health and inflammation. Modern diets contain excessive omega-6, leading to inflammation imbalances." },
    { buildAlphaLinolenicAcid, "ALA/Omega-3", CAT_FATS, "Plant omega-3",
      "ALA is an essential omega-3 fat from plants, found in flaxseed and walnuts. Your body converts ALA to EPA and DHA, though inefficiently. Getting both omega-3 and omega-6 in balance is crucial for health." },
    // === METAL COMPOUNDS (190-198) ===
    { buildRust, "Rust/Fe2O3", CAT_METALS, "Iron oxide",
      "Rust (iron oxide) forms when iron oxidizes, causing corrosion and weakness. Rust is electrochemically complex, requiring both oxygen and moisture. Preventing rust costs billions annually through painting, galvanizing, and stainless steel." },
    { buildMagnetite, "Magnetite/Fe3O4", CAT_METALS, "Magnetic iron",
      "Magnetite is a magnetic iron oxide naturally occurring as a mineral. Iron-loving bacteria use magnetite for navigation. Magnetite is red, compared to rust's red-brown, and found in ancient sediments." },
    { buildIronPentacarbonyl, "Iron Pentacarbonyl", CAT_METALS, "Fe(CO)5",
      "Iron pentacarbonyl is a volatile iron compound used in metal plating and synthesis. It's extremely toxic and flammable, requiring extreme caution. Despite dangers, it's valuable in nanomaterial production." },
    { buildCopperSulfate, "Copper Sulfate/CuSO4", CAT_METALS, "Blue vitriol",
      "Copper sulfate is bright blue, used as a fungicide and algicide in agriculture. It's toxic to aquatic life, making it effective for water treatment. Medieval times saw copper sulfate used as a pigment for blue paints." },
    { buildVerdigris, "Verdigris/Patina", CAT_METALS, "Copper carbonate",
      "Verdigris (copper carbonate patina) forms when copper oxidizes, giving it a distinctive green color. Statue of Liberty is green because of verdigris forming on its copper surface. Renaissance artists used verdigris as a green pigment in paintings." },
    { buildFerrocene, "Ferrocene", CAT_METALS, "Sandwich compound",
      "Ferrocene is an organometallic compound with an iron atom sandwiched between two aromatic rings. Its unusual structure revolutionized organometallic chemistry. It's used as a fuel additive and in research for advanced materials." },
    { buildCisplatin, "Cisplatin", CAT_METALS, "Cancer drug",
      "Cisplatin is a platinum-based chemotherapy drug that binds DNA, killing cancer cells. It's one of the most effective cancer drugs but causes severe side effects like hearing loss. Its geometry makes it work; the trans isomer is inactive." },
    { buildAluminumOxide, "Alumina/Al2O3", CAT_METALS, "Aluminum oxide",
      "Aluminum oxide (alumina) is extremely hard, used as sandpaper and polishing compound. It's also the main component of rubies and sapphires (colored by trace elements). Alumina protects aluminum from rusting by forming a protective layer." },
    { buildTitaniumDioxide, "Titania/TiO2", CAT_METALS, "White pigment",
      "Titanium dioxide is bright white, used in paints, cosmetics, and sunscreen. It's one of the whitest substances known, scattering all visible light. Nanoparticle TiO2 in sunscreen is controversial for potential toxicity." },
    // === NEUROTRANSMITTERS (199-202) ===
    { buildAcetylcholine, "Acetylcholine", CAT_NEUROTRANS, "Muscle/memory",
      "Acetylcholine enables muscle contraction and memory formation at synapses. Your brain uses acetylcholine to consolidate memories; Alzheimer's involves acetylcholine deficiency. Some pesticides work by blocking acetylcholine breakdown." },
    { buildGABA, "GABA", CAT_NEUROTRANS, "Inhibitory NT",
      "GABA is your brain's main inhibitory neurotransmitter, calming excitement. Benzodiazepine tranquilizers work by enhancing GABA function. Alcohol also enhances GABA, explaining its sedative effects but also addiction potential." },
    { buildGlutamate, "Glutamate", CAT_NEUROTRANS, "Excitatory NT",
      "Glutamate is your brain's main excitatory neurotransmitter, essential for learning and memory. Excessive glutamate damages neurons (excitotoxicity), linked to Alzheimer's and stroke. MSG in food provides glutamate, though it doesn't cross the blood-brain barrier." },
    { buildNorepinephrine, "Norepinephrine", CAT_NEUROTRANS, "Fight-or-flight",
      "Norepinephrine is your fight-or-flight hormone released by the sympathetic nervous system. It increases heart rate, blood pressure, and alertness during stress. Stimulant drugs like cocaine and amphetamines increase norepinephrine levels." },
    // === ENERGY MOLECULES (203-205) ===
    { buildATP, "ATP", CAT_ENERGY, "Energy currency",
      "ATP (adenosine triphosphate) is your cells' energy currency, released when cells burn glucose. You produce about your body weight in ATP daily, recycling it constantly. Every muscle contraction, nerve impulse, and heartbeat requires ATP hydrolysis." },
    { buildADP, "ADP", CAT_ENERGY, "Adenosine diphosphate",
      "ADP (adenosine diphosphate) is ATP with one phosphate removed, storing energy release. It's quickly recharged back to ATP by mitochondria. The ATP-ADP cycle occurs millions of times per second in your body." },
    { buildNADH, "NADH", CAT_ENERGY, "Electron carrier",
      "NADH shuttles electrons through cellular respiration, carrying energy from glucose breakdown. When NADH donates electrons, it becomes NAD+, recycled endlessly. Performance athletes sometimes supplement with NAD+ precursors, though benefits are unproven." },
    // === ANTIBIOTICS & IMMUNOSUPPRESSANTS (206-210) ===
    { buildPenicillinG, "Penicillin G", CAT_PHARMA, "Beta-lactam",
      "Penicillin G was the first widely used antibiotic, discovered accidentally by Fleming in 1928. It revolutionized medicine, saving millions from bacterial infections. It works by breaking bacterial cell walls, which human cells lack." },
    { buildAmoxicillin, "Amoxicillin", CAT_PHARMA, "Broad-spectrum",
      "Amoxicillin is a penicillin derivative with broader bacterial spectrum and better absorption. It's one of the most prescribed antibiotics worldwide. Amoxicillin remains highly effective because most bacteria haven't developed resistance." },
    { buildVancomycin, "Vancomycin", CAT_PHARMA, "Glycopeptide",
      "Vancomycin is a powerful antibiotic used for serious infections and antibiotic-resistant bacteria. It works differently from penicillin, making it effective against resistant strains. IV vancomycin is the last resort for some life-threatening infections." },
    { buildMethotrexate, "Methotrexate", CAT_PHARMA, "Immunosuppressant",
      "Methotrexate blocks folic acid metabolism, slowing rapidly dividing cells. Doctors use it for cancer, autoimmune diseases, and severe rheumatoid arthritis. It's toxic at high doses but effective at lower doses with careful monitoring." },
    { buildHydroxychloroquine, "Hydroxychloroquine", CAT_PHARMA, "Antimalarial/Immunosuppressant",
      "Hydroxychloroquine reduces inflammation by interfering with immune cell function. It's used for lupus, rheumatoid arthritis, and historically for malaria prevention. It gained attention during COVID-19 despite unproven benefits in treating the virus." },
    { buildFluconazole, "Fluconazole", CAT_PHARMA, "Antifungal triazole",
      "Fluconazole is a triazole antifungal that inhibits fungal cell membrane synthesis by blocking lanosterol 14α-demethylase. It effectively treats candida infections (thrush), cryptococcal meningitis, and other systemic fungal infections. Fluconazole is well-absorbed orally and penetrates body tissues well, making it a first-line treatment for many fungal infections. Its broad spectrum and favorable pharmacokinetics have made it one of the most prescribed antifungals worldwide. It works by disrupting ergosterol synthesis in fungal cell membranes, which is absent in human cells, providing good selectivity. Fluconazole has been used successfully for over 30 years since its introduction in 1990, and it remains effective despite some emerging resistant strains." },
    // === EXPLOSIVES (209-211) ===
    { buildTNT, "TNT", CAT_OTHER, "Trinitrotoluene",
      "TNT (trinitrotoluene) is a stable explosive used in military applications and mining. It requires a detonator to explode, making it safer to handle than other explosives. The yellow color of TNT can stain skin and clothing permanently." },
    { buildNitroglycerin, "Nitroglycerin", CAT_PHARMA, "Explosive/angina",
      "Nitroglycerin relieves angina (heart pain) by dilating blood vessels, improving blood flow. It was first synthesized in 1847 but took until the 1870s to recognize medical use. Ironically, Nobel Prize creator Alfred Nobel discovered this use treating his chest pain." },
    { buildRDX, "RDX", CAT_OTHER, "Military explosive",
      "RDX is a military explosive more powerful than TNT, used in detonators and shaped charges. It's less stable than TNT, requiring careful handling. Modern ammunition often uses RDX in primer formulations." },
    // === SWEETENERS (212-218) ===
    { buildSucrose, "Sucrose/Table Sugar", CAT_SUGARS, "Disaccharide",
      "Sucrose is table sugar made by plants combining glucose and fructose. Plants make it for energy storage; animals can't synthesize it. Sucrose was once rare and precious; now it's leading contributor to obesity and diabetes." },
    { buildLactose, "Lactose", CAT_SUGARS, "Milk sugar disaccharide",
      "Lactose is the natural sugar in milk, found in mammals from birth. It's made of glucose and galactose joined together. Many humans lose lactase enzyme after childhood, causing lactose intolerance. Lactose-free milk is made by breaking lactose into simple sugars." },
    { buildMaltose, "Maltose", CAT_SUGARS, "Malt sugar disaccharide",
      "Maltose is created when enzymes break down starch in grains and potatoes. Brewers use it to ferment beer; your saliva starts breaking starch into maltose when chewing. Malt sugar is slightly less sweet than table sugar but has a pleasant, subtle flavor." },
    { buildGalactose, "Galactose", CAT_SUGARS, "C6H12O6 - Milk monosaccharide",
      "Galactose is a monosaccharide derived from lactose in milk. It differs from glucose by the position of one hydroxyl group, yet this tiny difference makes it crucial for brain development. Your body uses galactose to make galactose-cerebroside, a critical component of myelin insulation." },
    { buildAspartame, "Aspartame", CAT_OTHER, "NutraSweet",
      "Aspartame is a sweetener 200 times sweeter than sugar with no calories. It's made from aspartic acid and phenylalanine, amino acids. People with phenylketonuria (PKU) cannot safely consume aspartame due to phenylalanine." },
    { buildSaccharin, "Saccharin", CAT_OTHER, "Sweet'N Low",
      "Saccharin is the oldest artificial sweetener, 300 times sweeter than sugar. Once suspected of causing cancer, decades of research proved it safe. It has a bitter aftertaste that makes it less popular than newer sweeteners." },
    { buildSucralose, "Sucralose/Splenda", CAT_OTHER, "Chlorinated sugar",
      "Sucralose is made by replacing hydroxyl groups in sucrose with chlorine. It's 600 times sweeter than sugar and completely calorie-free. Splenda is the most common brand, heat-stable enough for baking." },
    // === PESTICIDES (216-218) ===
    { buildDDT, "DDT", CAT_OTHER, "Pesticide",
      "DDT is an insecticide that revolutionized malaria prevention but later banned in most countries. It persists in the environment for decades, accumulating in animals. It causes eggshell thinning in birds, nearly driving eagles extinct." },
    { buildGlyphosate, "Glyphosate/Roundup", CAT_OTHER, "Herbicide",
      "Glyphosate (Roundup) kills weeds by inhibiting plant protein synthesis. It's the world's most used herbicide, used with genetically modified crops. The WHO classified it as 'probably carcinogenic,' sparking controversy and lawsuits." },
    { buildMalathion, "Malathion", CAT_OTHER, "Insecticide",
      "Malathion is an insecticide used against mosquitoes and agricultural pests. It's less toxic than older organophosphates but still dangerous to humans. California uses aerial malathion spraying against fruit flies, causing public concern." },
    // === FRAGRANCES (219-222) ===
    { buildVanillin, "Vanillin/Vanilla", CAT_FLAVORS, "Vanilla flavor",
      "Vanillin is the compound giving vanilla its distinctive flavor and aroma. It's also found in turmeric, clove oil, and wood smoke. Most 'vanilla' flavoring is synthetic vanillin, cheaper and faster than extracting from vanilla beans." },
    { buildLimonene, "Limonene/Citrus", CAT_FLAVORS, "Citrus scent",
      "Limonene is the terpene giving citrus fruits their fresh smell. It's used in cleaning products and has potential cancer-fighting properties. Your body absorbs limonene through skin, which is why lemon scent feels refreshing." },
    { buildMenthol, "Menthol/Mint", CAT_FLAVORS, "Cooling sensation",
      "Menthol gives peppermint its cooling sensation by activating temperature-sensing nerve receptors. It tricks your nerves into thinking your mouth is cold. Menthol appears in cough drops, toothpaste, and sports creams." },
    { buildCinnamaldehyde, "Cinnamaldehyde/Cinnamon", CAT_FLAVORS, "Cinnamon",
      "Cinnamaldehyde is the compound giving cinnamon its spicy warmth. It activates pain receptors, which is why cinnamon burns slightly. Pure cinnamaldehyde can cause allergic reactions; true cinnamon is safer than cassia." },
    // === LACHRYMATORY (223) ===
    { buildPropanethialSoxide, "Onion Factor/Tears", CAT_FLAVORS, "Makes you cry",
      "Propanethial S-oxide is released when cutting onions, irritating your eyes. Your body interprets this as pain, triggering tears. Chilling onions before cutting reduces this compound's release." },
    // === ENVIRONMENTAL (224-225) ===
    { buildCFC12, "CFC-12/Freon", CAT_OTHER, "Ozone depleter",
      "CFC-12 (Freon) was used in refrigerators and air conditioners for decades. It's extremely stable in the atmosphere, persisting for 100+ years. CFCs destroy ozone high in the atmosphere, creating the ozone hole." },
    { buildSF6, "SF6/Greenhouse", CAT_OTHER, "Potent GHG",
      "SF6 (sulfur hexafluoride) is a potent greenhouse gas used in electrical equipment. It traps 23,500 times more heat than CO2 over 100 years. SF6 lasts 3,200 years in the atmosphere, making it extremely concerning." },
    // === ANESTHETICS (226-228) ===
    { buildLidocaine, "Lidocaine/Xylocaine", CAT_PHARMA, "Local anesthetic",
      "Lidocaine is a local anesthetic used in dentistry, minor surgery, and topical creams. It blocks nerve signals near the injection site without affecting your consciousness. Lidocaine patches provide pain relief for arthritis and shingles." },
    { buildKetamine, "Ketamine", CAT_CONTROLLED, "Dissociative",
      "Ketamine is a dissociative anesthetic used in surgery and emergency rooms. It produces dreamlike detachment where patients feel pain-free and unaware. Recently approved for depression treatment despite abuse potential." },
    { buildSevoflurane, "Sevoflurane", CAT_PHARMA, "Inhalation anesthetic",
      "Sevoflurane is a modern inhalation anesthetic with rapid onset and recovery. It's gentler than older anesthetics and rarely causes liver problems. Children wake up faster from sevoflurane, reducing post-operative confusion." },
    // === PSYCHEDELICS & HALLUCINOGENS (229-231) ===
    { buildLSD, "LSD", CAT_CONTROLLED, "Psychedelic",
      "LSD (lysergic acid diethylamide) is a powerful psychedelic altering perception and consciousness. It's not addictive but causes psychological dependence and flashbacks. Some researchers study LSD for treating depression and PTSD." },
    { buildPsilocybin, "Psilocybin/Shrooms", CAT_CONTROLLED, "Magic mushroom",
      "Psilocybin is the active compound in magic mushrooms, producing hallucinations and altered thinking. Research suggests psilocybin-assisted therapy helps depression and cluster headaches. It's slowly becoming legal for therapeutic research." },
    { buildMescaline, "Mescaline/Peyote", CAT_CONTROLLED, "Peyote cactus",
      "Mescaline is a psychedelic in peyote cactus and San Pedro cactus, used in indigenous ceremonies. It causes profound hallucinations lasting 12 hours. Indigenous peoples use mescaline for spiritual experiences and healing." },
    // === SYNTHETIC DRUGS & NOVEL PSYCHOACTIVE SUBSTANCES (232-234) ===
    { buildPCP, "PCP/Phencyclidine", CAT_CONTROLLED, "Dissociative hallucinogen",
      "PCP (phencyclidine) is a powerful dissociative drug causing detachment from reality and pain. It's particularly dangerous because it triggers bizarre, aggressive behavior. PCP users sometimes harm themselves without realizing it." },
    { buildMDMA, "MDMA/Ecstasy", CAT_CONTROLLED, "Party drug",
      "MDMA (Ecstasy) increases serotonin and dopamine, creating euphoria and empathy. Users experience increased heart rate and body temperature. Repeated use damages serotonin neurons, causing lasting mood problems." },
    { buildMethylone, "Methylone/Bath Salts", CAT_CONTROLLED, "Synthetic cathinone",
      "Methylone (Bath Salts) is a synthetic cathinone similar to MDMA but more dangerous. It causes extreme paranoia, violent behavior, and hallucinations. Emergency rooms report severe complications including rapid heart rate and seizures." },
    // === PRESCRIPTION DRUGS (235-238) ===
    { buildValium, "Valium/Diazepam", CAT_CONTROLLED, "Benzodiazepine tranquilizer",
      "Valium (diazepam) is a benzodiazepine tranquilizer used for anxiety and muscle spasms. It was wildly overprescribed in the 1960s-70s, creating addiction epidemics. Long-term use causes tolerance, making higher doses necessary." },
    { buildMethamphetamine, "Methamphetamine/Crystal Meth", CAT_CONTROLLED, "Powerful stimulant",
      "Methamphetamine (crystal meth) is a powerful stimulant causing euphoria and extreme wakefulness. It's highly addictive and devastates teeth, skin, and brain dopamine neurons. Chronic use causes paranoia, violence, and psychosis." },
    { buildCodeine, "Codeine", CAT_PHARMA, "Opioid pain reliever",
      "Codeine is a mild opioid pain reliever and cough suppressant in many medicines. Your body converts codeine to morphine, the active form. Some people lack the enzyme for conversion, making it ineffective." },
    { buildMethadone, "Methadone", CAT_PHARMA, "Synthetic opioid agonist",
      "Methadone is a synthetic opioid used to treat opioid addiction. It prevents withdrawal symptoms and blocks heroin's euphoria. Methadone maintenance allows addicts to function and rebuild lives." },
    // === ANTIHISTAMINES & DECONGESTANTS (239-244) ===
    { buildDiphenhydramine, "Diphenhydramine/Benadryl", CAT_PHARMA, "1st gen antihistamine",
      "Diphenhydramine (Benadryl) is a first-generation antihistamine causing drowsiness. It's used for allergies and as a sleep aid, but newer antihistamines are safer. Anticholinergic effects increase dementia risk in elderly people." },
    { buildChlorpheniramine, "Chlorpheniramine/Chlor-Trimeton", CAT_PHARMA, "1st gen antihistamine",
      "Chlorpheniramine is an older antihistamine causing significant drowsiness. It's used in allergy medications and cough syrups combined with other drugs. Newer antihistamines have fewer side effects." },
    { buildLoratadine, "Loratadine/Claritin", CAT_PHARMA, "2nd gen antihistamine",
      "Loratadine (Claritin) is a non-drowsy second-generation antihistamine. It targets histamine receptors without crossing the blood-brain barrier. Most effective when taken regularly before allergy season starts." },
    { buildCetirizine, "Cetirizine/Zyrtec", CAT_PHARMA, "2nd gen antihistamine",
      "Cetirizine (Zyrtec) is another non-drowsy antihistamine for allergies. It's long-acting, effective for 24 hours from one dose. Some people experience mild drowsiness despite 'non-drowsy' labeling." },
    { buildFexofenadine, "Fexofenadine/Allegra", CAT_PHARMA, "2nd gen antihistamine",
      "Fexofenadine (Allegra) is a non-drowsy antihistamine that doesn't cross the blood-brain barrier. Food reduces absorption significantly, so it's best taken on empty stomach. It's effective within 30 minutes." },
    { buildPseudoephedrine, "Pseudoephedrine", CAT_PHARMA, "Decongestant",
      "Pseudoephedrine is a decongestant that shrinks nasal blood vessels, reducing congestion. It's restricted in many stores because it's used to make methamphetamine. Phenylephrine is a safer alternative but less effective." },
    // === BLOOD PRESSURE MEDICATIONS (245-263) ===
    // ACE Inhibitors
    { buildLisinopril, "Lisinopril", CAT_PHARMA, "ACE inhibitor",
      "Lisinopril (ACE inhibitor) lowers blood pressure by relaxing blood vessels. It prevents angiotensin II from constricting vessels. ACE inhibitors are first-line treatment for hypertension and heart failure." },
    { buildEnalapril, "Enalapril", CAT_PHARMA, "ACE inhibitor",
      "Enalapril is another ACE inhibitor that must be converted to its active form by the liver. It's longer-acting than some ACE inhibitors, effective for 24 hours. It's used after heart attacks to prevent further damage." },
    { buildRamipril, "Ramipril", CAT_PHARMA, "ACE inhibitor",
      "Ramipril (Altace) is an ACE inhibitor used for hypertension and heart disease prevention. Studies show it reduces stroke and heart attack risk in high-risk patients. It's taken once or twice daily." },
    // ARBs
    { buildLosartan, "Losartan", CAT_PHARMA, "ARB",
      "Losartan (ARB) blocks angiotensin receptors instead of reducing angiotensin production. ARBs work differently than ACE inhibitors and have fewer side effects. They're alternative first-line treatments for hypertension." },
    { buildValsartan, "Valsartan", CAT_PHARMA, "ARB",
      "Valsartan is another ARB used for hypertension and heart failure. It's effective when combined with other medications. Some batches were recalled for contamination concerns in recent years." },
    { buildTelmisartan, "Telmisartan", CAT_PHARMA, "ARB",
      "Telmisartan is a long-acting ARB requiring once-daily dosing. It's effective for 24 hours after a single dose. Patients appreciate the convenience of once-daily treatment." },
    // Beta Blockers
    { buildMetoprolol, "Metoprolol", CAT_PHARMA, "Beta blocker",
      "Metoprolol is a beta blocker slowing heart rate and reducing blood pressure. It's used after heart attacks to protect the heart. Metoprolol reduces anxiety by calming the sympathetic nervous system." },
    { buildAtenolol, "Atenolol", CAT_PHARMA, "Beta blocker",
      "Atenolol is a cardioselective beta blocker primarily affecting the heart. It's popular for hypertension and migraine prevention. Atenolol can cause fatigue and reduced sexual function." },
    { buildCarvedilol, "Carvedilol", CAT_PHARMA, "Beta blocker",
      "Carvedilol is a combined alpha/beta blocker, more effective than pure beta blockers. It also provides additional protection for the heart. Carvedilol improves survival after heart attacks." },
    // Calcium Channel Blockers
    { buildAmlodipine, "Amlodipine", CAT_PHARMA, "Calcium channel blocker",
      "Amlodipine (Norvasc) is a calcium channel blocker that relaxes blood vessels. It's effective with minimal side effects and once-daily dosing. Ankle swelling is a common side effect affecting about 10% of users." },
    { buildDiltiazem, "Diltiazem", CAT_PHARMA, "Calcium channel blocker",
      "Diltiazem is a calcium channel blocker and heart rate slower used for hypertension and angina. It's particularly useful in patients who can't tolerate beta blockers. Diltiazem also has antiarrhythmic effects." },
    { buildNifedipine, "Nifedipine", CAT_PHARMA, "Calcium channel blocker",
      "Nifedipine is a potent calcium channel blocker used for severe hypertension and angina. Immediate-release forms work rapidly; extended-release forms provide stable control. Nifedipine causes reflex tachycardia if you drop blood pressure too fast." },
    // Diuretics
    { buildHydrochlorothiazide, "Hydrochlorothiazide/HCTZ", CAT_PHARMA, "Diuretic",
      "HCTZ is a thiazide diuretic, one of the oldest blood pressure drugs. It works by reducing body fluid volume. HCTZ can cause gout and worsen diabetes due to glucose effects." },
    { buildFurosemide, "Furosemide/Lasix", CAT_PHARMA, "Diuretic",
      "Furosemide (Lasix) is a loop diuretic stronger than HCTZ, used for fluid overload. It works on the kidney's loop of Henle. Furosemide can cause severe electrolyte imbalances if not monitored." },
    { buildSpironolactone, "Spironolactone", CAT_PHARMA, "Potassium-sparing diuretic",
      "Spironolactone is a potassium-sparing diuretic blocking aldosterone. Unlike other diuretics, it retains potassium, requiring careful monitoring. It's used in heart failure and resistant hypertension." },
    // === ANTIVIRALS (264-269) ===
    { buildOseltamivir, "Oseltamivir/Tamiflu", CAT_PHARMA, "Influenza antiviral",
      "Oseltamivir (Tamiflu) stops flu virus replication if taken early. It reduces symptoms by 1-2 days if started within 48 hours. Its effectiveness against new flu strains varies yearly." },
    { buildNirmatrelvir, "Nirmatrelvir (Paxlovid)", CAT_PHARMA, "COVID protease inhibitor",
      "Nirmatrelvir is Paxlovid's main component, a protease inhibitor blocking COVID-19 virus multiplication. Taken with ritonavir as a booster, it's highly effective for early COVID-19. It must be started within 5 days of symptom onset." },
    { buildRitonavir, "Ritonavir (Paxlovid component)", CAT_PHARMA, "Protease inhibitor",
      "Ritonavir is a protease inhibitor from the 1990s that boosts other drugs' levels. In Paxlovid, it increases nirmatrelvir concentration. Ritonavir alone was an HIV medication before new better drugs emerged." },
    { buildMolnupiravir, "Molnupiravir/Lagevrio", CAT_PHARMA, "COVID antiviral",
      "Molnupiravir (Lagevrio) is an antiviral that causes viral mutations making replication impossible. Taken early, it reduces hospitalization and death. It's less effective than Paxlovid but works on variants better." },
    { buildRemdesivir, "Remdesivir/Veklury", CAT_PHARMA, "COVID/Ebola antiviral",
      "Remdesivir (Veklury) stops viral replication by incorporating into viral RNA. Developed for Ebola, it became famous treating severe COVID-19. Early studies showed modest benefit, reducing hospitalization time." },
    { buildDexamethasone, "Dexamethasone", CAT_PHARMA, "Corticosteroid",
      "Dexamethasone is a corticosteroid that reduces severe COVID-19 inflammation. It's most effective when oxygen is needed but doesn't help mild cases. It's also used for allergies, autoimmune diseases, and cancer." },
    // === RANDOM (270) ===
    { buildRandomMolecule, "Random", CAT_OTHER, "Random structure",
      "Random structures are computer-generated molecules for educational visualization. They help students understand chemical bonding diversity. Random molecules demonstrate how atoms can arrange in countless configurations." },
    // === EXOTIC STRUCTURES (233-249) ===
    { buildCubane, "Cubane", CAT_ORGANIC, "C8H8 - Cube-shaped",
      "Cubane is a remarkable cage structure where eight carbons form a perfect cube. It's highly strained but surprisingly stable. Scientists study cubane derivatives for potential drug development." },
    { buildBuckminsterfullerene, "Buckyball/C60", CAT_ORGANIC, "Carbon soccer ball",
      "Buckyball (C60) is 60 carbons arranged like a soccer ball, discovered in 1985. It's a single sheet of graphite rolled into a sphere. Scientists have explored buckyballs for drug delivery and materials." },
    { buildAdamantane, "Adamantane", CAT_ORGANIC, "C10H16 - Diamond cage",
      "Adamantane is a cage structure mimicking diamond's carbon arrangement. It's used in antivirals and as a building block for complex molecules. Its symmetry makes it useful in medicinal chemistry." },
    { buildChlorophyll, "Chlorophyll a", CAT_ENERGY, "Plant pigment core",
      "Chlorophyll absorbs light for photosynthesis, giving plants their green color. Its central magnesium atom captures photons. Chlorophyll is nearly identical to hemoglobin except hemoglobin has iron." },
    { buildXenonDifluoride, "Xenon Difluoride (XeF2)", CAT_OTHER, "Linear noble gas",
      "Xenon difluoride (XeF2) has a linear geometry with xenon bonded to fluorines. Despite noble gas inertness, XeF2 forms via intense UV. It's used in nuclear fuel reprocessing." },
    { buildRheniumChlorideDimer, "[Re2Cl8]^2- Dimer", CAT_OTHER, "Quadruple metal bond",
      "Rhenium chloride dimers have quadruple bonds between rhenium atoms. This was shocking to chemists; metal-metal quadruple bonds were thought impossible. It's a landmark in inorganic chemistry." },
    { buildTungstenHexacarbonyl, "W(CO)6 Hexacarbonyl", CAT_OTHER, "Octahedral carbonyl",
      "Tungsten hexacarbonyl (W(CO)6) is octahedral with carbons surrounding tungsten. Carbonyls are fundamental in organometallic chemistry. It's used in catalysis and organic synthesis." },
    { buildOCarborane, "o-Carborane (C2B10H12)", CAT_OTHER, "Icosahedral cage",
      "o-Carborane is an icosahedral cage with boron and carbon atoms. It's incredibly stable and used in boron neutron capture therapy for cancer. Its geometry resembles an icosahedron." },
    { buildDodecaborate, "[B12H12]2- Dodecaborate", CAT_OTHER, "Perfect icosahedron",
      "Dodecaborate ([B12H12]2-) is a perfect icosahedron of 12 boron atoms. It's one of the most symmetrical molecules in chemistry. Dodecaborate salts are used in boron research." },
    { buildCatenane, "[2]Catenane", CAT_OTHER, "Interlocked rings",
      "Catenane molecules consist of interlocked rings that cannot pass through each other. Creating interlocked rings requires careful synthesis. Catenanes represent breakthrough in topological chemistry." },
    { buildRotaxane, "Rotaxane", CAT_OTHER, "Ring on axle",
      "Rotaxane has a ring threaded on an axle, prevented from sliding off by bulky stoppers. Like a bead on a rod, the ring can rotate freely. Rotaxanes demonstrate mechanical properties at molecular scale." },
    { buildTrefoilKnot, "Molecular Trefoil Knot", CAT_OTHER, "Knotted loop",
      "Molecular trefoil knots are knotted loops resembling a trinity knot symbol. Scientists synthesized the first molecular knot in 1989, stunning chemists. Knotted molecules explore topology in chemistry." },
    { buildPropellane, "[1.1.1]Propellane", CAT_ORGANIC, "Inverted C-C bond",
      "Propellane ([1.1.1]) has an unusual inverted C-C bond between bridgeheads. It's incredibly strained, showing chemistry's limits. Propellanes are fascinating for theoretical chemistry." },
    { buildPrismane, "Prismane (C6H6)", CAT_ORGANIC, "Triangular prism",
      "Prismane is a benzene isomer with triangular prism geometry. It converts to benzene over time at room temperature. Prismane demonstrates how aromatic systems can rearrange." },
    { buildTetrahedrane, "Tetrahedrane (C4H4)", CAT_ORGANIC, "Carbon tetrahedron",
      "Tetrahedrane (C4H4) is a tetrahedron of four carbons. It's one of the smallest possible ring hydrocarbon. Tetrahedrane is extremely strained but has been synthesized." },
    { buildDewarBenzene, "Dewar Benzene", CAT_ORGANIC, "Bent benzene isomer",
      "Dewar benzene is a bent isomer of benzene with unusual bonding. It converts to benzene when heated. Dewar benzene challenges our understanding of aromaticity and resonance." },
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

const char* molecule_get_long_description(int index) {
    if (index < 0 || index >= NUM_MOLECULES) return "";
    const char* desc = molecules[index].longDescription;
    return desc ? desc : "No detailed description available yet.";
}

void molecule_build(int index, Molecule* mol) {
    if (index < 0 || index >= NUM_MOLECULES) index = 0;
    molecules[index].builder(mol);
}

void molecule_build_random(Molecule* mol) {
    buildRandomMolecule(mol);
}
