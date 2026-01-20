/*
 * MolVis - Molecule Database Implementation
 *
 * Contains molecule presets including organic compounds, drugs, and more.
 * Molecules are organized into categories for easy browsing.
 */

#include "molecule_db.h"
#include <string.h>
#include <stdlib.h>
#include <math.h>

// ============== Category Names ==============

static const char* categoryNames[] = {
    "Simple Molecules",
    "Organic Compounds",
    "Biochemistry",
    "Pharmaceuticals",
    "Household Chemicals",
    "Other"
};

const char* molecule_get_category_name(int category) {
    if (category < 0 || category >= CAT_COUNT) return "Unknown";
    return categoryNames[category];
}

// ============== Helper Functions ==============

static void addAtom(Molecule* mol, float x, float y, float z, int type) {
    if (mol->numAtoms >= MAX_ATOMS) return;
    Atom* a = &mol->atoms[mol->numAtoms];
    a->x = x;
    a->y = y;
    a->z = z;
    a->type = type;

    // Set radius based on type
    float radii[] = {
        0.25f,  // H
        0.40f,  // C
        0.38f,  // N
        0.35f,  // O
        0.45f,  // P
        0.45f,  // S
        0.45f,  // Cl
        0.50f,  // Br
        0.35f,  // F
        0.55f,  // I
        0.55f,  // Na
        0.48f,  // Si
        0.42f,  // B
        0.55f,  // Fe
        0.50f,  // Cu
        0.50f,  // Al
        0.52f,  // Ti
        0.55f,  // Pt
    };
    a->radius = radii[type < 18 ? type : 0];
    mol->numAtoms++;
}

static void addBond(Molecule* mol, int a1, int a2, int order) {
    if (mol->numBonds >= MAX_BONDS) return;
    Bond* b = &mol->bonds[mol->numBonds];
    b->atom1 = a1;
    b->atom2 = a2;
    b->order = order;
    mol->numBonds++;
}

static void initMolecule(Molecule* mol, const char* name) {
    mol->numAtoms = 0;
    mol->numBonds = 0;
    strncpy(mol->name, name, 63);
    mol->name[63] = '\0';
}

// ============== Molecule Builders ==============

static void buildWater(Molecule* mol) {
    initMolecule(mol, "Water (H2O)");
    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_O);
    addAtom(mol, 0.76f, 0.59f, 0.0f, ATOM_H);
    addAtom(mol, -0.76f, 0.59f, 0.0f, ATOM_H);
    addBond(mol, 0, 1, 1);
    addBond(mol, 0, 2, 1);
}

static void buildMethane(Molecule* mol) {
    initMolecule(mol, "Methane (CH4)");
    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_C);
    addAtom(mol, 0.63f, 0.63f, 0.63f, ATOM_H);
    addAtom(mol, -0.63f, -0.63f, 0.63f, ATOM_H);
    addAtom(mol, -0.63f, 0.63f, -0.63f, ATOM_H);
    addAtom(mol, 0.63f, -0.63f, -0.63f, ATOM_H);
    addBond(mol, 0, 1, 1);
    addBond(mol, 0, 2, 1);
    addBond(mol, 0, 3, 1);
    addBond(mol, 0, 4, 1);
}

static void buildEthanol(Molecule* mol) {
    initMolecule(mol, "Ethanol (C2H5OH)");
    // Carbon chain
    addAtom(mol, -0.7f, 0.0f, 0.0f, ATOM_C);  // 0: CH3
    addAtom(mol, 0.7f, 0.0f, 0.0f, ATOM_C);   // 1: CH2
    addAtom(mol, 1.8f, 0.0f, 0.0f, ATOM_O);   // 2: OH
    // Hydrogens on C1
    addAtom(mol, -1.1f, 1.0f, 0.0f, ATOM_H);
    addAtom(mol, -1.1f, -0.5f, 0.87f, ATOM_H);
    addAtom(mol, -1.1f, -0.5f, -0.87f, ATOM_H);
    // Hydrogens on C2
    addAtom(mol, 0.7f, 0.5f, 0.87f, ATOM_H);
    addAtom(mol, 0.7f, 0.5f, -0.87f, ATOM_H);
    // Hydrogen on O
    addAtom(mol, 2.3f, 0.8f, 0.0f, ATOM_H);
    // Bonds
    addBond(mol, 0, 1, 1);
    addBond(mol, 1, 2, 1);
    addBond(mol, 0, 3, 1);
    addBond(mol, 0, 4, 1);
    addBond(mol, 0, 5, 1);
    addBond(mol, 1, 6, 1);
    addBond(mol, 1, 7, 1);
    addBond(mol, 2, 8, 1);
}

static void buildBenzene(Molecule* mol) {
    initMolecule(mol, "Benzene (C6H6)");
    float r = 1.4f;
    for (int i = 0; i < 6; i++) {
        float angle = i * 3.14159f / 3.0f;
        addAtom(mol, r * cosf(angle), r * sinf(angle), 0.0f, ATOM_C);
    }
    for (int i = 0; i < 6; i++) {
        float angle = i * 3.14159f / 3.0f;
        addAtom(mol, 2.4f * cosf(angle), 2.4f * sinf(angle), 0.0f, ATOM_H);
    }
    // Alternating single/double bonds
    for (int i = 0; i < 6; i++) {
        addBond(mol, i, (i + 1) % 6, (i % 2 == 0) ? 2 : 1);
        addBond(mol, i, i + 6, 1);
    }
}

static void buildCaffeine(Molecule* mol) {
    initMolecule(mol, "Caffeine (C8H10N4O2)");
    // Purine ring system (simplified coordinates)
    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_C);      // 0
    addAtom(mol, 1.2f, 0.7f, 0.0f, ATOM_N);      // 1
    addAtom(mol, 2.4f, 0.0f, 0.0f, ATOM_C);      // 2
    addAtom(mol, 2.4f, -1.4f, 0.0f, ATOM_N);     // 3
    addAtom(mol, 1.2f, -2.1f, 0.0f, ATOM_C);     // 4
    addAtom(mol, 0.0f, -1.4f, 0.0f, ATOM_N);     // 5
    addAtom(mol, 3.6f, 0.7f, 0.0f, ATOM_N);      // 6
    addAtom(mol, 3.6f, -2.1f, 0.0f, ATOM_C);     // 7
    // Oxygens (carbonyl)
    addAtom(mol, -1.2f, 0.7f, 0.0f, ATOM_O);     // 8
    addAtom(mol, 1.2f, -3.3f, 0.0f, ATOM_O);     // 9
    // Methyl groups (simplified)
    addAtom(mol, 1.2f, 2.0f, 0.0f, ATOM_C);      // 10: N-CH3
    addAtom(mol, -1.2f, -2.1f, 0.0f, ATOM_C);    // 11: N-CH3
    addAtom(mol, 3.6f, -3.3f, 0.0f, ATOM_C);     // 12: N-CH3
    // Bonds
    addBond(mol, 0, 1, 1);
    addBond(mol, 1, 2, 1);
    addBond(mol, 2, 3, 2);
    addBond(mol, 3, 4, 1);
    addBond(mol, 4, 5, 1);
    addBond(mol, 5, 0, 1);
    addBond(mol, 2, 6, 1);
    addBond(mol, 3, 7, 1);
    addBond(mol, 6, 7, 2);
    addBond(mol, 0, 8, 2);
    addBond(mol, 4, 9, 2);
    addBond(mol, 1, 10, 1);
    addBond(mol, 5, 11, 1);
    addBond(mol, 7, 12, 1);
}

static void buildAspirin(Molecule* mol) {
    initMolecule(mol, "Aspirin (C9H8O4)");
    // Benzene ring
    float r = 1.4f;
    for (int i = 0; i < 6; i++) {
        float angle = i * 3.14159f / 3.0f;
        addAtom(mol, r * cosf(angle), r * sinf(angle), 0.0f, ATOM_C);
    }
    // Carboxylic acid group
    addAtom(mol, 2.5f, 0.0f, 0.0f, ATOM_C);      // 6
    addAtom(mol, 3.2f, 1.0f, 0.0f, ATOM_O);      // 7 (=O)
    addAtom(mol, 3.0f, -1.2f, 0.0f, ATOM_O);     // 8 (-OH)
    // Acetyl group
    addAtom(mol, -0.7f, 2.4f, 0.0f, ATOM_O);     // 9 (ester O)
    addAtom(mol, -0.7f, 3.6f, 0.0f, ATOM_C);     // 10 (acetyl C)
    addAtom(mol, -1.9f, 4.3f, 0.0f, ATOM_C);     // 11 (CH3)
    addAtom(mol, 0.3f, 4.3f, 0.0f, ATOM_O);      // 12 (=O)
    // Bonds
    for (int i = 0; i < 6; i++) {
        addBond(mol, i, (i + 1) % 6, (i % 2 == 0) ? 2 : 1);
    }
    addBond(mol, 0, 6, 1);
    addBond(mol, 6, 7, 2);
    addBond(mol, 6, 8, 1);
    addBond(mol, 2, 9, 1);
    addBond(mol, 9, 10, 1);
    addBond(mol, 10, 11, 1);
    addBond(mol, 10, 12, 2);
}

static void buildDopamine(Molecule* mol) {
    initMolecule(mol, "Dopamine (C8H11NO2)");
    // Catechol ring
    float r = 1.4f;
    for (int i = 0; i < 6; i++) {
        float angle = i * 3.14159f / 3.0f;
        addAtom(mol, r * cosf(angle), r * sinf(angle), 0.0f, ATOM_C);
    }
    // Hydroxyl groups
    addAtom(mol, 2.0f, 1.4f, 0.0f, ATOM_O);      // 6
    addAtom(mol, 2.0f, -1.4f, 0.0f, ATOM_O);     // 7
    // Ethylamine chain
    addAtom(mol, -2.5f, 0.0f, 0.0f, ATOM_C);     // 8
    addAtom(mol, -3.5f, 1.0f, 0.0f, ATOM_C);     // 9
    addAtom(mol, -4.5f, 1.0f, 0.0f, ATOM_N);     // 10
    // Bonds
    for (int i = 0; i < 6; i++) {
        addBond(mol, i, (i + 1) % 6, (i % 2 == 0) ? 2 : 1);
    }
    addBond(mol, 1, 6, 1);
    addBond(mol, 5, 7, 1);
    addBond(mol, 3, 8, 1);
    addBond(mol, 8, 9, 1);
    addBond(mol, 9, 10, 1);
}

static void buildGlucose(Molecule* mol) {
    initMolecule(mol, "Glucose (C6H12O6)");
    // Pyranose ring (chair conformation simplified)
    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_C);      // 0 C1
    addAtom(mol, 1.5f, 0.3f, 0.0f, ATOM_C);      // 1 C2
    addAtom(mol, 2.3f, -0.9f, 0.0f, ATOM_C);     // 2 C3
    addAtom(mol, 1.5f, -2.1f, 0.0f, ATOM_C);     // 3 C4
    addAtom(mol, 0.0f, -1.8f, 0.0f, ATOM_C);     // 4 C5
    addAtom(mol, -0.8f, -0.9f, 0.0f, ATOM_O);    // 5 ring O
    addAtom(mol, -0.8f, -3.0f, 0.0f, ATOM_C);    // 6 C6 (CH2OH)
    // Hydroxyl oxygens
    addAtom(mol, -0.5f, 1.2f, 0.0f, ATOM_O);     // 7 OH on C1
    addAtom(mol, 2.0f, 1.5f, 0.0f, ATOM_O);      // 8 OH on C2
    addAtom(mol, 3.5f, -0.9f, 0.0f, ATOM_O);     // 9 OH on C3
    addAtom(mol, 2.0f, -3.3f, 0.0f, ATOM_O);     // 10 OH on C4
    addAtom(mol, -2.0f, -3.0f, 0.0f, ATOM_O);    // 11 OH on C6
    // Bonds
    addBond(mol, 0, 1, 1);
    addBond(mol, 1, 2, 1);
    addBond(mol, 2, 3, 1);
    addBond(mol, 3, 4, 1);
    addBond(mol, 4, 5, 1);
    addBond(mol, 5, 0, 1);
    addBond(mol, 4, 6, 1);
    addBond(mol, 0, 7, 1);
    addBond(mol, 1, 8, 1);
    addBond(mol, 2, 9, 1);
    addBond(mol, 3, 10, 1);
    addBond(mol, 6, 11, 1);
}

static void buildATP(Molecule* mol) {
    initMolecule(mol, "ATP (Adenosine Triphosphate)");
    // Adenine base (simplified)
    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_N);      // 0
    addAtom(mol, 1.2f, 0.7f, 0.0f, ATOM_C);      // 1
    addAtom(mol, 2.4f, 0.0f, 0.0f, ATOM_N);      // 2
    addAtom(mol, 2.4f, -1.4f, 0.0f, ATOM_C);     // 3
    addAtom(mol, 1.2f, -2.1f, 0.0f, ATOM_C);     // 4
    addAtom(mol, 0.0f, -1.4f, 0.0f, ATOM_C);     // 5
    addAtom(mol, -1.2f, -2.1f, 0.0f, ATOM_N);    // 6 (NH2)
    // Ribose sugar
    addAtom(mol, 3.6f, -2.1f, 0.0f, ATOM_C);     // 7 C1'
    addAtom(mol, 4.3f, -3.3f, 0.0f, ATOM_C);     // 8 C2'
    addAtom(mol, 5.5f, -2.6f, 0.0f, ATOM_C);     // 9 C3'
    addAtom(mol, 5.5f, -1.2f, 0.0f, ATOM_C);     // 10 C4'
    addAtom(mol, 4.3f, -0.5f, 0.0f, ATOM_O);     // 11 O (ring)
    addAtom(mol, 6.7f, -0.5f, 0.0f, ATOM_C);     // 12 C5'
    // Triphosphate
    addAtom(mol, 7.9f, -1.2f, 0.0f, ATOM_O);     // 13
    addAtom(mol, 9.1f, -0.5f, 0.0f, ATOM_P);     // 14 P-alpha
    addAtom(mol, 10.3f, -1.2f, 0.0f, ATOM_O);    // 15
    addAtom(mol, 11.5f, -0.5f, 0.0f, ATOM_P);    // 16 P-beta
    addAtom(mol, 12.7f, -1.2f, 0.0f, ATOM_O);    // 17
    addAtom(mol, 13.9f, -0.5f, 0.0f, ATOM_P);    // 18 P-gamma
    // Phosphate oxygens (simplified)
    addAtom(mol, 9.1f, 0.9f, 0.0f, ATOM_O);      // 19
    addAtom(mol, 11.5f, 0.9f, 0.0f, ATOM_O);     // 20
    addAtom(mol, 13.9f, 0.9f, 0.0f, ATOM_O);     // 21
    addAtom(mol, 15.1f, -1.2f, 0.0f, ATOM_O);    // 22
    // Bonds
    addBond(mol, 0, 1, 2);
    addBond(mol, 1, 2, 1);
    addBond(mol, 2, 3, 2);
    addBond(mol, 3, 4, 1);
    addBond(mol, 4, 5, 2);
    addBond(mol, 5, 0, 1);
    addBond(mol, 5, 6, 1);
    addBond(mol, 3, 7, 1);
    addBond(mol, 7, 8, 1);
    addBond(mol, 8, 9, 1);
    addBond(mol, 9, 10, 1);
    addBond(mol, 10, 11, 1);
    addBond(mol, 11, 7, 1);
    addBond(mol, 10, 12, 1);
    addBond(mol, 12, 13, 1);
    addBond(mol, 13, 14, 1);
    addBond(mol, 14, 15, 1);
    addBond(mol, 15, 16, 1);
    addBond(mol, 16, 17, 1);
    addBond(mol, 17, 18, 1);
    addBond(mol, 14, 19, 2);
    addBond(mol, 16, 20, 2);
    addBond(mol, 18, 21, 2);
    addBond(mol, 18, 22, 1);
}

static void buildRandom(Molecule* mol) {
    initMolecule(mol, "Random Molecule");
    int numAtoms = 8 + rand() % 20;
    int types[] = {ATOM_C, ATOM_C, ATOM_C, ATOM_N, ATOM_O, ATOM_H, ATOM_H, ATOM_S};

    for (int i = 0; i < numAtoms; i++) {
        float x = (rand() % 100 - 50) / 20.0f;
        float y = (rand() % 100 - 50) / 20.0f;
        float z = (rand() % 100 - 50) / 20.0f;
        int type = types[rand() % 8];
        addAtom(mol, x, y, z, type);
    }

    // Add random bonds
    for (int i = 0; i < mol->numAtoms - 1; i++) {
        if (rand() % 3 != 0) {
            int j = i + 1 + rand() % (mol->numAtoms - i - 1);
            if (j < mol->numAtoms) {
                addBond(mol, i, j, 1 + rand() % 2);
            }
        }
    }
}

// ============== Additional Simple Molecules ==============

static void buildAmmonia(Molecule* mol) {
    initMolecule(mol, "Ammonia (NH3)");
    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_N);
    addAtom(mol, 0.94f, 0.0f, 0.34f, ATOM_H);
    addAtom(mol, -0.47f, 0.81f, 0.34f, ATOM_H);
    addAtom(mol, -0.47f, -0.81f, 0.34f, ATOM_H);
    addBond(mol, 0, 1, 1);
    addBond(mol, 0, 2, 1);
    addBond(mol, 0, 3, 1);
}

static void buildCarbonDioxide(Molecule* mol) {
    initMolecule(mol, "Carbon Dioxide (CO2)");
    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_C);
    addAtom(mol, -1.16f, 0.0f, 0.0f, ATOM_O);
    addAtom(mol, 1.16f, 0.0f, 0.0f, ATOM_O);
    addBond(mol, 0, 1, 2);
    addBond(mol, 0, 2, 2);
}

static void buildHydrogenPeroxide(Molecule* mol) {
    initMolecule(mol, "Hydrogen Peroxide (H2O2)");
    addAtom(mol, -0.7f, 0.0f, 0.0f, ATOM_O);
    addAtom(mol, 0.7f, 0.0f, 0.0f, ATOM_O);
    addAtom(mol, -1.2f, 0.8f, 0.0f, ATOM_H);
    addAtom(mol, 1.2f, -0.8f, 0.0f, ATOM_H);
    addBond(mol, 0, 1, 1);
    addBond(mol, 0, 2, 1);
    addBond(mol, 1, 3, 1);
}

// ============== Additional Organic Molecules ==============

static void buildAcetone(Molecule* mol) {
    initMolecule(mol, "Acetone (C3H6O)");
    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_C);      // Central C
    addAtom(mol, -1.5f, 0.0f, 0.0f, ATOM_C);     // CH3
    addAtom(mol, 1.5f, 0.0f, 0.0f, ATOM_C);      // CH3
    addAtom(mol, 0.0f, 1.2f, 0.0f, ATOM_O);      // =O
    // Hydrogens
    addAtom(mol, -1.9f, 1.0f, 0.0f, ATOM_H);
    addAtom(mol, -1.9f, -0.5f, 0.87f, ATOM_H);
    addAtom(mol, -1.9f, -0.5f, -0.87f, ATOM_H);
    addAtom(mol, 1.9f, 1.0f, 0.0f, ATOM_H);
    addAtom(mol, 1.9f, -0.5f, 0.87f, ATOM_H);
    addAtom(mol, 1.9f, -0.5f, -0.87f, ATOM_H);
    addBond(mol, 0, 1, 1);
    addBond(mol, 0, 2, 1);
    addBond(mol, 0, 3, 2);
    addBond(mol, 1, 4, 1);
    addBond(mol, 1, 5, 1);
    addBond(mol, 1, 6, 1);
    addBond(mol, 2, 7, 1);
    addBond(mol, 2, 8, 1);
    addBond(mol, 2, 9, 1);
}

static void buildFormaldehyde(Molecule* mol) {
    initMolecule(mol, "Formaldehyde (CH2O)");
    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_C);
    addAtom(mol, 0.0f, 1.2f, 0.0f, ATOM_O);
    addAtom(mol, -0.9f, -0.5f, 0.0f, ATOM_H);
    addAtom(mol, 0.9f, -0.5f, 0.0f, ATOM_H);
    addBond(mol, 0, 1, 2);
    addBond(mol, 0, 2, 1);
    addBond(mol, 0, 3, 1);
}

static void buildAceticAcid(Molecule* mol) {
    initMolecule(mol, "Acetic Acid (CH3COOH)");
    addAtom(mol, -0.7f, 0.0f, 0.0f, ATOM_C);     // CH3
    addAtom(mol, 0.7f, 0.0f, 0.0f, ATOM_C);      // COOH
    addAtom(mol, 1.3f, 1.0f, 0.0f, ATOM_O);      // =O
    addAtom(mol, 1.3f, -1.0f, 0.0f, ATOM_O);     // -OH
    addAtom(mol, 2.0f, -1.5f, 0.0f, ATOM_H);     // H
    addAtom(mol, -1.1f, 1.0f, 0.0f, ATOM_H);
    addAtom(mol, -1.1f, -0.5f, 0.87f, ATOM_H);
    addAtom(mol, -1.1f, -0.5f, -0.87f, ATOM_H);
    addBond(mol, 0, 1, 1);
    addBond(mol, 1, 2, 2);
    addBond(mol, 1, 3, 1);
    addBond(mol, 3, 4, 1);
    addBond(mol, 0, 5, 1);
    addBond(mol, 0, 6, 1);
    addBond(mol, 0, 7, 1);
}

// ============== Household Chemicals ==============

static void buildTableSalt(Molecule* mol) {
    initMolecule(mol, "Table Salt (NaCl)");
    addAtom(mol, -1.2f, 0.0f, 0.0f, ATOM_NA);
    addAtom(mol, 1.2f, 0.0f, 0.0f, ATOM_CL);
    addBond(mol, 0, 1, 1);
}

static void buildBakingSoda(Molecule* mol) {
    initMolecule(mol, "Baking Soda (NaHCO3)");
    addAtom(mol, -2.0f, 0.0f, 0.0f, ATOM_NA);
    addAtom(mol, 0.0f, 0.0f, 0.0f, ATOM_C);
    addAtom(mol, 1.2f, 0.0f, 0.0f, ATOM_O);
    addAtom(mol, -0.6f, 1.0f, 0.0f, ATOM_O);
    addAtom(mol, -0.6f, -1.0f, 0.0f, ATOM_O);
    addAtom(mol, -0.6f, -2.0f, 0.0f, ATOM_H);
    addBond(mol, 1, 2, 2);
    addBond(mol, 1, 3, 1);
    addBond(mol, 1, 4, 1);
    addBond(mol, 4, 5, 1);
    addBond(mol, 0, 3, 1);
}

static void buildBleach(Molecule* mol) {
    initMolecule(mol, "Bleach (NaOCl)");
    addAtom(mol, -1.5f, 0.0f, 0.0f, ATOM_NA);
    addAtom(mol, 0.5f, 0.0f, 0.0f, ATOM_O);
    addAtom(mol, 2.0f, 0.0f, 0.0f, ATOM_CL);
    addBond(mol, 1, 2, 1);
    addBond(mol, 0, 1, 1);
}

// ============== Molecule Registry ==============

// Structure to hold molecule metadata
struct MoleculeInfo {
    MoleculeBuilder builder;
    const char* name;
    int category;
    const char* description;
};

static MoleculeInfo molecules[] = {
    // Simple Molecules
    { buildWater,           "Water (H2O)",              CAT_SIMPLE,    "Essential for life, universal solvent" },
    { buildAmmonia,         "Ammonia (NH3)",            CAT_SIMPLE,    "Pungent gas, used in fertilizers" },
    { buildCarbonDioxide,   "Carbon Dioxide (CO2)",     CAT_SIMPLE,    "Greenhouse gas, product of respiration" },
    { buildHydrogenPeroxide,"Hydrogen Peroxide (H2O2)", CAT_SIMPLE,    "Antiseptic, bleaching agent" },
    { buildMethane,         "Methane (CH4)",            CAT_SIMPLE,    "Natural gas, simplest hydrocarbon" },

    // Organic Compounds
    { buildEthanol,         "Ethanol (C2H5OH)",         CAT_ORGANIC,   "Alcohol in beverages, fuel additive" },
    { buildBenzene,         "Benzene (C6H6)",           CAT_ORGANIC,   "Aromatic ring, industrial solvent" },
    { buildAcetone,         "Acetone (C3H6O)",          CAT_ORGANIC,   "Nail polish remover, industrial solvent" },
    { buildFormaldehyde,    "Formaldehyde (CH2O)",      CAT_ORGANIC,   "Preservative, building block chemical" },
    { buildAceticAcid,      "Acetic Acid (CH3COOH)",    CAT_ORGANIC,   "Vinegar, food preservative" },

    // Biochemistry
    { buildGlucose,         "Glucose (C6H12O6)",        CAT_BIOCHEM,   "Blood sugar, primary energy source" },
    { buildATP,             "ATP",                      CAT_BIOCHEM,   "Energy currency of cells" },

    // Pharmaceuticals
    { buildCaffeine,        "Caffeine (C8H10N4O2)",     CAT_DRUGS,     "Stimulant in coffee and tea" },
    { buildAspirin,         "Aspirin (C9H8O4)",         CAT_DRUGS,     "Pain reliever, anti-inflammatory" },
    { buildDopamine,        "Dopamine (C8H11NO2)",      CAT_DRUGS,     "Neurotransmitter, pleasure chemical" },

    // Household Chemicals
    { buildTableSalt,       "Table Salt (NaCl)",        CAT_HOUSEHOLD, "Seasoning, preservative" },
    { buildBakingSoda,      "Baking Soda (NaHCO3)",     CAT_HOUSEHOLD, "Leavening agent, cleaning" },
    { buildBleach,          "Bleach (NaOCl)",           CAT_HOUSEHOLD, "Disinfectant, whitening agent" },

    // Other
    { buildRandom,          "Random Molecule",          CAT_OTHER,     "Procedurally generated structure" },
};

static const int NUM_MOLECULES = sizeof(molecules) / sizeof(molecules[0]);

// ============== Public API ==============

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
    buildRandom(mol);
}
