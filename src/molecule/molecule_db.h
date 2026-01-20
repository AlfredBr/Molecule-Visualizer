/*
 * MolVis - Molecule Database
 *
 * Contains molecule presets and generation functions.
 */

#ifndef MOLECULE_DB_H
#define MOLECULE_DB_H

// Use platform-agnostic molecule types
#include "molecule/molecule_types.h"

// Molecule categories for UI organization
enum MoleculeCategory {
    CAT_SIMPLE = 0,     // Simple molecules (Water, Ammonia, gases)
    CAT_ORGANIC,        // Organic compounds (Alcohols, Ethers, etc.)
    CAT_AMINO_ACIDS,    // Amino acids
    CAT_SUGARS,         // Sugars and carbohydrates
    CAT_VITAMINS,       // Vitamins
    CAT_NEUROTRANS,     // Neurotransmitters
    CAT_HORMONES,       // Steroid hormones
    CAT_PHARMA,         // Pharmaceuticals & drugs
    CAT_HOUSEHOLD,      // Household chemicals
    CAT_ACIDS,          // Acids
    CAT_PLASTICS,       // Plastic monomers/polymers
    CAT_FATS,           // Fatty acids
    CAT_METALS,         // Metal compounds
    CAT_ENERGY,         // Energy molecules (ATP, etc.)
    CAT_FLAVORS,        // Fragrances and flavors
    CAT_OTHER,          // Other / Random
    CAT_COUNT
};

// Get category name
const char* molecule_get_category_name(int category);

// Molecule builder function type
typedef void (*MoleculeBuilder)(Molecule* mol);

// Get the number of available molecule presets
int molecule_get_count();

// Get the name of a molecule preset by index
const char* molecule_get_name(int index);

// Get the category of a molecule preset by index
int molecule_get_category(int index);

// Get the description of a molecule preset by index
const char* molecule_get_description(int index);

// Build a molecule by index
void molecule_build(int index, Molecule* mol);

// Build a random molecule
void molecule_build_random(Molecule* mol);

#endif // MOLECULE_DB_H
