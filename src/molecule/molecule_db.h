/*
 * MolVis - Molecule Database
 *
 * Contains molecule presets and generation functions.
 */

#ifndef MOLECULE_DB_H
#define MOLECULE_DB_H

#include "renderer/cuda_renderer.h"

// Molecule categories for UI organization
enum MoleculeCategory {
    CAT_SIMPLE = 0,     // Simple molecules (Water, Ammonia, etc.)
    CAT_ORGANIC,        // Organic compounds (Ethanol, Benzene, etc.)
    CAT_BIOCHEM,        // Biochemistry (ATP, Glucose, etc.)
    CAT_DRUGS,          // Pharmaceuticals (Aspirin, Caffeine, etc.)
    CAT_HOUSEHOLD,      // Household chemicals (Salt, Baking Soda, etc.)
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
