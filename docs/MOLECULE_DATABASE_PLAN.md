# Molecule Database Modernization Plan

## Goal

Keep MolVis self-contained with compiled presets while allowing users to add or explicitly override molecules after compilation through versioned JSON files.

## Design

1. Put the existing compiled registry behind a runtime database while preserving the current `molecule_get_*` and `molecule_build` API.
2. Give every entry a stable ID. Preserve selection by ID when the database reloads.
3. Load one molecule per JSON file from an application-adjacent `molecules/` directory and a platform user-data directory.
4. Merge in this order: compiled presets, application files, user files. Require `"override": true` to replace an existing ID. A bad override leaves the lower-precedence entry available.
5. Validate schema version, IDs, elements, finite coordinates, fixed atom/bond limits, bond indices, duplicate/self bonds, and bond orders. Reject one bad file without blocking startup.
6. Add host-only tests for valid parsing, malformed records, limits, overrides, ordering, and fallback behavior.
7. Add UI commands to reload the database and open the user molecule folder, plus source/error reporting.

## Schema v1

Required fields are `schemaVersion`, `id`, `name`, `category`, and `atoms`. Optional fields include `formula`, `description`, `longDescription`, `structureType`, `override`, and `bonds`. Positions are Cartesian coordinates in angstroms; bond indices are zero-based.

See `docs/MOLECULE_FORMAT.md` for the user-facing specification and examples.

## Semiconductor material presets

Add a **Materials & Semiconductors** category with finite, explicitly labeled crystal models:

- 2H-MoS2, 2H-WS2, and 2H-WSe2 monolayer fragments
- diamond-cubic silicon fragment
- representative diamond-lattice SiGe alloy fragment
- 3C-SiC fragment
- wurtzite GaN fragment

These are extended solids, not isolated molecules. Their descriptions and `structureType` must make the finite-model limitation clear. Supporting them requires Mo, W, Se, Ge, and Ga in shared atom constants, radii, both GPU color tables, both periodic-table UIs, and formula parsing.

## Delivery sequence

1. Runtime registry, stable IDs, loader, validation, and tests
2. Element expansion across shared data, CUDA, Metal, and UI
3. Crystal-building helpers and the seven built-in presets
4. Reload/open-folder UI, source diagnostics, examples, and documentation
5. Windows build, host tests, and macOS build/visual verification when each toolchain is available

Publishing, committing, pushing, tagging, and releasing remain human-only operations.
