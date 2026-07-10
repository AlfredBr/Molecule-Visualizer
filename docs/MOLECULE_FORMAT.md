# External Molecule Files

MolVis keeps its built-in database compiled into the application and also loads user-editable JSON files. Put one molecule in each `.json` file and choose **File > Reload External Molecules** after editing.

Use **File > Open Molecules Folder** to open the platform user directory:

- Windows: `%LOCALAPPDATA%\MolVis\molecules\`
- macOS: `~/Library/Application Support/MolVis/molecules/`

MolVis also scans a `molecules/` folder beside the executable before scanning the user directory. User-directory files therefore have higher precedence.

## Schema version 1

```json
{
  "schemaVersion": 1,
  "id": "water-example",
  "name": "Water Example",
  "formula": "H2O",
  "category": "Simple Molecules",
  "structureType": "molecule",
  "description": "An editable water model.",
  "longDescription": "An example external molecule file.",
  "atoms": [
    { "element": "O", "position": [0.0, 0.0, 0.0] },
    { "element": "H", "position": [0.76, 0.59, 0.0] },
    { "element": "H", "position": [-0.76, 0.59, 0.0] }
  ],
  "bonds": [
    { "atoms": [0, 1], "order": 1 },
    { "atoms": [0, 2], "order": 1 }
  ]
}
```

Required fields are `schemaVersion`, `id`, `name`, `category`, and a nonempty `atoms` array. IDs may contain lowercase ASCII letters, digits, and hyphens. Positions are Cartesian coordinates in angstroms. Bond indices are zero-based and orders are 1, 2, or 3.

Optional fields are `formula`, `description`, `longDescription`, `structureType`, `override`, and `bonds`. Set `"override": true` to intentionally replace a lower-precedence record having the same ID. An invalid override is rejected and leaves the earlier record available.

Supported categories must exactly match a category shown in MolVis. Supported element symbols are H, C, N, O, P, S, Cl, Br, F, I, Na, Si, B, Fe, Cu, Al, Ti, Pt, Re, Xe, Mo, W, Se, Ge, and Ga.

Files exceeding 200 atoms or 250 bonds are rejected. MolVis also rejects unknown elements, invalid coordinates, bad bond indices or orders, self-bonds, and duplicate bonds. One bad file does not stop other files or compiled presets from loading.
