# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

MolVis is a GPU-accelerated 3D molecular visualization application that renders ball-and-stick molecular models in real-time. It supports 300+ molecule presets across 17 categories with CPK-convention atomic coloring.

**Platform Support:**
- **Windows**: NVIDIA CUDA + DirectX 11 + Win32
- **macOS**: Apple Metal + SDL2

## Build Commands

### Windows
```powershell
.\build.bat           # Build (auto-detects Visual Studio)
.\build.bat run       # Build and run
.\build.bat clean     # Clean build artifacts
nmake                 # Manual build from VS Developer Command Prompt
```

### macOS
```bash
./build_mac.sh        # Build
./build_mac.sh run    # Build and run
./build_mac.sh clean  # Clean build artifacts
open build_mac/bin/MolVis.app  # Run the app
```

**Requirements:**
- Windows: Visual Studio 2019/2022/2025, CUDA Toolkit 12.0+, NVIDIA GPU (Compute Capability 5.0+)
- macOS: Xcode Command Line Tools, CMake, SDL2 (via Homebrew)

## Architecture

### Platform Abstraction Pattern

The renderer uses a Strategy pattern with platform-specific implementations:

```
IRenderer (src/renderer/renderer_interface.h)
├── CudaRenderer (src/renderer/cuda_renderer.cu)  - Windows implementation
└── MetalRenderer (src/renderer/metal_renderer.mm) - macOS implementation
```

### Key Data Structures (src/molecule/molecule_types.h)

```cpp
struct Atom { float x, y, z; int type; float radius; };
struct Bond { int atom1, atom2; int order; };  // order: 1=single, 2=double, 3=triple
struct Molecule { Atom atoms[200]; Bond bonds[250]; int numAtoms, numBonds; char name[64]; };
```

20 supported elements with CPK coloring: H, C, N, O, P, S, Cl, Br, F, I, Na, Si, B, Fe, Cu, Al, Ti, Pt, Re, Xe

### Source Layout

| Path | Description |
|------|-------------|
| `src/gui/main.cpp` | Windows entry point (Win32 + DX11 + ImGui) |
| `src/gui/main_mac.mm` | macOS entry point (SDL2 + Metal + ImGui) |
| `src/renderer/cuda_renderer.cu` | CUDA ray-tracing kernels (~25KB) |
| `src/renderer/metal_renderer.mm` | Metal implementation |
| `src/renderer/shaders/MoleculeShaders.metal` | Metal compute shaders |
| `src/molecule/molecule_db.cpp` | 300+ molecule presets (~16K lines) |
| `third_party/imgui/` | Dear ImGui library with all backends |

### Rendering Pipeline

1. GUI receives user input (rotation, zoom, molecule selection)
2. Camera parameters passed to platform renderer
3. GPU kernel performs per-pixel ray tracing:
   - Ray-sphere intersection for atoms
   - Ray-cylinder intersection for bonds (with offset for double/triple bonds)
   - Phong shading with Fresnel rim lighting
4. Rendered texture displayed in ImGui viewport

## Coding Conventions

### CUDA Kernels
- `__global__` for kernel entry points, `__device__` for helpers
- `__device__ __constant__` for read-only data (colors, radii)
- Prefix kernels with verbs: `render_`, `compute_`, `update_`

### Naming
- Functions: camelCase (host), snake_case (device)
- Structs: PascalCase (`Atom`, `Bond`, `Molecule`)
- Constants: UPPER_SNAKE_CASE, macros with `#define`

### Vector Math
- Use `float3` for 3D vectors
- Helper prefixes: `normalize3`, `dot3`, `reflect3`, `length3`

## Adding New Molecules

In `src/molecule/molecule_db.cpp`:

1. Create a builder function using the pattern:
```cpp
Molecule buildMoleculeName() {
    Molecule mol = {};
    strcpy(mol.name, "Molecule Name");
    // addAtom(mol, x, y, z, ATOM_TYPE);
    // addBond(mol, atom1_idx, atom2_idx, order);
    centerMolecule(mol);
    return mol;
}
```

2. Add to `getMoleculePresets()` in the appropriate category
3. Increment the molecule count in README.md if applicable

## Notes

- The Makefile targets sm_86 (RTX 3080) by default; adjust `NVCCFLAGS` for other GPUs
- Always update RELEASE_NOTES.md after changes
- ImGui docking is enabled for flexible panel layout
- Config is persisted: Windows uses `molvis.ini`, macOS uses `~/.config/molvis.ini`
