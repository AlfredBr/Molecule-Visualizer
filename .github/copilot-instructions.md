# MolVis - GitHub Copilot Instructions

## Project Overview

MolVis is a GPU-accelerated molecular visualization application built with CUDA and targeting Windows. The application renders 3D ball-and-stick molecular models with realistic atomic colors (CPK convention), metallic shading, and interactive camera controls.

## Technology Stack

- **Language**: CUDA C++ (.cu files)
- **GPU Framework**: NVIDIA CUDA Runtime API
- **Target GPU**: RTX 3080 (Compute Capability 8.6) - configurable via Makefile
- **Platform**: Windows (Win32 API for windowing)
- **Planned GUI**: Dear ImGui
- **Build System**: GNU Make with NVCC + MSVC host compiler

## Architecture

### Current Structure
- `cuda_molecule.cu` - Main application with CUDA kernels and rendering logic
- `win32_display.h` - Win32 abstraction layer for window management and input
- `Makefile` - Build configuration
- `build.bat` - Automated build script (sets up MSVC environment)
- `setup_env.ps1` - PowerShell environment setup script

### Planned Structure (with Dear ImGui)
```
MolVis/
├── src/
│   ├── main.cu              # Application entry point
│   ├── renderer/
│   │   ├── cuda_renderer.cu # CUDA rendering kernels
│   │   └── cuda_renderer.h
│   ├── molecule/
│   │   ├── molecule.h       # Molecule data structures
│   │   └── molecule_db.cu   # Molecule presets and generation
│   └── gui/
│       ├── gui.cpp          # Dear ImGui integration
│       └── gui.h
├── third_party/
│   └── imgui/               # Dear ImGui library
├── include/
│   └── win32_display.h
└── assets/
    └── molecules/           # Molecule data files (future)
```

## Coding Conventions

### CUDA Kernels
- Use `__global__` for kernel entry points
- Use `__device__` for device-only helper functions
- Use `__device__ __constant__` for read-only data (e.g., atom colors, radii)
- Prefix kernel names with descriptive verbs: `render_`, `compute_`, `update_`
- Always check CUDA errors in host code

### Naming Conventions
- **Functions**: camelCase for host, snake_case for device helpers
- **Structs**: PascalCase (e.g., `Atom`, `Bond`, `Molecule`)
- **Constants**: UPPER_SNAKE_CASE with `#define`
- **Device constants**: Include meaningful prefixes (e.g., `atomColors`, `atomRadii`)

### Memory Management
- Prefer `cudaMalloc`/`cudaFree` for device memory
- Use pinned memory (`cudaMallocHost`) for frequent host-device transfers
- Always pair allocations with deallocations
- Document memory ownership in comments

### Math Helpers
- Use `float3` for 3D vectors
- Prefix vector math functions: `normalize3`, `dot3`, `reflect3`, `length3`
- Use `__device__` qualifier for all GPU math helpers

## Dear ImGui Integration Guidelines

When integrating Dear ImGui:

1. **Backend**: Use `imgui_impl_win32` + `imgui_impl_dx11` backends
2. **Rendering**: Keep CUDA rendering separate, blit to ImGui texture
3. **Layout**: Use dockable windows for flexibility
4. **Panels to implement**:
   - Molecule selector/browser
   - Atom/bond property inspector
   - Rendering settings (lighting, colors, quality)
   - Camera controls
   - Performance metrics (FPS, GPU memory)

## Build Instructions

```powershell
# Option 1: Use build script (recommended)
.\build.bat

# Option 2: Set up environment manually
. .\setup_env.ps1
nmake

# Run the application
.\molvis.exe
```

## Key Constants

| Constant | Value | Description |
|----------|-------|-------------|
| `WIDTH` | 1024 | Default window width |
| `HEIGHT` | 768 | Default window height |
| `MAX_ATOMS` | 200 | Maximum atoms per molecule |
| `MAX_BONDS` | 250 | Maximum bonds per molecule |

## Atom Types (CPK Convention)

| ID | Element | Color |
|----|---------|-------|
| 0 | Hydrogen (H) | White |
| 1 | Carbon (C) | Dark Gray |
| 2 | Nitrogen (N) | Blue |
| 3 | Oxygen (O) | Red |
| 4 | Phosphorus (P) | Orange |
| 5 | Sulfur (S) | Yellow |
| ... | ... | ... |

## Performance Considerations

- Target 60 FPS at 1024x768 resolution
- Use shared memory for frequently accessed data in kernels
- Minimize host-device memory transfers
- Consider using CUDA streams for async operations when adding GUI

## Testing

- Test on multiple GPU architectures when possible
- Verify molecule rendering accuracy against reference images
- Profile with NVIDIA Nsight for performance optimization
