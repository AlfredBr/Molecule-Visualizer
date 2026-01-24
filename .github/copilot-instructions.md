# MolVis - GitHub Copilot Instructions

## Project Overview

MolVis is a GPU-accelerated molecular visualization application built with CUDA and targeting Windows. The application renders 3D ball-and-stick molecular models with realistic atomic colors (CPK convention), metallic shading, and interactive camera controls.

## Claude API Extended Thinking Requirements

When using Claude API with extended thinking enabled:

1. **Assistant message structure**: Every assistant message MUST start with a thinking block before any other content
2. **Valid block order**: `thinking` → `text` → `tool_use` (if applicable)
3. **Include previous thinking**: When continuing conversations, include thinking blocks from previous turns
4. **Alternative**: If you don't need extended thinking, disable it in your API request

### Code Implementation Notes:
- Check if `thinking` feature is enabled in your request
- If enabled, ensure assistant messages have structure: `[{type: "thinking", ...}, {type: "text", ...}]`
- If you're getting this error, either:
  - Restructure messages to include thinking blocks first, OR
  - Remove the thinking parameter from your API request


## Technology Stack

- **Language**: CUDA C++ (.cu files)
- **GPU Framework**: NVIDIA CUDA Runtime API
- **Target GPU**: RTX 3080 (Compute Capability 8.6) - configurable via Makefile
- **Platform**: Windows only (Win32 API + DirectX 11)
- **GUI**: Dear ImGui (Win32 + DX11 backends)
- **Build System**: nmake with NVCC + MSVC host compiler

## Architecture

### Project Structure
```
MolVis/
├── src/
│   ├── main_windows.cpp      # Application entry point (ImGui + Win32 + DX11)
│   ├── renderer/
│   │   ├── cuda_renderer.cu  # CUDA rendering kernels
│   │   └── cuda_renderer.h   # Renderer interface
│   └── molecule/
│       ├── molecule_db.h     # Molecule data structures
│       └── molecule_db.cpp   # Molecule presets and generation
├── third_party/
│   └── imgui/                # Dear ImGui library
├── legacy/
│   ├── cuda_molecule.cu      # Original CUDA-only version
│   └── win32_display.h       # Legacy Win32 abstraction
├── build/                    # Build artifacts (.obj files)
├── include/                  # Shared headers (future)
├── Makefile                  # nmake build configuration
├── build.bat                 # Automated build script (sets up MSVC)
└── setup_env.ps1             # PowerShell environment setup
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
- Always update the RELEASE_NOTES.md after every AI update

## Testing

- Test on various NVIDIA GPU architectures (sm_75, sm_86, sm_89)
- Verify molecule rendering accuracy against reference images
- Profile with NVIDIA Nsight for performance optimization
- Target platform: Windows 10/11 only (no Linux/macOS support)
