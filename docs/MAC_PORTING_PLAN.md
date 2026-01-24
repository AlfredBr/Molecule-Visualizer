# MolVis macOS Porting Plan

## Executive Summary

This document outlines a comprehensive plan to port MolVis from Windows (CUDA + DirectX 11 + Win32) to macOS. The primary challenge is that **CUDA is not supported on macOS** (Apple dropped NVIDIA support in 2019), requiring a complete replacement of the GPU compute and rendering stack.

---

## Current Architecture Analysis

### Windows Version Components

| Component | Technology | macOS Equivalent |
|-----------|------------|------------------|
| GPU Compute | NVIDIA CUDA | **Metal Compute Shaders** |
| Graphics API | DirectX 11 | **Metal** |
| Windowing | Win32 API | **Cocoa (AppKit)** or SDL2/GLFW |
| GUI Framework | Dear ImGui (Win32 + DX11 backends) | Dear ImGui (**macOS + Metal** backends) |
| Build System | nmake + MSVC + NVCC | **CMake + Clang** or Xcode |

### Files Requiring Changes

| File | Change Level | Notes |
|------|--------------|-------|
| `src/gui/main_windows.cpp` | **Major rewrite** | Replace Win32/DX11 with Cocoa/Metal |
| `src/renderer/cuda_renderer.cu` | **Full rewrite** | Convert CUDA kernels to Metal compute shaders |
| `src/renderer/cuda_renderer.h` | **Moderate** | Update API to use Metal types |
| `src/molecule/molecule_db.cpp` | **None** | Pure C++, fully portable |
| `src/molecule/molecule_db.h` | **Minor** | Remove DX11 include dependency |
| `Makefile` | **Replace** | Need CMake or Makefile for macOS |
| `third_party/imgui/` | **Update** | Add Metal + Cocoa backends |

---

## Porting Strategy

### Recommended Approach: Metal + SDL2/GLFW

I recommend using **SDL2 or GLFW** instead of raw Cocoa/AppKit because:
1. ImGui has excellent SDL2/GLFW backends that work on macOS
2. Simpler integration than raw Objective-C
3. Potential for cross-platform support later (Linux)
4. Well-documented and widely used

### Alternative: Pure Cocoa/Metal

If you want a more "native" macOS feel, you could use:
- Cocoa (AppKit) for windowing
- Metal for both rendering and compute
- ImGui Cocoa + Metal backends

---

## Phase 1: Project Setup & Build System (Week 1)

### 1.1 Create macOS Directory Structure

```
MolVis/
├── src/
│   ├── main_windows.cpp            # Windows entry point
│   ├── main_mac.mm                 # NEW: macOS entry point (Obj-C++)
│   ├── renderer/
│   │   ├── cuda_renderer.cu        # Keep for reference
│   │   ├── cuda_renderer.h         # Keep for reference
│   │   ├── metal_renderer.h        # NEW: Metal renderer interface
│   │   ├── metal_renderer.mm       # NEW: Metal implementation
│   │   └── shaders/
│   │       └── MoleculeShaders.metal  # NEW: Metal compute shaders
│   └── molecule/
│       ├── molecule_db.cpp         # Portable - no changes
│       └── molecule_db.h           # Minor changes
├── platform/
│   ├── windows/                    # Move Windows-specific files here
│   └── macos/                      # NEW: macOS-specific files
│       └── Info.plist              # NEW: macOS app bundle info
├── third_party/
│   └── imgui/                      # Update with Metal/SDL backends
├── CMakeLists.txt                  # NEW: Cross-platform build
└── build_mac.sh                    # NEW: macOS build script
```

### 1.2 Set Up CMake Build System

Create `CMakeLists.txt` for cross-platform builds:

```cmake
cmake_minimum_required(VERSION 3.20)
project(MolVis)

set(CMAKE_CXX_STANDARD 17)

if(APPLE)
    # macOS-specific settings
    set(CMAKE_OSX_DEPLOYMENT_TARGET "11.0")
    find_library(METAL_FRAMEWORK Metal)
    find_library(METALKIT_FRAMEWORK MetalKit)
    find_library(COCOA_FRAMEWORK Cocoa)
    find_library(QUARTZCORE_FRAMEWORK QuartzCore)

    # Option: Use SDL2 for windowing
    find_package(SDL2 REQUIRED)
endif()
```

### 1.3 Update Dear ImGui

Download/update ImGui to include these backend files:
- `imgui_impl_sdl2.cpp/.h` (or `imgui_impl_osx.mm` for Cocoa)
- `imgui_impl_metal.mm/.h`

---

## Phase 2: Platform Abstraction Layer (Week 1-2)

### 2.1 Create Renderer Interface

Create a platform-agnostic renderer interface that both CUDA and Metal can implement:

**`src/renderer/renderer_interface.h`**:
```cpp
#ifndef RENDERER_INTERFACE_H
#define RENDERER_INTERFACE_H

#include "molecule/molecule_db.h"

// Abstract renderer interface
class IRenderer {
public:
    virtual ~IRenderer() = default;

    virtual bool init(void* windowHandle, int width, int height) = 0;
    virtual void cleanup() = 0;
    virtual bool resize(int width, int height) = 0;
    virtual void render(const Molecule* mol, float rotX, float rotY, float zoom,
                        float offsetX, float offsetY) = 0;
    virtual void renderText(const char* text, int x, int y, int scale) = 0;
    virtual void* getTexture() = 0;  // Returns platform-specific texture handle
    virtual void getSize(int* width, int* height) = 0;
};

// Factory function (implemented per platform)
IRenderer* createRenderer();
```

### 2.2 Update molecule_db.h

Remove the DirectX dependency:

```cpp
// Before (Windows):
#include "renderer/cuda_renderer.h"  // Has <d3d11.h>

// After (Portable):
// Move Atom, Bond, Molecule structs to a shared header
#include "molecule_types.h"
```

---

## Phase 3: Metal Compute Shader Implementation (Week 2-3)

### 3.1 Convert CUDA Kernels to Metal

The CUDA renderer has these main kernels that need conversion:

| CUDA Kernel | Metal Equivalent | Complexity |
|-------------|------------------|------------|
| `clearKernel` | `kernel void clearPixels()` | Simple |
| `renderMoleculeKernel` | `kernel void renderMolecule()` | Complex |
| `renderTextKernel` | `kernel void renderText()` | Medium |

### 3.2 Example Metal Shader Conversion

**CUDA (`cuda_renderer.cu`):**
```cuda
__device__ float3 normalize3(float3 v) {
    float len = sqrtf(v.x * v.x + v.y * v.y + v.z * v.z);
    if (len > 0.0001f) {
        return make_float3(v.x / len, v.y / len, v.z / len);
    }
    return make_float3(0, 1, 0);
}

__global__ void clearKernel(unsigned char* pixels, int width, int height) {
    int x = blockIdx.x * blockDim.x + threadIdx.x;
    int y = blockIdx.y * blockDim.y + threadIdx.y;
    if (x >= width || y >= height) return;

    int idx = (y * width + x) * 4;
    float gy = (float)y / height;
    pixels[idx + 0] = (unsigned char)(20 + gy * 30);
    pixels[idx + 1] = (unsigned char)(25 + gy * 35);
    pixels[idx + 2] = (unsigned char)(35 + gy * 40);
    pixels[idx + 3] = 255;
}
```

**Metal (`MoleculeShaders.metal`):**
```metal
#include <metal_stdlib>
using namespace metal;

// Metal has built-in normalize() function

kernel void clearPixels(
    texture2d<half, access::write> output [[texture(0)]],
    uint2 gid [[thread_position_in_grid]])
{
    if (gid.x >= output.get_width() || gid.y >= output.get_height()) return;

    float gy = float(gid.y) / float(output.get_height());
    half4 color = half4(
        (20 + gy * 30) / 255.0h,
        (25 + gy * 35) / 255.0h,
        (35 + gy * 40) / 255.0h,
        1.0h
    );
    output.write(color, gid);
}
```

### 3.3 Metal Renderer Implementation Outline

**`src/renderer/metal_renderer.h`:**
```objc
#ifndef METAL_RENDERER_H
#define METAL_RENDERER_H

#include "renderer_interface.h"

#ifdef __OBJC__
#import <Metal/Metal.h>
#import <MetalKit/MetalKit.h>
#endif

class MetalRenderer : public IRenderer {
public:
    MetalRenderer();
    ~MetalRenderer() override;

    bool init(void* windowHandle, int width, int height) override;
    void cleanup() override;
    bool resize(int width, int height) override;
    void render(const Molecule* mol, float rotX, float rotY, float zoom,
                float offsetX, float offsetY) override;
    void renderText(const char* text, int x, int y, int scale) override;
    void* getTexture() override;
    void getSize(int* width, int* height) override;

private:
    struct Impl;
    Impl* pImpl;  // Pimpl pattern to hide Obj-C types from C++ headers
};

#endif
```

---

## Phase 4: Application Layer (Week 3-4)

### 4.1 Create macOS Entry Point

**`src/main_mac.mm`:**
```objc
#import <Cocoa/Cocoa.h>
#import <Metal/Metal.h>

#include "imgui.h"
#include "imgui_impl_osx.h"
#include "imgui_impl_metal.h"

#include "renderer/metal_renderer.h"
#include "molecule/molecule_db.h"

// ... (similar structure to main_windows.cpp but using Cocoa/Metal)
```

### 4.2 Alternative: SDL2 Entry Point

Using SDL2 makes the code more portable and easier to maintain:

**`src/main_sdl_metal.cpp`:**
```cpp
#include <SDL2/SDL.h>
#include <SDL2/SDL_metal.h>

#include "imgui.h"
#include "imgui_impl_sdl2.h"
#include "imgui_impl_metal.h"

#include "renderer/metal_renderer.h"
#include "molecule/molecule_db.h"

int main(int argc, char* argv[]) {
    SDL_Init(SDL_INIT_VIDEO);

    SDL_Window* window = SDL_CreateWindow(
        "MolVis - Molecule Visualizer",
        SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED,
        1280, 800,
        SDL_WINDOW_METAL | SDL_WINDOW_RESIZABLE | SDL_WINDOW_ALLOW_HIGHDPI
    );

    SDL_MetalView metalView = SDL_Metal_CreateView(window);
    // ... setup Metal device, ImGui, main loop
}
```

---

## Phase 5: Testing & Polish (Week 4-5)

### 5.1 Verify Feature Parity

| Feature | Windows | macOS | Status |
|---------|---------|-------|--------|
| Molecule rendering | ✓ | | To implement |
| Ball-and-stick model | ✓ | | To implement |
| Phong shading | ✓ | | To implement |
| Text overlay | ✓ | | To implement |
| Mouse rotation | ✓ | | To implement |
| Zoom control | ✓ | | To implement |
| Auto-rotation | ✓ | | To implement |
| Molecule browser | ✓ | | To implement |
| Window state save | ✓ | | To implement |

### 5.2 Performance Optimization

- Profile Metal compute shader performance
- Test on various Mac hardware (Intel, M1, M2, M3)
- Optimize threadgroup sizes for Apple Silicon

### 5.3 App Bundle Creation

Create a proper macOS `.app` bundle:
```
MolVis.app/
├── Contents/
│   ├── Info.plist
│   ├── MacOS/
│   │   └── MolVis
│   └── Resources/
│       └── AppIcon.icns
```

---

## Technical Considerations

### CUDA to Metal Mapping

| CUDA Concept | Metal Equivalent |
|--------------|------------------|
| `__global__` kernel | `kernel` function |
| `blockIdx.x * blockDim.x + threadIdx.x` | `thread_position_in_grid` |
| `__device__` function | Regular function in Metal shader |
| `__constant__` memory | `constant` address space |
| `cudaMalloc` / `cudaFree` | `MTLDevice.makeBuffer()` |
| `cudaMemcpy` | `MTLBuffer.contents()` + memcpy |
| Texture2D (DX11) | `MTLTexture` |
| SRV (Shader Resource View) | `MTLTexture` directly usable |

### Key Differences

1. **Thread Organization:**
   - CUDA: Blocks and threads (`dim3 blockSize(16, 16)`)
   - Metal: Threadgroups and threads (`MTLSize(16, 16, 1)`)

2. **Memory Model:**
   - CUDA: Explicit host/device separation
   - Metal: Unified memory on Apple Silicon, shared memory modes

3. **Synchronization:**
   - CUDA: `cudaDeviceSynchronize()`
   - Metal: Command buffer completion handlers or `waitUntilCompleted`

---

## Estimated Timeline

| Phase | Duration | Deliverable |
|-------|----------|-------------|
| Phase 1: Setup | 1 week | Build system, project structure |
| Phase 2: Abstraction | 1 week | Platform-agnostic interfaces |
| Phase 3: Metal Shaders | 2 weeks | Compute shaders, renderer |
| Phase 4: Application | 1 week | Working macOS app |
| Phase 5: Polish | 1 week | Testing, optimization, packaging |

**Total: 5-6 weeks**

---

## Dependencies to Install on macOS

```bash
# Xcode Command Line Tools (required)
xcode-select --install

# CMake
brew install cmake

# SDL2 (optional, recommended)
brew install sdl2

# Download ImGui
# git clone https://github.com/ocornut/imgui.git third_party/imgui
```

---

## Recommended First Steps

1. **Set up the build environment:**
   ```bash
   xcode-select --install
   brew install cmake sdl2
   ```

2. **Create CMakeLists.txt** with macOS support

3. **Download/update ImGui** with Metal and SDL2 backends

4. **Start with a minimal Metal "hello world"** before porting the full renderer

5. **Port the molecule database first** (it's already portable)

6. **Incrementally port the CUDA kernels** to Metal, testing each one

---

## Questions to Consider

1. **Do you want Apple Silicon (M1/M2/M3) optimization?**
   - Metal is already optimized for Apple Silicon
   - Consider using `MTLGPUFamily` for feature detection

2. **Do you need to support Intel Macs?**
   - Minimum macOS version affects Metal features available

3. **Should this become truly cross-platform (Linux too)?**
   - Consider Vulkan + compute shaders for maximum portability
   - Or use SDL2 abstraction throughout

4. **App Store distribution?**
   - Requires proper code signing, sandboxing, notarization

---

## References

- [Metal Best Practices Guide](https://developer.apple.com/library/archive/documentation/3DDrawing/Conceptual/MTLBestPracticesGuide/)
- [Metal Shading Language Specification](https://developer.apple.com/metal/Metal-Shading-Language-Specification.pdf)
- [Dear ImGui Metal Backend](https://github.com/ocornut/imgui/blob/master/backends/imgui_impl_metal.mm)
- [SDL2 Metal Support](https://wiki.libsdl.org/SDL2/CategoryMetal)
