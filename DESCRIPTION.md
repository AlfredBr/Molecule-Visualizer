 MolVis is a cross-platform, GPU-accelerated 3D molecular visualization application that renders interactive ball-and-stick
   molecular models in real-time.

  Purpose

  - Scientific education and visualization of molecular structures
  - Supports 233 molecules across 17 categories (simple gases, pharmaceuticals, amino acids, vitamins, neurotransmitters,
  and more)

  Architecture

  src/
  ├── gui/
  │   ├── main_windows.cpp  # Windows entry point (Win32 + DirectX 11)
  │   └── main_mac.mm       # macOS entry point (SDL2 + Metal)
  ├── molecule/
  │   ├── molecule_types.h  # Platform-agnostic data structures
  │   ├── molecule_db.h     # Database API
  │   └── molecule_db.cpp   # 233 molecule presets (~15,800 lines)
  └── renderer/
      ├── renderer_interface.h   # Abstract renderer interface
      ├── cuda_renderer.cu       # CUDA ray-tracing (Windows)
      ├── metal_renderer.mm      # Metal implementation (macOS)
      └── shaders/
          └── MoleculeShaders.metal

  Technologies
  ┌──────────┬─────────────┬───────────┬──────────────┐
  │ Platform │ GPU Compute │ Windowing │ Graphics API │
  ├──────────┼─────────────┼───────────┼──────────────┤
  │ Windows  │ CUDA 12.0+  │ Win32     │ DirectX 11   │
  ├──────────┼─────────────┼───────────┼──────────────┤
  │ macOS    │ Metal       │ SDL2      │ Metal        │
  └──────────┴─────────────┴───────────┴──────────────┘
  - GUI: Dear ImGui with docking support
  - Build: CMake, Makefile, and platform-specific scripts
  - Language: C++17 (Objective-C++ for macOS)

  Key Features

  - GPU Ray Tracing: Per-pixel parallel computation for atoms (spheres) and bonds (cylinders)
  - CPK Color Convention: Industry-standard atomic coloring
  - Multi-bond rendering: Single, double, and triple bonds with visual separation
  - Phong shading with Fresnel rim lighting and dual light sources
  - Interactive controls: Real-time rotation, zoom, panning
  - Periodic table panel: Visual element reference with CPK colors
  - Persistent config: Window state saved between sessions

  Molecule Categories

  Amino acids, nucleobases, sugars, lipids, vitamins, neurotransmitters, hormones, pharmaceuticals, household chemicals,
  acids, metal compounds, polymers, ATP/energy molecules, flavors, and exotic structures (C₆₀ fullerene, catenanes,
  molecular knots).

  Current Status

  - Version: v0.3.8
  - License: MIT
  - Active development with recent macOS improvements and molecule additions
  