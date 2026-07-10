# CLAUDE.md

This file gives Claude Code repository-specific guidance for MolVis.

## Non-negotiable workflow rules

- Never commit, push, tag, open a pull request, publish a release, or upload artifacts. A human always performs those steps.
- Do not run scripts under `scripts/release/` or any packaging command that uploads to GitHub unless a human explicitly asks for that exact operation.
- Preserve unrelated local changes. Inspect `git status` before editing and do not overwrite or clean files you did not create.
- Make the smallest focused change that solves the request. Do not modify vendored code in `third_party/` unless the task specifically requires it.
- Update `RELEASE_NOTES.md` for user-visible code or behavior changes. Documentation-only and agent-instruction changes do not need a release-note entry unless requested.
- Report what changed and what validation was run; leave committing and publishing to the human.

## Project overview

MolVis is a C++17, real-time ball-and-stick molecular visualizer with a shared molecule database and two native GPU rendering paths:

- Windows: Win32 + DirectX 11 + Dear ImGui, with CUDA ray tracing.
- macOS: SDL2 + Metal + Dear ImGui, packaged as a native app bundle.

The application contains roughly 300 compiled presets plus optional external JSON records. Treat exact preset counts as derived data: use `molecule_get_count()` rather than hard-coding a count.

## Build and validation

### Windows

Requirements: Visual Studio with the C++ workload, CUDA Toolkit 12+, and an NVIDIA CUDA-capable GPU.

```powershell
.\build.bat          # Configure the MSVC environment and build
.\build.bat run      # Build and run
.\build.bat clean    # Remove Windows build outputs
nmake                # Build from an initialized VS Developer Prompt
```

The normal Windows output is `molvis.exe`. The Makefile and CMake configuration currently target CUDA architecture `sm_86`; changing GPU compatibility should be deliberate and coordinated in both build systems.

### macOS

Requirements: Xcode command-line tools and Metal toolchain, CMake 3.20+, and SDL2.

```bash
./build_mac.sh          # Configure and build
./build_mac.sh run      # Build if needed, then run
./build_mac.sh clean    # Remove build_mac/
```

The current bundle output is `build_mac/bin/Molecule Visualizer.app` (note the space). CMake compiles the Metal shader and bundles it and SDL2 into the application.

There is no standalone unit-test suite. Validate changes with the narrowest relevant checks and then build the affected platform when its toolchain is available. If the current machine cannot build a platform-specific path, say so clearly rather than claiming it was validated.

## Architecture and source map

```text
src/
  gui/
    main_windows.cpp                 Windows entry point and ImGui UI
    main_mac.mm                      macOS entry point and ImGui UI
  molecule/
    molecule_types.h                 Shared Atom, Bond, and Molecule structs
    molecule_db.h                    Public preset/category API
    molecule_db.cpp                  Preset builders, metadata, and registry
    molecule_loader.{h,cpp}          External JSON parsing, validation, and paths
  renderer/
    renderer_interface.h             Platform-neutral IRenderer contract
    cuda_renderer.{h,cu}             Windows CUDA renderer and kernels
    metal_renderer.{h,mm}            macOS Metal renderer
    shaders/MoleculeShaders.metal    Metal compute shader
platform/macos/                      Info.plist and app icon
scripts/package/                     Local distribution packaging
scripts/release/                     Release/upload automation; do not run autonomously
third_party/imgui/                   Vendored Dear ImGui sources and backends
legacy/                              Original implementation; not the normal app path
```

The frame flow is: platform GUI gathers input and selects a `Molecule`; the platform renderer receives molecule and camera state; GPU code intersects rays with atom spheres and bond cylinders, applies lighting, and exposes a texture for the ImGui viewport.

`IRenderer` documents the common lifecycle (`init`, `resize`, `render`, `cleanup`, texture access), but each GUI currently instantiates its platform renderer directly. Changes to shared rendering behavior usually require corresponding CUDA and Metal shader/renderer changes. Keep the two UIs behaviorally aligned when changing controls or panels.

## Shared data constraints

`src/molecule/molecule_types.h` is consumed by ordinary C++, CUDA, Objective-C++, and Metal-side data transfer code. Preserve field order and compatible layouts unless every consumer is updated together.

Important limits and conventions:

- `MAX_ATOMS` is 200 and `MAX_BONDS` is 250.
- Bond order is `1` for single, `2` for double, and `3` for triple.
- `Molecule` includes `name[64]` and `formula[32]` as well as atom/bond arrays and counts.
- Atom type IDs are shared constants from `ATOM_H` through `ATOM_GA`; `ATOM_TYPE_COUNT` is 25.
- Validate indices and capacity implications when adding atoms or bonds. The local builder helpers silently stop adding data at the fixed limits.

## Adding or changing molecule presets

Preset work belongs in `src/molecule/molecule_db.cpp`:

1. Follow a nearby `static void buildX(Molecule* mol)` builder as the template.
2. Initialize the molecule, set its name/formula as the surrounding code does, add atoms before bonds, and finish with `centerMolecule(mol)`.
3. Add or update the entry in the `molecules[]` registry with the correct category and description metadata.
4. Keep bond indices valid and atom/bond totals below the fixed limits.
5. Update user-facing counts or descriptions only when they actually changed; avoid introducing another independently maintained hard-coded count.
6. Build and visually inspect the preset when possible, checking centering, bond orders, clipping, and both hydrogen-visible and hydrogen-hidden views where relevant.

## Coding conventions

- Use the style already present in the file being changed; this codebase does not use an enforced formatter.
- Host-side functions generally use `camelCase`; CUDA/Metal math and device helpers commonly use `snake_case`.
- Types use `PascalCase`; constants and macros use `UPPER_SNAKE_CASE`.
- CUDA entry kernels use `__global__`, helpers use `__device__`, and read-only device tables use `__device__ __constant__` where appropriate.
- Use existing vector helpers (`normalize3`, `dot3`, `reflect3`, `length3`, and platform equivalents) instead of creating subtly different math paths.
- Keep platform-neutral types free of Win32, DirectX, SDL, Objective-C, CUDA-runtime, and Metal-specific dependencies.
- Avoid broad cleanup in the very large molecule database. Keep preset edits localized and easy to review.

## Configuration and generated artifacts

- Windows runtime settings use `molvis.ini` in the working directory.
- macOS runtime settings use `~/.config/molvis.ini`.
- `build/`, `build_mac/`, executables, app bundles, DMGs, ZIPs, and generated shader outputs are build artifacts; do not treat them as source changes.
- Versions and release metadata appear in multiple places, including `CMakeLists.txt`, `platform/macos/Info.plist`, packaging scripts, and release notes. If asked to bump a version, search the repository and update all relevant locations consistently.

For extended build, signing, packaging, and troubleshooting details, consult `docs/BUILDING.md`. Packaging or deployment instructions are reference material, not authorization to publish anything.

External molecule files are documented in `docs/MOLECULE_FORMAT.md`. Preserve stable IDs and compiled fallback behavior when changing the runtime registry. Crystal presets represent finite fragments of extended materials and must be described as such.

Always save plans, specs and bugfixes in the repository. Do not rely on local-only notes or external documents to convey design intent or implementation details.
