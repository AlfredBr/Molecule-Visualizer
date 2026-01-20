# MolVis 🧬

**GPU-Accelerated Molecular Visualization**

MolVis is a real-time 3D molecular visualization application powered by NVIDIA CUDA and Dear ImGui. It renders beautiful ball-and-stick molecular models with realistic atomic colors, metallic shading, and interactive controls through a modern GUI.

![MolVis Screenshot](docs/screenshot.png)

## ✨ Features

- **Real-time GPU Rendering** — CUDA-powered parallel ray-sphere intersection for smooth performance
- **Dear ImGui Interface** — Modern, responsive GUI with dockable panels
- **249 Molecule Library** — From simple gases to exotic structures like Buckminsterfullerene (C₆₀)
- **16 Categories** — Organized by type: Organic, Pharma, Amino Acids, Sugars, Fats, Metals, and more
- **CPK Color Convention** — Industry-standard atomic coloring (Carbon=gray, Oxygen=red, Nitrogen=blue, etc.)
- **Ball-and-Stick Models** — Clear visualization of molecular structure with single, double, and triple bonds
- **Interactive Controls** — Rotation speed, direction, position offset, zoom—all from the GUI
- **Metallic Shading** — Specular highlights and realistic lighting for depth perception
- **Window Persistence** — Remembers your window size and position between sessions

## 🖥️ Interface

### Panels

| Panel | Description |
|-------|-------------|
| **Molecule Selector** | Browse 249 molecules by category or search by name |
| **View Controls** | Adjust rotation speed, direction, and position offset |
| **Atom Inspector** | View atom count, bond count, and molecular details |
| **Render Settings** | Adjust lighting and display options |
| **Performance** | Real-time FPS and GPU statistics |

### Categories

- Simple Molecules (H₂O, CO₂, NH₃...)
- Organic Compounds (Benzene, Ethanol, Caffeine...)
- Amino Acids (all 20 standard)
- Nucleobases (DNA/RNA bases)
- Sugars & Carbohydrates
- Lipids & Fatty Acids
- Pharmaceuticals (Aspirin, Ibuprofen, Morphine...)
- Neurotransmitters (Dopamine, Serotonin, GABA...)
- Metal Compounds (Ferrocene, Cisplatin...)
- Plastics & Polymers
- Energy Molecules (ATP, NADH...)
- Exotic Structures (C₆₀, Cubane, Catenanes, Molecular Knots...)

## 🛠️ Requirements

- **OS**: Windows 10/11
- **GPU**: NVIDIA GPU with CUDA support (Compute Capability 5.0+)
- **CUDA Toolkit**: 12.0 or newer
- **Compiler**: Visual Studio 2019/2022/2025 with C++ workload

## 🧪 Supported Elements

18 elements with CPK coloring:

| Element | Color | Element | Color |
|---------|-------|---------|-------|
| Hydrogen (H) | White | Sodium (Na) | Purple |
| Carbon (C) | Dark Gray | Silicon (Si) | Tan |
| Nitrogen (N) | Blue | Boron (B) | Salmon |
| Oxygen (O) | Red | Iron (Fe) | Orange |
| Phosphorus (P) | Orange | Copper (Cu) | Copper |
| Sulfur (S) | Yellow | Aluminum (Al) | Silver |
| Chlorine (Cl) | Green | Titanium (Ti) | Gray |
| Bromine (Br) | Dark Red | Platinum (Pt) | Silver |
| Fluorine (F) | Light Green | Iodine (I) | Purple |

## 📁 Project Structure

```
MolVis/
├── src/
│   ├── main.cpp              # ImGui application + Win32/DX11
│   ├── renderer/
│   │   ├── cuda_renderer.h   # Renderer interface & atom types
│   │   └── cuda_renderer.cu  # CUDA ray-tracing kernels
│   └── molecule/
│       ├── molecule_db.h     # Molecule database API
│       └── molecule_db.cpp   # 249 molecule definitions
├── third_party/
│   └── imgui/                # Dear ImGui (Win32 + DX11 backends)
├── Makefile                  # nmake build configuration
└── build.bat                 # One-click build script
```

## 🗺️ Roadmap

- [x] Dear ImGui integration with dockable panels
- [x] 249 molecule library across 16 categories
- [x] Search and filter functionality
- [x] View controls (rotation, position, zoom)
- [x] Window state persistence
- [ ] Molecule file import (PDB, MOL2, XYZ)
- [ ] Measurement tools (bond lengths, angles)
- [ ] Screenshot/video export
- [ ] Advanced rendering (ambient occlusion, depth of field)

## 📄 License

MIT License — see [LICENSE](LICENSE) for details.

## 🙏 Acknowledgments

- **CPK coloring** — Corey, Pauling, and Koltun convention
- **Dear ImGui** — Omar Cornut's immediate-mode GUI
- **NVIDIA CUDA** — GPU computing platform

---

Made with 💜 CUDA and Claude Opus 4.5
