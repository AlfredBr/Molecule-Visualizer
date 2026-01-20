# MolVis 🧬

**GPU-Accelerated Molecular Visualization**

MolVis is a real-time 3D molecular visualization application powered by NVIDIA CUDA. It renders beautiful ball-and-stick molecular models with realistic atomic colors, metallic shading, and smooth camera controls.

![MolVis Screenshot](docs/screenshot.png) <!-- TODO: Add screenshot -->

## ✨ Features

- **Real-time GPU Rendering** - Leverages CUDA for fast, parallel ray-sphere intersection
- **CPK Color Convention** - Industry-standard atomic coloring (Carbon=gray, Oxygen=red, etc.)
- **Ball-and-Stick Models** - Clear visualization of molecular structure and bonds
- **Multiple Molecule Presets** - Organic compounds, DNA bases, and more
- **Random Molecule Generation** - Explore procedurally generated structures
- **Interactive Camera** - Rotate, zoom, and auto-rotate views
- **Metallic Shading** - Specular highlights for depth perception

## 🎮 Controls

| Key | Action |
|-----|--------|
| `1-9` | Select molecule preset |
| `R` | Generate random molecule |
| `Arrow Keys` | Rotate view |
| `W/S` | Zoom in/out |
| `A` | Toggle auto-rotate |
| `Space` | Pause/resume rotation |
| `Q` / `Escape` | Quit |

## 🛠️ Requirements

### Hardware
- NVIDIA GPU with Compute Capability 5.0+ (GTX 900 series or newer)
- Recommended: RTX 3080 or equivalent

### Software
- Windows 10/11
- [NVIDIA CUDA Toolkit](https://developer.nvidia.com/cuda-downloads) (12.0+)
- [Visual Studio](https://visualstudio.microsoft.com/) 2019/2022/2025 with "Desktop development with C++" workload
- GNU Make (included with Git for Windows or install separately)

## 🚀 Quick Start

### Building

```powershell
# Clone the repository
git clone https://github.com/yourusername/MolVis.git
cd MolVis

# Build using the automated script (recommended)
.\build.bat

# Or manually set up environment and build
. .\setup_env.ps1
nmake
```

### Running

```powershell
.\molvis.exe
```

## 📁 Project Structure

```
MolVis/
├── cuda_molecule.cu    # Main CUDA application
├── win32_display.h     # Win32 windowing abstraction
├── Makefile            # Build configuration
├── build.bat           # Automated build script
├── setup_env.ps1       # PowerShell environment setup
└── README.md           # This file
```

## 🔧 Configuration

### GPU Architecture

The default target is RTX 3080 (sm_86). To change for your GPU, edit the `Makefile`:

```makefile
# Common architectures:
# sm_50 - GTX 900 series (Maxwell)
# sm_61 - GTX 1000 series (Pascal)
# sm_75 - RTX 2000 series (Turing)
# sm_86 - RTX 3000 series (Ampere)
# sm_89 - RTX 4000 series (Ada Lovelace)
NVCCFLAGS = -O3 -arch=sm_86 -allow-unsupported-compiler
```

## 🗺️ Roadmap

- [ ] **Dear ImGui Integration** - Modern GUI with dockable panels
- [ ] **Molecule File Import** - Load PDB, MOL2, XYZ formats
- [ ] **Advanced Rendering** - Ambient occlusion, depth of field
- [ ] **Measurement Tools** - Bond lengths, angles, dihedral angles
- [ ] **Animation** - Molecular dynamics playback
- [ ] **Export** - Screenshot and video capture

## 🧪 Supported Elements

MolVis supports visualization of 18 elements with CPK coloring:

| Element | Symbol | Color |
|---------|--------|-------|
| Hydrogen | H | White |
| Carbon | C | Dark Gray |
| Nitrogen | N | Blue |
| Oxygen | O | Red |
| Phosphorus | P | Orange |
| Sulfur | S | Yellow |
| Chlorine | Cl | Green |
| Bromine | Br | Dark Red |
| Fluorine | F | Light Green |
| Iodine | I | Purple |
| Sodium | Na | Metallic Purple |
| Silicon | Si | Tan |
| Boron | B | Salmon |
| Iron | Fe | Orange/Brown |
| Copper | Cu | Copper |
| Aluminum | Al | Silver |
| Titanium | Ti | Gray |
| Platinum | Pt | White/Silver |

## 🤝 Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

1. Fork the repository
2. Create your feature branch (`git checkout -b feature/AmazingFeature`)
3. Commit your changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to the branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 🙏 Acknowledgments

- CPK coloring convention by Corey, Pauling, and Koltun
- NVIDIA for the CUDA toolkit and documentation
- Dear ImGui by Omar Cornut (upcoming integration)

---

Made with 💜 and CUDA and Claude Opus 4.5
