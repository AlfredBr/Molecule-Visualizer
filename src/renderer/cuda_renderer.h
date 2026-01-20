/*
 * MolVis - CUDA Renderer Header
 *
 * Defines the interface between the CUDA renderer and the ImGui application.
 * Uses DX11-CUDA interop for efficient texture sharing.
 */

#ifndef CUDA_RENDERER_H
#define CUDA_RENDERER_H

#include <d3d11.h>
#include <cuda_runtime.h>

// Molecule limits
#define MAX_ATOMS 200
#define MAX_BONDS 250

// Atom types (CPK coloring)
#define ATOM_H  0   // Hydrogen - white
#define ATOM_C  1   // Carbon - dark gray
#define ATOM_N  2   // Nitrogen - blue
#define ATOM_O  3   // Oxygen - red
#define ATOM_P  4   // Phosphorus - orange
#define ATOM_S  5   // Sulfur - yellow
#define ATOM_CL 6   // Chlorine - green
#define ATOM_BR 7   // Bromine - dark red
#define ATOM_F  8   // Fluorine - light green
#define ATOM_I  9   // Iodine - purple
#define ATOM_NA 10  // Sodium
#define ATOM_SI 11  // Silicon
#define ATOM_B  12  // Boron
#define ATOM_FE 13  // Iron
#define ATOM_CU 14  // Copper
#define ATOM_AL 15  // Aluminum
#define ATOM_TI 16  // Titanium
#define ATOM_PT 17  // Platinum

// Atom structure
struct Atom {
    float x, y, z;
    int type;
    float radius;
};

// Bond structure
struct Bond {
    int atom1, atom2;
    int order;  // 1=single, 2=double, 3=triple
};

// Molecule structure
struct Molecule {
    Atom atoms[MAX_ATOMS];
    Bond bonds[MAX_BONDS];
    int numAtoms;
    int numBonds;
    char name[64];
};

// Renderer state (opaque)
struct CudaRenderer;

// Initialize the CUDA renderer with DX11 interop
// Returns nullptr on failure
CudaRenderer* renderer_init(ID3D11Device* device, int width, int height);

// Cleanup renderer resources
void renderer_cleanup(CudaRenderer* renderer);

// Resize the render target
bool renderer_resize(CudaRenderer* renderer, int width, int height);

// Render a molecule to the internal texture
// rotX, rotY: camera rotation angles
// zoom: camera distance
// offsetX, offsetY: horizontal/vertical offset (default 0)
void renderer_render(CudaRenderer* renderer,
                     const Molecule* molecule,
                     float rotX, float rotY, float zoom,
                     float offsetX = 0.0f, float offsetY = 0.0f);

// Render text overlay on the current frame
// text: null-terminated string to render
// x, y: position in pixels from top-left
// scale: font scale (1 = 6x8 pixels per character)
void renderer_render_text(CudaRenderer* renderer,
                          const char* text, int x, int y, int scale);

// Get the DX11 shader resource view for ImGui display
ID3D11ShaderResourceView* renderer_get_texture(CudaRenderer* renderer);

// Get current render dimensions
void renderer_get_size(CudaRenderer* renderer, int* width, int* height);

#endif // CUDA_RENDERER_H
