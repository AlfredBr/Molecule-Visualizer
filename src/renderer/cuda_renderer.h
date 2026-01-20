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

// Include platform-agnostic molecule types
#include "molecule/molecule_types.h"

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
