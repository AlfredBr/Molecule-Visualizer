/*
 * MolVis - Renderer Interface
 *
 * Platform-agnostic renderer interface.
 * Implemented by CudaRenderer (Windows) and MetalRenderer (macOS).
 */

#ifndef RENDERER_INTERFACE_H
#define RENDERER_INTERFACE_H

#include "molecule/molecule_types.h"

// Abstract renderer interface
class IRenderer {
public:
    virtual ~IRenderer() = default;

    // Initialize the renderer with the given window and dimensions
    virtual bool init(void* windowHandle, int width, int height) = 0;

    // Cleanup all renderer resources
    virtual void cleanup() = 0;

    // Resize the render target
    virtual bool resize(int width, int height) = 0;

    // Render a molecule with the given camera parameters
    virtual void render(const Molecule* mol, float rotX, float rotY, float zoom,
                        float offsetX = 0.0f, float offsetY = 0.0f) = 0;

    // Render text overlay on the current frame
    virtual void renderText(const char* text, int x, int y, int scale) = 0;

    // Get the texture for ImGui display (platform-specific type)
    virtual void* getTexture() = 0;

    // Get current render dimensions
    virtual void getSize(int* width, int* height) = 0;
};

#endif // RENDERER_INTERFACE_H
