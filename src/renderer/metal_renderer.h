/*
 * MolVis - Metal Renderer Header
 *
 * Defines the interface for the Metal-based molecular renderer (macOS).
 * Uses Metal compute shaders for GPU-accelerated ray tracing.
 */

#ifndef METAL_RENDERER_H
#define METAL_RENDERER_H

#include "renderer_interface.h"

// Forward declarations for Metal types (avoid importing Objective-C in headers)
#ifdef __OBJC__
@protocol MTLDevice;
@protocol MTLTexture;
#endif

// Opaque renderer state
struct MetalRendererImpl;

class MetalRenderer : public IRenderer {
public:
    MetalRenderer();
    ~MetalRenderer() override;

    // IRenderer interface implementation
    bool init(void* windowHandle, int width, int height) override;
    void cleanup() override;
    bool resize(int width, int height) override;
    void render(const Molecule* mol, float rotX, float rotY, float zoom,
                float offsetX = 0.0f, float offsetY = 0.0f) override;
    void renderText(const char* text, int x, int y, int scale) override;
    void* getTexture() override;
    void getSize(int* width, int* height) override;

    // Metal-specific accessors
#ifdef __OBJC__
    id<MTLDevice> getDevice();
    id<MTLTexture> getMTLTexture();
#endif

private:
    MetalRendererImpl* pImpl;
};

// Factory function
IRenderer* createMetalRenderer();

#endif // METAL_RENDERER_H
