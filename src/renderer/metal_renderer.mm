/*
 * MolVis - Metal Renderer Implementation
 *
 * GPU-accelerated molecular rendering using Metal compute shaders.
 * This is the macOS equivalent of cuda_renderer.cu.
 */

#import <Metal/Metal.h>
#import <MetalKit/MetalKit.h>
#import <Foundation/Foundation.h>
#include <stdio.h>
#include <string.h>

#include "metal_renderer.h"

// ============== Shader Parameter Structures ==============
// Must match the structures defined in MoleculeShaders.metal

struct RenderParams {
    int width;
    int height;
    int numAtoms;
    int numBonds;
    float rotX;
    float rotY;
    float zoom;
    float offsetX;
    float offsetY;
};

struct TextParams {
    int textLen;
    int startX;
    int startY;
    int scale;
    int width;
    int height;
};

// ============== Metal Renderer Implementation ==============

struct MetalRendererImpl {
    // Metal device and command queue
    id<MTLDevice> device;
    id<MTLCommandQueue> commandQueue;

    // Compute pipelines
    id<MTLComputePipelineState> clearPipeline;
    id<MTLComputePipelineState> renderMoleculePipeline;
    id<MTLComputePipelineState> renderTextPipeline;

    // Render target texture
    id<MTLTexture> renderTexture;

    // Buffers for molecule data
    id<MTLBuffer> atomBuffer;
    id<MTLBuffer> bondBuffer;
    id<MTLBuffer> paramsBuffer;
    id<MTLBuffer> textBuffer;
    id<MTLBuffer> textParamsBuffer;

    // Dimensions
    int width;
    int height;

    // Thread group size for compute shaders
    MTLSize threadGroupSize;
};

// ============== Helper Functions ==============

static bool createComputePipelines(MetalRendererImpl* impl) {
    NSError* error = nil;

    // Try to load the default Metal library from the app bundle
    id<MTLLibrary> library = [impl->device newDefaultLibrary];

    if (!library) {
        // If running outside a bundle, try to load from a metallib file
        NSString* libraryPath = [[NSBundle mainBundle] pathForResource:@"default" ofType:@"metallib"];
        if (libraryPath) {
            NSURL* libraryURL = [NSURL fileURLWithPath:libraryPath];
            library = [impl->device newLibraryWithURL:libraryURL error:&error];
        }

        if (!library) {
            // Last resort: compile from source
            NSString* shaderPath = @"src/renderer/shaders/MoleculeShaders.metal";

            // Try multiple paths
            NSArray* searchPaths = @[
                shaderPath,
                @"../src/renderer/shaders/MoleculeShaders.metal",
                @"MoleculeShaders.metal",
                [[NSBundle mainBundle] pathForResource:@"MoleculeShaders" ofType:@"metal"] ?: @""
            ];

            NSString* shaderSource = nil;
            for (NSString* path in searchPaths) {
                if (path.length > 0) {
                    shaderSource = [NSString stringWithContentsOfFile:path
                                                             encoding:NSUTF8StringEncoding
                                                                error:&error];
                    if (shaderSource) {
                        printf("Loaded Metal shaders from: %s\n", [path UTF8String]);
                        break;
                    }
                }
            }

            if (!shaderSource) {
                printf("Failed to load Metal shader source\n");
                return false;
            }

            MTLCompileOptions* options = [[MTLCompileOptions alloc] init];
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
            options.fastMathEnabled = YES;
#pragma clang diagnostic pop

            library = [impl->device newLibraryWithSource:shaderSource options:options error:&error];
            if (!library) {
                printf("Failed to compile Metal shaders: %s\n", [[error localizedDescription] UTF8String]);
                return false;
            }
        }
    }

    // Create compute pipeline for clear kernel
    id<MTLFunction> clearFunction = [library newFunctionWithName:@"clearKernel"];
    if (!clearFunction) {
        printf("Failed to find clearKernel function\n");
        return false;
    }
    impl->clearPipeline = [impl->device newComputePipelineStateWithFunction:clearFunction error:&error];
    if (!impl->clearPipeline) {
        printf("Failed to create clear pipeline: %s\n", [[error localizedDescription] UTF8String]);
        return false;
    }

    // Create compute pipeline for molecule rendering
    id<MTLFunction> renderFunction = [library newFunctionWithName:@"renderMoleculeKernel"];
    if (!renderFunction) {
        printf("Failed to find renderMoleculeKernel function\n");
        return false;
    }
    impl->renderMoleculePipeline = [impl->device newComputePipelineStateWithFunction:renderFunction error:&error];
    if (!impl->renderMoleculePipeline) {
        printf("Failed to create render pipeline: %s\n", [[error localizedDescription] UTF8String]);
        return false;
    }

    // Create compute pipeline for text rendering
    id<MTLFunction> textFunction = [library newFunctionWithName:@"renderTextKernel"];
    if (!textFunction) {
        printf("Failed to find renderTextKernel function\n");
        return false;
    }
    impl->renderTextPipeline = [impl->device newComputePipelineStateWithFunction:textFunction error:&error];
    if (!impl->renderTextPipeline) {
        printf("Failed to create text pipeline: %s\n", [[error localizedDescription] UTF8String]);
        return false;
    }

    printf("Metal compute pipelines created successfully\n");
    return true;
}

static bool createRenderTexture(MetalRendererImpl* impl, int width, int height) {
    MTLTextureDescriptor* textureDesc = [[MTLTextureDescriptor alloc] init];
    textureDesc.pixelFormat = MTLPixelFormatRGBA8Unorm;
    textureDesc.width = width;
    textureDesc.height = height;
    textureDesc.usage = MTLTextureUsageShaderRead | MTLTextureUsageShaderWrite;
    textureDesc.storageMode = MTLStorageModePrivate;

    // On macOS, use managed storage for better CPU access if needed
#if TARGET_OS_OSX
    textureDesc.storageMode = MTLStorageModeManaged;
#endif

    impl->renderTexture = [impl->device newTextureWithDescriptor:textureDesc];
    if (!impl->renderTexture) {
        printf("Failed to create render texture\n");
        return false;
    }

    impl->width = width;
    impl->height = height;

    printf("Render texture created: %dx%d\n", width, height);
    return true;
}

static bool createBuffers(MetalRendererImpl* impl) {
    // Atom buffer (MAX_ATOMS atoms)
    impl->atomBuffer = [impl->device newBufferWithLength:MAX_ATOMS * sizeof(Atom)
                                                  options:MTLResourceStorageModeShared];
    if (!impl->atomBuffer) {
        printf("Failed to create atom buffer\n");
        return false;
    }

    // Bond buffer (MAX_BONDS bonds)
    impl->bondBuffer = [impl->device newBufferWithLength:MAX_BONDS * sizeof(Bond)
                                                  options:MTLResourceStorageModeShared];
    if (!impl->bondBuffer) {
        printf("Failed to create bond buffer\n");
        return false;
    }

    // Render params buffer
    impl->paramsBuffer = [impl->device newBufferWithLength:sizeof(RenderParams)
                                                    options:MTLResourceStorageModeShared];
    if (!impl->paramsBuffer) {
        printf("Failed to create params buffer\n");
        return false;
    }

    // Text buffer (256 bytes max)
    impl->textBuffer = [impl->device newBufferWithLength:256
                                                  options:MTLResourceStorageModeShared];
    if (!impl->textBuffer) {
        printf("Failed to create text buffer\n");
        return false;
    }

    // Text params buffer
    impl->textParamsBuffer = [impl->device newBufferWithLength:sizeof(TextParams)
                                                        options:MTLResourceStorageModeShared];
    if (!impl->textParamsBuffer) {
        printf("Failed to create text params buffer\n");
        return false;
    }

    printf("Metal buffers created\n");
    return true;
}

// ============== MetalRenderer Class Implementation ==============

MetalRenderer::MetalRenderer() : pImpl(nullptr) {
}

MetalRenderer::~MetalRenderer() {
    cleanup();
}

bool MetalRenderer::init(void* windowHandle, int width, int height) {
    printf("Initializing Metal renderer (%dx%d)...\n", width, height);

    // Allocate implementation
    pImpl = new MetalRendererImpl();
    memset(pImpl, 0, sizeof(MetalRendererImpl));

    // Get default Metal device
    pImpl->device = MTLCreateSystemDefaultDevice();
    if (!pImpl->device) {
        printf("Failed to create Metal device\n");
        delete pImpl;
        pImpl = nullptr;
        return false;
    }
    printf("Metal device: %s\n", [[pImpl->device name] UTF8String]);

    // Create command queue
    pImpl->commandQueue = [pImpl->device newCommandQueue];
    if (!pImpl->commandQueue) {
        printf("Failed to create command queue\n");
        delete pImpl;
        pImpl = nullptr;
        return false;
    }

    // Determine optimal thread group size
    pImpl->threadGroupSize = MTLSizeMake(16, 16, 1);

    // Create compute pipelines
    if (!createComputePipelines(pImpl)) {
        delete pImpl;
        pImpl = nullptr;
        return false;
    }

    // Create render texture
    if (!createRenderTexture(pImpl, width, height)) {
        delete pImpl;
        pImpl = nullptr;
        return false;
    }

    // Create buffers
    if (!createBuffers(pImpl)) {
        delete pImpl;
        pImpl = nullptr;
        return false;
    }

    printf("Metal renderer initialized successfully\n");
    return true;
}

void MetalRenderer::cleanup() {
    if (!pImpl) return;

    // ARC will handle releasing Metal objects
    pImpl->renderTexture = nil;
    pImpl->atomBuffer = nil;
    pImpl->bondBuffer = nil;
    pImpl->paramsBuffer = nil;
    pImpl->textBuffer = nil;
    pImpl->textParamsBuffer = nil;
    pImpl->clearPipeline = nil;
    pImpl->renderMoleculePipeline = nil;
    pImpl->renderTextPipeline = nil;
    pImpl->commandQueue = nil;
    pImpl->device = nil;

    delete pImpl;
    pImpl = nullptr;

    printf("Metal renderer cleaned up\n");
}

bool MetalRenderer::resize(int width, int height) {
    if (!pImpl) return false;
    if (pImpl->width == width && pImpl->height == height) return true;

    // Release old texture
    pImpl->renderTexture = nil;

    // Create new texture
    return createRenderTexture(pImpl, width, height);
}

void MetalRenderer::render(const Molecule* mol, float rotX, float rotY, float zoom,
                           float offsetX, float offsetY) {
    if (!pImpl || !mol) return;

    // Copy molecule data to buffers
    memcpy([pImpl->atomBuffer contents], mol->atoms, mol->numAtoms * sizeof(Atom));
    memcpy([pImpl->bondBuffer contents], mol->bonds, mol->numBonds * sizeof(Bond));

    // Set render params
    RenderParams* params = (RenderParams*)[pImpl->paramsBuffer contents];
    params->width = pImpl->width;
    params->height = pImpl->height;
    params->numAtoms = mol->numAtoms;
    params->numBonds = mol->numBonds;
    params->rotX = rotX;
    params->rotY = rotY;
    params->zoom = zoom;
    params->offsetX = offsetX;
    params->offsetY = offsetY;

    // Create command buffer
    id<MTLCommandBuffer> commandBuffer = [pImpl->commandQueue commandBuffer];

    // Calculate grid size
    MTLSize gridSize = MTLSizeMake(pImpl->width, pImpl->height, 1);

    // Clear kernel
    {
        id<MTLComputeCommandEncoder> encoder = [commandBuffer computeCommandEncoder];
        [encoder setComputePipelineState:pImpl->clearPipeline];
        [encoder setTexture:pImpl->renderTexture atIndex:0];
        [encoder dispatchThreads:gridSize threadsPerThreadgroup:pImpl->threadGroupSize];
        [encoder endEncoding];
    }

    // Render molecule kernel
    {
        id<MTLComputeCommandEncoder> encoder = [commandBuffer computeCommandEncoder];
        [encoder setComputePipelineState:pImpl->renderMoleculePipeline];
        [encoder setTexture:pImpl->renderTexture atIndex:0];
        [encoder setBuffer:pImpl->atomBuffer offset:0 atIndex:0];
        [encoder setBuffer:pImpl->bondBuffer offset:0 atIndex:1];
        [encoder setBuffer:pImpl->paramsBuffer offset:0 atIndex:2];
        [encoder dispatchThreads:gridSize threadsPerThreadgroup:pImpl->threadGroupSize];
        [encoder endEncoding];
    }

    // Commit and wait
    [commandBuffer commit];
    [commandBuffer waitUntilCompleted];
}

void MetalRenderer::renderText(const char* text, int x, int y, int scale) {
    if (!pImpl || !text) return;

    int textLen = (int)strlen(text);
    if (textLen == 0 || textLen > 255) return;

    // Copy text to buffer
    memcpy([pImpl->textBuffer contents], text, textLen + 1);

    // Set text params
    TextParams* params = (TextParams*)[pImpl->textParamsBuffer contents];
    params->textLen = textLen;
    params->startX = x;
    params->startY = y;
    params->scale = scale;
    params->width = pImpl->width;
    params->height = pImpl->height;

    // Create command buffer
    id<MTLCommandBuffer> commandBuffer = [pImpl->commandQueue commandBuffer];

    // Calculate grid size
    MTLSize gridSize = MTLSizeMake(pImpl->width, pImpl->height, 1);

    // Render text kernel
    {
        id<MTLComputeCommandEncoder> encoder = [commandBuffer computeCommandEncoder];
        [encoder setComputePipelineState:pImpl->renderTextPipeline];
        [encoder setTexture:pImpl->renderTexture atIndex:0];
        [encoder setBuffer:pImpl->textBuffer offset:0 atIndex:0];
        [encoder setBuffer:pImpl->textParamsBuffer offset:0 atIndex:1];
        [encoder dispatchThreads:gridSize threadsPerThreadgroup:pImpl->threadGroupSize];
        [encoder endEncoding];
    }

    // Commit and wait
    [commandBuffer commit];
    [commandBuffer waitUntilCompleted];
}

void* MetalRenderer::getTexture() {
    if (!pImpl) return nullptr;
    return (__bridge void*)pImpl->renderTexture;
}

void MetalRenderer::getSize(int* width, int* height) {
    if (!pImpl) {
        if (width) *width = 0;
        if (height) *height = 0;
        return;
    }
    if (width) *width = pImpl->width;
    if (height) *height = pImpl->height;
}

#ifdef __OBJC__
id<MTLDevice> MetalRenderer::getDevice() {
    return pImpl ? pImpl->device : nil;
}

id<MTLTexture> MetalRenderer::getMTLTexture() {
    return pImpl ? pImpl->renderTexture : nil;
}
#endif

// Factory function
IRenderer* createMetalRenderer() {
    return new MetalRenderer();
}
