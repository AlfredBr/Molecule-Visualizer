/*
 * MolVis - macOS Entry Point
 *
 * A GPU-accelerated molecular visualization application.
 * Uses SDL2 + Metal + Dear ImGui on macOS.
 *
 * This is the macOS equivalent of main.cpp (Windows version).
 */

#import <Foundation/Foundation.h>
#import <Metal/Metal.h>
#import <QuartzCore/CAMetalLayer.h>

#include <SDL2/SDL.h>

#include "imgui.h"
#include "imgui_impl_sdl2.h"
#include "imgui_impl_metal.h"

#include "renderer/metal_renderer.h"
#include "molecule/molecule_db.h"

#include <stdio.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>

// ============== Configuration ==============

static const char* CONFIG_FILENAME = "molvis.ini";

struct WindowConfig {
    int x = SDL_WINDOWPOS_CENTERED;
    int y = SDL_WINDOWPOS_CENTERED;
    int width = 1280;
    int height = 800;
};

static char g_configPath[1024] = "";
static char g_imguiIniPath[1024] = "";

void InitConfigPath() {
    // Get user's home directory for config file
    const char* home = getenv("HOME");
    if (home) {
        snprintf(g_configPath, sizeof(g_configPath), "%s/.config/molvis.ini", home);
        snprintf(g_imguiIniPath, sizeof(g_imguiIniPath), "%s/.config/molvis_imgui.ini", home);

        // Create .config directory if it doesn't exist
        char configDir[1024];
        snprintf(configDir, sizeof(configDir), "%s/.config", home);
        mkdir(configDir, 0755);
    } else {
        strcpy(g_configPath, CONFIG_FILENAME);
        strcpy(g_imguiIniPath, "molvis_imgui.ini");
    }
}

void SaveWindowConfig(SDL_Window* window) {
    if (!window) return;

    int x, y, w, h;
    SDL_GetWindowPosition(window, &x, &y);
    SDL_GetWindowSize(window, &w, &h);

    FILE* f = fopen(g_configPath, "w");
    if (f) {
        fprintf(f, "[MainWindow]\n");
        fprintf(f, "Pos=%d,%d\n", x, y);
        fprintf(f, "Size=%d,%d\n", w, h);
        fflush(f);
        fclose(f);
        printf("Saved window config to %s\n", g_configPath);
    }
}

WindowConfig LoadWindowConfig() {
    WindowConfig cfg;
    FILE* f = fopen(g_configPath, "r");
    if (f) {
        char line[256];
        while (fgets(line, sizeof(line), f)) {
            sscanf(line, "Pos=%d,%d", &cfg.x, &cfg.y);
            sscanf(line, "Size=%d,%d", &cfg.width, &cfg.height);
        }
        fclose(f);
        printf("Loaded window config from %s\n", g_configPath);
    }
    return cfg;
}

// ============== Main Application ==============

int main(int argc, char* argv[]) {
    @autoreleasepool {
        printf("MolVis starting (macOS)...\n");

        // Initialize config path
        InitConfigPath();

        // Load window config
        WindowConfig winCfg = LoadWindowConfig();

        // Initialize SDL
        if (SDL_Init(SDL_INIT_VIDEO | SDL_INIT_TIMER) != 0) {
            printf("SDL_Init failed: %s\n", SDL_GetError());
            return 1;
        }

        // Create SDL window with Metal support
        SDL_Window* window = SDL_CreateWindow(
            "MolVis - Molecule Visualizer",
            winCfg.x, winCfg.y,
            winCfg.width, winCfg.height,
            SDL_WINDOW_METAL | SDL_WINDOW_RESIZABLE | SDL_WINDOW_ALLOW_HIGHDPI
        );
        if (!window) {
            printf("SDL_CreateWindow failed: %s\n", SDL_GetError());
            SDL_Quit();
            return 1;
        }
        printf("SDL window created\n");

        // Create Metal layer
        SDL_MetalView metalView = SDL_Metal_CreateView(window);
        CAMetalLayer* metalLayer = (__bridge CAMetalLayer*)SDL_Metal_GetLayer(metalView);

        // Get Metal device
        id<MTLDevice> device = MTLCreateSystemDefaultDevice();
        if (!device) {
            printf("Failed to create Metal device\n");
            SDL_DestroyWindow(window);
            SDL_Quit();
            return 1;
        }
        printf("Metal device: %s\n", [[device name] UTF8String]);

        metalLayer.device = device;
        metalLayer.pixelFormat = MTLPixelFormatBGRA8Unorm;

        // Create command queue for ImGui rendering
        id<MTLCommandQueue> commandQueue = [device newCommandQueue];

        // Initialize CUDA/Metal renderer
        printf("Initializing Metal renderer...\n");
        MetalRenderer* renderer = new MetalRenderer();
        if (!renderer->init((__bridge void*)metalLayer, 800, 600)) {
            printf("Failed to initialize Metal renderer\n");
            delete renderer;
            SDL_Metal_DestroyView(metalView);
            SDL_DestroyWindow(window);
            SDL_Quit();
            return 1;
        }
        printf("Metal renderer initialized\n");

        // Application state (matching Windows version)
        Molecule molecule = {};
        int currentMolecule = 0;
        float rotX = 0.3f;
        float rotY = 0.0f;
        float zoom = 10.0f;
        float offsetX = 0.0f;
        float offsetY = 0.0f;
        bool autoRotate = true;
        float rotSpeed = 1.0f;
        bool rotateX = false;
        bool rotateY = true;

        // Load initial molecule
        molecule_build(currentMolecule, &molecule);

        // Setup Dear ImGui context
        IMGUI_CHECKVERSION();
        ImGui::CreateContext();
        ImGuiIO& io = ImGui::GetIO();
        io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;

        // Set ImGui ini file path to persist window positions
        io.IniFilename = g_imguiIniPath;
        printf("ImGui config path: %s\n", g_imguiIniPath);

        // Setup Dear ImGui style
        ImGui::StyleColorsDark();

        // Load default font at larger size
        ImFontConfig fontConfig;
        fontConfig.SizePixels = 17.0f;
        io.Fonts->AddFontDefault(&fontConfig);

        // Setup Platform/Renderer backends
        ImGui_ImplSDL2_InitForMetal(window);
        ImGui_ImplMetal_Init(device);

        printf("ImGui initialized\n");

        // Main loop
        bool running = true;
        while (running) {
            @autoreleasepool {
                // Poll and handle events
                SDL_Event event;
                while (SDL_PollEvent(&event)) {
                    ImGui_ImplSDL2_ProcessEvent(&event);

                    if (event.type == SDL_QUIT) {
                        running = false;
                    }
                    if (event.type == SDL_WINDOWEVENT &&
                        event.window.event == SDL_WINDOWEVENT_CLOSE &&
                        event.window.windowID == SDL_GetWindowID(window)) {
                        running = false;
                    }
                }

                if (!running) break;

                // Get drawable size (handles Retina displays)
                int drawableWidth, drawableHeight;
                SDL_Metal_GetDrawableSize(window, &drawableWidth, &drawableHeight);
                metalLayer.drawableSize = CGSizeMake(drawableWidth, drawableHeight);

                // Get next drawable
                id<CAMetalDrawable> drawable = [metalLayer nextDrawable];
                if (!drawable) continue;

                // Create render pass descriptor for ImGui
                MTLRenderPassDescriptor* renderPassDesc = [MTLRenderPassDescriptor renderPassDescriptor];
                renderPassDesc.colorAttachments[0].texture = drawable.texture;
                renderPassDesc.colorAttachments[0].loadAction = MTLLoadActionClear;
                renderPassDesc.colorAttachments[0].storeAction = MTLStoreActionStore;
                renderPassDesc.colorAttachments[0].clearColor = MTLClearColorMake(0.1, 0.1, 0.15, 1.0);

                // Start the Dear ImGui frame
                ImGui_ImplMetal_NewFrame(renderPassDesc);
                ImGui_ImplSDL2_NewFrame();
                ImGui::NewFrame();

                // === MolVis UI ===

                // Main menu bar
                if (ImGui::BeginMainMenuBar()) {
                    if (ImGui::BeginMenu("File")) {
                        if (ImGui::MenuItem("Quit", "Cmd+Q")) {
                            running = false;
                        }
                        ImGui::EndMenu();
                    }
                    if (ImGui::BeginMenu("Molecule")) {
                        // Organize by category
                        for (int cat = 0; cat < CAT_COUNT; cat++) {
                            if (ImGui::BeginMenu(molecule_get_category_name(cat))) {
                                for (int i = 0; i < molecule_get_count(); i++) {
                                    if (molecule_get_category(i) == cat) {
                                        if (ImGui::MenuItem(molecule_get_name(i), nullptr, currentMolecule == i)) {
                                            currentMolecule = i;
                                            molecule_build(currentMolecule, &molecule);
                                        }
                                    }
                                }
                                ImGui::EndMenu();
                            }
                        }
                        ImGui::EndMenu();
                    }
                    if (ImGui::BeginMenu("View")) {
                        ImGui::MenuItem("Auto-Rotate", nullptr, &autoRotate);
                        if (ImGui::MenuItem("Reset View")) {
                            rotX = 0.3f;
                            rotY = 0.0f;
                            zoom = 10.0f;
                        }
                        ImGui::EndMenu();
                    }
                    ImGui::EndMainMenuBar();
                }

                // Auto-rotate
                if (autoRotate) {
                    if (rotateX) rotX += 0.01f * rotSpeed;
                    if (rotateY) rotY += 0.01f * rotSpeed;
                }

                // Resize renderer if needed
                int curW, curH;
                renderer->getSize(&curW, &curH);

                // Render molecule to Metal texture
                renderer->render(&molecule, rotX, rotY, zoom, offsetX, offsetY);

                // Render molecule name overlay
                renderer->renderText(molecule.name, 10, 10, 3);

                // Molecule Viewport Panel
                ImGui::SetNextWindowSize(ImVec2(820, 640), ImGuiCond_FirstUseEver);
                ImGui::Begin("Viewport", nullptr, ImGuiWindowFlags_NoScrollbar);
                {
                    ImVec2 viewportSize = ImGui::GetContentRegionAvail();
                    int vpWidth = (int)viewportSize.x;
                    int vpHeight = (int)viewportSize.y;

                    // Resize renderer if needed
                    if (vpWidth > 0 && vpHeight > 0 && (curW != vpWidth || curH != vpHeight)) {
                        renderer->resize(vpWidth, vpHeight);
                    }

                    // Display the Metal-rendered texture
                    id<MTLTexture> texture = (__bridge id<MTLTexture>)renderer->getTexture();
                    if (texture) {
                        ImGui::Image((ImTextureID)texture, viewportSize);
                    }

                    // Handle mouse interaction in viewport
                    if (ImGui::IsWindowHovered()) {
                        if (ImGui::IsMouseDragging(ImGuiMouseButton_Left)) {
                            ImVec2 delta = ImGui::GetMouseDragDelta(ImGuiMouseButton_Left);
                            ImGui::ResetMouseDragDelta(ImGuiMouseButton_Left);
                            rotY += delta.x * 0.01f;
                            rotX += delta.y * 0.01f;
                        }
                        float wheel = io.MouseWheel;
                        if (wheel != 0.0f) {
                            zoom -= wheel * 0.5f;
                            if (zoom < 3.0f) zoom = 3.0f;
                            if (zoom > 30.0f) zoom = 30.0f;
                        }
                    }
                }
                ImGui::End();

                // Molecule Info Panel
                ImGui::SetNextWindowSize(ImVec2(320, 350), ImGuiCond_FirstUseEver);
                ImGui::Begin("Molecule Info");
                {
                    ImGui::Text("Name: %s", molecule.name);
                    ImGui::TextWrapped("Description: %s", molecule_get_description(currentMolecule));
                    ImGui::Separator();

                    ImGui::Text("Atoms: %d", molecule.numAtoms);
                    ImGui::Text("Bonds: %d", molecule.numBonds);
                    ImGui::Text("Category: %s", molecule_get_category_name(molecule_get_category(currentMolecule)));
                    ImGui::Separator();

                    // Molecule browser by category
                    ImGui::Text("Browse Molecules:");
                    for (int cat = 0; cat < CAT_COUNT; cat++) {
                        if (ImGui::TreeNode(molecule_get_category_name(cat))) {
                            for (int i = 0; i < molecule_get_count(); i++) {
                                if (molecule_get_category(i) == cat) {
                                    bool isSelected = (currentMolecule == i);
                                    if (ImGui::Selectable(molecule_get_name(i), isSelected)) {
                                        currentMolecule = i;
                                        molecule_build(currentMolecule, &molecule);
                                    }
                                    if (ImGui::IsItemHovered()) {
                                        ImGui::SetTooltip("%s", molecule_get_description(i));
                                    }
                                }
                            }
                            ImGui::TreePop();
                        }
                    }
                }
                ImGui::End();

                // View Controls Panel
                ImGui::SetNextWindowSize(ImVec2(250, 220), ImGuiCond_FirstUseEver);
                ImGui::Begin("View Controls");
                {
                    ImGui::Checkbox("Auto-Rotate", &autoRotate);
                    if (autoRotate) {
                        ImGui::SliderFloat("Speed", &rotSpeed, -3.0f, 3.0f, "%.1fx");
                        ImGui::Checkbox("Rotate X", &rotateX);
                        ImGui::SameLine();
                        ImGui::Checkbox("Rotate Y", &rotateY);
                    }
                    ImGui::Separator();

                    ImGui::Text("Camera:");
                    ImGui::SliderFloat("Zoom", &zoom, 3.0f, 30.0f);
                    ImGui::SliderFloat("Rotation X", &rotX, -3.14159f, 3.14159f);
                    ImGui::SliderFloat("Rotation Y", &rotY, -3.14159f, 3.14159f);

                    ImGui::Separator();
                    ImGui::Text("Position:");
                    ImGui::SliderFloat("Horizontal", &offsetX, -2.0f, 2.0f);
                    ImGui::SliderFloat("Vertical", &offsetY, -2.0f, 2.0f);

                    ImGui::Separator();
                    if (ImGui::Button("Reset View")) {
                        rotX = 0.3f;
                        rotY = 0.0f;
                        zoom = 10.0f;
                        offsetX = 0.0f;
                        offsetY = 0.0f;
                    }
                    ImGui::SameLine();
                    if (ImGui::Button("Top View")) {
                        rotX = 1.57f;
                        rotY = 0.0f;
                    }
                    ImGui::SameLine();
                    if (ImGui::Button("Side View")) {
                        rotX = 0.0f;
                        rotY = 1.57f;
                    }
                }
                ImGui::End();

                // Performance Panel
                ImGui::SetNextWindowSize(ImVec2(200, 100), ImGuiCond_FirstUseEver);
                ImGui::Begin("Performance");
                {
                    ImGui::Text("%.3f ms/frame", 1000.0f / io.Framerate);
                    ImGui::Text("%.1f FPS", io.Framerate);
                }
                ImGui::End();

                // Rendering
                ImGui::Render();

                // Create command buffer for ImGui rendering
                id<MTLCommandBuffer> commandBuffer = [commandQueue commandBuffer];

                // Render ImGui
                id<MTLRenderCommandEncoder> renderEncoder = [commandBuffer renderCommandEncoderWithDescriptor:renderPassDesc];
                ImGui_ImplMetal_RenderDrawData(ImGui::GetDrawData(), commandBuffer, renderEncoder);
                [renderEncoder endEncoding];

                // Present
                [commandBuffer presentDrawable:drawable];
                [commandBuffer commit];
            }
        }

        // Save window config
        SaveWindowConfig(window);

        // Cleanup
        delete renderer;

        ImGui_ImplMetal_Shutdown();
        ImGui_ImplSDL2_Shutdown();
        ImGui::DestroyContext();

        SDL_Metal_DestroyView(metalView);
        SDL_DestroyWindow(window);
        SDL_Quit();

        printf("MolVis terminated\n");
        return 0;
    }
}
