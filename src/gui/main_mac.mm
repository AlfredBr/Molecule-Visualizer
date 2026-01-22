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
#include <cctype>

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
static char g_searchBuffer[256] = "";  // Search buffer for molecule browser
static int g_currentMolecule = 0;       // Current molecule index

// Forward declaration
static bool StringContains(const char* str, const char* search);

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
        fprintf(f, "[UserSettings]\n");
        fprintf(f, "LastMolecule=%d\n", g_currentMolecule);
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
            sscanf(line, "LastMolecule=%d", &g_currentMolecule);
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
        Molecule renderMolecule = {};
        // Note: g_currentMolecule is initialized from saved config via LoadWindowConfig()
        float rotX = 0.3f;
        float rotY = 0.0f;
        float zoom = 10.0f;
        float offsetX = 0.0f;
        float offsetY = 0.0f;
        bool autoRotate = true;
        float rotSpeed = 1.0f;
        bool rotateX = false;
        bool rotateY = true;
        static bool hideHydrogen = false;

        // Helper: Build filtered molecule (optionally hide hydrogens)
        auto BuildRenderMolecule = [](const Molecule& src, Molecule& dst, bool hideH) {
            dst.numAtoms = 0; dst.numBonds = 0;
            strncpy(dst.name, src.name, sizeof(dst.name));
            int mapOldToNew[MAX_ATOMS];
            for (int i = 0; i < src.numAtoms; ++i) mapOldToNew[i] = -1;
            for (int i = 0; i < src.numAtoms; ++i) {
                const Atom& a = src.atoms[i];
                if (hideH && a.type == ATOM_H) { mapOldToNew[i] = -1; continue; }
                if (dst.numAtoms >= MAX_ATOMS) break;
                mapOldToNew[i] = dst.numAtoms;
                dst.atoms[dst.numAtoms] = a;
                dst.numAtoms++;
            }
            for (int i = 0; i < src.numBonds; ++i) {
                int n1 = mapOldToNew[src.bonds[i].atom1];
                int n2 = mapOldToNew[src.bonds[i].atom2];
                if (n1 < 0 || n2 < 0) continue;
                if (dst.numBonds >= MAX_BONDS) break;
                dst.bonds[dst.numBonds].atom1 = n1;
                dst.bonds[dst.numBonds].atom2 = n2;
                dst.bonds[dst.numBonds].order = src.bonds[i].order;
                dst.numBonds++;
            }
        };

        // Load initial molecule
        molecule_build(g_currentMolecule, &molecule);
        BuildRenderMolecule(molecule, renderMolecule, hideHydrogen);

        // Setup Dear ImGui context
        IMGUI_CHECKVERSION();
        ImGui::CreateContext();
        ImGuiIO& io = ImGui::GetIO();
        io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;
        io.ConfigFlags |= ImGuiConfigFlags_DockingEnable;  // Uncomment if using docking branch

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

                // Enable Docking
		        ImGui::DockSpaceOverViewport();

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
                                        if (ImGui::MenuItem(molecule_get_name(i), nullptr, g_currentMolecule == i)) {
                                            g_currentMolecule = i;
                                            molecule_build(g_currentMolecule, &molecule);
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

                // Render filtered molecule to Metal texture
                renderer->render(&renderMolecule, rotX, rotY, zoom, offsetX, offsetY);

                // Render molecule name overlay
                renderer->renderText(renderMolecule.name, 10, 10, 3);

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
                    if (ImGui::IsItemHovered()) {
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

                // Browse Molecules Panel
                ImGui::SetNextWindowSize(ImVec2(280, 400), ImGuiCond_FirstUseEver);
                ImGui::Begin("Browse Molecules");
                {
                    // Search input field with clear button
                    ImGui::InputTextWithHint("##search", "Search molecules...", g_searchBuffer, sizeof(g_searchBuffer));
                    ImGui::SameLine();
                    if (ImGui::Button("Clear", ImVec2(50, 0)))
                    {
                        g_searchBuffer[0] = '\0';
                    }
                    ImGui::Separator();

                    // Browse molecules with search filtering
                    bool hasSearchResults = false;

                    if (g_searchBuffer[0] == '\0')
                    {
                        // No search: show categories
                        for (int cat = 0; cat < CAT_COUNT; cat++)
                        {
                            if (ImGui::TreeNode(molecule_get_category_name(cat)))
                            {
                                for (int i = 0; i < molecule_get_count(); i++)
                                {
                                    if (molecule_get_category(i) == cat)
                                    {
                                        bool isSelected = (g_currentMolecule == i);
                                        if (ImGui::Selectable(molecule_get_name(i), isSelected))
                                        {
                                            g_currentMolecule = i;
                                            molecule_build(g_currentMolecule, &molecule);
                                        }
                                        if (ImGui::IsItemHovered())
                                        {
                                            ImGui::SetTooltip("%s", molecule_get_description(i));
                                        }
                                    }
                                }
                                ImGui::TreePop();
                            }
                        }
                    }
                    else
                    {
                        // Search: show matching molecules across all categories
                        for (int i = 0; i < molecule_get_count(); i++)
                        {
                            const char* name = molecule_get_name(i);
                            const char* desc = molecule_get_description(i);

                            if (StringContains(name, g_searchBuffer) || StringContains(desc, g_searchBuffer))
                            {
                                hasSearchResults = true;
                                bool isSelected = (g_currentMolecule == i);
                                if (ImGui::Selectable(name, isSelected))
                                {
                                    g_currentMolecule = i;
                                    molecule_build(g_currentMolecule, &molecule);
                                }
                                if (ImGui::IsItemHovered())
                                {
                                    ImGui::SetTooltip("%s", desc);
                                }
                            }
                        }

                        if (!hasSearchResults)
                        {
                            ImGui::TextDisabled("No matches found");
                        }
                    }
                }
                ImGui::End();


                // Molecule Options Panel
                ImGui::SetNextWindowSize(ImVec2(240, 80), ImGuiCond_FirstUseEver);
                ImGui::Begin("Molecule Options");
                {
                    ImGui::TextDisabled("Rendering options");
                    ImGui::Separator();
                    if (ImGui::Checkbox("Hide Hydrogen (H)", &hideHydrogen)) {
                        BuildRenderMolecule(molecule, renderMolecule, hideHydrogen);
                    }
                    ImGui::TextDisabled("Does not affect Periodic Table");
                }
                ImGui::End();

                // Periodic Table Panel
                ImGui::SetNextWindowSize(ImVec2(860, 360), ImGuiCond_FirstUseEver);
                ImGui::Begin("Periodic Table");
                {
                    struct ElementInfo { const char* symbol; const char* name; int period; int group; int atomType; };
                    static const ElementInfo elements[] = {
                        {"H",  "Hydrogen",        1,  1, ATOM_H},     {"He", "Helium",          1, 18, -1},
                        {"Li", "Lithium",         2,  1, -1},         {"Be", "Beryllium",       2,  2, -1},
                        {"B",  "Boron",           2, 13, ATOM_B},     {"C",  "Carbon",          2, 14, ATOM_C},
                        {"N",  "Nitrogen",        2, 15, ATOM_N},     {"O",  "Oxygen",          2, 16, ATOM_O},
                        {"F",  "Fluorine",        2, 17, ATOM_F},     {"Ne", "Neon",            2, 18, -1},
                        {"Na", "Sodium",          3,  1, ATOM_NA},    {"Mg", "Magnesium",       3,  2, -1},
                        {"Al", "Aluminum",        3, 13, ATOM_AL},    {"Si", "Silicon",         3, 14, ATOM_SI},
                        {"P",  "Phosphorus",      3, 15, ATOM_P},     {"S",  "Sulfur",          3, 16, ATOM_S},
                        {"Cl", "Chlorine",        3, 17, ATOM_CL},    {"Ar", "Argon",           3, 18, -1},
                        {"K",  "Potassium",       4,  1, -1},         {"Ca", "Calcium",         4,  2, -1},
                        {"Sc", "Scandium",        4,  3, -1},         {"Ti", "Titanium",        4,  4, ATOM_TI},
                        {"V",  "Vanadium",        4,  5, -1},         {"Cr", "Chromium",        4,  6, -1},
                        {"Mn", "Manganese",       4,  7, -1},         {"Fe", "Iron",            4,  8, ATOM_FE},
                        {"Co", "Cobalt",          4,  9, -1},         {"Ni", "Nickel",          4, 10, -1},
                        {"Cu", "Copper",          4, 11, ATOM_CU},    {"Zn", "Zinc",            4, 12, -1},
                        {"Ga", "Gallium",         4, 13, -1},         {"Ge", "Germanium",       4, 14, -1},
                        {"As", "Arsenic",         4, 15, -1},         {"Se", "Selenium",        4, 16, -1},
                        {"Br", "Bromine",         4, 17, ATOM_BR},    {"Kr", "Krypton",         4, 18, -1},
                        {"Rb", "Rubidium",        5,  1, -1},         {"Sr", "Strontium",       5,  2, -1},
                        {"Y",  "Yttrium",         5,  3, -1},         {"Zr", "Zirconium",       5,  4, -1},
                        {"Nb", "Niobium",         5,  5, -1},         {"Mo", "Molybdenum",      5,  6, -1},
                        {"Tc", "Technetium",      5,  7, -1},         {"Ru", "Ruthenium",       5,  8, -1},
                        {"Rh", "Rhodium",         5,  9, -1},         {"Pd", "Palladium",       5, 10, -1},
                        {"Ag", "Silver",          5, 11, -1},         {"Cd", "Cadmium",         5, 12, -1},
                        {"In", "Indium",          5, 13, -1},         {"Sn", "Tin",             5, 14, -1},
                        {"Sb", "Antimony",        5, 15, -1},         {"Te", "Tellurium",       5, 16, -1},
                        {"I",  "Iodine",          5, 17, ATOM_I},     {"Xe", "Xenon",           5, 18, -1},
                        {"Cs", "Cesium",          6,  1, -1},         {"Ba", "Barium",          6,  2, -1},
                        {"La", "Lanthanum",       6,  3, -1},         {"Ce", "Cerium",          6,  4, -1},
                        {"Pr", "Praseodymium",    6,  5, -1},         {"Nd", "Neodymium",       6,  6, -1},
                        {"Pm", "Promethium",      6,  7, -1},         {"Sm", "Samarium",        6,  8, -1},
                        {"Eu", "Europium",        6,  9, -1},         {"Gd", "Gadolinium",      6, 10, -1},
                        {"Tb", "Terbium",         6, 11, -1},         {"Dy", "Dysprosium",      6, 12, -1},
                        {"Ho", "Holmium",         6, 13, -1},         {"Er", "Erbium",          6, 14, -1},
                        {"Tm", "Thulium",         6, 15, -1},         {"Yb", "Ytterbium",      6, 16, -1},
                        {"Lu", "Lutetium",        6, 17, -1},         {"Hf", "Hafnium",         6,  4, -1},
                        {"Ta", "Tantalum",        6,  5, -1},         {"W",  "Tungsten",        6,  6, -1},
                        {"Re", "Rhenium",         6,  7, ATOM_RE},    {"Os", "Osmium",          6,  8, -1},
                        {"Ir", "Iridium",         6,  9, -1},         {"Pt", "Platinum",        6, 10, ATOM_PT},
                        {"Au", "Gold",            6, 11, -1},         {"Hg", "Mercury",         6, 12, -1},
                        {"Tl", "Thallium",        6, 13, -1},         {"Pb", "Lead",            6, 14, -1},
                        {"Bi", "Bismuth",         6, 15, -1},         {"Po", "Polonium",        6, 16, -1},
                        {"At", "Astatine",        6, 17, -1},         {"Rn", "Radon",           6, 18, -1},
                        {"Fr", "Francium",        7,  1, -1},         {"Ra", "Radium",          7,  2, -1},
                        {"Ac", "Actinium",        7,  3, -1},         {"Th", "Thorium",         7,  4, -1},
                        {"Pa", "Protactinium",    7,  5, -1},         {"U",  "Uranium",         7,  6, -1},
                        {"Np", "Neptunium",       7,  7, -1},         {"Pu", "Plutonium",       7,  8, -1},
                        {"Am", "Americium",       7,  9, -1},         {"Cm", "Curium",          7, 10, -1},
                        {"Bk", "Berkelium",       7, 11, -1},         {"Cf", "Californium",     7, 12, -1},
                        {"Es", "Einsteinium",     7, 13, -1},         {"Fm", "Fermium",         7, 14, -1},
                        {"Md", "Mendelevium",     7, 15, -1},         {"No", "Nobelium",        7, 16, -1},
                        {"Lr", "Lawrencium",      7, 17, -1},         {"Rf", "Rutherfordium",   7,  4, -1},
                        {"Db", "Dubnium",         7,  5, -1},         {"Sg", "Seaborgium",      7,  6, -1},
                        {"Bh", "Bohrium",         7,  7, -1},         {"Hs", "Hassium",         7,  8, -1},
                        {"Mt", "Meitnerium",      7,  9, -1},         {"Ds", "Darmstadtium",    7, 10, -1},
                        {"Rg", "Roentgenium",     7, 11, -1},         {"Cn", "Copernicium",     7, 12, -1},
                        {"Nh", "Nihonium",        7, 13, -1},         {"Fl", "Flerovium",       7, 14, -1},
                        {"Mc", "Moscovium",       7, 15, -1},         {"Lv", "Livermorium",     7, 16, -1},
                        {"Ts", "Tennessine",      7, 17, -1},         {"Og", "Oganesson",       7, 18, -1},
                    };

                    // CPK color lookup for supported internal atom types
                    auto get_cpk_color = [](int atomType) -> ImVec4 {
                        switch (atomType) {
                        case ATOM_H:  return ImVec4(0.95f, 0.95f, 0.95f, 1.0f);
                        case ATOM_C:  return ImVec4(0.20f, 0.20f, 0.20f, 1.0f);
                        case ATOM_N:  return ImVec4(0.20f, 0.30f, 0.90f, 1.0f);
                        case ATOM_O:  return ImVec4(0.90f, 0.20f, 0.20f, 1.0f);
                        case ATOM_P:  return ImVec4(1.00f, 0.50f, 0.00f, 1.0f);
                        case ATOM_S:  return ImVec4(0.90f, 0.80f, 0.20f, 1.0f);
                        case ATOM_CL: return ImVec4(0.20f, 0.90f, 0.20f, 1.0f);
                        case ATOM_BR: return ImVec4(0.60f, 0.10f, 0.10f, 1.0f);
                        case ATOM_F:  return ImVec4(0.50f, 0.90f, 0.50f, 1.0f);
                        case ATOM_I:  return ImVec4(0.50f, 0.10f, 0.50f, 1.0f);
                        case ATOM_NA: return ImVec4(0.70f, 0.50f, 0.90f, 1.0f);
                        case ATOM_SI: return ImVec4(0.85f, 0.75f, 0.55f, 1.0f);
                        case ATOM_B:  return ImVec4(1.00f, 0.65f, 0.65f, 1.0f);
                        case ATOM_FE: return ImVec4(0.88f, 0.40f, 0.20f, 1.0f);
                        case ATOM_CU: return ImVec4(0.85f, 0.55f, 0.20f, 1.0f);
                        case ATOM_AL: return ImVec4(0.75f, 0.75f, 0.80f, 1.0f);
                        case ATOM_TI: return ImVec4(0.60f, 0.60f, 0.65f, 1.0f);
                        case ATOM_PT: return ImVec4(0.85f, 0.85f, 0.88f, 1.0f);
                        case ATOM_RE: return ImVec4(0.51f, 0.51f, 0.56f, 1.0f);
                        default:      return ImVec4(0.32f, 0.32f, 0.36f, 1.0f);
                        }
                    };

                    // Compute presence of each supported type in current molecule (not affected by hideHydrogen)
                    bool present[ATOM_TYPE_COUNT] = {false};
                    for (int i = 0; i < molecule.numAtoms; ++i) {
                        int t = molecule.atoms[i].type;
                        if (t >= 0 && t < ATOM_TYPE_COUNT) present[t] = true;
                    }

                    // Layout constants
                    const ImVec2 cellSize(36, 30);
                    const ImVec2 spacing(6, 6);
                    const float topMargin = 18.0f;
                    const float leftMargin = 18.0f;
                    ImVec2 origin = ImGui::GetCursorScreenPos();

                    // Group number row
                    for (int c = 1; c <= 18; ++c) {
                        ImVec2 pos(origin.x + leftMargin + (c-1)*(cellSize.x + spacing.x) + 10.0f, origin.y);
                        ImGui::SetCursorScreenPos(pos);
                        ImGui::Text("%d", c);
                    }

                    // Draw each element as a colored button placed by (period, group)
                    for (const auto& e : elements) {
                        float x = leftMargin + (e.group - 1) * (cellSize.x + spacing.x);
                        float y = topMargin + (e.period - 1) * (cellSize.y + spacing.y);
                        ImGui::SetCursorScreenPos(ImVec2(origin.x + x, origin.y + y));
                        bool isPresent = (e.atomType >= 0 && e.atomType < ATOM_TYPE_COUNT) ? present[e.atomType] : false;
                        ImVec4 neutral = ImVec4(0.32f, 0.32f, 0.36f, 1.0f);
                        ImVec4 base = isPresent ? get_cpk_color(e.atomType) : neutral;
                        ImVec4 hovered = ImVec4(base.x * 0.9f, base.y * 0.9f, base.z * 0.9f, 1.0f);
                        ImVec4 active  = ImVec4(base.x * 0.8f, base.y * 0.8f, base.z * 0.8f, 1.0f);
                        ImGui::PushID(e.symbol);
                        ImGui::PushStyleVar(ImGuiStyleVar_FrameRounding, 3.0f);
                        ImGui::PushStyleVar(ImGuiStyleVar_FrameBorderSize, isPresent ? 2.0f : 0.0f);
                        ImGui::PushStyleColor(ImGuiCol_Border, isPresent ? ImVec4(1.0f,1.0f,0.2f,1.0f) : ImVec4(0,0,0,0));
                        ImGui::PushStyleColor(ImGuiCol_Button, base);
                        ImGui::PushStyleColor(ImGuiCol_ButtonHovered, hovered);
                        ImGui::PushStyleColor(ImGuiCol_ButtonActive, active);
                        ImVec4 textCol = (base.x*0.299f + base.y*0.587f + base.z*0.114f) < 0.5f ? ImVec4(1,1,1,1) : ImVec4(0,0,0,1);
                        ImGui::PushStyleColor(ImGuiCol_Text, textCol);
                        ImGui::Button(e.symbol, cellSize);
                        if (ImGui::IsItemHovered()) {
                            ImGui::BeginTooltip();
                            ImGui::PushStyleColor(ImGuiCol_Text, ImVec4(1,1,1,1));
                            ImGui::TextUnformatted(e.name);
                            if (isPresent) {
                                ImGui::Spacing();
                                ImGui::TextColored(ImVec4(1.0f, 0.95f, 0.3f, 1.0f), "Present in molecule");
                            }
                            ImGui::PopStyleColor();
                            ImGui::EndTooltip();
                        }
                        ImGui::PopStyleColor(5);
                        ImGui::PopStyleVar(2);
                        ImGui::PopID();
                    }
                    ImGui::Dummy(ImVec2(0, topMargin + 7*(cellSize.y + spacing.y)));
                    ImGui::Separator();
                    ImGui::TextDisabled("Highlighted elements are present in the selected molecule.");
                }
                ImGui::End();

                // Description Panel - Educational info about the current molecule
                ImGui::SetNextWindowSize(ImVec2(350, 280), ImGuiCond_FirstUseEver);
                ImGui::Begin("Description");
                {
                    ImGui::PushStyleColor(ImGuiCol_Text, ImVec4(1.0f, 0.9f, 0.7f, 1.0f)); // Warm yellow text
                    ImGui::TextWrapped("%s", renderMolecule.name);
                    ImGui::PopStyleColor();
                    ImGui::Separator();
                    ImGui::Spacing();
                    ImGui::TextWrapped("%s", molecule_get_long_description(g_currentMolecule));
                    ImGui::Spacing();
                    ImGui::Separator();
                    ImGui::Text("Atoms: %d", renderMolecule.numAtoms);
                    ImGui::Text("Bonds: %d", renderMolecule.numBonds);
                    ImGui::Text("Category: %s", molecule_get_category_name(molecule_get_category(g_currentMolecule)));
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

                // CPK Color Reference Panel
                ImGui::SetNextWindowSize(ImVec2(280, 450), ImGuiCond_FirstUseEver);
                ImGui::Begin("CPK Color Reference");
                {
                    ImGui::Text("CPK Atom Colors:");
                    ImGui::Separator();

                    // Define colors matching the CUDA/Windows renderer
                    ImVec4 cpkColors[] = {
                        {0.95f, 0.95f, 0.95f, 1.0f},  // H - white
                        {0.2f,  0.2f,  0.2f,  1.0f},  // C - dark gray
                        {0.2f,  0.3f,  0.9f,  1.0f},  // N - blue
                        {0.9f,  0.2f,  0.2f,  1.0f},  // O - red
                        {1.0f,  0.5f,  0.0f,  1.0f},  // P - orange
                        {0.9f,  0.8f,  0.2f,  1.0f},  // S - yellow
                        {0.2f,  0.9f,  0.2f,  1.0f},  // Cl - green
                        {0.6f,  0.1f,  0.1f,  1.0f},  // Br - dark red
                        {0.5f,  0.9f,  0.5f,  1.0f},  // F - light green
                        {0.5f,  0.1f,  0.5f,  1.0f},  // I - purple
                        {0.7f,  0.5f,  0.9f,  1.0f},  // Na - metallic purple
                        {0.85f, 0.75f, 0.55f, 1.0f},  // Si - tan
                        {1.0f,  0.65f, 0.65f, 1.0f},  // B - salmon
                        {0.88f, 0.4f,  0.2f,  1.0f},  // Fe - orange/brown
                        {0.85f, 0.55f, 0.2f,  1.0f},  // Cu - copper
                        {0.75f, 0.75f, 0.8f,  1.0f},  // Al - silver
                        {0.6f,  0.6f,  0.65f, 1.0f},  // Ti - gray
                        {0.85f, 0.85f, 0.88f, 1.0f},  // Pt - white
                        {0.51f, 0.51f, 0.56f, 1.0f},  // Re - metallic gray
                    };

                    const char* atomNames[] = {
                        "H - Hydrogen (White)",
                        "C - Carbon (Dark Gray)",
                        "N - Nitrogen (Blue)",
                        "O - Oxygen (Red)",
                        "P - Phosphorus (Orange)",
                        "S - Sulfur (Yellow)",
                        "Cl - Chlorine (Green)",
                        "Br - Bromine (Dark Red)",
                        "F - Fluorine (Light Green)",
                        "I - Iodine (Purple)",
                        "Na - Sodium (Metallic Purple)",
                        "Si - Silicon (Tan)",
                        "B - Boron (Salmon)",
                        "Fe - Iron (Orange/Brown)",
                        "Cu - Copper (Copper)",
                        "Al - Aluminum (Silver)",
                        "Ti - Titanium (Gray)",
                        "Pt - Platinum (White)",
                        "Re - Rhenium (Gray)",
                    };

                    // Display color swatches in a grid
                    for (int i = 0; i < ATOM_TYPE_COUNT; i++) {
                        char id[32];
                        snprintf(id, sizeof(id), "##cpk_color_%d", i);
                        ImGui::ColorButton(id, cpkColors[i], ImGuiColorEditFlags_NoBorder | ImGuiColorEditFlags_NoTooltip, ImVec2(20, 20));
                        ImGui::SameLine();
                        ImGui::Text("%s", atomNames[i]);
                    }
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

// === Helper Functions ===

// Case-insensitive string search
static bool StringContains(const char* str, const char* search)
{
    if (!str || !search || search[0] == '\0')
        return true;  // Empty search matches everything

    // Simple case-insensitive substring search
    for (size_t i = 0; str[i]; ++i)
    {
        bool match = true;
        for (size_t j = 0; search[j]; ++j)
        {
            char c1 = tolower((unsigned char)str[i + j]);
            char c2 = tolower((unsigned char)search[j]);
            if (c1 != c2)
            {
                match = false;
                break;
            }
        }
        if (match)
            return true;
    }
    return false;
}
