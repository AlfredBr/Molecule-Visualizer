/*
 * MolVis - CUDA Molecular Visualization with Dear ImGui
 *
 * A GPU-accelerated molecular visualization application.
 * Uses Win32 + DirectX 11 + Dear ImGui + CUDA interop.
 */

#include "imgui.h"
#include "imgui_impl_win32.h"
#include "imgui_impl_dx11.h"
#include <d3d11.h>
#include <tchar.h>
#include <stdio.h>
#include <cctype>
#include <string.h>

#include "renderer/cuda_renderer.h"
#include "molecule/molecule_db.h"
#include "include/windowicon.h"

// App config file (stores main window position, last molecule, etc.)
static char g_configPath[MAX_PATH] = "molvis.ini";

// Global window handle for config save on close
static HWND g_hwnd = nullptr;

// Current molecule index
static int g_currentMolecule = 0;

struct WindowConfig {
    int x = 100;
    int y = 100;
    int width = 1280;
    int height = 800;
};

void InitConfigPath() {
    // Get executable directory for config file
    char exePath[MAX_PATH];
    GetModuleFileNameA(NULL, exePath, MAX_PATH);
    char* lastSlash = strrchr(exePath, '\\');
    if (lastSlash) {
        *(lastSlash + 1) = '\0';
        snprintf(g_configPath, MAX_PATH, "%smolvis.ini", exePath);
    }
}

void SaveWindowConfig() {
    if (!g_hwnd) return;
    RECT rect;
    if (GetWindowRect(g_hwnd, &rect)) {
        FILE* f = fopen(g_configPath, "w");
        if (f) {
            fprintf(f, "[MainWindow]\n");
            fprintf(f, "Pos=%d,%d\n", rect.left, rect.top);
            fprintf(f, "Size=%d,%d\n", rect.right - rect.left, rect.bottom - rect.top);
            fprintf(f, "[UserSettings]\n");
			fprintf(f, "LastMolecule=%d\n", g_currentMolecule);
            fflush(f);
            fclose(f);
        }
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
    }
    return cfg;
}

// DirectX 11 globals
static ID3D11Device*            g_pd3dDevice = nullptr;
static ID3D11DeviceContext*     g_pd3dDeviceContext = nullptr;
static IDXGISwapChain*          g_pSwapChain = nullptr;
static ID3D11RenderTargetView*  g_mainRenderTargetView = nullptr;

// CUDA Renderer
static CudaRenderer*            g_pRenderer = nullptr;
static Molecule                 g_molecule = {};
static float                    g_rotX = 0.3f;
static float                    g_rotY = 0.0f;
static float                    g_zoom = 10.0f;
static float                    g_offsetX = 0.0f;
static float                    g_offsetY = 0.0f;
static bool                     g_autoRotate = true;
static float                    g_rotSpeed = 1.0f;
static bool                     g_rotateX = false;
static bool                     g_rotateY = true;

// Molecule Browser Search
static char                     g_searchBuffer[256] = "";
// UI fonts
static ImFont*                  g_fontTooltip = nullptr;
// Molecule options
static bool                     g_hideHydrogen = false; // affects rendering only
// Docking options
static bool                     g_lockDocking = true;   // prevent undocking windows
static bool                     g_lockResize = true;    // prevent resizing docked windows

// Forward declarations
bool StringContains(const char* str, const char* search);
bool CreateDeviceD3D(HWND hWnd);
void CleanupDeviceD3D();
void CreateRenderTarget();
void CleanupRenderTarget();
LRESULT WINAPI WndProc(HWND hWnd, UINT msg, WPARAM wParam, LPARAM lParam);

// Build a copy of molecule used for rendering, optionally hiding hydrogens
static void BuildRenderMolecule(const Molecule& src, Molecule& dst, bool hideHydrogen)
{
    dst.numAtoms = 0; dst.numBonds = 0;
    strncpy(dst.name, src.name, sizeof(dst.name));
    strncpy(dst.formula, src.formula, sizeof(dst.formula));
    int mapOldToNew[MAX_ATOMS];
    for (int i = 0; i < src.numAtoms; ++i) mapOldToNew[i] = -1;
    // Copy atoms
    for (int i = 0; i < src.numAtoms; ++i) {
        const Atom& a = src.atoms[i];
        if (hideHydrogen && a.type == ATOM_H) { mapOldToNew[i] = -1; continue; }
        if (dst.numAtoms >= MAX_ATOMS) break;
        mapOldToNew[i] = dst.numAtoms;
        dst.atoms[dst.numAtoms] = a;
        dst.numAtoms++;
    }
    // Copy bonds, skipping any that reference hidden atoms; remap indices
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
}

// Compute element presence flags for our internal 18 types.
// Uses actual atoms plus a best-effort formula parse from the display name to include elements like H
// even if the model omits them. Rendering options (e.g., hide H) do not affect this.
static void ComputePresenceFlags(const Molecule& mol, int moleculeIndex, bool present[ATOM_TYPE_COUNT])
{
    for (int i = 0; i < ATOM_TYPE_COUNT; ++i) present[i] = false;
    // From atoms in the model
    for (int i = 0; i < mol.numAtoms; ++i) {
        int t = mol.atoms[i].type; if (t >= 0 && t < ATOM_TYPE_COUNT) present[t] = true;
    }
    // From formula in name (if available)
    auto markBySymbol = [&](const char* sym){
        if (strcmp(sym, "H") == 0) present[ATOM_H] = true;
        else if (strcmp(sym, "C") == 0) present[ATOM_C] = true;
        else if (strcmp(sym, "N") == 0) present[ATOM_N] = true;
        else if (strcmp(sym, "O") == 0) present[ATOM_O] = true;
        else if (strcmp(sym, "P") == 0) present[ATOM_P] = true;
        else if (strcmp(sym, "S") == 0) present[ATOM_S] = true;
        else if (strcmp(sym, "F") == 0) present[ATOM_F] = true;
        else if (strcmp(sym, "Cl") == 0) present[ATOM_CL] = true;
        else if (strcmp(sym, "Br") == 0) present[ATOM_BR] = true;
        else if (strcmp(sym, "I") == 0) present[ATOM_I] = true;
        else if (strcmp(sym, "Na") == 0) present[ATOM_NA] = true;
        else if (strcmp(sym, "Si") == 0) present[ATOM_SI] = true;
        else if (strcmp(sym, "B") == 0) present[ATOM_B] = true;
        else if (strcmp(sym, "Fe") == 0) present[ATOM_FE] = true;
        else if (strcmp(sym, "Cu") == 0) present[ATOM_CU] = true;
        else if (strcmp(sym, "Al") == 0) present[ATOM_AL] = true;
        else if (strcmp(sym, "Ti") == 0) present[ATOM_TI] = true;
        else if (strcmp(sym, "Pt") == 0) present[ATOM_PT] = true;
        else if (strcmp(sym, "Re") == 0) present[ATOM_RE] = true;
        else if (strcmp(sym, "Xe") == 0) present[ATOM_XE] = true;
    };

    const char* formula = molecule_get_formula(moleculeIndex);
    if (formula && formula[0]) {
        // Parse formula string directly (e.g., "H2O", "C6H12O6", "NaCl")
        const char* end = formula + strlen(formula);
        for (const char* p = formula; p < end; ) {
            if (*p >= 'A' && *p <= 'Z') {
                char sym[3] = {0,0,0}; sym[0] = *p; ++p;
                if (p < end && *p >= 'a' && *p <= 'z') { sym[1] = *p; ++p; }
                markBySymbol(sym);
                // skip optional digits
                while (p < end && *p >= '0' && *p <= '9') ++p;
            } else {
                ++p; // skip parentheses, charges, etc.
            }
        }
    }
}

int WINAPI WinMain(HINSTANCE hInstance, HINSTANCE hPrevInstance, LPSTR lpCmdLine, int nCmdShow)
{
    // Allocate console for debug output
    AllocConsole();
    freopen("CONOUT$", "w", stdout);
    freopen("CONOUT$", "w", stderr);
    printf("MolVis starting...\n");

    // Initialize config path to executable directory
    InitConfigPath();

    // Load window config from molvis.ini
    WindowConfig winCfg = LoadWindowConfig();

    // Create application window
    WNDCLASSEXW wc = {
        sizeof(wc), CS_CLASSDC, WndProc, 0L, 0L,
        GetModuleHandle(nullptr), nullptr, nullptr, nullptr, nullptr,
        L"MolVis", nullptr
    };
    RegisterClassExW(&wc);

    g_hwnd = CreateWindowW(
        wc.lpszClassName, L"MolVis - Molecule Visualizer",
        WS_OVERLAPPEDWINDOW,
        winCfg.x, winCfg.y, winCfg.width, winCfg.height,
        nullptr, nullptr, wc.hInstance, nullptr
    );
    HWND hwnd = g_hwnd;  // Local alias for convenience

    // Set window icon for title bar and taskbar
    SetWindowIconFromResource(hwnd);

    // Initialize Direct3D
    printf("Creating D3D device...\n");
    if (!CreateDeviceD3D(hwnd))
    {
        printf("Failed to create D3D device\n");
        CleanupDeviceD3D();
        UnregisterClassW(wc.lpszClassName, wc.hInstance);
        return 1;
    }
    printf("D3D device created\n");

    ShowWindow(hwnd, SW_SHOWDEFAULT);
    UpdateWindow(hwnd);

    // Initialize CUDA renderer
    printf("Initializing CUDA renderer...\n");
    g_pRenderer = renderer_init(g_pd3dDevice, 800, 600);
    if (!g_pRenderer)
    {
        printf("Failed to initialize CUDA renderer\n");
        CleanupDeviceD3D();
        UnregisterClassW(wc.lpszClassName, wc.hInstance);
        return 1;
    }
    printf("CUDA renderer initialized\n");

    // Load initial molecule
    molecule_build(g_currentMolecule, &g_molecule);

    // Setup Dear ImGui context
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGuiIO& io = ImGui::GetIO(); (void)io;
    io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;
    io.ConfigFlags |= ImGuiConfigFlags_DockingEnable;  // Uncomment if using docking branch

    // Setup Dear ImGui style
    ImGui::StyleColorsDark();

    // Load default font at larger size (13px * 1.3 = ~17px) for crisp text
    ImFontConfig fontConfig;
    fontConfig.SizePixels = 17.0f;
    io.Fonts->AddFontDefault(&fontConfig);
    // Larger font for tooltips/readability
    ImFontConfig tipCfg;
    tipCfg.SizePixels = 20.0f;
    g_fontTooltip = io.Fonts->AddFontDefault(&tipCfg);

    // Setup Platform/Renderer backends
    ImGui_ImplWin32_Init(hwnd);
    ImGui_ImplDX11_Init(g_pd3dDevice, g_pd3dDeviceContext);

    // Main loop
    bool running = true;
    while (running)
    {
        // Poll and handle messages
        MSG msg;
        while (PeekMessage(&msg, nullptr, 0U, 0U, PM_REMOVE))
        {
            TranslateMessage(&msg);
            DispatchMessage(&msg);
            if (msg.message == WM_QUIT)
                running = false;
        }
        if (!running)
            break;

        // Start the Dear ImGui frame
        ImGui_ImplDX11_NewFrame();
        ImGui_ImplWin32_NewFrame();
        ImGui::NewFrame();

		// Enable Docking with configurable lock options
		ImGuiDockNodeFlags dockFlags = 0;
		if (g_lockDocking) dockFlags |= ImGuiDockNodeFlags_NoUndocking;
		if (g_lockResize) dockFlags |= ImGuiDockNodeFlags_NoResize;
		ImGui::DockSpaceOverViewport(0, nullptr, dockFlags);

        // === MolVis UI ===

        // Main menu bar
        if (ImGui::BeginMainMenuBar())
        {
            if (ImGui::BeginMenu("File"))
            {
                if (ImGui::MenuItem("Exit", "Alt+F4"))
                    running = false;
                ImGui::EndMenu();
            }
            if (ImGui::BeginMenu("Molecule"))
            {
                // Organize by category
                for (int cat = 0; cat < CAT_COUNT; cat++)
                {
                    if (ImGui::BeginMenu(molecule_get_category_name(cat)))
                    {
                        for (int i = 0; i < molecule_get_count(); i++)
                        {
                            if (molecule_get_category(i) == cat)
                            {
                                if (ImGui::MenuItem(molecule_get_name(i), nullptr, g_currentMolecule == i))
                                {
                                    g_currentMolecule = i;
                                    molecule_build(g_currentMolecule, &g_molecule);
                                }
                            }
                        }
                        ImGui::EndMenu();
                    }
                }
                ImGui::EndMenu();
            }
            if (ImGui::BeginMenu("View"))
            {
                ImGui::MenuItem("Auto-Rotate", nullptr, &g_autoRotate);
                if (ImGui::MenuItem("Reset View"))
                {
                    g_rotX = 0.3f;
                    g_rotY = 0.0f;
                    g_zoom = 10.0f;
                }
                ImGui::EndMenu();
            }
            ImGui::EndMainMenuBar();
        }

        // Auto-rotate
        if (g_autoRotate) {
            if (g_rotateX) g_rotX += 0.01f * g_rotSpeed;
            if (g_rotateY) g_rotY += 0.01f * g_rotSpeed;
        }

        // Prepare molecule for rendering (apply view options without touching source data)
        Molecule renderMol;
        BuildRenderMolecule(g_molecule, renderMol, g_hideHydrogen);
        // Render molecule to CUDA texture
        renderer_render(g_pRenderer, &renderMol, g_rotX, g_rotY, g_zoom, g_offsetX, g_offsetY);

        // Render molecule name overlay using GPU text kernel
        // Format as "Name (Formula)" with subscript digits
        char displayText[128];
        if (g_molecule.formula[0]) {
            snprintf(displayText, sizeof(displayText), "%s (%s)", g_molecule.name, g_molecule.formula);
        } else {
            snprintf(displayText, sizeof(displayText), "%s", g_molecule.name);
        }
        renderer_render_text(g_pRenderer, displayText, 10, 10, 3);

        // Molecule Viewport Panel
        ImGui::SetNextWindowSize(ImVec2(820, 640), ImGuiCond_FirstUseEver);
        ImGui::Begin("Viewport", nullptr, ImGuiWindowFlags_NoScrollbar);
        {
            ImVec2 viewportSize = ImGui::GetContentRegionAvail();
            int vpWidth = (int)viewportSize.x;
            int vpHeight = (int)viewportSize.y;

            // Resize renderer if needed
            int curW, curH;
            renderer_get_size(g_pRenderer, &curW, &curH);
            if (vpWidth > 0 && vpHeight > 0 && (curW != vpWidth || curH != vpHeight))
            {
                renderer_resize(g_pRenderer, vpWidth, vpHeight);
            }

            // Display the CUDA-rendered texture
            ID3D11ShaderResourceView* srv = renderer_get_texture(g_pRenderer);
            if (srv)
            {
                ImGui::Image((ImTextureID)srv, viewportSize);
            }

            // Handle mouse interaction in viewport - only when hovering over the image
            if (ImGui::IsItemHovered())
            {
                if (ImGui::IsMouseDragging(ImGuiMouseButton_Left))
                {
                    ImVec2 delta = ImGui::GetMouseDragDelta(ImGuiMouseButton_Left);
                    ImGui::ResetMouseDragDelta(ImGuiMouseButton_Left);
                    g_rotY += delta.x * 0.01f;
                    g_rotX += delta.y * 0.01f;
                }
                float wheel = io.MouseWheel;
                if (wheel != 0.0f)
                {
                    g_zoom -= wheel * 0.5f;
                    if (g_zoom < 3.0f) g_zoom = 3.0f;
                    if (g_zoom > 30.0f) g_zoom = 30.0f;
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
                                    molecule_build(g_currentMolecule, &g_molecule);
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
                            molecule_build(g_currentMolecule, &g_molecule);
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

        // View Controls Panel
        ImGui::SetNextWindowSize(ImVec2(250, 220), ImGuiCond_FirstUseEver);
        ImGui::Begin("View Controls");
        {
            ImGui::Checkbox("Auto-Rotate", &g_autoRotate);
            if (g_autoRotate) {
                ImGui::SliderFloat("Speed", &g_rotSpeed, -3.0f, 3.0f, "%.1fx");
                ImGui::Checkbox("Rotate X", &g_rotateX);
                ImGui::SameLine();
                ImGui::Checkbox("Rotate Y", &g_rotateY);
            }
            ImGui::Separator();

            ImGui::Text("Camera:");
            ImGui::SliderFloat("Zoom", &g_zoom, 3.0f, 30.0f);
            ImGui::SliderFloat("Rotation X", &g_rotX, -3.14159f, 3.14159f);
            ImGui::SliderFloat("Rotation Y", &g_rotY, -3.14159f, 3.14159f);

            ImGui::Separator();
            ImGui::Text("Position:");
            ImGui::SliderFloat("Horizontal", &g_offsetX, -2.0f, 2.0f);
            ImGui::SliderFloat("Vertical", &g_offsetY, -2.0f, 2.0f);

            ImGui::Separator();
            if (ImGui::Button("Reset View"))
            {
                g_rotX = 0.3f;
                g_rotY = 0.0f;
                g_zoom = 10.0f;
                g_offsetX = 0.0f;
                g_offsetY = 0.0f;
            }
            ImGui::SameLine();
            if (ImGui::Button("Top View"))
            {
                g_rotX = 1.57f;
                g_rotY = 0.0f;
            }
            ImGui::SameLine();
            if (ImGui::Button("Side View"))
            {
                g_rotX = 0.0f;
                g_rotY = 1.57f;
            }

            ImGui::Separator();
            ImGui::Text("Layout:");
            ImGui::Checkbox("Lock Docking", &g_lockDocking);
            ImGui::Checkbox("Lock Resize", &g_lockResize);
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

        // Description Panel - Educational info about the current molecule
        ImGui::SetNextWindowSize(ImVec2(350, 280), ImGuiCond_FirstUseEver);
        ImGui::Begin("Description");
        {
            ImGui::PushStyleColor(ImGuiCol_Text, ImVec4(1.0f, 0.9f, 0.7f, 1.0f)); // Warm yellow text
            ImGui::TextWrapped("%s", g_molecule.name);
            ImGui::PopStyleColor();
            // Show formula below name if available
            if (g_molecule.formula[0]) {
                ImGui::PushStyleColor(ImGuiCol_Text, ImVec4(0.7f, 0.9f, 1.0f, 1.0f)); // Light blue text
                ImGui::Text("Formula: %s", g_molecule.formula);
                ImGui::PopStyleColor();
            }
            ImGui::Separator();
            ImGui::Spacing();
            ImGui::TextWrapped("%s", molecule_get_long_description(g_currentMolecule));
            ImGui::Spacing();
            ImGui::Separator();

            // Molecule Info
            ImGui::Text("Atoms: %d", g_molecule.numAtoms);
            ImGui::Text("Bonds: %d", g_molecule.numBonds);
            ImGui::Text("Category: %s", molecule_get_category_name(molecule_get_category(g_currentMolecule)));
        }
        ImGui::End();

        // Molecule Options Panel
        ImGui::SetNextWindowSize(ImVec2(240, 80), ImGuiCond_FirstUseEver);
        ImGui::Begin("Molecule Options");
        {
            ImGui::TextDisabled("Rendering options");
            ImGui::Separator();
            ImGui::Checkbox("Hide Hydrogen (H)", &g_hideHydrogen);
            ImGui::TextDisabled("Does not affect Periodic Table");
        }
        ImGui::End();

        // Periodic Table Panel
        ImGui::SetNextWindowSize(ImVec2(860, 360), ImGuiCond_FirstUseEver);
        ImGui::Begin("Periodic Table");
        {
            // Define all 118 elements with positions (period, group). Map to our internal types when applicable; otherwise -1.
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
                {"I",  "Iodine",          5, 17, ATOM_I},     {"Xe", "Xenon",           5, 18, ATOM_XE},
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
                case ATOM_RE: return ImVec4(0.51f, 0.51f, 0.56f, 1.0f); // metallic gray
                case ATOM_XE: return ImVec4(0.26f, 0.62f, 0.69f, 1.0f); // cyan (noble gas)
                default:      return ImVec4(0.32f, 0.32f, 0.36f, 1.0f);
                }
            };

            // Compute presence of each supported type in current molecule (independent of rendering options)
            bool present[ATOM_TYPE_COUNT];
            ComputePresenceFlags(g_molecule, g_currentMolecule, present);

            // Layout constants
            const ImVec2 cellSize(36, 30);
            const ImVec2 spacing(6, 6);
            const float topMargin = 18.0f; // space for group numbers
            const float leftMargin = 18.0f;

            ImVec2 origin = ImGui::GetCursorScreenPos();

            // Group number row (no grid)
            for (int c = 1; c <= 18; ++c) {
                ImVec2 pos(origin.x + leftMargin + (c-1)*(cellSize.x + spacing.x) + 10.0f,
                           origin.y);
                ImGui::SetCursorScreenPos(pos);
                ImGui::Text("%d", c);
            }

            // Draw each element as a colored button placed by (period, group)
            for (const auto& e : elements) {
                float x = leftMargin + (e.group - 1) * (cellSize.x + spacing.x);
                float y = topMargin + (e.period - 1) * (cellSize.y + spacing.y);
                ImGui::SetCursorScreenPos(ImVec2(origin.x + x, origin.y + y));

                        bool isPresent = (e.atomType >= 0 && e.atomType < ATOM_TYPE_COUNT) ? present[e.atomType] : false;
                // Gray for all elements by default; use CPK color only when present
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
                    if (g_fontTooltip) ImGui::PushFont(g_fontTooltip);
                    ImGui::TextUnformatted(e.name);
                    if (isPresent) {
                        ImGui::Spacing();
                        ImGui::TextColored(ImVec4(1.0f, 0.95f, 0.3f, 1.0f), "Present in molecule");
                    }
                    if (g_fontTooltip) ImGui::PopFont();
                    ImGui::PopStyleColor();
                    ImGui::EndTooltip();
                }

                ImGui::PopStyleColor(5); // Text, ButtonActive, ButtonHovered, Button, Border
                ImGui::PopStyleVar(2);
                ImGui::PopID();
            }

            ImGui::Dummy(ImVec2(0, topMargin + 7*(cellSize.y + spacing.y))); // ensure window content height
            ImGui::Separator();
            ImGui::TextDisabled("Highlighted elements are present in the selected molecule.");
        }
        ImGui::End();

        // CPK Color Reference Panel
        ImGui::SetNextWindowSize(ImVec2(280, 450), ImGuiCond_FirstUseEver);
        ImGui::Begin("CPK Color Reference");
        {
            ImGui::Text("CPK Atom Colors:");
            ImGui::Separator();

            // Define colors matching the CUDA renderer
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
                {0.26f, 0.62f, 0.69f, 1.0f},  // Xe - cyan (noble gas)
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
                "Xe - Xenon (Cyan)",
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
        const float clear_color[4] = { 0.1f, 0.1f, 0.15f, 1.0f };
        g_pd3dDeviceContext->OMSetRenderTargets(1, &g_mainRenderTargetView, nullptr);
        g_pd3dDeviceContext->ClearRenderTargetView(g_mainRenderTargetView, clear_color);
        ImGui_ImplDX11_RenderDrawData(ImGui::GetDrawData());

        g_pSwapChain->Present(1, 0); // Present with vsync
    }

    // Save window config to molvis.ini
    SaveWindowConfig();

    // Cleanup
    renderer_cleanup(g_pRenderer);
    ImGui_ImplDX11_Shutdown();
    ImGui_ImplWin32_Shutdown();
    ImGui::DestroyContext();

    CleanupDeviceD3D();
    DestroyWindow(hwnd);
    UnregisterClassW(wc.lpszClassName, wc.hInstance);

    return 0;
}

// === Molecule Browser Helper Functions ===

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

// === DirectX 11 Helper Functions ===

bool CreateDeviceD3D(HWND hWnd)
{
    DXGI_SWAP_CHAIN_DESC sd;
    ZeroMemory(&sd, sizeof(sd));
    sd.BufferCount = 2;
    sd.BufferDesc.Width = 0;
    sd.BufferDesc.Height = 0;
    sd.BufferDesc.Format = DXGI_FORMAT_R8G8B8A8_UNORM;
    sd.BufferDesc.RefreshRate.Numerator = 60;
    sd.BufferDesc.RefreshRate.Denominator = 1;
    sd.Flags = DXGI_SWAP_CHAIN_FLAG_ALLOW_MODE_SWITCH;
    sd.BufferUsage = DXGI_USAGE_RENDER_TARGET_OUTPUT;
    sd.OutputWindow = hWnd;
    sd.SampleDesc.Count = 1;
    sd.SampleDesc.Quality = 0;
    sd.Windowed = TRUE;
    sd.SwapEffect = DXGI_SWAP_EFFECT_DISCARD;

    UINT createDeviceFlags = 0;
    D3D_FEATURE_LEVEL featureLevel;
    const D3D_FEATURE_LEVEL featureLevelArray[2] = { D3D_FEATURE_LEVEL_11_0, D3D_FEATURE_LEVEL_10_0 };

    HRESULT res = D3D11CreateDeviceAndSwapChain(
        nullptr, D3D_DRIVER_TYPE_HARDWARE, nullptr, createDeviceFlags,
        featureLevelArray, 2, D3D11_SDK_VERSION,
        &sd, &g_pSwapChain, &g_pd3dDevice, &featureLevel, &g_pd3dDeviceContext
    );
    if (res != S_OK)
        return false;

    CreateRenderTarget();
    return true;
}

void CleanupDeviceD3D()
{
    CleanupRenderTarget();
    if (g_pSwapChain)       { g_pSwapChain->Release(); g_pSwapChain = nullptr; }
    if (g_pd3dDeviceContext){ g_pd3dDeviceContext->Release(); g_pd3dDeviceContext = nullptr; }
    if (g_pd3dDevice)       { g_pd3dDevice->Release(); g_pd3dDevice = nullptr; }
}

void CreateRenderTarget()
{
    ID3D11Texture2D* pBackBuffer;
    g_pSwapChain->GetBuffer(0, IID_PPV_ARGS(&pBackBuffer));
    g_pd3dDevice->CreateRenderTargetView(pBackBuffer, nullptr, &g_mainRenderTargetView);
    pBackBuffer->Release();
}

void CleanupRenderTarget()
{
    if (g_mainRenderTargetView) { g_mainRenderTargetView->Release(); g_mainRenderTargetView = nullptr; }
}

// Forward declare message handler from imgui_impl_win32.cpp
extern IMGUI_IMPL_API LRESULT ImGui_ImplWin32_WndProcHandler(HWND hWnd, UINT msg, WPARAM wParam, LPARAM lParam);

LRESULT WINAPI WndProc(HWND hWnd, UINT msg, WPARAM wParam, LPARAM lParam)
{
    if (ImGui_ImplWin32_WndProcHandler(hWnd, msg, wParam, lParam))
        return true;

    switch (msg)
    {
    case WM_SIZE:
        if (g_pd3dDevice != nullptr && wParam != SIZE_MINIMIZED)
        {
            CleanupRenderTarget();
            g_pSwapChain->ResizeBuffers(0, (UINT)LOWORD(lParam), (UINT)HIWORD(lParam), DXGI_FORMAT_UNKNOWN, 0);
            CreateRenderTarget();
        }
        return 0;
    case WM_SYSCOMMAND:
        if ((wParam & 0xfff0) == SC_KEYMENU) // Disable ALT application menu
            return 0;
        break;
    case WM_CLOSE:
        SaveWindowConfig();  // Save before closing
        break;
    case WM_DESTROY:
        CleanupWindowIcons();  // Clean up icon handles
        PostQuitMessage(0);
        return 0;
    }
    return DefWindowProcW(hWnd, msg, wParam, lParam);
}
