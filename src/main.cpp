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

#include "renderer/cuda_renderer.h"
#include "molecule/molecule_db.h"

// DirectX 11 globals
static ID3D11Device*            g_pd3dDevice = nullptr;
static ID3D11DeviceContext*     g_pd3dDeviceContext = nullptr;
static IDXGISwapChain*          g_pSwapChain = nullptr;
static ID3D11RenderTargetView*  g_mainRenderTargetView = nullptr;

// CUDA Renderer
static CudaRenderer*            g_pRenderer = nullptr;
static Molecule                 g_molecule = {};
static int                      g_currentMolecule = 0;
static float                    g_rotX = 0.3f;
static float                    g_rotY = 0.0f;
static float                    g_zoom = 10.0f;
static float                    g_offsetX = 0.0f;
static float                    g_offsetY = 0.0f;
static bool                     g_autoRotate = true;
static float                    g_rotSpeed = 1.0f;
static bool                     g_rotateX = false;
static bool                     g_rotateY = true;

// Forward declarations
bool CreateDeviceD3D(HWND hWnd);
void CleanupDeviceD3D();
void CreateRenderTarget();
void CleanupRenderTarget();
LRESULT WINAPI WndProc(HWND hWnd, UINT msg, WPARAM wParam, LPARAM lParam);

int WINAPI WinMain(HINSTANCE hInstance, HINSTANCE hPrevInstance, LPSTR lpCmdLine, int nCmdShow)
{
    // Allocate console for debug output
    AllocConsole();
    freopen("CONOUT$", "w", stdout);
    freopen("CONOUT$", "w", stderr);
    printf("MolVis starting...\n");

    // Create application window
    WNDCLASSEXW wc = {
        sizeof(wc), CS_CLASSDC, WndProc, 0L, 0L,
        GetModuleHandle(nullptr), nullptr, nullptr, nullptr, nullptr,
        L"MolVis", nullptr
    };
    RegisterClassExW(&wc);

    HWND hwnd = CreateWindowW(
        wc.lpszClassName, L"MolVis - Molecular Visualization",
        WS_OVERLAPPEDWINDOW,
        100, 100, 1280, 800,
        nullptr, nullptr, wc.hInstance, nullptr
    );

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
    // io.ConfigFlags |= ImGuiConfigFlags_DockingEnable;  // Uncomment if using docking branch

    // Setup Dear ImGui style
    ImGui::StyleColorsDark();

    // Load default font at larger size (13px * 1.3 = ~17px) for crisp text
    ImFontConfig fontConfig;
    fontConfig.SizePixels = 17.0f;
    io.Fonts->AddFontDefault(&fontConfig);

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

        // Render molecule to CUDA texture
        renderer_render(g_pRenderer, &g_molecule, g_rotX, g_rotY, g_zoom, g_offsetX, g_offsetY);

        // Render molecule name overlay using GPU text kernel
        renderer_render_text(g_pRenderer, g_molecule.name, 10, 10, 3);

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

            // Handle mouse interaction in viewport
            if (ImGui::IsWindowHovered())
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

        // Molecule Info Panel
        ImGui::SetNextWindowSize(ImVec2(320, 350), ImGuiCond_FirstUseEver);
        ImGui::Begin("Molecule Info");
        {
            ImGui::Text("Name: %s", g_molecule.name);
            ImGui::TextWrapped("Description: %s", molecule_get_description(g_currentMolecule));
            ImGui::Separator();

            ImGui::Text("Atoms: %d", g_molecule.numAtoms);
            ImGui::Text("Bonds: %d", g_molecule.numBonds);
            ImGui::Text("Category: %s", molecule_get_category_name(molecule_get_category(g_currentMolecule)));
            ImGui::Separator();

            // Molecule browser by category
            ImGui::Text("Browse Molecules:");
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
        const float clear_color[4] = { 0.1f, 0.1f, 0.15f, 1.0f };
        g_pd3dDeviceContext->OMSetRenderTargets(1, &g_mainRenderTargetView, nullptr);
        g_pd3dDeviceContext->ClearRenderTargetView(g_mainRenderTargetView, clear_color);
        ImGui_ImplDX11_RenderDrawData(ImGui::GetDrawData());

        g_pSwapChain->Present(1, 0); // Present with vsync
    }

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
    case WM_DESTROY:
        PostQuitMessage(0);
        return 0;
    }
    return DefWindowProcW(hWnd, msg, wParam, lParam);
}
