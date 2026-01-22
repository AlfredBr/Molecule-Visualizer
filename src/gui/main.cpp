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
// ...existing code...
