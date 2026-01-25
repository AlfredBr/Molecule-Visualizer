#pragma once
// Window Icon Management
// Sets custom window and taskbar icons from embedded resources on Windows

#ifdef _WIN32
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#include <windows.h>

// Global icon handles for cleanup
static HICON g_hIconBig = NULL;
static HICON g_hIconSmall = NULL;

// Sets the window icon for title bar and taskbar
// Must be called after CreateWindow()
// Uses the embedded icon resource (ID 1) from resource.rc
inline void SetWindowIconFromResource(HWND hwnd)
{
	HINSTANCE hInstance = GetModuleHandle(NULL);

	// Load icons from embedded resource
	g_hIconBig = (HICON)LoadImageA(hInstance, MAKEINTRESOURCEA(1), IMAGE_ICON, 32, 32, 0);
	g_hIconSmall = (HICON)LoadImageA(hInstance, MAKEINTRESOURCEA(1), IMAGE_ICON, 16, 16, 0);

	// Set icon on window instance
	if (g_hIconBig)
		SendMessage(hwnd, WM_SETICON, ICON_BIG, (LPARAM)g_hIconBig);
	if (g_hIconSmall)
		SendMessage(hwnd, WM_SETICON, ICON_SMALL, (LPARAM)g_hIconSmall);

	// Also set icon on window class (needed for taskbar on some Windows versions)
	if (g_hIconBig)
		SetClassLongPtr(hwnd, GCLP_HICON, (LONG_PTR)g_hIconBig);
	if (g_hIconSmall)
		SetClassLongPtr(hwnd, GCLP_HICONSM, (LONG_PTR)g_hIconSmall);
}

// Cleanup function to destroy icon handles
// Should be called during window destruction
inline void CleanupWindowIcons()
{
	if (g_hIconBig)
	{
		DestroyIcon(g_hIconBig);
		g_hIconBig = NULL;
	}
	if (g_hIconSmall)
	{
		DestroyIcon(g_hIconSmall);
		g_hIconSmall = NULL;
	}
}

#else
// Non-Windows stubs
inline void SetWindowIconFromResource(void* windowHandle) { (void)windowHandle; }
inline void CleanupWindowIcons() {}
#endif
