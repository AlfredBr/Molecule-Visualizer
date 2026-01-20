# MolVis - CUDA Molecular Visualization
# Makefile for Windows with NVCC and MSVC host compiler
#
# Target: RTX 3080 (Compute Capability 8.6)
#
# IMPORTANT: This Makefile must be run from a Visual Studio Developer Command Prompt
# or after running vcvarsall.bat to set up the MSVC environment.
# Use the provided build.bat or setup_env.ps1 scripts for automatic setup.

#==============================================================================
# Configuration
#==============================================================================

# Application name
APP_NAME = molvis

# Compiler
NVCC = nvcc

# NVCC flags
# -arch: Target GPU architecture (sm_86 for RTX 3080, adjust for your GPU)
# -allow-unsupported-compiler: Allow newer VS versions not yet officially supported
NVCCFLAGS = -O3 -arch=sm_86 -allow-unsupported-compiler

# Debug flags (uncomment for debug builds)
# NVCCFLAGS = -g -G -arch=sm_86 -allow-unsupported-compiler

# Windows libraries
LIBS = gdi32.lib user32.lib

# Source files
CUDA_SOURCES = cuda_molecule.cu
HEADERS = win32_display.h

# Output
TARGET = $(APP_NAME).exe

#==============================================================================
# Dear ImGui Integration (future)
#==============================================================================

# IMGUI_DIR = third_party/imgui
# IMGUI_SOURCES = $(IMGUI_DIR)/imgui.cpp \
#                 $(IMGUI_DIR)/imgui_draw.cpp \
#                 $(IMGUI_DIR)/imgui_tables.cpp \
#                 $(IMGUI_DIR)/imgui_widgets.cpp \
#                 $(IMGUI_DIR)/backends/imgui_impl_win32.cpp \
#                 $(IMGUI_DIR)/backends/imgui_impl_dx11.cpp

#==============================================================================
# Build Targets
#==============================================================================

.PHONY: all clean help run

all: $(TARGET)
	@echo.
	@echo === MolVis built successfully! ===
	@echo Run with: $(TARGET)

$(TARGET): $(CUDA_SOURCES) $(HEADERS)
	$(NVCC) $(NVCCFLAGS) -o $(TARGET) $(CUDA_SOURCES) $(LIBS)

run: $(TARGET)
	$(TARGET)

clean:
	@if exist $(TARGET) del /Q $(TARGET)
	@if exist *.obj del /Q *.obj
	@if exist *.exp del /Q *.exp
	@if exist *.lib del /Q *.lib
	@if exist *.pdb del /Q *.pdb
	@echo Cleaned build artifacts.

help:
	@echo.
	@echo MolVis Build System
	@echo ====================
	@echo.
	@echo Targets:
	@echo   all     - Build the application (default)
	@echo   run     - Build and run the application
	@echo   clean   - Remove build artifacts
	@echo   help    - Show this help message
	@echo.
	@echo Configuration:
	@echo   GPU Architecture: sm_86 (RTX 3080)
	@echo   Compiler: NVCC with MSVC host
	@echo.
	@echo Prerequisites:
	@echo   - NVIDIA CUDA Toolkit
	@echo   - Visual Studio with C++ workload
	@echo   - Run from VS Developer Command Prompt or use build.bat
	@echo.
