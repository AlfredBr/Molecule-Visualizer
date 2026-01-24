# MolVis - CUDA Molecular Visualization
# Makefile for Windows with MSVC and NVCC
#
# Usage:
#   build.bat         - Build the application (recommended)
#   build.bat clean   - Clean build artifacts
#   build.bat run     - Build and run
#
# Or from VS Developer Command Prompt:
#   nmake
#   nmake clean
#   nmake run

#==============================================================================
# Configuration
#==============================================================================

APP_NAME     = molvis
TARGET       = $(APP_NAME).exe

# Compilers
CXX          = cl
NVCC         = nvcc

# Directories
SRC_DIR      = src
BUILD_DIR    = build
IMGUI_DIR    = third_party/imgui
LEGACY_DIR   = legacy

# CUDA paths
CUDA_INC     = $(CUDA_PATH)/include
CUDA_LIB     = $(CUDA_PATH)/lib/x64

#==============================================================================
# Compiler Flags
#==============================================================================

# MSVC flags
CXXFLAGS     = /nologo /W3 /O2 /EHsc /std:c++17 /D_CRT_SECURE_NO_WARNINGS

# NVCC flags (sm_86 = RTX 3080, adjust for your GPU)
NVCCFLAGS    = -O3 -arch=sm_86 -allow-unsupported-compiler

# Include paths
INCLUDES     = /I. /I$(SRC_DIR) /Iinclude /I$(IMGUI_DIR) /I$(IMGUI_DIR)/backends /I"$(CUDA_INC)"
CUDA_INC_FLAGS = -I. -I$(SRC_DIR) -Iinclude -I$(IMGUI_DIR) -I$(IMGUI_DIR)/backends

# Libraries
LIBS         = user32.lib gdi32.lib d3d11.lib dxgi.lib d3dcompiler.lib cudart.lib
LIB_PATHS    = /LIBPATH:"$(CUDA_LIB)"

#==============================================================================
# Source Files
#==============================================================================

# Application sources
APP_SOURCES  = \
	$(SRC_DIR)/gui/main_windows.cpp \
	$(SRC_DIR)/molecule/molecule_db.cpp

# CUDA sources
CUDA_SOURCES = $(SRC_DIR)/renderer/cuda_renderer.cu

# Dear ImGui sources
IMGUI_SOURCES = \
	$(IMGUI_DIR)/imgui.cpp \
	$(IMGUI_DIR)/imgui_draw.cpp \
	$(IMGUI_DIR)/imgui_tables.cpp \
	$(IMGUI_DIR)/imgui_widgets.cpp \
	$(IMGUI_DIR)/backends/imgui_impl_win32.cpp \
	$(IMGUI_DIR)/backends/imgui_impl_dx11.cpp

# Object files
CUDA_OBJ     = $(BUILD_DIR)/cuda_renderer.obj
RES_OBJ      = $(BUILD_DIR)/resource.res

#==============================================================================
# Build Targets
#==============================================================================

.PHONY: all clean run legacy help dirs

all: dirs $(TARGET)
	@echo.
	@echo ========================================
	@echo   MolVis built successfully!
	@echo   Run with: $(TARGET)
	@echo ========================================

dirs:
	@if not exist $(BUILD_DIR) mkdir $(BUILD_DIR)

# Compile CUDA sources
$(CUDA_OBJ): $(CUDA_SOURCES)
	$(NVCC) $(NVCCFLAGS) $(CUDA_INC_FLAGS) -c $** -o $@

# Compile Windows resource file (for app icon)
$(RES_OBJ): resource.rc
	rc /nologo /fo $@ $**

# Link everything
$(TARGET): $(CUDA_OBJ) $(RES_OBJ) $(APP_SOURCES) $(IMGUI_SOURCES)
	$(CXX) $(CXXFLAGS) $(INCLUDES) /Fe:$@ /Fo:$(BUILD_DIR)/ \
		$(APP_SOURCES) $(IMGUI_SOURCES) $(CUDA_OBJ) \
		/link $(LIB_PATHS) $(LIBS) $(RES_OBJ)

run: all
	@echo Starting MolVis...
	@$(TARGET)

#==============================================================================
# Legacy Build (original CUDA-only version)
#==============================================================================

legacy: $(LEGACY_DIR)/cuda_molecule.cu $(LEGACY_DIR)/win32_display.h
	$(NVCC) $(NVCCFLAGS) -o $(LEGACY_DIR)/molvis_legacy.exe \
		$(LEGACY_DIR)/cuda_molecule.cu gdi32.lib user32.lib
	@echo Built: $(LEGACY_DIR)/molvis_legacy.exe

#==============================================================================
# Clean
#==============================================================================

clean:
	@echo Cleaning build artifacts...
	@if exist $(TARGET) del /Q $(TARGET)
	@if exist $(BUILD_DIR)\*.obj del /Q $(BUILD_DIR)\*.obj
	@if exist $(BUILD_DIR)\*.res del /Q $(BUILD_DIR)\*.res
	@if exist *.obj del /Q *.obj
	@if exist *.exp del /Q *.exp
	@if exist *.lib del /Q *.lib
	@if exist *.pdb del /Q *.pdb
	@if exist *.ilk del /Q *.ilk
	@if exist imgui.ini del /Q imgui.ini
	@echo Done.

#==============================================================================
# Help
#==============================================================================

help:
	@echo.
	@echo MolVis Build System
	@echo ========================================
	@echo.
	@echo Targets:
	@echo   all      Build the application (default)
	@echo   run      Build and run
	@echo   clean    Remove build artifacts
	@echo   legacy   Build original CUDA-only version
	@echo   help     Show this message
	@echo.
	@echo Configuration:
	@echo   GPU:     sm_86 (RTX 3080)
	@echo   Compiler: MSVC + NVCC
	@echo.
	@echo Requirements:
	@echo   - Visual Studio with C++ workload
	@echo   - NVIDIA CUDA Toolkit
	@echo   - Run from VS Developer Command Prompt
	@echo     or use build.bat
	@echo.
