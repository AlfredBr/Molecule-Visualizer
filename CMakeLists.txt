# MolVis - Cross-Platform Molecular Visualization
# CMake build system for Windows and macOS
#
# Usage:
#   mkdir build && cd build
#   cmake .. && cmake --build .

cmake_minimum_required(VERSION 3.20)
project(MolVis VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Output directory
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

#==============================================================================
# Options
#==============================================================================
option(MOLVIS_BUILD_WINDOWS "Build Windows version (CUDA + DX11)" OFF)
option(MOLVIS_BUILD_MACOS "Build macOS version (Metal)" OFF)

# Auto-detect platform
if(NOT MOLVIS_BUILD_WINDOWS AND NOT MOLVIS_BUILD_MACOS)
    if(APPLE)
        set(MOLVIS_BUILD_MACOS ON)
    elseif(WIN32)
        set(MOLVIS_BUILD_WINDOWS ON)
    endif()
endif()

message(STATUS "MolVis Build Configuration:")
message(STATUS "  Platform: ${CMAKE_SYSTEM_NAME}")
message(STATUS "  Build Windows: ${MOLVIS_BUILD_WINDOWS}")
message(STATUS "  Build macOS: ${MOLVIS_BUILD_MACOS}")

#==============================================================================
# Dear ImGui
#==============================================================================
set(IMGUI_DIR ${CMAKE_SOURCE_DIR}/third_party/imgui)

set(IMGUI_SOURCES
    ${IMGUI_DIR}/imgui.cpp
    ${IMGUI_DIR}/imgui_draw.cpp
    ${IMGUI_DIR}/imgui_tables.cpp
    ${IMGUI_DIR}/imgui_widgets.cpp
    ${IMGUI_DIR}/imgui_demo.cpp
)

#==============================================================================
# Common Sources (Platform-agnostic)
#==============================================================================
set(COMMON_SOURCES
    ${CMAKE_SOURCE_DIR}/src/molecule/molecule_db.cpp
)

set(COMMON_HEADERS
    ${CMAKE_SOURCE_DIR}/src/molecule/molecule_db.h
    ${CMAKE_SOURCE_DIR}/src/molecule/molecule_types.h
    ${CMAKE_SOURCE_DIR}/src/renderer/renderer_interface.h
)

#==============================================================================
# macOS Build (Metal + SDL2)
#==============================================================================
if(MOLVIS_BUILD_MACOS)
    message(STATUS "Configuring macOS build with Metal + SDL2...")

    # macOS deployment target
    set(CMAKE_OSX_DEPLOYMENT_TARGET "11.0" CACHE STRING "Minimum macOS version")

    # Enable Objective-C++
    enable_language(OBJCXX)
    set(CMAKE_OBJCXX_STANDARD 17)

    # Find required frameworks
    find_library(METAL_FRAMEWORK Metal REQUIRED)
    find_library(METALKIT_FRAMEWORK MetalKit REQUIRED)
    find_library(COCOA_FRAMEWORK Cocoa REQUIRED)
    find_library(QUARTZCORE_FRAMEWORK QuartzCore REQUIRED)
    find_library(IOKIT_FRAMEWORK IOKit REQUIRED)
    find_library(GAMECONTROLLER_FRAMEWORK GameController REQUIRED)

    # Find SDL2
    find_package(SDL2 REQUIRED)

    # ImGui backends for macOS
    set(IMGUI_BACKEND_SOURCES
        ${IMGUI_DIR}/backends/imgui_impl_sdl2.cpp
        ${IMGUI_DIR}/backends/imgui_impl_metal.mm
    )

    # macOS-specific sources
    set(MACOS_SOURCES
        ${CMAKE_SOURCE_DIR}/src/gui/main_mac.mm
        ${CMAKE_SOURCE_DIR}/src/renderer/metal_renderer.mm
    )

    # Metal shader source
    set(METAL_SHADER_SOURCE ${CMAKE_SOURCE_DIR}/src/renderer/shaders/MoleculeShaders.metal)
    set(METAL_SHADER_OUTPUT ${CMAKE_BINARY_DIR}/MoleculeShaders.metallib)

    # Compile Metal shaders to metallib
    add_custom_command(
        OUTPUT ${METAL_SHADER_OUTPUT}
        COMMAND xcrun -sdk macosx metal -c ${METAL_SHADER_SOURCE} -o ${CMAKE_BINARY_DIR}/MoleculeShaders.air
        COMMAND xcrun -sdk macosx metallib ${CMAKE_BINARY_DIR}/MoleculeShaders.air -o ${METAL_SHADER_OUTPUT}
        DEPENDS ${METAL_SHADER_SOURCE}
        COMMENT "Compiling Metal shaders..."
    )

    # Custom target for shaders
    add_custom_target(MolVisShaders DEPENDS ${METAL_SHADER_OUTPUT})

    # Create executable
    add_executable(MolVis MACOSX_BUNDLE
        ${MACOS_SOURCES}
        ${COMMON_SOURCES}
        ${IMGUI_SOURCES}
        ${IMGUI_BACKEND_SOURCES}
    )

    # Ensure shaders are built before the app
    add_dependencies(MolVis MolVisShaders)

    # Include directories
    target_include_directories(MolVis PRIVATE
        ${CMAKE_SOURCE_DIR}
        ${CMAKE_SOURCE_DIR}/src
        ${CMAKE_SOURCE_DIR}/include
        ${IMGUI_DIR}
        ${IMGUI_DIR}/backends
        ${SDL2_INCLUDE_DIRS}
    )

    # Link frameworks and libraries
    target_link_libraries(MolVis PRIVATE
        ${METAL_FRAMEWORK}
        ${METALKIT_FRAMEWORK}
        ${COCOA_FRAMEWORK}
        ${QUARTZCORE_FRAMEWORK}
        ${IOKIT_FRAMEWORK}
        ${GAMECONTROLLER_FRAMEWORK}
        ${SDL2_LIBRARIES}
    )

    # macOS bundle properties
    set_target_properties(MolVis PROPERTIES
        MACOSX_BUNDLE_INFO_PLIST ${CMAKE_SOURCE_DIR}/platform/macos/Info.plist
        MACOSX_BUNDLE_BUNDLE_NAME "MolVis"
        MACOSX_BUNDLE_BUNDLE_VERSION "${PROJECT_VERSION}"
        MACOSX_BUNDLE_SHORT_VERSION_STRING "${PROJECT_VERSION}"
        MACOSX_BUNDLE_GUI_IDENTIFIER "com.molvis.app"
    )

    # Copy Metal shader library and source to bundle Resources
    add_custom_command(TARGET MolVis POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E make_directory "$<TARGET_BUNDLE_DIR:MolVis>/Contents/Resources"
        COMMAND ${CMAKE_COMMAND} -E copy ${METAL_SHADER_OUTPUT} "$<TARGET_BUNDLE_DIR:MolVis>/Contents/Resources/default.metallib"
        COMMAND ${CMAKE_COMMAND} -E copy ${METAL_SHADER_SOURCE} "$<TARGET_BUNDLE_DIR:MolVis>/Contents/Resources/MoleculeShaders.metal"
        COMMENT "Copying Metal shaders to app bundle..."
    )

    # Copy application icon to bundle Resources (if it exists)
    set(MACOS_ICON_FILE ${CMAKE_SOURCE_DIR}/platform/macos/AppIcon.icns)
    if(EXISTS ${MACOS_ICON_FILE})
        add_custom_command(TARGET MolVis POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy ${MACOS_ICON_FILE} "$<TARGET_BUNDLE_DIR:MolVis>/Contents/Resources/AppIcon.icns"
            COMMENT "Copying application icon to app bundle..."
        )
    endif()

    # Bundle SDL2 dylib for distribution (so users don't need Homebrew)
    add_custom_command(TARGET MolVis POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E make_directory "$<TARGET_BUNDLE_DIR:MolVis>/Contents/Frameworks"
        COMMENT "Creating Frameworks directory..."
    )

    # Find and copy SDL2 dylib, then fix the library path
    add_custom_command(TARGET MolVis POST_BUILD
        COMMAND ${CMAKE_SOURCE_DIR}/scripts/bundle_sdl2.sh "$<TARGET_BUNDLE_DIR:MolVis>"
        COMMENT "Bundling SDL2 library..."
    )

    message(STATUS "macOS build configured successfully")
endif()

#==============================================================================
# Windows Build (CUDA + DX11) - Keep existing functionality
#==============================================================================
if(MOLVIS_BUILD_WINDOWS)
    message(STATUS "Configuring Windows build with CUDA + DX11...")

    # Find CUDA
    find_package(CUDAToolkit REQUIRED)
    enable_language(CUDA)

    # ImGui backends for Windows
    set(IMGUI_BACKEND_SOURCES
        ${IMGUI_DIR}/backends/imgui_impl_win32.cpp
        ${IMGUI_DIR}/backends/imgui_impl_dx11.cpp
    )

    # Windows-specific sources
    set(WINDOWS_SOURCES
        ${CMAKE_SOURCE_DIR}/src/gui/main.cpp
        ${CMAKE_SOURCE_DIR}/src/renderer/cuda_renderer.cu
    )

    # Create executable
    add_executable(MolVis WIN32
        ${WINDOWS_SOURCES}
        ${COMMON_SOURCES}
        ${IMGUI_SOURCES}
        ${IMGUI_BACKEND_SOURCES}
        ${CMAKE_SOURCE_DIR}/resource.rc
    )

    # CUDA architecture
    set_target_properties(MolVis PROPERTIES
        CUDA_ARCHITECTURES "86"
    )

    # Include directories
    target_include_directories(MolVis PRIVATE
        ${CMAKE_SOURCE_DIR}
        ${CMAKE_SOURCE_DIR}/src
        ${CMAKE_SOURCE_DIR}/include
        ${IMGUI_DIR}
        ${IMGUI_DIR}/backends
        ${CUDAToolkit_INCLUDE_DIRS}
    )

    # Link libraries
    target_link_libraries(MolVis PRIVATE
        CUDA::cudart
        d3d11
        dxgi
        d3dcompiler
        user32
        gdi32
    )

    message(STATUS "Windows build configured successfully")
endif()

#==============================================================================
# Installation
#==============================================================================
if(APPLE)
    install(TARGETS MolVis BUNDLE DESTINATION .)
else()
    install(TARGETS MolVis RUNTIME DESTINATION bin)
endif()
