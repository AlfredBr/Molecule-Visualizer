#!/bin/bash
# MolVis - macOS Build Script
#
# Usage:
#   ./build_mac.sh          Build the application
#   ./build_mac.sh clean    Clean build artifacts
#   ./build_mac.sh run      Build and run
#
# Requirements:
#   - Xcode Command Line Tools: xcode-select --install
#   - CMake: brew install cmake
#   - SDL2: brew install sdl2

set -e  # Exit on error

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BUILD_DIR="$SCRIPT_DIR/build_mac"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

print_status() {
    echo -e "${GREEN}==>${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}Warning:${NC} $1"
}

print_error() {
    echo -e "${RED}Error:${NC} $1"
}

# Check for required tools
check_dependencies() {
    print_status "Checking dependencies..."

    # Check for Xcode Command Line Tools
    if ! xcode-select -p &>/dev/null; then
        print_error "Xcode Command Line Tools not found"
        echo "Install with: xcode-select --install"
        exit 1
    fi

    # Check for Metal compiler (part of Xcode toolchain)
    if ! xcrun -find metal &>/dev/null; then
        print_error "Metal compiler not found"
        echo "Install the Metal toolchain with:"
        echo "  xcodebuild -downloadComponent MetalToolchain"
        echo "  xcodebuild -runFirstLaunch"
        exit 1
    fi

    # Check for CMake
    if ! command -v cmake &>/dev/null; then
        print_error "CMake not found"
        echo "Install with: brew install cmake"
        exit 1
    fi

    # Check for SDL2
    if ! pkg-config --exists sdl2 2>/dev/null; then
        # Try alternate location
        if [ ! -d "/usr/local/include/SDL2" ] && [ ! -d "/opt/homebrew/include/SDL2" ]; then
            print_warning "SDL2 not found via pkg-config"
            echo "Install with: brew install sdl2"
            echo "Continuing anyway - CMake may find it..."
        fi
    fi

    print_status "Dependencies OK"
}

# Clean build
clean() {
    print_status "Cleaning build directory..."
    rm -rf "$BUILD_DIR"
    print_status "Clean complete"
}

# Build
build() {
    check_dependencies

    # Generate macOS icon from PNG if available
    if [ -f "$SCRIPT_DIR/molvis-icon.png" ]; then
        if [ ! -f "$SCRIPT_DIR/platform/macos/AppIcon.icns" ] || \
           [ "$SCRIPT_DIR/molvis-icon.png" -nt "$SCRIPT_DIR/platform/macos/AppIcon.icns" ]; then
            print_status "Generating macOS icon from molvis-icon.png..."
            "$SCRIPT_DIR/scripts/generate_macos_icon.sh"
        fi
    fi

    print_status "Creating build directory..."
    mkdir -p "$BUILD_DIR"
    cd "$BUILD_DIR"

    print_status "Running CMake..."
    cmake .. -DCMAKE_BUILD_TYPE=Release

    print_status "Building..."
    cmake --build . --config Release -j$(sysctl -n hw.ncpu)

    print_status "Build complete!"
    echo ""
    echo "Application built at: $BUILD_DIR/bin/Molecule Visualizer.app"
    echo "Run with: open $BUILD_DIR/bin/Molecule Visualizer.app"
    echo "Or: $BUILD_DIR/bin/Molecule Visualizer.app/Contents/MacOS/Molecule Visualizer"
}

# Run
run() {
    if [ ! -d "$BUILD_DIR/bin/Molecule Visualizer.app" ]; then
        build
    fi

    print_status "Running Molecule Visualizer..."
    "$BUILD_DIR/bin/Molecule Visualizer.app/Contents/MacOS/Molecule Visualizer"
}

# Main
case "${1:-}" in
    clean)
        clean
        ;;
    run)
        run
        ;;
    help|--help|-h)
        echo "MolVis macOS Build Script"
        echo ""
        echo "Usage: $0 [command]"
        echo ""
        echo "Commands:"
        echo "  (none)    Build the application"
        echo "  clean     Remove build artifacts"
        echo "  run       Build and run the application"
        echo "  help      Show this help message"
        ;;
    *)
        build
        ;;
esac
