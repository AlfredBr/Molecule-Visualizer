#!/bin/bash
# Bundle SDL2 dylib into the app bundle for distribution
#
# This script:
# 1. Finds the SDL2 dylib (Homebrew location)
# 2. Copies it to the app's Frameworks folder
# 3. Updates the executable to look for SDL2 in @executable_path/../Frameworks
#
# Usage: ./bundle_sdl2.sh /path/to/MolVis.app

set -e

APP_BUNDLE="$1"

if [ -z "$APP_BUNDLE" ]; then
    echo "Usage: $0 /path/to/MolVis.app"
    exit 1
fi

FRAMEWORKS_DIR="$APP_BUNDLE/Contents/Frameworks"
EXECUTABLE="$APP_BUNDLE/Contents/MacOS/MolVis"

# Create Frameworks directory if it doesn't exist
mkdir -p "$FRAMEWORKS_DIR"

# Find SDL2 dylib
SDL2_DYLIB=""

# Check common Homebrew locations
if [ -f "/opt/homebrew/lib/libSDL2-2.0.0.dylib" ]; then
    SDL2_DYLIB="/opt/homebrew/lib/libSDL2-2.0.0.dylib"
elif [ -f "/usr/local/lib/libSDL2-2.0.0.dylib" ]; then
    SDL2_DYLIB="/usr/local/lib/libSDL2-2.0.0.dylib"
else
    # Try to find it via otool
    SDL2_PATH=$(otool -L "$EXECUTABLE" | grep -o '/.*libSDL2.*\.dylib' | head -1)
    if [ -f "$SDL2_PATH" ]; then
        SDL2_DYLIB="$SDL2_PATH"
    fi
fi

if [ -z "$SDL2_DYLIB" ] || [ ! -f "$SDL2_DYLIB" ]; then
    echo "Error: Could not find libSDL2-2.0.0.dylib"
    echo "Please install SDL2: brew install sdl2"
    exit 1
fi

echo "Found SDL2 at: $SDL2_DYLIB"

# Copy SDL2 to Frameworks
cp "$SDL2_DYLIB" "$FRAMEWORKS_DIR/"
chmod 755 "$FRAMEWORKS_DIR/libSDL2-2.0.0.dylib"

# Get the original install name from the dylib
ORIGINAL_ID=$(otool -D "$SDL2_DYLIB" | tail -1)

echo "Original SDL2 install name: $ORIGINAL_ID"

# Remove signature from the copied dylib before modifying
codesign --remove-signature "$FRAMEWORKS_DIR/libSDL2-2.0.0.dylib" 2>/dev/null || true

# Change the install name in the copied dylib
install_name_tool -id "@executable_path/../Frameworks/libSDL2-2.0.0.dylib" \
    "$FRAMEWORKS_DIR/libSDL2-2.0.0.dylib"

# Update the executable to look for SDL2 in Frameworks
# First, get the current reference path from the executable
SDL2_REF=$(otool -L "$EXECUTABLE" | grep SDL2 | awk '{print $1}')

if [ -n "$SDL2_REF" ]; then
    echo "Updating executable reference from: $SDL2_REF"
    
    # Remove code signature first (install_name_tool can't modify signed binaries)
    codesign --remove-signature "$EXECUTABLE" 2>/dev/null || true
    
    install_name_tool -change "$SDL2_REF" \
        "@executable_path/../Frameworks/libSDL2-2.0.0.dylib" \
        "$EXECUTABLE"
fi

# Re-sign everything with ad-hoc signature (order matters: frameworks first, then executable)
echo "Re-signing SDL2 dylib..."
codesign --force --sign - "$FRAMEWORKS_DIR/libSDL2-2.0.0.dylib"

echo "Re-signing executable..."
codesign --force --sign - "$EXECUTABLE"

# Verify the change
echo ""
echo "Verifying SDL2 linkage:"
otool -L "$EXECUTABLE" | grep -i sdl2 || echo "  (no SDL2 reference found - may be statically linked)"

echo ""
echo "SDL2 bundled successfully!"
