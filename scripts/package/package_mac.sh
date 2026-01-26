#!/bin/bash
# MolVis - macOS Packaging Script
#
# Creates a distributable DMG file for easy installation.
#
# Usage:
#   ./scripts/package/package_mac.sh              Create DMG package
#   ./scripts/package/package_mac.sh --notarize   Create and notarize (requires Apple Developer account)
#
# The resulting DMG will be in the dist/ directory.

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$SCRIPT_DIR/../.."
BUILD_DIR="$PROJECT_ROOT/build_mac"
DIST_DIR="$PROJECT_ROOT/dist"
APP_NAME="MolVis"
VERSION="0.4.1"
DMG_NAME="${APP_NAME}-${VERSION}-macOS"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

print_status() {
    echo -e "${GREEN}==>${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}Warning:${NC} $1"
}

print_error() {
    echo -e "${RED}Error:${NC} $1"
}

# Check if app is built
check_app() {
    if [ ! -d "$BUILD_DIR/bin/${APP_NAME}.app" ]; then
        print_status "App not found. Building first..."
        "$PROJECT_ROOT/build_mac.sh"
    fi
}

# Generate icon if needed
generate_icon() {
    if [ -f "$PROJECT_ROOT/molvis-icon.png" ]; then
        if [ ! -f "$PROJECT_ROOT/platform/macos/AppIcon.icns" ] || \
           [ "$PROJECT_ROOT/molvis-icon.png" -nt "$PROJECT_ROOT/platform/macos/AppIcon.icns" ]; then
            print_status "Generating macOS icon..."
            "$PROJECT_ROOT/scripts/generate_macos_icon.sh"
        fi
    fi
}

# Create DMG
create_dmg() {
    print_status "Creating DMG package..."

    # Create dist directory
    mkdir -p "$DIST_DIR"

    # Create temporary directory for DMG contents
    DMG_TEMP="$DIST_DIR/dmg_temp"
    rm -rf "$DMG_TEMP"
    mkdir -p "$DMG_TEMP"

    # Copy app to temp directory
    print_status "Copying application..."
    cp -R "$BUILD_DIR/bin/${APP_NAME}.app" "$DMG_TEMP/"

    # Create symlink to Applications folder
    ln -s /Applications "$DMG_TEMP/Applications"

    # Remove old DMG if exists
    rm -f "$DIST_DIR/${DMG_NAME}.dmg"

    # Create DMG
    print_status "Building DMG file..."

    # Check if create-dmg is available (prettier DMGs)
    if command -v create-dmg &>/dev/null; then
        create-dmg \
            --volname "$APP_NAME" \
            --volicon "$BUILD_DIR/bin/${APP_NAME}.app/Contents/Resources/AppIcon.icns" \
            --window-pos 200 120 \
            --window-size 600 400 \
            --icon-size 100 \
            --icon "${APP_NAME}.app" 150 190 \
            --hide-extension "${APP_NAME}.app" \
            --app-drop-link 450 185 \
            "$DIST_DIR/${DMG_NAME}.dmg" \
            "$DMG_TEMP" 2>/dev/null || {
                # Fallback to hdiutil if create-dmg fails
                print_warning "create-dmg failed, using hdiutil..."
                hdiutil create -volname "$APP_NAME" \
                    -srcfolder "$DMG_TEMP" \
                    -ov -format UDZO \
                    "$DIST_DIR/${DMG_NAME}.dmg"
            }
    else
        # Use built-in hdiutil
        hdiutil create -volname "$APP_NAME" \
            -srcfolder "$DMG_TEMP" \
            -ov -format UDZO \
            "$DIST_DIR/${DMG_NAME}.dmg"
    fi

    # Cleanup
    rm -rf "$DMG_TEMP"

    print_status "DMG created: $DIST_DIR/${DMG_NAME}.dmg"
}

# Calculate checksum
create_checksum() {
    print_status "Calculating SHA256 checksum..."

    cd "$DIST_DIR"
    shasum -a 256 "${DMG_NAME}.dmg" > "${DMG_NAME}.dmg.sha256"

    echo ""
    echo "Checksum:"
    cat "${DMG_NAME}.dmg.sha256"
    echo ""
}

# Get DMG size
show_info() {
    DMG_PATH="$DIST_DIR/${DMG_NAME}.dmg"
    DMG_SIZE=$(ls -lh "$DMG_PATH" | awk '{print $5}')

    echo ""
    echo "================================================"
    echo "  Package created successfully!"
    echo "================================================"
    echo ""
    echo "  File: ${DMG_NAME}.dmg"
    echo "  Size: ${DMG_SIZE}"
    echo "  Path: $DMG_PATH"
    echo ""
    echo "  To install:"
    echo "    1. Open the DMG file"
    echo "    2. Drag MolVis to Applications"
    echo "    3. Eject the DMG"
    echo ""
    echo "  To distribute:"
    echo "    Upload ${DMG_NAME}.dmg to GitHub Releases"
    echo "    Include the .sha256 checksum file"
    echo ""
}

# Optional: Code sign the app (requires Apple Developer account)
codesign_app() {
    if [ -n "$APPLE_IDENTITY" ]; then
        print_status "Code signing application..."
        codesign --force --deep --sign "$APPLE_IDENTITY" \
            "$BUILD_DIR/bin/${APP_NAME}.app"
        print_status "Code signing complete"
    else
        print_warning "APPLE_IDENTITY not set - skipping code signing"
        echo "  To code sign, set: export APPLE_IDENTITY='Developer ID Application: Your Name'"
    fi
}

# Optional: Notarize the DMG (requires Apple Developer account)
notarize_dmg() {
    if [ -z "$APPLE_ID" ] || [ -z "$APPLE_TEAM_ID" ]; then
        print_error "Notarization requires APPLE_ID and APPLE_TEAM_ID environment variables"
        echo "  export APPLE_ID='your@email.com'"
        echo "  export APPLE_TEAM_ID='XXXXXXXXXX'"
        exit 1
    fi

    print_status "Submitting for notarization..."
    xcrun notarytool submit "$DIST_DIR/${DMG_NAME}.dmg" \
        --apple-id "$APPLE_ID" \
        --team-id "$APPLE_TEAM_ID" \
        --wait

    print_status "Stapling notarization ticket..."
    xcrun stapler staple "$DIST_DIR/${DMG_NAME}.dmg"

    print_status "Notarization complete"
}

# Main
main() {
    echo ""
    echo "MolVis macOS Packager"
    echo "====================="
    echo ""

    generate_icon
    check_app

    # Optional code signing
    if [ -n "$APPLE_IDENTITY" ]; then
        codesign_app
    fi

    create_dmg
    create_checksum

    # Optional notarization
    if [ "$1" == "--notarize" ]; then
        notarize_dmg
    fi

    show_info
}

main "$@"
