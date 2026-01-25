#!/bin/bash
# MolVis - macOS Release Script
#
# Creates a local build, packages it, and uploads to GitHub Releases.
#
# Usage:
#   ./scripts/release/release_mac.sh                    Build, package, and release
#   ./scripts/release/release_mac.sh --version 0.3.2    Specify version explicitly
#   ./scripts/release/release_mac.sh --draft            Create as draft release
#   ./scripts/release/release_mac.sh --build-only       Build and package, don't upload
#   ./scripts/release/release_mac.sh --upload-only      Upload existing package to GitHub
#
# Requirements:
#   - GitHub CLI: brew install gh
#   - Authenticated: gh auth login

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$SCRIPT_DIR/../.."
BUILD_DIR="$PROJECT_ROOT/build_mac"
DIST_DIR="$PROJECT_ROOT/dist"
APP_NAME="MolVis"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
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

print_info() {
    echo -e "${BLUE}Info:${NC} $1"
}

# Default options
VERSION=""
DRAFT=false
BUILD_ONLY=false
UPLOAD_ONLY=false

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --version|-v)
            VERSION="$2"
            shift 2
            ;;
        --draft|-d)
            DRAFT=true
            shift
            ;;
        --build-only|-b)
            BUILD_ONLY=true
            shift
            ;;
        --upload-only|-u)
            UPLOAD_ONLY=true
            shift
            ;;
        --help|-h)
            echo "MolVis macOS Release Script"
            echo ""
            echo "Usage: $0 [options]"
            echo ""
            echo "Options:"
            echo "  --version, -v VERSION  Specify version (default: from RELEASE_NOTES.md)"
            echo "  --draft, -d            Create as draft release"
            echo "  --build-only, -b       Build and package only, don't upload"
            echo "  --upload-only, -u      Upload existing package to GitHub"
            echo "  --help, -h             Show this help message"
            echo ""
            echo "Examples:"
            echo "  $0                      # Full release with auto-detected version"
            echo "  $0 --version 0.3.2      # Release specific version"
            echo "  $0 --draft              # Create draft release for review"
            echo "  $0 --build-only         # Just build, upload later"
            exit 0
            ;;
        *)
            print_error "Unknown option: $1"
            exit 1
            ;;
    esac
done

# Get version from RELEASE_NOTES.md if not specified
get_version() {
    if [ -z "$VERSION" ]; then
        # Extract version from first ## v line in RELEASE_NOTES.md
        VERSION=$(grep -m1 "^## v" "$PROJECT_ROOT/RELEASE_NOTES.md" | sed 's/## v\([0-9.]*\).*/\1/')
        if [ -z "$VERSION" ]; then
            print_error "Could not detect version from RELEASE_NOTES.md"
            echo "Please specify version with --version"
            exit 1
        fi
    fi
    print_info "Version: $VERSION"
}

# Check for GitHub CLI
check_gh_cli() {
    if ! command -v gh &>/dev/null; then
        print_error "GitHub CLI (gh) not found"
        echo "Install with: brew install gh"
        exit 1
    fi

    # Check if authenticated
    if ! gh auth status &>/dev/null; then
        print_error "GitHub CLI not authenticated"
        echo "Run: gh auth login"
        exit 1
    fi

    print_status "GitHub CLI authenticated"
}

# Build the application
build_app() {
    print_status "Building MolVis..."
    "$PROJECT_ROOT/build_mac.sh"
}

# Package as DMG
package_app() {
        print_status "Packaging as DMG..."

        DMG_NAME="${APP_NAME}-${VERSION}-macOS"
        mkdir -p "$DIST_DIR"
        DMG_TEMP="$DIST_DIR/dmg_temp"
        rm -rf "$DMG_TEMP"
        mkdir -p "$DMG_TEMP"
        cp -R "$BUILD_DIR/bin/${APP_NAME}.app" "$DMG_TEMP/"
        ln -s /Applications "$DMG_TEMP/Applications"

        cat > "$DMG_TEMP/README.txt" << 'EOF'
MolVis - Molecule Visualizer
============================

Installation:
    1. Drag MolVis.app to the Applications folder
    2. If you see a warning, right-click and choose Open, or run:
         xattr -cr /Applications/MolVis.app
    3. Open MolVis from Applications

Requirements:
    - macOS 11.0 (Big Sur) or later
    - Any Mac with Metal support (Intel or Apple Silicon)

Usage:
    - Browse molecules by category in the sidebar
    - Drag to rotate, scroll to zoom
    - Adjust view settings in the View Controls panel

For more information, visit:
    https://github.com/AlfredBr/Molecule-Visualizer
EOF

        # Remove old DMG if exists
        rm -f "$DIST_DIR/${DMG_NAME}.dmg"

        # --- SIGNING AND NOTARIZATION ---
        # Use Developer ID Application if available
        if [ -n "$APPLE_IDENTITY" ]; then
                print_status "Signing app and SDL2 dylib with Developer ID..."
                codesign --force --deep --options runtime --sign "$APPLE_IDENTITY" "$DMG_TEMP/${APP_NAME}.app/Contents/Frameworks/libSDL2-2.0.0.dylib"
                codesign --force --deep --options runtime --sign "$APPLE_IDENTITY" "$DMG_TEMP/${APP_NAME}.app"
        fi

        # Create DMG
        hdiutil create -volname "$APP_NAME" \
                -srcfolder "$DMG_TEMP" \
                -ov -format UDZO \
                "$DIST_DIR/${DMG_NAME}.dmg"

        # Cleanup
        rm -rf "$DMG_TEMP"

        # Sign DMG if identity is set
        if [ -n "$APPLE_IDENTITY" ]; then
                print_status "Signing DMG with Developer ID..."
                codesign --force --sign "$APPLE_IDENTITY" "$DIST_DIR/${DMG_NAME}.dmg"
        fi

        # Notarize if credentials are set
        if [ -n "$APPLE_KEYCHAIN_PROFILE" ]; then
                print_status "Submitting DMG for notarization..."
                xcrun notarytool submit "$DIST_DIR/${DMG_NAME}.dmg" --keychain-profile "$APPLE_KEYCHAIN_PROFILE" --wait
                print_status "Stapling notarization ticket..."
                xcrun stapler staple "$DIST_DIR/${DMG_NAME}.dmg"
        fi

        # Calculate checksum
        cd "$DIST_DIR"
        shasum -a 256 "${DMG_NAME}.dmg" > "${DMG_NAME}.dmg.sha256"

        print_status "Package created: $DIST_DIR/${DMG_NAME}.dmg"
        echo "SHA256: $(cat ${DMG_NAME}.dmg.sha256)"
}

# Create GitHub release
create_release() {
    print_status "Creating GitHub release v${VERSION}..."

    DMG_NAME="${APP_NAME}-${VERSION}-macOS"
    TAG="v${VERSION}"

    # Check if tag exists
    if git rev-parse "$TAG" &>/dev/null; then
        print_warning "Tag $TAG already exists locally"
    fi

    # Check if release already exists
    if gh release view "$TAG" &>/dev/null; then
        print_warning "Release $TAG already exists on GitHub"
        echo ""
        read -p "Do you want to upload files to the existing release? (y/N) " -n 1 -r
        echo ""
        if [[ $REPLY =~ ^[Yy]$ ]]; then
            upload_to_existing "$TAG"
            return
        else
            print_error "Release already exists. Delete it first or use a different version."
            echo "  gh release delete $TAG --yes"
            exit 1
        fi
    fi

    # Extract release notes for this version from RELEASE_NOTES.md
    # Use sed instead of head -n -1 (which isn't supported on macOS)
    RELEASE_NOTES=$(awk "/^## v${VERSION}/,/^## v[0-9]/" "$PROJECT_ROOT/RELEASE_NOTES.md" | sed '$ d')

    # Create tag if it doesn't exist
    if ! git rev-parse "$TAG" &>/dev/null; then
        print_status "Creating git tag $TAG..."
        git tag -a "$TAG" -m "MolVis $TAG"
        git push origin "$TAG"
    fi

    # Build release command
    RELEASE_CMD="gh release create $TAG"
    RELEASE_CMD="$RELEASE_CMD --title \"MolVis $TAG\""

    if [ "$DRAFT" = true ]; then
        RELEASE_CMD="$RELEASE_CMD --draft"
    fi


    # Create release with files
    if [ "$DRAFT" = true ]; then
        gh release create "$TAG" \
            --title "MolVis $TAG" \
            --notes "$RELEASE_NOTES" \
            --draft \
            "$DIST_DIR/${DMG_NAME}.dmg" \
            "$DIST_DIR/${DMG_NAME}.dmg.sha256"
        print_status "Draft release created: https://github.com/AlfredBr/Molecule-Visualizer/releases/tag/$TAG"
        print_info "Edit and publish the draft on GitHub when ready"
    else
        gh release create "$TAG" \
            --title "MolVis $TAG" \
            --notes "$RELEASE_NOTES" \
            "$DIST_DIR/${DMG_NAME}.dmg" \
            "$DIST_DIR/${DMG_NAME}.dmg.sha256"
        print_status "Release published: https://github.com/AlfredBr/Molecule-Visualizer/releases/tag/$TAG"
    fi
}

# Upload to existing release
upload_to_existing() {
    local TAG="$1"
    DMG_NAME="${APP_NAME}-${VERSION}-macOS"

    print_status "Uploading files to existing release $TAG..."

    gh release upload "$TAG" \
        "$DIST_DIR/${DMG_NAME}.dmg" \
        "$DIST_DIR/${DMG_NAME}.dmg.sha256" \
        --clobber

    print_status "Files uploaded to: https://github.com/AlfredBr/Molecule-Visualizer/releases/tag/$TAG"
}

# Main
main() {
    echo ""
    echo "=========================================="
    echo "  MolVis macOS Release Script"
    echo "=========================================="
    echo ""

    get_version

    if [ "$UPLOAD_ONLY" = true ]; then
        check_gh_cli

        DMG_NAME="${APP_NAME}-${VERSION}-macOS"
        if [ ! -f "$DIST_DIR/${DMG_NAME}.dmg" ]; then
            print_error "Package not found: $DIST_DIR/${DMG_NAME}.dmg"
            echo "Run without --upload-only to build first"
            exit 1
        fi

        create_release
    elif [ "$BUILD_ONLY" = true ]; then
        build_app
        package_app

        echo ""
        echo "=========================================="
        echo "  Build complete!"
        echo "=========================================="
        echo ""
        echo "  Package: $DIST_DIR/${APP_NAME}-${VERSION}-macOS.dmg"
        echo ""
        echo "  To upload to GitHub, run:"
        echo "    $0 --upload-only --version $VERSION"
        echo ""
    else
        check_gh_cli
        build_app
        package_app
        create_release

        echo ""
        echo "=========================================="
        echo "  Release complete!"
        echo "=========================================="
        echo ""
    fi
}

main
