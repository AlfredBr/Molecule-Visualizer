#!/bin/bash
# Generate macOS .icns file from a PNG image
#
# Usage: ./scripts/generate_macos_icon.sh [input.png] [output.icns]
#
# If no arguments provided, uses molvis-icon.png and outputs to platform/macos/AppIcon.icns

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"

# Default input/output
INPUT_PNG="${1:-$PROJECT_DIR/molvis-icon.png}"
OUTPUT_ICNS="${2:-$PROJECT_DIR/platform/macos/AppIcon.icns}"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

print_status() {
    echo -e "${GREEN}==>${NC} $1"
}

print_error() {
    echo -e "${RED}Error:${NC} $1"
}

# Check if input file exists
if [ ! -f "$INPUT_PNG" ]; then
    print_error "Input file not found: $INPUT_PNG"
    echo "Please provide a PNG file (ideally 1024x1024 pixels)"
    exit 1
fi

# Check if sips is available (macOS built-in)
if ! command -v sips &>/dev/null; then
    print_error "sips command not found. This script requires macOS."
    exit 1
fi

print_status "Generating macOS icon from: $INPUT_PNG"

# Create temporary iconset directory
ICONSET_DIR="$PROJECT_DIR/AppIcon.iconset"
rm -rf "$ICONSET_DIR"
mkdir -p "$ICONSET_DIR"

# Generate all required icon sizes
# macOS requires specific sizes for the iconset
declare -a SIZES=(
    "16:icon_16x16.png"
    "32:icon_16x16@2x.png"
    "32:icon_32x32.png"
    "64:icon_32x32@2x.png"
    "128:icon_128x128.png"
    "256:icon_128x128@2x.png"
    "256:icon_256x256.png"
    "512:icon_256x256@2x.png"
    "512:icon_512x512.png"
    "1024:icon_512x512@2x.png"
)

for SIZE_NAME in "${SIZES[@]}"; do
    SIZE="${SIZE_NAME%%:*}"
    NAME="${SIZE_NAME##*:}"
    print_status "Generating ${SIZE}x${SIZE} -> $NAME"
    sips -z "$SIZE" "$SIZE" "$INPUT_PNG" --out "$ICONSET_DIR/$NAME" >/dev/null 2>&1
done

# Convert iconset to icns
print_status "Creating .icns file..."
iconutil -c icns "$ICONSET_DIR" -o "$OUTPUT_ICNS"

# Cleanup
rm -rf "$ICONSET_DIR"

print_status "Icon created successfully: $OUTPUT_ICNS"

# Show file info
ls -lh "$OUTPUT_ICNS"
