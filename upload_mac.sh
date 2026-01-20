#!/bin/bash
# upload_mac.sh - Upload a DMG to an existing GitHub release
# Usage: ./upload_mac.sh <version> <dmg_path>
# Example: ./upload_mac.sh 0.3.6 dist/MolVis-0.3.6-macOS.dmg

set -e

if [ $# -lt 2 ]; then
    echo "Usage: $0 <version> <dmg_path>"
    echo "Example: $0 0.3.6 dist/MolVis-0.3.6-macOS.dmg"
    exit 1
fi

VERSION="$1"
DMG_PATH="$2"
SHA256_PATH="${DMG_PATH}.sha256"
TAG="v${VERSION}"

if ! command -v gh &>/dev/null; then
    echo "Error: GitHub CLI (gh) not found. Install with: brew install gh"
    exit 1
fi

if [ ! -f "$DMG_PATH" ]; then
    echo "Error: DMG not found: $DMG_PATH"
    exit 1
fi

if [ ! -f "$SHA256_PATH" ]; then
    echo "Calculating SHA256..."
    shasum -a 256 "$DMG_PATH" > "$SHA256_PATH"
fi

echo "Uploading $DMG_PATH and $SHA256_PATH to release $TAG..."

gh release upload "$TAG" "$DMG_PATH" "$SHA256_PATH" --clobber

echo "Done. See: https://github.com/AlfredBr/Molecule-Visualizer/releases/tag/$TAG"
