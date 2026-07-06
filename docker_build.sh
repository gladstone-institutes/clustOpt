#!/usr/bin/env bash
set -euo pipefail

IMAGE="natalie23gill/clustopt"
PLATFORM="linux/amd64"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Read version from DESCRIPTION (e.g. "Version: 1.2.3")
VERSION=$(grep -m1 '^Version:' "$SCRIPT_DIR/DESCRIPTION" | sed 's/^Version:[[:space:]]*//')
if [ -z "$VERSION" ]; then
    echo "Error: could not read version from DESCRIPTION" >&2
    exit 1
fi

PUSH=false
NO_CACHE=""
PLAIN=""
for arg in "$@"; do
    case "$arg" in
        --push) PUSH=true ;;
        --no-cache) NO_CACHE="--no-cache" ;;
        --plain) PLAIN="--progress=plain" ;;
        *) echo "Usage: $0 [--push] [--no-cache] [--plain]" >&2; exit 1 ;;
    esac
done

TAG="${IMAGE}:${VERSION}"

# Warn if the tag already exists locally
if docker image inspect "$TAG" &>/dev/null; then
    echo "Warning: ${TAG} already exists locally and will be overwritten."
    read -rp "Continue? [y/N] " confirm
    if [[ ! "$confirm" =~ ^[Yy]$ ]]; then
        echo "Aborted." >&2
        exit 1
    fi
fi

echo "Building ${TAG} for ${PLATFORM}"

if [ "$PUSH" = true ]; then
    docker buildx build \
        --platform "$PLATFORM" \
        --tag "$TAG" \
        --tag "${IMAGE}:latest" \
        $NO_CACHE $PLAIN \
        --push \
        "$SCRIPT_DIR"
else
    docker buildx build \
        --platform "$PLATFORM" \
        --tag "$TAG" \
        --tag "${IMAGE}:latest" \
        $NO_CACHE $PLAIN \
        --load \
        "$SCRIPT_DIR"
fi

echo "Done: ${TAG}"
