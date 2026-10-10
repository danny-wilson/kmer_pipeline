#!/usr/bin/env bash
# Build a test SIF from kmer_pipeline's Dockerfile at COMMIT: the Dockerfile is
# translated by dockerfile2def.py and built with fakeroot from a git archive of
# COMMIT, as `docker build` would from a clean checkout. No Docker is used.
# Usage: build_test_sif.sh COMMIT [VERSION]
# Output: $KMER_E2E_ROOT/images/test-<sha7>.sif (and the .def next to it)
set -euo pipefail
{
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
# shellcheck source=config.sh
. "$HERE/config.sh"
REPO=${KMER_E2E_REPO:-$(cd "$HERE/../.." && git rev-parse --show-toplevel)}

SHA=$(git -C "$REPO" rev-parse --verify "${1:?usage: build_test_sif.sh COMMIT [VERSION]}^{commit}")
VERSION=${2:-test-${SHA:0:7}}
OUT=$KMER_E2E_ROOT/images/test-${SHA:0:7}.sif
[ -e "$OUT" ] && { echo "$OUT already exists" >&2; exit 1; }
mkdir -p "$KMER_E2E_ROOT/images"
TMP=$(mktemp -d "$KMER_E2E_ROOT/images/.tmp.XXXXXX")
trap 'rm -rf "$TMP"' EXIT

mkdir "$TMP/context"
git -C "$REPO" archive "$SHA" | tar -x -C "$TMP/context"
# .dockerignore: .git (git archive leaves it out anyway), docs, pycache/
# pytest-cache, and the e2e local config/SIF patterns. If this changes,
# update this check deliberately rather than silently building a different
# context from what a real `docker build` would see.
WANT=$'.git\ndocs\n**/__pycache__\n**/.pytest_cache\ntests/e2e/local.conf\n*.sif'
[ "$(cat "$TMP/context/.dockerignore")" = "$WANT" ] || { echo ".dockerignore changed: update this script" >&2; exit 1; }
rm -rf "$TMP/context/docs"
python3 "$HERE/dockerfile2def.py" "$TMP/context/Dockerfile" "$TMP/context" "$VERSION" > "${OUT%.sif}.def"
# Unpack on local disk: the rootfs is ~6 GB of small files, far too slow on a
# shared/networked filesystem
BUILD_TMP=$(mktemp -d "${KMER_E2E_LOCAL_TMP:-/tmp}/kmer_sif.XXXXXX")
trap 'rm -rf "$TMP" "$BUILD_TMP"' EXIT
export APPTAINER_TMPDIR=$BUILD_TMP APPTAINER_CACHEDIR=$KMER_E2E_ROOT/images/cache
apptainer build --fakeroot --ignore-fakeroot-command "$TMP/image.sif" "${OUT%.sif}.def"
mv "$TMP/image.sif" "$OUT"
echo "$OUT"
exit
}
