#!/usr/bin/env bash
# Build a code-stripped derivative of an existing SIF: everything the old code-full image
# installed into /usr/local/bin (the pipeline's own scripts) and /usr/share/kmer_pipeline is
# deleted, so a test running against this derivative cannot silently fall back to the image's own
# stale copy of a script instead of this checkout's bin/ (guards P2) - the exact gap that let B4's
# work (countkmers.py and friends) pass tests that were actually still running the old image's
# version. PLAN.md Part B, B9a: built before B8, since B8's own tests/ fix needs this derivative,
# not the plain base SIF, to actually prove anything about bin/ resolution.
# Usage: build_codefree_derivative.sh BASE_SIF [OUTPUT_NAME]
# Output: $KMER_E2E_ROOT/images/OUTPUT_NAME.sif (and the .def next to it)
set -euo pipefail
{
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
# shellcheck source=config.sh
. "$HERE/config.sh"

BASE_SIF=${1:?usage: build_codefree_derivative.sh BASE_SIF [OUTPUT_NAME]}
[ -f "$BASE_SIF" ] || { echo "BASE_SIF '$BASE_SIF' does not exist" >&2; exit 1; }
BASE_SIF=$(cd "$(dirname "$BASE_SIF")" && pwd)/$(basename "$BASE_SIF")
NAME=${2:-$(basename "$BASE_SIF" .sif)-codefree}
OUT=$KMER_E2E_ROOT/images/$NAME.sif
[ -e "$OUT" ] && { echo "$OUT already exists" >&2; exit 1; }
mkdir -p "$KMER_E2E_ROOT/images"
TMP=$(mktemp -d "$KMER_E2E_ROOT/images/.tmp.XXXXXX")
trap 'rm -rf "$TMP"' EXIT

cat > "$TMP/image.def" <<EOF
Bootstrap: localimage
From: $BASE_SIF

%post
	rm -f /usr/local/bin/*.py /usr/local/bin/*.R /usr/local/bin/report.js /usr/local/bin/report.css
	rm -rf /usr/share/kmer_pipeline
EOF
cp "$TMP/image.def" "${OUT%.sif}.def"

# Unpack/repack on local disk, not a shared/networked filesystem's own /tmp: the rootfs is
# several GB of small files, far too slow there (matching build_test_sif.sh's own pattern).
BUILD_TMP=$(mktemp -d "${KMER_E2E_LOCAL_TMP:-/tmp}/kmer_codefree.XXXXXX")
trap 'rm -rf "$TMP" "$BUILD_TMP"' EXIT
export APPTAINER_TMPDIR=$BUILD_TMP APPTAINER_CACHEDIR=$KMER_E2E_ROOT/images/cache
apptainer build --fakeroot --ignore-fakeroot-command "$TMP/image.sif" "$TMP/image.def"
mv "$TMP/image.sif" "$OUT"
echo "$OUT"
exit
}
