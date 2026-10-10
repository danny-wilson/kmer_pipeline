#!/usr/bin/env bash
# Build a staging directory from a kmer_pipeline commit. Usage: stage.sh COMMIT
# Output: $KMER_E2E_ROOT/staging/SHA/ containing the commit's pipeline scripts
# (mode 755), report css/js, kmer_pipeline.nf, symlinks to the image's C++
# tools, a software-file template and a STAGED_SHA stamp. Bind it into the
# container at the same path.
set -euo pipefail
{ # read whole script before running, so edits can't affect a running copy

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
# shellcheck source=config.sh
. "$HERE/config.sh"
REPO=${KMER_E2E_REPO:-$(cd "$HERE/../.." && git rev-parse --show-toplevel)}
CPP_TOOLS="kmerlist2pattern pattern2kinship patterncounts patternmerge sort_strings stringlist2count stringlist2pattern"

SHA=$(git -C "$REPO" rev-parse --verify "${1:?usage: stage.sh COMMIT}^{commit}")
DEST=$KMER_E2E_ROOT/staging/$SHA
[ -e "$DEST" ] && { echo "$DEST already exists" >&2; exit 1; }

# Drift warning (not a hard stop): the toolchain image was built from a
# different commit's C++/Dockerfile/dependencies than the one being staged.
if ! git -C "$REPO" diff --quiet "$KMER_E2E_IMAGE_COMMIT" "$SHA" -- C++ Makefile Dockerfile 2>/dev/null; then
  echo "warning: C++/Makefile/Dockerfile differ between the image's commit ($KMER_E2E_IMAGE_COMMIT)" >&2
  echo "         and $SHA -- the toolchain in the image may not match what is being staged" >&2
fi

mkdir -p "$KMER_E2E_ROOT/staging"
TMP=$(mktemp -d "$KMER_E2E_ROOT/staging/.tmp.XXXXXX")
trap 'rm -rf "$TMP"' EXIT

git -C "$REPO" archive "$SHA" | tar -x -C "$TMP"
mkdir -p "$DEST"
# Scripts are installed flat, as the Dockerfile installs them in /usr/local/bin
shopt -s nullglob
for f in "$TMP"/*.R "$TMP"/*.Rscript "$TMP"/python/*.py; do install -m 755 "$f" "$DEST/"; done
install -m 644 "$TMP/report.css" "$TMP/report.js" "$TMP/kmer_pipeline.nf" "$DEST/"
install -m 644 "$TMP/example/pipeline_software_location.txt" "$DEST/"
for t in $CPP_TOOLS; do ln -s "/usr/local/bin/$t" "$DEST/$t"; done
echo "$SHA" > "$DEST/STAGED_SHA"
chmod -R a-w "$DEST"
echo "$DEST"
exit
}
