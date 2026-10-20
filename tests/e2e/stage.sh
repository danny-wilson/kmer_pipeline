#!/usr/bin/env bash
# Build a staging directory from a kmer_pipeline commit: a plain git archive,
# nothing flattened or renamed. Usage: stage.sh COMMIT
# Output: $KMER_E2E_ROOT/staging/SHA/, the commit's own tree (bin/,
# kmer_pipeline.nf, nextflow.config, conf/, example/, a STAGED_SHA stamp).
# Run it with "nextflow run $DEST/kmer_pipeline.nf ..." directly -- Nextflow
# puts $DEST/bin on PATH automatically (it's next to the main script), and the
# committed conf/software_location.txt's toolpath still finds the compiled
# C++ tools in the container image itself, unchanged by which commit is
# staged (B9: bin/ is now native, so there is nothing left to symlink or
# flatten here, unlike the pre-refactor layout this script used to mimic).
set -euo pipefail
{ # read whole script before running, so edits can't affect a running copy

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
# shellcheck source=config.sh
. "$HERE/config.sh"
REPO=${KMER_E2E_REPO:-$(cd "$HERE/../.." && git rev-parse --show-toplevel)}

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
mv "$TMP" "$DEST"
chmod -R u+w,go-w "$DEST"
echo "$SHA" > "$DEST/STAGED_SHA"
chmod -R a-w "$DEST"
echo "$DEST"
exit
}
