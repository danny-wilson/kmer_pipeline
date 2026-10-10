#!/usr/bin/env bash
# Check a release image against a kmer_pipeline commit.
# Usage: check_image.sh IMAGE.sif COMMIT [EXPECTED_SHA256]
#   IMAGE.sif: the test SIF (build_test_sif.sh) or a pulled release tag
#   EXPECTED_SHA256: optional; if given, IMAGE.sif's own sha256 must match it
#     (the local-file hash, checked against the digest recorded when the image
#     was pulled -- see KMER_E2E_IMAGE_COMMIT and the project's own provenance
#     record -- not re-derived here by re-pulling the registry).
# Checks: every installed pipeline file equals the commit's, and no other
# pipeline files remain; /usr/share/kmer_pipeline equals the commit's tree
# less the installed files; the environment (PYTHONNOUSERSITE, LC_ALL);
# Biopython and pytest import; every script runs --help.
# If $KMER_E2E_BASE_SIF is set, also checks: pip freeze and conda list
# unchanged apart from the added packages; R's installed packages unchanged;
# the C++ tools, external tools and Nextflow unchanged -- against that base
# image. This comparison is skipped (not failed) if the key is unset, since a
# "previous production image" baseline is a locally-kept asset, not something
# every clone of this harness has.
#
# NOTE: pip/conda/R package-version lists and tool lists here (and the
# $ADDED list below) are a snapshot taken when this script was ported; they
# must be maintained by hand, or deliberately regenerated, whenever the image
# changes -- this script cannot discover that on its own.
set -euo pipefail
{
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
# shellcheck source=config.sh
. "$HERE/config.sh"
REPO=${KMER_E2E_REPO:-$(cd "$HERE/../.." && git rev-parse --show-toplevel)}
ADDED="biopython==1.83 exceptiongroup==1.2.2 iniconfig==2.0.0 pluggy==1.5.0 pytest==7.4.4 tomli==2.0.1"

IMAGE=${1:?usage: check_image.sh IMAGE.sif COMMIT [EXPECTED_SHA256]}
SHA=$(git -C "$REPO" rev-parse --verify "${2:?usage: check_image.sh IMAGE.sif COMMIT [EXPECTED_SHA256]}^{commit}")
EXPECTED_SHA256=${3:-}
TMP=$(mktemp -d "$KMER_E2E_ROOT/runs/.check_image.XXXXXX")
trap 'rm -rf "$TMP"' EXIT
problems=0
bad() { echo "FAIL: $*"; problems=$((problems + 1)); }
ok() { echo "ok: $*"; }
run() { apptainer exec --containall --cleanenv "$@"; }

if [ -n "$EXPECTED_SHA256" ]; then
  got=$(sha256sum "$IMAGE" | cut -d' ' -f1)
  [ "$got" = "$EXPECTED_SHA256" ] && ok "sha256 matches the recorded digest ($got)" \
    || bad "sha256 $got does not match the recorded digest $EXPECTED_SHA256"
else
  echo "note: no EXPECTED_SHA256 given -- skipping the digest check"
fi

inventory() {  # environment fingerprint, one file per item
	local sif=$1 out=$2
	mkdir -p "$out"
	run "$sif" pip freeze --all 2>/dev/null | sort > "$out/pip"
	run "$sif" conda list 2>/dev/null | grep -v '^#' | awk '{print $1, $2, $4}' | sort > "$out/conda"
	run "$sif" Rscript -e 'ip = installed.packages(); writeLines(sort(paste(ip[, "Package"], ip[, "Version"])))' > "$out/R"
	run "$sif" bash -c 'cd /usr/local/bin && md5sum dsk dsk2ascii gemma nextflow kmerlist2pattern pattern2kinship patterncounts patternmerge sort_strings stringlist2count stringlist2pattern; md5sum /usr/bin/bowtie2-align-s /usr/bin/nucmer /usr/bin/blastn /opt/conda/bin/samtools' > "$out/tools"
}

if [ -n "${KMER_E2E_BASE_SIF:-}" ]; then
  inventory "$KMER_E2E_BASE_SIF" "$TMP/base"
  inventory "$IMAGE" "$TMP/new"
  # pip: the new list is the base list plus exactly the added packages
  printf '%s\n' $ADDED | cat - "$TMP/base/pip" | sort > "$TMP/pip.want"
  if diff "$TMP/pip.want" "$TMP/new/pip" > "$TMP/pip.diff"; then ok "pip freeze = base + $ADDED"; else bad "pip freeze"; cat "$TMP/pip.diff"; fi
  # conda list shows pip-installed packages too (channel pypi): allow only those
  diff "$TMP/base/conda" "$TMP/new/conda" | grep '^[<>]' > "$TMP/conda.diff" || true
  if ! grep -q '^<' "$TMP/conda.diff" && ! grep '^>' "$TMP/conda.diff" | grep -v -E ' pypi$' -q; then
	  ok "conda list = base + $(grep -c '^>' "$TMP/conda.diff") pypi packages"
  else bad "conda list"; cat "$TMP/conda.diff"; fi
  if diff "$TMP/base/R" "$TMP/new/R" > /dev/null; then ok "R packages unchanged ($(wc -l < "$TMP/new/R"))"; else bad "R packages"; diff "$TMP/base/R" "$TMP/new/R"; fi
  if diff "$TMP/base/tools" "$TMP/new/tools" > /dev/null; then ok "C++ tools, external tools and Nextflow unchanged"; else bad "tools"; diff "$TMP/base/tools" "$TMP/new/tools"; fi
else
  echo "note: KMER_E2E_BASE_SIF not set -- skipping the pip/conda/R/tools comparison against a baseline image"
fi

# Installed pipeline files equal the commit's
mkdir "$TMP/src"
git -C "$REPO" archive "$SHA" | tar -x -C "$TMP/src"
( shopt -s nullglob; cd "$TMP/src" && md5sum *.R *.Rscript kmer_pipeline.nf report.js report.css && cd python && md5sum *.py ) | sort -k2 > "$TMP/want.md5"
run "$IMAGE" bash -c 'shopt -s nullglob; cd /usr/local/bin && md5sum *.R *.Rscript *.py kmer_pipeline.nf report.js report.css' | sort -k2 > "$TMP/got.md5"
if diff "$TMP/want.md5" "$TMP/got.md5" > "$TMP/md5.diff"; then ok "installed scripts = $SHA ($(wc -l < "$TMP/got.md5") files)"; else bad "installed scripts"; cat "$TMP/md5.diff"; fi
run "$IMAGE" bash -c 'cd /usr/local/bin && for f in *; do [ -x "$f" ] || echo "not executable: $f"; done' > "$TMP/exec"
[ -s "$TMP/exec" ] && { bad "permissions"; cat "$TMP/exec"; } || ok "all of /usr/local/bin executable"

# /usr/share/kmer_pipeline = commit tree less the installed files and docs/ (.dockerignore)
( shopt -s nullglob; cd "$TMP/src" && rm -rf *.R *.Rscript python kmer_pipeline.nf report.js report.css docs && find . -type f | sort | xargs md5sum ) > "$TMP/share.want"
run "$IMAGE" bash -c 'cd /usr/share/kmer_pipeline && find . -type f | sort | xargs md5sum' > "$TMP/share.got"
if diff "$TMP/share.want" "$TMP/share.got" > "$TMP/share.diff"; then ok "/usr/share/kmer_pipeline = $SHA tree less installed files"; else bad "/usr/share/kmer_pipeline"; head -20 "$TMP/share.diff"; fi

# Environment and imports (as Nextflow runs it: --containall --cleanenv)
env_out=$(run "$IMAGE" bash -c 'echo "PYTHONNOUSERSITE=$PYTHONNOUSERSITE LC_ALL=$LC_ALL"')
[ "$env_out" = "PYTHONNOUSERSITE=1 LC_ALL=en_US.UTF-8" ] && ok "environment: $env_out" || bad "environment: $env_out"
imp=$(run "$IMAGE" python3 -c 'import Bio, pytest, numpy, pandas, scipy, jinja2; print(Bio.__version__, pytest.__version__, numpy.__version__, pandas.__version__, scipy.__version__, jinja2.__version__)')
[ "$imp" = "1.83 7.4.4 1.21.6 1.4.2 1.8.1 3.1.2" ] && ok "imports: $imp" || bad "imports: $imp"
before=$problems
for name in $(git -C "$REPO" ls-tree --name-only "$SHA" python/ | xargs -n1 basename); do
	case $name in rcompat.py|sequence_functions.py|Manhattan_functions.py|alignmentfunctions.py|inventory.py|reference.py|report_assets.py) continue;; esac
	run "$IMAGE" "/usr/local/bin/$name" --help > /dev/null 2>&1 || bad "$name --help"
done
[ "$problems" -eq "$before" ] && ok "every script runs --help"

echo "$problems problems"
[ "$problems" -eq 0 ]
exit
}
