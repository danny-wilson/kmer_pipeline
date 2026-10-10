#!/usr/bin/env bash
# Full-run before/after check between any two commits (renamed from the private
# harness's check_fix.sh, which only ever compared against a hardcoded "r-fixes"
# label). Stages AFTER_SHA (if not already staged), generates nucleotide31 and
# protein11 references at maxp 2 under AFTER_LABEL, and compares every stage
# with the BEFORE references. Prints the differences; intended ones are named
# in the caller's own check, not decided here.
#
# Usage: compare_commits.sh AFTER_LABEL BEFORE_LABEL BEFORE_SHA AFTER_SHA
set -uo pipefail
{ # read whole script before running, so edits can't affect a running copy
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
# shellcheck source=config.sh
. "$HERE/config.sh"
REPO=${KMER_E2E_REPO:-$(cd "$HERE/../.." && git rev-parse --show-toplevel)}
AFTER_LABEL=${1:?usage: compare_commits.sh AFTER_LABEL BEFORE_LABEL BEFORE_SHA AFTER_SHA}
BEFORE_LABEL=${2:?}; BEFORE=${3:?}; AFTER=${4:?}
AFTER=$(git -C "$REPO" rev-parse "$AFTER"); BEFORE=$(git -C "$REPO" rev-parse "$BEFORE")
IMAGE_DIGEST=$(sha256sum "$KMER_E2E_SIF" | cut -d' ' -f1)

S=$KMER_E2E_ROOT/staging/$AFTER
[ -d "$S" ] || "$HERE/stage.sh" "$AFTER" > /dev/null

LOGS=$KMER_E2E_ROOT/runs/compare_commits-${AFTER:0:7}; mkdir -p "$LOGS"
pids=(); names=()
for spec in "nucleotide31 nucleotide 31" "protein11 protein 11"; do
  set -- $spec
  run=$1
  marker=$KMER_E2E_ROOT/goldens/$AFTER_LABEL-${AFTER:0:7}/maxp2/$run/.complete
  dir=$KMER_E2E_ROOT/goldens/$AFTER_LABEL-${AFTER:0:7}/maxp2/$run
  if [ -d "$dir" ]; then
    if [ -f "$marker" ] && [ "$(cat "$marker")" = "$IMAGE_DIGEST" ]; then
      echo "reusing existing $dir (image digest matches)"
      continue
    fi
    echo "$dir exists but has no completion marker matching the current image ($IMAGE_DIGEST) -- remove it to regenerate" >&2
    exit 1
  fi
  "$HERE/make_reference.sh" --label "$AFTER_LABEL" --sha "$AFTER" --run "$run" --kmer-type "$2" \
    --kmer-length "$3" --maxp 2 --staging "$S" > "$LOGS/$run.log" 2>&1 &
  pids+=("$!"); names+=("$run")
done
status=0
for i in "${!pids[@]}"; do
  if wait "${pids[$i]}"; then
    echo "$IMAGE_DIGEST" > "$KMER_E2E_ROOT/goldens/$AFTER_LABEL-${AFTER:0:7}/maxp2/${names[$i]}/.complete"
  else
    echo "make_reference.sh failed for ${names[$i]} (rc $?): see $LOGS/${names[$i]}.log" >&2
    status=1
  fi
done
[ "$status" -ne 0 ] && exit "$status"

for run in nucleotide31 protein11; do
  for k in 1 2 3 4 5 6 7; do
    out=$(python3 "$HERE/compare.py" "$KMER_E2E_ROOT/goldens/$BEFORE_LABEL-${BEFORE:0:7}/maxp2/$run/stage$k" \
      "$KMER_E2E_ROOT/goldens/$AFTER_LABEL-${AFTER:0:7}/maxp2/$run/stage$k") || status=1
    echo "== $run stage$k: $(echo "$out" | tail -n 1)"
    echo "$out" | sed '$d'
  done
done
exit $status
exit
}
