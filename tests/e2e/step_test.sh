#!/usr/bin/env bash
# Step test: run one Nextflow stage with a staged commit, seeded with a
# reference run's snapshot of the previous stage, and compare the result with
# the reference's snapshot of that stage. The run's other settings (k-mer type
# and length, maxp, covariates, ID file, copied inputs, extra parameters) are
# read from the reference's metadata.json. The container mounts base_dir at
# the same path for every run, so absolute paths in outputs agree.
#
# Usage: step_test.sh STAGING_DIR REF_RUN_DIR STAGE [PYTHONHASHSEED]
#   e.g. step_test.sh $KMER_E2E_ROOT/staging/<sha> $KMER_E2E_ROOT/goldens/<label>/maxp2/protein11 1 0
# Output: $KMER_E2E_ROOT/goldens/step<STAGE>h<SEED>-<sha7>/maxp<M>/<run>/stage<STAGE>
# Exit status 0 only if the comparison passes.
set -euo pipefail
{ # read whole script before running, so edits can't affect a running copy
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
# shellcheck source=config.sh
. "$HERE/config.sh"
STAGING=${1:?usage}; REF=${2:?usage}; STAGE=${3:?usage}; HSEED=${4:-0}
SHA=$(cat "$STAGING/STAGED_SHA")

ARGS=()
while IFS= read -r line; do ARGS+=("$line"); done < <(python3 - "$REF/metadata.json" "$STAGE" <<'EOF'
import json, sys
d = json.load(open(sys.argv[1])); k = int(sys.argv[2])
out = ["--run", d["run"], "--kmer-type", d["kmer_type"], "--kmer-length", str(d["kmer_length"]),
       "--maxp", str(d["maxp"])]
if d.get("covariates"): out += ["--covariates", d["covariates"]]
if d.get("id_file"): out += ["--id-file", d["id_file"]]
for c in (d.get("copies") or "").split():
    out += ["--copy", c]
for p in (d.get("extra_params") or "").split("\n"):
    if p.strip(): out += ["--param", p.strip()]
if k > 1: out += ["--seed", sys.argv[1].rsplit("/", 1)[0] + f"/stage{k - 1}"]
print("\n".join(out))
EOF
)
MAXP=$(python3 -c "import json; print(json.load(open('$REF/metadata.json'))['maxp'])")
RUN=$(python3 -c "import json; print(json.load(open('$REF/metadata.json'))['run'])")
LABEL=step${STAGE}h$HSEED
OUT=$KMER_E2E_ROOT/goldens/$LABEL-${SHA:0:7}/maxp$MAXP/$RUN
mkdir -p "$(dirname "$OUT")"
"$HERE/make_reference.sh" --label "$LABEL" --sha "$SHA" --staging "$STAGING" --stages "$STAGE" \
  --container-env "PYTHONHASHSEED=$HSEED" "${ARGS[@]}" > "$OUT.make.log" 2>&1 \
  || echo "make_reference failed: $OUT.make.log" >&2
python3 "$HERE/compare.py" "$REF/stage$STAGE" "$OUT/stage$STAGE"
exit
}
