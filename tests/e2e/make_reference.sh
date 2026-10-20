#!/usr/bin/env bash
# Generate a reference run: the example pipeline run one Nextflow stage at a time
# (skip1..skip7), with a snapshot of analysis_dir after each stage.
#
# Usage: make_reference.sh --label LABEL --sha SHA --run NAME --kmer-type TYPE
#                          --kmer-length K --maxp M --staging DIR [--sif IMAGE]
#                          [--covariates FILE] [--stages "1 2 3 4 5 6 7"] [--param "name = value"]...
#                          [--id-file FILE (paths may use @BASE@ for the base directory)] [--copy FILE]...
#                          [--seed DIR]
#
# --staging is the commit's own tree (stage.sh's output): kmer_pipeline.nf is run
# directly from there, by path, so bin/ and nextflow.config -- and therefore the
# -profile apptainer/apptainer-vs-docker container setup, the identity base_dir/
# projectDir bind, and the pinned OPENBLAS_NUM_THREADS/OMP_NUM_THREADS env{} scope
# (B9) -- all come from the staged commit itself, not this script. --seed copies a
# stage snapshot (logs and nextflow.out excluded) into analysis_dir before the
# first stage, for step tests.
# Output: $KMER_E2E_ROOT/goldens/LABEL-SHA7/maxpM/NAME/{stage1..stage7,metadata.json}
set -euo pipefail
{ # read whole script before running, so edits can't affect a running copy

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
# shellcheck source=config.sh
. "$HERE/config.sh"
SIF=$KMER_E2E_SIF
NEXTFLOW=$KMER_E2E_NEXTFLOW
STAGING=""; COVARIATES=""; STAGES="1 2 3 4 5 6 7"; EXTRA_PARAMS=""; ID_FILE=""; COPIES=(); SEED=""; EXTRA_ENV=""
while [ $# -gt 0 ]; do
  case "$1" in
    --label) LABEL=$2; shift 2;;
    --sha) SHA=$2; shift 2;;
    --run) RUN=$2; shift 2;;
    --kmer-type) KTYPE=$2; shift 2;;
    --kmer-length) KLEN=$2; shift 2;;
    --maxp) MAXP=$2; shift 2;;
    --staging) STAGING=$2; shift 2;;
    --sif) SIF=$2; shift 2;;
    --covariates) COVARIATES=$2; shift 2;;
    --stages) STAGES=$2; shift 2;;
    --id-file) ID_FILE=$2; shift 2;;
    --copy) COPIES+=("$2"); shift 2;;
    --seed) SEED=$2; shift 2;;
    --param) EXTRA_PARAMS+="	$2"$'\n'; shift 2;;
    --container-env) IFS='=' read -r k v <<< "$2"; EXTRA_ENV+="	$k = \"$v\""$'\n'; shift 2;;
    *) echo "unknown argument $1" >&2; exit 2;;
  esac
done
: "${LABEL:?} ${SHA:?} ${RUN:?} ${KTYPE:?} ${KLEN:?} ${MAXP:?} ${STAGING:?}"
[ "$(cat "$STAGING/STAGED_SHA")" = "$SHA" ] || { echo "--sha does not match $STAGING/STAGED_SHA" >&2; exit 1; }
NF=$STAGING/kmer_pipeline.nf

# example/ comes from the staged commit's own tree now, not extracted from the
# image -- the image's own copy still exists (code-full images only) and is
# checked byte-identical, scoped to the data files only (B9, guards P18): the
# id_file and software-location file are NOT checked here, since B6/B2 already
# rewrote both deliberately (bare paths; split scriptpath/toolpath).
if [ -n "$SIF" ] && apptainer exec --containall --cleanenv "$SIF" test -d /usr/share/kmer_pipeline/example 2>/dev/null; then
  for f in "$STAGING"/example/*.fa.gz "$STAGING"/example/*.fasta "$STAGING"/example/*.gb; do
    [ -e "$f" ] || continue
    b=$(basename "$f")
    img_sum=$(apptainer exec --containall --cleanenv "$SIF" sha256sum "/usr/share/kmer_pipeline/example/$b" 2>/dev/null | cut -d' ' -f1) || img_sum=""
    repo_sum=$(sha256sum "$f" | cut -d' ' -f1)
    if [ -n "$img_sum" ] && [ "$img_sum" != "$repo_sum" ]; then
      echo "example/$b differs between the staged commit and the image's own bundled copy" >&2
      exit 1
    fi
  done
fi

OUT=$KMER_E2E_ROOT/goldens/$LABEL-${SHA:0:7}/maxp$MAXP/$RUN
[ -e "$OUT" ] && { echo "$OUT already exists" >&2; exit 1; }
BASE=$KMER_E2E_ROOT/runs/$LABEL-${SHA:0:7}-maxp$MAXP-$RUN-$(date +%Y%m%d-%H%M%S)
mkdir -p "$BASE/tb20" "$OUT"
PREFIX=tb20_$KTYPE$KLEN
ANALYSIS=$BASE/tb20/kmergwas

cp -r "$STAGING/example/." "$BASE/tb20/"
chmod -R u+w "$BASE/tb20"
[ -n "$ID_FILE" ] && sed "s,@BASE@,$BASE,g" "$ID_FILE" > "$BASE/tb20/id_file.txt"
for f in "${COPIES[@]}"; do cp "$f" "$BASE/tb20/"; done

if [ -n "$SEED" ]; then
  mkdir -p "$ANALYSIS"
  rsync -a --exclude 'log.*' --exclude nextflow.out "$SEED/" "$ANALYSIS/"
fi

COV_PARAMS=""
if [ -n "$COVARIATES" ]; then
  cp "$COVARIATES" "$BASE/tb20/covariates.txt"
  COV_PARAMS="covariate_file = \"$BASE/tb20/covariates.txt\""
fi

write_config() {  # $1 = stage to run (1..7)
  local skips="" k
  for k in 1 2 3 4 5 6 7; do
    if [ "$k" = "$1" ]; then skips+="	skip$k = false"$'\n'; else skips+="	skip$k = true"$'\n'; fi
  done
  cat > "$BASE/nextflow.config" <<EOF
params {
	base_dir = "$BASE"
	output_prefix = "tb20"
	analysis_dir = "\$base_dir/\$output_prefix/kmergwas"
	kmer_type = "$KTYPE"
	kmer_length = $KLEN
	id_file = "\$base_dir/\$output_prefix/id_file.txt"
	ref_fa = "\$base_dir/\$output_prefix/Mtub_H37Rv_NC000962.3.fasta"
	ref_gb = "\$base_dir/\$output_prefix/Mtub_H37Rv_NC000962.3.gb"
	maxp = $MAXP
	$COV_PARAMS
$EXTRA_PARAMS$skips}
executor.queueSize = params.maxp
executor.cpus = params.maxp
EOF
  if [ -n "$EXTRA_ENV" ]; then
    # Merges with (does not replace) the staged commit's own env{} scope -- OPENBLAS_NUM_THREADS/
    # OMP_NUM_THREADS stay pinned regardless (B9); this is for a per-run override such as
    # PYTHONHASHSEED (step_test.sh's determinism checks).
    printf 'env {\n%s}\n' "$EXTRA_ENV" >> "$BASE/nextflow.config"
  fi
}

check_stage() {  # hidden failures: before a fix, R errors could exit 0
  local bad=0
  if find "$BASE/work" -name Rcoredump.rda | grep -q .; then echo "Rcoredump.rda found" >&2; bad=1; fi
  if grep -l -E 'Execution halted|^Error in |^Error:' "$ANALYSIS"/log.*/*.log 2>/dev/null | grep -q .; then
    echo "R errors in logs:" >&2; grep -l -E 'Execution halted|^Error in |^Error:' "$ANALYSIS"/log.*/*.log >&2; bad=1
  fi
  return $bad
}

snapshot() {  # $1 = stage; logs dereferenced, work.* link dir excluded
  rsync -a --copy-links --exclude "work.$PREFIX" "$ANALYSIS/" "$OUT/stage$1/"
}

export NXF_HOME=$KMER_E2E_ROOT/nxf_home NXF_OPTS="-Dnxf.ansi.log=false"
for stage in $STAGES; do
  write_config "$stage"
  echo "== stage $stage $(date)"
  (cd "$BASE" && nice -n "$KMER_E2E_NICE" "$NEXTFLOW" run "$NF" -profile apptainer --container "$SIF" -ansi-log false) > "$BASE/stage$stage.out" 2>&1
  check_stage
  snapshot "$stage"
  cp "$BASE/stage$stage.out" "$OUT/stage$stage/nextflow.out"
done

n=$(($(wc -l < "$BASE/tb20/id_file.txt") - 1))
# OPENBLAS_NUM_THREADS/OMP_NUM_THREADS are pinned unconditionally by the staged commit's own
# nextflow.config env{} scope (B9) now, not a per-run --container-env flag -- recorded as a fixed
# note here for provenance (guards P3), not a value this script chooses any more.
python3 - "$OUT/metadata.json" <<EOF
import json, math, sys
n, maxp = $n, $MAXP
json.dump({"label": "$LABEL", "sha": "$SHA", "run": "$RUN", "kmer_type": "$KTYPE",
  "kmer_length": $KLEN, "maxp": maxp, "n": n,
  "p": max(1, min(math.ceil(n/2), maxp)), "p5": max(1, min(math.ceil(n/5), maxp)),
  "staging": "$STAGING", "covariates": "$COVARIATES", "stages": "$STAGES", "id_file": "$ID_FILE", "copies": "${COPIES[*]}", "seed": "$SEED",
  "thread_pinning": "OPENBLAS_NUM_THREADS=1,OMP_NUM_THREADS=1 (nextflow.config env{}, always on)",
  "extra_params": $(printf %s "$EXTRA_PARAMS" | python3 -c "import json,sys; print(json.dumps(sys.stdin.read().strip()))"), "host": "$(hostname)", "nextflow": "$KMER_E2E_NEXTFLOW_VERSION",
  "container": "$SIF", "run_dir": "$BASE", "date": "$(date -Iseconds)",
  "stage_flags": {"1": "skip1 countkmers", "2": "skip2 createfullkmerlist",
    "3": "skip3 stringlist2patternandkinship", "4": "skip4 rungemma",
    "5": "skip5 kmercontigalign + kmercontigalignmerge", "6": "skip6 plotManhattan",
    "7": "skip7 gen*Report"}}, open(sys.argv[1], "w"), indent=1)
EOF
echo "== done: $OUT"
exit
}
