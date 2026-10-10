#!/usr/bin/env bash
# Generate a reference run: the example pipeline run one Nextflow stage at a time
# (skip1..skip7), with a snapshot of analysis_dir after each stage.
#
# Usage: make_reference.sh --label LABEL --sha SHA --run NAME --kmer-type TYPE
#                          --kmer-length K --maxp M [--staging DIR] [--nf FILE]
#                          [--covariates FILE] [--stages "1 2 3 4 5 6 7"] [--param "name = value"]...
#                          [--id-file FILE (paths may use @BASE@ for the base directory)] [--copy FILE]...
#                          [--seed DIR] [--container-env VAR=VALUE]... [--sif IMAGE]
#
# --seed copies a stage snapshot (logs and nextflow.out excluded) into analysis_dir
# before the first stage, for step tests. --container-env adds to the environment
# of the container (e.g. PYTHONHASHSEED=0); needs --staging. Without --staging,
# the image's own scripts and kmer_pipeline.nf are used, so the run only exercises
# whatever commit the image itself was built from (check with --sif a release
# image, whose scripts check_image.sh has checked against that commit).
# Output: $KMER_E2E_ROOT/goldens/LABEL-SHA7/maxpM/NAME/{stage1..stage7,metadata.json}
set -euo pipefail
{ # read whole script before running, so edits can't affect a running copy

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
# shellcheck source=config.sh
. "$HERE/config.sh"
SIF=$KMER_E2E_SIF
NEXTFLOW=$KMER_E2E_NEXTFLOW
STAGING=""; NF=""; COVARIATES=""; STAGES="1 2 3 4 5 6 7"; EXTRA_PARAMS=""; ID_FILE=""; COPIES=(); SEED=""; CENV=""
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
    --nf) NF=$2; shift 2;;
    --covariates) COVARIATES=$2; shift 2;;
    --stages) STAGES=$2; shift 2;;
    --id-file) ID_FILE=$2; shift 2;;
    --copy) COPIES+=("$2"); shift 2;;
    --seed) SEED=$2; shift 2;;
    --container-env) CENV+=",$2"; shift 2;;
    --param) EXTRA_PARAMS+="	$2"$'\n'; shift 2;;
    *) echo "unknown argument $1" >&2; exit 2;;
  esac
done
: "${LABEL:?} ${SHA:?} ${RUN:?} ${KTYPE:?} ${KLEN:?} ${MAXP:?}"
if [ -n "$STAGING" ]; then
  [ "$(cat "$STAGING/STAGED_SHA")" = "$SHA" ] || { echo "--sha does not match $STAGING/STAGED_SHA" >&2; exit 1; }
  [ -z "$NF" ] && NF=$STAGING/kmer_pipeline.nf
fi

OUT=$KMER_E2E_ROOT/goldens/$LABEL-${SHA:0:7}/maxp$MAXP/$RUN
[ -e "$OUT" ] && { echo "$OUT already exists" >&2; exit 1; }
BASE=$KMER_E2E_ROOT/runs/$LABEL-${SHA:0:7}-maxp$MAXP-$RUN-$(date +%Y%m%d-%H%M%S)
mkdir -p "$BASE/tb20" "$OUT"
PREFIX=tb20_$KTYPE$KLEN
ANALYSIS=$BASE/tb20/kmergwas

# Inputs: the example files from the image, paths rewritten to the user file system
(cd "$BASE/tb20" && apptainer exec --containall --cleanenv "$SIF" \
  bash -c 'cd /usr/share/kmer_pipeline/example && tar -c .' | tar -x)
sed -i "s,/usr/share/kmer_pipeline/example/,$BASE/tb20/,g" "$BASE/tb20/id_file.txt"
[ -n "$ID_FILE" ] && sed "s,@BASE@,$BASE,g" "$ID_FILE" > "$BASE/tb20/id_file.txt"
for f in "${COPIES[@]}"; do cp "$f" "$BASE/tb20/"; done
if [ -n "$NF" ]; then cp "$NF" "$BASE/kmer_pipeline.nf"
else apptainer exec --containall --cleanenv "$SIF" cat /usr/local/bin/kmer_pipeline.nf > "$BASE/kmer_pipeline.nf"; fi

if [ -n "$SEED" ]; then
  mkdir -p "$ANALYSIS"
  rsync -a --exclude 'log.*' --exclude nextflow.out "$SEED/" "$ANALYSIS/"
fi

COV_PARAMS=""
if [ -n "$COVARIATES" ]; then
  cp "$COVARIATES" "$BASE/tb20/covariates.txt"
  COV_PARAMS="covariate_file = \"$BASE/tb20/covariates.txt\""
fi

STAGING_PARAMS=""
if [ -n "$STAGING" ]; then
  # Software file must sit under base_dir; scriptpath is the staging dir, bound at the same path
  sed "s,^scriptpath\t.*,scriptpath\t$STAGING," "$STAGING/pipeline_software_location.txt" \
    > "$BASE/pipeline_software_location.txt"
  STAGING_PARAMS="software_file = \"$BASE/pipeline_software_location.txt\"
	container_args = \"--bind $STAGING --env PYTHONNOUSERSITE=1$CENV\""
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
	container_type = "singularity"
	container_file = "$SIF"
	$STAGING_PARAMS
	$COV_PARAMS
$EXTRA_PARAMS$skips}
executor.queueSize = params.maxp
executor.cpus = params.maxp
EOF
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
  (cd "$BASE" && nice -n "$KMER_E2E_NICE" "$NEXTFLOW" run kmer_pipeline.nf -ansi-log false) > "$BASE/stage$stage.out" 2>&1
  check_stage
  snapshot "$stage"
  cp "$BASE/stage$stage.out" "$OUT/stage$stage/nextflow.out"
done

n=$(($(wc -l < "$BASE/tb20/id_file.txt") - 1))
python3 - "$OUT/metadata.json" <<EOF
import json, math, sys
n, maxp = $n, $MAXP
json.dump({"label": "$LABEL", "sha": "$SHA", "run": "$RUN", "kmer_type": "$KTYPE",
  "kmer_length": $KLEN, "maxp": maxp, "n": n,
  "p": max(1, min(math.ceil(n/2), maxp)), "p5": max(1, min(math.ceil(n/5), maxp)),
  "staging": "$STAGING", "covariates": "$COVARIATES", "stages": "$STAGES", "id_file": "$ID_FILE", "copies": "${COPIES[*]}", "seed": "$SEED", "container_env": "$CENV", "extra_params": $(printf %s "$EXTRA_PARAMS" | python3 -c "import json,sys; print(json.dumps(sys.stdin.read().strip()))"), "host": "$(hostname)", "nextflow": "$KMER_E2E_NEXTFLOW_VERSION",
  "container": "$SIF", "run_dir": "$BASE", "date": "$(date -Iseconds)",
  "stage_flags": {"1": "skip1 countkmers", "2": "skip2 createfullkmerlist",
    "3": "skip3 stringlist2patternandkinship", "4": "skip4 rungemma",
    "5": "skip5 kmercontigalign + kmercontigalignmerge", "6": "skip6 plotManhattan",
    "7": "skip7 gen*Report"}}, open(sys.argv[1], "w"), indent=1)
EOF
echo "== done: $OUT"
exit
}
