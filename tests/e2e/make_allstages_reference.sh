#!/usr/bin/env bash
# A single-invocation, all-stages run (no skip* flags at all), unlike
# make_reference.sh's stage-by-stage loop -- for comparing the decomposed,
# stage-by-stage DAG against the pipeline's own normal, non-decomposed
# default. Only writes a stage7 snapshot (the only one golden_hashes.py
# needs). Ported from a scratch-only helper used at C12 (B9, guards P44: that
# version had no committed generating script).
#
# Usage: make_allstages_reference.sh --label LABEL --sha SHA --run NAME
#                                     --kmer-type TYPE --kmer-length K --maxp M --staging DIR [--sif IMAGE]
# Output: $KMER_E2E_ROOT/goldens/LABEL-SHA7/maxpM/NAME/stage7
set -euo pipefail
{ # read whole script before running, so edits can't affect a running copy

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
# shellcheck source=config.sh
. "$HERE/config.sh"
SIF=$KMER_E2E_SIF
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
    *) echo "unknown argument $1" >&2; exit 2;;
  esac
done
: "${LABEL:?} ${SHA:?} ${RUN:?} ${KTYPE:?} ${KLEN:?} ${MAXP:?} ${STAGING:?}"
[ "$(cat "$STAGING/STAGED_SHA")" = "$SHA" ] || { echo "--sha does not match $STAGING/STAGED_SHA" >&2; exit 1; }
NF=$STAGING/kmer_pipeline.nf

OUT=$KMER_E2E_ROOT/goldens/$LABEL-${SHA:0:7}/maxp$MAXP/$RUN
[ -e "$OUT" ] && { echo "$OUT already exists" >&2; exit 1; }
BASE=$KMER_E2E_ROOT/runs/$LABEL-${SHA:0:7}-maxp$MAXP-$RUN-$(date +%Y%m%d-%H%M%S)
mkdir -p "$BASE/tb20" "$OUT"
ANALYSIS=$BASE/tb20/kmergwas
PREFIX=tb20_$KTYPE$KLEN

cp -r "$STAGING/example/." "$BASE/tb20/"
chmod -R u+w "$BASE/tb20"

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
}
executor.queueSize = params.maxp
executor.cpus = params.maxp
EOF

export NXF_HOME=$KMER_E2E_ROOT/nxf_home NXF_OPTS="-Dnxf.ansi.log=false"
(cd "$BASE" && nice -n "$KMER_E2E_NICE" "$KMER_E2E_NEXTFLOW" run "$NF" -profile apptainer --container "$SIF" \
  -ansi-log false -with-trace trace.txt) > "$BASE/run.out" 2>&1
rsync -a --copy-links --exclude "work.$PREFIX" "$ANALYSIS/" "$OUT/stage7/"
cp "$BASE/run.out" "$OUT/stage7/nextflow.out"
echo "== done: $OUT"
exit
}
