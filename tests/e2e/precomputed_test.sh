#!/usr/bin/env bash
# precomputed_dir end-to-end test (renamed from the private harness's n5_test.sh,
# which named it after a private plan item). Steps 4, 6 and 7 run with
# precomputed_dir (the stage-5 snapshot of a full run of the example) and a
# phenotype from pheno_file must give the same step 4, 6 and 7 outputs as a
# full run with that phenotype; precomputed_dir must be left unchanged.
#
# Kept tied to the tb20/nucleotide/31 triple, since it is the only example
# dataset this harness has.
#
# Usage: precomputed_test.sh NAME SHA PRECOMPUTED_STAGE5 PHENO_ID_FILE FULL_STAGE7 [COVARIATES]
#   PRECOMPUTED_STAGE5  stage-5 snapshot of a full run (same SHA) with any phenotype
#   PHENO_ID_FILE       an id file (as make_reference.sh --id-file) whose id and pheno columns
#                       become pheno_file; id_file stays the example's
#   FULL_STAGE7         stage-7 snapshot of a full run with PHENO_ID_FILE (same SHA)
# Output: $KMER_E2E_ROOT/runs/precomputed/NAME/{run.out,compare.txt,result.txt}
set -uo pipefail
{
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
# shellcheck source=config.sh
. "$HERE/config.sh"
NAME=${1:?}; SHA=${2:?}; PRE=${3:?}; PHENO_ID=${4:?}; FULL=${5:?}; COV=${6:-}
STAGING=$KMER_E2E_ROOT/staging/$SHA
BASE=$KMER_E2E_ROOT/runs/precomputed/$NAME
[ -e "$BASE" ] && { echo "$BASE exists" >&2; exit 1; }
NF=$STAGING/kmer_pipeline.nf
mkdir -p "$BASE/tb20" "$BASE/pre"
cp -r "$STAGING/example/." "$BASE/tb20/"
chmod -R u+w "$BASE/tb20"
cut -f1,3 "$PHENO_ID" > "$BASE/tb20/pheno_file.txt"
rsync -a --exclude 'log.*' --exclude nextflow.out "$PRE/" "$BASE/pre/kmergwas/"
(cd "$BASE/pre/kmergwas" && find . -type f | sort | xargs md5sum) > "$BASE/pre_before.md5"
COV_LINE=""
if [ -n "$COV" ]; then cp "$COV" "$BASE/tb20/covariates.txt"; COV_LINE="covariate_file = \"$BASE/tb20/covariates.txt\""; fi
cat > "$BASE/nextflow.config" <<EOF
params {
	base_dir = "$BASE"
	output_prefix = "tb20"
	analysis_dir = "\$base_dir/\$output_prefix/kmergwas"
	kmer_type = "nucleotide"
	kmer_length = 31
	id_file = "\$base_dir/\$output_prefix/id_file.txt"
	pheno_file = "\$base_dir/\$output_prefix/pheno_file.txt"
	precomputed_dir = "\$base_dir/pre/kmergwas"
	ref_fa = "\$base_dir/\$output_prefix/Mtub_H37Rv_NC000962.3.fasta"
	ref_gb = "\$base_dir/\$output_prefix/Mtub_H37Rv_NC000962.3.gb"
	maxp = 2
	$COV_LINE
}
executor.queueSize = params.maxp
executor.cpus = params.maxp
EOF
export NXF_HOME=$KMER_E2E_ROOT/nxf_home NXF_OPTS="-Dnxf.ansi.log=false"
cd "$BASE" && "$KMER_E2E_NEXTFLOW" run "$NF" -profile apptainer --container "$KMER_E2E_SIF" -ansi-log false > run.out 2>&1
echo "run exit $?" > result.txt
(cd "$BASE/pre/kmergwas" && find . -type f | sort | xargs md5sum) > "$BASE/pre_after.md5"
cmp -s "$BASE/pre_before.md5" "$BASE/pre_after.md5" && echo "precomputed_dir unchanged" >> result.txt \
  || echo "PROBLEM: precomputed_dir changed" >> result.txt
# Compare the files of steps 4, 6 and 7 with the full run's
mkdir -p "$BASE/full_467"
(cd "$FULL" && find . -type f | sed 's|^\./||') | \
  apptainer exec --bind "$STAGING,$FULL" "$KMER_E2E_SIF" python3 -c "
import sys; sys.path.insert(0, '$STAGING/bin')
import inventory
inv = inventory.Inventory('tb20', 'nucleotide', 31)
for line in sys.stdin:
    p = line.strip()
    if inv.owner(p) in (4, 6, 7):
        print(p)" > "$BASE/full_467.txt"
rsync -a --files-from="$BASE/full_467.txt" "$FULL/" "$BASE/full_467/"
rsync -a --copy-links --exclude 'work.*' --exclude 'log.*' --exclude '*.container_id_file.txt' --exclude '*.analysis_file.txt' "$BASE/tb20/kmergwas/" "$BASE/new/"
# FULL's own run_dir (base_dir), for the manifest's "path" class (B9) -- FULL and this run are
# different run directories, so without this the path-bearing files would always DIFFER.
full_meta="$(dirname "$FULL")/metadata.json"
full_prefix=$([ -f "$full_meta" ] && python3 -c "import json; print(json.load(open('$full_meta')).get('run_dir',''))" || true)
python3 "$HERE/compare.py" "$BASE/full_467" "$BASE/new" --prefix-a "$full_prefix" --prefix-b "$BASE" > compare.txt 2>&1
echo "compare exit $? ($(tail -n 1 compare.txt))" >> result.txt
cat result.txt
exit
}
