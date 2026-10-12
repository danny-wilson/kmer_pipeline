#!/usr/bin/env bash
# -resume rule tests, on the local executor with a staged commit. One full run, then:
#   f1  -resume with a phenotype changed in id_file            -> refused (inputs changed)
#   f2  -resume with minor_allele_threshold changed            -> refused (parameter changed)
#   f3  -resume with merge_wait_minutes changed (operational)   -> allowed, everything cached
#   g1  resume = true in the config, nothing changed            -> allowed, everything cached
#   g2  resume = true in the config, minor_allele_threshold changed -> refused
#   h   -resume with another analysis_dir                       -> refused (nothing to continue)
#   n   no -resume, same analysis_dir                           -> stops: overwrite = false
# Usage: resume_rules_test.sh SHA
# Output: $KMER_E2E_ROOT/runs/resume_rules/SHA7/{*.out,result.txt}
set -uo pipefail
{
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
# shellcheck source=config.sh
. "$HERE/config.sh"
SHA=${1:?usage: resume_rules_test.sh SHA}
STAGING=$KMER_E2E_ROOT/staging/$SHA
[ -d "$STAGING" ] || "$HERE/stage.sh" "$SHA" > /dev/null
NF=$STAGING/kmer_pipeline.nf
BASE=$KMER_E2E_ROOT/runs/resume_rules/${SHA:0:7}
[ -e "$BASE" ] && { echo "$BASE exists" >&2; exit 1; }
mkdir -p "$BASE/tb20"
cp -r "$STAGING/example/." "$BASE/tb20/"
chmod -R u+w "$BASE/tb20"
config() {  # $1 analysis folder name, $2 extra lines
	cat > "$BASE/nextflow.config" <<EOF
params {
	base_dir = "$BASE"
	output_prefix = "tb20"
	analysis_dir = "\$base_dir/\$output_prefix/$1"
	kmer_type = "nucleotide"
	kmer_length = 31
	id_file = "\$base_dir/\$output_prefix/id_file.txt"
	ref_fa = "\$base_dir/\$output_prefix/Mtub_H37Rv_NC000962.3.fasta"
	ref_gb = "\$base_dir/\$output_prefix/Mtub_H37Rv_NC000962.3.gb"
	maxp = 2
$2
}
executor.queueSize = params.maxp
executor.cpus = params.maxp
EOF
}
export NXF_HOME=$KMER_E2E_ROOT/nxf_home NXF_OPTS="-Dnxf.ansi.log=false"
cd "$BASE"
# PLAN.md Part B, B5 (guards P11): f3/g1 take an optional 4th/5th argument (expect_submitted,
# expect_cached) asserting the exact submitted/cached task counts, not just exit code and a text
# pattern - a content check (does the trace show the full remaining task count actually cached,
# not merely that *something* was cached) rather than today's weaker presence check.
run() {  # $1 name, $2 expectation (ok|refused), $3 grep pattern for the expected message,
         # $4 expect_submitted (optional), $5 expect_cached (optional), rest: nextflow args
	local name=$1 want=$2 pat=$3 expect_submitted=${4:-} expect_cached=${5:-}
	shift 5
	"$KMER_E2E_NEXTFLOW" run "$NF" -profile apptainer --container "$KMER_E2E_SIF" -ansi-log false "$@" > "$name.out" 2>&1; local rc=$?
	local submitted=$(grep -c "Submitted process" "$name.out")
	local cached=$(grep -c "Cached process" "$name.out")
	local got=ok; [ $rc -ne 0 ] && got=refused
	local verdict=PASS
	[ "$got" != "$want" ] && verdict=FAIL
	[ -n "$pat" ] && ! grep -q -- "$pat" "$name.out" && verdict=FAIL
	[ -n "$expect_submitted" ] && [ "$submitted" != "$expect_submitted" ] && verdict=FAIL
	[ -n "$expect_cached" ] && [ "$cached" != "$expect_cached" ] && verdict=FAIL
	echo "$verdict $name: exit $rc, submitted $submitted, cached $cached (expected $want${pat:+, \"$pat\"}${expect_submitted:+, submitted=$expect_submitted}${expect_cached:+, cached=$expect_cached})" | tee -a result.txt
}
config kmergwas ""
run full ok "" "" ""
FULL_TASK_COUNT=$(grep -c "Submitted process" full.out)
cp tb20/id_file.txt tb20/id_file.orig
awk 'BEGIN{FS=OFS="\t"} NR==2{$3=$3+1} {print}' tb20/id_file.orig > tb20/id_file.txt
run f1 refused "changed since that run: id_file" "" "" -resume
cp tb20/id_file.orig tb20/id_file.txt
config kmergwas "	minor_allele_threshold = 0.05"
run f2 refused "changed since that run: minor_allele_threshold" "" "" -resume
config kmergwas "	merge_wait_minutes = 200"
# preflight only (operational parameter, not a result-affecting one, so everything else is
# unchanged and fully reused): submitted==1, cached==every other task from the full run.
run f3 ok "" 1 "$((FULL_TASK_COUNT - 1))" -resume
config kmergwas ""; echo "resume = true" >> nextflow.config
run g1 ok "" 1 "$((FULL_TASK_COUNT - 1))"
config kmergwas "	minor_allele_threshold = 0.05"; echo "resume = true" >> nextflow.config
run g2 refused "changed since that run" "" ""
config kmergwas_other ""
run h refused "no run to continue" "" "" -resume
config kmergwas ""
run n refused "set overwrite = true" "" ""
cat result.txt
exit
}
