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
BASE=$KMER_E2E_ROOT/runs/resume_rules/${SHA:0:7}
[ -e "$BASE" ] && { echo "$BASE exists" >&2; exit 1; }
mkdir -p "$BASE/tb20"
(cd "$BASE/tb20" && apptainer exec --containall --cleanenv "$KMER_E2E_SIF" \
  bash -c 'cd /usr/share/kmer_pipeline/example && tar -c .' | tar -x)
sed -i "s,/usr/share/kmer_pipeline/example/,$BASE/tb20/,g" "$BASE/tb20/id_file.txt"
cp "$STAGING/kmer_pipeline.nf" "$BASE/"
sed "s,^scriptpath\t.*,scriptpath\t$STAGING," "$STAGING/pipeline_software_location.txt" > "$BASE/software.txt"
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
	container_type = "singularity"
	container_file = "$KMER_E2E_SIF"
	software_file = "$BASE/software.txt"
	container_args = "--bind $STAGING --env PYTHONNOUSERSITE=1"
$2
}
executor.queueSize = params.maxp
executor.cpus = params.maxp
EOF
}
export NXF_HOME=$KMER_E2E_ROOT/nxf_home NXF_OPTS="-Dnxf.ansi.log=false"
cd "$BASE"
run() {  # $1 name, $2 expectation (ok|refused), $3 grep pattern for the expected message, rest: nextflow args
	local name=$1 want=$2 pat=$3; shift 3
	"$KMER_E2E_NEXTFLOW" run kmer_pipeline.nf -ansi-log false "$@" > "$name.out" 2>&1; local rc=$?
	local submitted=$(grep -c "Submitted process" "$name.out")
	local cached=$(grep -c "Cached process" "$name.out")
	local got=ok; [ $rc -ne 0 ] && got=refused
	local verdict=PASS
	[ "$got" != "$want" ] && verdict=FAIL
	[ -n "$pat" ] && ! grep -q -- "$pat" "$name.out" && verdict=FAIL
	echo "$verdict $name: exit $rc, submitted $submitted, cached $cached (expected $want${pat:+, \"$pat\"})" | tee -a result.txt
}
config kmergwas ""
run full ok ""
cp tb20/id_file.txt tb20/id_file.orig
awk 'BEGIN{FS=OFS="\t"} NR==2{$3=$3+1} {print}' tb20/id_file.orig > tb20/id_file.txt
run f1 refused "changed since that run: id_file" -resume
cp tb20/id_file.orig tb20/id_file.txt
config kmergwas "	minor_allele_threshold = 0.05"
run f2 refused "changed since that run: minor_allele_threshold" -resume
config kmergwas "	merge_wait_minutes = 200"
run f3 ok "" -resume
config kmergwas ""; echo "resume = true" >> nextflow.config
run g1 ok ""
config kmergwas "	minor_allele_threshold = 0.05"; echo "resume = true" >> nextflow.config
run g2 refused "changed since that run"
config kmergwas_other ""
run h refused "no run to continue" -resume
config kmergwas ""
run n refused "set overwrite = true"
cat result.txt
exit
}
