#!/usr/bin/env bash
# Run step_test.sh for one stage on every reference run under REF_LABEL (not
# bowtie), with PYTHONHASHSEED 0 and 1, JOBS at a time (default 4). Prints one
# summary line per test.
# Usage: run_step_tests.sh STAGING_DIR STAGE REF_LABEL
set -uo pipefail
{ # read whole script before running, so edits can't affect a running copy
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
# shellcheck source=config.sh
. "$HERE/config.sh"
STAGING=${1:?usage: run_step_tests.sh STAGING_DIR STAGE REF_LABEL}
STAGE=${2:?usage: run_step_tests.sh STAGING_DIR STAGE REF_LABEL}
LABEL=${3:?usage: run_step_tests.sh STAGING_DIR STAGE REF_LABEL}
export HERE STAGING STAGE
ls -d "$KMER_E2E_ROOT/goldens/$LABEL"/maxp*/* | grep -v -e bowtie -e '\.make\.log$' | while read -r ref; do
  for h in 0 1; do echo "$ref $h"; done
done | xargs -P "${JOBS:-4}" -L 1 bash -c '
  out=$("$HERE/step_test.sh" "$STAGING" "$0" "$STAGE" "$1" 2>&1); rc=$?
  echo "$(basename "$(dirname "$0")")/$(basename "$0") h$1 rc=$rc: $(echo "$out" | tail -n 1)"
  [ $rc -eq 0 ] || echo "$out" | head -n 20 | sed "s/^/    /"'
exit
}
