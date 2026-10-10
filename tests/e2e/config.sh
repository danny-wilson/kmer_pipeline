#!/usr/bin/env bash
# Load tests/e2e/local.conf (gitignored; see local.conf.example). Environment
# variables of the same name override the file. Source this file; don't run it.
#   . "$(dirname "$0")/config.sh"
# Required keys end up exported and non-empty, or this exits with a clear
# message. KMER_E2E_PRE (if set) is eval'd here, before Nextflow's version is
# checked, so any environment setup it does takes effect for the rest of the script.
E2E_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
CONF=${KMER_E2E_CONF:-$E2E_DIR/local.conf}

if [ -f "$CONF" ]; then
  set -a
  # shellcheck disable=SC1090
  source "$CONF"
  set +a
fi

: "${KMER_E2E_LOCAL_TMP:=/tmp}"
: "${KMER_E2E_NICE:=10}"

for k in KMER_E2E_ROOT KMER_E2E_SIF KMER_E2E_IMAGE_COMMIT KMER_E2E_NEXTFLOW KMER_E2E_NEXTFLOW_VERSION; do
  if [ -z "${!k:-}" ]; then
    echo "e2e config: $k is not set (see tests/e2e/local.conf.example, or set it as an environment variable)" >&2
    return 1 2>/dev/null || exit 1
  fi
done

if [ -n "${KMER_E2E_PRE:-}" ]; then eval "$KMER_E2E_PRE"; fi

got_version=$("$KMER_E2E_NEXTFLOW" -version 2>&1 | grep -o '[0-9][0-9.]*' | head -n 1 || true)
if [ "$got_version" != "$KMER_E2E_NEXTFLOW_VERSION" ]; then
  echo "e2e config: KMER_E2E_NEXTFLOW -version reports '$got_version', expected '$KMER_E2E_NEXTFLOW_VERSION'" >&2
  return 1 2>/dev/null || exit 1
fi

export KMER_E2E_ROOT KMER_E2E_SIF KMER_E2E_IMAGE_COMMIT KMER_E2E_NEXTFLOW KMER_E2E_NEXTFLOW_VERSION \
       KMER_E2E_PRE KMER_E2E_LOCAL_TMP KMER_E2E_NICE
