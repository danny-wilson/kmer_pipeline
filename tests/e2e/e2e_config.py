#!/usr/bin/env python3
"""Read tests/e2e/local.conf (gitignored; see local.conf.example): KEY=value
lines, one per line, optionally quoted, '#' comments. Environment variables of
the same name override the file. Used by compare.py and golden_hashes.py so
the shell scripts (config.sh) and the Python tools agree on where things are.

Usage: e2e_config.py            print every resolved key=value
       load() / require(key)    from Python
"""
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REQUIRED = ("KMER_E2E_ROOT", "KMER_E2E_SIF", "KMER_E2E_IMAGE_COMMIT",
            "KMER_E2E_NEXTFLOW", "KMER_E2E_NEXTFLOW_VERSION")
DEFAULTS = {"KMER_E2E_LOCAL_TMP": "/tmp", "KMER_E2E_NICE": "10", "KMER_E2E_PRE": ""}
_LINE = re.compile(r'^([A-Za-z_][A-Za-z0-9_]*)=(.*)$')


def _unquote(value):
    value = value.strip()
    if len(value) >= 2 and value[0] == value[-1] and value[0] in "\"'":
        return value[1:-1]
    return value


def load(conf_path=None):
    """Every resolved key, as a dict. Exits with a clear message if a required
    key is still unset after the file and the environment are applied."""
    conf_path = conf_path or os.environ.get("KMER_E2E_CONF") or os.path.join(HERE, "local.conf")
    values = dict(DEFAULTS)
    if os.path.exists(conf_path):
        with open(conf_path) as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                m = _LINE.match(line)
                if m:
                    values[m.group(1)] = _unquote(m.group(2))
    for k in set(values) | set(REQUIRED):
        if k in os.environ:
            values[k] = os.environ[k]
    missing = [k for k in REQUIRED if not values.get(k)]
    if missing:
        sys.exit(f"e2e config: missing {', '.join(missing)} "
                  f"(see tests/e2e/local.conf.example, or set them as environment variables)")
    return values


def require(key, conf_path=None):
    return load(conf_path)[key]


if __name__ == "__main__":
    for k, v in sorted(load().items()):
        print(f"{k}={v}")
