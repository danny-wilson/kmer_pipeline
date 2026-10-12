"""nextflow.config (PLAN.md Part B, B7): USER_KEYS excludes the repo's own config keys
(container, extra_binds), and preflight.py's KNOWN_PARAMS learns both names - the plan asks for
both, belt-and-suspenders. Checked directly (not assumed): with the current feature set, the two
fixes are redundant - removing the USER_KEYS exclusion alone does not make either test below fail,
because KNOWN_PARAMS already covers both names regardless of what USER_KEYS contains. The
exclusion is kept anyway, as the plan specifies, since nothing here guarantees KNOWN_PARAMS will
always be the only thing that reads USER_KEYS. Recorded as a genuine finding, not claimed as a
mutant either test guards."""
import os

import pytest

from conftest import REPO_DIR
from test_parameters import NEXTFLOW_RUNS
from test_workflow_start import make_base, run

NF = os.path.join(REPO_DIR, "kmer_pipeline.nf")
CONFIG = os.path.join(REPO_DIR, "nextflow.config")
pytestmark = pytest.mark.skipif(not NEXTFLOW_RUNS or not os.path.exists(NF), reason="needs Nextflow and the .nf")


def _make_base_with_repo_config(tmp_path):
    """Like test_workflow_start.make_base(), but the real nextflow.config (not a from-scratch
    one) holds the analysis params too - the repo's own params.container/extra_binds defaults
    must survive alongside them."""
    base = make_base(tmp_path)
    cfg = base / "nextflow.config"
    analysis_params = cfg.read_text()
    cfg.write_text(open(CONFIG).read())
    # Nextflow's own precedence (a launch-directory nextflow.config, parsed as one file here):
    # later params{} assignments of the same key override earlier ones. Appending the test's own
    # analysis params after the repo's config matches "the repo's config, then the launch
    # directory's" in spirit, since here they are the same file.
    with open(cfg, "a") as fh:
        fh.write(analysis_params)
    return base


def test_container_is_not_reported_as_an_unknown_parameter(tmp_path):
    base = _make_base_with_repo_config(tmp_path)
    with open(base / "nextflow.config", "a") as fh:
        for s in range(1, 8):
            fh.write(f"params.skip{s} = true\n")
    rc, out = run(base)
    assert "'container' is not used" not in out, out[-3000:]
    assert "unknown parameter 'container'" not in out, out[-3000:]


def test_a_typo_still_errors(tmp_path):
    """--kmer-lenght (a typo for --kmer-length) must still be caught, confirming the repo's own
    config keys were excluded from USER_KEYS without breaking real typo detection."""
    base = _make_base_with_repo_config(tmp_path)
    with open(base / "nextflow.config", "a") as fh:
        fh.write('params.kmer_lenght = 31\n')
    rc, out = run(base)
    assert rc != 0 and "did you mean 'kmer_length'" in out, out[-3000:]
