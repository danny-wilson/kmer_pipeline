# End-to-end test harness

Two layers of tests cover the pipeline:

- **`tests/`** holds self-checking assertions: no goldens, no stored outputs to compare against. Most run
  with the image's own Python/R directly (no container at all, "launcher mode" — `kmer_pipeline.nf` run by
  path, with the image's own `bin` on `PATH`); a few (`test_nextflow_config.py`) exercise a real
  `-profile apptainer --container` run.
- **`tests/e2e/`** (this directory) holds comparisons against golden outputs, with the toolchain image as
  the actual Nextflow container (`-profile apptainer --container <SIF>`, the repo's own `nextflow.config`).
  These need a local configuration (see below) and are skipped by a plain `pytest tests` run unless that
  configuration is present.

## Prerequisites

- [Apptainer](https://apptainer.org/) (or Singularity).
- Nextflow 22.04.5 (the version the image and every golden was made with; `local.conf`'s
  `KMER_E2E_NEXTFLOW_VERSION` is checked against it).
- Java 11 or newer to run that Nextflow version — if the host's default Java is older, set
  `KMER_E2E_PRE` in `local.conf` to whatever makes a newer one available before
  Nextflow runs.
- A host Python 3.8+ for `compare.py`, `golden_hashes.py`, `e2e_config.py`, `make_split_reference.py` and
  `make_stress_inputs.py`: these use only the standard library plus this project's own `e2e_config`, so
  they do not need the image. `split_compare.py` is the exception — it imports pipeline modules
  (`reference.py`) and the image's own packages, so it must run under `apptainer exec` with `PYTHONPATH`
  set to a staged commit's directory (see its own usage comment). `test_golden_hashes.py` (under `tests/`,
  not here) needs `pytest`, so it is always run the same way the rest of `tests/` is, under
  `apptainer exec`.

## `local.conf`

Copy `local.conf.example` to `local.conf` (gitignored) and fill in the required keys, or set the same
names as environment variables, which take precedence over the file:

- `KMER_E2E_ROOT` — scratch root holding `goldens/`, `runs/`, `staging/`, `images/`, `nxf_home/`, `inputs/`.
- `KMER_E2E_SIF` — the toolchain image (code-full or code-free; `make_reference.sh`'s `example/`
  byte-identity check against the image's own bundled copy is skipped gracefully against a code-free one,
  which has no `/usr/share/kmer_pipeline`).
- `KMER_E2E_CODEFREE_SIF` (optional) — a code-stripped derivative of `KMER_E2E_SIF`
  (`build_codefree_derivative.sh`), for checks that must not be able to silently fall back to a
  baked-in script copy instead of the checkout's own `bin/` (B9/B10's Step 1 validation; B8's
  `tests/` invocation).
- `KMER_E2E_IMAGE_COMMIT` — the commit the image's C++ tools/Dockerfile/dependencies were built from.
  `stage.sh` warns (not fails) if a staged commit's `C++`/`Makefile`/`Dockerfile` have since diverged from
  this commit, since that drift is exactly what this harness cannot detect (see below).
- `KMER_E2E_NEXTFLOW`, `KMER_E2E_NEXTFLOW_VERSION`.
- `KMER_E2E_PRE` (optional) — a shell snippet `config.sh` evaluates before checking the Nextflow version,
  e.g. to put a newer Java on `PATH`. **Quote it if it contains spaces** — `local.conf` is sourced by bash.
- `KMER_E2E_LOCAL_TMP` (default `/tmp`) — fast local disk for unpacking a SIF build.
- `KMER_E2E_NICE` (default `10`).
- `KMER_E2E_BASE_SIF` (optional, `check_image.sh` only) — a previous production image to diff
  pip/conda/R package lists against. Skipped, not failed, if unset.

## How to

- **Stage a commit:** `stage.sh COMMIT` → `$KMER_E2E_ROOT/staging/<sha>`, a plain `git archive` of that
  commit's own tree (`bin/`, `kmer_pipeline.nf`, `nextflow.config`, `conf/`, `example/`) — nothing
  flattened or symlinked: `bin/` is native, so Nextflow's own `projectDir`-relative `PATH` entry finds it,
  and the committed `conf/software_location.txt`'s `toolpath` finds the compiled C++ tools in the image
  itself, unchanged by which commit is staged.
- **Make a golden run:** `make_reference.sh --label L --sha SHA --run NAME --kmer-type T --kmer-length K
  --maxp M --staging $KMER_E2E_ROOT/staging/SHA` → `$KMER_E2E_ROOT/goldens/L-sha7/maxpM/NAME/stage{1..7}`,
  one `nextflow run $STAGING/kmer_pipeline.nf -profile apptainer --container <SIF>` invocation per stage
  (`skip1..skip7`) with a snapshot taken after each. `example/` comes from the staged commit's own tree.
- **Make an all-stages golden run:** `make_allstages_reference.sh` — the same invocation with no `skip*`
  flags at all (the pipeline's own default, non-decomposed DAG), one `stage7` snapshot only.
- **Compare two commits:** `compare_commits.sh AFTER_LABEL BEFORE_LABEL BEFORE_SHA AFTER_SHA` stages
  `AFTER_SHA` if needed, runs both example datasets at `maxp 2` under `AFTER_LABEL`, and prints
  `compare.py`'s stage-by-stage diff against the `BEFORE_LABEL` goldens. A reused `AFTER_LABEL` output
  directory is only trusted if its `.complete` marker's recorded image sha256 matches the current image.
- **Verify hashes:** `golden_hashes.py verify golden_hashes.tsv --run NAME DIR [--run NAME DIR ...]
  [--prefix NAME=BASE_DIR ...]` (`--prefix` is that run's own `base_dir`, for the `path` class;
  `--visual-advisory` downgrades a PNG hash mismatch to a warning; a listed-but-missing PNG still fails).
- **Regenerate hashes:** a deliberate, separate commit — `golden_hashes.py generate --run NAME DIR
  [--run NAME DIR ...] [--prefix NAME=BASE_DIR ...] > tests/e2e/golden_hashes.tsv`, reviewed like any
  other diff, with a commit message saying why outputs changed.
- **Run step tests:** `run_step_tests.sh STAGING_DIR STAGE REF_LABEL` (wraps `step_test.sh`, which reruns
  one stage seeded from a reference run's previous-stage snapshot and compares the result).
- **Run resume rules:** `resume_rules_test.sh SHA` — the seven `-resume`/`resume=true` scenarios this
  project's fixes depend on.
- **Run the precomputed-dir test:** `precomputed_test.sh NAME SHA PRECOMPUTED_STAGE5 PHENO_ID_FILE
  FULL_STAGE7 [COVARIATES]`.

## What this harness does not cover, and why

- **No Slurm-interrupt test.** The private harness's version used `squeue`/`scancel`; a local-process-kill
  equivalent would be a future addition.
- **`dockerfile2def.py` is not a general Dockerfile-to-Apptainer transpiler** — it understands only `FROM`
  (by digest), `ARG`, `LABEL`, `USER`, `WORKDIR`, `RUN`, `COPY . <dir>` and `ENV`, and exits on anything
  else. It also inherits a genuine Apptainer behaviour difference from Docker: Apptainer's `%labels` does
  **not** override a label key the base image already sets, so a locally fakeroot-built image's `version`
  label may still show the base image's value rather than the one just built — checked directly (C4): the
  actual script and image content were otherwise correct (installed scripts, environment, imports, every
  `--help` all matched the target commit exactly), only the `version` **label** was stale.
- **`check_image.sh`'s package/version lists are a point-in-time snapshot** (the `$ADDED` pip list, the
  exact `imports`/`environment` strings) that must be maintained, or deliberately regenerated, whenever the
  image changes; the script cannot discover drift on its own. Its pip/conda/R baseline comparison is
  skipped entirely unless `KMER_E2E_BASE_SIF` is set, since "a previous production image" is a locally-kept
  comparison asset, not something every clone of this harness has.
- **PNG hashes change with any change to the R graphics stack**, not only a change to this pipeline's own
  code — `golden_hashes.py --visual-advisory` exists for exactly this.
- **A regression confined to the C++ tools, the Dockerfile, or a Python/R dependency pin is not
  detected.** The toolchain (compiled binaries, Python/R packages) comes from the image, not from the
  staged `git archive` (only `bin/`, `nextflow.config`, `conf/`, `example/` do). `stage.sh`'s drift warning
  against `KMER_E2E_IMAGE_COMMIT` is a partial mitigation, not a fix. `example/` itself **is** now taken
  from the staged commit (B9), with `make_reference.sh` asserting its data files (`*.fa.gz`/`*.fasta`/
  `*.gb`) are byte-identical to a code-full image's own bundled copy, to surface drift there rather than
  silently diverging once the published image eventually stops carrying that copy at all.
- **The shipped example data cannot exercise two things even with full goldens:** `annotateGeneFile`
  reaching steps 6/7, and `min_contig_length > 0` actually dropping contigs inside Nextflow — `tests/`' own
  `test_workflow_run.py` (B12, B13) are the only guard for those. Low-complexity k-mer markup has unit-test
  coverage only (no example k-mer is low-complexity).
- **The two SIFs this project has used (`kmer_pipeline_2026-10-06.sif` and `kmer_pipeline_test-98e242a.sif`)
  are one image at one public registry digest, not two** (checked directly: identical sha256). This does
  **not** extend to every historical reference run cited in this project's own records — some predate this
  image and were made with a different one (the original `2022-10-26` image), on different hosts.
- **Sub-project 3's container refactor removed the fixed mount point** (`/home/jovyan`) in favour of
  identity host=container paths (`base_dir`/`projectDir` bound at their own, real spelling). 5
  path-bearing files per dataset (`manifest.tsv`'s `path` class) now embed the run's own `base_dir`, which
  differs between any two runs — `compare.py --prefix-a/--prefix-b` and `golden_hashes.py --prefix
  NAME=BASE_DIR` normalize that one known prefix to `@BASE_DIR@` before comparing/hashing. If `base_dir`
  is itself reached through a symlink (routine on a shared HPC filesystem), the prefix that is actually
  bound into the container — and that these path-bearing files must therefore be given relative to — is
  the *canonical* (symlink-resolved) spelling, not whatever raw spelling the run happened to be launched
  with; `kmer_pipeline.nf`'s `absPath()` handles this rewrite internally (B9).
