#!/usr/bin/env python3
"""preflight.py: checks run by kmer_pipeline.nf before the workflow writes anything.

- Parameters: unknown names are reported; names close to a known one (a typo, or kebab-case
  such as --kmer-min-count, which Nextflow stores as kmerMinCount) are errors.
- Reruns: a step that will run must not find its own files from an earlier run in analysis_dir
  (they would be overwritten, or mixed with the new results). With overwrite = true they are
  deleted first, together with the outputs of later skipped steps that the rerun makes out of
  date; a step that reads such out-of-date outputs is an error. Folders are checked against
  inventory.py.
- Reuse: when skipped steps' outputs will be used, the genomes and their order must be those
  the step-1 outputs were made with.
- -resume continues only the run recorded in the analysis's run manifest, with unchanged
  result-affecting parameters and input files.

Prints one JSON object on stdout: {"errors": [...], "warnings": [...], "deleted": [...]}; the
workflow stops if there are errors. Exit status is 0 unless the checks themselves fail.
With --finish, records the end of the run in the manifest instead."""
import argparse
import difflib
import hashlib
import json
import os
import re
import socket
import sys
import time

import inventory
import rcompat

# Every parameter kmer_pipeline.nf reads (including implied ones the manual allows to override)
KNOWN_PARAMS = {
    "analysis_dir", "analysis_file", "annotateGeneFile", "base_dir", "blastident", "bowtie_parameters",
    "container_analysis_dir", "container_analysis_file", "container_args", "container_cmd",
    "container_covariate_file", "container_file", "container_id_file", "container_logdir", "container_mount",
    "container_ref_fa", "container_ref_gb", "container_script_dir", "container_software_file", "container_type",
    "container_user_id_file",
    "covariate_file", "default_script_dir", "default_software_file", "gene_lookup_file", "id_file",
    "kmerFilePrefix", "kmergenecombination", "kmer_length", "kmer_min_count", "kmer_type", "logdir", "maxp",
    "merge_wait_minutes", "min_count", "minor_allele_threshold", "n", "ntopgenes", "nucmerident", "output_prefix",
    "overwrite", "override_signif", "p", "p5", "plot_min_genomes", "ref_fa", "ref_gb", "ref_name",
    "samtools_filter", "skip1", "skip2", "skip3", "skip4", "skip5", "skip6", "skip7", "software_file", "workdir",
}

MANIFEST_VERSION = 1
SHOW = 20  # files listed per message


def near_miss(name):
    """The known parameter name is close to (a likely misspelling of), or None."""
    snake = re.sub(r"(?<=[a-z0-9])([A-Z])", r"_\1", name).replace("-", "_").lower()
    if snake in KNOWN_PARAMS:
        return snake
    close = difflib.get_close_matches(snake, KNOWN_PARAMS, n=1, cutoff=0.8)
    return close[0] if close else None


def check_params(user_params, errors, warnings):
    for name in sorted(user_params):
        if name in KNOWN_PARAMS:
            continue
        guess = near_miss(name)
        if guess:
            errors.append(f"unknown parameter '{name}': did you mean '{guess}'?")
        else:
            warnings.append(f"parameter '{name}' is not used by the workflow")


def md5(path):
    if not path or not os.path.isfile(path):
        return None
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def read_manifest(path):
    try:
        with open(path) as fh:
            return json.load(fh)
    except (OSError, ValueError):
        return None


def write_manifest(path, manifest):
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(manifest, fh, indent=1, sort_keys=True)
        fh.write("\n")
    os.replace(tmp, path)


def summarise(paths):
    shown = ", ".join(paths[:SHOW])
    return shown + (f" and {len(paths) - SHOW} more" if len(paths) > SHOW else "")


def genomes_of_step1(filepaths_file, kmer_length):
    """Genome IDs, in order, that step 1 counted (from <prefix>_<type><k>_kmers_filepaths.txt)."""
    suffix = f".kmer{kmer_length}.txt.gz"
    with open(filepaths_file) as fh:
        return [os.path.basename(line.rstrip("\n"))[:-len(suffix)] for line in fh if line.strip()]


def genomes_of_id_file(id_file):
    """Genome IDs as the scripts read them (R's read.table type conversion, as in countkmers)."""
    table = rcompat.r_read_table(id_file, header=True, sep="\t")
    return [rcompat.r_as_character(i) for i in table["id"]]


def live_run(manifest):
    """A description of the run the manifest records as still running, or None."""
    if not manifest or manifest.get("status") != "running":
        return None
    host, pid = manifest.get("host"), manifest.get("pid")
    if host == socket.gethostname() and pid:
        try:
            os.kill(int(pid), 0)
        except (ProcessLookupError, ValueError):
            return None
        except PermissionError:
            pass
    return f"started {manifest.get('started')} on {host} (process {pid})"


def delete(analysis_dir, relpaths, errors):
    """Two phases: check every file can be removed, then remove them all (symbolic links are
    removed, never followed), then any folders left empty. Returns the deleted paths."""
    blocked = [p for p in relpaths
               if not os.access(os.path.dirname(os.path.join(analysis_dir, p)) or analysis_dir, os.W_OK)]
    if blocked:
        errors.append("overwrite = true, but these files cannot be removed (permissions): " + summarise(blocked))
        return []
    folders = set()
    for p in relpaths:
        full = os.path.join(analysis_dir, p)
        os.unlink(full)
        d = os.path.dirname(p)
        while d:
            folders.add(d)
            d = os.path.dirname(d)
    for d in sorted(folders, key=lambda x: -x.count("/")):
        full = os.path.join(analysis_dir, d)
        if os.path.isdir(full) and not os.path.islink(full) and not os.listdir(full):
            os.rmdir(full)
    return relpaths


def check(args):
    errors, warnings, deleted = [], [], []
    run_steps = sorted(int(s) for s in args.run_steps.split(",") if s)
    skipped = [s for s in inventory.DEPENDS if s not in run_steps]
    params = json.loads(args.params_json)
    inputs = {k: md5(v) for k, v in json.loads(args.input_files).items()}
    prefix_full = f"{args.output_prefix}_{args.kmer_type}{args.kmer_length}"
    manifest_path = os.path.join(args.analysis_dir, prefix_full + ".run_manifest.json")
    manifest = read_manifest(manifest_path)

    check_params([p for p in args.user_params.split(",") if p], errors, warnings)
    resume = args.resume == "true"
    overwrite = args.overwrite == "true"

    if resume:
        if overwrite:
            errors.append("overwrite = true cannot be used with -resume (it would delete the outputs of the "
                          "tasks Nextflow reuses)")
        elif manifest is None:
            errors.append("-resume: there is no run to continue in " + args.analysis_dir)
        elif manifest.get("session") != args.session_id:
            errors.append("-resume continues the last run started from this launch folder (session "
                          f"{args.session_id}), but {args.analysis_dir} holds the outputs of another run "
                          f"(session {manifest.get('session')})")
        else:
            changed = sorted(k for k in set(params) | set(manifest.get("params", {}))
                             if params.get(k) != manifest.get("params", {}).get(k))
            changed += sorted(k for k in set(inputs) | set(manifest.get("inputs", {}))
                              if inputs.get(k) != manifest.get("inputs", {}).get(k))
            if changed:
                errors.append("-resume only continues an interrupted run with the same inputs; changed since "
                              "that run: " + ", ".join(changed) + ". Run without -resume (with overwrite = true "
                              "to replace the earlier outputs)")
        return finish_check(errors, warnings, deleted, manifest_path, manifest, args, params, inputs, run_steps,
                            write=not errors, resume=True)

    inv = inventory.Inventory(args.output_prefix, args.kmer_type, args.kmer_length)
    present = inv.files(args.analysis_dir)
    stale = inventory.downstream(run_steps) & set(skipped)
    for s in run_steps:
        for d in inventory.DEPENDS[s]:
            if d in stale:
                cause = [r for r in run_steps if d in inventory.downstream([r])]
                errors.append(f"step {s} would read the outputs of step {d}, which are out of date because step "
                              f"{min(cause)} runs: run step {d} too (skip{d} = false)")

    targets = [p for s in run_steps + sorted(stale) for p in present.get(s, [])]
    if targets:
        live = live_run(manifest)
        if live:
            errors.append(f"another run of this analysis may still be running ({live}); if it is not, delete "
                          f"{manifest_path} and try again")
        elif not overwrite:
            per_step = ", ".join(f"step {s}: {len(present.get(s, []))}" for s in run_steps + sorted(stale)
                                 if present.get(s))
            errors.append("outputs of an earlier run are in " + args.analysis_dir + " (" + per_step + " files, e.g. "
                          + summarise(targets) + "). To replace them set overwrite = true; to keep them, use "
                          "another analysis_dir (precomputed_dir can reuse its steps 1-3 and 5)")

    # Reuse: genomes and their order as step 1 counted them
    filepaths = os.path.join(args.analysis_dir, prefix_full + "_kmers_filepaths.txt")
    if 1 not in run_steps and run_steps and os.path.isfile(filepaths) and os.path.isfile(args.id_file):
        if genomes_of_step1(filepaths, args.kmer_length) != genomes_of_id_file(args.id_file):
            errors.append("the genomes in id_file, or their order, differ from those the existing step-1 outputs "
                          f"were made with ({filepaths}): the earlier steps must be rerun (skip1 = false, "
                          "overwrite = true)")
    if manifest is None and skipped and run_steps and any(present.get(s) for s in skipped):
        warnings.append("this analysis_dir has no run manifest (made by an earlier release), so the reference "
                        "and k-mer settings of the outputs being reused cannot be checked")
    elif manifest is not None:
        prov = manifest.get("steps", {})
        for s in skipped:
            p = prov.get(str(s))
            if p is None or s in stale or not any(s in inventory.DEPENDS[r] for r in run_steps):
                continue
            keys = {3: ["kmer_min_count"], 5: ["nucmerident", "ref_fa", "ref_gb"]}.get(s, [])
            diff = [k for k in keys if p.get("params", {}).get(k, p.get("inputs", {}).get(k)) !=
                    params.get(k, inputs.get(k))]
            if diff:
                errors.append(f"the step-{s} outputs being reused were made with different "
                              + ", ".join(diff) + f": rerun step {s} (skip{s} = false, overwrite = true)")

    if not errors and targets and overwrite:
        deleted = delete(args.analysis_dir, targets, errors)
        if deleted:
            listing = os.path.join(args.analysis_dir, "log." + prefix_full,
                                   time.strftime("overwrite_deleted_%Y%m%d-%H%M%S.txt"))
            os.makedirs(os.path.dirname(listing), exist_ok=True)
            with open(listing, "w") as fh:
                fh.write("".join(p + "\n" for p in deleted))
            warnings.append(f"overwrite = true: deleted {len(deleted)} files of an earlier run from "
                            f"{args.analysis_dir}, listed in {listing}")
    return finish_check(errors, warnings, deleted, manifest_path, manifest, args, params, inputs, run_steps,
                        write=not errors, resume=False, stale=stale)


def finish_check(errors, warnings, deleted, manifest_path, manifest, args, params, inputs, run_steps, write,
                 resume, stale=()):
    if write:
        os.makedirs(args.analysis_dir, exist_ok=True)
        steps = dict((manifest or {}).get("steps", {}))
        if not resume:
            for s in list(steps):
                if int(s) in run_steps or int(s) in stale:
                    del steps[s]
            for s in run_steps:
                steps[str(s)] = {"session": args.session_id, "params": params, "inputs": inputs}
        write_manifest(manifest_path, {
            "version": MANIFEST_VERSION, "session": args.session_id, "status": "running",
            "started": time.strftime("%Y-%m-%dT%H:%M:%S%z"), "host": socket.gethostname(), "pid": args.pid,
            "params": params, "inputs": inputs, "steps": steps})
    return {"errors": errors, "warnings": warnings, "deleted": deleted}


def finish(args):
    prefix_full = f"{args.output_prefix}_{args.kmer_type}{args.kmer_length}"
    path = os.path.join(args.analysis_dir, prefix_full + ".run_manifest.json")
    manifest = read_manifest(path)
    if manifest is not None and manifest.get("session") == args.session_id:
        manifest["status"] = args.finish
        manifest["finished"] = time.strftime("%Y-%m-%dT%H:%M:%S%z")
        write_manifest(path, manifest)
    return {"errors": [], "warnings": [], "deleted": []}


def main():
    parser = argparse.ArgumentParser(description="preflight.py checks before kmer_pipeline.nf runs",
                                     allow_abbrev=False)
    parser.add_argument("--analysis-dir", required=True)
    parser.add_argument("--output-prefix", required=True)
    parser.add_argument("--kmer-type", required=True)
    parser.add_argument("--kmer-length", required=True)
    parser.add_argument("--session-id", required=True)
    parser.add_argument("--run-steps", default="", help="steps that will run, e.g. 1,2,3")
    parser.add_argument("--overwrite", default="false")
    parser.add_argument("--resume", default="false")
    parser.add_argument("--pid", default="")
    parser.add_argument("--user-params", default="", help="names of the parameters the user set, comma-separated")
    parser.add_argument("--params-json", default="{}", help="result-affecting parameters (JSON)")
    parser.add_argument("--input-files", default="{}", help="input files to checksum (JSON name: path)")
    parser.add_argument("--id-file", default="")
    parser.add_argument("--finish", default=None, help="record the end of the run: finished or failed")
    args = parser.parse_args()
    result = finish(args) if args.finish else check(args)
    sys.stdout.write(json.dumps(result) + "\n")


if __name__ == "__main__":
    main()
