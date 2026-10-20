"""inventory.py: the files each workflow step creates under analysis_dir, and how the steps
depend on each other. Used by preflight.py for the overwrite check (a step that will run must
not find its own files from an earlier run) and to find outputs made out of date by a rerun.

One analysis per (analysis_dir, output_prefix, kmer_type, kmer_length): the per-genome k-mer
counts and the step folders are named by kmer type and length only, so two analyses with the
same type and length need separate analysis_dir folders (precomputed_dir shares steps 1-3 and 5).

Patterns are relative to analysis_dir. {P} is <output_prefix>_<kmer_type><kmer_length>,
{O} the output prefix, {T} the kmer type, {K} the kmer length, {TK} <kmer_type>kmer<kmer_length>.
"**" matches anything (including "/"), "*" anything within one path component."""
import os
import re

# Steps and the steps whose outputs they read (kmer_pipeline.nf)
DEPENDS = {1: [], 2: [1], 3: [1, 2], 4: [3], 5: [1, 2], 6: [3, 4, 5], 7: [4, 5, 6]}

STEP_FILES = {
    1: ["{TK}/**", "translated_contigs/**", "{P}_kmers_filepaths.txt"],
    2: ["{P}.kmermerge.txt.gz", "{O}.{T}{K}.j.*"],
    3: ["{TK}_patternbatches/**", "{P}.patternmerge.patternKey.txt.gz", "{P}.patternmerge.patternKeySize.txt",
        "{P}.patternmerge.patternIndex.txt.gz", "{P}.kinshipmerge.kinship.txt.gz", "{P}.kinshipmerge.kinshipWeight.txt"],
    # The presence counts depend on the phenotype: made by step 4 (by step 3 before Phase 4, N4)
    4: ["{TK}_gemma/**", "{P}.patternmerge.presenceCount.txt.gz", "{P}.patternmerge.presenceCount.txt.gz.tmp"],
    5: ["{TK}_kmergenealign/**", "{P}.*.kmeralignmerge.txt.gz", "{P}.*.kmeralignmerge.count.txt.gz"],
    6: ["{TK}_kmergenealign_figures/**", "{P}.summary.json"],
    7: ["{P}.report.html", "{P}.report_*.html", "report.css", "report.js"],
}

# Written by the workflow itself on every run, not by a step
RUN_FILES = ["{P}.container_id_file.txt", "{P}.analysis_file.txt", "{P}.run_manifest.json", "log.{P}/**",
             "work.{P}", "work.{P}/**"]


def _regex(pattern, prefix, kmer_type, kmer_length):
    names = {"P": f"{prefix}_{kmer_type}{kmer_length}", "O": prefix, "T": kmer_type, "K": str(kmer_length),
             "TK": f"{kmer_type}kmer{kmer_length}"}
    out, i = [], 0
    while i < len(pattern):
        if pattern.startswith("/**", i) and i + 3 == len(pattern):  # the folder itself (or a link to one) too
            out.append("(/.*)?")
            i += 3
        elif pattern.startswith("**", i):
            out.append(".*")
            i += 2
        elif pattern[i] == "*":
            out.append("[^/]*")
            i += 1
        elif pattern[i] == "{":
            j = pattern.index("}", i)
            out.append(re.escape(names[pattern[i + 1:j]]))
            i = j + 1
        else:
            out.append(re.escape(pattern[i]))
            i += 1
    return re.compile("".join(out) + r"\Z")


class Inventory:
    def __init__(self, prefix, kmer_type, kmer_length):
        self.steps = {s: [_regex(p, prefix, kmer_type, kmer_length) for p in pats] for s, pats in STEP_FILES.items()}
        self.run = [_regex(p, prefix, kmer_type, kmer_length) for p in RUN_FILES]

    def owner(self, relpath):
        """The step that creates relpath, "run" for workflow files, or None (not this analysis's)."""
        for s, regexes in self.steps.items():
            if any(r.match(relpath) for r in regexes):
                return s
        if any(r.match(relpath) for r in self.run):
            return "run"
        return None

    def files(self, analysis_dir):
        """{step: [relative paths]} of the files (and symbolic links) present under analysis_dir;
        "run" and None as in owner()."""
        found = {}
        if not os.path.isdir(analysis_dir):
            return found
        for root, dirs, names in os.walk(analysis_dir):
            rel_root = os.path.relpath(root, analysis_dir)
            for d in list(dirs):  # symbolic links to folders (work.*) are listed, not followed
                if os.path.islink(os.path.join(root, d)):
                    names.append(d)
                    dirs.remove(d)
            for n in names:
                rel = n if rel_root == "." else rel_root + "/" + n
                found.setdefault(self.owner(rel), []).append(rel)
        for v in found.values():
            v.sort()
        return found


def downstream(steps):
    """Every step that depends, directly or not, on any of steps."""
    out = set()
    changed = True
    while changed:
        changed = False
        for s, deps in DEPENDS.items():
            if s not in out and any(d in steps or d in out for d in deps):
                out.add(s)
                changed = True
    return out - set(steps)
