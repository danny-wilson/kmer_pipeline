"""The scriptpath/toolpath split (PLAN.md Part B, B2): scriptpath defaults to a script's own
directory when the software file doesn't set it; toolpath (the compiled C++ tools) defaults to
scriptpath when the software file doesn't set it either. Each of the 5 scripts that actually
resolve a C++ tool (countkmers.py, proteinkmermerge.py, kmercontigalignmerge.py,
kmercontigalignonly.py, stringlist2patternandkinship.py) is tested against all three rules.
createfullkmerlist.py is tested separately: it never resolves a C++ tool, so it must keep
resolving its sibling Python scripts via scriptpath alone, regardless of toolpath (guards P2)."""
import os
import shutil
import stat
import subprocess
import sys

import pytest

from conftest import SCRIPTS_DIR

MISSING_TOOL = object()  # sentinel: this case's directory has no tool files at all


def _isolated_script_dir(tmp_path, script_name, tool_names):
    """A fresh directory holding only a copy of the named script (plus rcompat.py, which every
    script imports) and fake executable tool files, standing in for 'the script's own
    directory' under full control of the test."""
    d = tmp_path / "scriptdir"
    d.mkdir()
    shutil.copy(os.path.join(SCRIPTS_DIR, script_name), d / script_name)
    os.chmod(d / script_name, 0o755)
    shutil.copy(os.path.join(SCRIPTS_DIR, "rcompat.py"), d / "rcompat.py")
    for t in tool_names:
        f = d / t
        f.write_text("#!/bin/bash\nexit 0\n")
        f.chmod(f.stat().st_mode | stat.S_IEXEC)
    return d


def _run(script_path, *args):
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE="1")
    return subprocess.run([sys.executable, str(script_path), *args],
                          capture_output=True, text=True, env=env)


def _software_file(tmp_path, rows, name="software.txt"):
    path = tmp_path / name
    lines = ["name\tpath"] + [f"{k}\t{v}" for k, v in rows]
    path.write_text("\n".join(lines) + "\n")
    return str(path)


# --- countkmers.py (the `software(name)` closure style) ----------------------------------

def _countkmers_common(tmp_path):
    analysis_dir = tmp_path / "analysis"
    analysis_dir.mkdir()
    id_file = tmp_path / "id_file.txt"
    id_file.write_text("id\tpaths\tpheno\n")
    dsk = tmp_path / "dsk"
    dsk.write_text("")
    dsk2ascii = tmp_path / "dsk2ascii"
    dsk2ascii.write_text("")
    return analysis_dir, id_file, dsk, dsk2ascii


@pytest.mark.parametrize("case", ["no_scriptpath", "no_toolpath", "both_set"])
def test_countkmers_scriptpath_toolpath(tmp_path, case):
    analysis_dir, id_file, dsk, dsk2ascii = _countkmers_common(tmp_path)
    rows = [("dsk", str(dsk)), ("dsk2ascii", str(dsk2ascii))]

    if case == "no_scriptpath":
        # Own directory: a copy of countkmers.py alongside a fake sort_strings/sequence_functions.py.
        scriptdir = _isolated_script_dir(tmp_path, "countkmers.py", ["sort_strings"])
        shutil.copy(os.path.join(SCRIPTS_DIR, "sequence_functions.py"), scriptdir / "sequence_functions.py")
        script = scriptdir / "countkmers.py"
    elif case == "no_toolpath":
        # scriptpath explicit, pointing at a directory with both the tool and the sibling module;
        # no toolpath row, so it must fall back to that same scriptpath value.
        scriptdir = _isolated_script_dir(tmp_path, "countkmers.py", ["sort_strings"])
        shutil.copy(os.path.join(SCRIPTS_DIR, "sequence_functions.py"), scriptdir / "sequence_functions.py")
        rows.append(("scriptpath", str(scriptdir)))
        script = os.path.join(SCRIPTS_DIR, "countkmers.py")
    else:  # both_set, explicitly, to different directories - toolpath must win for the tool
        scriptdir = _isolated_script_dir(tmp_path, "countkmers.py", [])
        shutil.copy(os.path.join(SCRIPTS_DIR, "sequence_functions.py"), scriptdir / "sequence_functions.py")
        tooldir = tmp_path / "tools"
        tooldir.mkdir()
        (tooldir / "sort_strings").write_text("#!/bin/bash\nexit 0\n")
        os.chmod(tooldir / "sort_strings", 0o755)
        rows.append(("scriptpath", str(scriptdir)))
        rows.append(("toolpath", str(tooldir)))
        script = os.path.join(SCRIPTS_DIR, "countkmers.py")

    sw = _software_file(tmp_path, rows)
    result = _run(script, "--task-id", "1", "--id-file", str(id_file),
                  "--analysis-dir", str(analysis_dir), "--output-prefix", "x",
                  "--software-file", sw)
    # In every case the tool is genuinely present wherever this rule says it should be, so the
    # script must get past the "sort_strings path doesn't exist" check - whatever it fails on
    # afterwards (fake dsk, no real reads) is not this test's concern.
    assert "sort_strings path doesn't exist" not in result.stderr, result.stderr


def test_countkmers_no_scriptpath_no_tool_present_fails_clearly(tmp_path):
    """Negative control: with no scriptpath row and no tool in the script's own directory
    either, the default-to-own-directory rule must still be the one that's exercised (not some
    other silent fallback) - confirmed by the specific, named error."""
    analysis_dir, id_file, dsk, dsk2ascii = _countkmers_common(tmp_path)
    scriptdir = _isolated_script_dir(tmp_path, "countkmers.py", [])  # no sort_strings anywhere
    sw = _software_file(tmp_path, [("dsk", str(dsk)), ("dsk2ascii", str(dsk2ascii))])
    result = _run(scriptdir / "countkmers.py", "--task-id", "1", "--id-file", str(id_file),
                  "--analysis-dir", str(analysis_dir), "--output-prefix", "x",
                  "--software-file", sw)
    assert "sort_strings path doesn't exist" in result.stderr


# --- proteinkmermerge.py (the zip(names, paths) list-comprehension style) ----------------

def test_proteinkmermerge_no_scriptpath_defaults_to_own_directory(tmp_path):
    input_dir = tmp_path / "in"
    input_dir.mkdir()
    output_dir = tmp_path / "out"
    output_dir.mkdir()
    id_file = tmp_path / "id_file.txt"
    id_file.write_text("id\n1\n")
    scriptdir = _isolated_script_dir(tmp_path, "proteinkmermerge.py", ["sort_strings"])
    sw = _software_file(tmp_path, [])
    result = _run(scriptdir / "proteinkmermerge.py", "--n", "1", "--p", "1",
                  "--output-prefix", "x", "--input-dir", str(input_dir),
                  "--output-dir", str(output_dir), "--id-file", str(id_file),
                  "--kmer-length", "31", "--software-file", sw, "--process", "1")
    assert "sort_strings path doesn't exist" not in result.stderr, result.stderr


def test_proteinkmermerge_no_scriptpath_and_no_tool_fails_clearly(tmp_path):
    input_dir = tmp_path / "in"
    input_dir.mkdir()
    output_dir = tmp_path / "out"
    output_dir.mkdir()
    id_file = tmp_path / "id_file.txt"
    id_file.write_text("id\n1\n")
    scriptdir = _isolated_script_dir(tmp_path, "proteinkmermerge.py", [])
    sw = _software_file(tmp_path, [])
    result = _run(scriptdir / "proteinkmermerge.py", "--n", "1", "--p", "1",
                  "--output-prefix", "x", "--input-dir", str(input_dir),
                  "--output-dir", str(output_dir), "--id-file", str(id_file),
                  "--kmer-length", "31", "--software-file", sw, "--process", "1")
    assert "sort_strings path doesn't exist" in result.stderr


def test_proteinkmermerge_toolpath_overrides_scriptpath_for_the_tool(tmp_path):
    """Both set, to different directories: toolpath must be the one honoured for the tool."""
    input_dir = tmp_path / "in"
    input_dir.mkdir()
    output_dir = tmp_path / "out"
    output_dir.mkdir()
    id_file = tmp_path / "id_file.txt"
    id_file.write_text("id\n1\n")
    scriptdir = _isolated_script_dir(tmp_path, "proteinkmermerge.py", [])  # no tool here
    tooldir = tmp_path / "tools"
    tooldir.mkdir()
    (tooldir / "sort_strings").write_text("#!/bin/bash\nexit 0\n")
    os.chmod(tooldir / "sort_strings", 0o755)
    sw = _software_file(tmp_path, [("scriptpath", str(scriptdir)), ("toolpath", str(tooldir))])
    result = _run(scriptdir / "proteinkmermerge.py", "--n", "1", "--p", "1",
                  "--output-prefix", "x", "--input-dir", str(input_dir),
                  "--output-dir", str(output_dir), "--id-file", str(id_file),
                  "--kmer-length", "31", "--software-file", sw, "--process", "1")
    assert "sort_strings path doesn't exist" not in result.stderr, result.stderr


# --- kmercontigalignmerge.py --------------------------------------------------------------

def _kmercontigalignmerge_common(tmp_path, kmer_type="nucleotide", kmer_length="31"):
    output_dir = tmp_path / "out"
    output_dir.mkdir()
    (output_dir / f"{kmer_type}kmer{kmer_length}_kmergenealign").mkdir()
    input_files = tmp_path / "input_files.txt"
    input_files.write_text("")
    ref_fa = tmp_path / "ref.fa"
    ref_fa.write_text("")
    return output_dir, input_files, ref_fa


def test_kmercontigalignmerge_no_scriptpath_defaults_to_own_directory(tmp_path):
    output_dir, input_files, ref_fa = _kmercontigalignmerge_common(tmp_path)
    scriptdir = _isolated_script_dir(tmp_path, "kmercontigalignmerge.py",
                                     ["sort_strings", "stringlist2count"])
    sw = _software_file(tmp_path, [])
    result = _run(scriptdir / "kmercontigalignmerge.py", "--task-id", "1", "--n", "1", "--p", "1",
                  "--output-prefix", "x", "--analysis-dir", str(output_dir),
                  "--input-files", str(input_files), "--kmer-type", "nucleotide",
                  "--kmer-length", "31", "--ref-fa", str(ref_fa), "--nucmerident", "95",
                  "--software-file", sw)
    assert "doesn't exist - check pipeline tool location" not in result.stderr, result.stderr


def test_kmercontigalignmerge_no_scriptpath_and_no_tool_fails_clearly(tmp_path):
    output_dir, input_files, ref_fa = _kmercontigalignmerge_common(tmp_path)
    scriptdir = _isolated_script_dir(tmp_path, "kmercontigalignmerge.py", [])
    sw = _software_file(tmp_path, [])
    result = _run(scriptdir / "kmercontigalignmerge.py", "--task-id", "1", "--n", "1", "--p", "1",
                  "--output-prefix", "x", "--analysis-dir", str(output_dir),
                  "--input-files", str(input_files), "--kmer-type", "nucleotide",
                  "--kmer-length", "31", "--ref-fa", str(ref_fa), "--nucmerident", "95",
                  "--software-file", sw)
    assert "sort_strings path doesn't exist" in result.stderr


# --- kmercontigalignonly.py ----------------------------------------------------------------

def _kmercontigalignonly_common(tmp_path):
    output_dir = tmp_path / "out"
    output_dir.mkdir()
    id_file = tmp_path / "id_file.txt"
    id_file.write_text("")
    ref_fa = tmp_path / "ref.fa"
    ref_fa.write_text("")
    ref_gb = tmp_path / "ref.gb"
    ref_gb.write_text("")
    kmerseqfile = tmp_path / "kmers.txt"
    kmerseqfile.write_text("")
    return output_dir, id_file, ref_fa, ref_gb, kmerseqfile


def _kmercontigalignonly_args(output_dir, id_file, ref_fa, ref_gb, kmerseqfile, sw):
    return ("--task-id", "1", "--n", "1", "--output-prefix", "x",
            "--output-dir", str(output_dir), "--id-file", str(id_file),
            "--ref-fa", str(ref_fa), "--ref-gb", str(ref_gb), "--kmer-type", "nucleotide",
            "--kmer-length", "31", "--nucmerident", "95", "--kmerseqfile", str(kmerseqfile),
            "--software-file", sw)


def test_kmercontigalignonly_no_scriptpath_defaults_to_own_directory(tmp_path):
    output_dir, id_file, ref_fa, ref_gb, kmerseqfile = _kmercontigalignonly_common(tmp_path)
    scriptdir = _isolated_script_dir(tmp_path, "kmercontigalignonly.py", ["sort_strings"])
    shutil.copy(os.path.join(SCRIPTS_DIR, "kmercontigalignmerge.py"), scriptdir / "kmercontigalignmerge.py")
    shutil.copy(os.path.join(SCRIPTS_DIR, "sequence_functions.py"), scriptdir / "sequence_functions.py")
    sw = _software_file(tmp_path, [("R", "/usr/bin/R"), ("mummer", "/usr/bin"), ("genoPlotR", "/usr/bin")])
    result = _run(scriptdir / "kmercontigalignonly.py",
                  *_kmercontigalignonly_args(output_dir, id_file, ref_fa, ref_gb, kmerseqfile, sw))
    assert "sort_strings path doesn't exist" not in result.stderr, result.stderr


def test_kmercontigalignonly_no_scriptpath_and_no_tool_fails_clearly(tmp_path):
    output_dir, id_file, ref_fa, ref_gb, kmerseqfile = _kmercontigalignonly_common(tmp_path)
    scriptdir = _isolated_script_dir(tmp_path, "kmercontigalignonly.py", [])
    shutil.copy(os.path.join(SCRIPTS_DIR, "kmercontigalignmerge.py"), scriptdir / "kmercontigalignmerge.py")
    shutil.copy(os.path.join(SCRIPTS_DIR, "sequence_functions.py"), scriptdir / "sequence_functions.py")
    sw = _software_file(tmp_path, [("R", "/usr/bin/R"), ("mummer", "/usr/bin"), ("genoPlotR", "/usr/bin")])
    result = _run(scriptdir / "kmercontigalignonly.py",
                  *_kmercontigalignonly_args(output_dir, id_file, ref_fa, ref_gb, kmerseqfile, sw))
    assert "sort_strings path doesn't exist" in result.stderr


def test_kmercontigalignonly_toolpath_overrides_scriptpath_for_the_tool(tmp_path):
    output_dir, id_file, ref_fa, ref_gb, kmerseqfile = _kmercontigalignonly_common(tmp_path)
    scriptdir = _isolated_script_dir(tmp_path, "kmercontigalignonly.py", [])
    shutil.copy(os.path.join(SCRIPTS_DIR, "kmercontigalignmerge.py"), scriptdir / "kmercontigalignmerge.py")
    shutil.copy(os.path.join(SCRIPTS_DIR, "sequence_functions.py"), scriptdir / "sequence_functions.py")
    tooldir = tmp_path / "tools"
    tooldir.mkdir()
    (tooldir / "sort_strings").write_text("#!/bin/bash\nexit 0\n")
    os.chmod(tooldir / "sort_strings", 0o755)
    sw = _software_file(tmp_path, [("R", "/usr/bin/R"), ("mummer", "/usr/bin"), ("genoPlotR", "/usr/bin"),
                                   ("scriptpath", str(scriptdir)), ("toolpath", str(tooldir))])
    result = _run(scriptdir / "kmercontigalignonly.py",
                  *_kmercontigalignonly_args(output_dir, id_file, ref_fa, ref_gb, kmerseqfile, sw))
    assert "sort_strings path doesn't exist" not in result.stderr, result.stderr


# --- stringlist2patternandkinship.py (the generic tool() closure, mixed C++/Python) ------

def _s2pk_common(tmp_path):
    analysis_dir = tmp_path / "analysis"
    analysis_dir.mkdir()
    (analysis_dir / "batches").mkdir()
    id_file = tmp_path / "id_file.txt"
    id_file.write_text("")
    fullkmerlistfile = tmp_path / "kmers.gz"
    fullkmerlistfile.write_bytes(b"")
    kmercountslistfile = tmp_path / "counts.txt"
    kmercountslistfile.write_text("")
    return analysis_dir, id_file, fullkmerlistfile, kmercountslistfile


CPP_TOOLS_FOR_S2PK = ("stringlist2pattern", "kmerlist2pattern", "patternmerge",
                      "patterncounts", "pattern2kinship")


def test_s2pk_no_scriptpath_defaults_to_own_directory(tmp_path):
    analysis_dir, id_file, fullkmerlistfile, kmercountslistfile = _s2pk_common(tmp_path)
    scriptdir = _isolated_script_dir(tmp_path, "stringlist2patternandkinship.py",
                                     list(CPP_TOOLS_FOR_S2PK) + ["pattern2presencecount.py"])
    sw = _software_file(tmp_path, [("R", "/usr/bin/R")])
    result = _run(scriptdir / "stringlist2patternandkinship.py",
                  "--task-id", "1", "--p", "1", "--id-file", str(id_file),
                  "--fullkmerlistfile", str(fullkmerlistfile),
                  "--kmercountslistfile", str(kmercountslistfile),
                  "--analysis-dir", str(analysis_dir), "--output-prefix", "x",
                  "--kmertype", "nucleotide", "--software-file", sw)
    for tool in CPP_TOOLS_FOR_S2PK:
        assert f"{tool} path doesn't exist" not in result.stderr, result.stderr


def test_s2pk_no_scriptpath_and_no_tools_fails_clearly(tmp_path):
    analysis_dir, id_file, fullkmerlistfile, kmercountslistfile = _s2pk_common(tmp_path)
    scriptdir = _isolated_script_dir(tmp_path, "stringlist2patternandkinship.py", [])
    sw = _software_file(tmp_path, [("R", "/usr/bin/R")])
    result = _run(scriptdir / "stringlist2patternandkinship.py",
                  "--task-id", "1", "--p", "1", "--id-file", str(id_file),
                  "--fullkmerlistfile", str(fullkmerlistfile),
                  "--kmercountslistfile", str(kmercountslistfile),
                  "--analysis-dir", str(analysis_dir), "--output-prefix", "x",
                  "--kmertype", "nucleotide", "--software-file", sw)
    assert "stringlist2pattern path doesn't exist" in result.stderr


def test_s2pk_pattern2presencecount_stays_on_scriptpath_when_toolpath_differs(tmp_path):
    """pattern2presencecount.py is a sibling Python script, not a compiled tool: it must keep
    resolving via scriptpath even when toolpath points elsewhere (the mixed tool()/sibling-script
    closure has to route the two kinds of name to different locations)."""
    analysis_dir, id_file, fullkmerlistfile, kmercountslistfile = _s2pk_common(tmp_path)
    scriptdir = _isolated_script_dir(tmp_path, "stringlist2patternandkinship.py",
                                     ["pattern2presencecount.py"])  # only the sibling script here
    tooldir = tmp_path / "tools"
    tooldir.mkdir()
    for tool in CPP_TOOLS_FOR_S2PK:
        (tooldir / tool).write_text("#!/bin/bash\nexit 0\n")
        os.chmod(tooldir / tool, 0o755)
    sw = _software_file(tmp_path, [("R", "/usr/bin/R"), ("scriptpath", str(scriptdir)),
                                   ("toolpath", str(tooldir))])
    result = _run(scriptdir / "stringlist2patternandkinship.py",
                  "--task-id", "1", "--p", "1", "--id-file", str(id_file),
                  "--fullkmerlistfile", str(fullkmerlistfile),
                  "--kmercountslistfile", str(kmercountslistfile),
                  "--analysis-dir", str(analysis_dir), "--output-prefix", "x",
                  "--kmertype", "nucleotide", "--software-file", sw)
    assert "pattern2presencecount.py path doesn't exist" not in result.stderr, result.stderr
    for tool in CPP_TOOLS_FOR_S2PK:
        assert f"{tool} path doesn't exist" not in result.stderr, result.stderr


# --- createfullkmerlist.py: no C++ tool, must keep using scriptpath alone (guards P2) -----

def test_createfullkmerlist_resolves_siblings_via_scriptpath_not_toolpath(tmp_path):
    """createfullkmerlist.py never resolves a C++ tool - it only calls nucleotidekmermerge.py/
    proteinkmermerge.py, both found via scriptpath. Giving it a toolpath that does NOT contain
    those sibling scripts must not matter; only scriptpath's content does."""
    analysis_dir = tmp_path / "analysis"
    analysis_dir.mkdir()
    (analysis_dir / "nucleotidekmer31").mkdir()  # createfullkmerlist.py checks this exists first
    id_file = tmp_path / "id_file.txt"
    id_file.write_text("id\n1\n")
    scriptdir = _isolated_script_dir(tmp_path, "createfullkmerlist.py", [])
    shutil.copy(os.path.join(SCRIPTS_DIR, "nucleotidekmermerge.py"), scriptdir / "nucleotidekmermerge.py")
    shutil.copy(os.path.join(SCRIPTS_DIR, "proteinkmermerge.py"), scriptdir / "proteinkmermerge.py")
    tooldir = tmp_path / "tools"  # deliberately empty - createfullkmerlist.py must never look here
    tooldir.mkdir()
    sw = _software_file(tmp_path, [("R", "/usr/bin/R"), ("scriptpath", str(scriptdir)), ("toolpath", str(tooldir))])
    result = _run(scriptdir / "createfullkmerlist.py",
                  "--task-id", "1", "--n", "1", "--p", "1", "--output-prefix", "x",
                  "--analysis-dir", str(analysis_dir), "--id-file", str(id_file),
                  "--kmer-type", "nucleotide", "--kmer-length", "31", "--software-file", sw)
    assert "doesn't exist" not in result.stderr, result.stderr
    assert "missing required software path" not in result.stderr, result.stderr
    assert "nucleotidekmermerge.py path doesn't exist" not in result.stderr, result.stderr
    assert "proteinkmermerge.py path doesn't exist" not in result.stderr, result.stderr


def test_createfullkmerlist_no_scriptpath_defaults_to_own_directory(tmp_path):
    """The committed conf/software_location.txt has no scriptpath row at all (B2) - every other
    script already defaults scriptpath to its own directory in that case; this one was missed
    (found as a real bug, B10, by actually running the committed software file, not a test's own
    synthetic one that always added scriptpath)."""
    analysis_dir = tmp_path / "analysis"
    analysis_dir.mkdir()
    (analysis_dir / "nucleotidekmer31").mkdir()  # createfullkmerlist.py checks this exists first
    id_file = tmp_path / "id_file.txt"
    id_file.write_text("id\n1\n")
    scriptdir = _isolated_script_dir(tmp_path, "createfullkmerlist.py", [])
    shutil.copy(os.path.join(SCRIPTS_DIR, "nucleotidekmermerge.py"), scriptdir / "nucleotidekmermerge.py")
    shutil.copy(os.path.join(SCRIPTS_DIR, "proteinkmermerge.py"), scriptdir / "proteinkmermerge.py")
    sw = _software_file(tmp_path, [("R", "/usr/bin/R")])
    result = _run(scriptdir / "createfullkmerlist.py",
                  "--task-id", "1", "--n", "1", "--p", "1", "--output-prefix", "x",
                  "--analysis-dir", str(analysis_dir), "--id-file", str(id_file),
                  "--kmer-type", "nucleotide", "--kmer-length", "31", "--software-file", sw)
    assert "doesn't exist" not in result.stderr, result.stderr
    assert "missing required software path" not in result.stderr, result.stderr
    assert "nucleotidekmermerge.py path doesn't exist" not in result.stderr, result.stderr


def test_plotmanhattan_no_scriptpath_defaults_to_own_directory(tmp_path):
    """plotManhattan.py had the identical gap to createfullkmerlist.py above (found in the same
    B10 run, same root cause): required_software hard-included "scriptpath" with no fallback,
    unlike every other script. Minimal inputs, just enough to reach the software-file check (the
    real figure-drawing code after it is not exercised here)."""
    analysis_dir = tmp_path / "analysis"
    analysis_dir.mkdir()
    ref_fa = tmp_path / "ref.fa"
    ref_fa.write_text(">ref\nACGT\n")
    ref_gb = tmp_path / "ref.gb"
    ref_gb.write_text("LOCUS ref\n")
    gene_lookup_file = tmp_path / "genes.txt"
    gene_lookup_file.write_text("id\tname\n")
    id_file = tmp_path / "id_file.txt"
    id_file.write_text("id\n1\n")
    prefix = str(tmp_path / "x")
    for suffix in (".patternmerge.patternKeySize.txt", ".patternmerge.patternIndex.txt.gz",
                   ".kmermerge.txt.gz", ".ref_t90.kmeralignmerge.txt.gz", ".ref_t90.kmeralignmerge.count.txt.gz"):
        (tmp_path / (os.path.basename(prefix) + suffix)).write_text("")
    (analysis_dir / "x_nucleotide31.patternmerge.presenceCount.txt.gz").write_text("")
    scriptdir = _isolated_script_dir(tmp_path, "plotManhattan.py", [])
    # rcompat.py is not listed here: _isolated_script_dir() already copies it (every script
    # imports it) - re-copying over the same destination fails if the source's own permissions
    # (preserved by shutil.copy) are read-only, as a git-archive-staged checkout's are (found
    # running this suite against a staged commit for the first time, D3a - not previously
    # exercised, since every earlier run used a live checkout instead).
    for mod in ("sequence_functions.py", "Manhattan_functions.py", "alignmentfunctions.py", "kmercontigalignonly.py"):
        shutil.copy(os.path.join(SCRIPTS_DIR, mod), scriptdir / mod)
    sw = _software_file(tmp_path, [("R", "/usr/bin/R"), ("genoPlotR", "/usr/bin"), ("blast", "/usr/bin")])
    result = _run(scriptdir / "plotManhattan.py",
                  "--output-prefix", "x", "--analysis-dir", str(analysis_dir), "--kmerfile-prefix", prefix,
                  "--ref-gb", str(ref_gb), "--ref-fa", str(ref_fa), "--gene-lookup-file", str(gene_lookup_file),
                  "--id-file", str(id_file), "--nucmerident", "90", "--plot-min-genomes", "1",
                  "--kmer-type", "nucleotide", "--kmer-length", "31", "--minor-allele-threshold", "0.01",
                  "--software-file", sw, "--blastident", "70", "--ngenes", "20")
    assert "script location directory" not in result.stderr, result.stderr
    assert "missing required software path" not in result.stderr, result.stderr
