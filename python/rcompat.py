"""R-compatibility helpers for the Python port of kmer_pipeline.

Each helper reproduces what R 4.1.3 does in the 2022-10-26 container image, so
that the Python scripts write the same files as the R scripts they replace:

- ordering:   r_order, r_unique, r_sort_strings, r_collate_key
- numbers:    r_format_num, r_as_character, r_cat, r_paste, r_paste0
- tables:     r_read_table, r_write_table
- commands:   r_system, r_system2, r_system_intern
- statistics: neg_log10_pchisq1

R's integers and doubles print differently (100000L is "100000", 1e5 is
"1e+05"), so callers pass int for R integers and float for R doubles.
Paths are built by string concatenation, as R's paste0 does; os.path.join would
remove the "//" that R keeps.
"""
import math
import re
import subprocess
import sys

import numpy as np
import pandas as pd
from scipy.special import log_ndtr

NA = None

# --------------------------------------------------------------------------
# Script start-up
# --------------------------------------------------------------------------


def script_setup(script_path, announce=True):
    """Common start-up for every pipeline script: line-buffered output, so
    progress lines, command output and tracebacks reach the log in order; refuse
    to run with user site-packages enabled, which can shadow the image's
    packages; and (if announce) print the script's path and the git SHA in the
    STAGED_SHA file beside it, if any, so a log shows which code ran."""
    import os
    import site
    sys.stdout.reconfigure(line_buffering=True)
    sys.stderr.reconfigure(line_buffering=True)
    if site.ENABLE_USER_SITE:
        raise RuntimeError("Python user site-packages are enabled; set PYTHONNOUSERSITE=1")
    if announce:
        path = os.path.abspath(script_path)
        stamp = os.path.join(os.path.dirname(path), "STAGED_SHA")
        sha = open(stamp).read().strip() if os.path.exists(stamp) else "not staged"
        print(f"Running {path} ({sha})")


# --------------------------------------------------------------------------
# Ordering
# --------------------------------------------------------------------------


def _is_na(v):
    return v is None or v is pd.NA or (isinstance(v, (float, np.floating)) and math.isnan(v))


def r_order(x, decreasing=False):
    """R's order(x, decreasing=) for one numeric vector: stable (ties keep their
    original order whichever the direction) with NA and NaN last. Returns
    0-based indices as a numpy int64 array."""
    x = np.asarray(x, dtype=float)
    na = np.isnan(x)
    idx = np.flatnonzero(~na)
    vals = x[idx]
    o = np.argsort(-vals if decreasing else vals, kind="stable")
    return np.concatenate([idx[o], np.flatnonzero(na)]).astype(np.int64)


def r_unique(x):
    """R's unique(): first occurrences, in order of appearance; all NAs count as one."""
    out, seen, seen_na = [], set(), False
    for v in x:
        if _is_na(v):
            if not seen_na:
                seen_na = True
                out.append(v)
        elif v not in seen:
            seen.add(v)
            out.append(v)
    return out


# R sorts strings with ICU root collation (capabilities("ICU") is TRUE in the
# image, LC_COLLATE en_US.UTF-8). Punctuation is not ignored: space and
# punctuation sort first, in this order, then digits, then letters; letters are
# compared without case first and lower case wins a tie.
_ICU_ORDER = " _-,;:!?.'\"()[]{}@*/\\&#%`^+<=>|~$0123456789abcdefghijklmnopqrstuvwxyz"
_PRIMARY = {c: i for i, c in enumerate(_ICU_ORDER)}
_PRIMARY.update({c.upper(): _PRIMARY[c] for c in "abcdefghijklmnopqrstuvwxyz"})


def r_collate_key(s):
    """Sort key giving R's (ICU) order for printable ASCII strings. Raises
    ValueError on any other character rather than risk a silent difference."""
    try:
        primary = tuple(_PRIMARY[c] for c in s)
    except KeyError as e:
        raise ValueError(f"r_collate_key: no R collation weight for character {e.args[0]!r} in {s!r}")
    return primary, tuple(c.isupper() for c in s)


def r_sort_strings(x, decreasing=False):
    """R's sort() for a character vector: NAs removed, ICU order, stable."""
    return sorted((s for s in x if not _is_na(s)), key=r_collate_key, reverse=decreasing)


# --------------------------------------------------------------------------
# Number formatting (R's formatReal for a single value)
# --------------------------------------------------------------------------

_KP_MAX = 22
_TBL = [np.longdouble(10) ** i for i in range(_KP_MAX + 1)]


def _scientific(x, digits):
    """Port of scientific() in R's src/main/format.c (long double branch).
    Returns (neg, kpower, nsig, roundingwidens)."""
    if x == 0.0:
        return 0, 0, 1, False
    neg = 1 if x < 0 else 0
    r = -x if neg else x
    if digits >= 16:
        # R uses sprintf("%#.*e") once digits exceed DBL_DIG
        buff = "%#.*e" % (digits - 1, r)
        kpower = int(buff[digits + 2:])
        j = digits
        while j > 0 and buff[j] == "0":
            j -= 1
        return neg, kpower, j, False
    kp = math.floor(math.log10(r)) - digits + 1
    r_prec = np.longdouble(r)
    if abs(kp) < 10:
        if kp > 0:
            r_prec /= _TBL[kp]
        elif kp < 0:
            r_prec *= _TBL[-kp]
    elif kp <= -308:  # R_dec_min_exponent
        r_prec = (np.longdouble(r) * np.longdouble(1e303)) / np.longdouble(10) ** (kp + 303)
    else:
        r_prec /= np.longdouble(10) ** kp
    if r_prec < _TBL[digits - 1]:
        r_prec *= 10
        kp -= 1
    alpha = float(np.rint(r_prec))
    nsig = digits
    for _ in range(digits):
        alpha /= 10.0
        if alpha == math.floor(alpha):
            nsig -= 1
        else:
            break
    if nsig == 0 and digits > 0:
        nsig = 1
        kp += 1
    kpower = kp + digits - 1
    rgt = min(max(digits - kpower, 0), _KP_MAX)
    fuzz = 0.5 / float(_TBL[rgt])
    roundingwidens = 0 < kpower <= _KP_MAX and r < float(_TBL[kpower]) - fuzz
    return neg, kpower, nsig, roundingwidens


def r_format_num(x, digits=7):
    """Format one R double as R does with `digits` significant digits
    (7: cat/print; 15: as.character, paste, write.table). Never gives "-0"."""
    x = float(x)
    if math.isnan(x):
        return "NaN"
    if math.isinf(x):
        return "Inf" if x > 0 else "-Inf"
    if x == 0:
        return "0"
    neg, kpower, nsig, roundingwidens = _scientific(x, digits)
    left = kpower + 1
    if roundingwidens:
        left -= 1
    sleft = neg + (1 if left <= 0 else left)
    rgt = max(nsig - left, 0)
    rgt = min(rgt, 350)
    e_digits = 2 if (kpower >= 100 or kpower <= -100) else 1
    d = nsig - 1
    width_e = neg + (d > 0) + d + 4 + e_digits
    width_f = sleft + rgt + (rgt != 0)
    if width_f <= width_e:  # options(scipen = 0)
        s = "%.*f" % (rgt, x)
    else:
        s = "%.*e" % (d, x)
    return s


def r_str(v, digits):
    """One R value as text: None/NaN -> "NA", bool -> TRUE/FALSE, int as is,
    float through r_format_num, anything else str()."""
    if v is None or v is pd.NA:
        return "NA"
    if isinstance(v, (bool, np.bool_)):
        return "TRUE" if v else "FALSE"
    if isinstance(v, (int, np.integer)):
        return str(int(v))
    if isinstance(v, (float, np.floating)):
        if math.isnan(v):
            return "NA"  # R's NA_real_; NaN proper is rare in this pipeline
        return r_format_num(v, digits)
    return str(v)


def r_as_character(v):
    """R's as.character() for one value (15 significant digits for doubles)."""
    return r_str(v, 15)


def _flatten(args):
    for a in args:
        if isinstance(a, (list, tuple, np.ndarray, pd.Series)):
            yield from a
        else:
            yield a


def r_cat(*args, sep=" ", file=None):
    """R's cat(): each element formatted separately (7 significant digits),
    separated by `sep`; no newline is added. Vectors (list, tuple, numpy,
    pandas) are expanded element by element. Writes to stdout and flushes."""
    text = sep.join(r_str(v, 7) for v in _flatten(args))
    out = sys.stdout if file is None else file
    out.write(text)
    out.flush()
    return text


def r_paste(*args, sep=" "):
    """R's paste() with scalar arguments (as.character on each, 15 digits)."""
    return sep.join(r_as_character(v) for v in args)


def r_paste0(*args):
    return r_paste(*args, sep="")


# --------------------------------------------------------------------------
# read.table / write.table
# --------------------------------------------------------------------------

_INT_RE = re.compile(r"^\s*[+-]?[0-9]+$")
_DBL_RE = re.compile(
    r"^\s*[+-]?(?:(?:[0-9]+\.?[0-9]*|\.[0-9]+)(?:[eE][+-]?[0-9]+)?"
    r"|0[xX][0-9a-fA-F]+(?:\.[0-9a-fA-F]*)?(?:[pP][+-]?[0-9]+)?"
    r"|Inf|inf|NaN|NA)\s*$")
_LOGICAL = {"T": True, "TRUE": True, "true": True, "True": True,
            "F": False, "FALSE": False, "false": False, "False": False}


def _records(text, sep, quote, comment_char):
    """Split text into records of fields as R's scan() does, skipping blank lines.
    Whitespace separated (sep None): a quote counts only at the start of a field
    and the closing quote ends the field. Explicit sep: .csv-style, a quote
    anywhere in a field opens a quoted section ("" inside it is a quote), which
    may span separators and newlines. A comment character outside quotes ends
    the line."""
    records, i, n = [], 0, len(text)

    def unterminated():
        raise ValueError("r_read_table: EOF within quoted string")

    while i < n:
        fields = []
        if sep is None:
            while True:
                while i < n and text[i] in " \t":
                    i += 1
                if i >= n or text[i] == "\n":
                    i += 1
                    break
                c = text[i]
                if comment_char and c == comment_char:
                    j = text.find("\n", i)
                    i = n if j < 0 else j
                elif c in quote:
                    j = text.find(c, i + 1)
                    if j < 0:
                        unterminated()
                    fields.append(text[i + 1:j])
                    i = j + 1
                else:
                    j = i
                    while j < n and text[j] not in " \t\n" and text[j] != comment_char:
                        j += 1
                    fields.append(text[i:j])
                    i = j
            if fields:
                records.append(fields)
        else:
            cur = []
            while True:
                if i >= n or text[i] == "\n":
                    i += 1
                    fields.append("".join(cur))
                    break
                c = text[i]
                if c == sep:
                    fields.append("".join(cur))
                    cur = []
                    i += 1
                elif comment_char and c == comment_char:
                    j = text.find("\n", i)
                    i = n if j < 0 else j
                elif c in quote:
                    i += 1
                    while True:
                        j = text.find(c, i)
                        if j < 0:
                            unterminated()
                        cur.append(text[i:j])
                        i = j + 1
                        if i < n and text[i] == c:
                            cur.append(c)
                            i += 1
                        else:
                            break
                else:
                    cur.append(c)
                    i += 1
            if not (len(fields) == 1 and fields[0].strip(" \t") == ""):
                records.append(fields)
    return records


_R_RESERVED = {"if", "else", "repeat", "while", "function", "for", "next", "break", "in", "TRUE", "FALSE",
               "NULL", "Inf", "NaN", "NA", "NA_integer_", "NA_real_", "NA_character_", "NA_complex_"}


def r_make_names(names):
    """R's make.names(names, unique = TRUE). As in R, names that make.names
    leaves unchanged keep priority when duplicates are suffixed."""
    out = []
    for s in names:
        if not re.match(r"[A-Za-z]|\.(?![0-9])", s):
            s = "X" + s
        s = re.sub(r"[^A-Za-z0-9._]", ".", s)
        if s in _R_RESERVED:
            s += "."
        out.append(s)
    o = [i for i in range(len(out)) if out[i] == names[i]] + [i for i in range(len(out)) if out[i] != names[i]]
    res = list(out)
    for i, u in zip(o, r_make_unique([out[i] for i in o])):
        res[i] = u
    return res


def r_make_unique(names):
    """R's make.unique(names): the second "a" becomes "a.1", the third "a.2",
    skipping any name already present."""
    taken = set(names)
    seen, cnt, res = set(), {}, []
    for s in names:
        if s in seen:
            k = cnt.get(s, 1)
            while f"{s}.{k}" in taken:
                k += 1
            new = f"{s}.{k}"
            cnt[s] = k + 1
            taken.add(new)
            res.append(new)
        else:
            seen.add(s)
            res.append(s)
    return res


def _type_convert(values):
    """R's type.convert(as.is = TRUE) on one column of strings (None = NA)."""
    def blank(v):
        return v is None or v.strip() == ""
    present = [v for v in values if not blank(v)]
    if all(v in _LOGICAL for v in present):
        return pd.array([None if blank(v) else _LOGICAL[v] for v in values], dtype="boolean")
    if all(_INT_RE.match(v) and -2147483647 <= int(v) <= 2147483647 for v in present):
        return pd.array([None if blank(v) else int(v) for v in values], dtype="Int64")
    if all(_DBL_RE.match(v) for v in present):
        def dbl(v):
            t = v.strip().lstrip("+")
            if t in ("NA", "-NA"):
                return np.nan
            if t.lower().lstrip("-").startswith("0x"):
                return float.fromhex(t)
            return float(t.replace("Inf", "inf").replace("NaN", "nan"))
        return np.array([np.nan if blank(v) else dbl(v) for v in values], dtype=float)
    return pd.array(values, dtype=object)


def r_read_table(file, header=None, sep="", quote="\"'", comment_char="#",
                 na_strings=("NA",), check_names=True):
    """R's read.table() with the arguments the pipeline uses (as.is is always
    TRUE, as stringsAsFactors is FALSE in R 4.1). header=None means R's missing
    header: TRUE if the first line has one field fewer than the data. As in R, a
    header one field short makes the first column the row names (the index), and
    the number of columns comes from the first five lines. Column types follow
    R: logical -> pandas "boolean", integer -> "Int64", double -> float64,
    character -> object (NA as None)."""
    with open(file, encoding="utf-8") as f:
        text = f.read()
    rows = _records(text, None if sep == "" else sep, quote, comment_char)
    if not rows:
        raise ValueError("no lines available in input")
    first = [s.strip(" \t") for s in rows[0]]  # header is read with strip.white
    cols = max(len(r) for r in rows[:5])
    row_names = cols - len(first) == 1
    if header is None:
        header = row_names
    if not header:
        row_names = False
    if header:
        names, rows = first, rows[1:]
        if len(names) + row_names < cols:
            raise ValueError("more columns than column names")
        cols = len(names) + row_names
        names = r_make_names(names) if check_names else names
    else:
        names = [f"V{j + 1}" for j in range(cols)]
    for i, r in enumerate(rows):
        if len(r) != cols:
            raise ValueError(f"line {i + 1} did not have {cols} elements")
    columns = [[None if r[j] in na_strings else r[j] for r in rows] for j in range(cols)]
    index = None
    if row_names:
        index, columns = columns[0], columns[1:]
    df = pd.DataFrame({n: _type_convert(c) for n, c in zip(names, columns)})
    if index is not None:
        df.index = index
    return df


def r_write_table(rows, file, col_names=None, sep="\t", na="NA"):
    """R's write.table(x, row.names = FALSE, quote = FALSE). `rows` is a
    DataFrame (columns kept with their R types) or a list of row tuples; values
    are written as R does (doubles to 15 significant digits). Column names are
    written when `col_names` is given (a list), or True to use a DataFrame's."""
    if isinstance(rows, pd.DataFrame):
        if col_names is True:
            col_names = list(rows.columns)
        data = rows.astype(object).itertuples(index=False, name=None)
    else:
        data = rows
    with open(file, "w", encoding="utf-8") as f:
        if col_names:
            f.write(sep.join(col_names) + "\n")
        for r in data:
            f.write(sep.join(na if _is_na(v) else r_str(v, 15) for v in r) + "\n")


# --------------------------------------------------------------------------
# External commands. R's three forms, reproduced exactly.
# --------------------------------------------------------------------------


def _flush():
    sys.stdout.flush()
    sys.stderr.flush()


def r_system(cmd):
    """stopifnot(system(cmd) == 0): run with /bin/sh, fail on non-zero exit."""
    _flush()
    rc = subprocess.run(cmd, shell=True, executable="/bin/sh").returncode
    if rc != 0:
        raise RuntimeError(f"command exited with status {rc}: {cmd}")
    return rc


def r_shquote(s):
    """R's shQuote(type = "sh")."""
    if "'" not in s:
        return "'" + s + "'"
    return '"' + re.sub(r'(["$`\\])', r"\\\1", s) + '"'


def r_system2(command, args):
    """stopifnot(system2(command, args) == 0). As in R, the command line is
    shQuote(command) followed by args, run by /bin/sh."""
    return r_system(r_shquote(command) + " " + args)


def r_system_intern(cmd):
    """system(cmd, intern = TRUE): /bin/sh, stdout as a list of lines. The exit
    status is ignored, as in R; a non-zero status is logged as R warns."""
    _flush()
    p = subprocess.run(cmd, shell=True, executable="/bin/sh", stdout=subprocess.PIPE)
    if p.returncode != 0:
        print(f"Warning: running command '{cmd}' had status {p.returncode}", file=sys.stderr, flush=True)
    out = p.stdout.decode("utf-8")
    lines = out.split("\n")
    if lines and lines[-1] == "":
        lines.pop()
    return lines


# --------------------------------------------------------------------------
# Statistics
# --------------------------------------------------------------------------


def neg_log10_pchisq1(D):
    """-log10 of the upper tail of chi-squared(1) at D, as R's
    pchisq(D, 1, lower = FALSE, log.p = TRUE) / -log(10), but finite for large D:
    P(chi2_1 > D) = 2 * Phi(-sqrt(D)). D <= 0 (and NaN) give R's result: 0 and NaN."""
    D = np.asarray(D, dtype=float)
    with np.errstate(invalid="ignore"):
        pos = D > 0
        out = np.where(pos, -(log_ndtr(-np.sqrt(np.where(pos, D, 1.0))) + math.log(2)) / math.log(10), 0.0)
    out = np.where(np.isnan(D), np.nan, out)
    return out + 0.0  # turns -0.0 into 0.0
