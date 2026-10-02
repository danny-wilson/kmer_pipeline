"""R-compatibility helpers for the Python port of kmer_pipeline.

Each helper reproduces what R 4.1.3 does in the 2022-10-26 container image, so
that the Python scripts write the same files as the R scripts they replace:

- ordering:   r_order, r_unique, r_sort_strings, r_collate_key
- numbers:    r_format_num, r_as_character, r_cat, r_paste, r_paste0
- files:      r_dir_create, r_open, r_scan_lines, r_cat_lines, r_read_table, r_write_table
- language:   r_stop, r_colon, r_index, r_seq, r_as_integer, r_as_numeric
- commands:   r_system, r_system2, r_system_intern
- statistics: neg_log10_pchisq1_r (R's pchisq, bit for bit), neg_log10_pchisq1

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
# Language
# --------------------------------------------------------------------------


class RError(RuntimeError):
    """An error raised where the R original calls stop()."""


def r_stop(*parts):
    """R's stop(...): the parts pasted together with no separator (as.character
    on each). The traceback gives the line number R could not."""
    raise RError("".join(r_as_character(p) for p in parts))


def r_colon(a, b):
    """R's a:b for integers: inclusive, and descending when a > b (so 1:0 is
    c(1, 0)). Returned as a Python list of ints."""
    return list(range(a, b + 1)) if a <= b else list(range(a, b - 1, -1))


def r_as_integer(s):
    """as.integer() of a command-line string: None (NA) if it isn't a number,
    otherwise truncated towards zero, as R does."""
    v = r_as_numeric(s)
    if v is None or math.isnan(v) or math.isinf(v) or abs(v) > 2147483647:
        return None
    return int(v)


def r_as_numeric(s):
    """as.numeric() of a string: None (NA) if R can't read it as a number."""
    if s is None:
        return None
    t = s.strip()
    if not _DBL_RE.match(t) or t in ("NA", "-NA", "+NA"):
        return None
    t = t.lstrip("+")
    if t.lower().lstrip("-").startswith("0x"):
        return float.fromhex(t)
    return float(t)  # Python reads nan, inf and infinity in any case, as R_strtod does


def r_as_numeric_value(v):
    """as.numeric() of one value from a data frame column (r_read_table types):
    logical TRUE/FALSE -> 1/0, numbers as floats, strings as as.numeric(); NA -> None."""
    if _is_na(v):
        return None
    if isinstance(v, (bool, np.bool_)):
        return 1.0 if v else 0.0
    if isinstance(v, (int, float, np.integer, np.floating)):
        return float(v)
    return r_as_numeric(str(v))


def r_dir_create(path):
    """R's dir.create(path): never fails. If the directory already exists (e.g.
    created at the same moment by a parallel task) or can't be made, R warns and
    returns FALSE."""
    import os
    try:
        os.mkdir(path)
        return True
    except FileExistsError:
        print(f"Warning message:\nIn dir.create({path!r}) : '{path}' already exists", file=sys.stderr, flush=True)
    except OSError as e:
        print(f"Warning message:\nIn dir.create({path!r}) : cannot create dir '{path}', reason '{e.strerror}'",
              file=sys.stderr, flush=True)
    return False


def r_pipe(cmd):
    """The text R reads from pipe(cmd) (e.g. scan(pipe(cmd))): /bin/sh, stdout."""
    _flush()
    return subprocess.run(cmd, shell=True, executable="/bin/sh", stdout=subprocess.PIPE).stdout.decode("utf-8")


def r_seq(from_, to, by):
    """seq(from, to, by) for doubles, as R's seq.default (with its 1e-10 fuzz)."""
    from_, to, by = float(from_), float(to), float(by)
    delta = to - from_
    if delta == 0 and to == 0:
        return [to]
    n = delta / by
    if n < 0:
        raise ValueError("wrong sign in 'by' argument")
    if abs(delta) / max(abs(to), abs(from_)) < 100 * 2.220446049250313e-16:
        return [from_]
    n = int(n + 1e-10)
    x = [from_ + k * by for k in range(n + 1)]
    return [min(v, to) for v in x] if by > 0 else [max(v, to) for v in x]


def r_index(x, idx):
    """x[idx] with R's positive 1-based indices: 0 is dropped and an index past
    the end gives NA (None)."""
    return [None if i > len(x) else x[i - 1] for i in idx if i != 0]


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


def r_dir(path, glob=None, full_names=False):
    """R's dir(path, pattern = glob2rx(glob), full.names =): file names in R's
    sort order. glob2rx makes * any characters and ? one, anchored at both ends."""
    import fnmatch
    import os
    names = [n for n in os.listdir(path) if not n.startswith(".")]  # all.files = FALSE
    if glob is not None:
        names = [n for n in names if fnmatch.fnmatchcase(n, glob.replace("[", "[[]"))]
    if full_names:  # as R: path "/" name, keeping any "//"
        names = [path + "/" + n for n in names]
    return r_sort_strings(names)


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


def r_signif(x, digits):
    """R's signif() (fprec in src/nmath/fprec.c): round half to even at
    `digits` significant digits."""
    x = float(x)
    if math.isnan(x) or math.isinf(x) or x == 0:
        return x
    dig = max(1, int(round(digits)))
    sgn = -1.0 if x < 0 else 1.0
    x = abs(x)
    l10 = math.log10(x)
    e10 = int(dig - 1 - math.floor(l10))
    if e10 > 0:
        p10 = 10.0 ** e10
        return sgn * (float(np.rint(x * p10)) / p10)
    p10 = 10.0 ** (-e10)
    return sgn * (float(np.rint(x / p10)) * p10)


def r_formatC_fg(x, digits, keep_trailing_zeros=True):
    """formatC(x, digits =, format = "fg", flag = "#") for one double (R's
    str_signif): fixed notation with `digits` significant digits."""
    x = float(x)
    if math.isnan(x):
        return "NA"
    if math.isinf(x):
        return "Inf" if x > 0 else "-Inf"
    if x == 0:
        return "0"
    dig = digits
    xxx = abs(x)
    iex = int(math.floor(math.log10(xxx) + 1e-12))
    X = round(xxx / 10.0 ** iex + 1e-12, dig - 1)
    xx = x
    if iex > 0 and X >= 10:
        xx = X * 10.0 ** iex
        iex += 1
    if iex == -4 and abs(xx) < 1e-4:
        iex = -5
    if iex < -4:
        return "%#.*f" % (dig - 1 - iex, xx)
    return "%#.*g" % (iex + 1 if iex >= dig else dig, xx)


def r_s3(x, digits=3):
    """The reports' s3(): gsub("\\.$", "", formatC(signif(x, 3), digits = 3, format = "fg", flag = "#"))."""
    if x is None:  # NA; formatC right-justifies non-finite values to width digits + 1
        return "NA".rjust(digits + 1)
    x = float(x)
    if math.isnan(x):
        return "NaN".rjust(digits + 1)
    if math.isinf(x):
        return ("Inf" if x > 0 else "-Inf").rjust(digits + 1)
    out = r_formatC_fg(r_signif(x, digits), digits)
    return out[:-1] if out.endswith(".") else out


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


def r_paste_collapse(x, collapse=""):
    """paste(x, collapse =) for a character vector; NA (None) becomes "NA"."""
    return collapse.join("NA" if v is None else v for v in x)


# --------------------------------------------------------------------------
# read.table / write.table
# --------------------------------------------------------------------------

_INT_RE = re.compile(r"^\s*[+-]?[0-9]+$")
_DBL_RE = re.compile(
    r"^\s*[+-]?(?:(?:[0-9]+\.?[0-9]*|\.[0-9]+)(?:[eE][+-]?[0-9]+)?"
    r"|0[xX][0-9a-fA-F]+(?:\.[0-9a-fA-F]*)?(?:[pP][+-]?[0-9]+)?"
    r"|(?i:inf(?:inity)?|nan)|NA)\s*$")
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
            return float(t)
        return np.array([np.nan if blank(v) else dbl(v) for v in values], dtype=float)
    return pd.array(values, dtype=object)


def r_open(file):
    """Open a file for reading text as R's file() does: gzip, bzip2 and xz
    compression are detected from the first bytes, and LF, CRLF and CR all end
    a line. Bytes that are not UTF-8 pass through unchanged (surrogateescape)."""
    import bz2
    import gzip
    import lzma
    with open(file, "rb") as f:
        magic = f.read(6)
    if magic[:2] == b"\x1f\x8b":
        opener = gzip.open
    elif magic[:3] == b"BZh":
        opener = bz2.open
    elif magic[:6] == b"\xfd7zXZ\x00":
        opener = lzma.open
    else:
        opener = open
    return opener(file, "rt", encoding="utf-8", errors="surrogateescape", newline=None)


def r_scan_lines(file, quiet=True):
    """scan(file, what = character(0), sep = "\\n", quiet =): the lines of the
    file (compressed or not; any line ending), skipping empty lines, with no
    quote or comment processing. Unless quiet, writes R's "Read N items" to stderr."""
    with r_open(file) as f:
        lines = [l for l in f.read().split("\n") if l != ""]
    if not quiet:
        print(f"Read {len(lines)} item{'' if len(lines) == 1 else 's'}", file=sys.stderr, flush=True)
    return lines


def r_cat_lines(x, file, append=False):
    """cat(x, file = file, sep = "\\n", append =) for a character vector: each
    element followed by a newline (R ends the output with the separator when it
    contains a newline). An empty vector gives a single newline, as in R."""
    with open(file, "a" if append else "w", encoding="utf-8", errors="surrogateescape") as f:
        f.write("".join(v + "\n" for v in x) if len(x) else "\n")


def r_read_table(file, header=None, sep="", quote="\"'", comment_char="#",
                 na_strings=("NA",), check_names=True):
    """R's read.table() with the arguments the pipeline uses (as.is is always
    TRUE, as stringsAsFactors is FALSE in R 4.1). header=None means R's missing
    header: TRUE if the first line has one field fewer than the data. As in R, a
    header one field short makes the first column the row names (the index), and
    the number of columns comes from the first five lines. Column types follow
    R: logical -> pandas "boolean", integer -> "Int64", double -> float64,
    character -> object (NA as None)."""
    with r_open(file) as f:
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


# --- R 4.1.3's pchisq(D, df = 1, lower.tail = FALSE, log.p = TRUE), i.e.
# pgamma(D/2, shape 0.5, lower.tail = FALSE, log.p = TRUE), ported from src/nmath
# (pgamma.c, dpois.c, bd0.c, stirlerr.c) for shape 0.5 only, so that -log10 p is
# bit-identical to R's. Constants R computes with lgammafn() are R's own values.
_DBL_EPSILON = 2.220446049250313e-16
_DBL_MIN = 2.2250738585072014e-308
_M_LN2 = math.log(2)
_M_2PI = 6.283185307179586476925286766559
_M_CUTOFF = _M_LN2 * 1024 / _DBL_EPSILON          # M_LN2 * DBL_MAX_EXP / DBL_EPSILON
_SCALEFACTOR = 4294967296.0 ** 8                   # SQR(SQR(SQR(2^32)))
_LGAMMA_1_5 = float.fromhex("-0x1.eeb95b094c18ep-4")   # R's lgammafn(1.5) = lgamma1p(0.5)
_LGAMMA_0_5 = float.fromhex("0x1.250d048e7a1bdp-1")    # R's lgammafn(0.5)
_STIRLERR_0_5 = 0.1534264097200273452913848         # sferr_halves[1]


def _log1_exp(x):
    """R_Log1_Exp: log(1 - exp(x)) for x <= 0."""
    return math.log(-math.expm1(x)) if x > -_M_LN2 else math.log1p(-math.exp(x))


def _bd0(x, np_):
    if not math.isfinite(x) or not math.isfinite(np_) or np_ == 0.0:
        return math.nan
    if abs(x - np_) < 0.1 * (x + np_):
        v = (x - np_) / (x + np_)
        s = (x - np_) * v
        if abs(s) < _DBL_MIN:
            return s
        ej = 2 * x * v
        v *= v
        for j in range(1, 1000):
            ej *= v
            s_ = s
            s += ej / ((j << 1) + 1)
            if s == s_:
                return s
    return x * math.log(x / np_) + np_ - x


def _dpois_raw_log(x, lam):
    """dpois_raw(x, lambda, give_log = TRUE) for x = 0.5."""
    if lam == 0:
        return 0.0 if x == 0 else -math.inf
    if not math.isfinite(lam):
        return -math.inf
    if x < 0:
        return -math.inf
    if x <= lam * _DBL_MIN:
        return -lam
    if lam < x * _DBL_MIN:
        return -lam + x * math.log(lam) - _LGAMMA_1_5  # lgammafn(x + 1) for x = 0.5
    return -0.5 * math.log(_M_2PI * x) + (-_STIRLERR_0_5 - _bd0(x, lam))


def _dpois_wrap_log(x_plus_1, lam):
    """dpois_wrap(0.5, lambda, give_log = TRUE)."""
    if not math.isfinite(lam):
        return -math.inf
    if lam > abs(x_plus_1 - 1) * _M_CUTOFF:
        return -lam - _LGAMMA_0_5
    return _dpois_raw_log(x_plus_1, lam) + math.log(x_plus_1 / lam)


def _pd_lower_cf(y, d):
    if y == 0:
        return 0.0
    f0 = y / d
    if abs(y - 1) < abs(d) * _DBL_EPSILON:
        return f0
    if f0 > 1.0:
        f0 = 1.0
    c2, c4 = y, d
    a1, b1, a2, b2 = 0.0, 1.0, y, d
    while b2 > _SCALEFACTOR:
        a1 /= _SCALEFACTOR
        b1 /= _SCALEFACTOR
        a2 /= _SCALEFACTOR
        b2 /= _SCALEFACTOR
    i, of, f = 0.0, -1.0, 0.0
    while i < 200000:
        i += 1
        c2 -= 1
        c3 = i * c2
        c4 += 2
        a1 = c4 * a2 + c3 * a1
        b1 = c4 * b2 + c3 * b1
        i += 1
        c2 -= 1
        c3 = i * c2
        c4 += 2
        a2 = c4 * a1 + c3 * a2
        b2 = c4 * b1 + c3 * b2
        if b2 > _SCALEFACTOR:
            a1 /= _SCALEFACTOR
            b1 /= _SCALEFACTOR
            a2 /= _SCALEFACTOR
            b2 /= _SCALEFACTOR
        if b2 != 0:
            f = a2 / b2
            af = abs(f)
            if abs(f - of) <= _DBL_EPSILON * (af if f0 < af else f0):
                return f
            of = f
    return f


def r_pchisq1_upper_log(D):
    """R's pchisq(D, 1, lower.tail = FALSE, log.p = TRUE), bit for bit."""
    x = float(D) / 2.0
    alph = 0.5
    if math.isnan(x):
        return x
    if x <= 0:
        return 0.0
    if math.isinf(x):
        return -math.inf
    if x < 1:  # pgamma_smallx
        s, c, n = 0.0, alph, 0.0
        while True:
            n += 1
            c *= -x / n
            term = c / (alph + n)
            s += term
            if not abs(term) > _DBL_EPSILON * abs(s):
                break
        lf2 = alph * math.log(x) - _LGAMMA_1_5
        return _log1_exp(math.log1p(s) + lf2)
    # alph - 1 < x and alph < 0.8 * (x + 50): always for x >= 1 and shape 0.5
    d = _dpois_wrap_log(alph, x)
    if x * _DBL_EPSILON > 1 - alph:
        s = 0.0
    else:
        f = _pd_lower_cf(alph, x - (alph - 1)) * x / alph
        s = math.log(f)
    return s + d


_NEG_LOG_10 = float.fromhex("-0x1.26bb1bbb55516p+1")  # R's -log(10)


def neg_log10_pchisq1_r(D):
    """-log10 p as R computes it: pchisq(D, 1, lower = FALSE, log.p = TRUE) / -log(10),
    bit-identical (D <= 0 gives 0, NaN gives NaN)."""
    D = np.asarray(D, dtype=float)
    out = np.array([r_pchisq1_upper_log(v) / _NEG_LOG_10 for v in D.ravel()], dtype=float).reshape(D.shape)
    return out + 0.0  # -0.0 -> 0.0


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
