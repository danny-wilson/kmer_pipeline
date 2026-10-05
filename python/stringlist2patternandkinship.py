#!/usr/bin/env python3
"""stringlist2patternandkinship.py: define presence/absence patterns of a k-mer list
among the k-mer count files and build the kinship matrix, in batches merged as a
pyramid between tasks. Port of stringlist2patternandkinship.Rscript.

R types are kept where they reach file names or files: R doubles are Python
floats and R integers Python ints, so r_paste0 and r_str format them as R does."""
import argparse
import math
import os
import sys
import time

import numpy as np

import rcompat
from rcompat import r_cat, r_paste0, r_stop

###################################################################################################
## Functions and software paths
###################################################################################################


def size_of(out):
    """as.numeric() of the output of system("ls -l <path> | cut -d ' ' -f5", intern = T):
    the size, or None (R's numeric(0)) when ls printed nothing."""
    if not out:
        return None
    v = rcompat.r_as_numeric(out[0])
    return math.nan if v is None else v


def size_is_zero(size):
    """if(as.numeric(system("ls -l ...", intern = T)) == 0): R stops when ls prints nothing."""
    if size is None:
        raise rcompat.RError("argument is of length zero")
    return size == 0


def r_any_zero(sizes):
    """any(sizes == 0) on R's sapply result; an element that is numeric(0)
    (missing file) makes it NA, which stops R in if()."""
    if any(s == 0 for s in sizes if s is not None):
        return True
    if any(s is None or s != s for s in sizes):
        raise rcompat.RError("missing value where TRUE/FALSE needed")
    return False


def create_pattern_batch(fullkmerlistfile, p, t, stringlist2patternpath, kmerlist2patternpath, kmercountslistfile,
                         kmerlen, mincount, kmertype, output_prefix, output_dir):
    # Pattern batches directory - create if doesn't exist
    batches_dir = output_dir + "/" + r_paste0(kmertype, "kmer", kmerlen, "_patternbatches/")  # file.path
    if not os.path.isdir(batches_dir):
        r_cat("Creating temp directory for pattern batches:", batches_dir, "\n")
        rcompat.r_dir_create(batches_dir)

    # Read kmers
    # Total number of kmers
    n = rcompat.r_as_integer(rcompat.r_pipe("zcat " + fullkmerlistfile + " | wc -l").split()[0])
    if n < 1:
        r_stop("No kmers found in", fullkmerlistfile)
    # Number of kmers (batch size) per process
    nkmersbatch = get_kmer_batch_numbers(n=n, p=p)
    beg = nkmersbatch[t - 1][0]
    end = nkmersbatch[t - 1][1]

    output_prefix_batch = r_paste0(output_prefix, "_", kmertype, kmerlen, ".", beg, "-", end)

    kmersublistfile = r_paste0(batches_dir, output_prefix_batch, ".", t, ".temp_kmerlist.txt.gz")
    r_cat("Creating temp file:", kmersublistfile, "\n")
    cmd = r_paste0("zcat ", fullkmerlistfile, " | head -n ", end, " | tail -n ", end - beg + 1, " | gzip -c > ",
                   kmersublistfile)
    rcompat.r_system(cmd)
    softwarepath = stringlist2patternpath if kmertype == "protein" else kmerlist2patternpath
    cmd = rcompat.r_paste(softwarepath, kmersublistfile, kmercountslistfile, batches_dir + output_prefix_batch,
                          kmerlen, mincount)
    rcompat.r_system(cmd)

    cmd = rcompat.r_paste("rm", kmersublistfile)
    rcompat.r_system(cmd)

    # Check that files have been created and are not empty
    outfiles = [batches_dir + output_prefix_batch + s for s in (".patternKey.txt.gz", ".patternIndex.txt.gz")]
    outfiles_size = [size_of(rcompat.r_system_intern("ls -l " + x + " | cut -d ' ' -f5")) for x in outfiles]
    if r_any_zero(outfiles_size):
        r_stop("One or more task_id ", t, " ", output_prefix_batch,
               " patternKey, patternKeySize or patternIndex files are empty")

    # Write output file prefix to output files
    r_cat("Output file prefix: " + output_prefix_batch, "\n")
    outfile_completed = batches_dir + output_prefix_batch + ".patternbatch.completed.txt"
    rcompat.r_system2("/bin/bash", "-c 'touch " + outfile_completed + "'")

    return {"nkmersbatch": nkmersbatch, "output_prefix_batch": output_prefix_batch, "batches_dir": batches_dir}


def get_kmer_batch_numbers(n, p):
    """Rows (beg, end) of k-mer numbers for each of p tasks, as integers, so they are
    written as plain integers in file names and commands (D1a; R wrote 100000 as "1e+05",
    which head -n rejects)."""
    # Number of kmers (batch size) per process
    b = math.ceil(n / p)
    return [((t - 1) * b + 1, min(t * b, n)) for t in range(1, int(p) + 1)]


def merge_pattern_batches_parameters(nkmersbatch, p, t, output_dir, output_prefix, kmertype, kmerlen):
    prefix = r_paste0(output_prefix, "_", kmertype, kmerlen)
    output_prefix_batches = [r_paste0(prefix, ".", b0, "-", b1) for b0, b1 in nkmersbatch]
    batches_dir = output_dir + "/" + r_paste0(kmertype, "kmer", kmerlen, "_patternbatches/")
    infiles = [batches_dir + o + ".patternKey.txt.gz" for o in output_prefix_batches]
    infiles_completed = [batches_dir + o + ".patternbatch.completed.txt" for o in output_prefix_batches]
    nattempts = 0
    while not all(os.path.exists(f) for f in infiles_completed) or not all(os.path.exists(f) for f in infiles):
        nattempts = nattempts + 1
        if nattempts > 100:
            r_stop("Could not find files", "".join(infiles_completed))
        time.sleep(60)

    n = len(infiles)
    if not all(1 + nkmersbatch[k][1] == nkmersbatch[k + 1][0] for k in range(n - 1)):
        r_stop("Input files do not have consecutive ranges")

    # Determine remaining parameters
    b = math.ceil(n / p)
    if n == 1:
        imax = 1
    else:
        if b == 1:
            r_stop("Pattern merging - cannot have batchsize = 1. Try p < n/2")
        if p != math.ceil(n / b):
            r_cat("Warning: pattern merging - adjusting number of processes to equal ceiling(n/b)\n")
            p = math.ceil(n / b)  # as in R, the caller's p is not changed
        imax = math.ceil(math.log(n) / math.log(b))

    return {"n": n, "b": b, "imax": imax, "infiles": infiles}


def key_size_name(f):
    return f.replace(".patternKey.txt.gz", ".patternKeySize.txt")


def merge_patterns(t, n, b, p, imax, prefix, files, patternmergepath, output_dir, kmertype, kmerlen, batches_dir):
    # Merge
    i = 0.0
    outfile_patternKey = outfile_patternKeySize = outfile_patternIndex = outfile_prefix = None
    while True:
        i = i + 1
        if (n == 1 and i > 1) or not ((t % b ** int(i - 1)) == 0 or (t == p and i <= imax)):
            break
        r_cat("t =", t, "i =", i, "\n")
        if i == 1:
            # First round: merge source files
            outfile_prefix = r_paste0(prefix, "_", kmertype, kmerlen, ".patternmerge.j.", i, ".", t)
            outfile_patternKey = batches_dir + outfile_prefix + ".patternKey.txt.gz"
            outfile_patternKeySize = batches_dir + outfile_prefix + ".patternKeySize.txt"
            outfile_patternIndex = batches_dir + outfile_prefix + ".patternIndex.txt.gz"
            outfile_completed = batches_dir + outfile_prefix + ".patternbatch.completed.txt"
            beg = b * (t - 1) + 1
            end = min(b * t, n)
            if end < beg:
                r_stop("Problem with input arguments, please check")
            infiles = ["NA" if f is None else f for f in rcompat.r_index(files, rcompat.r_colon(int(beg), end))]
            nattempts = 0
            while not all(os.path.exists(f) for f in infiles):
                nattempts = nattempts + 1
                if nattempts > 1000:
                    r_stop("Could not find files", "".join(infiles))
                time.sleep(1)

            # Check files aren't empty
            infiles_size = [size_of(rcompat.r_system_intern("ls -l " + x + " | cut -d ' ' -f5")) for x in infiles]
            if r_any_zero(infiles_size):
                r_stop("One or more file size is zero ", "".join(infiles), "\n")
            # Create patternKeySize files
            for f in infiles:
                rcompat.r_system("zcat " + f + " | wc -l > " + key_size_name(f))
                r_cat("Created ", key_size_name(f), "\n")

            # Merge the patternKey files and redefine the patternKeySize and patternIndex files
            if len(infiles) == 1:
                rcompat.r_system("cp " + infiles[0] + " " + outfile_patternKey)
                rcompat.r_system("cp " + key_size_name(infiles[0]) + " " + outfile_patternKeySize)
                rcompat.r_system("cp " + infiles[0].replace(".patternKey.txt.gz", ".patternIndex.txt.gz") + " "
                                 + outfile_patternIndex)
                rcompat.r_system2("/bin/bash", "-c 'touch " + outfile_completed + "'")
            elif len(infiles) == 2:
                prefixA = infiles[0].replace(".patternKey.txt.gz", "")
                prefixB = infiles[1].replace(".patternKey.txt.gz", "")
                rcompat.r_system(patternmergepath + " " + prefixA + " " + prefixB + " " + batches_dir + outfile_prefix)
                rcompat.r_system2("/bin/bash", "-c 'touch " + outfile_completed + "'")
            else:
                # Create temporary storage filenames
                tmpfile_prefix = [batches_dir + "tmpfile" + str(k) + "." + outfile_prefix for k in (1, 2)]
                tmpfile_patternKey = [x + ".patternKey.txt.gz" for x in tmpfile_prefix]
                tmpfile_patternKeySize = [x + ".patternKeySize.txt" for x in tmpfile_prefix]
                tmpfile_patternIndex = [x + ".patternIndex.txt.gz" for x in tmpfile_prefix]

                # Do the merges (j is 1-based, as in R; tmpfile_prefix[1+(j%%2)] is [j % 2] here)
                for j in range(2, len(infiles) + 1):
                    if j == 2:
                        prefixA = infiles[0].replace(".patternKey.txt.gz", "")
                        prefixB = infiles[1].replace(".patternKey.txt.gz", "")
                        rcompat.r_system(patternmergepath + " " + prefixA + " " + prefixB + " " + tmpfile_prefix[j % 2])
                    else:
                        prefixA = tmpfile_prefix[(j - 1) % 2]
                        prefixB = infiles[j - 1].replace(".patternKey.txt.gz", "")
                        rcompat.r_system(patternmergepath + " " + prefixA + " " + prefixB + " " + tmpfile_prefix[j % 2])
                # Finally move the tmpfiles into place
                # Behaviour: j equals the last value in the list
                rcompat.r_system("mv " + tmpfile_patternKey[j % 2] + " " + outfile_patternKey)
                rcompat.r_system("mv " + tmpfile_patternKeySize[j % 2] + " " + outfile_patternKeySize)
                rcompat.r_system("mv " + tmpfile_patternIndex[j % 2] + " " + outfile_patternIndex)
                # Delete the other tmpfiles
                for f in (tmpfile_patternKey[(j - 1) % 2], tmpfile_patternKeySize[(j - 1) % 2],
                          tmpfile_patternIndex[(j - 1) % 2]):
                    if os.path.lexists(f):
                        os.remove(f)
                # Create file stating that the round is completed
                rcompat.r_system2("/bin/bash", "-c 'touch " + outfile_completed + "'")
        else:
            # Subsequent rounds: merge merged files
            te = math.ceil(t / (b ** (i - 1))) * int(b ** (i - 1))
            outfile_prefix = r_paste0(prefix, "_", kmertype, kmerlen, ".patternmerge.j.", i, ".", te)
            outfile_patternKey = batches_dir + outfile_prefix + ".patternKey.txt.gz"
            outfile_patternKeySize = batches_dir + outfile_prefix + ".patternKeySize.txt"
            outfile_patternIndex = batches_dir + outfile_prefix + ".patternIndex.txt.gz"
            outfile_completed = batches_dir + outfile_prefix + ".patternbatch.completed.txt"
            beg = te - float(b) ** (i - 1) + float(b) ** (i - 2)
            end = min(te, math.ceil(t / b ** (i - 2)) * int(b ** (i - 2)))
            if end < beg:
                r_stop("Problem with input arguments, please check")
            inc = float(b) ** (i - 2)
            js = rcompat.r_seq(beg, end, inc)
            infiles = [r_paste0(batches_dir, prefix, "_", kmertype, kmerlen, ".patternmerge.j.", i - 1, ".", x,
                                ".patternKey.txt.gz") for x in js]
            infiles_completed = [r_paste0(batches_dir, prefix, "_", kmertype, kmerlen, ".patternmerge.j.", i - 1, ".",
                                          x, ".patternbatch.completed.txt") for x in js]
            nattempts = 0
            while not all(os.path.exists(f) for f in infiles_completed) or not all(os.path.exists(f) for f in infiles):
                nattempts = nattempts + 1
                if nattempts > 3600:
                    r_stop("Could not find files", "".join(infiles_completed))
                time.sleep(1)

            infiles_size = [size_of(rcompat.r_system_intern("ls -l " + x + " | cut -d ' ' -f5")) for x in infiles]
            if r_any_zero(infiles_size):
                r_stop("One or more file size is zero ", " ".join(infiles), "\n")
            # Merge the patternKey files and redefine the patternKeySize and patternIndex files
            if len(infiles) == 1:
                rcompat.r_system("cp " + infiles[0] + " " + outfile_patternKey)
                rcompat.r_system("cp " + key_size_name(infiles[0]) + " " + outfile_patternKeySize)
                rcompat.r_system("cp " + infiles[0].replace(".patternKey.txt.gz", ".patternIndex.txt.gz") + " "
                                 + outfile_patternIndex)
                rcompat.r_system2("/bin/bash", "-c 'touch " + outfile_completed + "'")
            elif len(infiles) == 2:
                prefixA = infiles[0].replace(".patternKey.txt.gz", "")
                prefixB = infiles[1].replace(".patternKey.txt.gz", "")
                rcompat.r_system(patternmergepath + " " + prefixA + " " + prefixB + " " + batches_dir + outfile_prefix)
                rcompat.r_system2("/bin/bash", "-c 'touch " + outfile_completed + "'")
            else:
                # Create temporary storage filenames
                ntmp = 2
                tmpfile_prefix = [batches_dir + "tmpfile" + str(k) + "." + outfile_prefix for k in range(1, ntmp + 1)]
                tmpfile_patternKey = [x + ".patternKey.txt.gz" for x in tmpfile_prefix]
                tmpfile_patternKeySize = [x + ".patternKeySize.txt" for x in tmpfile_prefix]
                tmpfile_patternIndex = [x + ".patternIndex.txt.gz" for x in tmpfile_prefix]
                # Do the merges
                for j in range(2, len(infiles) + 1):
                    if j == 2:
                        prefixA = infiles[0].replace(".patternKey.txt.gz", "")
                        prefixB = infiles[1].replace(".patternKey.txt.gz", "")
                        rcompat.r_system(patternmergepath + " " + prefixA + " " + prefixB + " "
                                         + tmpfile_prefix[(j - 1) % ntmp])
                    else:
                        prefixA = tmpfile_prefix[(j - 2) % ntmp]
                        prefixB = infiles[j - 1].replace(".patternKey.txt.gz", "")
                        rcompat.r_system(patternmergepath + " " + prefixA + " " + prefixB + " "
                                         + tmpfile_prefix[(j - 1) % ntmp])
                        # Delete the previous tmpfile
                        for f in (tmpfile_patternKey[(j - 2) % ntmp], tmpfile_patternKeySize[(j - 2) % ntmp],
                                  tmpfile_patternIndex[(j - 2) % ntmp]):
                            if os.path.lexists(f):
                                os.remove(f)
                # Finally move the tmpfiles into place
                # Behaviour: j equals the last value in the list
                rcompat.r_system("mv " + tmpfile_patternKey[(j - 1) % ntmp] + " " + outfile_patternKey)
                rcompat.r_system("mv " + tmpfile_patternKeySize[(j - 1) % ntmp] + " " + outfile_patternKeySize)
                rcompat.r_system("mv " + tmpfile_patternIndex[(j - 1) % ntmp] + " " + outfile_patternIndex)
                rcompat.r_system2("/bin/bash", "-c 'touch " + outfile_completed + "'")

    if t == p:
        create_final_merged_pattern_files(outfile_patternKey=outfile_patternKey,
                                          outfile_patternKeySize=outfile_patternKeySize,
                                          outfile_patternIndex=outfile_patternIndex, outfile_prefix=outfile_prefix,
                                          prefix=prefix, output_dir=output_dir, kmertype=kmertype, kmerlen=kmerlen,
                                          files=files, batches_dir=batches_dir, t=t)


def create_final_merged_pattern_files(outfile_patternKey, outfile_patternKeySize, outfile_patternIndex, outfile_prefix,
                                      prefix, output_dir, kmertype, kmerlen, files, batches_dir, t):
    # Check final files aren't empty
    outfiles_size = [size_of(rcompat.r_system_intern("ls -l " + x + " | cut -d ' ' -f5")) for x in (outfile_patternKey, outfile_patternKeySize, outfile_patternIndex)]
    if r_any_zero(outfiles_size):
        r_stop("One or more task_id ", t, " final ", outfile_prefix, " patternKey, patternKeySize or patternIndex files are empty")
    final_file_prefix = r_paste0(output_dir, prefix, "_", kmertype, kmerlen)
    rcompat.r_system("mv " + outfile_patternKey + " " + final_file_prefix + ".patternmerge.patternKey.txt.gz")
    rcompat.r_system("mv " + outfile_patternKeySize + " " + final_file_prefix + ".patternmerge.patternKeySize.txt")
    rcompat.r_system("mv " + outfile_patternIndex + " " + final_file_prefix + ".patternmerge.patternIndex.txt.gz")

    r_cat("Final pattern files: " + final_file_prefix + ".patternmerge.patternKey.txt.gz " + final_file_prefix
          + ".patternmerge.patternKeySize.txt " + final_file_prefix + ".patternmerge.patternIndex.txt.gz ", "\n")

    # Remove all completed files
    stem = r_paste0(batches_dir, prefix, "_", kmertype, kmerlen)
    completed_files = rcompat.r_system_intern("ls " + stem + "*.patternbatch.completed.txt")
    # Remove all intermediate pattern files
    r_cat("Removing intermediate files", "\n")
    cmd = " ".join(["rm", " ".join(completed_files)])
    rcompat.r_system2("/bin/bash", "-c '" + cmd + "'")

    pattern_files = rcompat.r_system_intern("ls " + stem + ".patternmerge.j*.patternKey.txt.gz")
    pattern_files = pattern_files + rcompat.r_system_intern("ls " + stem + ".patternmerge.j*.patternKeySize.txt")
    pattern_files = pattern_files + rcompat.r_system_intern("ls " + stem + ".patternmerge.j*.patternIndex.txt.gz")
    pattern_files = pattern_files + [key_size_name(f) for f in files]
    cmd = " ".join(["rm", " ".join(pattern_files)])
    rcompat.r_system2("/bin/bash", "-c '" + cmd + "'")


def create_kinship_batch(batch_prefix, batches_dir, nkmersbatch, patterncountspath, pattern2kinshippath, kmertype,
                         kmerlen):
    # Calculate patternCounts
    cmd = patterncountspath + " " + batches_dir + batch_prefix
    rcompat.r_system(cmd)

    # Compute kinship matrix
    cmd = pattern2kinshippath + " " + batches_dir + batch_prefix
    rcompat.r_system(cmd)

    # Check kinship file not empty
    if size_is_zero(size_of(rcompat.r_system_intern("ls -l " + batches_dir + batch_prefix + ".kinship.txt.gz | cut -d ' ' -f5"))):
        r_stop(batches_dir + batch_prefix + " kinship matrix file empty", "\n")

    # Output kinship matrix weight (i.e. total count), a plain integer (D1a)
    with open(batches_dir + batch_prefix + ".kinshipWeight.txt", "w") as f:
        f.write("%d\n" % (nkmersbatch[1] - nkmersbatch[0] + 1))
    # Check kinship weight file not empty
    if size_is_zero(size_of(rcompat.r_system_intern("ls -l " + batches_dir + batch_prefix + ".kinshipWeight.txt | cut -d ' ' -f5"))):
        r_stop(batches_dir + batch_prefix + " kinship weight file empty", "\n")

    r_cat("Created kinship matrix batch:", batches_dir + batch_prefix + ".kinship.txt.gz", "\n")
    # Write output file prefix to output files
    outfile_completed = batches_dir + batch_prefix + ".kinship.completed.txt"
    rcompat.r_system2("/bin/bash", "-c 'touch " + outfile_completed + "'")


def scan_doubles(path, nlines):
    """scan(path, what = double(0), nlines =, quiet = TRUE): numbers on the first
    nlines lines (blank lines skipped); "NA" is NA."""
    vals, k = [], 0
    with rcompat.r_open(path) as f:
        for line in f:
            if k >= nlines:
                break
            toks = line.split()
            if not toks:
                continue
            k += 1
            for tok in toks:
                v = rcompat.r_as_numeric(tok)
                if v is None and tok != "NA":
                    r_stop("scan() expected 'a real', got '", tok, "'")
                vals.append(math.nan if v is None else v)
    return vals


def scan_integers(path):
    """scan(path, what = integer(0), quiet = TRUE)."""
    with rcompat.r_open(path) as f:
        toks = f.read().split()
    out = []
    for tok in toks:
        if tok == "NA":
            out.append(None)
            continue
        try:
            out.append(rcompat.parse_index(tok))  # also reads "1e+05" from older runs (D1a)
        except ValueError:
            r_stop("scan() expected 'an integer', got '", tok, "'")
    return out


def read_kinship(path, nsamp):
    """matrix(scan(path, what = double(0), nlines = nsamp), nsamp, byrow = TRUE)."""
    vals = scan_doubles(path, nsamp)
    if len(vals) % nsamp != 0:
        r_stop("data length [", len(vals), "] is not a sub-multiple or multiple of the number of rows [", nsamp, "]")
    return np.array(vals, dtype=float).reshape(nsamp, len(vals) // nsamp)


def write_kinship(kinship, weight, kinship_file, weight_file):
    """write(sprintf("%.17g", kinship), gzfile(...), ncol = nsamp): column-major
    order, nsamp values per line; and the weight with sprintf("%d")."""
    import gzip
    nsamp = kinship.shape[0]
    vals = ["%.17g" % v for v in kinship.flatten(order="F")]
    text = "".join(" ".join(vals[k:k + nsamp]) + "\n" for k in range(0, len(vals), nsamp))
    with gzip.open(kinship_file, "wt") as f:
        f.write(text)
    with open(weight_file, "w") as f:
        f.write("NA\n" if weight is None else "%d\n" % weight)


def add_weight(total, w):
    """Integer addition as R: NA on overflow."""
    if total is None or w is None:
        return None
    s = total + w
    return s if abs(s) <= 2147483647 else None


def merge_kinship_matrices(t, n, b, p, imax, files, prefix, nkmersbatch, output_dir, kmertype, kmerlen, batches_dir):
    output_prefix_batches = [r_paste0(batches_dir, prefix, "_", kmertype, kmerlen, ".", b0, "-", b1)
                             for b0, b1 in nkmersbatch]

    # Merge
    i = 0.0
    nsamp = 0
    while True:
        i = i + 1
        if (n == 1 and i > 1) or not ((t % b ** int(i - 1)) == 0 or (t == p and i <= imax)):
            break
        r_cat("t=", t, "i=", i, "\n")
        if i == 1:
            # First round: merge initial kinship matrices
            kinship_weight = 0
            # Define input/output files
            outfile_prefix = r_paste0(batches_dir, prefix, "_", kmertype, kmerlen, ".kinshipmerge.j.", i, ".", t)
            beg = b * (t - 1) + 1
            end = min(b * t, n)
            if end < beg:
                r_stop("Problem with input arguments, please check")
            idx = rcompat.r_colon(int(beg), end)
            infiles = ["NA" if f is None else f for f in rcompat.r_index(files, idx)]
            infiles_completed = ["NA" if f is None else f + ".kinship.completed.txt"
                                 for f in rcompat.r_index(output_prefix_batches, idx)]
            nattempts = 0
            # Wait for every input and its completion marker, then check the sizes (D1b: R
            # measured the sizes once, before waiting, so a file that appeared late stopped the run)
            while not all(os.path.exists(f) for f in infiles) or not all(os.path.exists(f) for f in infiles_completed):
                nattempts = nattempts + 1
                if nattempts > 100:
                    r_stop("Could not find files", "".join(infiles))
                time.sleep(60)
            infiles_size = [size_of(rcompat.r_system_intern("ls -l " + x + " | cut -d ' ' -f5")) for x in infiles]
            if any(s is None or s != s or not s > 0 for s in infiles_size):
                r_stop("One or more kinship matrix files are empty: ", " ".join(infiles), "\n")
            # Read one line from the first kinship matrix
            kin1 = scan_doubles(infiles[0], 1)
            nsamp = len(kin1)
            # Allocate memory
            kinship = np.zeros((nsamp, nsamp))
            # Read each file and augment the kinship matrix
            for f in infiles:
                kinship_j = read_kinship(f, nsamp)
                if kinship_j.shape[1] != nsamp:
                    r_stop("Expected ", nsamp, " columns in ", f)
                if np.isnan(kinship_j).any():
                    r_stop("Found NAs in ", f)
                weightfile = f.replace(".kinship.txt.gz", "") + ".kinshipWeight.txt"
                kinship_weight_j = scan_integers(weightfile)
                if len(kinship_weight_j) != 1:
                    r_stop("Expected one entry in kinship weight file ", weightfile)
                w = kinship_weight_j[0]
                # Augment the kinship matrix and its weights
                kinship = kinship + (math.nan if w is None else float(w)) * kinship_j
                kinship_weight = add_weight(kinship_weight, w)
            # Save the kinship matrix and its weight
            kinship = kinship / (math.nan if kinship_weight is None else float(kinship_weight))
            write_kinship(kinship, kinship_weight, outfile_prefix + ".kinship.txt.gz", outfile_prefix + ".kinshipWeight.txt")
            rcompat.r_system2("/bin/bash", "-c 'touch " + outfile_prefix + ".kinship.completed.txt'")
        else:
            # Subsequent rounds: merge merged kinship matrices
            kinship = np.zeros((nsamp, nsamp))
            kinship_weight = 0
            # Define input/output files
            te = math.ceil(t / (b ** (i - 1))) * int(b ** (i - 1))
            outfile_prefix = r_paste0(batches_dir, prefix, "_", kmertype, kmerlen, ".kinshipmerge.j.", i, ".", te)
            outfile_kinship = outfile_prefix + ".kinship.txt.gz"
            outfile_kinshipWeight = outfile_prefix + ".kinshipWeight.txt"
            outfile_completed = outfile_prefix + ".kinship.completed.txt"
            beg = te - float(b) ** (i - 1) + float(b) ** (i - 2)
            end = min(te, math.ceil(t / b ** (i - 2)) * int(b ** (i - 2)))
            if end < beg:
                r_stop("Problem with input arguments, please check")
            inc = float(b) ** (i - 2)
            js = rcompat.r_seq(beg, end, inc)
            stem = r_paste0(batches_dir, prefix, "_", kmertype, kmerlen, ".kinshipmerge.j.", i - 1, ".")
            infiles = [stem + rcompat.r_as_character(x) + ".kinship.txt.gz" for x in js]
            infiles_weights = [stem + rcompat.r_as_character(x) + ".kinshipWeight.txt" for x in js]
            infiles_completed = [stem + rcompat.r_as_character(x) + ".kinship.completed.txt" for x in js]
            nattempts = 0
            while (not all(os.path.exists(f) for f in infiles) or not all(os.path.exists(f) for f in infiles_weights)
                   or not all(os.path.exists(f) for f in infiles_completed)):
                nattempts = nattempts + 1
                if nattempts > 3600:
                    r_stop("Could not find files", "".join(infiles))
                time.sleep(1)
            infiles_size = [size_of(rcompat.r_system_intern("ls -l " + x + " | cut -d ' ' -f5")) for x in infiles]
            if r_any_zero(infiles_size):
                r_stop("One or more file size is zero ", " ".join(infiles), "\n")
            infiles_weights_size = [size_of(rcompat.r_system_intern("ls -l " + x + " | cut -d ' ' -f5")) for x in infiles_weights]
            if r_any_zero(infiles_weights_size):
                r_stop("One or more file size is zero ", " ".join(infiles_weights), "\n")
            # Read each file and augment the kinship matrix
            for f, wf in zip(infiles, infiles_weights):
                kinship_j = read_kinship(f, nsamp)
                if kinship_j.shape[1] != nsamp:
                    r_stop("Expected ", nsamp, " columns in ", f)
                if np.isnan(kinship_j).any():
                    r_stop("Found NAs in ", f)
                kinship_weight_j = scan_integers(wf)
                if len(kinship_weight_j) != 1:
                    r_stop("Expected one entry in kinship weight file ", wf)
                w = kinship_weight_j[0]
                # Augment the kinship matrix and its weights
                kinship = kinship + (math.nan if w is None else float(w)) * kinship_j
                kinship_weight = add_weight(kinship_weight, w)
            # Save the kinship matrix and its weight
            kinship = kinship / (math.nan if kinship_weight is None else float(kinship_weight))
            write_kinship(kinship, kinship_weight, outfile_kinship, outfile_kinshipWeight)
            rcompat.r_system2("/bin/bash", "-c 'touch " + outfile_completed + "'")

    if t == p:
        # Define input/output files: bug fix DJW 20220522
        i = i - 1
        te = math.ceil(t / (b ** (i - 1))) * int(b ** (i - 1))
        outfile_prefix = r_paste0(batches_dir, prefix, "_", kmertype, kmerlen, ".kinshipmerge.j.", i, ".", te)
        outfile_kinship = outfile_prefix + ".kinship.txt.gz"
        outfile_kinshipWeight = outfile_prefix + ".kinshipWeight.txt"
        create_final_kinship_file(outfile_kinship=outfile_kinship, outfile_kinshipWeight=outfile_kinshipWeight, t=t,
                                  outfile_prefix=outfile_prefix, prefix=prefix, output_dir=output_dir,
                                  kmertype=kmertype, kmerlen=kmerlen, files=files, batches_dir=batches_dir)


def create_final_kinship_file(outfile_kinship, outfile_kinshipWeight, t, outfile_prefix, output_dir, prefix, kmertype,
                              kmerlen, files, batches_dir):
    outfiles_size = [size_of(rcompat.r_system_intern("ls -l " + x + " | cut -d ' ' -f5")) for x in (outfile_kinship, outfile_kinshipWeight)]
    if r_any_zero(outfiles_size):
        r_stop("One or more task_id ", t, " ", outfile_prefix, " kinship or kinship weight files are empty")
    final_file_prefix = r_paste0(output_dir, prefix, "_", kmertype, kmerlen)
    rcompat.r_system("mv " + outfile_kinship + " " + final_file_prefix + ".kinshipmerge.kinship.txt.gz")
    rcompat.r_system("mv " + outfile_kinshipWeight + " " + final_file_prefix + ".kinshipmerge.kinshipWeight.txt")

    r_cat("Final kinship files: " + final_file_prefix + ".kinshipmerge.kinship.txt.gz " + final_file_prefix
          + ".kinshipmerge.kinshipWeight.txt", "\n")

    # Remove all completed files
    stem = r_paste0(batches_dir, prefix, "_", kmertype, kmerlen)
    completed_files = rcompat.r_system_intern("ls " + stem + "*.kinship.completed.txt")
    cmd = " ".join(["rm", " ".join(completed_files)])
    rcompat.r_system2("/bin/bash", "-c '" + cmd + "'")

    # Remove all intermediate kinship files
    kinship_files = rcompat.r_system_intern("ls " + stem + ".kinshipmerge.j*.kinship.txt.gz")
    kinship_files = kinship_files + rcompat.r_system_intern("ls " + stem + ".kinshipmerge.j*.kinshipWeight.txt")
    # Remove those in pattern subdirectory
    kinship_files = kinship_files + list(files)
    kinship_files = kinship_files + [f.replace(".kinship.txt.gz", ".patternCounts.txt.gz") for f in files]
    kinship_files = kinship_files + [f.replace(".kinship.txt.gz", ".kinshipWeight.txt") for f in files]
    r_cat("Removing intermediate files", "\n")
    cmd = " ".join(["rm", " ".join(kinship_files)])
    rcompat.r_system2("/bin/bash", "-c '" + cmd + "'")


###################################################################################################


def main():
    rcompat.script_setup(__file__)
    start_time = time.monotonic()
    parser = argparse.ArgumentParser(
        description="stringlist2patternandkinship.py define presence/absence of a kmer list among kmer count files "
                    "and build kinship matrix. Adapted from Daniel Wilson scripts (2018) kmerlist2pattern.Rscript, "
                    "patternmerge.Rscript, pattern2kinship.Rscript, kinshipmerge.Rscript", allow_abbrev=False)
    parser.add_argument("--task-id", required=True)
    parser.add_argument("--p", required=True, help="number of parallel tasks")
    parser.add_argument("--id-file", required=True)
    parser.add_argument("--fullkmerlistfile", required=True, help="merged list of all k-mers (gzipped)")
    parser.add_argument("--kmercountslistfile", required=True, help="file listing the per-sample k-mer count files")
    parser.add_argument("--analysis-dir", required=True)
    parser.add_argument("--output-prefix", required=True)
    parser.add_argument("--kmertype", required=True, help="protein or nucleotide")
    parser.add_argument("--software-file", required=True)
    parser.add_argument("--kmer-length", default="31")
    parser.add_argument("--mincount", default="5")
    args = parser.parse_args()

    # Initialize variables
    t = rcompat.r_as_integer(args.task_id)
    p = rcompat.r_as_integer(args.p)
    id_file = args.id_file
    fullkmerlistfile = args.fullkmerlistfile
    kmercountslistfile = args.kmercountslistfile
    output_dir = args.analysis_dir
    output_prefix = args.output_prefix
    kmertype = args.kmertype.lower()
    software_file = args.software_file
    kmerlen = rcompat.r_as_integer(args.kmer_length)
    mincount = rcompat.r_as_integer(args.mincount)

    # Check input arguments
    if p is None:
        r_stop("Error: p must be an integer", "\n")
    if not os.path.exists(id_file):
        r_stop("Error: sample ID file doesn't exist", "\n")
    if not os.path.exists(fullkmerlistfile):
        r_stop("Error: full kmer list file doesn't exist", "\n")
    if not os.path.exists(kmercountslistfile):
        r_stop("Error: kmer counts list file doesn't exist", "\n")
    if not os.path.exists(output_dir):
        r_stop("Error: output directory doesn't exist", "\n")
    if not output_dir.endswith("/"):
        output_dir = output_dir + "/"
    if kmertype != "protein" and kmertype != "nucleotide":
        r_stop("Error: kmer type must be either 'protein' or 'nucleotide'", "\n")
    if kmerlen is None:
        r_stop("Error: kmer length must be an integer", "\n")
    if mincount is None:
        r_stop("Error: min count must be an integer", "\n")

    # Read in software file
    software_paths = rcompat.r_read_table(software_file, header=True, sep="\t", quote="")
    names = [rcompat.r_as_character(v) for v in software_paths["name"]]
    paths = [rcompat.r_as_character(v) for v in software_paths["path"]]
    # Children are Python scripts run with this interpreter; software files keep the "R" entry
    required_software = ["R", "scriptpath"]
    if any(r not in names for r in required_software):
        r_stop("Error: missing required software path in the software file - requires " + ", ".join(required_software), "\n")
    python_path = sys.executable

    script_location = [pth for nm, pth in zip(names, paths) if nm.lower() == "scriptpath"][0]
    if not os.path.isdir(script_location):
        r_stop("Error: script location directory specified in the software paths file doesn't exist", "\n")

    def tool(name, what):
        path = script_location + "/" + name
        if not os.path.exists(path):
            r_stop("Error: " + what + " path doesn't exist - check pipeline script location in the software file", "\n")
        return path
    stringlist2patternpath = tool("stringlist2pattern", "stringlist2pattern")
    kmerlist2patternpath = tool("kmerlist2pattern", "kmerlist2pattern")
    patternmergepath = tool("patternmerge", "patternmerge")
    patterncountspath = tool("patterncounts", "patterncounts")
    pattern2kinshippath = tool("pattern2kinship", "pattern2kinship")
    pattern2presencecountscript = tool("pattern2presencecount.py", "pattern2presencecount.py")

    # Report variables
    r_cat("#############################################", "\n")
    r_cat("Running on host: ", rcompat.r_system_intern("hostname"), "\n")
    r_cat("Command line arguments", "\n")
    r_cat(sys.argv[1:], "\n\n")
    r_cat("Parameters:", "\n")
    r_cat("task_id:", t, "\n")
    r_cat("p:", p, "\n")
    r_cat("ID file path:", id_file, "\n")
    r_cat("Full kmer list file:", fullkmerlistfile, "\n")
    r_cat("Kmer counts list file:", kmercountslistfile, "\n")
    r_cat("Analysis directory:", output_dir, "\n")
    r_cat("Output prefix:", output_prefix, "\n")
    r_cat("Kmer type:", kmertype, "\n")
    r_cat("Software file:", software_file, "\n")
    r_cat("Script location:", script_location, "\n")
    r_cat("Python path:", python_path, "\n")
    r_cat("Kmer length:", kmerlen, "\n")
    r_cat("Min count:", mincount, "\n")
    r_cat("#############################################", "\n\n")

    if t > p:
        r_cat("Task", t, "not required\n")
        return

    # Create pattern batch
    batch = create_pattern_batch(fullkmerlistfile=fullkmerlistfile, p=p, t=t,
                                 stringlist2patternpath=stringlist2patternpath,
                                 kmerlist2patternpath=kmerlist2patternpath, kmercountslistfile=kmercountslistfile,
                                 kmerlen=kmerlen, mincount=mincount, kmertype=kmertype, output_prefix=output_prefix,
                                 output_dir=output_dir)

    r_cat("\n")

    # Create kinship matrix for the batch
    create_kinship_batch(batch_prefix=batch["output_prefix_batch"], batches_dir=batch["batches_dir"],
                         nkmersbatch=batch["nkmersbatch"][t - 1], patterncountspath=patterncountspath,
                         pattern2kinshippath=pattern2kinshippath, kmertype=kmertype, kmerlen=kmerlen)

    r_cat("\n")

    # For a subset of the processes, use to merge the batches created by the other processes
    # Which processes to keep (floor() gives an R double)
    p = float(math.floor(p / 5))
    p = 1.0 if p == 0 else p

    if t > p:
        r_cat("Created pattern batch", t, "\n")
        r_cat("Task", t, "not required for pattern and kinship matrix merging", "\n")
        r_cat("Finished in", (time.monotonic() - start_time) / 3600, "hours\n")
        return
    else:
        r_cat("Created pattern and kinship batch", t, "\n")
        r_cat("Running kmer pattern merge", "\n")
        # Get merge parameters
        merge_parameters = merge_pattern_batches_parameters(nkmersbatch=batch["nkmersbatch"], p=p, t=t,
                                                            output_dir=output_dir, output_prefix=output_prefix,
                                                            kmertype=kmertype, kmerlen=kmerlen)

        # Merge patterns
        merge_patterns(t=t, n=merge_parameters["n"], b=merge_parameters["b"], p=p, imax=merge_parameters["imax"],
                       prefix=output_prefix, files=merge_parameters["infiles"], patternmergepath=patternmergepath,
                       output_dir=output_dir, kmertype=kmertype, kmerlen=kmerlen, batches_dir=batch["batches_dir"])
        r_cat("\n")
        r_cat("Running kinship matrix merge", "\n")
        # Merge kinship matrices
        kinship_matrix_files = [f.replace(".patternKey.txt.gz", ".kinship.txt.gz") for f in merge_parameters["infiles"]]
        merge_kinship_matrices(t=t, n=merge_parameters["n"], b=merge_parameters["b"], p=p,
                               imax=merge_parameters["imax"], files=kinship_matrix_files, prefix=output_prefix,
                               nkmersbatch=batch["nkmersbatch"], output_dir=output_dir, kmertype=kmertype,
                               kmerlen=kmerlen, batches_dir=batch["batches_dir"])
        r_cat("\n")

    # Keep the last process to calculate the MACs for the patterns
    if t == p:
        r_cat("Calculating number of genomes each pattern is present in (of the genomes with a non NA phenotype)", "\n")
        rcompat.r_system(python_path + " " + pattern2presencecountscript + " --kmerfile-prefix "
                         + r_paste0(output_dir, output_prefix, "_", kmertype, kmerlen) + " --output-dir " + output_dir
                         + " --id-file " + id_file + " --include-na FALSE")
        r_cat("\n")

    r_cat("Finished in", (time.monotonic() - start_time) / 3600, "hours\n")


if __name__ == "__main__":
    main()
