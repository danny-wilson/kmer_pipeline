#!/usr/bin/env python3
"""gen-unmapped-report.py: generate the HTML report on k-mers that did not align
to the reference. Port of gen-unmapped-report.Rscript (Daniel Wilson, 2022)."""
import argparse
import math
import os
import sys
import time

import rcompat
from rcompat import r_s3 as s3

NL = "\n"
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
gr = __import__("gen-report")
ggr = __import__("gen-gene-report")


def main():
    rcompat.script_setup(__file__)
    parser = argparse.ArgumentParser(description="gen-unmapped-report.py Generate a kmer GWAS report for unmapped kmers. "
                                                 "Daniel Wilson (2022)", allow_abbrev=False)
    for name in ("prefix", "anatype", "k", "refname", "ref-gb", "maf", "alignident", "mincount", "srcdir", "outdir",
                 "logdir"):
        parser.add_argument("--" + name, required=True)
    args = parser.parse_args()
    PREFIX, ANATYPE, K, REFNAME = args.prefix, args.anatype, args.k, args.refname
    MAF, ALIGNIDENT, MINCOUNT, PWD = args.maf, args.alignident, args.mincount, args.outdir
    MACORMAF = "maf" if gr.r_lt(MAF, "1") else "mac"
    FIGDIR = ANATYPE + "kmer" + K + "_kmergenealign_figures/"

    os.chdir(PWD)
    is_maf = True if MAF == "0" else gr.r_lt(MAF, "1")
    is_maf_text3 = "MAF" if is_maf else "MAC"

    html_head = NL.join(["<!DOCTYPE html>", "<html>", "<head>", "  <title>Kmer GWAS report: unmapped kmers</title>",
                         "  <link rel='stylesheet' href='report.css'>", NL])
    html_body = NL.join([
        "</head>", "<body>", "  <h1>Kmer GWAS report: unmapped kmers</h1>",
        "  <div><p class='timestamp'><code>Prefix: " + PREFIX + "; KmerType: " + ANATYPE + "; K:",
        "  " + K + "; ReferenceGenome: " + REFNAME + "; " + is_maf_text3 + ": " + MAF + "; MinCount:",
        "  " + MINCOUNT + "; AlignIdent: " + ALIGNIDENT + "; ReportTimeStamp:",
        "  " + time.ctime() + ".</code></p></div>", NL])
    html_foot = NL.join(["<script src='report.js'></script>", "</body>", "</html>", ""])

    outfile_html = PREFIX + "_" + ANATYPE + K + "." + "report_unmapped.html"
    filename_unmapped = (FIGDIR + PREFIX + "_" + ANATYPE + K + "_" + REFNAME + "_" + MACORMAF + "_" + MAF + "_alignIdent_"
                         + ALIGNIDENT + "_alignPosMinCount_" + MINCOUNT + "_unaligned_kmersandpvals.txt")
    rows = ggr.read_delim(filename_unmapped)[1]
    neg = [r["negLog10"] for r in rows]
    max_signif_unmapped = max(neg) if neg else -math.inf
    summary = gr.read_summary_json(PREFIX + "_" + ANATYPE + K + ".summary.json")
    thr_signif = float(summary["bonferroni_threshold"])

    html_body = NL.join([
        html_body,
        "  <p>Not all kmers were aligned to the user-supplied reference genome. The strongest",
        "  significance among the unmapped kmers was " + s3(max_signif_unmapped) + ","])
    if max_signif_unmapped >= thr_signif:
        html_body = NL.join([html_body, "  which was genome-wide significant.</p>", "", NL])
        tb = [{"kmer": r["kmer"], "Signif": s3(r["negLog10"]), "beta": s3(r["beta"]), "MAC": r["mac"]}
              for r in rows if r["negLog10"] >= thr_signif]
        lines = [html_body, "  <h2>Significant unmapped kmers</h2>", "  <div class='divkmertab'>",
                 "  <table class='kmertab center'>", "    <tr>",
                 "      <th>" + "</th><th>".join(["kmer", "Signif", "beta", "MAC"]) + "</th>", "    </tr>"]
        for r in tb:
            lines += ["    <tr>", "      <td>" + "</td><td>".join(ggr.rstr(r[c]) for c in ("kmer", "Signif", "beta", "MAC"))
                      + "</tr>", "    </tr>"]
        lines += ["  </table>", "  </div>", NL]
        html_body = NL.join(lines)
    else:
        html_body = NL.join([html_body, "  which was not genome-wide significant.</p>", NL])

    html_body = NL.join([
        html_body,
        "  <p>For a summary of the associations at all unmapped kmers,",
        "  follow the link to <a href=" + filename_unmapped + " target='_blank' rel='noopener noreferrer'>this file</a>.</p>",
        NL])
    with open(outfile_html, "w") as f:
        f.write(html_head + " " + html_body + " " + html_foot)


if __name__ == "__main__":
    main()
