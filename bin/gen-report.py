#!/usr/bin/env python3
"""gen-report.py: generate the k-mer GWAS HTML report. Port of gen-report.Rscript
(Daniel Wilson, 2022). The HTML is built line by line as the R script builds it,
so the report is the same apart from its time stamp."""
import argparse
import math
import os
import sys
import time

import shutil

import rcompat
from rcompat import r_s3 as s3

NL = "\n"


def r_lt(a, b):
    """a < b for two R strings (collation order, as R compares character values)."""
    return rcompat.r_collate_key(a) < rcompat.r_collate_key(b)


def read_summary_json(filename):
    """The flat JSON object written by plotManhattan, as R's read_summary_json: values as text."""
    import re
    out = {}
    for l in open(filename).read().split("\n"):
        if ":" not in l:
            continue
        key = re.sub(r'^\s*"([^"]+)".*$', r"\1", l)
        val = re.sub(r",\s*$", "", re.sub(r"^[^:]*:\s*", "", l, count=1), count=1).replace('"', "")
        out[key] = val
    return out


def zgrep_last(regex, filename):
    """The last word of the lines zgrep finds, as a number (NA on failure)."""
    try:
        text = rcompat.r_pipe("zgrep '" + regex + "' " + filename)
        words = [w for line in text.split("\n")[:-1] for w in line.split(" ")]
        return rcompat.r_as_numeric(words[-1]) if words else None
    except Exception:  # noqa: BLE001  (R: tryCatch(..., error = function(e) NA))
        return None


def main():
    rcompat.script_setup(__file__)
    parser = argparse.ArgumentParser(description="gen-report.py Generate a kmer GWAS report. Daniel Wilson (2022)",
                                     allow_abbrev=False)
    for name in ("prefix", "anatype", "k", "refname", "ref-gb", "maf", "alignident", "mincount", "ngenes", "srcdir",
                 "outdir", "logdir"):
        if name == "mincount":  # D4: --plot-min-genomes; --mincount kept as an alias
            parser.add_argument("--plot-min-genomes", "--mincount", dest="mincount", required=True,
                                help="genomes a k-mer/gene combination must be seen in to be plotted (as step 6)")
        else:
            parser.add_argument("--" + name, required=True)
    args = parser.parse_args()

    PREFIX, ANATYPE, K, REFNAME, REF_GB = args.prefix, args.anatype, args.k, args.refname, args.ref_gb
    MAF, ALIGNIDENT, MINCOUNT, NGENES = args.maf, args.alignident, args.mincount, args.ngenes
    SRC, PWD, LOGDIR = args.srcdir, args.outdir, args.logdir
    MACORMAF = "maf" if r_lt(MAF, "1") else "mac"  # R compares the text of MAF with 1

    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import sequence_functions

    # MAF vs MAC text
    is_maf = True if MAF == "0" else r_lt(MAF, "1")
    is_maf_text = ("minor allele frequency (MAF) threshold of " + MAF) if is_maf else \
        ("minor allele count (MAC) threshold of " + MAF)
    is_maf_text2 = "minor allele frequency (MAF)" if is_maf else "minor allele count (MAC)"
    is_maf_text3 = "MAF" if is_maf else "MAC"

    # Copy (overwrite if necessary) the CSS stylesheet and JavaScript code
    for f in ("report.css", "report.js"):
        shutil.copyfile(SRC + "/" + f, PWD + "/" + f)

    outfile_prefix = PREFIX + "_" + ANATYPE + K + "."
    outfile_html = outfile_prefix + "report.html"
    FIGDIR = ANATYPE + "kmer" + K + "_kmergenealign_figures/"
    stem = PREFIX + "_" + ANATYPE + K + "_" + REFNAME
    filename_topgenes = (FIGDIR + stem + "_top20genes_toppvals_" + MACORMAF + "_" + MAF + "_nucmerAlign_alignIdent_"
                         + ALIGNIDENT + "_alignPosMinCount_" + MINCOUNT + ".txt")
    filename_heritability = ANATYPE + "kmer" + K + "_gemma/output/" + PREFIX + "_" + ANATYPE + K + ".1-*.log.txt.gz"
    man = FIGDIR + stem + "_LMM_kmergenealign_ct" + MINCOUNT + "_Manhattan_"
    filename_Manhattan_maf = man + "mafCOL_" + MACORMAF + MAF + ".png"
    filename_Manhattan_beta = man + "betaCOL_" + MACORMAF + MAF + ".png"
    filename_Manhattan_align = man + "alignCOL_" + MACORMAF + MAF + ".png"
    filename_Manhattan_maf0 = man + "mafCOL_" + MACORMAF + "0.png"
    filename_qq_maf = FIGDIR + PREFIX + "_" + ANATYPE + K + "_QQplot_" + MACORMAF + MAF + ".png"
    filename_qq_maf0 = FIGDIR + PREFIX + "_" + ANATYPE + K + "_QQplot_allkmers.png"
    filename_unmapped = (FIGDIR + stem + "_" + MACORMAF + "_" + MAF + "_alignIdent_" + ALIGNIDENT + "_alignPosMinCount_"
                         + MINCOUNT + "_unaligned_kmersandpvals.txt")

    # Preliminaries
    os.chdir(PWD)

    html_head = NL.join(["<!DOCTYPE html>", "<html>", "<head>", "  <title>Kmer GWAS report</title>",
                         "  <link rel='stylesheet' href='report.css'>", NL])
    html_body = NL.join([
        "</head>", "<body>", "  <h1>Kmer GWAS report</h1>",
        "  <div><p class='timestamp'><code>Prefix: " + PREFIX + "; KmerType: " + ANATYPE + "; K:",
        "  " + K + "; ReferenceGenome: " + REFNAME + "; " + is_maf_text3 + ": " + MAF + "; MinCount:",
        "  " + MINCOUNT + "; AlignIdent: " + ALIGNIDENT + "; ReportTimeStamp:",
        "  " + time.ctime() + ".</code></p></div>", NL])
    html_foot = NL.join(["<script src='report.js'></script>", "</body>", "</html>", ""])

    # Genbank file
    gbk = sequence_functions.read_dna_seg_from_file(REF_GB, tagsToParse=("CDS",))
    gbk_names, gbk_product = list(gbk["name"]), list(gbk["product"])
    with rcompat.r_open(REF_GB) as f:
        first = f.readline().rstrip("\n")
    print("Read 1 item", file=sys.stderr, flush=True)
    toks = [t for t in first.split(" ") if t != ""]
    if len(toks) < 3 or rcompat.r_as_numeric(toks[2]) is None:
        rcompat.r_stop("Error retrieving the reference genome length from the genbank file", "\n")

    # Heritability
    d = os.path.dirname(filename_heritability)
    hits = rcompat.r_dir(d, glob=os.path.basename(filename_heritability), full_names=True) if os.path.isdir(d) else []
    filename_heritability_1 = hits[0] if hits else "NA"
    pve = zgrep_last("pve estimate in the null (linear mixed) model", filename_heritability_1)
    se_pve = zgrep_last("se(pve) in the null (linear mixed) model", filename_heritability_1)

    def qnorm(p, mean, sd):
        if mean is None or sd is None:
            return None  # NA
        if math.isnan(mean) or math.isnan(sd):
            return math.nan
        from scipy.stats import norm
        return float(norm.ppf(p, mean, sd))
    lo = qnorm(0.025, pve, se_pve)
    hi = qnorm(0.975, pve, se_pve)
    html_body = NL.join([
        html_body,
        "  <h2>Heritability</h2>",
        "  <p>The sample heritability (proportion of variance explained) under the null",
        "  linear mixed model (LMM) was " + s3(pve) + " with a standard error of " + s3(se_pve) + ",",
        "  which implies a 95% confidence interval of",
        "  (" + s3(None if lo is None else (lo if math.isnan(lo) else max(0.0, lo))) + ", "
        + s3(None if hi is None else (hi if math.isnan(hi) else min(1.0, hi))) + ").",
        "  </p>", NL])

    # Obtain the actual p-value threshold used, after filtering samples with
    # no phenotypes and applying the MAF or MAC filter
    summary = read_summary_json(PREFIX + "_" + ANATYPE + K + ".summary.json")
    thr_signif = float(summary["bonferroni_threshold"])
    total_impliedtests = int(float(summary["n_tests"]))
    thr_p = 0.05 / total_impliedtests
    total_npatterns = int(float(summary["n_patterns"]))
    total_nkmers = int(float(summary["n_kmers"]))

    html_body = NL.join([
        html_body,
        "  <h2>Significance threshold</h2>",
        "  <p>A total of " + str(total_nkmers) + " distinct kmers were observed, of which there",
        "  were " + str(total_npatterns) + " unique phylopatterns (patterns of presence or absence)",
        "  across the sample.",
        "  After filtering any individuals lacking phenotype information, and",
        "  applying a " + is_maf_text + ", there were " + str(total_impliedtests) + " " +
        "  unique phylopatterns to be tested. Assuming a",
        "  familywise error rate of 5%, this implied a Bonferroni-corrected",
        "  <i>p</i>-value threshold of " + s3(thr_p) + ", or 10<sup>-" + s3(thr_signif) + "</sup>.</p>", NL])

    # N9: the genomes and patterns the results describe (summaries from earlier releases lack them)
    if "n_genomes_analysed" in summary:
        n_untested = int(float(summary["n_untested_patterns"]))
        n_nan = int(float(summary["n_patterns_nan"]))
        groups = ""
        if summary.get("pheno_type") == "binary" and "n_cases" in summary:
            groups = (": " + summary["n_cases"] + " with phenotype " + summary["case_value"] + " and "
                      + summary["n_controls"] + " with phenotype " + summary["control_value"])
        html_body = NL.join([
            html_body,
            "  <h2>Genomes and patterns analysed</h2>",
            "  <p>Of the " + summary["n_genomes"] + " genomes, " + summary["n_genomes_analysed"] + " were analysed"
            + groups + ". Genomes are analysed if they have a phenotype and, when covariates are used, a value for"
            + " every covariate.",
            "  GEMMA tested " + str(total_npatterns - n_untested - n_nan) + " of the " + str(total_npatterns)
            + " phylopatterns; " + str(n_untested) + " were not tested because they do not vary among the analysed"
            + " genomes" + (", and " + str(n_nan) + " could not be fitted (no result)" if n_nan else "") + ".</p>",
            NL])

    # Top regions by min p-value
    table_topgenes = rcompat.r_read_table(filename_topgenes)
    top_names = [rcompat.r_as_character(v) for v in table_topgenes.iloc[:, 0]]
    top_p = [rcompat.r_as_numeric_value(v) for v in table_topgenes.iloc[:, 1]]
    ngenes_signif = sum(1 for v in top_p if v is not None and v >= thr_signif)

    first_index = {}
    records = list(gbk["record"]) if "record" in gbk.columns else None
    for k, n in enumerate(gbk_names):
        first_index.setdefault(n, k)
        if records is not None:  # D6: names used in several records appear as name@record
            first_index.setdefault(n + "@" + records[k], k)
    import reference
    single_record = reference.records(REF_GB)[0].name

    def record(region):
        """The reference record of a gene or intergenic region (its first gene's)."""
        if records is None:
            return single_record
        k = first_index.get(region.split(":")[0])
        return "NA" if k is None else records[k]

    def product(n):
        k = first_index.get(n)
        return "NA" if k is None or gbk_product[k] is None else gbk_product[k]
    propro = [" :<br> ".join(product(g) for g in sequence_functions.r_strsplit(genes, ":")) for genes in top_names]

    lines = [html_body,
             "  <h2>Most significant regions</h2>",
             "  <p>The " + NGENES + " most significant genes or intergenic regions are summarized",
             "  in the Table below. Of those, " + str(ngenes_signif) + " were",
             "  genome-wide significant. The gene or (if an intergenic region)",
             "  flanking genes are named for each region, alongside its significance.",
             "  In what follows, <i>significance</i> is defined as the -log<sub>10</sub>",
             "  <i>p</i>-value. The significance of each region was based on the",
             "  smallest <i>p</i>-value in that region.</p>",
             NL,
             "  <table>",
             "    <tr>",
             "      <th>Region</th><th>Record</th><th>Significance</th><th>Product</th>",
             "    </tr>"]
    for g, p, prod in zip(top_names, top_p, propro):
        tag = ("<b>", "</b>") if not (p < thr_signif) else ("", "")
        report_filename = outfile_prefix + "report_" + g.replace(":", "_") + ".html"
        lines += ["    <tr>",
                  "      <td><a href='" + report_filename + "' target='_blank' rel='noopener noreferrer'>" + g
                  + "</a></td><td>" + record(g) + "</td><td>" + tag[0] + s3(p) + tag[1] + "</td><td>" + prod + "</td>",
                  "    </tr>"]
    lines += ["  </table>", NL]
    html_body = NL.join(lines)

    # Unmapped kmers
    table_unmapped = rcompat.r_read_table(filename_unmapped, header=True, sep="\t", quote="\"", comment_char="")
    table_unmapped = table_unmapped.iloc[:1]
    if len(table_unmapped) > 0:
        max_signif_unmapped = max(rcompat.r_as_numeric_value(v) for v in table_unmapped["negLog10"])
        html_body = NL.join([
            html_body,
            "  <h2>Unmapped kmers</h2>",
            "  <p>Not all kmers were aligned to the user-supplied reference genome. The strongest",
            "  significance among the unmapped kmers was " + s3(max_signif_unmapped) + ","])
        if max_signif_unmapped >= thr_signif:
            html_body = NL.join([html_body, "  which was genome-wide significant.</p>"])
        else:
            html_body = NL.join([html_body, "  which was not genome-wide significant.</p>"])
        html_body = NL.join([
            html_body,
            "  <p>Follow this link for the <a href='" + outfile_prefix + "report_unmapped.html' target='_blank' "
            "rel='noopener noreferrer'>Report on unmapped kmers</a>.</p>", NL])

    # Manhattan plot
    html_body = NL.join([
        html_body,
        "  <h2>Manhattan plot</h2>",
        "  <p>The Figure displays the significance of each kmer against the position in the",
        "  reference genome to which it mapped. Kmers that did not map are shown at the far",
        "  right hand side. The Bonferroni-corrected significance threshold is shown as a horizontal black dashed line.",
        "  The names of significant regions are plotted above. Points are colour-coded in",
        "  an adjustable manner to display " + is_maf_text2 + ", <i>&beta;</i> (direction of effect) or uniqueness",
        "  of mapping. The " + is_maf_text3 + " threshold can also be removed (although the significance",
        "  threshold is not updated since we do not recommend reporting low-" + is_maf_text3 + "  kmers",
        "  as significant).</p>", NL])

    filenames = [filename_Manhattan_maf, filename_Manhattan_beta, filename_Manhattan_align, filename_Manhattan_maf0]
    descriptions = ["Kmers colour-coded by minor allele frequency.",
                    "Kmers colour-coded by direction of effect. When <i>&beta;</i>&nbsp;>&nbsp;0, the presence of the "
                    "kmer is associated with larger values of the phenotype.",
                    "Kmers colour-coded by mapping uniqueness.",
                    "Kmers colour-coded by minor allele frequency. No " + is_maf_text3 + " filter."]
    html_body = slideshow(html_body, filenames, descriptions, "")
    html_body = NL.join([html_body, NL])

    # QQ plots
    filenames = [filename_qq_maf, filename_qq_maf0]
    descriptions = ["QQ plot with " + MACORMAF + " filter of " + MAF + ".", "QQ plot with no " + MACORMAF + " filter."]
    html_body = NL.join([
        html_body,
        "  <h2>QQ plots</h2>",
        "  <p>The QQ plots in the Figure below allow an assessment of whether",
        "  there were any problems with inflation of significance in the analysis.",
        "  Inflation is detected by an elevation of the black solid line above",
        "  the red dashed line at relatively small -log<sub>10</sub> <i>p</i>-values.",
        "  An elevation of the black solid line above the red dashed line only",
        "  at relatively large values (e.g. above the significance threshold) is evidence of association, rather",
        "  than inflation.</p>",
        "  <p>If the black solid line falls below the red dashed line, that may",
        "  provide evidence of deflation, which occurs when the analysis is",
        "  under-powered. The removal of low " + is_maf_text3 + " variants is one measure aimed",
        "  at avoiding deflation by avoiding under-powered tests. Note that the QQ plot",
        "  is noisier at larger -log<sub>10</sub> <i>p</i>-values.",
        ""])
    # In R the expression adding the closing tags continues into the cat() below, which
    # writes the file before they are added: the report ends without them
    html_body = slideshow(html_body, filenames, descriptions, '  style="width:60%"', close=False)

    # Final: cat(head, body, foot) separates them with spaces
    with open(outfile_html, "w") as f:
        f.write(html_head + " " + html_body + " " + html_foot)

    if False:  # as in R (if(FALSE)): the other reports are run by Nextflow instead
        print("Generating top hit gene reports:")

        def fp(s):
            return s.replace(" ", "\\ ")
        for i in range(1, int(NGENES) + 1):
            CMD = " ".join([sys.executable, fp(SRC + "/gen-gene-report.py"), "--hit-num", str(i), "--prefix", PREFIX,
                            "--anatype", ANATYPE, "--k", K, "--refname", REFNAME, "--ref-gb", fp(REF_GB), "--maf", MAF,
                            "--alignident", ALIGNIDENT, "--plot-min-genomes", MINCOUNT, "--srcdir", fp(SRC), "--outdir", fp(PWD),
                            "--logdir", fp(LOGDIR)])
            rcompat.r_system(CMD)
            print("Done", i, "of", NGENES)
        # And for unmapped kmers
        CMD = " ".join([sys.executable, fp(SRC + "/gen-unmapped-report.py"), "--prefix", PREFIX, "--anatype", ANATYPE,
                        "--k", K, "--refname", REFNAME, "--ref-gb", fp(REF_GB), "--maf", MAF, "--alignident", ALIGNIDENT,
                        "--plot-min-genomes", MINCOUNT, "--srcdir", fp(SRC), "--outdir", fp(PWD), "--logdir", fp(LOGDIR)])
        rcompat.r_system(CMD)
        print("Done unmapped kmer report")


def slideshow(html_body, filenames, descriptions, img_style, close=True):
    """The slide-show block (R builds it inline in gen-report and the gene reports)."""
    n = len(filenames)
    lines = [html_body, '  <div class="slideshow-container">']
    for i, (f, d) in enumerate(zip(filenames, descriptions), start=1):
        lines += ['    <div class="mySlides fade">',
                  '      <div class="numbertext">' + str(i) + ' / ' + str(n) + '</div>',
                  '      <img src="' + f + '" class="center"' + img_style + '>',
                  '      <div class="text">' + d + '</div>',
                  '    </div>']
    lines += ['    <a class="prev" onclick="plusSlides(-1,this)">&#10094;</a>',
              '    <a class="next" onclick="plusSlides(1,this)">&#10095;</a>',
              '    <br>',
              '    <div style="text-align:center">']
    lines += ['    <span class="dot" onclick="currentSlide(' + str(i) + ',this)"></span>' for i in range(1, n + 1)]
    if close:
        lines += ['    </div>', '  </div>']
    return NL.join(lines)


if __name__ == "__main__":
    main()
