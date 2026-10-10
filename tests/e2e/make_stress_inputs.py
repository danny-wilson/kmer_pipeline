#!/usr/bin/env python3
"""Build inputs for the stress variants from the example.

Writes to $KMER_E2E_ROOT/inputs/stress/:
  id_binary_na.txt   phenotype 1 if > 0 else 0, with three samples NA
  id_with_Ns.txt     sample 702 replaced by an assembly containing runs of N
  702_with_Ns.velvet_assembly.contigs.fa.gz
  id_truefalse.txt   phenotype TRUE if > 0 else FALSE (logical column)
  id_nan_pheno.txt   example phenotypes with sample row 7 set to NaN
  covariates_na.txt  covariates_example.txt with row 9's second column NA
ID files use @BASE@ for the base directory (make_reference.sh --id-file).

The example id_file/covariates themselves must already exist under
$KMER_E2E_ROOT/inputs/example (extracted from the image's own example/, as
make_reference.sh does) before this script is run.
"""
import gzip
import os

from e2e_config import load

ROOT = load()["KMER_E2E_ROOT"]
EXAMPLE = os.path.join(ROOT, "inputs", "example")
OUT = os.path.join(ROOT, "inputs", "stress")
NA_ROWS = {5, 12, 20}  # 1-based sample rows set to NA


def example_rows():
    with open(os.path.join(EXAMPLE, "id_file.txt")) as fh:
        header, *rows = [l.rstrip("\n").split("\t") for l in fh if l.strip()]
    return header, rows


def path(sample_file):
    return f"@BASE@/tb20/{os.path.basename(sample_file)}"


def main():
    os.makedirs(OUT, exist_ok=True)
    header, rows = example_rows()

    with open(os.path.join(OUT, "id_binary_na.txt"), "w") as fh:
        fh.write("\t".join(header) + "\n")
        for i, (sid, p, pheno) in enumerate(rows, start=1):
            value = "NA" if i in NA_ROWS else ("1" if float(pheno) > 0 else "0")
            fh.write(f"{sid}\t{path(p)}\t{value}\n")

    with gzip.open(os.path.join(EXAMPLE, "702.velvet_assembly.contigs.fa.gz"), "rt") as fh:
        lines = fh.read().split("\n")
    out, contig = [], 0
    for line in lines:
        if line.startswith(">"):
            contig += 1
        elif contig <= 3 and len(line) >= 150:
            line = line[:100] + "N" * 50 + line[150:]          # a run of N in each of the first 3 contigs
        elif contig == 4 and len(line) > 60:
            line = line[:60] + "N" + line[61:]                 # a single N
        out.append(line)
    n_file = "702_with_Ns.velvet_assembly.contigs.fa.gz"
    with gzip.open(os.path.join(OUT, n_file), "wt") as fh:
        fh.write("\n".join(out))

    with open(os.path.join(OUT, "id_with_Ns.txt"), "w") as fh:
        fh.write("\t".join(header) + "\n")
        for sid, p, pheno in rows:
            fh.write(f"{sid}\t{path(n_file) if sid == '702' else path(p)}\t{pheno}\n")

    with open(os.path.join(OUT, "id_truefalse.txt"), "w") as fh:
        fh.write("\t".join(header) + "\n")
        for sid, p, pheno in rows:
            fh.write(f"{sid}\t{path(p)}\t{'TRUE' if float(pheno) > 0 else 'FALSE'}\n")

    with open(os.path.join(OUT, "id_nan_pheno.txt"), "w") as fh:
        fh.write("\t".join(header) + "\n")
        for i, (sid, p, pheno) in enumerate(rows, start=1):
            fh.write(f"{sid}\t{path(p)}\t{'NaN' if i == 7 else pheno}\n")

    with open(os.path.join(ROOT, "inputs", "covariates_example.txt")) as fh:
        cov = [l.rstrip("\n").split("\t") for l in fh if l.strip()]
    cov[8][1] = "NA"
    with open(os.path.join(OUT, "covariates_na.txt"), "w") as fh:
        fh.write("".join("\t".join(r) + "\n" for r in cov))


if __name__ == "__main__":
    main()
