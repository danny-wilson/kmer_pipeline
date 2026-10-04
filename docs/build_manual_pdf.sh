#!/usr/bin/env bash
# Build kmer_pipeline_manual.pdf from docs/manual.md, with the pandoc (2.18) and XeLaTeX in the
# kmer_pipeline image. Run from the repository root, inside the container, e.g.
#   singularity exec --containall --cleanenv -B "$PWD" --pwd "$PWD" kmer_pipeline_2026-10-04.sif \
#     docs/build_manual_pdf.sh
# Output: kmer_pipeline_manual.pdf in the current directory.
set -euo pipefail
pandoc docs/manual.md -f gfm -o kmer_pipeline_manual.pdf \
	--lua-filter=docs/manual_pdf.lua --resource-path=docs \
	--include-in-header=docs/manual_pdf.tex --highlight-style=monochrome \
	--pdf-engine=xelatex -V fontfamily=fontspec \
	-V mainfont="Nimbus Sans" -V sansfont="Nimbus Sans" -V monofont="DejaVu Sans Mono" \
	-V fontsize=10pt -V papersize=a4 -V geometry:margin=2cm -V colorlinks=true
echo "Written kmer_pipeline_manual.pdf"
