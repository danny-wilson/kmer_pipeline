# The kmer_pipeline image: the 2022-10-26 image, which provides R, genoPlotR, GEMMA, dsk,
# bowtie2, BLAST, MUMmer, samtools, Nextflow 22.04.5 and the compiled C++ tools (unchanged
# since then), with the pipeline scripts replaced from this commit and Biopython and pytest
# added. The base image is amd64 (x86-64) only. Build from a clean checkout of a release tag:
#   docker buildx build --platform linux/amd64 -t dannywilson/kmer_pipeline:2026-10-04 .
# The Dockerfile that built the base image is in this repository at tag 2022-10-26.
FROM dannywilson/kmer_pipeline:2022-10-26@sha256:d38900db59b92128dc7fb1118d71482452b37361253fc7348af72bf425bfb7ad
LABEL app="kmer_pipeline"
LABEL description="Pipeline for kmer (oligo)-based genome-wide association studies"
LABEL maintainer="Daniel Wilson"
LABEL version="2026-10-04"

# Set user and working directory
USER root
WORKDIR /tmp

# Install the Python packages missing from the base image, without upgrading any it has
RUN pip install --no-cache-dir --no-deps \
	biopython==1.83 \
	pytest==7.4.4 \
	pluggy==1.5.0 \
	iniconfig==2.0.0 \
	tomli==2.0.1 \
	exceptiongroup==1.2.2 \
	&& pip check \
	&& fix-permissions "${CONDA_DIR}"

# Replace the pipeline installed in the base image
COPY . /usr/share/kmer_pipeline.new
# (the base image's R workflow scripts are removed: only the two R files below remain)
RUN cd /usr/local/bin \
	&& rm -f *.R *.Rscript kmer_pipeline.nf report.js report.css \
	&& cd /usr/share/kmer_pipeline.new \
	&& install plot_figures.R Rscript_launcher.R python/*.py kmer_pipeline.nf report.js report.css /usr/local/bin \
	&& rm -r plot_figures.R Rscript_launcher.R python kmer_pipeline.nf report.js report.css \
	&& rm -r /usr/share/kmer_pipeline \
	&& mv /usr/share/kmer_pipeline.new /usr/share/kmer_pipeline \
	&& chmod -R a+rX,go-w /usr/share/kmer_pipeline

# Ignore any Python packages in the user's home directory
ENV PYTHONNOUSERSITE=1

# Set user, home and working directory
USER jovyan
ENV HOME=/home/jovyan
WORKDIR /home/jovyan
