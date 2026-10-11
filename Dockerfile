# The kmer_pipeline toolchain image: code-free. It holds only the slow-changing environment (apt
# packages, conda/samtools, R packages, dsk, gemma, Nextflow, and the compiled C++ tools) — no pipeline
# scripts, no kmer_pipeline.nf, no nextflow.config. Pipeline code is delivered through git (bin/,
# kmer_pipeline.nf, nextflow.config, pinned by `process.container`'s digest and a release tag), so a code
# change never needs an image rebuild or a slow Apptainer re-conversion.
#
# Closely follows the toolchain recipe at git tag 2022-10-26 (`git show 2022-10-26:Dockerfile`), with
# three deliberate changes: every download is checksum-verified; the C++ tools are built from this
# commit's own source in the build context, not downloaded as a tarball of a previously-pushed tag; and
# the build happens at the same path the 2022-10-26 recipe used, so the compiled binaries' debug info is
# directly comparable, not just functionally equivalent.
#
# Build (GitHub Actions only — Docker is not available on the maintainer's HPC):
#   docker buildx build --platform linux/amd64 \
#     --build-arg KMER_PIPELINE_VERSION=img-YYYY-MM-DD \
#     -t dannywilson/kmer_pipeline:img-YYYY-MM-DD .
FROM jupyter/datascience-notebook:2022-05-31@sha256:227d569a52d91cacea8f10c858f14b1099475b0435b1edba7372cbf898a8b6fb

ARG KMER_PIPELINE_VERSION=dev
LABEL app="kmer_pipeline"
LABEL description="Pipeline for kmer (oligo)-based genome-wide association studies: toolchain, code-free"
LABEL maintainer="Daniel Wilson"
LABEL version="${KMER_PIPELINE_VERSION}"

USER root
WORKDIR /tmp

# openjdk-11 and gcc-9/g++-9 at the exact versions the 2026-10-06 image has are no longer served by the
# live Ubuntu archive (superseded by later security updates). Fetched instead from a snapshot.ubuntu.com
# mirror by exact version, checksum-verified, and installed with dpkg before apt resolves anything else
# — so the meta-packages below (default-jdk, g++) bind to these exact versions, not whatever is currently
# the live-archive candidate.
RUN set -eux; \
	for spec in \
		"o/openjdk-lts openjdk-11-jdk_11.0.16+8-0ubuntu1~20.04_amd64.deb 1e03e1e2c0136efe480b7d41aa5041fae054986faa1a15a78d455ba4a391ab11" \
		"o/openjdk-lts openjdk-11-jdk-headless_11.0.16+8-0ubuntu1~20.04_amd64.deb 2a4d8d3c1b5636fae870e98b511cd3f3f7ada511f7b503e60bd99c2c424392df" \
		"o/openjdk-lts openjdk-11-jre_11.0.16+8-0ubuntu1~20.04_amd64.deb 5e7bc3c6ea5d974bfa72040a2bb0cf8a790befcecadf421d8dd2118766aaa129" \
		"o/openjdk-lts openjdk-11-jre-headless_11.0.16+8-0ubuntu1~20.04_amd64.deb dd7653014d564b62dc7c0a585d2b9aab557162b9fece6fbe5b4e37e5049b72f6" \
		"g/gcc-9 gcc-9-base_9.4.0-1ubuntu1~20.04.1_amd64.deb 0c9fae85529c13c9162dc81c5cbf9c4159f39df2898a9d9f94c404ce648e97a9" \
		"g/gcc-9 gcc-9_9.4.0-1ubuntu1~20.04.1_amd64.deb f2c456cd79618795fbb2ebb5748e3f9b6aada8e7795462d217b2abafe3eeaa90" \
		"g/gcc-9 g++-9_9.4.0-1ubuntu1~20.04.1_amd64.deb 819e3aad1e0dfeef9a27317a20a960f9df516d1da2cc874f1a0474e2af91dd07" \
	; do \
		pool=$(echo "$spec" | cut -d' ' -f1); \
		fname=$(echo "$spec" | cut -d' ' -f2); \
		sha=$(echo "$spec" | cut -d' ' -f3); \
		wget -q "https://snapshot.ubuntu.com/ubuntu/20221026T000000Z/pool/main/${pool}/${fname}"; \
		echo "${sha}  ${fname}" | sha256sum -c -; \
	done; \
	apt-get update --fix-missing; \
	dpkg -i gcc-9-base_*.deb gcc-9_*.deb g++-9_*.deb \
		openjdk-11-jre-headless_*.deb openjdk-11-jre_*.deb \
		openjdk-11-jdk-headless_*.deb openjdk-11-jdk_*.deb \
		|| true; \
	apt-get install -yf; \
	rm -f /tmp/*.deb

# The rest of the toolchain, pinned exactly as in the 2022-10-26 recipe. default-jdk and g++ are thin
# meta-packages that resolve against the exact openjdk-11/gcc-9/g++-9 versions already installed above,
# not a newer one the live archive might otherwise offer.
RUN apt-get update --fix-missing && apt-get -yqq install -f \
	bowtie2=2.3.5.1-6build1 \
	default-jdk=2:1.11-72 \
	g++=4:9.3.0-1ubuntu2 \
	libatlas-base-dev=3.10.3-8ubuntu7 \
	libgsl-dev=2.5+dfsg-6build1 \
	make=4.2.1-1.2 \
	mummer=3.23+dfsg-4build1 \
	ncbi-blast+=2.9.0-2 \
	zlib1g-dev=1:1.2.11.dfsg-2ubuntu1.5 \
	&& rm -rf /var/lib/apt/lists/*

# samtools via conda, pinned by an explicit spec captured from the existing, already-published
# 2026-10-06 image (not a fresh solve against today's conda channels, which could silently change
# numpy/pandas or other base-environment packages the pipeline's Python scripts rely on).
COPY docker/env/conda-explicit.txt /tmp/conda-explicit.txt
RUN conda install --file /tmp/conda-explicit.txt && rm /tmp/conda-explicit.txt

# dsk, from a checksummed precompiled binary release.
RUN wget -q https://github.com/GATB/dsk/releases/download/v2.3.3/dsk-v2.3.3-bin-Linux.tar.gz \
	&& echo "626a4663f70323a80834c013f1d9dcb4a3b3a73a1c9ab5da95a013783d81b0f6  dsk-v2.3.3-bin-Linux.tar.gz" | sha256sum -c - \
	&& tar -xzf dsk-v2.3.3-bin-Linux.tar.gz \
	&& install dsk-v2.3.3-bin-Linux/bin/* /usr/local/bin/ \
	&& rm -rf dsk-v2.3.3-bin-Linux.tar.gz dsk-v2.3.3-bin-Linux

# gemma0.93b, built from a checksummed source release.
RUN wget -q -O v0.1.tar.gz https://github.com/danny-wilson/gemma0.93b/archive/refs/tags/v0.1.tar.gz \
	&& echo "7368e51dd36a8988b7174024bb1bb5635860f7b414675a25fb254c7a2a54230b  v0.1.tar.gz" | sha256sum -c - \
	&& tar -xzf v0.1.tar.gz \
	&& cd gemma0.93b-0.1 \
	&& mkdir bin \
	&& make \
	&& install bin/gemma /usr/local/bin \
	&& cd /tmp && rm -rf v0.1.tar.gz gemma0.93b-0.1

# R packages, version-pinned via remotes (2.4.2, matching the 2026-10-06 image exactly). genoPlotR/ade4
# are kept for Tier-1 equivalence even though the pipeline's own R scripts no longer call them directly
# (OQ-S3-7: removing them is deferred to a later sub-project).
RUN wget -q https://cran.r-project.org/src/contrib/Archive/remotes/remotes_2.4.2.tar.gz \
	&& echo "f2ef875f24a485bf4f55a8c830f87cdd5db868f9a8cdb624dc452d0bf66ba516  remotes_2.4.2.tar.gz" | sha256sum -c - \
	&& R CMD INSTALL remotes_2.4.2.tar.gz \
	&& rm -f remotes_2.4.2.tar.gz \
	&& R -e "if (as.character(packageVersion('remotes')) != '2.4.2') stop('unexpected remotes version: ', packageVersion('remotes'))" \
	&& R -e "remotes::install_version('sp', version = '1.5-0', upgrade = FALSE, repos='https://cloud.r-project.org')" \
	&& R -e "remotes::install_version('pixmap', version = '0.4-12', upgrade = FALSE, repos='https://cloud.r-project.org')" \
	&& R -e "remotes::install_version('ade4', version = '1.7-19', upgrade = FALSE, repos='https://cloud.r-project.org')" \
	&& R -e "remotes::install_version('genoPlotR', version = '0.8.11', upgrade = FALSE, repos='https://cloud.r-project.org')"

# The compiled C++ tools: built from this commit's own source in the build context (never a downloaded
# tarball of a pushed git tag), at the same path the 2022-10-26 recipe used so the binaries' embedded
# debug-info build path matches, not just their behaviour. The compiled tools stay in the image, built by
# this toolchain stage — they are not interpreted scripts, so they do not move to bin/ (PI decision).
RUN mkdir -p /tmp/kmer_pipeline
COPY C++ /tmp/kmer_pipeline/C++
COPY Makefile /tmp/kmer_pipeline/Makefile
RUN cd /tmp/kmer_pipeline \
	&& make \
	&& install kmerlist2pattern pattern2kinship patterncounts patternmerge sort_strings stringlist2count stringlist2pattern /usr/local/bin \
	&& cd / && rm -rf /tmp/kmer_pipeline

# The Python packages missing from the base image, without upgrading any it has (unchanged from the
# 2026-10-06 recipe). pytest is environment, used by tests/, not pipeline code.
RUN pip install --no-cache-dir --no-deps \
	biopython==1.83 \
	pytest==7.4.4 \
	pluggy==1.5.0 \
	iniconfig==2.0.0 \
	tomli==2.0.1 \
	exceptiongroup==1.2.2 \
	&& pip check \
	&& fix-permissions "${CONDA_DIR}"

# Ignore any Python packages in the user's home directory
ENV PYTHONNOUSERSITE=1

# Set user, home and working directory
USER jovyan
ENV HOME=/home/jovyan
WORKDIR /home/jovyan
