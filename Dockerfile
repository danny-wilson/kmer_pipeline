FROM jupyter/datascience-notebook:2022-05-31
LABEL app="kmer_pipeline"
LABEL description="Pipeline for kmer (oligo)-based genome-wide association studies"
LABEL maintainer="Daniel Wilson"
LABEL version="2022-10-26"

# Set user and working directory
USER root
WORKDIR /tmp

# Install packages
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

# Install samtools with conda
RUN conda install -c bioconda -y \
	samtools=1.15.1

# Install dsk from precompiled binary
RUN wget https://github.com/GATB/dsk/releases/download/v2.3.3/dsk-v2.3.3-bin-Linux.tar.gz \
	&& tar -xvzf dsk-v2.3.3-bin-Linux.tar.gz \
	&& install dsk-v2.3.3-bin-Linux/bin/* /usr/local/bin/

# Install gemma0.93b from source
RUN wget https://github.com/danny-wilson/gemma0.93b/archive/refs/tags/v0.1.tar.gz \
	&& tar -xvzf v0.1.tar.gz \
	&& cd gemma0.93b-0.1 \
	&& mkdir bin \
	&& make \
	&& install bin/gemma /usr/local/bin

# Install R package genoPlotR and dependency (version control using remotes)
RUN R -e "install.packages('remotes', repos = c(CRAN = 'https://cloud.r-project.org'))" \
	&& R -e "remotes::install_version('ade4', version = '1.7-19', upgrade = FALSE, repos='https://cloud.r-project.org')" \
	&& R -e "remotes::install_version('genoPlotR', version = '0.8.11', upgrade = FALSE, repos='https://cloud.r-project.org')"

# Install kmer_pipeline
RUN wget -O kmer_pipeline.tgz https://github.com/danny-wilson/kmer_pipeline/archive/refs/tags/2022-10-26.tar.gz \
	&& mkdir kmer_pipeline \
	&& tar -xvzf kmer_pipeline.tgz --directory kmer_pipeline --strip-components 1 \
	&& cd kmer_pipeline \
	&& make \
	&& mkdir /usr/bin/kmer_pipeline \
	&& ls \
	&& install *.R *.Rscript *.nf report.js report.css kmerlist2pattern pattern2kinship patterncounts patternmerge sort_strings stringlist2count stringlist2pattern /usr/local/bin \
	&& rm *.R *.Rscript *.nf report.js report.css kmerlist2pattern pattern2kinship patterncounts patternmerge sort_strings stringlist2count stringlist2pattern \
	&& cd .. \
	&& mv kmer_pipeline /usr/share/

# Install Nextflow
RUN wget -O nextflow https://github.com/nextflow-io/nextflow/releases/download/v22.04.5/nextflow-22.04.5-all \
	&& chmod a+x nextflow \
	&& mv nextflow /usr/local/bin/

# Set user, home and working directory
USER jovyan
ENV HOME /home/jovyan
WORKDIR /home/jovyan
