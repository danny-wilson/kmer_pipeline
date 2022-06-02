FROM jupyter/datascience-notebook:2022-05-31
LABEL app="kmer_pipeline"
LABEL description="Pipeline for kmer (oligo)-based genome-wide association studies"
LABEL maintainer="Daniel Wilson"
LABEL version="2022-06-01"

# Set user and working directory
USER root
WORKDIR /tmp

# Install packages
RUN apt-get update && apt-get -yqq install ncbi-blast+ mummer bowtie2 \
	libgsl0-dev libatlas-base-dev make g++ zlib1g-dev

# Install samtools with conda
RUN conda install -c bioconda -y samtools=1.15.1

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

# Install R package genoPlotR
RUN Rscript -e 'install.packages("genoPlotR", repos="https://www.stats.bris.ac.uk/R")'

# Install kmer_pipeline
RUN wget -O kmer_pipeline.tgz https://github.com/danny-wilson/kmer_pipeline/archive/refs/tags/2022-06-02.tar.gz \
	&& mkdir kmer_pipeline \
	&& tar -xvzf kmer_pipeline.tgz --directory kmer_pipeline --strip-components 1 \
	&& cd kmer_pipeline \
	&& make \
	&& mkdir /usr/bin/kmer_pipeline \
	&& install *.R *.Rscript kmerlist2pattern pattern2kinship patterncounts patternmerge sort_strings stringlist2count stringlist2pattern /usr/local/bin \
	&& rm *.R *.Rscript kmerlist2pattern pattern2kinship patterncounts patternmerge sort_strings stringlist2count stringlist2pattern \
	&& cd .. \
	&& mv kmer_pipeline /usr/share/

# Set user, home and working directory
USER jovyan
ENV HOME /home/jovyan
WORKDIR /home/jovyan
