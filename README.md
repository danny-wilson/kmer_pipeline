# kmer_pipeline
Pipeline for kmer (oligo)-based genome-wide association studies

Implementing methods described in

**Genome-wide association studies of global *Mycobacterium tuberculosis* resistance to thirteen antimicrobials in 10,228 genomes**
The CRyPTIC Consortium (2021)
*PLOS Biology* 20: e3001755 ([article](https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.3001755))

**Identifying lineage effects when controlling for population structure improves power in bacterial association studies.**
Earle, S. G., Wu, C.-H., Charlesworth, J., Stoesser, N., Gordon, N. C., Walker, T. M., Spencer, C. C. A., Iqbal, Z., Clifton, D. A., Hopkins, K. L., Woodford, N., Smith, E. G., Ismail, N., Llewelyn, M. J., Peto, T. E., Crook, D. W., McVean, G., Walker, A. S. and D. J. Wilson (2016)
*Nature Microbiology* 1: 16041 ([preprint](http://arxiv.org/abs/1510.06863))

## Authors
Written by [Sarah G Earle](https://github.com/sgearle) and [Daniel J Wilson](https://github.com/danny-wilson) at the [University of Oxford](https://www.ox.ac.uk), [Big Data Institute](https://www.bdi.ox.ac.uk).

## Funding
This project was supported by the [Wellcome Trust](https://wellcome.org), the [Royal Society](https://royalsociety.org) and the [Robertson Foundation](https://robertsonfoundation.org).

## Dependencies
The pipeline utilizes software including [GEMMA](https://github.com/genetics-statistics/GEMMA), [DSK](https://github.com/GATB/dsk), [NCBI BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi), [MUMmer](http://mummer.sourceforge.net), [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml), [Samtools](http://www.htslib.org), the [GNU Scientific Library](https://www.gnu.org/software/gsl/), the [Automatically Tuned Linear Algebra Software](http://math-atlas.sourceforge.net)

The container implements a [Nextflow](https://www.nextflow.io) pipeline which runs on top of [Apache Groovy](https://groovy-lang.org) and [Java](https://www.java.com).

Original code and scripts were also written in [R](https://www.r-project.org) and [C++](https://isocpp.org), using the [genoPlotR library](https://cran.r-project.org/web/packages/genoPlotR/index.html) and the [zstr library](https://github.com/mateidavid/zstr), a C++ wrapper for the [zlib library](https://github.com/madler/zlib). C++ code was compiled with the [GNU Compiler Collection](https://gcc.gnu.org).

The container was written using [Docker](https://www.docker.com) and based on the [Jupyter Data Science Notebook](https://jupyter-docker-stacks.readthedocs.io/en/latest/index.html).

## License
The zstr headers are licensed under the MIT license. The myutils headers are licensed under the GNU Lesser General Public License Version 3. All other code is licensed under the GNU General Public License v3.0.

## Installation
To download a prebuilt [Docker](https://www.docker.com) image

    docker pull dannywilson/kmer_pipeline:2022-10-26

To build a [Singularity](https://sylabs.io/guides/3.3/user-guide/index.html) container

    singularity pull -F docker://dannywilson/kmer_pipeline:2022-10-26

## Running the Nextflow pipeline
For instructions on running the Nextflow pipeline, including the *Mycobacterium tuberculosis* example, download [the manual](https://github.com/danny-wilson/docs/blob/main/kmer_pipeline_nf.pdf).

## Running the Jupyter Data Science Notebook
To launch as a Jupyter Data Science Notebook using Docker

    docker container run --name kmer -p 8888:8888 -v $MNT_DIR:/home/jovyan dannywilson/kmer_pipeline:latest
  
where $MNT_DIR represents the local directory outside the container you wish to access inside the container from /home/jovyan. Having launched the container, navigate to the user interface in a web browser by following one of the URLs provided at the command line.
