# kmer_pipeline
Pipeline for kmer (oligo)-based genome-wide association studies

Implementing methods described in

**Genome-wide association studies of global Mycobacterium tuberculosis resistance to thirteen antimicrobials in 10,228 genomes**
The CRyPTIC Consortium (2021)
*bioRxiv* doi: 10.1101/2021.09.14.460272 ([preprint](https://www.biorxiv.org/content/10.1101/2021.09.14.460272))

**Identifying lineage effects when controlling for population structure improves power in bacterial association studies.**
Earle, S. G., Wu, C.-H., Charlesworth, J., Stoesser, N., Gordon, N. C., Walker, T. M., Spencer, C. C. A., Iqbal, Z., Clifton, D. A., Hopkins, K. L., Woodford, N., Smith, E. G., Ismail, N., Llewelyn, M. J., Peto, T. E., Crook, D. W., McVean, G., Walker, A. S. and D. J. Wilson (2016)
*Nature Microbiology* 1: 16041 ([preprint](http://arxiv.org/abs/1510.06863))

## Authors
Written by [Sarah G Earle](https://github.com/sgearle) and [Daniel J Wilson](https://github.com/danny-wilson) at the [University of Oxford](https://www.ox.ac.uk), [Big Data Institute](https://www.bdi.ox.ac.uk).

## Funding
This project was supported by the [Wellcome Trust](https://wellcome.org), the [Royal Society](https://royalsociety.org) and the [Robertson Foundation](https://robertsonfoundation.org).

## Dependencies
The pipeline utilizes software including [GEMMA](https://github.com/genetics-statistics/GEMMA), [DSK](https://github.com/GATB/dsk), [NCBI BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi), [MUMmer](http://mummer.sourceforge.net), [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml), [Samtools](http://www.htslib.org), the [GNU Scientific Library](https://www.gnu.org/software/gsl/), the [Automatically Tuned Linear Algebra Software](http://math-atlas.sourceforge.net)

Original code and scripts were written in [R](https://www.r-project.org) and [C++](https://isocpp.org), using the [genoPlotR library](https://cran.r-project.org/web/packages/genoPlotR/index.html) and the [zstr library](https://github.com/mateidavid/zstr), a C++ wrapper for the [zlib library](https://github.com/madler/zlib). C++ code was compiled with the [GNU Compiler Collection](https://gcc.gnu.org).

The container was written for [Docker](https://www.docker.com) based on the [Jupyter Data Science Notebook](https://jupyter-docker-stacks.readthedocs.io/en/latest/index.html).

The example data comprise genomes of *Staphylococcus aureus* including the [MRSA252 reference genome](https://pubmed.ncbi.nlm.nih.gov/15213324/) (NCBI accession number [BX571856.1](https://www.ncbi.nlm.nih.gov/nuccore/49240382)) and twelve genomes with fusidic acid sensitivity originally reported by [Gordon and colleagues](https://pubmed.ncbi.nlm.nih.gov/24501024/) (NCBI accession number [PRJNA308283](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA308283)).

## License
The zstr headers are licensed under the MIT license. The myutils headers are licensed under the GNU Lesser General Public License Version 3. All other code is licensed under the GNU General Public License v3.0.

## Installation
To build from source
  docker build github.com/danny-wilson/kmer_pipeline

To download a prebuilt docker image
  docker pull danny-wilson/kmer_pipeline:latest

To build a [Singularity](https://sylabs.io/guides/3.3/user-guide/index.html) container
  singularity build kmer.sif docker://dannywilson/kmer_pipeline:latest

## Running
To launch the Jupyter Data Science Notebook using Docker
  docker container run --name kmer -p 8888:8888 -v $USERMOUNT:/home/jovyan/ext/ dannywilson/kmer_pipeline:latest
where $USERMOUNT represents the local directory outside the container you wish to access inside the container from /home/jovyan/ext. Having launched the container, navigate to the user interface in a web browser by following one of the URLs provided at the command line.

To launch a shell within the Docker container
  docker exec -it kmer bash

To launch a shell within the Singularity container
  singularity exec --containall --cleanenv --home $USERMOUNT:/home/jovyan --workdir /tmp kmer.sif /bin/bash

Detailed information on the example data to follow.

