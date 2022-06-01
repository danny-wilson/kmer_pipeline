
help = c("createcommands.R",
			"Usage: Rscript createcommands.R scriptDir analysisDir [output_prefix id_file analysis_file ref_fa ref_gb software_file] [nucmerident=90 min_count=1 minor_allele_threshold=0.01 bowtie_parameters=--very-sensitive samtools_filter=10 blastident=70]")

args = commandArgs(trailingOnly=TRUE)

if(length(args!=0)){
	if(args[1]=="-help" | args[1]=="-h"){
		cat(help,sep="\n")
		q("no")
	}
}

if(length(args)!=2 & length(args)!=8 & length(args)!=14){
	cat(help,sep="\n")
	cat("Received arguments: ", args, "\n")
	stop("\nIncorrect usage\n")
}

scriptDir = as.character(args[1])
analysisDir = as.character(args[2])

if(!dir.exists(scriptDir)) stop("Error: required argument scriptDir directory doesn't exist","\n")
if(unlist(strsplit(scriptDir,""))[length(unlist(strsplit(scriptDir,"")))]!="/") scriptDir = paste0(scriptDir,"/")
if(!dir.exists(analysisDir)) dir.create(analysisDir)
if(!file.exists(analysisDir)) stop("Error: analysis directory doesn't exist","\n")
analysisDir = normalizePath(analysisDir)
if(unlist(strsplit(analysisDir,""))[length(unlist(strsplit(analysisDir,"")))]!="/") analysisDir = paste0(analysisDir,"/")

if(length(args)>2){
	cat("Using command line input files","\n")
	prefix = as.character(args[3])
	id_file = as.character(args[4])
	analysis_file = as.character(args[5])
	ref_fa = as.character(args[6])
	ref_gb = as.character(args[7])
	software_file = as.character(args[8])
} else {
	cat("Using default example files","\n")
	prefix="saur12_example"
	id_file = paste0(scriptDir, "example/saur12_example_id_path_pheno.txt")
	analysis_file = paste0(scriptDir, "example/analysistype_nucleotide31_protein11.txt")
	ref_fa = paste0(scriptDir, "example/R00000022.fa")
	ref_gb = paste0(scriptDir, "example/R00000022.gbk")
	software_file = paste0(scriptDir, "example/pipeline_software_location.txt")
}
if(length(args)>8){
	cat("Using command line input parameters","\n")
	nucmerident = as.integer(args[9])
	min_count = as.integer(args[10])
	minor_allele_threshold = as.numeric(args[11])
	bowtie_parameters = as.character(args[12])
	samtools_filter = as.integer(args[13])
	blastident = as.integer(args[14])
} else {
	nucmerident = 90
	min_count = 1
	minor_allele_threshold = 0.01
	bowtie_parameters = "--very-sensitive"
	samtools_filter = 10
	blastident = 70
}


if(!file.exists(id_file)) stop("Error: id file doesn't exist","\n")
if(!file.exists(analysis_file)) stop("Error: analysis file doesn't exist","\n")
if(!file.exists(ref_fa)) stop("Error: reference fasta file doesn't exist","\n")
if(!file.exists(ref_gb)) stop("Error: reference genbank file doesn't exist","\n")
if(!file.exists(software_file)) stop("Error: software paths file doesn't exist","\n")
if(minor_allele_threshold>0.5 & minor_allele_threshold<1) stop("Error: minor allele threshold cannot be between 0.5-1","\n")
if(bowtie_parameters!="--very-sensitive"){
	if(!file.exists(bowtie_parameters)) stop("Error: bowtie_parameters file doesn't exist","\n")
}
if(is.na(samtools_filter)) stop("Error: samtools filter for bowtie2 results must be an integer","\n")


# If running the example data, create phenotype file with expanded file paths for the contigs
if(length(args)==2){
	id_file_example = read.table(id_file, h = T, sep = "\t")
	contig_paths = paste0(scriptDir, "example/", as.character(id_file_example[,2]))
	if(!all(file.exists(contig_paths))) stop("Error: not all example contig files exist", "\n")
	id_file_newpaths = id_file_example
	id_file_newpaths[,2] = contig_paths
	id_file = paste0(analysisDir, "saur12_example_id_path_pheno_expandedpaths.txt")
	write.table(id_file_newpaths, file = id_file, row = F, col = T, sep = "\t", quote = F)
	if(!file.exists(id_file)) stop("Error: created example id_file with expanded file paths does not exist","\n")
}


# Read in sample number
n = system(paste0("wc -l ",id_file), intern = T)
n = as.numeric(unlist(strsplit(n," "))[1])-1
cat("Total number of samples:",n,"\n")

# Read in kmer type and lengths to create commands for
analysis_type = read.table(analysis_file, h = T, sep = "\t")
kmer_type = as.character(analysis_type$kmertype)
kmer_length = as.numeric(analysis_type$kmerlength)

cat("Analyses:\n")
for(i in 1:length(kmer_type)) cat(kmer_type[i]," ",kmer_length[i],"\n")
cat("\n")

# Get directory for std out and err files
stdouterr = paste0(analysisDir, "runfiles/")

cat("Step 1 - count kmers","\n")
cat("########################################################################\n\n")

tc = n; if(tc>200) tc = 200

cat(paste0("qsub -t 1:", n, " -tc ", tc, " -o ", stdouterr, " -e ", stdouterr, " ", scriptDir, "countkmers.Rscript ", id_file, " ", analysisDir, " ", prefix, " ", software_file, " ", analysis_file),"\n\n")

for(i in 1:length(kmer_type)){
	
	cat("\n")
	cat("########################################################################\n")	
	cat("Steps for kmer type",kmer_type[i],"and kmer length",kmer_length[i],"\n")
	cat("########################################################################\n\n")
	
	extended_path_prefix = paste0(analysisDir, prefix, "_", kmer_type[i], kmer_length[i])
	
	cat("Step 2 - create unique kmer list","\n")
	cat("########################################################################\n\n")
	
	p = floor(n/5); p = ifelse(p==0 | p==1, 2, p)
	if(p>20) p = 20
	
	cat(paste0("qsub -t 1:", p, " -tc ", p, " -o ", stdouterr, " -e ", stdouterr, " ", scriptDir, "createfullkmerlist.Rscript ", n, " ", p, " ", prefix, " ", analysisDir, " ", id_file, " ", kmer_type[i], " ", kmer_length[i], " ", software_file),"\n\n")

	cat("Step 3 - create kmer presence/absence matrix and kinship matrix","\n")
	cat("########################################################################\n\n")
	
	p = 20
	
	cat(paste0("qsub -t 1:", p, " -tc ", p, " -o ", stdouterr, " -e ", stdouterr, " ", scriptDir, "stringlist2patternandkinship.Rscript ", p, " ", id_file, " ", extended_path_prefix, ".kmermerge.txt.gz ", extended_path_prefix, "_kmers_filepaths.txt ", analysisDir, " ", prefix, " ", kmer_type[i], " ", software_file, " ", kmer_length[i], " 1"),"\n\n")
	
	cat("Step 4 - run GEMMA","\n")
	cat("########################################################################\n\n")
	
	p = 10
	
	cat(paste0("qsub -t 1:",p," -tc ",p, " -o ", stdouterr, " -e ", stdouterr," ", scriptDir, "rungemma.Rscript ", p, " ", extended_path_prefix, " ", id_file, " ", prefix, " ", analysisDir, " ", kmer_type[i], " ", kmer_length[i], " ", software_file),"\n\n")
	
	cat("Step 5A - run contig alignment","\n")
	cat("########################################################################\n\n")
	
	tc = n; if(tc>200) tc = 200
	
	cat(paste0("qsub -t 1:",n," -tc ",tc, " -o ", stdouterr, " -e ", stdouterr," ", scriptDir, "kmercontigalign.Rscript ", n, " ", prefix, " ", analysisDir, " ", id_file, " ", ref_fa, " ", ref_gb, " ", kmer_type[i], " ", kmer_length[i], " ", nucmerident, " ", extended_path_prefix, ".kmermerge.txt.gz ", software_file),"\n\n")
	
	if(kmer_type[i]=="nucleotide"){
		cat("Step 5B - or run bowtie2 mapping","\n")
		cat("########################################################################\n\n")
	
		cat(paste0("qsub -o ", stdouterr, " -e ", stdouterr, " ", scriptDir, "runbowtie.Rscript ", prefix, " ", analysisDir, " ", extended_path_prefix, " ", ref_fa, " ", kmer_type[i], " ", kmer_length[i], " ", software_file, " ", bowtie_parameters, " ", samtools_filter),"\n\n")
	}
	
	cat("Step 6A - plot Manhattan figures using contig alignment positions","\n")
	cat("########################################################################\n\n")
	
	# Get the reference name
	ref.name = scan(ref_fa, what = character(0), sep = "\n", nlines = 1, quiet = TRUE)
	ref.name = unlist(strsplit(ref.name, " "))[1]
	if(substr(ref.name,1,1)!=">") stop("Error: reference fasta file does not begin with a name starting with '>'","\n")
	ref.name = substr(ref.name,2,1e6)

	cat(paste0("qsub -o ", stdouterr, " -e ", stdouterr, " ", scriptDir, "plotManhattan.Rscript ", prefix, " ", analysisDir, " ", extended_path_prefix, " ", ref_gb, " ", ref_fa, " ", analysisDir, kmer_type[i], "kmer", kmer_length[i], "_kmergenealign/", prefix, "_", kmer_type[i], kmer_length[i], "_", ref.name, "_gene_id_name_lookup.txt ", id_file, " ", nucmerident, " ", min_count, " ", kmer_type[i], " ", kmer_length[i], " ", minor_allele_threshold, " ", software_file, " ", blastident),"\n\n")
	
	if(kmer_type[i]=="nucleotide"){
		cat("Step 6B - plot Manhattan figures using bowtie2 mapping positions","\n")
		cat("########################################################################\n\n")
	
		cat(paste0("qsub -o ", stdouterr, " -e ", stdouterr, " ", scriptDir, "plotManhattanbowtie.Rscript ", prefix, " ", analysisDir, " ", extended_path_prefix, " ", ref_gb, " ", ref_fa, " ", id_file, " ", kmer_type[i], " ", kmer_length[i], " ", minor_allele_threshold, " ", samtools_filter, " ", software_file, " ", blastident),"\n\n")
	}
	
}







