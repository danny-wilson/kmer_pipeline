
# Start time
start.time = proc.time()[3]
args = commandArgs(trailingOnly=TRUE)

# Process arguments - SGE added output_prefix kmerpath_file output_dir to usage
help = paste(
"kmermerge.Rscript merge kmer files",
"Daniel Wilson (2018)",
"",
"Usage: nucleotidekmermerge.R n p output_prefix input_dir output_dir id_file kmer_length process",
sep="\n")
if(length(args)!=8) { # SGE changed to 8 as added output_prefix, input_dir, output_dir, id_file, kmer_length, process
	cat(help,sep="\n")
	cat("Received arguments: ",args,"\n")
	stop("\nIncorrect usage\n")
}


###################################################################################################
## Functions and software paths
###################################################################################################



create_final_file = function(outfile = NULL, output_dir = NULL, output_prefix = NULL, kmer_length = NULL){
	
	# Check outfile isn't empty # SGE added
	outfile_size = system(paste0("ls -l ",outfile," | cut -d ' ' -f5"), intern = T) # SGE added
	if(outfile_size==0) stop(outfile, " file is empty","\n") # SGE added
	cmd = paste0("mv ",outfile," ",output_dir, output_prefix,"_nucleotide", kmer_length, ".kmermerge.txt") # SGE added output_prefix and _ before kmermerge.txt
	stopifnot(system2("/bin/bash",paste0("-c '",cmd,"'"), wait = TRUE)==0)
	cmd=paste0("gzip ",output_dir, output_prefix,"_nucleotide", kmer_length, ".kmermerge.txt") # Added in as next script takes in a gzipped file
	stopifnot(system2("/bin/bash",paste0("-c '",cmd,"'"), wait = TRUE)==0)

	# Remove all completed files
	completed_files = dir(output_dir,pattern=glob2rx(paste0( output_prefix,".nucleotide",kmer_length,"*.completed.txt")),full.names=TRUE)
	cmd = paste("rm",paste(completed_files,collapse=" "))
	stopifnot(system2("/bin/bash",paste0("-c '",cmd,"'"), wait = TRUE)==0)

}

###################################################################################################

# Initialize variables
n = as.integer(args[1])
p = as.integer(args[2])
output_prefix = as.character(args[3]) # SGE added
input_dir = as.character(args[4]) # SGE added
output_dir = as.character(args[5]) # SGE added
id_file = as.character(args[6]) # SGE added
kmer_length = as.numeric(args[7]) # SGE added
process = as.integer(args[8]) # SGE added


if(is.na(n)) stop("Error: n must be an integer","\n") # SGE added
if(is.na(p)) stop("Error: p must be an integer","\n") # SGE added
if(!file.exists(input_dir)) stop("Error: input directory doesn't exist","\n") # SGE added
if(unlist(strsplit(input_dir,""))[length(unlist(strsplit(input_dir,"")))]!="/") input_dir = paste0(input_dir, "/") # SGE added
if(!file.exists(output_dir)) stop("Error: output directory doesn't exist","\n") # SGE added
if(unlist(strsplit(output_dir,""))[length(unlist(strsplit(output_dir,"")))]!="/") output_dir = paste0(output_dir, "/") # SGE added
if(!file.exists(id_file)) stop("Error: sample ID file doesn't exist","\n") # SGE added
if(is.na(process)) stop("Error: process must be an integer","\n") # SGE added


# Read in sample IDs
id_file = read.table(id_file, h = T, sep = "\t")
sample_id = as.character(id_file$id)
kmerpaths = paste0(input_dir,sample_id,".kmer", kmer_length, ".txt.gz")

b = as.integer(ceiling(n/p))
if(b==1) stop("Cannot have batchsize = 1. Try p < n/2")
p = as.integer(ceiling(n/b))
# t = as.integer(Sys.getenv("SGE_TASK_ID"))
t = process # SGE added
imax = as.integer(ceiling(log(n)/log(b)))

# Report variables
# cat("SGE job name: ",Sys.getenv("JOB_NAME"),"\n")
# cat("SGE job ID: ",Sys.getenv("JOB_ID"),"\n")
# cat("SGE task ID: ",Sys.getenv("SGE_TASK_ID"),"\n")
# cat("Running on host: ",Sys.getenv("HOST"),"\n")
# cat("Arguments [string]: ",args,"\n")
# cat("Parameters [n batchsize processes taskid imax]: ",n,b,p,t,imax,"\n")
if(t>p) {
cat("Task",t,"not required\n")
quit("no")
}

# Merge
i = 0
while(TRUE) {
  i=i+1
  if(!((t %% b^(i-1))==0 | (t==p & i<=imax))) break()
  cat("t=",t,"i=",i,"\n")
  if(i==1) {
    # First round: merge source files
    outfile = paste0(output_dir, output_prefix,".nucleotide",kmer_length,".j.",i,".",t,".txt") # SGE added output_dir output_prefix nucleotide and kmer_length and . before j
    cat("Creating temp file:",outfile,"\n")
	outfile_completed = paste0(output_dir, output_prefix,".nucleotide",kmer_length,".j.",i,".",t,".completed.txt") # SGE added output_dir output_prefix nucleotide and kmer_length and . before j
    beg = b*(t-1)+1
    end = min(b*t, n)
	cat("Beg:",beg,"End:",end,"\n")
	if(end<beg) stop("Problem with input arguments, please check")
    infiles = kmerpaths[beg:end] # SGE changed from paste0(input_dir,[beg:end],".txt.gz")
	infiles_size = sapply(infiles, function(x) as.numeric(system(paste0("zcat ",x," | wc -c"), intern = T)), USE.NAMES = F)
	# nattempts = 0
	# while(!all(file.exists(infiles)) & !all(infiles_size>0)) { # SGE added check for non zero file size
		# nattempts=nattempts+1
		# if(nattempts>1000) stop("Could not find files",infiles)
		# Sys.sleep(1)
	# }
	if(!all(file.exists(infiles)) | !all(infiles_size>0)) stop("Could not find files or files empty",paste(infiles, collapse = " ")) # SGE added replacing above as input files should not be being created as the script is running - could mean partially written files are used

	# Unzip the input files and remove the counts
	onestep = FALSE
	if(!onestep) {
		tmpinfiles = gsub(".txt.gz",paste0(".sorted.j.",i,".",t,".txt"),infiles)
		subcmds = paste0("zcat ",infiles," | cut -d \" \" -f1 > ",tmpinfiles)
		for(subcmd in subcmds) stopifnot(system2("/bin/bash",paste0("-c '",subcmd,"'"), wait = TRUE)==0)
	}
	
    if(length(infiles)==1) {
		if(onestep) cmd = paste0("zcat ",infiles[1]," | cut -d \" \" -f1 > ",outfile) # Changed cut delimiter
		if(!onestep) cmd = paste0("cat ",tmpinfiles[1]," > ",outfile) # Changed cut delimiter
    } else if(length(infiles)==2) {
		if(onestep) cmd = paste0("LC_ALL=C sort -um <(zcat ",infiles[1]," | cut -d \" \" -f1) <(zcat ",infiles[2]," | cut -d \" \" -f1) > ",outfile) # Changed cut delimiter
		if(!onestep) cmd = paste0("LC_ALL=C sort -um ",tmpinfiles[1]," ",tmpinfiles[2]," > ",outfile) # Changed cut delimiter
    } else {
		if(onestep) cmd = paste0("LC_ALL=C sort -um <(zcat ",infiles[1]," | cut -d \" \" -f1) <(zcat ",infiles[2]," | cut -d \" \" -f1)",
			paste0(" | LC_ALL=C sort -um - <(zcat ",infiles[3:length(infiles)]," | cut -d \" \" -f1)",collapse=""),
			" > ",outfile) # Changed cut delimiter
		if(!onestep) cmd = paste0("LC_ALL=C sort -um ",tmpinfiles[1]," ",tmpinfiles[2],
			paste0(" | LC_ALL=C sort -um - ",tmpinfiles[3:length(tmpinfiles)],collapse=""),
			" > ",outfile) # Changed cut delimiter
    }
	stopifnot(system2("/bin/bash",paste0("-c '",cmd,"'"), wait = TRUE)==0)
	stopifnot(system2("/bin/bash",paste0("-c 'touch ",outfile_completed,"'"), wait = TRUE)==0)
	
	# Remove the temporary input files
	if(!onestep) unlink(tmpinfiles)
  } else {
    # Subsequent rounds: merge merged files
    te = as.integer(ceiling(t/(b^(i-1)))*b^(i-1))
    outfile = paste0(output_dir, output_prefix,".nucleotide",kmer_length,".j.",i,".",te,".txt") # SGE added output_dir output_prefix nucleotide and kmer_length and . before j
    cat("Creating temp file:",outfile,"\n")
	outfile_completed = paste0(output_dir, output_prefix,".nucleotide",kmer_length,".j.",i,".",te,".completed.txt") # SGE added output_dir output_prefix nucleotide and kmer_length and . before j
    beg = te-b^(i-1)+b^(i-2)
    end = min(te,as.integer(ceiling(t/b^(i-2))*b^(i-2)))
	cat(paste0("Beg",i,":"),beg,paste0("End",i,":"),end,"\n") # SGE changed from 2 to i
	if(end<beg) stop("Problem with input arguments, please check")
    inc = b^(i-2)
    infiles = paste0(output_dir, output_prefix,".nucleotide",kmer_length,".j.",i-1,".",seq(from=beg,to=end,by=inc),".txt") # SGE added output_dir output_prefix and . before j
	infiles_completed = paste0(output_dir, output_prefix,".nucleotide",kmer_length,".j.",i-1,".",seq(from=beg,to=end,by=inc),".completed.txt") # SGE added
	nattempts = 0
	while(!all(file.exists(infiles_completed)) | !all(file.exists(infiles))) {
		nattempts=nattempts+1
		if(nattempts>100) stop("Could not find files",infiles_completed)
		Sys.sleep(60) # SGE changed from 1 to 60
	}

	infiles_size = sapply(infiles, function(x) as.numeric(system(paste0("ls -l ",x," | cut -d ' ' -f5"), intern = T)), USE.NAMES = F) # SGE added this and in below while statement - could set higher than zero, an appropriate number for completed
	if(any(infiles_size==0)) stop("One or more file size is zero ",infiles,"\n") # SGE added

    if(length(infiles)==1) {
      cmd = paste0("mv ",infiles[1]," ",outfile)
    } else if(length(infiles)==2) {
      cmd = paste0("LC_ALL=C sort -um ",infiles[1]," ",infiles[2]," > ",outfile)
    } else {
      cmd = paste0("LC_ALL=C sort -um ",infiles[1]," ",infiles[2],
        paste0(" | LC_ALL=C sort -um - ",infiles[3:length(infiles)],collapse=""),
        " > ",outfile)
    }
	stopifnot(system2("/bin/bash",paste0("-c '",cmd,"'"), wait = TRUE)==0)
	# If length of infiles is one, it has been moved and doesn't exist. If more than one, delete the temp files. # SGE added
	if(length(infiles)>1){
		cmd = paste("rm",paste(infiles,collapse=" "))
		stopifnot(system2("/bin/bash",paste0("-c '",cmd,"'"), wait = TRUE)==0)
	}
	stopifnot(system2("/bin/bash",paste0("-c 'touch ",outfile_completed,"'"), wait = TRUE)==0)
  }
}
if(t==p) {
	## Copy final file to end location # SGE changed to function
	create_final_file(outfile = outfile, output_dir = output_dir, output_prefix = output_prefix, kmer_length = kmer_length)
}

cat("Finished in",(proc.time()[3]-start.time)/60,"minutes\n")
