# Start time
start.time = proc.time()[3]
args = commandArgs(trailingOnly=TRUE)

# Process arguments
help = paste(
"proteinkmermerge.Rscript merge protein kmer files names 1.txt.gz 2.txt.gz etc",
"Daniel Wilson (2018) kmermerge.Rscript modified to proteinkmermerge.Rscript by Sarah Earle (2019)",
"",
"Usage: qsub -t 1:p proteinkmermerge.Rscript n p output_prefix input_dir output_dir id_file kmer_length software_file process",
sep="\n")
if(length(args)!=9) { # SGE changed to 9 as added output_prefix, input_dir, output_dir, id_file, kmer_length, software_file, process
	cat(help,sep="\n")
	cat("Received arguments: ",args,"\n")
	stop("\nIncorrect usage\n")
}



###################################################################################################
## Functions and software paths
###################################################################################################

# sort_strings = "/well/bag/earle/scripts/software/sort_strings"

sort_final_file = function(output_dir = NULL, output_prefix = NULL, sort_strings = NULL, kmer_length = NULL){
	
	# Read in kmers
	kmerFile = paste0(output_dir, output_prefix,"_protein", kmer_length , ".kmermerge.unsorted.txt.gz")
	kmers = system(paste0("zcat ", kmerFile), intern = T)
	kmers = cbind(kmers, rep(1, length(kmers)))
	kmer_output_file = paste0(output_dir, output_prefix, "_protein", kmer_length, ".kmermerge.wdummycount.txt")
	write.table(kmers, file = kmer_output_file, row = F, col = F, sep = "\t", quote = F)
	# Gzip file
	system(paste0("gzip ", kmer_output_file))
	# Run sort strings
	kmer_output_file_gz = paste0(kmer_output_file, ".gz")
	final.kmer.txt.gz = paste0(output_dir, output_prefix, "_protein", kmer_length, ".kmermerge.sorted.wdummycount.txt.gz")
	sortCommand = paste(c(sort_strings, kmer_output_file_gz, "| gzip -c >", final.kmer.txt.gz), collapse=" ")
	system(sortCommand)
	# Remove column of dummy counts
	final.kmer.txt.gz.sorted = paste0(output_dir, output_prefix, "_protein", kmer_length, ".kmermerge.txt.gz")
	system(paste0("zcat ", final.kmer.txt.gz, " | cut -f1 | gzip -c > ", final.kmer.txt.gz.sorted))
	# Remove temp dummy count files
	system(paste0("rm ", kmer_output_file_gz))
	system(paste0("rm ", final.kmer.txt.gz))
	# Tests
	kmers_new = system(paste0("zcat ", final.kmer.txt.gz.sorted), intern = T)
	# cat("Length of input kmer file:", nrow(kmers), "\n")
	# cat("Length of output kmer file:", length(kmers_new), "\n")
	# cat("All kmers match:", all(!is.na(match(kmers[,1], kmers_new))),"\n")
	# cat("All kmers are in the same order:", all(kmers[,1]==kmers_new),"\n")
	if(nrow(kmers)!=length(kmers_new)) stop("Error: issue when sorting the final kmer file, number of kmers does not match","\n")
	if(!all(!is.na(match(kmers[,1], kmers_new)))) stop("Error: issue when sorting the final kmer file, the kmers do not match","\n")
	if(!all(kmers[,1]==kmers_new)) stop("Error: issue when sorting the final kmer file, the kmers are not in the same order","\n")
	# Remove the unsorted file
	system(paste0("rm ", output_dir, output_prefix,"_protein", kmer_length , ".kmermerge.unsorted.txt.gz"))
	
	
}

create_final_file = function(outfile = NULL, output_dir = NULL, output_prefix = NULL, kmer_length = NULL){
	
	# Check outfile isn't empty # SGE added
	outfile_size = system(paste0("ls -l ",outfile," | cut -d ' ' -f5"), intern = T) # SGE adde
	if(outfile_size==0) stop(outfile, " file is empty","\n") # SGE adde
	cmd = paste0("mv ",outfile," ", output_dir, output_prefix,"_protein", kmer_length, ".kmermerge.unsorted.txt") # SGE added output_dir output_prefix and _ before kmermerge.txt
	system2("/bin/bash",paste0("-c '",cmd,"'"), wait = TRUE) # SGE added wait = TRUE
	cmd=paste0("gzip ", output_dir, output_prefix,"_protein", kmer_length, ".kmermerge.unsorted.txt") # Added in as next script takes in a gzipped file
	system2("/bin/bash",paste0("-c '",cmd,"'"), wait = TRUE) # Added in as next script takes in a gzipped file # SGE added wait = TRUE

	# Remove all completed files
	# completed_files = dir(pattern=glob2rx(paste0(output_dir, output_prefix,".protein",kmer_length,"*.completed.txt")))
	completed_files = system(paste0("ls ", paste0(output_dir, output_prefix,".protein",kmer_length,"*.completed.txt")), intern = T)
	cmd = paste("rm",paste(completed_files,collapse=" "))
	system2("/bin/bash",paste0("-c '",cmd,"'"), wait = TRUE) # SGE added wait = TRUE

}

###################################################################################################




# Initialize variables
n = as.integer(args[1])
p = as.integer(args[2])
output_prefix = as.character(args[3]) # SGE added
input_dir = as.character(args[4]) # SGE added
output_dir = as.character(args[5]) # SGE added
id_file = as.character(args[6]) # SGE added
kmer_length = as.integer(args[7]) # SGE added
software_file = as.character(args[8]) # SGE added
process = as.integer(args[9]) # SGE added


if(is.na(n)) stop("Error: n must be an integer","\n") # SGE added
if(is.na(p)) stop("Error: p must be an integer","\n") # SGE added
if(!file.exists(input_dir)) stop("Error: input directory doesn't exist","\n") # SGE added
if(unlist(strsplit(input_dir,""))[length(unlist(strsplit(input_dir,"")))]!="/") input_dir = paste0(input_dir, "/") # SGE added
if(!file.exists(output_dir)) stop("Error: output directory doesn't exist","\n") # SGE added
if(unlist(strsplit(output_dir,""))[length(unlist(strsplit(output_dir,"")))]!="/") output_dir = paste0(output_dir, "/") # SGE added
if(!file.exists(id_file)) stop("Error: sample ID file doesn't exist","\n") # SGE added
if(!file.exists(software_file)) stop("Error: software file doesn't exist","\n") # SGE added
if(is.na(process)) stop("Error: process must be an integer","\n") # SGE added


# Read in software file
software_paths = read.table(software_file, h = T, sep = "\t", quote = "")
# Required software and script paths
# Begin with software
required_software = c("scriptpath")
if(any(is.na(match(required_software, as.character(software_paths$name))))) stop(paste0("Error: missing required software path in the software file - requires ",paste(required_software, collapse = ", ")),"\n")

script_location = as.character(software_paths$path)[which(tolower(as.character(software_paths$name))=="scriptpath")]
if(!dir.exists(script_location)) stop("Error: script location directory specified in the software paths file doesn't exist","\n")
sort_strings = file.path(script_location, "sort_strings")
if(!file.exists(sort_strings)) stop("Error: sort_strings path doesn't exist - check pipeline script location in the software file","\n")


# Read in sample IDs
id_file = read.table(id_file, h = T, sep = "\t")
sample_id = as.character(id_file$id)

if(n!=length(sample_id)) stop("n does not equal the number of samples in the sample ID file","\n")

b = as.integer(ceiling(n/p))
if(b==1) stop("Cannot have batchsize = 1. Try p < n/2")
p = as.integer(ceiling(n/b))
# t = as.integer(Sys.getenv("SGE_TASK_ID"))
t = process # SGE added
imax = as.integer(ceiling(log(n)/log(b)))

# Report variables
# cat("SGE job name: ",Sys.getenv("JOB_NAME"),"\n")
# cat("SGE job ID: ",Sys.getenv("JOB_ID"),"\n")
# cat("SGE task ID: ",Sys.getenv("SGE_TASK_ID"),"\n")
# cat("Arguments [string]: ",args,"\n")
# cat("Parameters [n batchsize processes taskid imax]: ",n,b,p,t,imax,"\n")

if(t>p) {
cat("Task",t,"not required\n")
quit("no")
}

# Input files are gzipped
gz = TRUE

# Merge
i = 0
while(TRUE) {
  i=i+1
  if(!((t %% b^(i-1))==0 | (t==p & i<=imax))) break()
  cat("t =",t,"i =",i,"\n")
  if(i==1) {
    # First round: merge source files
    outfile = paste0(output_dir, output_prefix,".protein",kmer_length,".j.",i,".",t,".txt") # SGE added output_dir output_prefix and . before j
    cat("Creating temp file:",outfile,"\n")
	outfile_completed = paste0(output_dir, output_prefix,".protein",kmer_length,".j.",i,".",t,".completed.txt") # SGE added output_dir output_prefix , protein and kmer_length and . before j
    beg = b*(t-1)+1
    end = min(b*t, n)
	cat("Beg:",beg," End:",end,"\n")
	if(end<beg) stop("Problem with input arguments, please check")
	if(gz){
		infiles = paste0(input_dir,sample_id[beg:end],".kmer", kmer_length, ".txt.gz") # SGE changed from beg:end to sample_id[beg:end] and removed .gz from file name and added kmer_length and kmers before .txt
	} else {
	    infiles = paste0(input_dir,sample_id[beg:end],".kmer", kmer_length, ".txt") # SGE changed from beg:end to sample_id[beg:end] and removed .gz from file name and added kmer_length and kmers before .txt
	}
	infiles_size = sapply(infiles, function(x) as.numeric(system(paste0("ls -l ",x," | cut -d ' ' -f5"), intern = T)), USE.NAMES = F) # SGE added - could make check greater than zero
	# nattempts = 0
	# while(!all(file.exists(infiles)) & !all(infiles_size>0)) { # SGE added check for non zero file size
		# nattempts=nattempts+1
		# if(nattempts>1000) stop("Could not find files",infiles)
		# Sys.sleep(1)
	# }
	if(!all(file.exists(infiles)) | !all(infiles_size>0)) stop("Could not find files or files empty ",paste(infiles, collapse = " ")) # SGE added replacing above as input files should not be being created as the script is running - could mean partially written files are used

    if(length(infiles)==1) {
    		if(gz){
	    		cmd = paste0("zcat ",infiles[1]," | cut -f1 > ",outfile) # SGE changed from zcat to cat
	    } else {
	    		cmd = paste0("cat ",infiles[1]," | cut -f1 > ",outfile) # SGE changed from zcat to cat
	    }
    } else if(length(infiles)==2) {
    		if(gz){
	    		cmd = paste0("sort -um <(zcat ",infiles[1]," | cut -f1) <(zcat ",infiles[2]," | cut -f1) > ",outfile) # SGE changed from zcat to cat
	    } else {
	    		cmd = paste0("sort -um <(cat ",infiles[1]," | cut -f1) <(cat ",infiles[2]," | cut -f1) > ",outfile) # SGE changed from zcat to cat
	    }
    } else {
    		if(gz){
      		cmd = paste0("sort -um <(zcat ",infiles[1]," | cut -f1) <(zcat ",infiles[2]," | cut -f1)",
        		paste0(" | sort -um - <(zcat ",infiles[3:length(infiles)]," | cut -f1)",collapse=""),
        		" > ",outfile) # SGE changed from zcat to cat
        	} else {
        		cmd = paste0("sort -um <(cat ",infiles[1]," | cut -f1) <(cat ",infiles[2]," | cut -f1)",
        	paste0(" | sort -um - <(cat ",infiles[3:length(infiles)]," | cut -f1)",collapse=""),
        	" > ",outfile) # SGE changed from zcat to cat
        	}
    }
	system2("/bin/bash",paste0("-c '",cmd,"'"), wait = TRUE) # SGE added wait = TRUE
	if(length(infiles)>1){
		kmers.i = scan(outfile, what = character(0), sep = "\n", quiet = TRUE) # SGE added as bash sort ignores *s so is keeping duplicate kmers
		kmers.i = unique(kmers.i) # SGE added
		cat(kmers.i, file = outfile, sep = "\n") # SGE added
	}
	system2("/bin/bash",paste0("-c 'touch ",outfile_completed,"'"), wait = TRUE) # SGE added wait = TRUE
  } else {
    # Subsequent rounds: merge merged files
    te = as.integer(ceiling(t/(b^(i-1)))*b^(i-1))
    outfile = paste0(output_dir, output_prefix,".protein",kmer_length,".j.",i,".",te,".txt") # SGE added output_dir output_prefix and . before j
    cat("Creating temp file:",outfile,"\n")
	outfile_completed = paste0(output_dir, output_prefix,".protein",kmer_length,".j.",i,".",te,".completed.txt") # SGE added output_dir output_prefix and . before j
    beg = te-b^(i-1)+b^(i-2)
    end = min(te,as.integer(ceiling(t/b^(i-2))*b^(i-2)))
	cat(paste0("Beg",i,":"),beg,paste0("End",i,":"),end,"\n") # SGE changed from 2 to i
	if(end<beg) stop("Problem with input arguments, please check")
    inc = b^(i-2)
    infiles = paste0(output_dir, output_prefix,".protein",kmer_length,".j.",i-1,".",seq(from=beg,to=end,by=inc),".txt") # SGE added output_dir output_prefix and . before j
	infiles_completed = paste0(output_dir, output_prefix,".protein",kmer_length,".j.",i-1,".",seq(from=beg,to=end,by=inc),".completed.txt") # SGE added
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
      cmd = paste0("sort -um ",infiles[1]," ",infiles[2]," > ",outfile)
    } else {
      cmd = paste0("sort -um ",infiles[1]," ",infiles[2],
        paste0(" | sort -um - ",infiles[3:length(infiles)],collapse=""),
        " > ",outfile)
    }
	system2("/bin/bash",paste0("-c '",cmd,"'"), wait = TRUE) # SGE added wait = TRUE
	if(length(infiles)>1){
		kmers.i = scan(outfile, what = character(0), sep = "\n", quiet = TRUE) # SGE added as bash sort ignores *s so is keeping duplicate kmers
		kmers.i = unique(kmers.i) # SGE added
		cat(kmers.i, file = outfile, sep = "\n") # SGE added
	}
	# If length of infiles is one, it has been moved and doesn't exist. If more than one, delete the temp files. # SGE added
	if(length(infiles)>1){
		cmd = paste("rm",paste(infiles,collapse=" "))
		# cat("Running remove command:",cmd,"\n")
		# cat("Infiles:",paste(infiles, collapse = " "), "\n")
		# cat("Infiles exist:",paste(file.exists(infiles), collapse = " "), "\n")
		system2("/bin/bash",paste0("-c '",cmd,"'"), wait = TRUE) # SGE added wait = TRUE
	}
	system2("/bin/bash",paste0("-c 'touch ",outfile_completed,"'"), wait = TRUE) # SGE added wait = TRUE
  }
}
if(t==p) {
	
	## Copy final file to end location # SGE changed to function
	create_final_file(outfile = outfile, output_dir = output_dir, output_prefix = output_prefix, kmer_length = kmer_length)
	
	## Sort final kmer file # SGE added function
	sort_final_file(output_dir = output_dir, output_prefix = output_prefix, sort_strings = sort_strings, kmer_length = kmer_length)
	
}

cat("Finished in",(proc.time()[3]-start.time)/60,"minutes\n")






