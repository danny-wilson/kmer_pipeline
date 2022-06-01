
read_reference = function(ref_file) {
	r = scan(ref_file,what=character(0),sep="\n")
	rcat = paste(r[2:length(r)],collapse="")
	return(unlist(strsplit(rcat,"")))
}

reorder_reference_gbk = function(ref_gb = NULL){
	
	# Read in reference genbank file
	ref = read_dna_seg_from_file(ref_gb)
	ref = ref[which(ref$feature=="CDS"),]
	# For each name, if there is more than one entry, label as 'gene_1', 'gene_2'
	ref_unique_names_multiple = table(as.character(ref$name))
	ref_unique_names_multiple = names(ref_unique_names_multiple[which(ref_unique_names_multiple>1)])
	for(i in 1:length(ref_unique_names_multiple)){
		w.i = which(ref$name==ref_unique_names_multiple[i])
		ref$name[w.i] = paste0(ref$name[w.i], "_", 1:length(w.i))
	}
	ref = ref[order(as.numeric(ref$start)),]
	return(ref)
	
}

rev_base = list("A" = "T", "T" = "A", "C" = "G", "G" = "C", "-" = "-")


transcribe = function(x) {
	y = t(sapply(1:nrow(x),function(i) totriplet(x[i,])))
	rownames(y) = rownames(x)
	return(y)
}

translate = function(x,oneLetter=FALSE) {
	x = toupper(x)
	tr = t(apply(x,1,function(y)sapply(y,function(i) {aa=geneticCode[[i]];ifelse(is.null(aa),"---",aa)} )))
	if(oneLetter) tr = t(apply(tr,1,function(y) oneLetterCodes[y]))
	rownames(tr) = rownames(x)
	return(tr)
}

revcompl = c("A"="T","C"="G","G"="C","T"="A")
revcompl_full = c("A"="T","C"="G","G"="C","T"="A","-"="-","N"="N")

rc = function(x) revcompl[rev(x)]
rc_full = function(x) revcompl_full[rev(x)]

totriplet = function(x) {
	L = floor(length(x)/3)*3
	paste(x[seq(1,L,by=3)],x[seq(2,L,by=3)],x[seq(3,L,by=3)],sep="")
}

geneticCode = list(
"TTT"="Phe","TTC"="Phe","TTA"="Leu","TTG"="Leu",
"TCT"="Ser","TCC"="Ser","TCA"="Ser","TCG"="Ser",
"TAT"="Tyr","TAC"="Tyr","TAA"="STO","TAG"="STO",
"TGT"="Cys","TGC"="Cys","TGA"="STO","TGG"="Trp",
"CTT"="Leu","CTC"="Leu","CTA"="Leu","CTG"="Leu",
"CCT"="Pro","CCC"="Pro","CCA"="Pro","CCG"="Pro",
"CAT"="His","CAC"="His","CAA"="Gln","CAG"="Gln",
"CGT"="Arg","CGC"="Arg","CGA"="Arg","CGG"="Arg",
"ATT"="Ile","ATC"="Ile","ATA"="Ile","ATG"="Met",
"ACT"="Thr","ACC"="Thr","ACA"="Thr","ACG"="Thr",
"AAT"="Asn","AAC"="Asn","AAA"="Lys","AAG"="Lys",
"AGT"="Ser","AGC"="Ser","AGA"="Arg","AGG"="Arg",
"GTT"="Val","GTC"="Val","GTA"="Val","GTG"="Val",
"GCT"="Ala","GCC"="Ala","GCA"="Ala","GCG"="Ala",
"GAT"="Asp","GAC"="Asp","GAA"="Glu","GAG"="Glu",
"GGT"="Gly","GGC"="Gly","GGA"="Gly","GGG"="Gly")

oneLetterCodes = list("Gly" = "G", "Ala" = "A", "Leu" = "L", "Met" = "M", "Phe" = "F", "Trp" = "W", "Lys" = "K", "Gln" = "Q", "Glu" = "E", "Ser" = "S", "Pro" = "P", "Val" = "V", "Ile" = "I", "Cys" = "C", "Tyr" = "Y", "His" = "H", "Arg" = "R", "Asn" = "N", "Asp" = "D", "Thr" = "T", "STO" = "*")

# Don't think this should be used now
# oneLetterCodesCorrected = list("Gly" = "G", "Ala" = "A", "Leu" = "L", "Met" = "M", "Phe" = "F", "Trp" = "W", "Lys" = "K", "Gln" = "Q", "Glu" = "E", "Ser" = "S", "Pro" = "P", "Val" = "V", "Ile" = "I", "Cys" = "C", "Tyr" = "Y", "His" = "H", "Arg" = "R", "Asn" = "N", "Asp" = "D", "Thr" = "T", "STO" = "*", "---" = "X")

### Translate functions


translate_function = function(contig = NULL, frame = NULL, oneLetterCodes = NULL, revcompl = NULL){
	contig = unlist(strsplit(contig,""))
	if(any(is.na(match(contig, c("A", "C", "G", "T"))))) stop("Error: contig contains non base characters", "\n")
	if(frame<=3){
		return(paste(unlist(oneLetterCodes[translate(matrix(totriplet(contig[frame:length(contig)]), nrow = 1))]), collapse = ""))
	} else {
		contig = rev(revcompl[contig])
		if(any(is.na(contig))) stop("Error: reverse complemented contig contains NAs", "\n")
		return(paste(unlist(oneLetterCodes[translate(matrix(totriplet(contig[(frame-3):length(contig)]), nrow = 1))]), collapse = ""))
	}
}

get_output_file = function(outDir, id){
	return(paste0(outDir, id, "_translated_all_reading_frames.fa"))
}

write_proteins_to_file = function(contig_names = NULL, frame = NULL, contigs = NULL, id = NULL, append = FALSE, outDir = NULL){
	out = paste(paste0(contig_names,"_rf",frame), contigs, sep = "\n")
	outfile = get_output_file(outDir, id)
	cat(out, file = outfile, sep = "\n", append = append)
}

remove_Ns = function(contig){
	contigs = unlist(strsplit(contig,"N"))
	return(contigs[which(contigs!="")])
}

read_contigs_removeNs = function(contig_file = NULL){
	
	if(substr(contig_file, (nchar(contig_file)-2), nchar(contig_file))==".gz"){
		contigs = scan(gzfile(contig_file), what = character(0), sep = "\n", quiet = TRUE)
	} else {
		contigs = scan(contig_file, what = character(0), sep = "\n", quiet = TRUE)
	}
	contig_name_positions = which(sapply(contigs, function(x) any(unlist(gregexpr(">",x))!=-1), USE.NAMES = F))
	cat("Number of contigs:", length(contig_name_positions),"\n")
	contig_names = sapply(contigs[contig_name_positions], function(x) unlist(strsplit(x," "))[1], USE.NAMES = F)
	contig_names_final = c()
	contig_seq = c()
	for(i in 1:length(contig_names)){
		
		if(i!=length(contig_names)){
			contig.i = paste(contigs[(contig_name_positions[i]+1):(contig_name_positions[i+1]-1)], collapse = "")
		} else {
			contig.i = paste(contigs[(contig_name_positions[i]+1):length(contigs)], collapse = "")
		}
		contig.i = toupper(contig.i)
		if(any(unlist(strsplit(contig.i,""))=="N")) contig.i = remove_Ns(contig.i)
		contig_seq = c(contig_seq, contig.i)
		contig_names_final = c(contig_names_final, paste0(contig_names[i], "_", 1:length(contig.i)))
		
	}
	return(list("contig_seq" = contig_seq, "contig_names_final" = contig_names_final))

}

translate_6_frames = function(contig_path = NULL, id = NULL, outDir = NULL, oneLetterCodes = NULL, revcompl = NULL){
	
	# if(substr(contig_path, nchar(contig_path)-2, nchar(contig_path))==".gz") gzipped = TRUE else gzipped = FALSE
	
	# if(gzipped){
		# contigs = scan(gzfile(contig_path), what = character(0), sep = "\n", quiet = TRUE)
	# } else {
		# contigs = scan(contig_path, what = character(0), sep = "\n", quiet = TRUE)
	# }
		
	# contig_name_positions = which(sapply(contigs, function(x) any(unlist(gregexpr(">",x))!=-1), USE.NAMES = F))
	# cat("Number of contigs:", length(contig_name_positions),"\n")
	# contig_names = sapply(contigs[contig_name_positions], function(x) unlist(strsplit(x," "))[1], USE.NAMES = F)
	# contig_names_final = c()
	# contig_seq = c()
	# for(i in 1:length(contig_names)){
		
		# if(i!=length(contig_names)){
			# contig.i = paste(contigs[(contig_name_positions[i]+1):(contig_name_positions[i+1]-1)], collapse = "")
		# } else {
			# contig.i = paste(contigs[(contig_name_positions[i]+1):length(contigs)], collapse = "")
		# }
		# contig.i = toupper(contig.i)
		# if(any(unlist(strsplit(contig.i,""))=="N")) contig.i = remove_Ns(contig.i)
		# contig_seq = c(contig_seq, contig.i)
		# contig_names_final = c(contig_names_final, paste0(contig_names[i], "_", 1:length(contig.i)))
		
	# }
	
	contig_seq = read_contigs_removeNs(contig_path)
	contig_names_final = contig_seq$contig_names_final
	contig_seq = contig_seq$contig_seq
	
	rf1 = sapply(contig_seq, translate_function, frame = 1, oneLetterCodes = oneLetterCodes, revcompl = revcompl, USE.NAMES = FALSE)
	rf2 = sapply(contig_seq, translate_function, frame = 2, oneLetterCodes = oneLetterCodes, revcompl = revcompl, USE.NAMES = FALSE)
	rf3 = sapply(contig_seq, translate_function, frame = 3, oneLetterCodes = oneLetterCodes, revcompl = revcompl, USE.NAMES = FALSE)
	rf4 = sapply(contig_seq, translate_function, frame = 4, oneLetterCodes = oneLetterCodes, revcompl = revcompl, USE.NAMES = FALSE)
	rf5 = sapply(contig_seq, translate_function, frame = 5, oneLetterCodes = oneLetterCodes, revcompl = revcompl, USE.NAMES = FALSE)
	rf6 = sapply(contig_seq, translate_function, frame = 6, oneLetterCodes = oneLetterCodes, revcompl = revcompl, USE.NAMES = FALSE)

	write_proteins_to_file(contig_names = contig_names_final, frame = 1, contigs = rf1, id = id, outDir = outDir)
	write_proteins_to_file(contig_names = contig_names_final, frame = 2, contigs = rf2, id = id, append = TRUE, outDir = outDir)
	write_proteins_to_file(contig_names = contig_names_final, frame = 3, contigs = rf3, id = id, append = TRUE, outDir = outDir)
	write_proteins_to_file(contig_names = contig_names_final, frame = 4, contigs = rf4, id = id, append = TRUE, outDir = outDir)
	write_proteins_to_file(contig_names = contig_names_final, frame = 5, contigs = rf5, id = id, append = TRUE, outDir = outDir)
	write_proteins_to_file(contig_names = contig_names_final, frame = 6, contigs = rf6, id = id, append = TRUE, outDir = outDir)
	
	# gzip final file
	final_file = get_output_file(outDir, id)
	system(paste0("gzip ", final_file))
	return(paste0(final_file,".gz"))
	
}



create_gene_lookup = function(ref = NULL, ref_length = NULL){

	gene_start = ref$start
	gene_end = ref$end
	gene_strand = ref$strand
	
	# Get a gene/IR ID for every position in the reference
	# A position can have multiple gene IDs as there are overlapping genes
	ref.pos.gene.id = vector(mode = "list", length = ref_length)
	for(i in 1:nrow(ref)){
		for(j in as.numeric(ref$start[i]):as.numeric(ref$end[i])){
			ref.pos.gene.id[[j]] = c(ref.pos.gene.id[[j]], i)
		}
	}

	# Intergenic ID only assigned if the start of one gene is after the end of the previous gene
	intergenic_names = c()
	for(i in 2:nrow(ref)){
		inter_start = as.numeric(ref$end[i-1])+1
		inter_end = as.numeric(ref$start[i])-1
		if(inter_start<=inter_end){
			intergenic_names = c(intergenic_names, paste(c(ref$name[i-1], ref$name[i]), collapse = ":"))
			gene_start = c(gene_start, inter_start)
			gene_end = c(gene_end, inter_end)
			for(j in (inter_start):(inter_end)){
				ref.pos.gene.id[[j]] = c(ref.pos.gene.id[[j]], nrow(ref)+length(intergenic_names))
			}
		}
		if(i==nrow(ref)){
			intergenic_names = c(intergenic_names, paste0(ref$name[length(ref$name)], ":"))
			gene_start = c(gene_start, (ref$end[nrow(ref)]+1))
			gene_end = c(gene_end, ref_length)
			for(j in (as.numeric(ref$end[nrow(ref)]):length(ref.pos.gene.id))){
				ref.pos.gene.id[[j]] = c(ref.pos.gene.id[[j]], nrow(ref)+length(intergenic_names))
			}
		}
	}
	# Create a lookup table to get a gene name/IR name for every ID assigned
	gene_lookup = as.list(c(ref$name, intergenic_names)); names(gene_lookup) = 1:length(c(ref$name, intergenic_names))
	gene_lookup = cbind(gene_lookup, names(gene_lookup))
	
	# Set the strand for all intergenic regions to 1 - there is no correct strand for these
	gene_strand = c(gene_strand, rep(1, length(intergenic_names)))
	gene_lookup = cbind(gene_lookup, gene_start, gene_end, gene_strand)
	
	return(list("gene_lookup" = gene_lookup, "ref.pos.gene.id" = ref.pos.gene.id))
	
}
