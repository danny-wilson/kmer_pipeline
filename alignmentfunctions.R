
######################################################################
## Plotting kmer/gene alignment functions
######################################################################


add_xaxis_top = function(reverse.xaxis.start = NULL, forward.xaxis.start = NULL, genestart = NULL, gene_db_i = NULL){

	# If the x-axis is to be plotted on the forward strand
	if(is.null(reverse.xaxis.start)){

		# If the x-axis is to be remapped to new positions change genestart
		if(!is.null(forward.xaxis.start)) genestart = forward.xaxis.start

		# Find the first position which is a multiple of 10 for clean numbering
		axis.start = which(c(genestart:(genestart+ncol(gene_db_i)-1))%%10==0)[1]

		# If none of the positions are a multiple of 10, label every position
		if(is.na(axis.start)){

			axis(3, at = 1:ncol(gene_db_i), labels = c(genestart:(genestart+ncol(gene_db_i)-1)),
				 cex.axis = 0.7, lwd.ticks = NA, line = 0.4, lwd = NA, xpd = T)
			axis(3, at = 1:ncol(gene_db_i), lwd.ticks = NA, label = NA, lwd = 0.7, line = 0.81, xpd = T)
			axis(3, at = 1:ncol(gene_db_i), labels = NA, cex.axis = 0.7, lwd = 0.7, line = 0.81, xpd = T)

		# Else label in multiples of 10
		} else {

			axis(3, at = seq(from = axis.start, by = 10, to = ncol(gene_db_i)), labels = c(seq(from = c(genestart+axis.start-1),
					 by = 10, to = c(genestart +ncol(gene_db_i)-1))), cex.axis = 0.7, lwd.ticks = NA, line = 0.4, lwd = NA, xpd = T)
			axis(3, at = c(1-0.5, ncol(gene_db_i)+0.5), lwd.ticks = NA, label = NA, lwd = 0.7, line = 0.81, xpd = T)
			axis(3, at = seq(from = axis.start, by = 10, to = ncol(gene_db_i)), labels = NA, cex.axis = 0.7, lwd = 0.7, line = 0.81, xpd = T)

		}


	# Else if it is to be plotted on the reverse strand
	} else if(!is.null(reverse.xaxis.start)){

		# Find the first position which is a multiple of 10 for clean numbering
		axis.start = which(c(reverse.xaxis.start:(reverse.xaxis.start-ncol(gene_db_i)+1))%%10==0)[1]

		# If none of the positions are a multiple of 10, label every position
		if(is.na(axis.start)){

			axis(3, at = 1:ncol(gene_db_i), labels = c(genestart:(genestart-ncol(gene_db_i)+1)),
				 cex.axis = 0.7, lwd.ticks = NA, line = 0.4, lwd = NA, xpd = T)
			axis(3, at = 1:ncol(gene_db_i), lwd.ticks = NA, label = NA, lwd = 0.7, line = 0.81, xpd = T)
			axis(3, at = 1:ncol(gene_db_i), labels = NA, cex.axis = 0.7, lwd = 0.7, line = 0.81, xpd = T)

		# Else label in multiples of 10
		} else {

			axis(3, at = seq(from = axis.start, by = 10, to = ncol(gene_db_i)), labels = c(seq(from = c(reverse.xaxis.start-axis.start+1),
					 by = -10, to = c(reverse.xaxis.start-ncol(gene_db_i)+1))), cex.axis = 0.7, lwd.ticks = NA, line = 0.4, lwd = NA, xpd = T)
			axis(3, at = c(1-0.5, ncol(gene_db_i)+0.5), lwd.ticks = NA, label = NA, lwd = 0.7, line = 0.81, xpd = T)
			axis(3, at = seq(from = axis.start, by = 10, to = ncol(gene_db_i)), labels = NA, cex.axis = 0.7, lwd = 0.7, line = 0.81, xpd = T)

		}


	}


}



get_kmer_cols_grad = function(odds){
	
	blue_black = colorRamp(c("#d3d3d3","#8181d7"))
	black_red = colorRamp(c("#d3d3d3","#d78181"))
	
	kmercols1 = rgb(blue_black(st(odds[which(odds<1 & odds!="Inf")])^0.6), maxColorValue = 256)
	kmercols2 = rgb(black_red(st(odds[which(odds>1 & odds!="Inf")])^0.6), maxColorValue = 256) 
	kmercols = rep("#000000",length(odds))
	kmercols[which(odds<1 & odds!="Inf")] = kmercols1
	kmercols[which(odds>1 & odds!="Inf")] = kmercols2
	return(kmercols)
}

get_kmer_cols_lm = function(odds, lm.col){
	
	blue_black = colorRamp(c("#767676","#8181d7"))
	black_red = colorRamp(c("#767676","#d78181"))
	
	lm = lm.col
	if(length(which(odds!="Inf"))!=0){
		if(length(which(odds<1 & odds!="Inf"))!=0){
			if(max(odds[which(odds<1 & odds!="Inf")])==0){
				kmercols1 = rep("blue",length(which(odds<1 & odds!="Inf")))
			} else {
				kmercols1 = rgb(lm*(1-st(odds[which(odds<1 & odds!="Inf")])^0.99), lm*(1-st(odds[which(odds<1 & odds!="Inf")])^0.99), lm+st(odds[which(odds<1 & odds!="Inf")])^0.99*(1-lm))
			}
		}
		if(length(which(odds>1 & odds!="Inf"))!=0){
			kmercols2 = rgb(lm+st(odds[which(odds>1 & odds!="Inf")])^0.99*(1-lm), lm*(1-st(odds[which(odds>1 & odds!="Inf")])^0.99), lm*(1-st(odds[which(odds>1 & odds!="Inf")])^0.99))
		}
		kmercols = rep("#d3d3d3",length(odds))
		if(length(which(odds<1 & odds!="Inf"))!=0){
			kmercols[which(odds<1 & odds!="Inf")] = kmercols1
		}
		if(length(which(odds>1 & odds!="Inf"))!=0){
			kmercols[which(odds>1 & odds!="Inf")] = kmercols2
		}
	}
	kmercols[which(odds=="Inf")] = "red"
	return(kmercols)
}

# Function to colour variants by the beta point estimate
get_betaCOL = function(beta, cols, se){
	beta.min = as.numeric(beta[2])-(as.numeric(beta[3])*se)
	beta.max = as.numeric(beta[2])+(as.numeric(beta[3])*se)
	beta.point = as.numeric(beta[2])
	
	if(beta.min<0 & beta.max>0){
		return("#d3d3d3")
	} else {
		if(beta.point<(-0.5)){
			return(cols[1])
		} else if(beta.point>=(-0.5) & beta.point<0){
			return(cols[2])
		} else if(beta.point>0 & beta.point<=0.5){
			return(cols[3])
		} else if(beta.point>0.5){
			return(cols[4])
		}
	}
}


# Function to colour variants by the beta point estimate
get_betaCOL_2cols = function(beta, cols, se){
	beta.min = as.numeric(beta[2])-(as.numeric(beta[3])*se)
	beta.max = as.numeric(beta[2])+(as.numeric(beta[3])*se)
	beta.point = as.numeric(beta[2])
	if(beta.min<0 & beta.max>0){
		return("#d3d3d3")
	} else {
		if(beta.point<(0)){
			return(cols[1])
		} else if(beta.point>0){
			return(cols[4])
		}
	}
}


get_alignment_col = function(seq = NULL, cols = NULL, len.col = NULL, snp_cols = NULL, which_translucent = NULL, col_lib = NULL){
	col = matrix("#d3d3d3",nrow=nrow(seq),ncol=ncol(seq))
	col[which(seq=="-")] = "#ffffff"
	if(!is.null(cols)){
		for(i in 1:(nrow(col)-len.col)){
			col[(i+as.numeric(len.col)),which(col[(i+as.numeric(len.col)),]!="#ffffff")] = cols[i]
		}
	}
	
	if(is.null(snp_cols)) snp_cols = matrix(rep("#d3d3d3",ncol(seq)),nrow=1)
	for(i in 1:ncol(col)){
		if(length(unique(seq[which(seq[,i]!="-"),i]))>1 | any(snp_cols[,i]!="#d3d3d3")){
			for(j in 1:nrow(col)){
				col[j,i] = as.character(col_lib[seq[j,i]])
			}
		}
	}
	col[which(col=="NULL")] = "#ffffff"
	
	if(!is.null(which_translucent)){
		if(length(which_translucent)>0){
			for(i in c(1:(nrow(col)-len.col))[which_translucent]){
				col[(i+as.numeric(len.col)),] = paste0(col[(i+as.numeric(len.col)),],"55")
			}
		}
	}
	
	return(col)
}

read_reference = function(ref_file) {
	r = scan(ref_file,what=character(0),sep="\n", quiet = TRUE)
	rcat = paste(r[2:length(r)],collapse="")
	return(unlist(strsplit(rcat,"")))
}


build_kmer_matrix = function(kmers = NULL, genestart = 0, ref = NULL, snps = NULL, kmerpos = NULL, snp_num = NULL){
	seq_all_reads = matrix("-", nrow = (length(kmers)+nrow(ref)+snp_num), ncol = ncol(ref))
	seq_all_reads[unique(1:nrow(ref)),] = as.character(ref)
	if(genestart!=0) genestart = -genestart+1
	for(i in 1:length(kmers)){
		pos = kmerpos[i] + genestart
		# Check now that some kmers have changed length due to indels that they still overlap the region
		if((pos+nchar(kmers[i])-1)>=1 & pos<=ncol(ref)){
		# If the starting position of the kmer, leftmost, is greater than 1 (start of figure) and ends before the end of the figure
		# if((pos+30) <= ncol(ref) & pos>=1){
		# Updating to take in kmers of length different to 31 - happens with indels
		if((pos+nchar(kmers[i])-1) <= ncol(ref) & pos>=1){
			# seq_all_reads[(i+nrow(ref)+ snp_num),pos:(pos+30)] = as.vector(unlist(strsplit(kmers[i],"")))
			seq_all_reads[(i+nrow(ref)+ snp_num),pos:(pos+nchar(kmers[i])-1)] = as.vector(unlist(strsplit(kmers[i],""))) # Updating to take in kmers of length different to 31
		# Else if it overlaps with the end of the figure but starts after the start of the figure
		} else if(pos>=1){
			seq_all_reads[(i+nrow(ref)+ snp_num),pos:ncol(seq_all_reads)] = as.vector(unlist(strsplit(kmers[i],"")))[1:length(pos:ncol(seq_all_reads))]
		# Else if it overlaps with the start of the figure
		} else {
			# Check that the lengths match up
			# if(length(1:(pos+30))!=length(as.vector(unlist(strsplit(kmers[i],"")))[-c(1:(31-length(c(1:(pos+30)))))])) stop("Error","\n")
			if(length(1:(pos+nchar(kmers[i])-1))!=length(as.vector(unlist(strsplit(kmers[i],"")))[-c(1:(nchar(kmers[i])-length(c(1:(pos+nchar(kmers[i])-1)))))])) stop("Error","\n")
			# seq_all_reads[(i+nrow(ref)+ snp_num),1:(pos+30)] = as.vector(unlist(strsplit(kmers[i],"")))[-c(1:(31-length(c(1:(pos+30)))))]
			seq_all_reads[(i+nrow(ref)+ snp_num),1:(pos+nchar(kmers[i])-1)] = as.vector(unlist(strsplit(kmers[i],"")))[-c(1:(nchar(kmers[i])-length(c(1:(pos+nchar(kmers[i])-1)))))]
		}
	}
	}
	# seq_all_reads = seq_all_reads[nrow(seq_all_reads):1,]
	return(seq_all_reads)
}




# col_lib = list("A" = "#008000", "C" = "blue","G" = "black", "T" = "#ff9000")
# col_lib = list("A" = "#009E73", "C" = "#0072B2","G" = "black", "T" = "#E69F00")
col_lib_nuc = list("A" = "#009E73", "C" = "#0072B2","G" = "black", "T" = "#E69F00")
col_lib_pro = list("*" = "#000000",
			   "A" = "#009E73",
			   "C" = "#0072B2",
			   "D" = "#E69F00",
			   "E" = "#A01FF0",
			   "F" = "#50FF00",
			   "G" = "#FAC0CB",
			   "H" = "#F8A503",
			   "I" = "#ADD8E6",
			   "K" = "#0C008B",
			   "L" = "#8B0000",
			   "M" = "#1A6400",
			   "N" = "#a52a2a",
			   "P" = "#ffbbff",
			   "Q" = "#F78C02",
			   "R" = "#df9797",
			   "S" = "#90EE90",
			   "T" = "#FDFF00",
			   "V" = "#9d9d00",
			   "W" = "#ff0f39",
			   "Y" = "#A52A29",
			   "-" = "#ffffff")

# rev_compl_col = list("#009E73" = "#E69F00", "#E69F00" = "#009E73", "#0072B2" = "black", "black" = "#0072B2")
# rev_compl_col = list("#008000" = "#ff9000", "#ff9000" = "#008000", "blue" = "black", "black" = "blue", "#d3d3d3" = "#d3d3d3")
rev_compl_col = list("#008000" = "#ff9000", "#ff9000" = "#008000", "blue" = "black", "black" = "blue", "#d3d3d3" = "#d3d3d3")
rev_base = list("A" = "T", "T" = "A", "C" = "G", "G" = "C", "-" = "-")
colour_selection = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00")
beta_cols_list = colour_selection[c(2, 3, 1, 6)]



plot_alignment = function(prefix = NULL, seq = NULL, align_col = NULL, ref_num = NULL, snp_num = NULL, labs.cex = 1, legend.txt = NULL, legend.fill = NULL, sep_lines = NULL, lm.col = 0.8, max.odds = NULL, genestart = NULL, lwd = NULL, bonferroni_threshold = NULL, n_snp_cols = NULL, plot.dim = c(100,40), legend.cex = 1, legend.xpos = NULL, legend.ypos = NULL, text.xpos = NULL, text.ypos = NULL, axis.cex = 1, legend.lty = NULL, legend.pch = NULL, legend.col = NULL, line.sep.lwd = NULL, legend.odds = NULL, ref.name = NULL, reverse.xaxis = NULL, reverse.xaxis.start = NULL, kmer.numbers = NULL, close.file = NULL, margins = NULL, plot_ref = NULL, forward.xaxis.start = NULL){
	png(paste0(prefix,"_alignment.png"), width = plot.dim[1], height = plot.dim[2], units = "cm", res = 600)
	par(mar=margins)
	image(x = c(1:ncol(seq)), y = c(1:nrow(seq)), z = t(matrix(c(1:(ncol(seq)*nrow(seq))), nrow = nrow(seq), ncol = ncol(seq))), col = align_col, bty = "n", xlab = "", ylab = "", xaxt = "n", yaxt = "n")
	# cat("Plotted image","\n")
	sapply(c(1:ncol(seq)), function(i) abline(v = (i - 0.5), col = "white", lwd = line.sep.lwd))
	#cat("Plotted first lines","\n")
	sapply(c((ref_num+snp_num+1):nrow(seq)), function(i) abline(h = (i - 0.5), col = "white", lwd = line.sep.lwd))
	#cat("Plotted second lines","\n")
	for(i in 1:length(sep_lines)){
		abline(h = sep_lines[i], col = "white",lwd= 1)
	}
	if(!is.null(bonferroni_threshold)){
		abline(h = bonferroni_threshold+ref_num+snp_num+0.5, col = "black", lwd = 0.5, lty = 2)
	}
	#if(!is.null(snp_num)) abline(h = ref_num+snp_num+0.5, col = "#808080",lwd=0.5)
	if(!is.null(kmer.numbers)){
		for(i in 1:nrow(kmer.numbers)){
			w = which(align_col[i+ref_num+snp_num,]!="white"); w = w[length(w)]+1
		
			text(x = w, y = (i+ref_num+snp_num), labels = paste(kmer.numbers[i,], collapse = ", "), cex = 0.1, adj = 0)
		
		}
	}
	
	if(is.null(reverse.xaxis.start)){
		if(!is.null(forward.xaxis.start)) genestart = forward.xaxis.start
		axis.start = which(c(genestart:(genestart+ncol(seq)-1))%%10==0)[1]
		if(is.na(axis.start)){
			axis(1, at = 1:ncol(seq), labels = c(genestart:(genestart+ncol(seq)-1)), cex.axis = axis.cex, lwd.ticks = NA, line = -0.4, lwd = NA)
		} else {
			axis(1, at = seq(from = axis.start, by = 10, to = ncol(seq)), labels = c(seq(from = c(genestart+axis.start-1), by = 10, to = c(genestart+ncol(seq)-1))), cex.axis = axis.cex, lwd.ticks = NA, line = -0.4, lwd = NA)
		}
	} else if(!is.null(reverse.xaxis.start)){
		axis.start = which(c(reverse.xaxis.start:(reverse.xaxis.start-ncol(seq)+1))%%10==0)[1]
		if(is.na(axis.start)){
			axis(1, at = 1:ncol(seq), labels = c(seq(from = c(reverse.xaxis.start), by = -1, to = c(reverse.xaxis.start-ncol(seq)+1))), cex.axis = axis.cex, lwd.ticks = NA, line = -0.4, lwd = NA)
		} else {
			axis(1, at = seq(from = axis.start, by = 10, to = ncol(seq)), labels = c(seq(from = c(reverse.xaxis.start-axis.start+1), by = -10, to = c(reverse.xaxis.start-ncol(seq)+1))), cex.axis = axis.cex, lwd.ticks = NA, line = -0.4, lwd = NA)
		}
	}
	
	axis(1, at = c(1-0.5, ncol(seq)+0.5), lwd.ticks = NA, label = NA, lwd = axis.cex, line = 0.01)
	if(is.na(axis.start)){
		axis(1, at = seq(from = 1, by = 1, to = ncol(seq)), labels = NA, cex.axis = axis.cex, lwd = axis.cex, line = 0.01)
	} else {
		axis(1, at = seq(from = axis.start, by = 10, to = ncol(seq)), labels = NA, cex.axis = axis.cex, lwd = axis.cex, line = 0.01)
	}
	
	
	#axis(1,at = seq(from=1,by=10,to=ncol(seq)), cex = 0.5)
	if(plot_ref) text(x = 1,y = 2.5,ref.name,xpd=TRUE, pos = 2, cex = labs.cex, adj = 1)
	if(n_snp_cols!=0){
		text(x = 1,y = 6,"SNPs",xpd=TRUE, pos = 2, cex = labs.cex, adj = 1)
		if(snp_num==5 | snp_num==6) text(x = 1, y = 5.8, "SNP TYPE", xpd = T, pos = 2, cex = labs.cex, adj = 1)
		if(snp_num==6) text(x = 1, y = 6.8, "FEATURES", xpd = T, pos = 2, cex = labs.cex, adj = 1)
	}
	
	
	if(!is.null(legend.txt)){
		par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
		plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
		if(!is.null(legend.lty)){
			if(any(!is.na(legend.lty))){
				legend("topright", legend = legend.txt, border = NA, bty = "n", cex = legend.cex, lty = legend.lty, pch = legend.pch, col = legend.col, text.col = "white", inset = c(-0.01,0), xpd = T)
			}
		}
		legend("topright", legend.txt, fill = legend.fill, border = NA, bty = "o", cex = legend.cex, col = legend.col, bg = NA)
		if(legend.odds){
			get_legend_col_alignment(lm.col = lm.col, max.odds = max.odds, legend.xpos = legend.xpos, legend.ypos = legend.ypos, text.xpos = text.xpos, text.ypos = text.ypos)
		}
	}
	
	if(close.file) dev.off()
}

run_plot_alignment = function(kmerseq = NULL, snp_cols = NULL, kmerpos = NULL, effect_size = NULL, prefix = NULL, ref_num = NULL, snp_num = 4, labs.cex = 1, ref = NULL, plot_subset = NULL, rev_compl = FALSE, snp_type = NULL, sep_lines = c(1.5, 4.5), legend.fill = c("#008000","blue","black","red","#d3d3d3"), legend.txt = c("A","C","G","T","Invariant"), genestart = 0, lm.col = 0.8, lwd = 0.5, rev_compl_sites = NULL, bonferroni_threshold = NULL, plot.dim = c(100,40), legend.cex = 1, legend.xpos = NULL, legend.ypos = NULL, text.xpos = NULL, text.ypos = NULL, axis.cex = 1, legend.lty = NULL, legend.pch = NULL, legend.col = NULL, line.sep.lwd = 0.5, legend.odds = TRUE, ref.name = "REF", reverse.xaxis = FALSE, reverse.xaxis.start = NULL, kmer.numbers = NULL, manhattan.order = FALSE, beta_estimate = NULL, se.d.kmers = NULL, beta_cols_list = c("blue","green","orange","red"), close.file = TRUE, margins = c(2,3,7.5,7), which_translucent = NULL, plot_ref = TRUE, forward.xaxis.start = NULL, col_lib = NULL){
	# Build matrix filled in by mapped kmers
	kmer_matrix = build_kmer_matrix(kmers = kmerseq, genestart = genestart, ref = ref, snps = snp_cols, kmerpos = kmerpos, snp_num = snp_num)
	# Get kmer colours by their odds ratio or beta point estimate
	if(!is.null(effect_size)){
		kmer_odds_COL = get_kmer_cols_lm(effect_size, lm.col)
	} else if(!is.null(beta_estimate)){
		kmer_odds_COL = apply(beta_estimate, 1, function(x, cols, se) get_betaCOL_2cols(x, cols, se), cols = beta_cols_list, se = se.d.kmers)
	} else {
		kmer_odds_COL = rep("#d3d3d3", length(kmerseq))
	}
	# Get SNP alignment colours
	# How many lines not to be coloured
	len.col = ref_num+snp_num
	#cat("len.col:",len.col,"dimkmermatrix:",dim(kmer_matrix),"lenkmeroddscol",length(kmer_odds_COL),"\n")
	alignment_reads_col = get_alignment_col(seq = kmer_matrix, cols = kmer_odds_COL, len.col = len.col, snp_cols = snp_cols, which_translucent = which_translucent, col_lib = col_lib)
	if(!is.null(plot_subset)){
		kmer_matrix = kmer_matrix[1:(ref_num+snp_num+plot_subset),]
		alignment_reads_col = alignment_reads_col[1:(ref_num+snp_num+plot_subset),]
		#kmer_odds_COL = kmer_odds_COL[1:plot_subset]
		if(length(which(effect_size[1:plot_subset]!="Inf"))!=0){
			max.odds = round(max(effect_size[(1:plot_subset)[which(effect_size[1:plot_subset]!="Inf")]]))
		} else {
			max.odds = "NA"
		}
	} else {
		if(length(which(effect_size!="Inf"))!=0){
			max.odds = round(max(effect_size[which(effect_size!="Inf")]))
		} else {
			max.odds = "NA"
		}
	}
	
	# Fill in SNP rows with mapped data SNP colours
	# New line
	if(!is.null(snp_cols)) alignment_reads_col[((1:nrow(snp_cols))+ref_num),] = snp_cols
	if(is.null(snp_cols)) n_snp_cols = 0 else n_snp_cols = nrow(snp_cols)
	if(!is.null(snp_type)) alignment_reads_col[(ref_num+n_snp_cols+1):(ref_num+snp_num),] = snp_type
	if(rev_compl & !is.null(snp_cols)){
		
		# alignment_reads_col[((1:nrow(snp_cols))+ref_num),] = matrix(unlist(rev_compl_col[snp_cols[,ncol(snp_cols):1]]),nrow = snp_num)
		alignment_reads_col[((1:nrow(snp_cols))+ref_num),c(rev_compl_sites)] = matrix(unlist(rev_compl_col[snp_cols[,c(rev_compl_sites)]]),nrow = nrow(snp_cols))
		alignment_reads_col[(1:ref_num),c(rev_compl_sites)] = matrix(unlist(rev_compl_col[alignment_reads_col[(1:ref_num),c(rev_compl_sites)]]),nrow = ref_num)
		# if(!is.null(snp_type)){
			# alignment_reads_col[ref_num+snp_num,] = rev(snp_type)
		# }
	} #else {
		#alignment_reads_col[((1:nrow(snp_cols))+ref_num),] = snp_cols
		#if(!is.null(snp_type)){
			#cat("refsumsnpnum:",ref_num+snp_num,"dimsnptype:",length(snp_type),"dimaligncol:",dim(alignment_reads_col),"\n")
			#alignment_reads_col[ref_num+snp_num,] = snp_type
		#}
	#}
	
	# If want to show x-axis decreasing rather than increasing
	if(reverse.xaxis){
		
		kmer_matrix = kmer_matrix[,c(ncol(kmer_matrix):1)]
		alignment_reads_col = alignment_reads_col[,c(ncol(alignment_reads_col):1)]
		
	}
	
	if(manhattan.order){
		kmer_matrix = kmer_matrix[c((1:(ref_num+snp_num)), (nrow(kmer_matrix):(ref_num+snp_num+1))),]
		alignment_reads_col = alignment_reads_col[c((1:(ref_num+snp_num)), (nrow(alignment_reads_col):(ref_num+snp_num+1))),]
	}
	
	alignment_reads_col_unadjusted = alignment_reads_col
	if(plot_ref==FALSE){
		alignment_reads_col[1:ref_num,] = "#ffffff"
	}
	# Plot alignment
	plot_alignment(prefix = prefix, seq = kmer_matrix, align_col = alignment_reads_col, ref_num = ref_num, snp_num = snp_num, labs.cex = labs.cex, legend.txt = legend.txt, legend.fill = legend.fill, sep_lines = sep_lines, lm.col = lm.col, max.odds = max.odds, genestart = genestart, bonferroni_threshold = bonferroni_threshold, n_snp_cols = n_snp_cols, plot.dim = plot.dim, legend.cex = legend.cex, legend.xpos = legend.xpos, legend.ypos = legend.ypos, text.xpos = text.xpos, text.ypos = text.ypos, axis.cex = axis.cex, legend.lty = legend.lty, legend.pch = legend.pch, legend.col = legend.col, line.sep.lwd = line.sep.lwd, legend.odds = legend.odds, ref.name = ref.name, reverse.xaxis = reverse.xaxis, reverse.xaxis.start = reverse.xaxis.start, kmer.numbers = kmer.numbers, close.file = close.file, margins = margins, plot_ref = plot_ref, forward.xaxis.start = forward.xaxis.start)
	return(alignment_reads_col_unadjusted)
	
}



st = function(x) (x-min(x,na.rm=TRUE))/diff(range(x,na.rm=TRUE))


# For if plotting with snp type legend labels
get_legend_col_alignment = function(lm.col = 0.8, max.odds = NULL, legend.xpos = NULL, legend.ypos = NULL, text.xpos = NULL, text.ypos = NULL){
	#lm = 0.8
	lm = lm.col
	testcol1 = seq(from = 1,by = -0.01, length.out=100)^0.99
	testcol1 = rgb(lm*(1-testcol1), lm*(1-testcol1), lm+testcol1*(1-lm))

	testcol2 = seq(from = 0.01,by = 0.01, length.out=100)^0.99
	testcol2 = rgb(lm+testcol2*(1-lm), lm*(1-testcol2), lm*(1-testcol2))
	
	# cat("xpos:",c(legend.xpos[1], ((legend.xpos[2]-legend.xpos[1])/2)+legend.xpos[1]),"\n")
	# cat("ypos:",c(legend.ypos[1], legend.ypos[2]),"\n")
	
	par(fig = c(legend.xpos[1], ((legend.xpos[2]-legend.xpos[1])/2)+legend.xpos[1], legend.ypos[1], legend.ypos[2]), mar=c(0,0,0,0), new=TRUE)
	image(c(1:100),1, (matrix(c(1:100),ncol=1,nrow=100)), col=testcol1, axes=FALSE)
	axis(1, at = c(1,100), labels = c(NA,NA), xpd = T, tck = -0.1, cex.axis = 0.5, lwd = 0.8)
	axis(1, at = c(1,100), labels = c(0,1), xpd = T, tck = 0, cex.axis = 0.5, lwd = 0, line = -1.3)

	# par(fig = c(0.9723, 0.9727, 0.949, 0.954), mar=c(0,0,0,0), new=TRUE)
	par(fig = c(0,1,0,1), mar=c(0,0,0,0), new = TRUE)
	# plot(c(1:10),c(1:10))
	plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
	# text(labels = c("0"),cex = 0.6, x = c(text.xpos[1]), y = text.ypos, xpd = T)
	# text(labels = c("1"),cex = 0.6, x = c(text.xpos[2]), y = text.ypos, xpd = T)
	# text(labels = c(max.odds),cex = 0.6, x = c(text.xpos[3]), y = text.ypos, xpd = T)
	
	par(fig = c(0.9762, 0.9767,0.949, 0.954), mar=c(0,0,0,0), new=TRUE)

	# mtext(side=1,"1", line = -1.1,adj = 0, cex = 0.7)
	
	
	# cat("xpos:",c(((legend.xpos[2]-legend.xpos[1])/2)+legend.xpos[1], legend.xpos[2]),"\n")
	# cat("ypos:",c(legend.ypos[1], legend.ypos[2]),"\n")

	par(fig = c(((legend.xpos[2]-legend.xpos[1])/2)+legend.xpos[1], legend.xpos[2],legend.ypos[1], legend.ypos[2]), mar=c(0,0,0,0), new=TRUE)
	image(c(1:100),1, (matrix(c(1:100),ncol=1,nrow=100)), col=testcol2, axes=FALSE)
	axis(1, at = c(1,100), labels = c(NA,NA), xpd = T, tck = -0.1, cex.axis = 0.5, lwd = 0.8)
	# cat("max.odds:",max.odds,"\n")
	# cat("is.na:",is.na(max.odds),"\n")
	if(max.odds=="NA") max.odds = "Inf"
	axis(1, at = c(1,100), labels = c(NA,max.odds), xpd = T, tck = 0, cex.axis = 0.5, lwd = 0, line = -1.3)
	# axis(1, at = c(1,100), labels = c(0,1), xpd = T, tck = 0, cex.axis = 0.5, lwd = 0, line = -1.3)
	par(fig = c(0.9802, 0.9807,0.949, 0.954), mar=c(0,0,0,0), new=TRUE)
	# mtext(side=1,max.odds, line = -1.1,adj = 0, cex = 0.7)
	
	
}

run_manhattan_allframes = function(gene_i_results_list = NULL, prefix = NULL, gene_name = NULL, ref_gene_i = NULL, which_kmers_no_result = NULL, bonferroni = NULL, minor_allele_threshold = NULL, macormaf = NULL, output_dir = NULL, kmer_type = NULL, kmer_length = NULL, ref.name = NULL, ref_gb_full = NULL){
	
	prefix = paste0(output_dir, prefix, "_", kmer_type, kmer_length, "_", ref.name)
	
	ymax = plot_allframes_manhattan(gene_i_results_list = gene_i_results_list, prefix = prefix, genes_name = gene_name, length_correct = ref_gene_i$length_correct, all_translations = ref_gene_i$all_translations, correct_frame = ref_gene_i$correct_frame, which_kmers_no_result = which_kmers_no_result, bonferroni = bonferroni, ylim = NULL, x.adjust = 333, gene_end_line = ref_gene_i$length_protein, start = ref_gene_i$ref_start_i, end = ref_gene_i$ref_end_i, macormaf = macormaf, ref_gb_full = ref_gb_full)
	if(ymax>100){
		ymax = plot_allframes_manhattan(gene_i_results_list = gene_i_results_list, prefix = prefix, genes_name = gene_name, length_correct = ref_gene_i$length_correct, all_translations = ref_gene_i$all_translations, correct_frame = ref_gene_i$correct_frame, which_kmers_no_result = which_kmers_no_result, bonferroni = bonferroni, ylim = 50, x.adjust = 333, gene_end_line = ref_gene_i$length_protein, start = ref_gene_i$ref_start_i, end = ref_gene_i$ref_end_i, macormaf = macormaf, ref_gb_full = ref_gb_full)
	}
	ymax = plot_allframes_manhattan(gene_i_results_list = gene_i_results_list, prefix = prefix, genes_name = gene_name, length_correct = ref_gene_i$length_correct, all_translations = ref_gene_i$all_translations, correct_frame = ref_gene_i$correct_frame, which_kmers_no_result = which_kmers_no_result, bonferroni = bonferroni, ylim = NULL, maname = paste0("_", macormaf, minor_allele_threshold), malim = minor_allele_threshold, x.adjust = 333, gene_end_line = ref_gene_i$length_protein, start = ref_gene_i$ref_start_i, end = ref_gene_i$ref_end_i, macormaf = macormaf, ref_gb_full = ref_gb_full)
	if(ymax>100){
		ymax = plot_allframes_manhattan(gene_i_results_list = gene_i_results_list, prefix = prefix, genes_name = gene_name, length_correct = ref_gene_i$length_correct, all_translations = ref_gene_i$all_translations, correct_frame = ref_gene_i$correct_frame, which_kmers_no_result = which_kmers_no_result, bonferroni = bonferroni, ylim = 50, maname = paste0("_", macormaf, minor_allele_threshold), malim = minor_allele_threshold, x.adjust = 333, gene_end_line = ref_gene_i$length_protein, start = ref_gene_i$ref_start_i, end = ref_gene_i$ref_end_i, macormaf = macormaf, ref_gb_full = ref_gb_full)
	}

}

run_manhattan_single_protein = function(which_kmers_no_result = NULL, res = NULL, ref_gene_i = NULL, kmer_length = NULL, prefix = NULL, gene_name = NULL, j = NULL, bonferroni = NULL, ref_gb_full = NULL, ref_length = NULL, kmer_type = NULL, nsamples = NULL, minor_allele_threshold = NULL, macormaf = NULL, output_dir = NULL, ref.name = NULL){

	if(j==ref_gene_i$correct_frame) correct_or_wrong = "correct_frame" else correct_or_wrong = "wrong_frame"
	
	# Plot Manhattan for the gene for each reading frame
	# Give position to unmapped kmers
	xpos = cbind(as.numeric(res$sstart), as.numeric(res$send))
	ypos = as.numeric(res$negLog10)
	beta = as.numeric(res$beta)
	whichMAthreshold = which(c(as.numeric(res[[macormaf]]))>=(minor_allele_threshold))
	if(!is.null(which_kmers_no_result)){
		ypos = c(ypos, as.numeric(which_kmers_no_result$negLog10))
		xpos_unmapped = sample(seq(from = ref_gene_i$length_correct+50, to = ref_gene_i$length_correct+100, length.out = nrow(which_kmers_no_result)), nrow(which_kmers_no_result), replace = F)
		xpos_unmapped = cbind(xpos_unmapped, xpos_unmapped+kmer_length-1)
		xpos = rbind(xpos, xpos_unmapped)
		beta = c(beta, as.numeric(which_kmers_no_result$beta))
		whichMAthreshold = which(c(as.numeric(res[[macormaf]]), as.numeric(which_kmers_no_result[[macormaf]]))>=(minor_allele_threshold))
	}
	beta_col = rep("#d3d3d3", length(beta)); beta_col[which(beta>0)] = "#838383"
	
	# Change colour of unaligned to red rather than their beta colour
	beta_col[(nrow(res)+1):length(beta_col)][which(beta_col[(nrow(res)+1):length(beta_col)]=="#d3d3d3")] = "#ffbc87"
	beta_col[(nrow(res)+1):length(beta_col)][which(beta_col[(nrow(res)+1):length(beta_col)]=="#838383")] = "#D55E00"
	
	prefix = paste0(output_dir, prefix, "_", kmer_type, kmer_length, "_", ref.name)
	
	# First with no limit on y-axis
	plot_singleframe_manhattan_protein(prefix = prefix, genes_name = gene_name, correct_or_wrong = correct_or_wrong, j = j, xpos = xpos, ypos = ypos, all_translations = ref_gene_i$all_translations, beta_col = beta_col, bonferroni = bonferroni, ylim_max = NULL, x.adjust = 999, ref_gb_full = ref_gb_full, start = ref_gene_i$ref_start_i, end = ref_gene_i$ref_end_i, ref_length = ref_length, kmer_length = kmer_length)
	plot_singleframe_manhattan_protein(prefix = prefix, genes_name = gene_name, correct_or_wrong = correct_or_wrong, j = j, xpos = xpos[whichMAthreshold,], ypos = ypos[whichMAthreshold], all_translations = ref_gene_i$all_translations, beta_col = beta_col[whichMAthreshold], bonferroni = bonferroni, ylim_max = NULL, maname = paste0("_",macormaf,minor_allele_threshold), x.adjust = 999, ref_gb_full = ref_gb_full, start = ref_gene_i$ref_start_i, end = ref_gene_i$ref_end_i, ref_length = ref_length, kmer_length = kmer_length)

	# Then limit to ylim = 50 for those with ylim > 100
	if(max(as.numeric(ypos))>100){
		plot_singleframe_manhattan_protein(prefix = prefix, genes_name = gene_name, correct_or_wrong = correct_or_wrong, j = j, xpos = xpos, ypos = ypos, all_translations = ref_gene_i$all_translations, beta_col = beta_col, bonferroni = bonferroni, ylim_max = 50, x.adjust = 999, ref_gb_full = ref_gb_full, start = ref_gene_i$ref_start_i, end = ref_gene_i$ref_end_i, ref_length = ref_length, kmer_length = kmer_length)
	}
	if(max(as.numeric(ypos)[whichMAthreshold])>100){
		plot_singleframe_manhattan_protein(prefix = prefix, genes_name = gene_name, correct_or_wrong = correct_or_wrong, j = j, xpos = xpos[whichMAthreshold,], ypos = ypos[whichMAthreshold], all_translations = ref_gene_i$all_translations, beta_col = beta_col[whichMAthreshold], bonferroni = bonferroni, ylim_max = 50, maname = paste0("_",macormaf,minor_allele_threshold), x.adjust = 999, ref_gb_full = ref_gb_full, start = ref_gene_i$ref_start_i, end = ref_gene_i$ref_end_i, ref_length = ref_length, kmer_length = kmer_length)
	}

}

run_manhattan_single_nucleotide = function(which_kmers_no_result = NULL, res = NULL, ref_gene_i = NULL, prefix = NULL, gene_name = NULL, bonferroni = NULL, ref_gb_full = NULL, ref_length = NULL, kmer_type = NULL, kmer_length = NULL, nsamples = NULL, minor_allele_threshold = NULL, macormaf = NULL, output_dir = NULL, ref.name = NULL){
	
	ref_start_i = ref_gene_i$ref_start_i
	ref_end_i = ref_gene_i$ref_end_i
	
	# Plot Manhattan for the gene for each reading frame
	# Give position to unmapped kmers and plot those on every reading frame
	xpos = cbind(as.numeric(res$sstart), as.numeric(res$send))
	ypos = as.numeric(res$negLog10)
	beta = as.numeric(res$beta)
	whichMAthreshold = which(c(as.numeric(res[[macormaf]]))>=(minor_allele_threshold))
	if(!is.null(which_kmers_no_result)){
		ypos = c(ypos, as.numeric(which_kmers_no_result$negLog10))
		xpos_unmapped = sample(seq(from = nchar(ref_gene_i$ref_gene_i)+50, to = nchar(ref_gene_i$ref_gene_i)+100, length.out = nrow(which_kmers_no_result)), nrow(which_kmers_no_result), replace = F)
		xpos_unmapped = cbind(xpos_unmapped, xpos_unmapped+(kmer_length-1))
		xpos = rbind(xpos, xpos_unmapped)
		beta = c(beta, as.numeric(which_kmers_no_result$beta))
		whichMAthreshold = which(c(as.numeric(res[[macormaf]]), as.numeric(which_kmers_no_result[[macormaf]]))>=(minor_allele_threshold))
	}
	beta_col = rep("#d3d3d3", length(beta)); beta_col[which(beta>0)] = "#838383"
	# Change colour of unaligned to red rather than their beta colour
	if(!is.null(which_kmers_no_result)){
		beta_col[(nrow(res)+1):length(beta_col)][which(beta_col[(nrow(res)+1):length(beta_col)]=="#d3d3d3")] = "#ffbc87"
		beta_col[(nrow(res)+1):length(beta_col)][which(beta_col[(nrow(res)+1):length(beta_col)]=="#838383")] = "#D55E00"
	}
	
	prefix = paste0(output_dir, prefix, "_", kmer_type, kmer_length, "_", ref.name)
	# First with no limit on y-axis
	plot_singleframe_manhattan_nucleotide(prefix = prefix, genes_name = gene_name, xpos = xpos, ypos = ypos, ref_fa = ref_gene_i$ref_gene_i, beta_col = beta_col, bonferroni = bonferroni, ylim_max = NULL, x.adjust = 999, ref_gb_full = ref_gb_full, start = ref_start_i, end = ref_end_i)
	plot_singleframe_manhattan_nucleotide(prefix = prefix, genes_name = gene_name, xpos = xpos[whichMAthreshold,], ypos = ypos[whichMAthreshold], ref_fa = ref_gene_i$ref_gene_i, beta_col = beta_col[whichMAthreshold], bonferroni = bonferroni, ylim_max = NULL, maname = paste0("_",macormaf,minor_allele_threshold), x.adjust = 999, ref_gb_full = ref_gb_full, start = ref_start_i, end = ref_end_i)
	
	if(max(as.numeric(ypos))>100){
		plot_singleframe_manhattan_nucleotide(prefix = prefix, genes_name = gene_name, xpos = xpos, ypos = ypos, ref_fa = ref_gene_i$ref_gene_i, beta_col = beta_col, bonferroni = bonferroni, ylim_max = 50, x.adjust = 999, ref_gb_full = ref_gb_full, start = ref_start_i, end = ref_end_i)
	}
	if(max(as.numeric(ypos)[whichMAthreshold])>100){
		plot_singleframe_manhattan_nucleotide(prefix = prefix, genes_name = gene_name, xpos = xpos[whichMAthreshold,], ypos = ypos[whichMAthreshold], ref_fa = ref_gene_i$ref_gene_i, beta_col = beta_col[whichMAthreshold], bonferroni = bonferroni, ylim_max = 50, maname = paste0("_",macormaf,minor_allele_threshold), x.adjust = 999, ref_gb_full = ref_gb_full, start = ref_start_i, end = ref_end_i)
	}


	
}



read_res_table = function(resfile = NULL, refseq = NULL, i = NULL, genes_names = NULL, perident = 70, kmer_type = NULL, kmer_length = NULL){

	# if(j==ref_gene_i$correct_frame) correct_or_wrong = "correct_frame" else correct_or_wrong = "wrong_frame"
	# res = read.table(paste0(prefix, "_", drug, "_blast_results_top_gene_",i,"_", drug_genes_names[i],"_", j, "_", correct_or_wrong, ".txt"), h = T, sep = "\t")
	res = read.table(resfile, h = T, sep = "\t")
	res = res[order(as.numeric(res$negLog10), decreasing = T),]
	# Fix res so that the full kmer is aligned, not just plotting the BLAST result
	if(kmer_type=="protein"){
		res = fix_res_table_protein(res = res, refseq = refseq)
	} else {
		res = fix_res_table_nucleotide(res = res, ref_fa = refseq, kmer_len = kmer_length)
	}
	# Now redo % identity
	nmatch = sapply(1:nrow(res), function(p, x, y) length(which(unlist(strsplit(as.character(x[p]),""))==unlist(strsplit(as.character(y[p]),"")))), x = res$sseq, y = res$qseq, USE.NAMES = F)
	if(nchar(refseq)>nchar(res$kmer[1])){
		low_match_results_j = res[which(nmatch<(nchar(res$kmer[1])*(perident/100))),]
		res = res[which(nmatch>=(nchar(res$kmer[1])*(perident/100))),]
	} else {
		low_match_results_j = res[which(nmatch<(nchar(refseq)*(perident/100))),]
		res = res[which(nmatch>=(nchar(refseq)*(perident/100))),]
	}

	return(list("res" = res, "low_match_results_j" = low_match_results_j))
}

run_alignment_nplots_protein = function(ref_gene_i = NULL, res = NULL, nsamples = NULL, bonferroni = NULL, prefix = NULL, gene_name = NULL, j = NULL, override_signif = FALSE, prange = NULL, plot_figures = TRUE, col_lib = NULL, minor_allele_threshold = NULL, macormaf = NULL, output_dir = NULL, kmer_type = NULL, kmer_length = NULL, ref.name = NULL){

	if(j==ref_gene_i$correct_frame) correct_or_wrong = "correct_frame" else correct_or_wrong = "wrong_frame"

	# Plot in a sliding window across the protein
	nplots = seq(from = 1, by = 20, to = nchar(ref_gene_i$all_translations[j]))


	nplots = cbind(nplots, nplots+39)
	if(any(nplots[,2]>nchar(ref_gene_i$all_translations[j]))){
		nplots[which(nplots[,2]>nchar(ref_gene_i$all_translations[j])),2] = nchar(ref_gene_i$all_translations[j])
	}

	if(nrow(nplots)>1){
		if(length(nplots[nrow(nplots),1]:nplots[nrow(nplots),2])<max(nchar(as.character(res$kmer))) & nplots[(nrow(nplots)-1),2]==nchar(ref_gene_i$all_translations[j])) nplots = nplots[-nrow(nplots),]
	}
	if(is.null(nrow(nplots))){
		nplots = matrix(nplots, ncol = 2)
	}

	if(is.null(prange)) prange = 1:nrow(nplots)

	out_results_correct_frame = list()

	for(p in prange){
		
		out_mafthreshold = plot_alignment_function_protein(genestart = nplots[p,1], geneend = nplots[p,2], minor_allele_threshold = minor_allele_threshold, res = res, nsamples = nsamples, maname = paste0("_", macormaf,minor_allele_threshold), bonferroni = bonferroni, translation = ref_gene_i$all_translations[j], prefix = prefix, gene_name = gene_name, correct_or_wrong = correct_or_wrong, j = j, p = p, plot_ref = FALSE, main = paste0("Protein kmers \u2265 ",macormaf," ", minor_allele_threshold), x.adjust = 333, nplots = nplots, override_signif = override_signif, plot_figures = plot_figures, col_lib = col_lib, macormaf = macormaf, output_dir = output_dir, kmer_type = kmer_type, kmer_length = kmer_length, ref.name = ref.name)
		if(!is.null(out_mafthreshold) | override_signif){
			out_allmaf = plot_alignment_function_protein(genestart = nplots[p,1], geneend = nplots[p,2], minor_allele_threshold = 0, res = res, nsamples = nsamples, maname = "", bonferroni = bonferroni, translation = ref_gene_i$all_translations[j], prefix = prefix, gene_name = gene_name, correct_or_wrong = correct_or_wrong, j = j, p = p, lowfreq = minor_allele_threshold, plot_ref = FALSE, main = "All protein kmers", x.adjust = 333, nplots = nplots, override_signif = override_signif, plot_figures = plot_figures, col_lib = col_lib, macormaf = macormaf, output_dir = output_dir, kmer_type = kmer_type, kmer_length = kmer_length, ref.name = ref.name)
			if(!is.null(out_allmaf)) out_results_correct_frame = rbind(out_results_correct_frame, out_allmaf)
		}
		

	}
	if(j==ref_gene_i$correct_frame) return(out_results_correct_frame)

}

run_alignment_nplots_nucleotide = function(ref_gene_i = NULL, res = NULL, nsamples = NULL, bonferroni = NULL, prefix = NULL, gene_name = NULL, override_signif = FALSE, prange = NULL, gene_lookup = NULL, wh_genelookup = NULL, col_lib = NULL, plot_figures = TRUE, minor_allele_threshold = NULL, macormaf = NULL, output_dir = NULL, kmer_type = NULL, kmer_length = NULL, ref.name = NULL){

	ref_start_i = ref_gene_i$ref_start_i
	ref_end_i = ref_gene_i$ref_end_i
	ref_gene_i = ref_gene_i$ref_gene_i
	
	nplots = seq(from = 1, by = 40, to = nchar(ref_gene_i))
	nplots = cbind(nplots, nplots+99)
	if(any(nplots[,2]>nchar(ref_gene_i))){
		nplots[which(nplots[,2]>nchar(ref_gene_i)),2] = nchar(ref_gene_i)
	}
	if(nrow(nplots)>1){
		if(length(nplots[nrow(nplots),1]:nplots[nrow(nplots),2])<max(nchar(as.character(res$kmer))) & nplots[(nrow(nplots)-1),2]==nchar(ref_gene_i)) nplots = nplots[-nrow(nplots),]
	}
	if(is.null(nrow(nplots))){
		nplots = matrix(nplots, ncol = 2)
	}

	if(is.null(prange)) prange = 1:nrow(nplots)
	
	allgenes_results_table_out = list()

	for(p in prange){
		
		rev_xaxis = as.numeric(gene_lookup[wh_genelookup,5])!=1
		if(rev_xaxis){
			reverse.xaxis.start = (ref_end_i:ref_start_i)[nplots[p,1]]
			forward.xaxis.start = NULL
		} else {
			reverse.xaxis.start = NULL
			forward.xaxis.start = c(ref_start_i:ref_end_i)[nplots[p,1]]
		}
		rev_xaxis = FALSE
		
		
		out_mafthreshold = plot_alignment_function_nucleotide(genestart = nplots[p,1], geneend = nplots[p,2],
									minor_allele_threshold = minor_allele_threshold, res = res, nsamples = nsamples, maname = paste0("_", macormaf,minor_allele_threshold),
									bonferroni = bonferroni, ref_fa = ref_gene_i, prefix = prefix,
									gene_name = gene_name, p = p, plot_ref = FALSE,
									main = paste0("Nucleotide kmers \u2265 ",macormaf," ", minor_allele_threshold), x.adjust = 999,
									reverse.xaxis = rev_xaxis, reverse.xaxis.start = reverse.xaxis.start,
									forward.xaxis.start = forward.xaxis.start, plot_figures = plot_figures,
									alignment_range = c(ref_start_i, ref_end_i), override_signif = override_signif,
									col_lib = col_lib, macormaf = macormaf, output_dir = output_dir,
									kmer_type = kmer_type, kmer_length = kmer_length, ref.name = ref.name)
		
		if(!is.null(out_mafthreshold) | override_signif){
		
			out_allmaf = plot_alignment_function_nucleotide(genestart = nplots[p,1], geneend = nplots[p,2],
									minor_allele_threshold = 0, res = res, nsamples = nsamples, maname = "",
									bonferroni = bonferroni, ref_fa = ref_gene_i,
									prefix = prefix, gene_name = gene_name,
									p = p, lowfreq = minor_allele_threshold, plot_ref = FALSE, main = "All nucleotide kmers",
									x.adjust = 999, reverse.xaxis = rev_xaxis,
									reverse.xaxis.start = reverse.xaxis.start,
									forward.xaxis.start = forward.xaxis.start, plot_figures = plot_figures,
									alignment_range = c(ref_start_i, ref_end_i), override_signif = override_signif,
									col_lib = col_lib, macormaf = macormaf, output_dir = output_dir,
									kmer_type = kmer_type, kmer_length = kmer_length, ref.name = ref.name)

			if(!is.null(out_allmaf)) allgenes_results_table_out = rbind(allgenes_results_table_out, out_allmaf)
		}

	}

	return(allgenes_results_table_out)

}

plot_alignment_function_nucleotide = function(genestart = NULL, geneend = NULL, minor_allele_threshold = NULL, res = NULL, nsamples = NULL, maname = NULL, bonferroni = NULL, ref_fa = NULL, prefix = NULL, gene_name = NULL, p = NULL, lowfreq = NULL, plot_ref = NULL, main = NULL, x.adjust = 0, reverse.xaxis = NULL, reverse.xaxis.start = NULL, forward.xaxis.start = NULL, override_signif = FALSE, plot_figures = TRUE, alignment_range = NULL, col_lib = NULL, macormaf = NULL, output_dir = NULL, kmer_type = NULL, kmer_length = NULL, ref.name = NULL){
	
	
	alignment_start = as.numeric(alignment_range[1]); alignment_end = as.numeric(alignment_range[2])

	legend.txt = c(names(col_lib_nuc),"Invariant","","\u03B2 < 0","\u03B2 > 0","","Significance","threshold")
	legend.txt[which(legend.txt=="-")] = "Gap"
	# legend.fill = c(unlist(col_lib_nuc),"#d3d3d3",NA, "#d3d3d3", "#838383", NA, NA, NA, NA, NA)[which_legend]
	# legend.col = c(rep(NA,length(col_lib_nuc)+5),"black",NA, "black", NA)[which_legend]
	legend.col = c(unlist(col_lib_nuc),"#d3d3d3",NA, "#d3d3d3", "#838383", NA, "black", NA)
	legend.pch = c(rep(15,length(col_lib_nuc)+5), NA,NA)
	legend.lty = c(rep(NA, length(col_lib_nuc)+5), 2,NA)
	legend.pch.cex = rep(1, length(legend.txt))
	sep_lines = c(0.5)
	ref_num = 1

	# Get the position in the reference genome of the plot
	if(is.null(reverse.xaxis.start)){
		plot_start_position = forward.xaxis.start
		plot_end_position = forward.xaxis.start+(length(genestart:geneend))-1
	}  else {
		plot_start_position = reverse.xaxis.start
		plot_end_position = reverse.xaxis.start-(length(genestart:geneend))+1
	}

	# Which of the BLAST results should be plotted (which are in the region and pass the minor allele threshold)
	which_to_align = which(as.numeric(res$send)>=genestart & as.numeric(res$sstart)<=geneend & as.numeric(res[[macormaf]])>=(minor_allele_threshold))

	# If wanting to show which are low frequency (if they are in the plot)
	# then find which are below the cutoff and feed them in to be coloured translucent
	if(!is.null(lowfreq)){
		which_low_freq = which(c(as.numeric(res$mac)<(nsamples*lowfreq))[which_to_align])
	} else {
		which_low_freq = NULL
	}

	if(length(which_to_align)>0){
		# How many kmers are above the bonferroni threshold
		bonferroni_lim = length(which(as.numeric(res$negLog10)[which_to_align]<bonferroni))
		if(bonferroni_lim!=length(which_to_align) | override_signif){
			# Pull out reference amino acid sequence for the region
			gene_db_i = matrix(unlist(strsplit(ref_fa, ""))[genestart:geneend], nrow = 1)
			# gene_db_i = rbind(gene_db_i, gene_db_i, gene_db_i, gene_db_i, gene_db_i, gene_db_i)

			snp_num = round((length(which_to_align)+ref_num)*0.05)

			beta_to_input = cbind(rep(0, nrow(res)), as.numeric(res$beta), rep(0, nrow(res)))[which_to_align,]
			if(is.null(dim(beta_to_input))) beta_to_input = matrix(beta_to_input, nrow = 1)

			if(plot_figures){
				# Plot name was "_pos_",(genestart-x.adjust),"_to_",(geneend-x.adjust)
				aligncols = run_plot_alignment(kmerseq = as.character(res$qseq)[which_to_align], snp_cols = NULL,
									   kmerpos = as.numeric(res$sstart)[which_to_align],
									   prefix = paste0(output_dir, prefix, "_", kmer_type, kmer_length, "_",
									   ref.name, "_", gene_name, "_plot_",p,"_pos_",
									   plot_start_position,"_to_", plot_end_position, maname),
									   ref_num = ref_num, snp_num = snp_num, labs.cex = 0.8, ref = gene_db_i,
								   	   plot_subset = NULL, snp_type = NULL, legend.txt = NULL, legend.fill = legend.fill,
								   	   sep_lines = sep_lines, genestart = genestart, rev_compl=FALSE,
								   	   bonferroni_threshold = bonferroni_lim, plot.dim = c(21,18), legend.cex = 0.7,
								       legend.xpos = c(0.868, 0.9), legend.ypos = c(0.703,0.717),
								       text.xpos = c(0.794,0.830, 0.869), text.ypos = c(0.42), axis.cex = 0.7,
								       legend.lty = legend.lty, legend.pch = legend.pch, legend.col= legend.col,
								       line.sep.lwd = 0.15, reverse.xaxis = reverse.xaxis, manhattan.order = TRUE,
								       reverse.xaxis.start = reverse.xaxis.start,
								       beta_estimate = beta_to_input,
								       se.d.kmers = 0, beta_cols_list = c("#d3d3d3", "#009E73", "#E69F00", "#838383"),
								       ref.name = "", legend.odds = FALSE, close.file = FALSE, margins = c(2,3,2.5,6),
								       which_translucent = NULL, plot_ref = plot_ref,
								       forward.xaxis.start = forward.xaxis.start,
								       col_lib = col_lib)

				# Plot reference colours along bottom
				refcols_ytop = ((nrow(aligncols))*0.04)+0.5
				for(k in 1:ncol(aligncols)){
					rect(xleft = k-0.5, xright = k+0.5, ybottom = 0.5, ytop = refcols_ytop, col = aligncols[1,k], border = NA, xpd = TRUE)
				}
				# Add lines and reference name
				sapply(c(1:ncol(aligncols)), function(k, ytop, lwd) lines(x = c(k-0.5,k-0.5), y = c(0.5, ytop), lwd = lwd, col = "white"), lwd = 0.15, ytop = refcols_ytop, USE.NAMES = F)
				text(x = 1,y = mean(c(refcols_ytop, 0.5)),"REF",xpd=TRUE, pos = 2, cex = 0.8, adj = 1)


				add_xaxis_top(reverse.xaxis.start = reverse.xaxis.start, forward.xaxis.start = forward.xaxis.start, genestart = genestart, gene_db_i = gene_db_i)

				axis.label.pos = axis(3, at = 1:ncol(gene_db_i), labels = rep("", ncol(gene_db_i)), cex.axis = 0.4, lwd = NA, line = -0.75, xpd = T)
				for(k in axis.label.pos) axis(3, at = k, labels = as.character(gene_db_i[1,k]), cex.axis = 0.4, lwd = NA, line = -0.75, xpd = T)



				# Plot legend
				par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
				plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
				if(!is.null(main)) text(x = -0.99, y = 1.03, main, pos = 2, xpd = TRUE, cex = 0.8, srt = 90)
				# if(!is.null(main)) legend(x = -1.10, y = 1.03, main, col = "black", bg = "white", bty = "o",xjust = 0, yjust = 0.5, cex = 0.8, xpd = TRUE, box.col = "white")
				legend("topright", legend.txt, border = NA, bty = "o", cex = 0.7, col = legend.col, bg = NA, lty = legend.lty, pch = legend.pch, pt.cex = legend.pch.cex)
				dev.off()
			}

			# Create table with pvals, beta, maf etc. for all kmers plotted in the alignment figure
			cols_to_keep = c("kmer","negLog10","beta","mac","maf","sstart")
			if(override_signif){
				out_table = cbind(res[,match(cols_to_keep, names(res))][which_to_align,])
			} else {
				out_table = cbind(res[,match(cols_to_keep, names(res))][which_to_align[1:length(which(as.numeric(res$negLog10)[which_to_align]>=bonferroni))],])
			}
			colnames(out_table)[which(colnames(out_table)=="sstart")] = "ps"
			if(is.null(reverse.xaxis.start)){
				out_table$ps = sapply(as.numeric(out_table$ps), function(x, start, end) c(start:end)[x], start = alignment_start, end = alignment_end, USE.NAMES = F)
			}  else {
				# plot_start_position = reverse.xaxis.start
				# plot_end_position = reverse.xaxis.start-(length(genestart:geneend))+1
				out_table$ps = sapply(as.numeric(out_table$ps), function(x, start, end) c(end:start)[x], start = alignment_start, end = alignment_end, USE.NAMES = F)
			}

			out_table = cbind("gene" = gene_name, "plot" = paste0(gene_name, "_plot_",p,"_ps_",plot_start_position,"_to_",plot_end_position), out_table)
			return(out_table)
		}
	}
}



get_ref_gene_i = function(gene_lookup = NULL, genename_i = NULL, ref_length = NULL, ref_fa = NULL, oneLetterCodes = NULL){

	# Which row in gene_lookup contains genename_i
	wh_genelookup = which(gene_lookup[,1]==genename_i)
	if(length(wh_genelookup)==0) stop("Error: no match for gene name", genename_i, "\n")

	# # Create reference file for gene
	# Create reference file for gene - extend either side of the gene
	ref_start_i = (as.numeric(gene_lookup[wh_genelookup,3])-999)
	if(ref_start_i<1) ref_start_i = 1
	ref_end_i = (as.numeric(gene_lookup[wh_genelookup,4])+999)
	if(ref_end_i>ref_length) ref_end_i = ref_length

	ref_gene_i = paste(ref_fa[ref_start_i:ref_end_i], collapse = "")
	length_protein = length(as.numeric(gene_lookup[wh_genelookup,3]):as.numeric(gene_lookup[wh_genelookup,4]))/3

	# Get all translations
	all_translations = translate_6_frames_alignment(contig = ref_gene_i, oneLetterCodes = oneLetterCodes)
	# Which is the correct frame
	if(as.numeric(gene_lookup[wh_genelookup,5])==1) correct_frame = 1 else correct_frame = 4

	length_correct = nchar(all_translations[correct_frame])

	return(list("ref_start_i" = ref_start_i, "ref_end_i" = ref_end_i, "ref_gene_i" = ref_gene_i, "length_protein" = length_protein, "all_translations" = all_translations, "correct_frame" = correct_frame, "length_correct" = length_correct))

}

read_kmer_file = function(kmerfile = NULL, i = NULL, prefix = NULL, kmer_type = NULL, kmer_length = NULL, output_dir = NULL){

	# Read in kmers
	kmers_gene_i = read.table(kmerfile, h = T, sep = "\t")
	kmers_gene_i = unique(kmers_gene_i)
	kmers_gene_i = kmers_gene_i[order(as.numeric(kmers_gene_i$negLog10), decreasing = T),]

	kmers_gene_i_write = as.vector(rbind(paste0(">kmer", 1:nrow(kmers_gene_i)), as.character(kmers_gene_i$kmer)))
	### *** CHANGE TO TEMP ***
	cat(kmers_gene_i_write, file = paste0(output_dir, prefix, "_", kmer_type, kmer_length, "_gene_", i, "_tmp_kmer_file.txt"), sep = "\n")
	return(kmers_gene_i)

}

get_kmers_noresult = function(gene_i_results_list = NULL, kmers_gene_i = NULL, prefix = NULL, i = NULL, genes_names = NULL, kmer_type = NULL, kmer_length = NULL, output_dir = NULL, ref.name = NULL){

	# Find the kmers which haven't mapped to any reading frame
	all_kmer_matches = as.vector(unlist(sapply(gene_i_results_list, function(x) return(x$origkmer), USE.NAMES = F)))
	all_kmer_matches = unique(all_kmer_matches)
	which_kmers_no_result = which(is.na(match(as.character(kmers_gene_i$kmer), all_kmer_matches)))
	if(length(which_kmers_no_result)>0){
		which_kmers_no_result = kmers_gene_i[which_kmers_no_result, ]
		new_colnames = c("kmer","negLog10", "beta","mac")
		if(is.null(nrow(which_kmers_no_result))){
			which_kmers_no_result = matrix(which_kmers_no_result, nrow=1)
			if(ncol(which_kmers_no_result)!=length(new_colnames)) stop("Error: dimensions of which_kmers_no_result in get_kmers_noresult is incorrect for gene",gene_names[i], "\n")
		}
		colnames(which_kmers_no_result) = new_colnames
		which_kmers_no_result = which_kmers_no_result[order(as.numeric(which_kmers_no_result$negLog10), decreasing = T),]
		write.table(which_kmers_no_result, file = paste0(output_dir, prefix, "_", kmer_type, kmer_length, "_", ref.name, "_top_gene_",i,"_", genes_names[i], "_no_blast_result_or_poor_alignment.txt"), row = F, col = T, sep = "\t", quote = F)
		return(which_kmers_no_result)
	} else {
		return(NULL)
	}

}


draw_arrow = function(start = NULL, end = NULL, arrow_length = 0.1, height1 = NULL,
					  height2 = NULL, arrowdiff = NULL, fillCOL = "grey50",
					  border = "black", name = NULL, text_adjust = 0, rev.compl = FALSE, text.cex = 0.5, lwd = 1.5,
					  text.col = "black"){
	# If the gene is reverse complemented switch the start and end positions
	if(rev.compl){
		start.new = end
		end.new = start
		start = start.new
		end = end.new
	}
	# Draw the initial rectangle (without the arrow head)
	rect(start, height1, (end+((start-end)*arrow_length)), height2, col = fillCOL, border = fillCOL, xpd = T)
	# Draw the arrow head
	polygon(x = c((end+((start-end)*arrow_length)), (end+((start-end)*arrow_length)), end),
			y = c(height1-arrowdiff, height2+arrowdiff, height1+abs(height1-height2)/2), xpd = T,
			col = fillCOL, border = NA)
	# Surround the rectangle and arrow head with a border
	lines(x = c(start, (end+((start-end)*arrow_length))), y = c(height1, height1),
		  lwd = lwd, xpd = T, col = border)
	lines(x = c(start, (end+((start-end)*arrow_length))), y = c(height2, height2),
		  lwd = lwd, xpd = T, col = border)
	lines(x = c(start, start), y = c(height1, height2), lwd = lwd, xpd = T, col = border)
	lines(x = c((end+((start-end)*arrow_length)), (end+((start-end)*arrow_length))),
		  y = c(height1-arrowdiff, height1), lwd = lwd, xpd = T, col = border)
	lines(x = c((end+((start-end)*arrow_length)), (end+((start-end)*arrow_length))),
		  y = c(height2+arrowdiff, height2), lwd = lwd, xpd = T, col = border)
	lines(x = c((end+((start-end)*arrow_length)), end),
		  y = c(height2+arrowdiff, height1+abs(height1-height2)/2),lwd = lwd, xpd = T, col = border)
	lines(x = c((end+((start-end)*arrow_length)), end),
		  y = c(height1-arrowdiff, height1+abs(height1-height2)/2),lwd = lwd, xpd = T, col = border)
	# Add text for the gene name
	if(start<end){
		text(x = ((start+((end-start)/2))+text_adjust), y = height1-((height1-height2)/2), labels = name, xpd = T, cex = text.cex, font = 3, adj = c(0.5, 0.5), col = text.col)
	} else {
		text(x = ((end+((start-end)/2))+text_adjust), y = height1-((height1-height2)/2), labels = name, xpd = T, cex = text.cex, font = 3, adj = c(0.5, 0.5), col = text.col)
	}
}


draw_gene_arrows = function(genes = NULL, ref = NULL, height1 = NULL, height2 = NULL, arrowdiff = NULL, arrow_length = NULL, text_adjust = 0, plot_names = TRUE, fillCOL = NULL, name_replace = NULL, text.cex = 0.5, lwd = 1.5, text.col = "black", draw_line = TRUE){
	if(length(text_adjust)==1) text_adjust = rep(text_adjust, length(genes))
	if(is.null(fillCOL)){
		fillCOL = rep("#BDBDBD", length(genes))
	} else {
		if(length(fillCOL)==1) fillCOL = rep(fillCOL, length(genes))
	}
	if(length(plot_names)==1) plot_names = rep(plot_names, length(genes))
	if(length(text.col)==1) text.col = rep(text.col, length(genes))
	for(i in 1:length(genes)){
		start = ref$start[which(ref$name==genes[i])]
		end = ref$end[which(ref$name==genes[i])]
		rev.compl = ref$strand[which(ref$name==genes[i])]==-1
		if(plot_names[i]) name = genes[i] else name = NA
		if(!is.null(name_replace)){
			if(plot_names[i] & !is.na(name_replace[i])) name = name_replace[i] else name = name
		}
		draw_arrow(start = start, end = end, arrow_length = arrow_length, height1 = height1, height2 = height2, arrowdiff = arrowdiff, fillCOL = fillCOL[i], border = "black", name = name, text_adjust = text_adjust[i], rev.compl = rev.compl, text.cex = text.cex, lwd = lwd, text.col = text.col[i])
	}
	for(i in 1:(length(genes)-1)){
		start = ref$end[which(ref$name==genes[i])]
		end = ref$start[which(ref$name==genes[i+1])]
		if(draw_line) lines(x = c(start, end), y = c(height1-((height1-height2)/2), height1-((height1-height2)/2)), lwd = lwd, xpd = T)
	}

}

# oneLetterCodes = list("Gly" = "G", "Ala" = "A", "Leu" = "L", "Met" = "M", "Phe" = "F", "Trp" = "W", "Lys" = "K", "Gln" = "Q", "Glu" = "E", "Ser" = "S", "Pro" = "P", "Val" = "V", "Ile" = "I", "Cys" = "C", "Tyr" = "Y", "His" = "H", "Arg" = "R", "Asn" = "N", "Asp" = "D", "Thr" = "T", "STO" = "*")

oneLetterCodes = list("Gly" = "G", "Ala" = "A", "Leu" = "L", "Met" = "M", "Phe" = "F", "Trp" = "W", "Lys" = "K", "Gln" = "Q", "Glu" = "E", "Ser" = "S", "Pro" = "P", "Val" = "V", "Ile" = "I", "Cys" = "C", "Tyr" = "Y", "His" = "H", "Arg" = "R", "Asn" = "N", "Asp" = "D", "Thr" = "T", "STO" = "*", "---" = "X")

translate_function_alignment = function(contig = NULL, frame = NULL, oneLetterCodes = NULL){
	contig = unlist(strsplit(contig,""))
	if(frame<=3){
		return(paste(unlist(oneLetterCodes[translate(matrix(totriplet(contig[frame:length(contig)]), nrow = 1))]), collapse = ""))
	} else {
		contig = rev(revcompl_full[contig])
		return(paste(unlist(oneLetterCodes[translate(matrix(totriplet(contig[(frame-3):length(contig)]), nrow = 1))]), collapse = ""))
	}
}

translate_6_frames_alignment = function(contig = NULL, oneLetterCodes = NULL){
	return(sapply(1:6, function(frame, contig) translate_function_alignment(contig, frame, oneLetterCodes), contig = contig, USE.NAMES = F))
}

rev_bases_list = list("A" = "T", "T" = "A", "G" = "C", "C" = "G", "-" = "-", "N" = "N")


read_gene_lookup = function(gene_lookup = NULL, ref_gb = NULL){


	# Read in gene lookup
	gene_lookup = read.table(gene_lookup, h = F, sep = "\t")

	# Read in reference genbank file
	# ref = read_dna_seg_from_file(ref_gb)
	ref = read.table(ref_gb, h = T, sep = "\t", as.is = T, quote = "")
	ref = ref[which(ref$feature=="CDS"),]
	# For each name, if there is more than one entry, label as 'gene_1', 'gene_2'
	ref_unique_names_multiple = table(as.character(ref$name))
	ref_unique_names_multiple = names(ref_unique_names_multiple[which(ref_unique_names_multiple>1)])
	for(i in 1:length(ref_unique_names_multiple)){
		w.i = which(ref$name==ref_unique_names_multiple[i])
		ref$name[w.i] = paste0(ref$name[w.i], "_", 1:length(w.i))
	}
	ref = ref[order(as.numeric(ref$start)),]

	gene_start = ref$start
	gene_end = ref$end
	gene_strand = ref$strand

	# Intergenic ID only assigned if the start of one gene is after the end of the previous gene
	intergenic_names = c()
	for(i in 2:nrow(ref)){
		inter_start = as.numeric(ref$end[i-1])+1
		inter_end = as.numeric(ref$start[i])-1
		if(inter_start<=inter_end){
			intergenic_names = c(intergenic_names, paste(c(ref$name[i-1], ref$name[i]), collapse = ":"))
			gene_start = c(gene_start, inter_start)
			gene_end = c(gene_end, inter_end)
		}
		if(i==nrow(ref)){
			intergenic_names = c(intergenic_names, paste0(ref$name[length(ref$name)], ":"))
			gene_start = c(gene_start, (ref$end[nrow(ref)]+1))
			gene_end = c(gene_end, ref_length)
		}
	}
	# Set the strand for all intergenic regions to 1 - there is no correct strand for these
	gene_strand = c(gene_strand, rep(1, length(intergenic_names)))
	gene_lookup = cbind(gene_lookup, gene_start, gene_end, gene_strand)
	return(gene_lookup)

}
plot_beta_legend = function(){
	par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
	plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
	legend.txt = c("\u03B2 < 0","\u03B2 > 0","","Significance","threshold")
	legend(x = 0.75, y = 0.95, legend.txt, fill = c("#d3d3d3", "#838383", rep("white", 3)), border = NA, bty = "n", cex = 0.7, bg = NA, lty = c(rep(NA, 3), 2, NA), col = "red", pch = NA)
}

plot_frame_legend = function(correct_frame = NULL, frame_cols = NULL){
	par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
	plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
	legend.txt = c(1:6); legend.txt[correct_frame] = paste0(legend.txt[correct_frame], " (Correct)")
	legend.txt = c(legend.txt, "Unaligned", "", "Significance", "threshold","","CDS","rRNA","tRNA","ncRNA","Repeat","Mobile element","Other")
	legend(x = 0.75, y = 0.95, legend.txt, fill = c(frame_cols, "#80808066", rep("white", 4),"#E3E3E3","#15abff","#ff7676","#9cbea6","#ffd370","#F0E442","#f1e7ff"), border = NA, bty = "o", cex = 0.65, bg = NA, lty = c(rep(NA, length(frame_cols)+2), 2, rep(NA,9)), col = "black", pch = NA)
}

plot_axes = function(translation){
	axis(1, at = pretty(0:nchar(translation)))
	axis(2)
	box()
}


plot_allframes_manhattan = function(gene_i_results_list = NULL, prefix = NULL, genes_name = NULL, length_correct = NULL, all_translations = NULL, correct_frame = NULL, which_kmers_no_result = NULL, bonferroni = NULL, ylim_max = NULL, maname = NULL, malim = NULL, kmer_len = NULL, x.adjust = 0, gene_end_line = NULL, start = NULL, end = NULL, macormaf = NULL, output_dir = NULL, ref_gb_full = NULL){
	
	
	whcol_sstart = which(colnames(gene_i_results_list[[1]])=="sstart")
	whcol_send = which(colnames(gene_i_results_list[[1]])=="send")
	if(is.null(kmer_len)) kmer_len = max(apply(gene_i_results_list[[1]][, whcol_sstart:whcol_send], 1, function(x) length(as.numeric(x[1]):as.numeric(x[2])))[which(as.numeric(gene_i_results_list[[1]][["gapopen"]])==0)])

	ymax = max(as.numeric(unlist(sapply(gene_i_results_list, function(x) return(x$negLog10), USE.NAMES = F))))
	if(!is.null(which_kmers_no_result)) ymax = max(c(ymax, as.numeric(which_kmers_no_result$negLog10)))

	if(is.null(maname)) maname = "_allkmers"

	if(!is.null(malim)){
		ymax = max(as.numeric(unlist(sapply(gene_i_results_list, function(x, malim, macormaf) return(x$negLog10[which(as.numeric(x[[macormaf]])>=malim)]), malim = malim, macormaf = macormaf, USE.NAMES = F))))
		if(!is.null(which_kmers_no_result)) ymax = max(c(ymax, as.numeric(which_kmers_no_result$negLog10)[which(as.numeric(which_kmers_no_result[[macormaf]])>=malim)]))
	}

	if(is.null(ylim_max)){
		ylim = c(0, ymax)
		plotname = paste0("_allframes_Manhattan", maname,".png")
	} else {
		ylim = c(0, ylim_max)
		plotname = paste0("_allframes_Manhattan_ylim", ylim_max, maname, ".png")
	}


	xpos_list = list()
	ypos_list = list()

	for(f in 1:6){

		xpos = cbind(as.numeric(as.character(gene_i_results_list[[f]]$sstart)), as.numeric(as.character(gene_i_results_list[[f]]$send)))
		# ypos = c(as.numeric(as.character(gene_i_results_list[[f]]$negLog10)), as.numeric(which_kmers_no_result[,2]))
		ypos = c(as.numeric(as.character(gene_i_results_list[[f]]$negLog10)))

		beta = c(as.numeric(gene_i_results_list[[f]]$beta))
		beta_col = rep("#d3d3d3", length(beta)); beta_col[which(beta>0)] = "#838383"


		# Adjust xpos for frame
		if(f<=3){
			seq_pos1 = seq(from = start+(f-1), by = 3, to = end)
			seq_pos2 = seq(from = (start+2+f-1), by = 3, to = end)
		} else {
			seq_pos1 = seq(from = end-(length(4:f)-1), by = -3, to = start)
			seq_pos2 = seq(from = end-(length(4:f)-1)-2, by = -3, to = start)
		}

		xpos_1_new_aligned = sapply(xpos[,1], function(x, s) s[x], s = seq_pos1, USE.NAMES = F)
		xpos_2_new_aligned = sapply(xpos[,2], function(x, s) s[x], s = seq_pos2, USE.NAMES = F)

		xpos[,1] = xpos_1_new_aligned
		xpos[,2] = xpos_2_new_aligned

		# # Flip positions for those that have been reversed
		# if((correct_frame==1 & f>=4) | (correct_frame==4 & f<=3)){
			# xpos[,1] = length_correct-xpos[,1]+1
			# xpos[,2] = length_correct-xpos[,2]+1
		# }

		# # Adjust for times where the x axis includes upstream and downstream regions
		# xpos[,1] = as.numeric(xpos[,1])-x.adjust
		# xpos[,2] = as.numeric(xpos[,2])-x.adjust


		if(length(ypos)!=nrow(xpos)) stop("Error in Multi manhattan plot","\n")

		xpos_list[[f]] = xpos
		ypos_list[[f]] = ypos


	}

	max_xpos = max(unlist(xpos_list))
	min_xpos = min(unlist(xpos_list))

	png(paste0(prefix, "_", genes_name, plotname), width = 20, height = 15, units = "cm", res = 600)
	par(mar = c(5.1, 4.1, 2, 6))
	plot(range(c(start, end+80)), c(0, ymax), type = "n", xlab = "Amino acid position in reference", ylab = expression(paste("Significance (-log"[10],italic(' p'),") LMM",collapse="")), main = paste0(genes_name," (Protein kmers)"), axes = F, ylim = ylim)
	# plot_axes(all_translations[correct_frame])
	frame_cols_incorrect = c("#009E73", "#0072B2", "#E69F00", "#800080", "#ADD8E6")
	frame_cols = rep("#FF0000", 6); frame_cols[-correct_frame] = frame_cols_incorrect
	frame_cols = paste0(frame_cols, "66")

	if(!is.null(gene_end_line)){
		abline(v = gene_end_line, col = "#80808066")
		abline(v = 1, col = "#80808066")
	}

	for(f in 1:6){

		if(is.null(malim)) which_to_plot = 1:length(ypos_list[[f]]) else which_to_plot = which(c(as.numeric(gene_i_results_list[[f]][[macormaf]]))>=(malim))

		for(k in which_to_plot){
			lines(x = c(as.numeric(xpos_list[[f]][k,1]), as.numeric(xpos_list[[f]][k,2])), y = rep(as.numeric(ypos_list[[f]][k]), 2), col = frame_cols[f])
		}
	}

	# Plot unmapped kmers
	if(!is.null(which_kmers_no_result)){
		xpos_unmapped = sample(seq(from = end+50, to = end+100, length.out = nrow(which_kmers_no_result)), nrow(which_kmers_no_result), replace = F)
		xpos_unmapped = cbind(xpos_unmapped, xpos_unmapped+(kmer_len-1))
		ypos_unmapped = as.numeric(which_kmers_no_result$negLog10)
		if(is.null(malim)) which_to_plot = 1:length(ypos_unmapped) else which_to_plot = which(as.numeric(which_kmers_no_result[[macormaf]])>=(malim))
		if(length(which_to_plot)>0){
			for(k in which_to_plot){
				lines(x = c(as.numeric(xpos_unmapped[k,1]), as.numeric(xpos_unmapped[k,2])), y = rep(as.numeric(ypos_unmapped[k]), 2), col = "#80808066")
			}
		}
	}

	abline(h = bonferroni, col = "black", lty = 2)


	plot_gene_arrows_manhattan(ref_gb_full = ref_gb_full, start = start, end = end, xpos = xpos)

	plot_frame_legend(correct_frame = correct_frame, frame_cols = frame_cols)
	dev.off()

	return(ymax)

}

plot_singleframe_manhattan_protein = function(prefix = NULL, genes_name = NULL, correct_or_wrong = NULL, j = NULL, xpos = NULL, ypos = NULL, all_translations = NULL, beta_col = NULL, bonferroni = NULL, ylim_max = NULL, maname = NULL, x.adjust = NULL, ref_gb_full = NULL, start = NULL, end = NULL, ref_length = NULL, kmer_length = NULL){

	if(is.null(maname)) maname = "_allkmers"


	# Adjust xpos for frame
	if(j<=3){
		seq_pos1 = seq(from = start+(j-1), by = 3, to = end)
		seq_pos2 = seq(from = (start+2+j-1), by = 3, to = end)
	} else {
		seq_pos1 = seq(from = end-(length(4:j)-1), by = -3, to = start)
		seq_pos2 = seq(from = end-(length(4:j)-1)-2, by = -3, to = start)
	}


	which_aligned = which(beta_col=="#d3d3d3" | beta_col=="#838383")
	which_unaligned = which(beta_col!="#d3d3d3" & beta_col!="#838383")
	xpos_1_new_aligned = sapply(xpos[which_aligned,1], function(x, s) s[x], s = seq_pos1, USE.NAMES = F)
	xpos_2_new_aligned = sapply(xpos[which_aligned,2], function(x, s) s[x], s = seq_pos2, USE.NAMES = F)
	xpos_1_new_unaligned = sample(seq(from = end+50, to = end+100, length.out = length(which_unaligned)), length(which_unaligned), replace = F)
	xpos_2_new_unaligned = xpos_1_new_unaligned+((kmer_length*3)-1)

	xpos_1_new = rep(0, nrow(xpos)); xpos_1_new[which_aligned] = xpos_1_new_aligned; xpos_1_new[which_unaligned] = xpos_1_new_unaligned
	xpos_2_new = rep(0, nrow(xpos)); xpos_2_new[which_aligned] = xpos_2_new_aligned; xpos_2_new[which_unaligned] = xpos_2_new_unaligned

	xpos[,1] = xpos_1_new
	xpos[,2] = xpos_2_new

	if(is.null(ylim_max)){
		ylim = c(0, max(as.numeric(ypos)))
		plotname = paste0("_Manhattan", maname,".png")
	} else {
		ylim = c(0, ylim_max)
		plotname = paste0("_Manhattan_ylim", ylim_max, maname, ".png")
	}

	png(paste0(prefix, "_", genes_name, "_", correct_or_wrong, "_", j, plotname), width = 20, height = 15, units = "cm", res = 600)
	par(mar = c(5.1, 4.1, 2, 6))
	plot(range(as.numeric(as.vector(xpos))), c(0, max(as.numeric(ypos))), type = "n", xlab = "Amino acid position in reference", ylab = expression(paste("Significance (-log"[10],italic(' p'),") LMM",collapse="")), main = paste0(genes_name, " (Protein kmers)"), axes = F, ylim = ylim)
	# plot_axes(all_translations[j])
	for(k in 1:nrow(xpos)){
		lines(x = c(as.numeric(xpos[k,1]), as.numeric(xpos[k,2])), y = rep(as.numeric(ypos[k]), 2), col = beta_col[k])
	}
	abline(h = bonferroni, col = "black", lty = 2)

	plot_gene_arrows_manhattan(ref_gb_full = ref_gb_full, start = start, end = end, xpos = xpos)


	plot_beta_and_feature_legend()
	dev.off()

}


plot_singleframe_manhattan_nucleotide = function(prefix = NULL, genes_name = NULL, xpos = NULL, ypos = NULL, ref_fa = NULL, beta_col = NULL, bonferroni = NULL, ylim_max = NULL, maname = NULL, gene_end_line = NULL, x.adjust = 0, ref_gb_full = NULL, start = NULL, end = NULL){

	if(is.null(maname)) maname = "_allkmers"

	which_gene_arrows = which(ref_gb_full$end>=(start-800) & ref_gb_full$start<=(end+800))
	# Just keep genes in the region
	ref_subset = ref_gb_full[which_gene_arrows,]
	ref_subset$name = as.character(ref_subset$name)
	ref_subset$name[which(is.na(ref_subset$name))] = paste0("NA",1:length(which(is.na(ref_subset$name))))
	# Rename duplicate names to avoid complication later
	ref_unique_names_multiple = table(as.character(ref_subset$name))
	ref_unique_names_multiple = names(ref_unique_names_multiple[which(ref_unique_names_multiple>1)])
	ref_unique_names_multiple = ref_unique_names_multiple[ref_unique_names_multiple!=""]
	if(length(ref_unique_names_multiple)>0){
		for(i in 1:length(ref_unique_names_multiple)){
			w.i = which(ref_subset$name==ref_unique_names_multiple[i] & ref_subset$feature!="CDS")
			ref_subset$name[w.i] = paste0(ref_subset$name[w.i], "_", 2:(length(w.i)+1))
		}
	}

	if(any(ref_subset[,1]==genes_name)){
		if(ref_subset$strand[which(ref_subset$name==genes_name & ref_subset$feature=="CDS")]==1){
			x.adjust = -(start-1)
			ref_subset = ref_subset[order(as.numeric(ref_subset$start)),]
			wh_genematch = which(ref_subset[,1]==genes_name)
			genematch_start = as.numeric(ref_subset$start[wh_genematch])
		} else {
			length_refregion = length(start:end)
			new_xpos1 = sapply(as.numeric(xpos[,1]), function(x,adjust) length(adjust:x), adjust = length_refregion, USE.NAMES = F)
			new_xpos2 = sapply(as.numeric(xpos[,2]), function(x,adjust) length(adjust:x), adjust = length_refregion, USE.NAMES = F)
			xpos[which(beta_col=="#d3d3d3" | beta_col=="#838383"),1] = new_xpos1[which(beta_col=="#d3d3d3" | beta_col=="#838383")]
			xpos[which(beta_col=="#d3d3d3" | beta_col=="#838383"),2] = new_xpos2[which(beta_col=="#d3d3d3" | beta_col=="#838383")]
			x.adjust = -(start-1)
			ref_subset = ref_subset[order(as.numeric(ref_subset$start)),]
			wh_genematch = which(ref_subset[,1]==genes_name & ref_subset$feature=="CDS")
		}
	} else {
		x.adjust = -(start-1)
		first_gene = unlist(strsplit(genes_name,":"))[1]
		wh_genematch = which(ref_subset[,1]==first_gene & ref_subset$feature=="CDS")
	}


	# Adjust for times where the x axis includes upstream and downstream regions
	xpos[,1] = as.numeric(xpos[,1])-x.adjust
	xpos[,2] = as.numeric(xpos[,2])-x.adjust


	if(is.null(ylim_max)){
		ylim = c(0, max(as.numeric(ypos)))
		plotname = paste0("_Manhattan", maname,".png")
	} else {
		ylim = c(0, ylim_max)
		plotname = paste0("_Manhattan_ylim", ylim_max, maname, ".png")
	}

	png(paste0(prefix, "_", genes_name, plotname), width = 20, height = 15, units = "cm", res = 600)
	par(mar = c(5.1, 4.1, 2, 6))
	plot(range(as.numeric(as.vector(xpos))), c(0, max(as.numeric(ypos))), type = "n", xlab = "Position in reference", ylab = expression(paste("Significance (-log"[10],italic(' p'),") LMM",collapse="")), main = paste0(genes_name, " (Nucleotide kmers)"), axes = F, ylim = ylim)

	if(!is.null(gene_end_line)){
		abline(v = gene_end_line, col = "#80808066")
		abline(v = 1, col = "#80808066")
	}


	for(k in 1:nrow(xpos)){
		lines(x = c(as.numeric(xpos[k,1]), as.numeric(xpos[k,2])), y = rep(as.numeric(ypos[k]), 2), col = beta_col[k])
	}
	abline(h = bonferroni, col = "black", lty = 2)


	plot_height = par("usr")[4]-par("usr")[3]
	height1 = -(plot_height*0.05)
	height2 = -(plot_height*0.02)
	arrowdiff = plot_height*0.01
	axis_height = mean(c(height1, height2))

	axis(1, at = pretty(range(xpos)), pos = axis_height, tck = -0.03)

	gene_cols = rep("#E3E3E3",nrow(ref_subset))
	gene_cols[which(ref_subset$feature=="tRNA")] = "#ff7676"
	gene_cols[which(ref_subset$feature=="rRNA")] = "#15abff"
	gene_cols[which(ref_subset$feature=="repeat_region" | ref_subset$feature=="repeat_region_pseudo")] = "#ffd370"
	gene_cols[which(ref_subset$feature=="ncRNA")] = "#9cbea6"
	gene_cols[which(ref_subset$feature=="mobile_element")] = "#F0E442"
	gene_cols[which(ref_subset$feature=="misc_feature" | ref_subset$feature=="misc_feature_intron_pseudo" | ref_subset$feature=="misc_feature_pseudo")] = "#f1e7ff"
	text.col = rep("black", nrow(ref_subset)) # ; text.col[which(gene_cols!="#E3E3E3")] = "white"

	which_na = which(substr(as.character(ref_subset$name),1,2)=="NA")
	if(length(which_na)>0){
		name_replace = as.character(ref_subset$name)
		name_replace[which_na] = ""
	} else {
		name_replace = NULL
	}

	draw_gene_arrows(genes = as.character(ref_subset$name), ref = ref_subset, arrow_length = 0.15, height1 = height1, height2 = height2, arrowdiff = arrowdiff, fillCOL = gene_cols, plot_names = TRUE, text.cex = 0.45, text.col = text.col, lwd = 0.8, draw_line = FALSE, name_replace = name_replace)
	xleft = par("usr")[1]
	xright = par("usr")[2]
	ytop = par("usr")[3]
	rect(xright = xleft, xleft = xleft-10000, ytop = 0, ybottom = (height1*3), col = "white", border = NA, xpd = TRUE)
	rect(xright = xright+10000, xleft = xright, ytop = 0, ybottom = (height1*3), col = "white", border = NA, xpd = TRUE)

	axis(2)

	plot_beta_and_feature_legend()
	dev.off()

}


plot_gene_arrows_manhattan = function(ref_gb_full = NULL, start = NULL, end = NULL, xpos = NULL){

	which_gene_arrows = which(ref_gb_full$end>=(start-999) & ref_gb_full$start<=(end+999))
	# Just keep genes in the region
	ref_subset = ref_gb_full[which_gene_arrows,]
	ref_subset$name = as.character(ref_subset$name)
	ref_subset$name[which(is.na(ref_subset$name))] = paste0("NA",1:length(which(is.na(ref_subset$name))))
	# Rename duplicate names to avoid complication later
	ref_unique_names_multiple = table(as.character(ref_subset$name))
	ref_unique_names_multiple = names(ref_unique_names_multiple[which(ref_unique_names_multiple>1)])
	ref_unique_names_multiple = ref_unique_names_multiple[ref_unique_names_multiple!=""]
	if(length(ref_unique_names_multiple)>0){
		for(i in 1:length(ref_unique_names_multiple)){
			w.i = which(ref_subset$name==ref_unique_names_multiple[i] & ref_subset$feature!="CDS")
			ref_subset$name[w.i] = paste0(ref_subset$name[w.i], "_", 2:(length(w.i)+1))
		}
	}

	plot_height = par("usr")[4]-par("usr")[3]
	height1 = -(plot_height*0.05)
	height2 = -(plot_height*0.02)
	arrowdiff = plot_height*0.01
	axis_height = mean(c(height1, height2))

	axis(1, at = pretty(range(xpos)), pos = axis_height, tck = -0.03)

	gene_cols = rep("#E3E3E3",nrow(ref_subset))
	gene_cols[which(ref_subset$feature=="tRNA")] = "#ff7676"
	gene_cols[which(ref_subset$feature=="rRNA")] = "#15abff"
	gene_cols[which(ref_subset$feature=="repeat_region" | ref_subset$feature=="repeat_region_pseudo")] = "#ffd370"
	gene_cols[which(ref_subset$feature=="ncRNA")] = "#9cbea6"
	gene_cols[which(ref_subset$feature=="mobile_element")] = "#F0E442"
	gene_cols[which(ref_subset$feature=="misc_feature" | ref_subset$feature=="misc_feature_intron_pseudo" | ref_subset$feature=="misc_feature_pseudo")] = "#f1e7ff"
	text.col = rep("black", nrow(ref_subset))

	which_na = which(substr(as.character(ref_subset$name),1,2)=="NA")
	if(length(which_na)>0){
		name_replace = as.character(ref_subset$name)
		name_replace[which_na] = ""
	} else {
		name_replace = NULL
	}

	draw_gene_arrows(genes = as.character(ref_subset$name), ref = ref_subset, arrow_length = 0.15, height1 = height1, height2 = height2, arrowdiff = arrowdiff, fillCOL = gene_cols, plot_names = TRUE, text.cex = 0.45, text.col = text.col, lwd = 0.8, draw_line = FALSE, name_replace = name_replace)
	xleft = par("usr")[1]
	xright = par("usr")[2]
	ytop = par("usr")[3]
	rect(xright = xleft, xleft = xleft-10000, ytop = 0, ybottom = (height1*3), col = "white", border = NA, xpd = TRUE)
	rect(xright = xright+10000, xleft = xright, ytop = 0, ybottom = (height1*3), col = "white", border = NA, xpd = TRUE)

	axis(2)


}


plot_beta_and_feature_legend = function(){
	par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
	plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
	legend.txt = c("\u03B2 < 0","\u03B2 > 0","Unaligned","","Significance","threshold","","CDS","rRNA","tRNA","ncRNA","Repeat","Mobile element","Other")
	legend(x = 0.75, y = 0.95, legend.txt, fill = c("#d3d3d3", "#838383", "#D55E00", rep("white", 4), "#E3E3E3","#15abff","#ff7676","#9cbea6","#ffd370","#F0E442","#f1e7ff"), border = NA, bty = "o", cex = 0.65, bg = NA, lty = c(rep(NA, 4), 2, rep(NA,9)), col = "black", pch = NA)
}







fix_res_table_protein = function(res = NULL, refseq = NULL){

	res$kmer = as.character(res$kmer)
	res$sseq = as.character(res$sseq)
	res$qseq = as.character(res$qseq)

	# Fix start and end positions so the whole kmer is aligned, not just plotting the BLAST result
	which_start_short = which(as.numeric(res$qstart)>1)
	new_start = as.numeric(res$sstart)[which_start_short]
	res$sstart[which_start_short] = sapply(1:length(new_start), function(p, s, q) return(c(s[p]-q[p]+1)), s = new_start, q = as.numeric(res$qstart)[which_start_short],USE.NAMES = F)
	# Now check if the start is less than 1 for any
	any_minus_start = which(as.numeric(res$sstart)<1)
	if(length(any_minus_start)>0){
		minus_start_seq = as.character(res$kmer)[any_minus_start]
		res$sstart[any_minus_start] = 1
		minus_start_seq = sapply(1:length(minus_start_seq), function(p, seq, pos) return(substr(seq[p], (length(pos[p]:1)), nchar(seq[p]))), seq = minus_start_seq, pos = as.numeric(res$sstart)[any_minus_start], USE.NAMES = F)
		# res$kmer[any_minus_start] = as.character(minus_start_seq)
	}
	# Update ref seq
	new_which_start_short = which(as.numeric(res$sstart[which_start_short])!=new_start)
	if(length(new_which_start_short)>0){
		new_start = new_start[new_which_start_short]
		which_start_short = which_start_short[new_which_start_short]
		res$qseq[which_start_short] = sapply(1:length(new_start), function(p, o, s, q, kmer, seq) paste0(substr(kmer[p], q[p]-length(s[p]:(o[p]-1)), q[p]-1), seq[p]),o = new_start, s = as.numeric(res$sstart[which_start_short]), q = as.numeric(res$qstart)[which_start_short], kmer = as.character(res$kmer)[which_start_short], seq = as.character(res$qseq)[which_start_short], USE.NAMES = F)
		res$sseq[which_start_short] = sapply(1:length(new_start), function(p, s, q, refseq, seq) paste0(substr(refseq, min((s[p]-q[p]+1):(s[p]-1)), max((s[p]-q[p]+1):(s[p]-1))), seq[p]),s = new_start, q = as.numeric(res$qstart)[which_start_short], refseq = refseq, seq = as.character(res$sseq)[which_start_short], USE.NAMES = F)
	}

	which_end_short = which(as.numeric(res$qend)<11)
	new_end = as.numeric(res$send)[which_end_short]
	res$send[which_end_short] = sapply(1:length(new_end), function(p, s, q) return(c(s[p]-q[p]+11)), s = new_end, q = as.numeric(res$qend)[which_end_short],USE.NAMES = F)
	# res$sseq = substring(refseq, as.numeric(res$sstart), as.numeric(res$send))
	# Now check if the end is after the end of the protein
	any_long_end = which(as.numeric(res$send)>nchar(refseq))
	if(length(any_long_end)>0){
		# long_end_seq = as.character(res$kmer)[any_long_end]
		res$send[any_long_end] = nchar(refseq)
		# long_end_seq = sapply(1:length(long_end_seq), function(p, seq, pos, len) return(substr(seq[p], 1, length(pos[p]:len))), seq = long_end_seq, pos = as.numeric(res$sstart)[any_long_end], len = nchar(refseq), USE.NAMES = F)
		# res$kmer[any_long_end] = long_end_seq
	}
	# Update ref seq
	new_which_end_short = which(as.numeric(res$send[which_end_short])!=new_end)
	if(length(new_which_end_short)>0){
		new_end = new_end[new_which_end_short]
		which_end_short = which_end_short[new_which_end_short]
		res$sseq[which_end_short] = sapply(1:length(new_end), function(p, s, q, refseq, seq) paste0(seq[p], substr(refseq, (s[p]+1), (c(s[p]-q[p]+11)))),s = new_end, q = as.numeric(res$qend)[which_end_short], refseq = refseq, seq = as.character(res$sseq)[which_end_short], USE.NAMES = F)
		res$qseq[which_end_short] = sapply(1:length(new_end), function(p, o, s, q, kmer, seq) paste0(seq[p], substr(kmer[p], (q[p]+1), (q[p]+length((o[p]+1):s[p])))), o = new_end, s = as.numeric(res$send)[which_end_short], q = as.numeric(res$qend)[which_end_short], kmer = as.character(res$kmer)[which_end_short], seq = as.character(res$qseq)[which_end_short], USE.NAMES = F)
	}
	res$origkmer = res$kmer
	return(res)
}

fix_res_table_nucleotide = function(res = NULL, ref_fa = NULL, kmer_len = NULL){

	res$kmer = as.character(res$kmer)
	res$sseq = as.character(res$sseq)
	res$qseq = as.character(res$qseq)
	res$origkmer = as.character(res$kmer)
	# For those which aligned to the reverse strand, reverse complement the sequences
	wh_revstrand = which(res$sstrand=="minus")
	new_sseq = sapply(as.character(res$sseq[wh_revstrand]), function(x) paste(rc_full(unlist(strsplit(x,""))),collapse = ""), USE.NAMES = F)
	new_qseq = sapply(as.character(res$qseq[wh_revstrand]), function(x) paste(rc_full(unlist(strsplit(x,""))),collapse = ""), USE.NAMES = F)
	new_sstart = as.numeric(res$send)[wh_revstrand]
	new_send = as.numeric(res$sstart)[wh_revstrand]
	new_kmer = sapply(as.character(res$kmer[wh_revstrand]), function(x) paste(rc_full(unlist(strsplit(x,""))),collapse = ""), USE.NAMES = F)
	new_qend = kmer_len-as.numeric(res$qstart)[wh_revstrand]+1
	new_qstart = 31-as.numeric(res$qend)[wh_revstrand]+1

	res$sstart[wh_revstrand] = new_sstart
	res$send[wh_revstrand] = new_send
	res$sseq[wh_revstrand] = new_sseq
	res$qseq[wh_revstrand] = new_qseq
	res$kmer[wh_revstrand] = new_kmer
	res$qstart[wh_revstrand] = new_qstart
	res$qend[wh_revstrand] = new_qend

	# Fix start and end positions so the whole kmer is aligned, not just plotting the BLAST result
	which_start_short = which(as.numeric(res$qstart)>1)
	new_start = as.numeric(res$sstart)[which_start_short]
	res$sstart[which_start_short] = sapply(1:length(new_start), function(p, s, q) return(c(s[p]-q[p]+1)), s = new_start, q = as.numeric(res$qstart)[which_start_short], USE.NAMES = F)
	# Now check if the start is less than 1 for any
	any_minus_start = which(as.numeric(res$sstart)<1)
	if(length(any_minus_start)>0){
		minus_start_seq = as.character(res$kmer)[any_minus_start]
		res$sstart[any_minus_start] = 1
		minus_start_seq = sapply(1:length(minus_start_seq), function(p, seq, pos) return(substr(seq[p], (length(pos[p]:1)), nchar(seq[p]))), seq = minus_start_seq, pos = as.numeric(res$sstart)[any_minus_start], USE.NAMES = F)
		# res$kmer[any_minus_start] = as.character(minus_start_seq)
	}

	# Update ref seq
	new_which_start_short = which(as.numeric(res$sstart[which_start_short])!=new_start)
	if(length(new_which_start_short)>0){
		new_start = new_start[new_which_start_short]
		which_start_short = which_start_short[new_which_start_short]
		res$qseq[which_start_short] = sapply(1:length(new_start), function(p, o, s, q, kmer, seq) paste0(substr(kmer[p], q[p]-length(s[p]:(o[p]-1)), q[p]-1), seq[p]),o = new_start, s = as.numeric(res$sstart[which_start_short]), q = as.numeric(res$qstart)[which_start_short], kmer = as.character(res$kmer)[which_start_short], seq = as.character(res$qseq)[which_start_short], USE.NAMES = F)
		res$sseq[which_start_short] = sapply(1:length(new_start), function(p, s, q, ref_fa, seq) paste0(substr(ref_fa, min((s[p]-q[p]+1):(s[p]-1)), max((s[p]-q[p]+1):(s[p]-1))), seq[p]),s = new_start, q = as.numeric(res$qstart)[which_start_short], ref_fa = ref_fa, seq = as.character(res$sseq)[which_start_short], USE.NAMES = F)
	}

	which_end_short = which(as.numeric(res$qend)<kmer_len)
	new_end = as.numeric(res$send)[which_end_short]
	res$send[which_end_short] = sapply(1:length(new_end), function(p, s, q) return(c(s[p]-q[p]+kmer_len)), s = new_end, q = as.numeric(res$qend)[which_end_short],USE.NAMES = F)
	# Now check if the end is after the end of the protein
	any_long_end = which(as.numeric(res$send)>nchar(ref_fa))
	if(length(any_long_end)>0){
		res$send[any_long_end] = nchar(ref_fa)
	}
	# Update ref seq
	new_which_end_short = which(as.numeric(res$send[which_end_short])!=new_end)
	if(length(new_which_end_short)>0){
		new_end = new_end[new_which_end_short]
		which_end_short = which_end_short[new_which_end_short]
		res$sseq[which_end_short] = sapply(1:length(new_end), function(p, s, q, ref_fa, seq) paste0(seq[p], substr(ref_fa, (s[p]+1), (c(s[p]-q[p]+kmer_len)))),s = new_end, q = as.numeric(res$qend)[which_end_short], ref_fa = ref_fa, seq = as.character(res$sseq)[which_end_short], USE.NAMES = F)
		res$qseq[which_end_short] = sapply(1:length(new_end), function(p, o, s, q, kmer, seq) paste0(seq[p], substr(kmer[p], (q[p]+1), (q[p]+length((o[p]+1):s[p])))), o = new_end, s = as.numeric(res$send)[which_end_short], q = as.numeric(res$qend)[which_end_short], kmer = as.character(res$kmer)[which_end_short], seq = as.character(res$qseq)[which_end_short], USE.NAMES = F)
	}


	return(res)
}


plot_alignment_function_protein = function(genestart = NULL, geneend = NULL, minor_allele_threshold = NULL, res = NULL, nsamples = NULL, maname = NULL, bonferroni = NULL, translation = NULL, prefix = NULL, gene_name = NULL, correct_or_wrong = NULL, j = NULL, p = NULL, lowfreq = NULL, plot_ref = NULL, main = NULL, x.adjust = 0, nplots = NULL, override_signif = FALSE, out_table = NULL, plot_figures = TRUE, col_lib = NULL, macormaf = NULL, output_dir = NULL, kmer_type = NULL, kmer_length = NULL, ref.name = NULL){

	legend.txt = c(names(col_lib_pro),"Invariant","","\u03B2 < 0","\u03B2 > 0","","Significance","threshold")
	legend.txt[which(legend.txt=="-")] = "Gap"
	# legend.fill = c(unlist(col_lib_pro),"#d3d3d3",NA, "#d3d3d3", "#838383", NA, NA, NA, NA, NA)[which_legend]
	# legend.col = c(rep(NA,length(col_lib_pro)+5),"black",NA, "black", NA)[which_legend]
	legend.col = c(unlist(col_lib_pro),"#d3d3d3",NA, "#d3d3d3", "#838383", NA, "black", NA)
	legend.pch = c(rep(15,length(col_lib_pro)+5), NA,NA)
	legend.lty = c(rep(NA, length(col_lib_pro)+5), 2,NA)
	legend.pch.cex = rep(1, length(legend.txt))
	sep_lines = c(0.5)
	ref_num = 1

	# Which of the BLAST results should be plotted (which are in the region and above minor allele threshold)
	which_to_align = which(as.numeric(res$send)>= genestart & as.numeric(res$sstart)<=geneend & as.numeric(res[[macormaf]])>=(minor_allele_threshold))
	# If wanting to show which are low frequency (if they are in the plot)
	# then find which are below the cutoff and feed them in to be coloured translucent
	if(!is.null(lowfreq)){
		which_low_freq = which(c(as.numeric(res$mac)<(nsamples*lowfreq))[which_to_align])
	} else {
		which_low_freq = NULL
	}
	if(length(which_to_align)>0){
		# How many kmers are above the bonferroni threshold
		bonferroni_lim = length(which(as.numeric(res$negLog10)[which_to_align]<bonferroni))
		if(bonferroni_lim!=length(which_to_align) | override_signif){
			# Pull out reference amino acid sequence for the region
			gene_db_i = matrix(unlist(strsplit(translation, ""))[genestart:geneend], nrow = 1)
			# gene_db_i = rbind(gene_db_i, gene_db_i, gene_db_i, gene_db_i, gene_db_i, gene_db_i)

			snp_num = round((length(which_to_align)+ref_num)*0.05)

			beta_to_input = cbind(rep(0, nrow(res)), as.numeric(res$beta), rep(0, nrow(res)))[which_to_align,]
			if(is.null(dim(beta_to_input))) beta_to_input = matrix(beta_to_input, nrow = 1)

			if(plot_figures){

				aligncols = run_plot_alignment(kmerseq = as.character(res[["qseq"]])[which_to_align], snp_cols = NULL,
									   kmerpos = as.numeric(res$sstart)[which_to_align],
									   prefix = paste0(output_dir, prefix, "_", kmer_type, kmer_length, "_", ref.name, "_", gene_name, "_", correct_or_wrong, "_", j,
									   "_plot_",p,"_aminoacids_",(genestart-x.adjust),"_to_",(geneend-x.adjust), maname),
									   ref_num = ref_num, snp_num = snp_num, labs.cex = 0.8, ref = gene_db_i,
								   	   plot_subset = NULL, snp_type = NULL, legend.txt = NULL, legend.fill = legend.fill,
								   	   sep_lines = sep_lines, genestart = genestart, rev_compl=FALSE,
								   	   bonferroni_threshold = bonferroni_lim, plot.dim = c(21,18), legend.cex = 0.7,
								       legend.xpos = c(0.868, 0.9), legend.ypos = c(0.703,0.717),
								       text.xpos = c(0.794,0.830, 0.869), text.ypos = c(0.42), axis.cex = 0.7,
								       legend.lty = legend.lty, legend.pch = legend.pch, legend.col= legend.col,
								       line.sep.lwd = 0.15, reverse.xaxis = FALSE, manhattan.order = TRUE,
								       reverse.xaxis.start = NULL,
								       beta_estimate = beta_to_input,
								       se.d.kmers = 0, beta_cols_list = c("#d3d3d3", "#009E73", "#E69F00", "#838383"),
								       ref.name = "", legend.odds = FALSE, close.file = FALSE, margins = c(2,3,2.5,6),
								       which_translucent = NULL, plot_ref = plot_ref, forward.xaxis.start = genestart-x.adjust,
								       col_lib = col_lib)

				# Plot reference colours along bottom
				refcols_ytop = ((nrow(aligncols))*0.04)+0.5
				for(k in 1:ncol(aligncols)){
					rect(xleft = k-0.5, xright = k+0.5, ybottom = 0.5, ytop = refcols_ytop, col = aligncols[1,k], border = NA, xpd = TRUE)
				}
				# Add lines and reference name
				sapply(c(1:ncol(aligncols)), function(k, ytop, lwd) lines(x = c(k-0.5,k-0.5), y = c(0.5, ytop), lwd = lwd, col = "white"), lwd = 0.15, ytop = refcols_ytop, USE.NAMES = F)
				text(x = 1,y = mean(c(refcols_ytop, 0.5)),"REF",xpd=TRUE, pos = 2, cex = 0.8, adj = 1)

				# Annotate the ref genome
				# text(x = 1,y = length(which_to_align)+ref_num+snp_num+10,"REF",xpd=TRUE, pos = 2, cex = 0.8, adj = 1)
				# Plot axes
				axis.start = which(c((genestart-x.adjust):((genestart-x.adjust)+ncol(gene_db_i)-1))%%10==0)[1]
				if(is.na(axis.start)){
					axis(3, at = 1:ncol(gene_db_i), labels = c((genestart-x.adjust):((genestart-x.adjust)+ncol(gene_db_i)-1)), cex.axis = 0.7, lwd.ticks = NA, line = 0.4, lwd = NA, xpd = T)
					axis(3, at = 1:ncol(gene_db_i), lwd.ticks = NA, label = NA, lwd = 0.7, line = 0.81, xpd = T)
					axis(3, at = 1:ncol(gene_db_i), labels = NA, cex.axis = 0.7, lwd = 0.7, line = 0.81, xpd = T)
				} else {
					axis(3, at = seq(from = axis.start, by = 10, to = ncol(gene_db_i)), labels = c(seq(from = c((genestart-x.adjust) +axis.start-1),
						 by = 10, to = c((genestart-x.adjust) +ncol(gene_db_i)-1))), cex.axis = 0.7, lwd.ticks = NA, line = 0.4, lwd = NA, xpd = T)
					axis(3, at = c(1-0.5, ncol(gene_db_i)+0.5), lwd.ticks = NA, label = NA, lwd = 0.7, line = 0.81, xpd = T)
					axis(3, at = seq(from = axis.start, by = 10, to = ncol(gene_db_i)), labels = NA, cex.axis = 0.7, lwd = 0.7, line = 0.81, xpd = T)
				}
				axis.label.pos = axis(3, at = 1:ncol(gene_db_i), labels = rep("", ncol(gene_db_i)), cex.axis = 0.4, lwd = NA, line = -0.75, xpd = T)
				for(k in axis.label.pos) axis(3, at = k, labels = as.character(gene_db_i[1,k]), cex.axis = 0.4, lwd = NA, line = -0.75, xpd = T)

				# Plot legend
				par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
				plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
				# if(!is.null(main)) text(x = -1.06, y = 1.03, main, pos = 4, xpd = TRUE, cex = 0.8)
				if(!is.null(main)) text(x = -0.99, y = 1.03, main, pos = 2, xpd = TRUE, cex = 0.8, srt = 90)
				legend("topright", legend.txt, border = NA, bty = "o", cex = 0.7, col = legend.col, bg = NA, lty = legend.lty, pch = legend.pch, pt.cex = legend.pch.cex)
				dev.off()
			}
			# Create table with pvals, beta, maf etc. for all kmers plotted in the alignment figure
			if(correct_or_wrong=="correct_frame"){
				cols_to_keep = c("kmer","negLog10","beta","mac","maf","sstart")
				if(override_signif){
					out_table = cbind(res[,match(cols_to_keep, names(res))][which_to_align,])
				} else {
					out_table = cbind(res[,match(cols_to_keep, names(res))][which_to_align[1:length(which(as.numeric(res$negLog10)[which_to_align]>=bonferroni))],])
				}
				colnames(out_table)[which(colnames(out_table)=="sstart")] = "ps"
				out_table$ps = as.numeric(out_table$ps)-x.adjust
				out_table = cbind("gene" = gene_name, "plot" = paste0(gene_name, "_plot_",p,"_aminoacids_",(genestart-x.adjust),"_to_",(geneend-x.adjust)), out_table)
			} else {
				out_table = NULL
			}
		}
	}
	if(!is.null(out_table)) return(out_table)
}



get_top_genes = function(input_dir = NULL, prefix = NULL, ngenes = NULL, nsamples = NULL, bonferroni = NULL, ref.name = NULL){
	genes = system(paste0("ls ", input_dir, prefix, "_", kmer_type, kmer_length, "_", ref.name, "*topgene_*kmersandpvals.txt"), intern = T)
	# Get the order (from 1-n) of the genes
	genes_order = sapply(genes, function(x) unlist(strsplit(x, "/"))[length(unlist(strsplit(x, "/")))], USE.NAMES = F)
	genes_order = sapply(genes_order, function(x, prefix) unlist(strsplit(x, paste0("*topgene_")))[2], prefix = prefix, USE.NAMES = F)
	genes_order = sapply(genes_order, function(x) as.numeric(unlist(strsplit(x,"_"))[1]), USE.NAMES = F)
	# Order the files
	genes = genes[order(genes_order)]
	# Pull out the gene names from the files
	genes_names = sapply(1:length(genes), function(x, genes) unlist(strsplit(unlist(strsplit(genes[x], paste0("_topgene_", x, "_")))[2], "_kmersandpvals.txt")), genes  = genes, USE.NAMES = F)
	genes_all = genes

	cat("Read in top gene names", "\n")
	return(list("genes" = genes, "genes_names" = genes_names, "genes_all" = genes_all))
}


process_blast_protein = function(prefix = NULL, i = NULL, kmers_gene_i = NULL, genes_names = NULL, j = NULL, blastPath = NULL, correct_frame = NULL, genename_i = NULL, all_translations = NULL, kmer_type = NULL, kmer_length = NULL, perident = 70, nsamples = NULL, output_dir = NULL, ref.name = NULL){

	# Write reference translation to file
	if(j==correct_frame) correct_or_wrong = "correct_frame" else correct_or_wrong = "wrong_frame"
	header = paste0(">", genename_i, "_", correct_or_wrong,"_", j)
	ref_trans_file = paste0(output_dir, prefix, "_", kmer_type, kmer_length, "_", genename_i, "_", correct_or_wrong,"_", j, "_pro.fa")
	cat(c(header, all_translations[j]), file = ref_trans_file, sep = "\n")
	
	kmer_sequence_file = paste0(output_dir, prefix,"_", kmer_type, kmer_length, "_gene_", i, "_tmp_kmer_file.txt")
	blast_output_file = paste0(output_dir, prefix,"_", kmer_type, kmer_length, "_blast_out_",i,".txt")
	system(paste0(blastPath," -query ", kmer_sequence_file, " -subject ", ref_trans_file," -max_hsps 1 -outfmt '6 qseqid sseqid sacc pident length mismatch gapopen qstart qend evalue sstart send sseq qseq' -evalue 200000 -word_size 2 -gapopen 9 -gapextend 1 -matrix PAM30 -out ", blast_output_file))


	# Read in blast results
	blast.search=scan(blast_output_file,what=character(0),quiet = TRUE)
	# Is there are results:
	if(length(blast.search)!=0){
		# Create a matrix of the blast results
		blast.search = matrix(blast.search,ncol=14,byrow=TRUE)
		blast.search = blast.search[which(as.numeric(blast.search[,4])>=perident),]
		kmer.index = as.numeric(substr(as.character(blast.search[,1]), 5, 100000000))
		blast.search = cbind(kmers_gene_i[kmer.index,], blast.search)
		colnames(blast.search) = c("kmer","negLog10","beta","mac","qseqid","sseqid", "sacc", "pident", "length", "mismatch", "gapopen", "qstart", "qend" ,"evalue", "sstart","send","sseq","qseq")
		blast.search$maf = as.numeric(blast.search$mac)/nsamples
		blast_output_file_processed = paste0(output_dir, prefix, "_", kmer_type, kmer_length, "_", ref.name, "_top_gene_",i,"_", genes_names[i],"_", j, "_", correct_or_wrong, "_blast_results.txt")
		write.table(blast.search, file = blast_output_file_processed, row = F, col = T, sep = "\t", quote = F)
		system(paste0("rm ", ref_trans_file))
		if(j==6) system(paste0("rm ", kmer_sequence_file))
		system(paste0("rm ", blast_output_file))


	}

	return(blast_output_file_processed)

}

process_blast_nucleotide = function(blastPath = NULL, prefix = NULL, kmer_type = NULL, kmer_length = NULL, i = NULL, perident = NULL, kmers_gene_i = NULL, genes_names = NULL, ref_gene_i = NULL, genename_i = NULL, nsamples = NULL, output_dir = NULL, ref.name = NULL){
	
	header = paste0(">", genename_i)
	ref_file = paste0(output_dir, prefix, "_", kmer_type, kmer_length, "_", genename_i, "_nuc.fa")
	cat(c(header, ref_gene_i), file = ref_file, sep = "\n")
	
	kmer_sequence_file = paste0(output_dir, prefix,"_", kmer_type, kmer_length, "_gene_", i, "_tmp_kmer_file.txt")
	blast_output_file = paste0(output_dir, prefix,"_", kmer_type, kmer_length, "_blast_out_",i,".txt")

	
	# Run blast
	system(paste0(blastPath," -query ", kmer_sequence_file, " -subject ", ref_file," -max_hsps 1 -outfmt '6 qseqid sseqid sacc pident length mismatch gapopen qstart qend evalue sstart send sseq qseq sstrand' -evalue 0.05 -gapopen 5 -gapextend 2 -penalty -3 -reward 2 -word_size 4 -out ", blast_output_file))

	# Read in blast results
	blast.search=scan(blast_output_file,what=character(0),quiet = TRUE)
	# Is there are results:
	if(length(blast.search)!=0){
		# Create a matrix of the blast results
		blast.search = matrix(blast.search,ncol=15,byrow=TRUE)
		blast.search = blast.search[which(as.numeric(blast.search[,4])>=perident),]
		kmer.index = as.numeric(substr(as.character(blast.search[,1]), 5, 100000000))
		blast.search = cbind(kmers_gene_i[kmer.index,], blast.search)
		colnames(blast.search) = c("kmer","negLog10","beta","mac","qseqid","sseqid", "sacc", "pident", "length", "mismatch", "gapopen", "qstart", "qend" ,"evalue", "sstart","send","sseq","qseq", "sstrand")
		blast.search$maf = as.numeric(blast.search$mac)/nsamples
		blast_output_file_processed = paste0(output_dir, prefix, "_",kmer_type, kmer_length,  "_", ref.name, "_top_gene_",i,"_", genes_names[i],"_blast_results.txt")
		write.table(blast.search, file = blast_output_file_processed, row = F, col = T, sep = "\t", quote = F)
		system(paste0("rm ", ref_file))
		system(paste0("rm ", kmer_sequence_file))
		system(paste0("rm ", blast_output_file))
	}
	return(blast_output_file_processed)
	
}



plot_closeup_alignments = function(ref = NULL, ref_length = NULL, ref_gb = NULL, ref_fa = NULL, figures_dir = NULL, output_prefix = NULL, ngenes = 20, nsamples = NULL, bonferroni = NULL, gene_lookup = NULL, oneLetterCodes = NULL, kmer_type = NULL, kmer_length = NULL, blastPath = NULL, perident = NULL, col_lib_nuc = NULL, col_lib_pro = NULL, ref.name = NULL, alignmenttype = NULL, genes_all = NULL, override_signif = NULL, correct_only = TRUE){
	
	
	# Create alignment directory
	###########################
	
	alignment_dir = file.path(figures_dir, "alignments/")
	if(!dir.exists(alignment_dir)) dir.create(alignment_dir)

	
	# Read in reference files
	###########################
	gene_lookup = create_gene_lookup(ref = ref, ref_length = ref_length)$gene_lookup
	ref_gb = read_dna_seg_from_file(ref_gb, tagsToParse=c("CDS","repeat_region","tRNA","rRNA","ncRNA"))
	ref_fa = read_reference(ref_fa)	
	
	# Pull out the files containing the top genes
	###########################

	genes = genes_all$genes
	genes_names = genes_all$genes_names
	genes_all = genes_all$genes

	allgenes_results_table_out = list()
	
	for(i in 1:length(genes_names)){
		
		# What is gene i
		genename_i = genes_names[i]
		# If there is a colon in the name, replace with an underscore for directory name
		genename_i_dir = genename_i
		if(any(unlist(strsplit(genename_i_dir,""))==":")){
			genename_i_dir = paste(unlist(strsplit(genename_i_dir,":")), collapse = "_")
		}
		
		# Create gene subdirectory
		gene_dir = file.path(alignment_dir, paste0(genename_i_dir,"/"))
		if(!dir.exists(gene_dir)) dir.create(gene_dir)
	
		# Which row in gene_lookup contains genename_i
		wh_genelookup = which(gene_lookup[,1]==genename_i)
		if(length(wh_genelookup)==0) stop("Error: no match for gene name", genename_i, "\n")


		# # Create reference file for gene
		ref_gene_i = get_ref_gene_i(gene_lookup = gene_lookup, genename_i = genename_i, ref_length = ref_length, ref_fa = ref_fa, oneLetterCodes = oneLetterCodes)
		if(kmer_type=="nucleotide" & ref_gene_i$correct_frame!=1){
			ref_gene_i$ref_gene_i = paste(rev(unlist(revcompl_full[unlist(strsplit(ref_gene_i$ref_gene_i,""))])), collapse = "")
		}

		# Read in kmers
		kmers_gene_i = read_kmer_file(kmerfile = genes[i], i = i, prefix = output_prefix,
									kmer_type = kmer_type, kmer_length = kmer_length, output_dir = gene_dir)
		
		if(kmer_type=="protein"){
			
			gene_i_results_list = list()

			for(j in 1:6){

				blast.search = process_blast_protein(prefix = output_prefix, i = i,
								kmers_gene_i = kmers_gene_i, genes_names = genes_names,
								j = j, blastPath = blastPath,
								correct_frame = ref_gene_i$correct_frame,
								genename_i = genename_i,
								all_translations = ref_gene_i$all_translations,
								kmer_type = kmer_type, kmer_length = kmer_length,
								perident = perident, nsamples = nsamples,
								output_dir = gene_dir, ref.name = ref.name)

				res = read_res_table(resfile = blast.search, refseq = ref_gene_i$all_translations[j],
								i = i, genes_names = genes_names, perident = perident, kmer_type = kmer_type)$res

				gene_i_results_list[[j]] = res

				if((correct_only & j==ref_gene_i$correct_frame) | correct_only==FALSE){
					# Plot the alignment figures
					out_results_correct_frame = run_alignment_nplots_protein(ref_gene_i = ref_gene_i, res = res,
									nsamples = nsamples, bonferroni = bonferroni, prefix = output_prefix,
									gene_name = genes_names[i], j = j, col_lib = col_lib_pro,
									minor_allele_threshold = minor_allele_threshold,
									macormaf = macormaf, output_dir = gene_dir,
									kmer_type = kmer_type, kmer_length = kmer_length,
									ref.name = ref.name, override_signif = override_signif)
					if(j==ref_gene_i$correct_frame) allgenes_results_table_out = rbind(allgenes_results_table_out, out_results_correct_frame)
				}	
			}

		} else {

			blast.search = process_blast_nucleotide(blastPath = blastPath, prefix = output_prefix,
												kmer_type = kmer_type, kmer_length = kmer_length,
												i = i, perident = perident, kmers_gene_i = kmers_gene_i,
												genes_names = genes_names, ref_gene_i = ref_gene_i$ref_gene_i,
												genename_i = genename_i, nsamples = nsamples,
												output_dir = gene_dir, ref.name = ref.name)

			res = read_res_table(resfile = blast.search, refseq = ref_gene_i$ref_gene_i, 
										i = i, genes_names = genes_names, perident = perident, 
										kmer_type = kmer_type, kmer_length = kmer_length)$res

			gene_i_results_list = list(); gene_i_results_list[[1]] = res
		
			out_results = run_alignment_nplots_nucleotide(ref_gene_i = ref_gene_i, res = gene_i_results_list[[1]],
										nsamples = nsamples, bonferroni = bonferroni, prefix = output_prefix,
										gene_name = genename_i, gene_lookup = gene_lookup,
										wh_genelookup = wh_genelookup, col_lib = col_lib_nuc,
										minor_allele_threshold = minor_allele_threshold,
										macormaf = macormaf, output_dir = gene_dir,
										kmer_type = kmer_type, kmer_length = kmer_length,
										ref.name = ref.name, override_signif = override_signif)

			allgenes_results_table_out = rbind(allgenes_results_table_out, out_results)

		

		
		}

		# Find the kmers which haven't aligned to any reading frame
		which_kmers_no_result = get_kmers_noresult(gene_i_results_list = gene_i_results_list,
								kmers_gene_i = kmers_gene_i, prefix = output_prefix,
								i = i, genes_names = genes_names,
								kmer_type = kmer_type, kmer_length = kmer_length,
								output_dir = gene_dir, ref.name = ref.name)
		
		# Plot Manhattans per frame
		if(kmer_type=="protein"){
			for(j in 1:6){
				if(((correct_only & j==ref_gene_i$correct_frame) | correct_only==FALSE)){
					run_manhattan_single_protein(which_kmers_no_result = which_kmers_no_result,
									res = gene_i_results_list[[j]], ref_gene_i = ref_gene_i,
									kmer_length = kmer_length, prefix = output_prefix,
									gene_name = genes_names[i], j = j,
									bonferroni = bonferroni,
									ref_gb_full = ref_gb, ref_length = ref_length,
									nsamples = nsamples, kmer_type = kmer_type,
									minor_allele_threshold = minor_allele_threshold, macormaf = macormaf,
									output_dir = gene_dir, ref.name = ref.name)
				}
			}
		
			run_manhattan_allframes(gene_i_results_list = gene_i_results_list,
								prefix = output_prefix, gene_name = genes_names[i],
								ref_gene_i = ref_gene_i,
								which_kmers_no_result = which_kmers_no_result,
								bonferroni = bonferroni,
								minor_allele_threshold = minor_allele_threshold, macormaf = macormaf,
								output_dir = gene_dir,
								kmer_type = kmer_type, kmer_length = kmer_length, ref.name = ref.name,
								ref_gb_full = ref_gb)
		
		} else {
			run_manhattan_single_nucleotide(which_kmers_no_result = which_kmers_no_result, res = gene_i_results_list[[1]], ref_gene_i = ref_gene_i, prefix = output_prefix, gene_name = genes_names[i], bonferroni = bonferroni, ref_gb_full = ref_gb, ref_length = ref_length, kmer_type = kmer_type, kmer_length = kmer_length, nsamples = nsamples, minor_allele_threshold = minor_allele_threshold, macormaf = macormaf, output_dir = gene_dir, ref.name = ref.name)
		}



	}
	
	write.table(allgenes_results_table_out, file = paste0(alignment_dir, output_prefix, "_", kmer_type, kmer_length, "_", ref.name, "_", alignmenttype, "_all_top_genes_significant_kmers_per_alignment_plot.txt"), row = F, col = T, sep = "\t", quote = F)


}
