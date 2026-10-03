# plot_figures.R: draw the pipeline's figures from the figure-data files written by
# plotManhattan.py (and plotManhattanbowtie.py). Python does all the computation; this
# script only draws (PLAN 5.5, amendment A1).
#
# Run it through Rscript_launcher.R, which reports errors with file:line for each call
# and exits with status 1:
#   Rscript --vanilla Rscript_launcher.R plot_figures.R --data-dir <figures_dir>/figure_data
#
# The drawing functions below are copied unchanged from alignmentfunctions.R and
# Manhattan_functions.R at r-fixed (ba418d6), except that set.seed(0) now precedes each
# sample() call (jitter of k-mers without a BLAST hit), so figures are reproducible.

###################################################################################################
## Drawing functions (from Manhattan_functions.R and alignmentfunctions.R, r-fixed)
###################################################################################################

colour_selection = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00")
get_genes_to_plot = function(gene_names = NULL, y = NULL, gene_conversion = NULL, ymax = NULL, gene_panel = NULL, ref = NULL, xadjust = NULL, ngenes = 20){
	if(!is.null(gene_names)){
		o = order(as.numeric(y), decreasing = T)
		if(is.null(y)) top_genes = gene_names else top_genes = head(unique(as.character(gene_names[o])), ngenes)
		if(!is.null(xadjust)) xadjust = xadjust[order(match(top_genes, ref[,"name"]))]
		top_genes = top_genes[order(match(top_genes, ref[,"name"]))]
		gene_name_conversion = sapply(top_genes, function(x) as.character(gene_conversion[x]), USE.NAMES = F)
		gene_name_conversion[which(is.na(gene_name_conversion))] = top_genes[which(is.na(gene_name_conversion))]
		gene_col = rep("black", length(top_genes))
		gene_col[which(!is.na(match(gene_name_conversion, gene_panel)))] = "#D55E00"
		wh_intergenicmatch = sapply(gene_name_conversion, function(x, gene_panel) length(which(!is.na(match(unlist(strsplit(x,":")), gene_panel)))), gene_panel = gene_panel, USE.NAMES = F)
		if(any(wh_intergenicmatch)>0){
			if(any(wh_intergenicmatch==1 & gene_col!="#D55E00")) gene_col[which(wh_intergenicmatch==1 & gene_col!="#D55E00")] = "#E69F00"
			if(any(wh_intergenicmatch==2)) gene_col[which(wh_intergenicmatch==2)] = "#D55E00"
		}
		if(is.null(xadjust)) xadjust = rep(0, length(top_genes))
		gene_lines_to_plot = cbind("genes" = as.character(top_genes),
			"ytop" = rep(c((ymax[1]+(ymax[2]/40)), (ymax[1]+(ymax[2]/12))), length(top_genes))[1:length(top_genes)],
			"xadjust" = xadjust,
			"replace_gene_name" = as.character(gene_name_conversion),
			"gene_col" = gene_col)
	}
	return(gene_lines_to_plot)
}
plot_gene_lines = function(genes = NULL, ytop = NULL, col = NULL, ref = NULL, rect = TRUE, ytext = 0, line.angle = 2, line.gap.y = 0.2, line.gap.x = 1000, xadjust = 0, ybottom = 0, replace_gene_name = NULL, line.length = 150000, gene_name_col = NULL, gene_name_cex = 0.4){
	if(length(ytop)==1) ytop = rep(ytop, length(genes))
	if(length(col)==1) col = rep(col, length(genes))
	if(length(xadjust)==1) xadjust = rep(xadjust, length(genes)) 
	if(length(ybottom)==1) ybottom = rep(ybottom, length(genes))
	if(length(line.length)==1) line.length = rep(line.length, length(genes))
	if(is.null(gene_name_col)) gene_name_col = rep("black", length(genes))
	# o = order(ytop, decreasing = T)
	# genes = genes[o]; ytop = ytop[o]; col = col[o]
	for(i in 1:length(genes)){
		if(!is.na(genes[i])){
			if(!any(unlist(strsplit(genes[i],""))==":")){
				pos1 = as.numeric(ref[,"start"][which(ref[,"name"]==genes[i])[1]])
				pos2 = as.numeric(ref[,"end"][which(ref[,"name"]==genes[i])[length(which(ref[,"name"]==genes[i]))]])
			} else {
				genes.i = as.character(unlist(strsplit(genes[i],":")))
				pos1 = as.numeric(ref[,"end"][which(ref[,"name"]==genes.i[1])[length(which(ref[,"name"]==genes.i[1]))]])+1
				pos2 = as.numeric(ref[,"start"][which(ref[,"name"]==genes.i[2])[1]])-1
			}
			if(rect){
				rect(xleft = pos1, xright = pos2, ybottom = 0, ytop = ytop[i], border = NA, col = col[i], xpd = T)
			} else {
				# cat("lines x:",c((pos1+(pos2-pos1)/2), (pos1+(pos2-pos1)/2)), "\n")
				# cat("lines y:", c(ybottom[i], ytop[i]), "\n")
				lines(x = c((pos1+(pos2-pos1)/2), (pos1+(pos2-pos1)/2)), y = c(ybottom[i], ytop[i]), lty = 3, col = col[i], xpd = T)
			}
			if(!is.null(replace_gene_name)){
				if(replace_gene_name[i]!=""){
					genes[i] = replace_gene_name[i]
				}
			}
			text(x = (pos1+(pos2-pos1)/2)+xadjust[i], y = ytop[i]+ytext, srt = 45, labels = genes[i], adj = 0, cex = gene_name_cex, xpd = T, col = gene_name_col[i])
		}
	}
}
get_legend_col_manhattan = function(beta = NULL, legend.xpos = NULL, legend.ypos = NULL, text.xpos = NULL, text.ypos = NULL){
	bluegrey = colorRamp(c(colour_selection[5], "grey50"))
	greyred = colorRamp(c("grey50", colour_selection[6]))
	testcol1 = seq(from = 0,by = 0.01, length.out=100)
	testcol1 = bluegrey(testcol1)
	testcol1 = rgb(testcol1, maxColorValue = 256)
	testcol2 = seq(from = 0,by = 0.01, length.out=100)
	testcol2 = greyred(testcol2)
	testcol2 = rgb(testcol2, maxColorValue = 256)
	par(fig = c(legend.xpos[1], ((legend.xpos[2]-legend.xpos[1])/2)+legend.xpos[1], legend.ypos[1], legend.ypos[2]), mar=c(0,0,0,0), new=TRUE)
	image(c(1:100),1, (matrix(c(1:100),ncol=1,nrow=100)), col=testcol1, axes=FALSE)
	axis(1, at = c(1,100), labels = c(NA,NA), xpd = T, tck = -0.1, cex.axis = 0.5, lwd = 0.8)
	axis(1, at = c(1,100), labels = c(round(min(beta, na.rm = T)),0), xpd = T, tck = 0, cex.axis = 0.5, lwd = 0, line = -1.3)
	par(fig = c(0,1,0,1), mar=c(0,0,0,0), new = TRUE)
	plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
	par(fig = c(0.9762, 0.9767,0.949, 0.954), mar=c(0,0,0,0), new=TRUE)
	par(fig = c(((legend.xpos[2]-legend.xpos[1])/2)+legend.xpos[1], legend.xpos[2],legend.ypos[1], legend.ypos[2]), mar=c(0,0,0,0), new=TRUE)
	image(c(1:100),1, (matrix(c(1:100),ncol=1,nrow=100)), col=testcol2, axes=FALSE)
	axis(1, at = c(1,100), labels = c(NA,NA), xpd = T, tck = -0.1, cex.axis = 0.5, lwd = 0.8)
	axis(1, at = c(1,100), labels = c(0,round(max(beta, na.rm = T))), xpd = T, tck = 0, cex.axis = 0.5, lwd = 0, line = -1.3)
	par(fig = c(0.9802, 0.9807,0.949, 0.954), mar=c(0,0,0,0), new=TRUE)	
}
plot_QQ = function(kmerIndex = NULL, assoc = NULL, output_dir = NULL, prefix = NULL, minor_allele_threshold = NULL, macormaf = NULL, mapatterns = NULL, kmer_type = NULL, kmer_length = NULL){
	
	# Get expected and empirical for LMM p-values
	if(minor_allele_threshold==0) which_kmers = which(mapatterns[unique(kmerIndex)]>0) else which_kmers = which(mapatterns[unique(kmerIndex)]>=minor_allele_threshold)
	qqplot.x = -log10((1:length(which_kmers))/length(which_kmers))
	qqplot.y = as.numeric(assoc[,6])[unique(kmerIndex)[which_kmers]]
	qqplot.y = qqplot.y[order(qqplot.y, decreasing = T)]
	
	if(minor_allele_threshold==0) file_suffix = "_QQplot_allkmers.png" else file_suffix = paste0("_QQplot_", macormaf, minor_allele_threshold, ".png")
	png(paste0(output_dir, prefix, "_", kmer_type, kmer_length, file_suffix), width = 12, height = 12, units = "cm", res = 600)
	par("mar" = c(5.1, 4.1, 1, 1))
	plot(x = qqplot.x, y = qqplot.y, xlab=expression(paste("Null distribution of -log"[10],italic(' p')," values",collapse="")), ylab = expression(paste("Empirical distribution of -log"[10],italic(' p')," values",collapse="")), cex.axis = 0.8, cex.lab = 0.8, type = "l", log = "")
	abline(0,1,col = "red", lty = 2)		
	dev.off()
	
}
get_Manhattan_colours = function(final_kmer_pos_index = NULL, assoc_patterns = NULL, kmerIndex = NULL, colour_selection = NULL, ypos = NULL, bonferroni = NULL, mafpatterns = NULL, pheno_type = NULL){
	
	## Colour by whether the kmer has mapped more than once
	multialignCOL = rep("grey50", length(final_kmer_pos_index))
	matchcount = table(as.numeric(final_kmer_pos_index))
	matchcount = matchcount[which(matchcount>1)]
	multialignCOL[which(!is.na(match(final_kmer_pos_index, as.numeric(names(matchcount)))))] = colour_selection[6]
	cat("Created multialignCOL","\n")
	
	
	# Colour by beta (kmers above significance threshold)
	cat("Range beta:", range(as.numeric(assoc_patterns[,2]), na.rm = T), "\n")
	beta = as.numeric(assoc_patterns[,2])[kmerIndex[final_kmer_pos_index]]
	betaCOL = rep("grey50", length(ypos))
	if(pheno_type=="binary"){
		betaCOL[which(beta>0)] = colour_selection[6]
		betaCOL[which(beta<0)] = colour_selection[5]
	} else {
		greyred = colorRamp(c("grey50", colour_selection[6]))
		bluegrey = colorRamp(c(colour_selection[5], "grey50"))
		betapos = beta[which(beta>0)]
		betapos = (betapos-min(betapos, na.rm = T))/(max(betapos, na.rm = T)-min(betapos, na.rm = T))
		betaneg = beta[which(beta<0)]
		betaneg = (betaneg-min(betaneg, na.rm = T))/(max(betaneg, na.rm = T)-min(betaneg, na.rm = T))
		bposcols = greyred(betapos); bposcols = rgb(bposcols, maxColorValue = 256)
		bnegcols = bluegrey(betaneg); bnegcols = rgb(bnegcols, maxColorValue = 256)
		betaCOL[which(beta>0)] = bposcols
		betaCOL[which(beta<0)] = bnegcols
	}
	betaCOL[which(ypos<bonferroni)] = "grey50"
	cat("Created betaCOL","\n")
	
	
	# Colour by MAF
	maf = mafpatterns[kmerIndex[final_kmer_pos_index]]
	mafCOL = rep("grey50", length(final_kmer_pos_index))
	mafCOL[which(maf<0.01)] = colour_selection[6]
	mafCOL[which(maf>=0.01 & maf<0.05)] = colour_selection[5]
	mafCOL[which(maf>=0.05)] = colour_selection[3]
	cat("Created mafCOL","\n")
	
	return(list("multialignCOL" = multialignCOL, "betaCOL" = betaCOL, "mafCOL" = mafCOL))
	
}
plot_manhattan = function(outfilename = NULL, xpos = NULL, ma_threshold_pass = NULL, ypos = NULL, ylims.i = NULL, annotateGeneFile = NULL, ref = NULL, which_genes_to_annotate.i = NULL, allCOLS = NULL, allPCH = NULL, i = NULL, bonferroni = NULL, legendtext = NULL, legendcol = NULL, legendpch = NULL, legendlty = NULL, beta = NULL, gene_names = NULL, gene_conversion = NULL, pheno_type = NULL, ref_length = NULL){
	
			
	
	png(outfilename, width = 22, height = 12, units = "cm", res = 600)
	par(mar=c(4.1,4.1,3,7.7))
	plot(x = xpos[ma_threshold_pass], y = ypos[ma_threshold_pass], col = "grey50", cex = 0.5, cex.lab = 0.8, cex.axis = 0.8, xlab = "", ylab = "", axes = F, type = "n", ylim = ylims.i)
	ymax = c(par("usr")[4], (par("usr")[4]-par("usr")[3]))
	
	if(!is.null(annotateGeneFile)){
		# annotateGene = read.table(annotateGeneFile, h = F, sep = "\t", as.is = T)
		# if(ncol(annotateGene)>1) annotateGeneXadjust = as.character(annotateGene[,2]) else annotateGeneXadjust = rep(0,nrow(annotateGene))
		# annotateGene = as.character(annotateGene[,2])
		# if(any(annotateGeneXadjust)=="") annotateGeneXadjust[which(annotateGeneXadjust=="")] = 0
		# annotateGeneXadjust = as.numeric(annotateGeneXadjust)
		annotateGene = scan(annotateGeneFile, what = character(0), sep = "\n", quiet = TRUE)
		annotateGeneXadjust = rep(0,length(annotateGene))
		annotateGeneConversion = annotateGene
		names(annotateGeneConversion) = annotateGeneConversion
		cat("Genes/IRs to annotate on the Manhattan plot:",paste(annotateGene, collapse = " "),"\n")
		gene_lines_to_plot = get_genes_to_plot(gene_names = annotateGene, y = NULL, gene_conversion = annotateGeneConversion, ymax = ymax, gene_panel = c(), ref = ref, xadjust = annotateGeneXadjust)
	} else {
		gene_lines_to_plot = get_genes_to_plot(gene_names = gene_names[which_genes_to_annotate.i], y = ypos[which_genes_to_annotate.i], gene_conversion = gene_conversion[which_genes_to_annotate.i], ymax = ymax, gene_panel = c(), ref = ref)
	}

	
	plot_gene_lines(genes = as.character(gene_lines_to_plot[,1]), col = "#cecece", ytop = as.numeric(gene_lines_to_plot[,2]), ref = ref, rect = FALSE, ytext = 0, line.angle = 1.62, line.gap.y = 0.165, line.gap.x = 15000, xadjust = as.numeric(gene_lines_to_plot[,3]), ybottom = 0, replace_gene_name = as.character(gene_lines_to_plot[,4]), gene_name_col = as.character(gene_lines_to_plot[,5]), gene_name_cex = 0.6)
	
	points(x = xpos[ma_threshold_pass], y = ypos[ma_threshold_pass], col = allCOLS[[i]][ma_threshold_pass], cex = 0.5, pch = allPCH[[i]][ma_threshold_pass])
	
	mtext("Position in reference genome (Mb)", side = 1, line = 2.5, cex = 0.8)
	mtext(expression(paste("Significance (-log"[10],italic(' p'),") LMM",collapse="")), side = 2, line = 2.8, cex = 0.8)
	axis(1, cex.axis = 0.8, at = c(0:floor(ref_length/1e6))*1e6, labels = as.character(0:floor(ref_length/1e6)))
	axis(2, cex.axis = 0.8)
	abline(h = bonferroni, col = "black", lty = 2)
	par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
	plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
	legend_pos = c(0.71, 0.99)
	if(filecol[i]=="betaCOL" & pheno_type=="continuous"){
		legend(legend_pos[1], legend_pos[2], legendtext[[i]][-c(3:5)], col = legendcol[[i]][-c(3:5)], pch = legendpch[-c(3:5)], lty = legendlty[-c(3:5)], bty = "n", cex = 0.65, xpd = TRUE, pt.bg = "#949494", lwd = 1)
		get_legend_col_manhattan(beta = c(min(beta, na.rm = T),max(beta, na.rm = T)), legend.xpos = c(0.828, 0.858), legend.ypos = c(0.84, 0.86), text.xpos = c(0.8,0.9), text.ypos = c(0.8,0.9))
		par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
		plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
		text(x = 0.797, y = 0.755, label = "\u03B2", cex = 0.65, xpd = T)
		rect(xleft = 0.68, xright = 1.07, ybottom = 0.64, ytop = 0.99)
	} else {
		legend(legend_pos[1], legend_pos[2], legendtext[[i]], col = legendcol[[i]], pch = legendpch, lty = legendlty, bty = "o", cex = 0.65, xpd = TRUE, pt.bg = "#949494", lwd = 1)
	}
	dev.off()

	
	
}
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
rev_compl_col = list("#008000" = "#ff9000", "#ff9000" = "#008000", "blue" = "black", "black" = "blue", "#d3d3d3" = "#d3d3d3")
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
		set.seed(0); xpos_unmapped = sample(seq(from = ref_gene_i$length_correct+50, to = ref_gene_i$length_correct+100, length.out = nrow(which_kmers_no_result)), nrow(which_kmers_no_result), replace = F)
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
	plot_singleframe_manhattan_protein(prefix = prefix, genes_name = gene_name, correct_or_wrong = correct_or_wrong, j = j, xpos = xpos[whichMAthreshold,,drop=FALSE], ypos = ypos[whichMAthreshold], all_translations = ref_gene_i$all_translations, beta_col = beta_col[whichMAthreshold], bonferroni = bonferroni, ylim_max = NULL, maname = paste0("_",macormaf,minor_allele_threshold), x.adjust = 999, ref_gb_full = ref_gb_full, start = ref_gene_i$ref_start_i, end = ref_gene_i$ref_end_i, ref_length = ref_length, kmer_length = kmer_length)

	# Then limit to ylim = 50 for those with ylim > 100
	if(max(as.numeric(ypos))>100){
		plot_singleframe_manhattan_protein(prefix = prefix, genes_name = gene_name, correct_or_wrong = correct_or_wrong, j = j, xpos = xpos, ypos = ypos, all_translations = ref_gene_i$all_translations, beta_col = beta_col, bonferroni = bonferroni, ylim_max = 50, x.adjust = 999, ref_gb_full = ref_gb_full, start = ref_gene_i$ref_start_i, end = ref_gene_i$ref_end_i, ref_length = ref_length, kmer_length = kmer_length)
	}
	if(max(as.numeric(ypos)[whichMAthreshold])>100){
		plot_singleframe_manhattan_protein(prefix = prefix, genes_name = gene_name, correct_or_wrong = correct_or_wrong, j = j, xpos = xpos[whichMAthreshold,,drop=FALSE], ypos = ypos[whichMAthreshold], all_translations = ref_gene_i$all_translations, beta_col = beta_col[whichMAthreshold], bonferroni = bonferroni, ylim_max = 50, maname = paste0("_",macormaf,minor_allele_threshold), x.adjust = 999, ref_gb_full = ref_gb_full, start = ref_gene_i$ref_start_i, end = ref_gene_i$ref_end_i, ref_length = ref_length, kmer_length = kmer_length)
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
		set.seed(0); xpos_unmapped = sample(seq(from = nchar(ref_gene_i$ref_gene_i)+50, to = nchar(ref_gene_i$ref_gene_i)+100, length.out = nrow(which_kmers_no_result)), nrow(which_kmers_no_result), replace = F)
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
	plot_singleframe_manhattan_nucleotide(prefix = prefix, genes_name = gene_name, xpos = xpos[whichMAthreshold,,drop=FALSE], ypos = ypos[whichMAthreshold], ref_fa = ref_gene_i$ref_gene_i, beta_col = beta_col[whichMAthreshold], bonferroni = bonferroni, ylim_max = NULL, maname = paste0("_",macormaf,minor_allele_threshold), x.adjust = 999, ref_gb_full = ref_gb_full, start = ref_start_i, end = ref_end_i)
	
	if(max(as.numeric(ypos))>100){
		plot_singleframe_manhattan_nucleotide(prefix = prefix, genes_name = gene_name, xpos = xpos, ypos = ypos, ref_fa = ref_gene_i$ref_gene_i, beta_col = beta_col, bonferroni = bonferroni, ylim_max = 50, x.adjust = 999, ref_gb_full = ref_gb_full, start = ref_start_i, end = ref_end_i)
	}
	if(max(as.numeric(ypos)[whichMAthreshold])>100){
		plot_singleframe_manhattan_nucleotide(prefix = prefix, genes_name = gene_name, xpos = xpos[whichMAthreshold,,drop=FALSE], ypos = ypos[whichMAthreshold], ref_fa = ref_gene_i$ref_gene_i, beta_col = beta_col[whichMAthreshold], bonferroni = bonferroni, ylim_max = 50, maname = paste0("_",macormaf,minor_allele_threshold), x.adjust = 999, ref_gb_full = ref_gb_full, start = ref_start_i, end = ref_end_i)
	}


	
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

		out_allmaf = plot_alignment_function_protein(genestart = nplots[p,1], geneend = nplots[p,2], minor_allele_threshold = 0, res = res, nsamples = nsamples, maname = "", bonferroni = bonferroni, translation = ref_gene_i$all_translations[j], prefix = prefix, gene_name = gene_name, correct_or_wrong = correct_or_wrong, j = j, p = p, lowfreq = minor_allele_threshold, plot_ref = FALSE, main = "All protein kmers", x.adjust = 333, nplots = nplots, override_signif = override_signif, plot_figures = plot_figures, col_lib = col_lib, macormaf = macormaf, output_dir = output_dir, kmer_type = kmer_type, kmer_length = kmer_length, ref.name = ref.name)

		out_mafthreshold = plot_alignment_function_protein(genestart = nplots[p,1], geneend = nplots[p,2], minor_allele_threshold = minor_allele_threshold, res = res, nsamples = nsamples, maname = paste0("_", macormaf,minor_allele_threshold), bonferroni = bonferroni, translation = ref_gene_i$all_translations[j], prefix = prefix, gene_name = gene_name, correct_or_wrong = correct_or_wrong, j = j, p = p, plot_ref = FALSE, main = paste0("Protein kmers \u2265 ",macormaf," ", minor_allele_threshold), x.adjust = 333, nplots = nplots, override_signif = override_signif, plot_figures = plot_figures, col_lib = col_lib, macormaf = macormaf, output_dir = output_dir, kmer_type = kmer_type, kmer_length = kmer_length, ref.name = ref.name)
		if(!is.null(out_allmaf)) out_results_correct_frame = rbind(out_results_correct_frame, out_allmaf)

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

		if(!is.null(out_allmaf)) allgenes_results_table_out = rbind(allgenes_results_table_out, out_allmaf)

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
	# Bug fix in indices in loop below DJW 31 May 2022
	for(i in 1:length(genes)){
		start = ref$start[i]
		end = ref$end[i]
		rev.compl = ref$strand[i]==-1
		if(plot_names[i]) name = genes[i] else name = NA
		if(!is.null(name_replace)){
			if(plot_names[i] & !is.na(name_replace[i])) name = name_replace[i] else name = name
		}
		draw_arrow(start = start, end = end, arrow_length = arrow_length, height1 = height1, height2 = height2, arrowdiff = arrowdiff, fillCOL = fillCOL[i], border = "black", name = name, text_adjust = text_adjust[i], rev.compl = rev.compl, text.cex = text.cex, lwd = lwd, text.col = text.col[i])
	}
	for(i in 1:(length(genes)-1)){
		start = ref$end[i]
		end = ref$start[i+1]
		if(draw_line) lines(x = c(start, end), y = c(height1-((height1-height2)/2), height1-((height1-height2)/2)), lwd = lwd, xpd = T)
	}

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

		if(is.null(malim)) which_to_plot = seq(1,by=1,len=length(ypos_list[[f]])) else which_to_plot = which(c(as.numeric(gene_i_results_list[[f]][[macormaf]]))>=(malim))

		for(k in which_to_plot){
			lines(x = c(as.numeric(xpos_list[[f]][k,1]), as.numeric(xpos_list[[f]][k,2])), y = rep(as.numeric(ypos_list[[f]][k]), 2), col = frame_cols[f])
		}
	}

	# Plot unmapped kmers
	if(!is.null(which_kmers_no_result)){
		set.seed(0); xpos_unmapped = sample(seq(from = end+50, to = end+100, length.out = nrow(which_kmers_no_result)), nrow(which_kmers_no_result), replace = F)
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
	if(length(which_aligned)==0) xpos_1_new_aligned = c()
	xpos_2_new_aligned = sapply(xpos[which_aligned,2], function(x, s) s[x], s = seq_pos2, USE.NAMES = F)
	if(length(which_aligned)==0) xpos_2_new_aligned = c()
	set.seed(0); xpos_1_new_unaligned = sample(seq(from = end+50, to = end+100, length.out = length(which_unaligned)), length(which_unaligned), replace = F)
	if(length(which_unaligned)==0) xpos_1_new_unaligned = c()
	xpos_2_new_unaligned = xpos_1_new_unaligned+((kmer_length*3)-1)
	if(length(which_unaligned)==0) xpos_2_new_unaligned = c()

	xpos_1_new = rep(0, nrow(xpos)); if(length(which_aligned)>0) xpos_1_new[which_aligned] = xpos_1_new_aligned; if(length(which_unaligned)>0) xpos_1_new[which_unaligned] = xpos_1_new_unaligned
	xpos_2_new = rep(0, nrow(xpos)); if(length(which_aligned)>0) xpos_2_new[which_aligned] = xpos_2_new_aligned; if(length(which_unaligned)>0) xpos_2_new[which_unaligned] = xpos_2_new_unaligned

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
	if(length(xpos)>0) {
		plot(range(as.numeric(as.vector(xpos))), c(0, max(as.numeric(ypos))), type = "n", xlab = "Amino acid position in reference", ylab = expression(paste("Significance (-log"[10],italic(' p'),") LMM",collapse="")), main = paste0(genes_name, " (Protein kmers)"), axes = F, ylim = ylim)
		# plot_axes(all_translations[j])
		for(k in 1:nrow(xpos)){
			lines(x = c(as.numeric(xpos[k,1]), as.numeric(xpos[k,2])), y = rep(as.numeric(ypos[k]), 2), col = beta_col[k])
		}
		abline(h = bonferroni, col = "black", lty = 2)

		plot_gene_arrows_manhattan(ref_gb_full = ref_gb_full, start = start, end = end, xpos = xpos)


		plot_beta_and_feature_legend()
	}
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


###################################################################################################
## Reading the figure data (written by plotManhattan.py / plotManhattanbowtie.py)
###################################################################################################

FIGURE_DATA_VERSION = "kmer_pipeline figure data v1"
NA_SENTINEL = "__NA__"

# Read the key/value parameters file
read_params = function(file) {
	if(!file.exists(file)) stop("Figure data file missing: ", file)
	l = readLines(file, warn = FALSE)
	if(length(l) == 0 || l[1] != paste("#", FIGURE_DATA_VERSION)) stop("Not a figure data file of format '", FIGURE_DATA_VERSION, "': ", file)
	l = l[-1]
	k = sub("\t.*$", "", l)
	v = sub("^[^\t]*\t", "", l)
	v[v == NA_SENTINEL] = NA
	setNames(v, k)
}

# Read one table: line 1 the format version, line 2 the columns as name:type
read_fd = function(dir, name) {
	file = file.path(dir, paste0(name, ".tsv.gz"))
	if(!file.exists(file)) stop("Figure data file missing: ", file)
	con = gzfile(file)
	head = readLines(con, n = 2, warn = FALSE)
	close(con)
	if(length(head) < 2 || head[1] != paste("#", FIGURE_DATA_VERSION)) stop("Not a figure data file of format '", FIGURE_DATA_VERSION, "': ", file)
	cols = strsplit(head[2], "\t", fixed = TRUE)[[1]]
	cname = sub(":[a-z]+$", "", cols)
	ctype = sub("^.*:", "", cols)
	if(!all(ctype %in% c("character", "numeric", "integer", "logical"))) stop("Unknown column type in ", file, ": ", paste(cols, collapse = " "))
	what = lapply(ctype, function(t) vector(t, 0))
	names(what) = cname
	con = gzfile(file)
	x = scan(con, what = what, sep = "\t", skip = 2, quote = "", na.strings = NA_SENTINEL, comment.char = "", quiet = TRUE, allowEscapes = FALSE, strip.white = FALSE, blank.lines.skip = TRUE)
	close(con)
	data.frame(x, stringsAsFactors = FALSE, check.names = FALSE)
}

need = function(p, keys, file) {
	missing = setdiff(keys, names(p))
	if(length(missing)) stop("Missing parameters in ", file, ": ", paste(missing, collapse = ", "))
}

###################################################################################################
## Figure families
###################################################################################################

# QQ plots: plot_QQ() as plotManhattan.Rscript calls it
draw_qq = function(dd, p) {
	patterns = read_fd(dd, "patterns")
	kmers = read_fd(dd, "kmers")
	assoc = matrix(NA_real_, nrow(patterns), 6)
	assoc[, 2] = patterns$beta
	assoc[, 6] = patterns$neglog10p
	minor_allele_threshold = as.numeric(p[["minor_allele_threshold"]])
	for(thr in c(0, minor_allele_threshold)) {
		plot_QQ(kmerIndex = kmers$kmer_index, assoc = assoc, output_dir = p[["figures_dir"]], prefix = p[["output_prefix"]], minor_allele_threshold = thr, macormaf = p[["macormaf"]], mapatterns = patterns$ma, kmer_type = p[["kmer_type"]], kmer_length = as.numeric(p[["kmer_length"]]))
	}
}

# Genome-wide Manhattan plots: the drawing part of plotManhattan.Rscript (r-fixed), unchanged
# apart from reading its inputs and seeding the subsample
draw_genome_manhattan = function(dd, p) {
	patterns = read_fd(dd, "patterns")
	kmers = read_fd(dd, "kmers")
	positions = read_fd(dd, "positions")
	ref = read_fd(dd, "reference_cds")
	assoc = matrix(NA_real_, nrow(patterns), 6)
	assoc[, 2] = patterns$beta
	assoc[, 6] = patterns$neglog10p
	mafpatterns = patterns$maf
	kmerIndex = kmers$kmer_index
	ma = kmers$ma
	final_kmer_pos_index = positions$kmer
	final_kmer_pos = positions$position
	final_kmer_genes = positions$gene
	bonferroni = as.numeric(p[["bonferroni"]])
	minor_allele_threshold = as.numeric(p[["minor_allele_threshold"]])
	macormaf = p[["macormaf"]]
	pheno_type = p[["pheno_type"]]
	ref_length = as.numeric(p[["ref_length"]])
	annotateGeneFile = if(file.exists(file.path(dd, "annotate_genes.txt"))) file.path(dd, "annotate_genes.txt") else NULL

	## Get y position
	ypos = as.numeric(assoc[,6])[kmerIndex[final_kmer_pos_index]]
	cat("Got ypos","\n")

	kmerCOLS = get_Manhattan_colours(final_kmer_pos_index = final_kmer_pos_index, assoc_patterns = assoc, kmerIndex = kmerIndex, colour_selection = colour_selection, ypos = ypos, bonferroni = bonferroni, mafpatterns = mafpatterns, pheno_type = pheno_type)
	multialignCOL = kmerCOLS$multialignCOL
	betaCOL = kmerCOLS$betaCOL
	mafCOL = kmerCOLS$mafCOL
	rm(kmerCOLS)

	# PCH by alignment count
	pch_standard = rep(1, length(final_kmer_pos))

	# Subsample for faster plotting (seeded so the figure is reproducible)
	if(length(which(!is.na(ypos)))<1e6){
		s = 1:length(ypos)
	} else {
		set.seed(0)
		s = union(sample(c(which(!is.na(ypos))),1e6),which(ypos>2))
		cat("Subsampling kmers below -log10(p)=2 for faster plotting, plotting",length(s),"kmers","\n")
	}

	xpos = final_kmer_pos[s]
	ypos = ypos[s]
	multialignCOL = multialignCOL[s]
	betaCOL = betaCOL[s]
	mafCOL = mafCOL[s]
	ma = ma[final_kmer_pos_index[s]]
	gene_names = final_kmer_genes[s]
	pch_standard = pch_standard[s]

	gene_conversion = gene_names
	names(gene_conversion) = gene_names

	# How to colour each figure by name (plot_manhattan() reads filecol as a global)
	filecol <<- c("alignCOL","betaCOL","mafCOL","mafCOL")
	# Set which MAF threshold to plot for each figure
	ma_threshold_all = c(minor_allele_threshold, minor_allele_threshold, 0, minor_allele_threshold)
	allCOLS = list(multialignCOL, betaCOL, mafCOL, mafCOL)
	allPCH = list(pch_standard, pch_standard, pch_standard, pch_standard)

	# Figure legend
	legendtext = c("Bonferroni-corrected","significance threshold","")
	legendcol = c("black","white","white")
	legendtext_align = c("Multiple alignments","Single alignments")
	legendtext_MAF = c("MAF < 0.01","0.01 ≤ MAF < 0.05", "MAF ≥ 0.05")
	legendtext_beta = c("β < 0", "β > 0")
	redgrey = c(colour_selection[6], "grey50")
	bluered = c(colour_selection[5], colour_selection[6])
	redbluegreengrey = c(colour_selection[6], colour_selection[5], colour_selection[3],"grey50")
	legendtext = list(c(legendtext, legendtext_align),
						c(legendtext, legendtext_beta),
						c(legendtext, legendtext_MAF),
						c(legendtext, legendtext_MAF))
	legendcol = list(c(legendcol, redgrey),
						c(legendcol, bluered),
						c(legendcol, redbluegreengrey),
						c(legendcol, redbluegreengrey))
	legendpch = c(rep(NA,3),rep(16,3))
	legendlty = c(2,rep(NA,5))

	# Alternate the y-axis limit between max in figure and ylimit of 50 (if above a threshold)
	ylims_options = c(NA, 50)

	for(i in seq_along(filecol)){
		outfilename_prefix = paste0(p[["manhattan_stem"]], "_Manhattan_", filecol[i],"_", macormaf, ma_threshold_all[i])
		for(j in 1:length(ylims_options)){
			if(is.na(ylims_options[j])){
				outfilename = paste0(outfilename_prefix, ".png")
				ylims.i = NULL
				which_genes_to_annotate.i = which(!is.na(gene_names) & ma>=ma_threshold_all[i])
				ma_threshold_pass = which(ma>=ma_threshold_all[i])
				plot.i = TRUE
			} else {
				outfilename = paste0(outfilename_prefix, "_ylim",ylims_options[j],".png")
				ylims.i = c(0, ylims_options[j])
				which_genes_to_annotate.i = which(!is.na(gene_names) & ypos<=max(ylims.i) & ma>=ma_threshold_all[i])
				ma_threshold_pass = which(ma>=ma_threshold_all[i])
				if(max(ypos, na.rm = T)<(max(ylims.i)+(max(ylims.i)/2))) plot.i = FALSE else plot.i = TRUE
			}
			if(plot.i){
				plot_manhattan(outfilename = outfilename, xpos = xpos, ma_threshold_pass = ma_threshold_pass, ypos = ypos, ylims.i = ylims.i, annotateGeneFile = annotateGeneFile, ref = ref, gene_names = gene_names, which_genes_to_annotate.i = which_genes_to_annotate.i, gene_conversion = gene_conversion, allCOLS = allCOLS, allPCH = allPCH, i = i, bonferroni = bonferroni, legendtext = legendtext, legendcol = legendcol, legendpch = legendpch, legendlty = legendlty, beta = as.numeric(assoc[,2]), pheno_type = pheno_type, ref_length = ref_length)
			}
		}
	}
}

# Close-up figures of the top genes: the drawing calls of plot_closeup_alignments()
# (alignmentfunctions.R, r-fixed), with the BLAST results already processed by Python
draw_closeups = function(dd, p) {
	genes = read_fd(dd, "genes")
	if(nrow(genes) == 0) return(invisible())
	ref_gb = read_fd(dd, "features")
	nsamples = as.numeric(p[["nsamples"]])
	bonferroni = as.numeric(p[["bonferroni"]])
	minor_allele_threshold = as.numeric(p[["minor_allele_threshold"]])
	macormaf = p[["macormaf"]]
	kmer_type = p[["kmer_type"]]
	kmer_length = as.numeric(p[["kmer_length"]])
	ref.name = p[["ref_name"]]
	ref_length = as.numeric(p[["ref_length"]])
	output_prefix = p[["output_prefix"]]
	figures_dir = p[["figures_dir"]]
	override_signif = as.logical(p[["override_signif"]])
	correct_only = TRUE

	for(r in seq_len(nrow(genes))) {
		i = genes$index[r]
		genename_i = genes$gene[r]
		gfile = file.path(dd, paste0("gene_", i, "_params.tsv"))
		gp = read_params(gfile)
		need(gp, c("ref_start_i", "ref_end_i", "length_protein", "correct_frame", "length_correct", "strand"), gfile)
		seqs = read_fd(dd, paste0("gene_", i, "_sequences"))
		ref_gene_i = list("ref_start_i" = as.numeric(gp[["ref_start_i"]]), "ref_end_i" = as.numeric(gp[["ref_end_i"]]),
			"ref_gene_i" = seqs$sequence[seqs$name == "region"],
			"length_protein" = as.numeric(gp[["length_protein"]]),
			"all_translations" = seqs$sequence[match(paste0("frame", 1:6), seqs$name)],
			"correct_frame" = as.numeric(gp[["correct_frame"]]), "length_correct" = as.numeric(gp[["length_correct"]]))
		# run_alignment_nplots_nucleotide() reads the strand from column 5 of the gene look-up
		gene_lookup = matrix(c(genename_i, i, ref_gene_i$ref_start_i, ref_gene_i$ref_end_i, gp[["strand"]]), nrow = 1)
		which_kmers_no_result = if(file.exists(file.path(dd, paste0("gene_", i, "_no_result.tsv.gz")))) read_fd(dd, paste0("gene_", i, "_no_result")) else NULL
		cat("Drawing figures for", genename_i, "\n")

		if(kmer_type=="protein"){
			gene_i_results_list = list()
			for(j in 1:6){
				res = read_fd(dd, paste0("gene_", i, "_res_", j))
				gene_i_results_list[[j]] = res
				if((correct_only & j==ref_gene_i$correct_frame) | correct_only==FALSE){
					run_alignment_nplots_protein(ref_gene_i = ref_gene_i, res = res,
									nsamples = nsamples, bonferroni = bonferroni, prefix = output_prefix,
									gene_name = genename_i, j = j, col_lib = col_lib_pro,
									minor_allele_threshold = minor_allele_threshold,
									macormaf = macormaf, output_dir = figures_dir,
									kmer_type = kmer_type, kmer_length = kmer_length,
									ref.name = ref.name, override_signif = override_signif)
				}
			}
			for(j in 1:6){
				if(((correct_only & j==ref_gene_i$correct_frame) | correct_only==FALSE)){
					run_manhattan_single_protein(which_kmers_no_result = which_kmers_no_result,
									res = gene_i_results_list[[j]], ref_gene_i = ref_gene_i,
									kmer_length = kmer_length, prefix = output_prefix,
									gene_name = genename_i, j = j,
									bonferroni = bonferroni,
									ref_gb_full = ref_gb, ref_length = ref_length,
									nsamples = nsamples, kmer_type = kmer_type,
									minor_allele_threshold = minor_allele_threshold, macormaf = macormaf,
									output_dir = figures_dir, ref.name = ref.name)
				}
			}
			run_manhattan_allframes(gene_i_results_list = gene_i_results_list,
								prefix = output_prefix, gene_name = genename_i,
								ref_gene_i = ref_gene_i,
								which_kmers_no_result = which_kmers_no_result,
								bonferroni = bonferroni,
								minor_allele_threshold = minor_allele_threshold, macormaf = macormaf,
								output_dir = figures_dir,
								kmer_type = kmer_type, kmer_length = kmer_length, ref.name = ref.name,
								ref_gb_full = ref_gb)
		} else {
			res = read_fd(dd, paste0("gene_", i, "_res_1"))
			run_alignment_nplots_nucleotide(ref_gene_i = ref_gene_i, res = res,
										nsamples = nsamples, bonferroni = bonferroni, prefix = output_prefix,
										gene_name = genename_i, gene_lookup = gene_lookup,
										wh_genelookup = 1, col_lib = col_lib_nuc,
										minor_allele_threshold = minor_allele_threshold,
										macormaf = macormaf, output_dir = figures_dir,
										kmer_type = kmer_type, kmer_length = kmer_length,
										ref.name = ref.name, override_signif = override_signif)
			run_manhattan_single_nucleotide(which_kmers_no_result = which_kmers_no_result, res = res, ref_gene_i = ref_gene_i, prefix = output_prefix, gene_name = genename_i, bonferroni = bonferroni, ref_gb_full = ref_gb, ref_length = ref_length, kmer_type = kmer_type, kmer_length = kmer_length, nsamples = nsamples, minor_allele_threshold = minor_allele_threshold, macormaf = macormaf, output_dir = figures_dir, ref.name = ref.name)
		}
	}
}

###################################################################################################
## Main
###################################################################################################

main = function(args) {
	usage = "Usage: plot_figures.R --data-dir DIR   (DIR: the figure_data directory written by plotManhattan.py)"
	if(length(args) != 2 || args[1] != "--data-dir") stop(usage)
	dd = args[2]
	if(!dir.exists(dd)) stop("Figure data directory doesn't exist: ", dd)
	pfile = file.path(dd, "params.tsv")
	p = read_params(pfile)
	need(p, c("figures_dir", "output_prefix", "kmer_type", "kmer_length", "ref_name", "ref_length", "macormaf",
		"minor_allele_threshold", "bonferroni", "pheno_type", "nsamples", "override_signif", "manhattan_stem"), pfile)
	if(!dir.exists(p[["figures_dir"]])) stop("Figures directory doesn't exist: ", p[["figures_dir"]])
	expected = readLines(file.path(dd, "expected_figures.txt"), warn = FALSE)
	start.time = Sys.time() - 1

	cat("Drawing QQ plots", "\n")
	draw_qq(dd, p)
	cat("Drawing genome-wide Manhattan plots", "\n")
	draw_genome_manhattan(dd, p)
	cat("Drawing close-up figures for the top genes", "\n")
	draw_closeups(dd, p)
	graphics.off()

	# Every expected figure was drawn, and nothing else
	sz = file.size(expected)
	missing = expected[is.na(sz) | sz == 0]
	drawn = list.files(p[["figures_dir"]], pattern = "\\.png$", full.names = TRUE)
	drawn = drawn[file.mtime(drawn) >= start.time]
	unexpected = setdiff(normalizePath(drawn), normalizePath(expected, mustWork = FALSE))
	if(length(missing)) stop(length(missing), " expected figure(s) not drawn or empty:\n  ", paste(missing, collapse = "\n  "))
	if(length(unexpected)) stop(length(unexpected), " figure(s) drawn that were not expected:\n  ", paste(unexpected, collapse = "\n  "))
	cat("Drew", length(expected), "figures", "\n")
}
