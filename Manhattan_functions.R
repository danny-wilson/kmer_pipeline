###################################################################################################
## Functions
###################################################################################################



### Get colours
colour_selection = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00")

create_figures_dir = function(dir = NULL, kmer_type = NULL, kmer_length = NULL, alignmenttype = NULL){
	figures_dir = file.path(dir,paste0(kmer_type,"kmer", kmer_length,"_", alignmenttype, "_figures/"))
	if(!dir.exists(figures_dir)) dir.create(figures_dir)
	return(figures_dir)
}



get_genes_to_plot = function(gene_names = NULL, y = NULL, gene_conversion = NULL, ymax = NULL, gene_panel = NULL, ref = NULL, xadjust = NULL, ngenes = 20){
	if(!is.null(gene_names)){
		o = order(as.numeric(y), decreasing = T)
		if(is.null(y)) top_genes = gene_names else top_genes = unique(as.character(gene_names[o]))[1:ngenes]
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

# Function to plot lines behind the points either as a rectangle covering a whole region
# or as a dashed line which is drawn at the  midpoint of the region 
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


get_new_table = function(t1, t2){
	n = intersect(names(t1), names(t2))
	return(list(c(t1[!(names(t1) %in% n)], t2[!(names(t2) %in% n)], t1[n] + t2[n])))
}


get_final_kmer_pos = function(x, index, p_threshold){
	if(length(index[[x]])==0){
		return(x)
	} else {
		t = table(index[[x]]); t = t[which(t>=p_threshold)]
		if(length(t)>0){
			return(rep(x, length(t)))
		} else {
			return(x)
		}
	}
}

get_final_kmer_genes = function(x, p_threshold){
	if(length(x)==0){
		return(NA)
	} else {
		t = table(x); t = t[which(t>=p_threshold)]
		if(length(t)>0){
			return(as.numeric(names(t)))
		} else {
			return(NA)
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


extract_lambda_lognull <- function(datafiles = NULL){
	a <- scan(gzfile(datafiles), what=character(0), sep="\n", quiet = TRUE)
	b <- unlist(strsplit(a[13], " "))
	b <- b[length(b)]
	c <- unlist(strsplit(a[17], " "))
	c <- c[length(c)]
	return(list("lambda" = b, "lognull" = c))
}

get_loglik <- function(LH1 = NULL, lognull = NULL){
	D <- as.numeric(2*(as.numeric(LH1) - as.numeric(lognull)))
	return(D)
}


read_gemma_files = function(input_dir = NULL, prefix = NULL, kmer_type = NULL, kmer_length = NULL, nPatterns = NULL){
	
	files = system(paste0("ls ", input_dir, prefix, "_", kmer_type, kmer_length, "*-*.assoc.txt.gz"), intern = T)
	file.range = gsub(".assoc.txt.gz","",gsub(paste0(input_dir, prefix, "_", kmer_type, kmer_length, "."),"",files,fixed=TRUE),fixed=TRUE)
	file.beg = as.integer(sapply(file.range,function(s)unlist(strsplit(s,"-",fixed=TRUE))[1]))
	file.end = as.integer(sapply(file.range,function(s)unlist(strsplit(s,"-",fixed=TRUE))[2]))
	if(max(file.end)!=nPatterns) stop("Error: max gemma pattern index does not equal total number of patterns","\n")
	gemma_file_index = unlist(apply(cbind(file.beg, file.end), 1, function(x) x[1]:x[2]))
	if(any(is.na(match(1:nPatterns, gemma_file_index)))) stop("Error: not all patterns are present in gemma files", "\n")
	files = files[order(file.beg)]

	for(i in 1:length(files)){
		header = scan(pipe(paste0("zcat ",files[i], " | cut -f2,5,6,10,12")), what = character(0), sep = "\t", nlines = 1, quiet = TRUE)
		gemma.i = scan(pipe(paste0("zcat ",files[i], " | cut -f2,5,6,10,12")), what = character(0), sep = "\t", skip = 1, quiet = TRUE)
		gemma.i = matrix(gemma.i, ncol = 5, byrow = T)
		colnames(gemma.i) = header
		if(i==1){
			assoc = gemma.i
		} else {
			assoc = rbind(assoc, gemma.i)
		}
		# cat("Read in gemma file",i,"of", length(files),"\n")
	}
	rm(gemma.i)

	gemma.log.file = system(paste0("ls ", input_dir, prefix, "_", kmer_type, kmer_length, ".1-*.log.txt.gz"), intern = T)
	l0 = extract_lambda_lognull(gemma.log.file)["lognull"]

	D <- sapply(as.numeric(assoc[,5]), get_loglik, lognull=l0, USE.NAMES=FALSE)
	pvals = pchisq(as.numeric(D), 1, low=F, log.p = TRUE)/-log(10)
	assoc = cbind(assoc, "negLog10" = pvals)
	rm(pvals)
	cat("Number of untested patterns:", length(which(is.na(match(1:nPatterns, as.numeric(assoc[,1]))))), "\n")
	cat("GEMMA range of pvalues:",range(as.numeric(assoc[,4])),"\n")
	cat("GEMMA range of -log10(pvalues):",range(as.numeric(assoc[,6])),"\n")

	# Match gemma results to patterns
	D = D[match(1:nPatterns, as.numeric(assoc[,1]))]
	assoc = assoc[match(1:nPatterns, as.numeric(assoc[,1])),]

	cat("Matched gemma results to all patterns", "\n")

	return(assoc)	
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

write_top_gene_kmers_to_file = function(wh.i = NULL, final_kmer_list = NULL, final_kmer_pos_index = NULL, assoc = NULL, kmerIndex = NULL, mac = NULL, output_file = NULL){
	
	gene.i.kmers = final_kmer_list[final_kmer_pos_index[wh.i]]
	gene.i.beta = as.numeric(assoc[,2])[kmerIndex[final_kmer_pos_index[wh.i]]]
	gene.i.mac = mac[final_kmer_pos_index[wh.i]]
	gene.i.signif = as.numeric(assoc[,6])[kmerIndex[final_kmer_pos_index[wh.i]]]
	# Remove kmers which have not been tested for this phenotype
	gene.i.kmers = gene.i.kmers[which(!is.na(gene.i.signif))]
	gene.i.beta = gene.i.beta[which(!is.na(gene.i.signif))]
	gene.i.mac = gene.i.mac[which(!is.na(gene.i.signif))]
	gene.i.signif = gene.i.signif[which(!is.na(gene.i.signif))]
	# Order from most significant to least
	gene.i.kmers = gene.i.kmers[order(gene.i.signif, decreasing = T)]
	gene.i.beta = gene.i.beta[order(gene.i.signif, decreasing = T)]
	gene.i.mac = gene.i.mac[order(gene.i.signif, decreasing = T)]
	gene.i.signif = gene.i.signif[order(gene.i.signif, decreasing = T)]
write.table(cbind("kmer" = gene.i.kmers, "negLog10" = gene.i.signif, "beta" = gene.i.beta, "mac" = gene.i.mac), file = output_file, row = F, col = T, sep = "\t", quote = F)
	
}


top20genes = function(gene_names = NULL, ma = NULL, minor_allele_threshold = NULL, ypos = NULL, macormaf = NULL, output_dir = NULL, prefix = NULL, min_count = NULL, ident_threshold = NULL, kmer_type = NULL, kmer_length = NULL, ref.name = NULL){
	
	gene_conversion = gene_names
	names(gene_conversion) = gene_names
	
	# Find the top 20 genes (for kmers with MAC/MAF above the threshold) by p-value and store the gene name plus the most significant p-value per gene
	which_genes_to_annotate = which(!is.na(gene_names) & ma>=minor_allele_threshold)
	top20genes = unique(gene_names[which_genes_to_annotate][order(ypos[which_genes_to_annotate], decreasing = T)])[1:20]
	top20genespvals = c()
	for(i in 1:length(top20genes)){
		top20genespvals[i] = max(ypos[which_genes_to_annotate][which(gene_names[which_genes_to_annotate]==top20genes[i])], na.rm = T)
	}
	cat(paste0("Top 20 genes above the ",macormaf," threshold ", minor_allele_threshold,":"),"\n")
	cat("Gene -log10pvalue","\n")
	printtop20 = apply(cbind(top20genes, top20genespvals), 1, function(x) cat(x[1],x[2], "\n"))

	outprefix = paste0(output_dir, prefix, "_", kmer_type, kmer_length, "_", ref.name, "_top20genes_toppvals_", macormaf, "_", minor_allele_threshold)
	if(!is.null(min_count)){
		outfile = paste0(outprefix, "_nucmerAlign_alignIdent_", ident_threshold,"_alignPosMinCount_", min_count, ".txt")
	} else {
		outfile = paste0(outprefix, "_bowtie2mapping.txt")
	}

	write.table(cbind(top20genes, top20genespvals), file = outfile, row = F, col = F, sep = "\t", quote = F)
	
	
}


plot_manhattan = function(outfilename = NULL, xpos = NULL, ma_threshold_pass = NULL, ypos = NULL, ylims.i = NULL, annotateGeneFile = NULL, ref = NULL, which_genes_to_annotate.i = NULL, allCOLS = NULL, allPCH = NULL, i = NULL, bonferroni = NULL, legendtext = NULL, legendcol = NULL, legendpch = NULL, legendlty = NULL, beta = NULL, gene_names = NULL, gene_conversion = NULL, pheno_type = NULL){
	
			
	
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
	axis(1, cex.axis = 0.8, at = c(0,1,2,3,4,5)*1e6, labels = c("0","1","2","3","4","5"))
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


get_pheno_type = function(pheno){
	
	if(length(table(pheno))==2) return("binary") else return("continuous")
	
}


