#setwd("E:\\Courses\\edX\\Harvard_PH525.4x\\Week4\\Spielman_et_al_2007")

##35 genes
#Table1_genes = c("UGT2B17","ROBO1","ATP8B1","ARPC4","DPYSL2",
#				"RGS20","FCER2","MEIS2","COPG","HSPB1",
#				"TCF7","HDGFRP3","STXBP2","D21S2056E","MAN2A1",
#				"STS","CTSS","PDE4B","OCIL","TMPRSS3",
#				"NFIL3","LEF1","DNAJB9","ARL7","PTPRO",
#				"ADM","IGF1","PSPHL","CD38","SEC10L1",
#				"WFDC2","P2RY5","SLC2A5","NK4","CLECSF2")
##COPG --> COPG1 --> 217749_at (HG-Focus.na35.annot.csv, from http://www.affymetrix.com/support/technical/byproduct.affx?product=focus)
##D21S2056E --> RRP1 --> 218758_s_at (HG-Focus.na35.annot.csv, from http://www.affymetrix.com/support/technical/byproduct.affx?product=focus)
##OCIL --> CLEC2D (NCBI Gene) --> 220132_s_at (HG-Focus.na35.annot.csv, from http://www.affymetrix.com/support/technical/byproduct.affx?product=focus)
##ARL7 --> ARL4C (NCBI Gene) --> 202207_at (HG-Focus.na35.annot.csv, from http://www.affymetrix.com/support/technical/byproduct.affx?product=focus)
##PSPHL --> PSPH --> 205048_s_at (HG-Focus.na35.annot.csv, from http://www.affymetrix.com/support/technical/byproduct.affx?product=focus)
##SEC10L1 --> EXOC5 --> 218748_s_at (HG-Focus.na35.annot.csv, from http://www.affymetrix.com/support/technical/byproduct.affx?product=focus)
##P2RY5 --> LPAR6 --> 218589_at (HG-Focus.na35.annot.csv, from http://www.affymetrix.com/support/technical/byproduct.affx?product=focus)
##NK4 --> IL32 --> 203828_s_at (HG-Focus.na35.annot.csv, from http://www.affymetrix.com/support/technical/byproduct.affx?product=focus, match with DEF)
##CLECSF2 --> CLEC2B (NCBI Gene) --> 209732_at (HG-Focus.na35.annot.csv, from http://www.affymetrix.com/support/technical/byproduct.affx?product=focus)

##If I found an alias in the Affy file, then I confirmed using NCBI Gene. If I couldn't find an alias in the Affy file, I looked for an official gene symbol in NCBI Gene (and then I looked for a probe)
##PSPHL is an exception in that NCBI gene did not list "PSPH" as the official gene symbol, and I think the "L" often means "like" (which would mean you would expect a different gene).  However, the result in NCBI gene for PSPHL is PSPHP1, and I could not find that in the annotation file.
##The other 5 symbols had an alternative mapping that matched NCBI Gene.
##For the 3 remaining symbols (OCIL, ARL7, CLECSF2), I looked for an alternative symbol in NCBI Gene first (and then I looked for a probe)
Table1_genes = c("UGT2B17","ROBO1","ATP8B1","ARPC4","DPYSL2",
				"RGS20","FCER2","MEIS2","COPG1","HSPB1",
				"TCF7","HDGFRP3","STXBP2","RRP1","MAN2A1",
				"STS","CTSS","PDE4B","CLEC2D","TMPRSS3",
				"NFIL3","LEF1","DNAJB9","ARL4C","PTPRO",
				"ADM","IGF1","PSPH","CD38","EXOC5",
				"WFDC2","LPAR6","SLC2A5","IL32","CLEC2B")
	
#use as a shortcut for annotations https://github.com/genomicsclass/GSE5859Subset
#load first because some variables share names
library(GSE5859Subset)
data(GSE5859Subset)
	
subset_probes = rownames(geneExpression)
print(dim(geneAnnotation))

library(Biobase)
#installed from https://github.com/genomicsclass/GSE5859		
library(GSE5859)
data(GSE5859)

geneExpression = exprs(e)
sampleInfo = pData(e)

print(dim(geneExpression))

matched_probeInfo = geneAnnotation[geneAnnotation$SYMBOL %in% Table1_genes,]
print(dim(matched_probeInfo))
print(length(unique(matched_probeInfo$SYMBOL)))

paper_group = rep("CHB+JPT",nrow(sampleInfo))
paper_group[sampleInfo$ethnicity == "CEU"]="CEU"
paper_group = factor(paper_group, levels = c("CHB+JPT","CEU"))

#defined batch
year = substr(sampleInfo$date,1,4)

lm_ethnicity_pvalue.adj_batch_year = c()
for (i in 1:nrow(geneExpression)){
	X = model.matrix(~paper_group+year)
	y = geneExpression[i,]
	fit = lm(y~X-1)
	results = summary(fit)$coef
	
	lm_ethnicity_pvalue.adj_batch_year[i]=results[2,4]
}#end for (i in 1:nrow(geneExpression))

lm_ethnicity_fdr.adj_batch_year = p.adjust(lm_ethnicity_pvalue.adj_batch_year,"fdr")

#1st 2 principal components (centered expression input for singular value decomposition)
y = geneExpression - rowMeans(geneExpression)

svd_obj = svd(y)
potential_hidden_factors = svd_obj$v[,1:2]

lm_ethnicity_pvalue.adj_PCA = c()
lm_PC1_pvalue = c()
lm_PC2_pvalue = c()
for (i in 1:nrow(geneExpression)){
	X = model.matrix(~paper_group+potential_hidden_factors)
	y = geneExpression[i,]
	fit = lm(y~X-1)
	results = summary(fit)$coef
	
	lm_ethnicity_pvalue.adj_PCA[i]=results[2,4]
	lm_PC1_pvalue[i]=results[3,4]
	lm_PC2_pvalue[i]=results[4,4]
}#end for (i in 1:nrow(geneExpression))

lm_ethnicity_fdr.adj_PCA = p.adjust(lm_ethnicity_pvalue.adj_PCA,"fdr")
lm_PC1_fdr = p.adjust(lm_PC1_pvalue,"fdr")
lm_PC2_fdr = p.adjust(lm_PC2_pvalue,"fdr")

output.summary = data.frame(Gene=matched_probeInfo$SYMBOL, Probe=matched_probeInfo$PROBEID,
							lm_ethnicity_pvalue.adj_batch_year=lm_ethnicity_pvalue.adj_batch_year[match(matched_probeInfo$PROBEID,rownames(geneExpression))],
							lm_ethnicity_fdr.adj_batch_year=lm_ethnicity_fdr.adj_batch_year[match(matched_probeInfo$PROBEID,rownames(geneExpression))],
							lm_ethnicity_pvalue.adj_PCA=lm_ethnicity_pvalue.adj_PCA[match(matched_probeInfo$PROBEID,rownames(geneExpression))],
							lm_ethnicity_fdr.adj_PCA=lm_ethnicity_fdr.adj_PCA[match(matched_probeInfo$PROBEID,rownames(geneExpression))],
							lm_PC1_pvalue=lm_PC1_pvalue[match(matched_probeInfo$PROBEID,rownames(geneExpression))],
							lm_PC1_fdr=lm_PC1_fdr[match(matched_probeInfo$PROBEID,rownames(geneExpression))],
							lm_PC2_pvalue=lm_PC2_pvalue[match(matched_probeInfo$PROBEID,rownames(geneExpression))],
							lm_PC2_fdr=lm_PC2_fdr[match(matched_probeInfo$PROBEID,rownames(geneExpression))])
write.table(output.summary, "lm_summary.txt",quote=F, row.names=F, sep="\t")

print(length(output.summary$lm_ethnicity_fdr.adj_batch_year[output.summary$lm_ethnicity_fdr.adj_batch_year < 0.05])/nrow(output.summary))
print(length(output.summary$lm_ethnicity_fdr.adj_PCA[output.summary$lm_ethnicity_fdr.adj_PCA < 0.05])/nrow(output.summary))
print(length(output.summary$lm_PC1_fdr[output.summary$lm_PC1_fdr < 0.05])/nrow(output.summary))
print(length(output.summary$lm_PC2_fdr[output.summary$lm_PC2_fdr < 0.05])/nrow(output.summary))

#dot plot
groupCol=rep("blue",length(paper_group))
groupCol[paper_group == "CEU"]="red"

png("UGT2B17_dot_plot-edX_provided.png")
geneAnnotation$SYMBOL[is.na(geneAnnotation$SYMBOL)]=""
boxplot(geneExpression[geneAnnotation$SYMBOL == "UGT2B17",]~year,
	ylab="Normalized Expression", col="grey90", outline=F)
points(jitter(as.numeric(as.factor(year))),geneExpression[geneAnnotation$SYMBOL == "UGT2B17",],
		col=groupCol, pch=16)
legend("top",legend=c("CHB+JPT","CEU"),col=c("blue","red"),
		xpd=T, inset=-0.1, ncol=2, pch=16)
dev.off()

#heatmap
library(gplots)

plot_mat = geneExpression[match(matched_probeInfo$PROBEID,rownames(geneExpression)),]
colnames(plot_mat)=paper_group
rownames(plot_mat)=matched_probeInfo$SYMBOL

yearCol = rep("black",ncol(plot_mat))
yearCol[year == 2002]=rainbow(5)[1]
yearCol[year == 2003]=rainbow(5)[2]
yearCol[year == 2004]=rainbow(5)[3]
yearCol[year == 2005]=rainbow(5)[4]
yearCol[year == 2006]=rainbow(5)[5]

png("lm_summary.png")
heatmap.2(plot_mat, trace="none", margins=c(10,10),
			scale="row",
			ColSideColors=yearCol)
legend("topright",legend=2002:2006, col=rainbow(5), pch=15, cex=0.8)
dev.off()
