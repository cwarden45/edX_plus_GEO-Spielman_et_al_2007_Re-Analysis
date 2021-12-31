#setwd("E:\\Courses\\edX\\Harvard_PH525.4x\\Week4\\Spielman_et_al_2007\\Re-Analysis")

celpath="Raw_Data"
metadata.in = "GEO_sample_infomation-v2.txt"
metadata.out = "GEO_sample_infomation-v3.txt"
expr.out = "GEO_RMA_expr.txt"
edX_analysis.in = "../lm_summary.txt"
GEO_analysis.out = "cwarden_GEO-lm_summary.txt"
GEO_analysis.plot1 = "cwarden_GEO-lm_summary.png"
GEO_analysis.plot2 = "UGT2B17_dot_plot-GEO_unique.png"

library(affy)

affy_data = ReadAffy(celfile.path=celpath)
sampleNames = colnames(affy_data)
probes = rownames(affy_data)

scan_date = protocolData(affy_data)$ScanDate
#Get scan date from file: https://www.biostars.org/p/61822/


##use to manually create a metadata file with sample information
#write.table(data.frame(Samples=sampleNames),"GEO_sample_infomation.txt", quote=F, sep="\t", row.names=F)

metadata.table = read.table(metadata.in, head=T, sep="\t")
metadata.table = data.frame(metadata.table, Scan_Date=scan_date)
write.table(metadata.table, metadata.out, quote=F, sep="\t", row.names=F)

eset = rma(affy_data)
write.exprs(eset, file=expr.out)

#degregation plot (from https://www.bioconductor.org/packages/release/bioc/vignettes/affy/inst/doc/affy.pdf)
degregation_stats = AffyRNAdeg(affy_data)
png("QC_RNAdegradation.png")
plotAffyRNAdeg(degregation_stats)
dev.off()

#use first listed sample, if there are duplicates
unique_samples = unique(metadata.table$Subject)
print(dim(metadata.table))
metadata_filtered = metadata.table[match(unique_samples, metadata.table$Subject),]
print(dim(metadata_filtered))

print(dim(eset))
expr_filtered = eset[,match(metadata_filtered$Samples,colnames(eset))]
print(dim(expr_filtered))

##copy and modify similar code as used for edX expression

paper_group = rep("CHB+JPT",nrow(metadata_filtered))
paper_group[metadata_filtered$Alt_Pop_ID == "CEU"]="CEU"
paper_group = factor(paper_group, levels = c("CHB+JPT","CEU"))

year = paste("20",substr(metadata_filtered$Scan_Date,7,8),sep="")
year = as.numeric(year)

geneExpression = exprs(expr_filtered)

lm_ethnicity_pvalue.adj_batch_year = c()
for (i in 1:nrow(geneExpression)){
	X = model.matrix(~paper_group+year)
	y = as.numeric(geneExpression[i,])
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


Table1.edX = read.table(edX_analysis.in, head=T, sep="\t")
output.summary = data.frame(Gene=Table1.edX$Gene, Probe=Table1.edX$Probe,
							lm_ethnicity_pvalue.adj_batch_year=lm_ethnicity_pvalue.adj_batch_year[match(Table1.edX$Probe,rownames(expr_filtered))],
							lm_ethnicity_fdr.adj_batch_year=lm_ethnicity_fdr.adj_batch_year[match(Table1.edX$Probe,rownames(expr_filtered))],
							lm_ethnicity_pvalue.adj_PCA=lm_ethnicity_pvalue.adj_PCA[match(Table1.edX$Probe,rownames(expr_filtered))],
							lm_ethnicity_fdr.adj_PCA=lm_ethnicity_fdr.adj_PCA[match(Table1.edX$Probe,rownames(expr_filtered))],
							lm_PC1_pvalue=lm_PC1_pvalue[match(Table1.edX$Probe,rownames(expr_filtered))],
							lm_PC1_fdr=lm_PC1_fdr[match(Table1.edX$Probe,rownames(expr_filtered))],
							lm_PC2_pvalue=lm_PC2_pvalue[match(Table1.edX$Probe,rownames(expr_filtered))],
							lm_PC2_fdr=lm_PC2_fdr[match(Table1.edX$Probe,rownames(expr_filtered))])
write.table(output.summary, GEO_analysis.out,quote=F, row.names=F, sep="\t")

print(length(output.summary$lm_ethnicity_fdr.adj_batch_year[output.summary$lm_ethnicity_fdr.adj_batch_year < 0.05])/nrow(output.summary))
print(length(output.summary$lm_ethnicity_fdr.adj_PCA[output.summary$lm_ethnicity_fdr.adj_PCA < 0.05])/nrow(output.summary))
print(length(output.summary$lm_PC1_fdr[output.summary$lm_PC1_fdr < 0.05])/nrow(output.summary))
print(length(output.summary$lm_PC2_fdr[output.summary$lm_PC2_fdr < 0.05])/nrow(output.summary))

plot_mat = geneExpression[match(Table1.edX$Probe,rownames(expr_filtered)),]
colnames(plot_mat)=paper_group
rownames(plot_mat)=Table1.edX$Gene

#dot plot
groupCol=rep("blue",length(paper_group))
groupCol[paper_group == "CEU"]="red"

png(GEO_analysis.plot2)
boxplot(plot_mat[Table1.edX$Gene == "UGT2B17",]~year,
	ylab="Normalized Expression", col="grey90", outline=F)
points(jitter(as.numeric(as.factor(year))),plot_mat[Table1.edX$Gene == "UGT2B17",],
		col=groupCol, pch=16)
legend("top",legend=c("CHB+JPT","CEU"),col=c("blue","red"),
		xpd=T, inset=-0.1, ncol=2, pch=16)
dev.off()

#heatmap
library(gplots)

yearCol = rep("black",ncol(plot_mat))
yearCol[year == 2002]=rainbow(5)[1]
yearCol[year == 2003]=rainbow(5)[2]
yearCol[year == 2004]=rainbow(5)[3]
yearCol[year == 2005]=rainbow(5)[4]
yearCol[year == 2006]=rainbow(5)[5]

png(GEO_analysis.plot1)
heatmap.2(plot_mat, trace="none", margins=c(10,10),
			scale="row",
			ColSideColors=yearCol)
legend("topright",legend=2002:2006, col=rainbow(5), pch=15, cex=0.8)
dev.off()
