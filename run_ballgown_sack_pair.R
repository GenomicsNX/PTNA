
library(ballgown)
# library(RSkittleBrewer)
library(genefilter)
library(dplyr)
library(devtools)
library(ggplot2)
library(ggfortify)
library("biomaRt")
library(edgeR)

setwd("/data/NHLBI_BCB/Sack_Kim-Han_RNAseq_v2/07-ballgown")

## Read phenotype sample data
pheno_data = read.csv(file="/data/NHLBI_BCB/Sack_Kim-Han_RNAseq_v2/07-ballgown/phenodata.csv", header=T, row.names=1)

## Read sample expression data (created by stringtie) into ballgown, fpkm, counts, trasncript, gene level, fpkm, ballgown object
# bg_samples = ballgown(dataDir = "/data/NHLBI_BCB/Sack_Kim-Han_RNAseq_v2/05-StringTie", samplePattern = "S")
## save the ballgown object created above 
# save(bg_samples, file="bg_samples.rda")
load("/data/NHLBI_BCB/Sack_Kim-Han_RNAseq_v2/07-ballgown/bg_samples.rda")

## Read sample expression data (created by stringtie) into ballgown, fpkm, counts, trasncript, gene level, counts, count matrix --> COUNT
#bg_samples = read.csv(file="/data/NHLBI_BCB/Sack_Kim-Han_RNAseq_v2/06-StringTie_count/sack.gencode.v25.transcript.count.txt", header=T, row.names=1)

#bg_samples = read.csv(file="/data/NHLBI_BCB/Sack_Kim-Han_RNAseq_v2/06-StringTie_count/sack.gencode.v25.gene.count.txt", header=T, row.names=1)

## Build phenotype matrix in the same order as bg_samples, ballgown object
samples = data.frame(sampleNames(bg_samples))
colnames(samples) = "SampleName"
pheno_data = merge(samples, pheno_data, by="SampleName")

# Run on ballgown object
pData(bg_samples) = data.frame(id=sampleNames(bg_samples), group=pheno_data$group, group_text=pheno_data$group_text, pair=pheno_data$pair ,  id_group=pheno_data$SampleName_group)


#mkdir -p pair/01-Baseline_Fasting/03-Gene_fpkm #0/1
#mkdir -p pair/02-Baseline_Refeeding/03-Gene_fpkm #0/2
#mkdir -p pair/03-Fasting_Refeeding/03-Gene_fpkm #1/2
  

## Filter by condition (change comparison groups below to whatever you are comparing in your pheno_data matrix) 
setwd("/data/NHLBI_BCB/Sack_Kim-Han_RNAseq_v2/07-ballgown/pair/01-Baseline_Fasting/03-Gene_fpkm")
group1 = "0"
group2 = "1"

setwd("/data/NHLBI_BCB/Sack_Kim-Han_RNAseq_v2/07-ballgown/pair/02-Baseline_Refeeding/03-Gene_fpkm")
group1 = "0"
group2 = "2"


setwd("/data/NHLBI_BCB/Sack_Kim-Han_RNAseq_v2/07-ballgown/pair/03-Fasting_Refeeding/03-Gene_fpkm")
group1 = "1"
group2 = "2"





pheno_data_comparison = subset(pheno_data, pheno_data$group %in% c(group1,group2))
pheno_data_comparison = pheno_data_comparison[order(pheno_data_comparison$group),]

# count matrix --> COUNT
#bg_samples_filt = subset(bg_samples, select=as.vector(pheno_data_comparison$SampleName))

# ballgown object  
bg_samples_filt = subset(bg_samples, "pData(bg_samples)$id %in% pheno_data_comparison$SampleName", genomesubset=FALSE)

## Filter low abundance transcripts, genes, fpkm, ballgown object
# bg_samples_filt = subset(bg_samples_filt, "rowVars(texpr(bg_samples_filt)) > 1", genomesubset=TRUE)
# bg_samples_filt_all_zeros_fpkm = subset(bg_samples_filt, "rowSums((texpr(bg_samples_filt)[,]>0))>=1", genomesubset=TRUE)
# bg_samples_filt = bg_samples_filt_all_zeros_fpkm

# MPMPMP transcript
#bg_samples_filt_all_less_than_one_fpkm = subset(bg_samples_filt, "rowSums((texpr(bg_samples_filt)[,]>1))>=1", genomesubset=TRUE)
#bg_samples_filt = bg_samples_filt_all_less_than_one_fpkm
# bg_samples_filt_all_less_than_one_half_samples_fpkm = subset(bg_samples_filt, "rowSums((texpr(bg_samples_filt)[,]>1))>=4", genomesubset=TRUE)

#run for genes MPMPMP
gexpr = subset(gexpr(bg_samples_filt), select=paste("FPKM.",as.vector(pheno_data_comparison$SampleName),sep="")) #run for genes
select = rowSums((gexpr)[,]>1)>=1 #run for genes
keep = which(select) #run for genes
bg_samples_filt_all_less_than_one_half_samples_fpkm = gexpr[keep,] #run for genes
bg_samples_filt = bg_samples_filt_all_less_than_one_half_samples_fpkm #run for genes

## Filter low abundance transcripts, genes, counts, count matrix on row counts
#select = rowSums((bg_samples_filt)[,]>5)>=1
#keep = which(select)
#bg_samples_filt = bg_samples_filt[keep,]

## Filter low abundance transcripts, genes, counts, count matrix on counts per million
# select = rowSums(cpm(bg_samples_filt)>1)>=1
# keep = which(select)
# bg_samples_filt_all_less_than_one_half_samples_cpm = bg_samples_filt[keep,] 
# bg_samples_filt = bg_samples_filt_all_less_than_one_half_samples_cpm 

## FPKM
## Differential Expression Analysis, transcripts, genes, fpkm
#feature='transcript'
feature='gene'

## create design matrices: ONLY RUN IF YOU HAVE COVARIATES

#run for transcript
# stat_results_transcripts = stattest(bg_samples_filt, feature=feature, meas='FPKM', covariate='group', getFC=TRUE, adjustvars=c('pair'))
#run for genes
stat_results_transcripts = stattest(gowntable=bg_samples_filt, feature=feature, meas='FPKM', pData=pheno_data_comparison, covariate='group', getFC=TRUE, adjustvars=c('pair')) 



# texpr = subset(texpr(bg_samples_filt), select=paste("FPKM.",as.vector(pheno_data_comparison$SampleName),sep="")) #sorted by comparison groups
# stat_results_transcripts = data.frame(geneName=ballgown::geneNames(bg_samples_filt), geneID=ballgown::geneIDs(bg_samples_filt), transcriptNames=ballgown::transcriptNames(bg_samples_filt), stat_results_transcripts, log2fc=log2(stat_results_transcripts$fc), texpr)
# trans_info = read.csv('/data/NHLBI_BCB/Sack_Kim-Han_RNA_seq/07-ballgown/transcript_annotations_genecodev25.csv')
# stat_results_transcripts = merge(trans_info, stat_results_transcripts, by.x="t_id", by.y="id", all.x=FALSE, all.y=TRUE)
# lncrna_trans_info = read.table('/data/NHLBI_BCB/Fayaz/02-RNA-seq-transcript-level/data_source/GENCODE_lncRNA_Mouse_Human/gencode.v25.long_noncoding_RNAs.transcript_ids.txt', header=T)
# lncrna_trans_info = data.frame(lncrna_transid = lncrna_trans_info$lncrna_transid, lncrna_transid = lncrna_trans_info$lncrna_transid)
# stat_results_transcripts = merge(x=stat_results_transcripts, y=lncrna_trans_info, by.x="t_name", by.y="lncrna_transid", all.x=TRUE)

###custom model
# mod = model.matrix(~ pData(bg_samples_filt)$pair + pData(bg_samples_filt)$group)
# mod0 = model.matrix(~ pData(bg_samples_filt)$group)
# stat_results_transcripts_for_fc = stattest(bg_samples_filt, feature=feature, meas='FPKM', covariate='group', getFC=TRUE, libadjust = FALSE)
# stat_results_transcripts = stattest(bg_samples_filt, feature=feature, meas='FPKM', mod0=mod0, mod=mod)
# texpr = subset(texpr(bg_samples_filt), select=paste("FPKM.",as.vector(pheno_data_comparison$SampleName),sep="")) #sorted by comparison groups
# stat_results_transcripts = data.frame(geneName=ballgown::geneNames(bg_samples_filt), geneID=ballgown::geneIDs(bg_samples_filt), transcriptNames=ballgown::transcriptNames(bg_samples_filt), stat_results_transcripts, texpr)
# stat_results_transcripts = data.frame(stat_results_transcripts, fc=stat_results_transcripts_for_fc$fc, log2fc=log2(stat_results_transcripts_for_fc$fc))

# stat_results_transcripts = stattest(gowntable=bg_samples_filt, feature=feature, meas='FPKM', pData=pheno_data_comparison, covariate='group', getFC=TRUE, adjustvars=c('pair')) #run for genes
# stat_results_transcripts = data.frame(stat_results_transcripts, log2fc=log2(stat_results_transcripts$fc), ensembl_gene_id=sub("[.].*","",stat_results_transcripts$id), bg_samples_filt) #run for genes
# gene_info = read.csv('/data/NHLBI_BCB/Sack_Kim-Han_RNA_seq/07-ballgown/gene_annotations_genecodev25.csv')
# stat_results_transcripts = merge(gene_info, stat_results_transcripts, by.x="gene_id", by.y="id", all.x=FALSE, all.y = TRUE)
# lncrna_gene_info = read.table('/data/NHLBI_BCB/Fayaz/02-RNA-seq-transcript-level/data_source/GENCODE_lncRNA_Mouse_Human/gencode.v25.long_noncoding_RNAs.gene_ids.txt', header=T)
# lncrna_gene_info = data.frame(lncrna_geneid = lncrna_gene_info$lncrna_geneid, lncrna_geneid = lncrna_gene_info$lncrna_geneid)
# stat_results_transcripts = merge(x=stat_results_transcripts, y=lncrna_gene_info, by.x="gene_id", by.y="lncrna_geneid", all.x=TRUE)

###custom model
# mod = model.matrix(~ pheno_data_comparison$pair + pheno_data_comparison$group) #run for genes
# mod0 = model.matrix(~ pheno_data_comparison$group) #run for genes
# stat_results_transcripts_for_fc = stattest(gowntable=bg_samples_filt, feature=feature, meas='FPKM', pData=pheno_data_comparison, covariate='group', getFC=TRUE, libadjust = FALSE) #run for genes
# stat_results_transcripts = stattest(gowntable=bg_samples_filt, feature=feature, meas='FPKM', pData=pheno_data_comparison, mod0=mod0, mod=mod) #run for genes
# stat_results_transcripts = data.frame(stat_results_transcripts, fc=stat_results_transcripts_for_fc$fc, log2fc=log2(stat_results_transcripts_for_fc$fc), ensembl_gene_id=sub("[.].*","",stat_results_transcripts$id), bg_samples_filt) #run for genes

## RUN THIS WITHOUT COVARIATES
#stat_results_transcripts = stattest(bg_samples_filt, feature=feature, meas='FPKM', covariate='group', getFC=TRUE) #run for transcripts



## Add gene name, gene ids, transcript ids, log2fc, FPKM expression values for comparison group
#texpr = subset(texpr(bg_samples_filt), select=paste("FPKM.",as.vector(pheno_data_comparison$SampleName),sep="")) #run for transcripts

# transcript
#colnames(texpr) = paste("FPKM.",as.vector(pheno_data_comparison$SampleName_group),sep="") 
#stat_results_transcripts = data.frame(geneName=ballgown::geneNames(bg_samples_filt), geneID=ballgown::geneIDs(bg_samples_filt), transcriptNames=ballgown::transcriptNames(bg_samples_filt), stat_results_transcripts, log2fc=log2(stat_results_transcripts$fc), texpr) #run for transcripts
#trans_info = read.csv('/data/NHLBI_BCB/Sack_Kim-Han_RNA_seq/07-ballgown/transcript_annotations_genecodev25.csv')
#stat_results_transcripts = merge(trans_info, stat_results_transcripts, by.x="t_name", by.y="transcriptNames", all.x=FALSE, all.y=TRUE)
#lncrna_trans_info = read.table('/data/NHLBI_BCB/Fayaz/02-RNA-seq-transcript-level/data_source/GENCODE_lncRNA_Mouse_Human/gencode.v25.long_noncoding_RNAs.transcript_ids.txt', header=T)
#lncrna_trans_info = data.frame(lncrna_transid = lncrna_trans_info$lncrna_transid, lncrna_transid = lncrna_trans_info$lncrna_transid)
#stat_results_transcripts = merge(x=stat_results_transcripts, y=lncrna_trans_info, by.x="t_name", by.y="lncrna_transid", all.x=TRUE)

# MPMPMP gene
colnames(gexpr) = paste("FPKM.",as.vector(pheno_data_comparison$SampleName_group),sep="") 
stat_results_transcripts = data.frame(stat_results_transcripts, log2fc=log2(stat_results_transcripts$fc), bg_samples_filt) #run for genes
gene_info = read.csv('/data/NHLBI_BCB/Sack_Kim-Han_RNA_seq/07-ballgown/gene_annotations_genecodev25.csv')
stat_results_transcripts = merge(gene_info, stat_results_transcripts, by.x="gene_id", by.y="id", all.x=FALSE, all.y = TRUE)
lncrna_gene_info = read.table('/data/NHLBI_BCB/Fayaz/02-RNA-seq-transcript-level/data_source/GENCODE_lncRNA_Mouse_Human/gencode.v25.long_noncoding_RNAs.gene_ids.txt', header=T)
lncrna_gene_info = data.frame(lncrna_geneid = lncrna_gene_info$lncrna_geneid, lncrna_geneid = lncrna_gene_info$lncrna_geneid)
stat_results_transcripts = merge(x=stat_results_transcripts, y=lncrna_gene_info, by.x="gene_id", by.y="lncrna_geneid", all.x=TRUE)

#annotating to ensembl
# listMarts()
ensembl = useMart("ensembl")
# listDatasets(ensembl)
ensembl = useMart("ensembl", dataset="hsapiens_gene_ensembl")
# listAttributes(ensembl)
#stat_results_transcripts_annotated = getBM(attributes=c("ensembl_gene_id","external_gene_name"), filters="ensembl_gene_id", values=sub("[.].*","",stat_results_transcripts$id), mart=ensembl)
#stat_results_transcripts = merge(stat_results_transcripts, stat_results_transcripts_annotated, by="ensembl_gene_id", all=TRUE)
head(stat_results_transcripts)
## Sort results from smallest p-value
stat_results_transcripts = arrange(stat_results_transcripts, pval)
write.csv(stat_results_transcripts, paste("stat_results_",feature,"_fpkm_",group1,"_",group2,".csv", sep=""), row.names=FALSE)

## Differential Expression Analysis (LM), transcripts, genes, fpkm
#texpr.fpkm = subset(texpr(bg_samples_filt), select=paste("FPKM.",as.vector(pheno_data_comparison$SampleName),sep=""))
#row.names(texpr.fpkm) = ballgown::transcriptNames(bg_samples_filt)
#lm_stats = apply(texpr.fpkm, 1, function(x) linear_model(log(x+1, base=2), grp=pheno_data_comparison))
#lm_stats = data.frame(t(lm_stats))
#colnames(lm_stats) = c("foldchange", "pvalue")
#lm_stats = lm_stats[order(lm_stats$pvalue), ]
#write.csv(lm_stats,file=paste("lm_stat_results_",feature,"_fpkm_",group1,"_",group2,".csv", sep=""), row.names=TRUE, quote=FALSE)

#linear_model = function(x,grp=NULL){
#  grp = model.matrix(~grp[,"group"] + grp[,"pair"])
#  lm = lm(x ~ grp)
#  beta_coefficient = tryCatch(coefficients(summary(lm))[2,1], error=function(e) NA)
#  if (is.na(beta_coefficient[1])==F) {foldchange = 2^beta_coefficient} else {foldchange="NA"}
# pval = tryCatch(coefficients(summary(lm))[2,4], error=function(e) NA)
#  if (is.na(pval[1])==F) {pval = pval} else {pval="NA"}
#  stats = cbind(foldchange, pval)
#}

# pca plot - using results with p-value < 0.05 or any other subset based on significant DE for e.g. fc, log2fc, transcripts, genes, fpkm
# stat_results_transcripts_subset = subset(stat_results_transcripts, subset=stat_results_transcripts$log2fc < -1 | stat_results_transcripts$log2fc > 1)
stat_results_transcripts_subset = subset(stat_results_transcripts, subset=stat_results_transcripts$pval < 0.05)
stat_results_transcripts_subset_fpkm = subset(stat_results_transcripts_subset, select = grep("FPKM", names(stat_results_transcripts_subset)))

# stat_results_transcripts_subset_fpkm = subset(stat_results_transcripts, select = grep("FPKM", names(stat_results_transcripts))) #all transcripts or all genes after filtering low FPKM transcripts

# pca = prcomp(t(stat_results_transcripts_subset_fpkm), scale = F)
pca = prcomp(t(log2(stat_results_transcripts_subset_fpkm+1)), scale = F)
pdf(paste("pca_plot_p_less_than_0.05_pairs_plot_first_five_PCs_",feature,"_",group1,"_",group2,".pdf", sep=""))
pairs(pca$x[,1:5], main=paste("scatter plot matrix of Principal Components_first_five_PCs_",feature,"_",group1,"_",group2, sep=""), pch=21, bg=c('green','red')[as.vector(unclass(as.factor(as.vector(pheno_data_comparison$group))))], cex.main=0.75)
dev.off()
pdf(paste("pca_plot_p_less_than_0.05_pc1_vs_pc2_",feature,"_",group1,"_",group2,".pdf", sep=""))
# Add extra space to right of plot area; change clipping to figure
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
plot(pca$x[,1], pca$x[,2], main=paste("scatter plot matrix of Principal Components_PC1.vs.PC2_",feature,"_",group1,"_",group2, sep=""), pch=16, col=c('green','red')[as.vector(unclass(as.factor(as.vector(pheno_data_comparison$group))))], cex=1.20, cex.main=0.90, xlab="PC1", ylab="PC2")
text(pca$x[,1], pca$x[,2], labels=pheno_data_comparison$SampleName, pos=3, cex=0.70)
legend("topright", inset=c(-0.2,0), legend=c(group1,group2), pch=c(16,16), col=levels(as.factor(c('green','red')[as.vector(unclass(as.factor(as.vector(pheno_data_comparison$group))))])), title="Group", cex=0.70)
dev.off()

pdf(paste("pca_plot_all_transcripts_pc2_vs_pc3_",feature,"_",group1,"_",group2,".pdf", sep=""))
# Add extra space to right of plot area; change clipping to figure
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
plot(pca$x[,2], pca$x[,3], main=paste("scatter plot matrix of Principal Components_PC2.vs.PC3_",feature,"_",group1,"_",group2, sep=""), pch=16, col=c('green','red')[as.vector(unclass(as.factor(as.vector(pheno_data_comparison$group))))], cex=1.20, cex.main=0.90, xlab="PC2", ylab="PC3")
text(pca$x[,2], pca$x[,3], labels=pheno_data_comparison$SampleName, pos=3, cex=0.70)
legend("topright", inset=c(-0.2,0), legend=c(group1,group2), pch=c(16,16), col=levels(as.factor(c('green','red')[as.vector(unclass(as.factor(as.vector(pheno_data_comparison$group))))])), title="Group", cex=0.70)
dev.off()

## Plotting abundance distribution, after applying filter, transcripts, genes, fpkm
#texpr_comparison = subset(texpr(bg_samples_filt), select=paste("FPKM.",as.vector(pheno_data_comparison$SampleName),sep=""))
#fpkm = texpr_comparison
gexpr_comparison = subset(stat_results_transcripts, select = grep("FPKM", names(stat_results_transcripts))) #run for genes
fpkm = gexpr_comparison #run for genes
colnames(fpkm) = pheno_data_comparison$SampleName_group
fpkm = log2(fpkm+1)
pdf(file=paste("fpkm_distribution_",feature,"_",group1,"_",group2,".pdf", sep=""))
boxplot(fpkm, col=c('green','red')[as.vector(unclass(as.factor(as.vector(pheno_data_comparison$group))))], las=2, ylab='log2(FPKM+1)', main=paste("FPKM distribution:",feature,"_",group1,"_",group2), xlab="sample", cex.main=0.80, cex.axis=0.50)
dev.off()

## Plotting volcano plot, Differential Expression Analysis-qval, transcripts, genes, fpkm
data = read.csv(file=paste("stat_results_",feature,"_fpkm_",group1,"_",group2,".csv", sep=""), header=T, row.names=1)
q.trans = -1*log(data$qval, base = 10)
logfold = log(data$fc, base = 2)
pdf(file=paste("volcano_qval_",feature,"_",group1,"_",group2,".pdf", sep=""))
plot(logfold,q.trans,type="n",ylab="-1*log10(q-value)",xlab="log2(fold change)",main=paste("FPKM:",feature,"_",group1,"_",group2), cex.main=0.80, xlim=range(logfold), ylim=range(q.trans))
points(logfold,q.trans,col="black",cex=0.65)
points(logfold[(q.trans>1.3&logfold>1.0)],q.trans[(q.trans>1.3&logfold>1.0)],col="red",pch=16,cex=0.65)
points(logfold[(q.trans>1.3&logfold<(-1.0))],q.trans[(q.trans>1.3&logfold<(-1.0))],col="green",pch=16,cex=0.65)
# text(logfold[(q.trans>1.3&logfold>1.0)],q.trans[(q.trans>1.3&logfold>1.0)],labels=as.character(data$geneName[(q.trans>1.3&logfold>1.0)]), cex=0.65)
# text(logfold[(q.trans>1.3&logfold<(-1.0))],q.trans[(q.trans>1.3&logfold<(-1.0))],labels=as.character(data$geneName[(q.trans>1.3&logfold<(-1.0))]), cex=0.65)
abline(h=1.3)
abline(v=-1.0)
abline(v=1.0)
dev.off()

## Plotting volcano plot, Differential Expression Analysis-pval, transcripts, genes, fpkm
p.trans = -1*log(data$pval, base = 10)
logfold = log(data$fc, base = 2)
pdf(file=paste("volcano_pval_",feature,"_",group1,"_",group2,".v2.pdf", sep=""))
plot(logfold,p.trans,type="n",ylab="-1*log10(p-value)",xlab="log2(fold change)",main=paste("FPKM:",feature,"_",group1,"_",group2), cex.main=0.80, xlim=range(logfold), ylim=range(p.trans))
points(logfold,p.trans,col="black",cex=0.65)
# points(logfold[(p.trans>1.3&logfold>1.0)],p.trans[(p.trans>1.3&logfold>1.0)],col="red",pch=16,cex=0.65)
# points(logfold[(p.trans>1.3&logfold<(-1.0))],p.trans[(p.trans>1.3&logfold<(-1.0))],col="green",pch=16,cex=0.65)
points(logfold[(p.trans>1.3&logfold>0.5)],p.trans[(p.trans>1.3&logfold>0.5)],col="red",pch=16,cex=0.65)
points(logfold[(p.trans>1.3&logfold<(-0.5))],p.trans[(p.trans>1.3&logfold<(-0.5))],col="green",pch=16,cex=0.65)
# text(logfold[(p.trans>1.3&logfold>1.0)],p.trans[(p.trans>1.3&logfold>1.0)],labels=as.character(data$geneName[(p.trans>1.3&logfold>1.0)]), cex=0.65)
# text(logfold[(p.trans>1.3&logfold<(-1.0))],p.trans[(p.trans>1.3&logfold<(-1.0))],labels=as.character(data$geneName[(p.trans>1.3&logfold<(-1.0))]), cex=0.65)
abline(h=1.3)
# abline(v=-1.0)
# abline(v=1.0)
abline(v=-0.5)
abline(v=0.5)
dev.off()

pdf(file=paste("qval_distribution_",feature,"_",group1,"_",group2,".pdf", sep=""))
hist(data$qval, col="grey", border="white", xlab="qval", ylab="", main=paste("FPKM:",feature,"_",group1,"_",group2))
dev.off()
pdf(file=paste("pval_distribution_",feature,"_",group1,"_",group2,".pdf", sep=""))
hist(data$pval, col="grey", border="white", xlab="pval", ylab="", main=paste("FPKM:",feature,"_",group1,"_",group2))
dev.off()

## COUNTS
## Differential Expression Analysis, transcripts, genes, counts

## create design matrices: ONLY RUN IF YOU HAVE COVARIATES
# mod = model.matrix(~ pData(bg_samples_filt)$age + pData(bg_samples_filt)$group)
# mod0 = model.matrix(~ pData(bg_samples_filt)$group)
# stat_results_transcripts_for_fc = stattest(bg_samples_filt, feature='transcript', meas='FPKM', covariate='group', getFC=TRUE)
# stat_results_transcripts = stattest(bg_samples_filt, feature='transcript', meas='FPKM', mod0=mod0, mod=mod)
# texpr_rh = subset(texpr(bg_samples_filt), select=paste("FPKM.",as.vector(pheno_data_comparison$SampleName),sep="")) #sorted by comparison groups
# stat_results_transcripts = data.frame(geneName=ballgown::geneNames(bg_samples_filt), geneID=ballgown::geneIDs(bg_samples_filt), transcriptNames=ballgown::transcriptNames(bg_samples_filt), stat_results_transcripts, texpr_rh)
# stat_results_transcripts = data.frame(stat_results_transcripts, fc=stat_results_transcripts_for_fc$fc, log2fc=log2(stat_results_transcripts_for_fc$fc))

## create design matrices: ONLY RUN IF YOU HAVE COVARIATES
##custom model
# mod = model.matrix(~ pheno_data_comparison$pair + pheno_data_comparison$group)
# mod0 = model.matrix(~ pheno_data_comparison$group)
# stat_results_transcripts_for_fc = stattest(gowntable=bg_samples_filt, feature=feature, pData=pheno_data_comparison, covariate='group', getFC=TRUE, libadjust = FALSE)
# stat_results_transcripts = stattest(gowntable=bg_samples_filt, pData=pheno_data_comparison, feature=feature, mod0=mod0, mod=mod)
# stat_results_transcripts = data.frame(stat_results_transcripts, fc=stat_results_transcripts_for_fc$fc, log2fc=log2(stat_results_transcripts_for_fc$fc), ensembl_transcript_id=sub("[.].*","",stat_results_transcripts$id), bg_samples_filt)

# stat_results_transcripts = stattest(gowntable=bg_samples_filt, feature=feature, pData=pheno_data_comparison, covariate='group', getFC=TRUE, adjustvars=c('pair')) 
# stat_results_transcripts = data.frame(stat_results_transcripts, log2fc=log2(stat_results_transcripts$fc), bg_samples_filt)
# trans_info = read.csv('/data/NHLBI_BCB/Sack_Kim-Han_RNA_seq/07-ballgown/transcript_annotations_genecodev25.csv')
# stat_results_transcripts = merge(trans_info, stat_results_transcripts, by.x="t_name", by.y="id", all.x=FALSE, all.y=TRUE)
# lncrna_trans_info = read.table('/data/NHLBI_BCB/Fayaz/02-RNA-seq-transcript-level/data_source/GENCODE_lncRNA_Mouse_Human/gencode.v25.long_noncoding_RNAs.transcript_ids.txt', header=T)
# lncrna_trans_info = data.frame(lncrna_transid = lncrna_trans_info$lncrna_transid, lncrna_transid = lncrna_trans_info$lncrna_transid)
# stat_results_transcripts = merge(x=stat_results_transcripts, y=lncrna_trans_info, by.x="t_name", by.y="lncrna_transid", all.x=TRUE)

##custom model
# mod = model.matrix(~ pheno_data_comparison$pair + pheno_data_comparison$group) #run for genes
# mod0 = model.matrix(~ pheno_data_comparison$group) #run for genes
# stat_results_transcripts_for_fc = stattest(gowntable=bg_samples_filt, feature=feature, pData=pheno_data_comparison, covariate='group', getFC=TRUE, libadjust = FALSE) #run for genes
# stat_results_transcripts = stattest(gowntable=bg_samples_filt, feature=feature, pData=pheno_data_comparison, mod0=mod0, mod=mod) #run for genes
# stat_results_transcripts = data.frame(stat_results_transcripts, fc=stat_results_transcripts_for_fc$fc, log2fc=log2(stat_results_transcripts_for_fc$fc), ensembl_gene_id=sub("[.].*","",stat_results_transcripts$id), bg_samples_filt) #run for genes

# stat_results_transcripts = stattest(gowntable=bg_samples_filt, feature=feature, pData=pheno_data_comparison, covariate='group', getFC=TRUE, adjustvars=c('pair')) 
# stat_results_transcripts = data.frame(stat_results_transcripts, log2fc=log2(stat_results_transcripts$fc), bg_samples_filt)
# gene_info = read.csv('/data/NHLBI_BCB/Sack_Kim-Han_RNA_seq/07-ballgown/gene_annotations_genecodev25.csv')
# stat_results_transcripts = merge(gene_info, stat_results_transcripts, by.x="gene_id", by.y="id", all.x=FALSE, all.y = TRUE)
# lncrna_gene_info = read.table('/data/NHLBI_BCB/Fayaz/02-RNA-seq-transcript-level/data_source/GENCODE_lncRNA_Mouse_Human/gencode.v25.long_noncoding_RNAs.gene_ids.txt', header=T)
# lncrna_gene_info = data.frame(lncrna_geneid = lncrna_gene_info$lncrna_geneid, lncrna_geneid = lncrna_gene_info$lncrna_geneid)
# stat_results_transcripts = merge(x=stat_results_transcripts, y=lncrna_gene_info, by.x="gene_id", by.y="lncrna_geneid", all.x=TRUE)

## RUN THIS WITHOUT COVARIATES
#stat_results_transcripts = stattest(gowntable=bg_samples_filt, pData=pheno_data_comparison, covariate='group', feature=feature, getFC=TRUE)
## Add gene name, gene ids, transcript names, count values for comparison groups (sorted by group)
#stat_results_transcripts = data.frame(stat_results_transcripts, log2fc=log2(stat_results_transcripts$fc), bg_samples_filt)
#trans_info = read.csv('/data/NHLBI_BCB/Sack_Kim-Han_RNA_seq/07-ballgown/transcript_annotations_genecodev25.csv')
#stat_results_transcripts = merge(trans_info, stat_results_transcripts, by.x="t_name", by.y="id", all.x=FALSE, all.y=TRUE)
#lncrna_trans_info = read.table('/data/NHLBI_BCB/Fayaz/02-RNA-seq-transcript-level/data_source/GENCODE_lncRNA_Mouse_Human/gencode.v25.long_noncoding_RNAs.transcript_ids.txt', header=T)
#lncrna_trans_info = data.frame(lncrna_transid = lncrna_trans_info$lncrna_transid, lncrna_transid = lncrna_trans_info$lncrna_transid)
#stat_results_transcripts = merge(x=stat_results_transcripts, y=lncrna_trans_info, by.x="t_name", by.y="lncrna_transid", all.x=TRUE)

#stat_results_transcripts = data.frame(stat_results_transcripts, log2fc=log2(stat_results_transcripts$fc), bg_samples_filt) #run for genes
#gene_info = read.csv('/data/NHLBI_BCB/Sack_Kim-Han_RNA_seq/07-ballgown/gene_annotations_genecodev25.csv')
#stat_results_transcripts = merge(gene_info, stat_results_transcripts, by.x="gene_id", by.y="id", all.x=FALSE, all.y = TRUE)
#lncrna_gene_info = read.table('/data/NHLBI_BCB/Fayaz/02-RNA-seq-transcript-level/data_source/GENCODE_lncRNA_Mouse_Human/gencode.v25.long_noncoding_RNAs.gene_ids.txt', header=T)
#lncrna_gene_info = data.frame(lncrna_geneid = lncrna_gene_info$lncrna_geneid, lncrna_geneid = lncrna_gene_info$lncrna_geneid)
#stat_results_transcripts = merge(x=stat_results_transcripts, y=lncrna_gene_info, by.x="gene_id", by.y="lncrna_geneid", all.x=TRUE)

# listMarts()
ensembl = useMart("ensembl")
# listDatasets(ensembl)
ensembl = useMart("ensembl", dataset="hsapiens_gene_ensembl")
# listAttributes(ensembl)
stat_results_transcripts_annotated = getBM(attributes=c("ensembl_transcript_id","ensembl_gene_id","external_gene_name"), filters="ensembl_transcript_id", values=sub("[.].*","",stat_results_transcripts$id), mart=ensembl)
stat_results_transcripts_annotated = getBM(attributes=c("ensembl_gene_id","external_gene_name"), filters="ensembl_gene_id", values=sub("[.].*","",stat_results_transcripts$id), mart=ensembl) #run for genes
stat_results_transcripts = merge(stat_results_transcripts, stat_results_transcripts_annotated, by="ensembl_transcript_id", all=TRUE)
stat_results_transcripts = merge(stat_results_transcripts, stat_results_transcripts_annotated, by="ensembl_gene_id", all=TRUE) #run for genes

## Sort results from smallest p-value
stat_results_transcripts = arrange(stat_results_transcripts, pval)
## Write results to CSV
write.csv(stat_results_transcripts, paste("stat_results_",feature,"_count_",group1,"_",group2,".csv", sep=""), row.names=FALSE)

# pca plot - using results with p-value < 0.05 or any other subset based on significant DE for e.g. fc, log2fc, transcripts, genes, counts

# stat_results_transcripts_subset = subset(stat_results_transcripts, subset=stat_results_transcripts$log2fc < -1 | stat_results_transcripts$log2fc > 1)
stat_results_transcripts_subset = subset(stat_results_transcripts, subset=stat_results_transcripts$pval < 0.05)
stat_results_transcripts_subset = subset(stat_results_transcripts_subset, select = grep("S", names(stat_results_transcripts_subset)))

# stat_results_transcripts_subset = subset(stat_results_transcripts, select = grep("S", names(stat_results_transcripts))) #all transcripts or all genes (after applying initial filter of expression)

pca = prcomp(t(log2(stat_results_transcripts_subset+1)), scale = F)
# pca = prcomp(t(stat_results_transcripts_subset), scale = F)
pdf(paste("pca_plot_p_less_than_0.05_pairs_plot_first_five_PCs_",feature,"_",group1,"_",group2,".pdf", sep=""))
pairs(pca$x[,1:5], main=paste("scatter plot matrix of Principal Components_first_five_PCs_",feature,"_",group1,"_",group2, sep=""), pch=21, bg=c('blue','orange')[as.vector(unclass(as.factor(as.vector(pheno_data_comparison$group))))], cex.main=0.75)
dev.off()
pdf(paste("pca_plot_p_less_than_0.05_pc1_vs_pc2_",feature,"_",group1,"_",group2,".pdf", sep=""))
# Add extra space to right of plot area; change clipping to figure
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
plot(pca$x[,1], pca$x[,2], main=paste("scatter plot matrix of Principal Components_PC1.vs.PC2_",feature,"_",group1,"_",group2, sep=""), pch=16, col=c('blue','orange')[as.vector(unclass(as.factor(as.vector(pheno_data_comparison$group))))], cex=1.20, cex.main=0.90, xlab="PC1", ylab="PC2")
text(pca$x[,1], pca$x[,2], labels=pheno_data_comparison$SampleName, pos=3, cex=0.70)
legend("topright", inset=c(-0.2,0), legend=c(group1,group2), pch=c(16,16), col=levels(as.factor(c('blue','orange')[as.vector(unclass(as.factor(as.vector(pheno_data_comparison$group))))])), title="Group", cex=0.70)
dev.off()

pca = prcomp(t(stat_results_transcripts_subset), scale = F)
pdf(paste("pca_plot_pairs_plot_first_five_PCs_",feature,"_",group1,"_",group2,".pdf", sep=""))
pairs(pca$x[,1:5], main=paste("scatter plot matrix of Principal Components_first_five_PCs_",feature,"_",group1,"_",group2, sep=""), pch=21, bg=c('green','red')[as.vector(unclass(as.factor(as.vector(pheno_data_comparison$group))))], cex.main=0.75)
dev.off()
pdf(paste("pca_plot_pc1_vs_pc2_",feature,"_",group1,"_",group2,".pdf", sep=""))
# Add extra space to right of plot area; change clipping to figure
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
plot(pca$x[,1], pca$x[,2], main=paste("scatter plot matrix of Principal Components_PC1.vs.PC2_",feature,"_",group1,"_",group2, sep=""), pch=16, col=c('green','red')[as.vector(unclass(as.factor(as.vector(pheno_data_comparison$group))))], cex=1.20, cex.main=0.90, xlab="PC1", ylab="PC2")
text(pca$x[,1], pca$x[,2], labels=pheno_data_comparison$SampleName, pos=3, cex=0.70)
legend("topright", inset=c(-0.2,0), legend=c(group1,group2), pch=c(16,16), col=levels(as.factor(c('green','red')[as.vector(unclass(as.factor(as.vector(pheno_data_comparison$group))))])), title="Group", cex=0.70)
dev.off()

pdf(paste("pca_plot_transcripts_pc2_vs_pc5_",feature,"_",group1,"_",group2,".pdf", sep=""))
# Add extra space to right of plot area; change clipping to figure
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
plot(pca$x[,2], pca$x[,5], main=paste("scatter plot matrix of Principal Components_PC2.vs.PC5_",feature,"_",group1,"_",group2, sep=""), pch=16, col=c('green','red')[as.vector(unclass(as.factor(as.vector(pheno_data_comparison$group))))], cex=1.20, cex.main=0.90, xlab="PC2", ylab="PC5")
text(pca$x[,2], pca$x[,5], labels=pheno_data_comparison$SampleName, pos=3, cex=0.70)
legend("topright", inset=c(-0.2,0), legend=c(group1,group2), pch=c(16,16), col=levels(as.factor(c('green','red')[as.vector(unclass(as.factor(as.vector(pheno_data_comparison$group))))])), title="Group", cex=0.70)
dev.off()

## Plotting abundance distribution, before and after applying filter, transcripts, genes, fpkm
counts_comparison = subset(stat_results_transcripts, select = grep("S", names(stat_results_transcripts)))
counts = counts_comparison
colnames(counts) = pheno_data_comparison$SampleName_group
# counts = log2(cpm(bg_samples_filt) +1)
# counts = log2(cpm(counts) +1)
counts = log2(counts+1)
pdf(file=paste("counts_distribution_",feature,"_",group1,"_",group2,".pdf", sep=""))
boxplot(counts, col=c('blue','orange')[as.vector(unclass(as.factor(as.vector(pheno_data_comparison$group))))], las=2, ylab='log2(counts+1)', main=paste("counts distribution:",feature,"_",group1,"_",group2), xlab="sample", cex.main=0.80, cex.axis=0.50)
dev.off()

## Plotting volcano plot, Differential Expression Analysis-qval, transcripts, genes, fpkm
data = read.csv(file=paste("stat_results_",feature,"_count_",group1,"_",group2,".csv", sep=""), header=T, row.names=1)
q.trans = -1*log(data$qval, base = 10)
logfold = log(data$fc, base = 2)
pdf(file=paste("volcano_qval_",feature,"_",group1,"_",group2,".pdf", sep=""))
plot(logfold,q.trans,type="n",ylab="-1*log10(q-value)",xlab="log2(fold change)",main=paste("COUNT:",feature,"_",group1,"_",group2), cex.main=0.80, xlim=range(logfold), ylim=range(q.trans))
points(logfold,q.trans,col="black",cex=0.65)
points(logfold[(q.trans>1.3&logfold>1.0)],q.trans[(q.trans>1.3&logfold>1.0)],col="red",pch=16,cex=0.65)
points(logfold[(q.trans>1.3&logfold<(-1.0))],q.trans[(q.trans>1.3&logfold<(-1.0))],col="green",pch=16,cex=0.65)
# text(logfold[(q.trans>1.3&logfold>1.0)],q.trans[(q.trans>1.3&logfold>1.0)],labels=as.character(data$geneName[(q.trans>1.3&logfold>1.0)]), cex=0.65)
# text(logfold[(q.trans>1.3&logfold<(-1.0))],q.trans[(q.trans>1.3&logfold<(-1.0))],labels=as.character(data$geneName[(q.trans>1.3&logfold<(-1.0))]), cex=0.65)
abline(h=1.3)
abline(v=-1.0)
abline(v=1.0)
dev.off()

## Plotting volcano plot, Differential Expression Analysis-pval, transcripts, genes, fpkm
p.trans = -1*log(data$pval, base = 10)
logfold = log(data$fc, base = 2)
pdf(file=paste("volcano_pval_",feature,"_",group1,"_",group2,".pdf", sep=""))
plot(logfold,p.trans,type="n",ylab="-1*log10(p-value)",xlab="log2(fold change)",main=paste("COUNT:",feature,"_",group1,"_",group2), cex.main=0.80, xlim=range(logfold), ylim=range(p.trans))
points(logfold,p.trans,col="black",cex=0.65)
points(logfold[(p.trans>1.3&logfold>1.0)],p.trans[(p.trans>1.3&logfold>1.0)],col="red",pch=16,cex=0.65)
points(logfold[(p.trans>1.3&logfold<(-1.0))],p.trans[(p.trans>1.3&logfold<(-1.0))],col="green",pch=16,cex=0.65)
# text(logfold[(p.trans>1.3&logfold>1.0)],p.trans[(p.trans>1.3&logfold>1.0)],labels=as.character(data$geneName[(p.trans>1.3&logfold>1.0)]), cex=0.65)
# text(logfold[(p.trans>1.3&logfold<(-1.0))],p.trans[(p.trans>1.3&logfold<(-1.0))],labels=as.character(data$geneName[(p.trans>1.3&logfold<(-1.0))]), cex=0.65)
abline(h=1.3)
abline(v=-1.0)
abline(v=1.0)
dev.off()

pdf(file=paste("qval_distribution_",feature,"_",group1,"_",group2,".pdf", sep=""))
hist(data$qval, col="grey", border="white", xlab="qval", ylab="", main=paste("COUNT:",feature,"_",group1,"_",group2))
dev.off()
pdf(file=paste("pval_distribution_",feature,"_",group1,"_",group2,".pdf", sep=""))
hist(data$pval, col="grey", border="white", xlab="pval", ylab="", main=paste("COUNT:",feature,"_",group1,"_",group2))
dev.off()

# transcript_fpkm = texpr(bg_samples_filt, 'FPKM')
# transcript_cov = texpr(bg_samples, 'cov')
# whole_tx_table = texpr(bg_samples, 'all')
# exon_mcov = eexpr(bg_samples, 'mcov')
# junction_rcount = iexpr(bg_samples)
# whole_intron_table = iexpr(bg_samples, 'all')
# gene_expression = gexpr(bg_samples)

# transcript_gene_table = indexes(bg_samples_filt)$t2g
# head(transcript_gene_table)

## Plotting setup
#tropical = c('darkorange', 'dodgerblue', 'hotpink', 'limegreen', 'yellow')
#palette(tropical)

## Plotting gene abundance distribution
# fpkm = texpr(bg_samples, meas='FPKM')
# fpkm = log2(fpkm +1)
# boxplot(fpkm, col=as.numeric(pheno_data$Gender), las=2,ylab='log2(FPKM+1)')

## Plot individual transcripts
#ballgown::transcriptNames(bg_chrX)[12]
#plot(fpkm[12,] ~ pheno_data$sex, border=c(1,2),
#     main=paste(ballgown::geneNames(bg_chrX)[12], ' : ',ballgown::transcriptNames(bg_chrX)[12]),
#     pch=19, xlab="Sex", ylab='log2(FPKM+1)')
#points(fpkm[12,] ~ jitter(as.numeric(pheno_data$sex)), col=as.numeric(pheno_data$sex))

## Plot gene of transcript 1729
#plotTranscripts(ballgown::geneIDs(bg_chrX)[1729], bg_chrX,
#                main=c('Gene XIST in sample ERR188234'), sample=c('ERR188234'))

## Plot average expression
# plotMeans(ballgown::geneIDs(bg_samples)[203], bg_samples, groupvar="sex", legend=FALSE)

transcripts_fpkm_bg = read.csv(file='/data/NHLBI_BCB/Fayaz/02-RNA-seq-transcript-level/02-HISAT2-StringTie-Ballgown/04-ballgown/human_only/stat_results_transcripts_untreated_vs_treated.csv', header= T, row.names=3)
transcripts_count_bg = read.csv(file='/data/NHLBI_BCB/Fayaz/02-RNA-seq-transcript-level/02-HISAT2-StringTie-Ballgown/04-ballgown/human_only_count/stat_results_transcripts_1_2.csv', header=T, row.names=3)
transcripts_fpkm_counts_bg = merge(transcripts_fpkm_bg, transcripts_count_bg, by="row.names")
# plot(x=-log(transcripts_fpkm_counts_bg$pval.x,10), y=-log(transcripts_fpkm_counts_bg$pval.y,10))
# dev.off()

library(ffpe)
CATplot(vec1=cbind(rownames(tfpkm_res), abs(tfpkm_res$log2fc)), vec2=cbind(rownames(tcount_res), abs(tcount_res$log2fc)), main="transcript-fpkm.vs.count-concordance.by.log2fc",maxrank=1000,xlab="Size of top-ranked gene lists",ylab="Concordance")
legend("topright",lty=1:2,legend=c("Actual concordance","Concordance expected by chance"), bty="n")
dev.off()

CATplot(vec1=cbind(rownames(gfpkm_res), abs(gfpkm_res$log2fc)), vec2=cbind(rownames(gcount_res), abs(gcount_res$log2fc)), main="gene-fpkm.vs.count-concordance.by.log2fc",maxrank=1000,xlab="Size of top-ranked gene lists",ylab="Concordance")
legend("topleft",lty=1:2,legend=c("Actual concordance","Concordance expected by chance"), bty="n")
dev.off()




