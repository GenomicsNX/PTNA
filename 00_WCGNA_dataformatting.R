rm(list=ls())
library(WGCNA)
library(dplyr)
library(dynamicTreeCut)
library(fastcluster)
library(stringi)
options(stringsAsFactors = FALSE)
allowWGCNAThreads()

###Data input formatting
#Base.data=read.csv("stat_results_gene_fpkm_FED_2_average_FPKM_obo2_log2.csv")
Base.data=read.csv("stat_Fast.csv")
names(Base.data)
#first part of data
Base.data=Base.data1[,c(1:14,22)]
#second part of data
Base.data=Base.data1[,c(1:7,15:21)]
Base.data.pre=Base.data[,c(8:14)]+1
Base.data.log=log(Base.data.pre,2)
Final.log.data=cbind(Base.data[,c(1:7)], Base.data.log)
write.csv(Final.log.data, file = "stat_Fast_90percentobo1_log2.csv")

##filtering low expressing data points

gsg=goodSamplesGenes(ExprData1, verbose = 3)
gsg$allOK

if (!gsg$allOK)
{
  if(sum(gsg$goodGenes)>0)
    printFlush(paste("Removing genes:",paste(names(ExprData1)[!gsg$goodGenes], collapse=",")))
  if(sum(gsg$goodSamples)>0)
    printFlush(paste("Removing samples:",paste(names(ExprData1)[!gsg$goodSamples], collapse=",")))
  
  ExprData1=ExprData1[gsg$goodSamples,gsg$goodGenes]
}


##ks edits...just do this step.
ExprDataFinal=ExprData1
nGenes=ncol(ExprDataFinal)
nSamples=nrow(ExprDataFinal)

save(ExprDataFinal,nGenes,nSamples, file="stat_fast_90percentobo2_log2_dataInput.RData")


##there is a need to make a good gene annotation file from here:
#after applying the gsg$goodGene filter, save the good genes order along with the gene_id and gene_name
#in a separate csv file as follows:
BaseDataWithGoodGenes1<-as.data.frame(Base.data[gsg$goodGenes,])
FAstData_goodgenes<-BaseDataWithGoodGenes1
FEdData_goodgenes<-BaseDataWithGoodGenes
FedGenes<-FEdData_goodgenes[,c(2,3)]

FastGSG<-left_join(FedGenes,FAstData_goodgenes)
FastGSG1<-FastGSG[!is.na(FastGSG$ensembl_gene_id),]
FastGSG2<-FastGSG1[!is.na(FastGSG1$hgnc_symbol),]
write.csv(FastGSG2, file="DEseq2_VarianceStabilizingTransformation_FAST1_geneNameID_34374.csv")
FastGenes<-FastGSG1[,c(1,2)]
FedGSG<-left_join(FastGenes,FEdData_goodgenes)
FedGSG1<-FedGSG[!is.na(FedGSG$ensembl_gene_id),]
FedGSG2<-FedGSG1[!is.na(FedGSG1$hgnc_symbol),]
write.csv(FedGSG2, file="DEseq2_VarianceStabilizingTransformation_FED2_geneNameID_34374.csv")
GeneIDSymbol=as.data.frame(BaseDataWithGoodGenes[gsg$goodGenes,][-c(3:32)])
GeneIDSymbol=FedGSG2[,c(1,2)]
GeneIDSymbol=Base.data[,c(1,2)]
write.csv(GeneIDSymbol, file="GeneAnnotKS_FAST_FED_VarianceStabilizingTransformation19168_oneCondition.csv")
#the third column of the csv file is gene_name


