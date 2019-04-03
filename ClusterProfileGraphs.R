rm(list=ls())
library("clusterProfiler")
library(ReactomePA)
library(igraph)
library(GO.db)
library(WGCNA)
#install.packages("org.Hs.eg.db", repos="http://bioconductor.org/packages/3.1/bioc")
library(org.Hs.eg.db)
library(DOSE)
#install.packages("data.table")
library(data.table)
library(IRanges)
library(S4Vectors)

aa<-read.delim("CytoscapeInput-nodes-OneCondition0brown TOP 200.9percentFED2ABStop.txt", 
               header = T)
annotid<-aa$nodeName

install.packages("devtools")
library(devtools)
BiocInstaller::biocLite('grimbough/biomaRt')
library("biomaRt")
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))

FullGene_list2 <- getBM(filters= "external_gene_name", 
                        attributes= c("ensembl_gene_id",
                                      "entrezgene"),values=annotid,mart= mart)
geneID<-FullGene_list2[!is.na(FullGene_list2$entrezgene),]


#eNRICHEMENT WITH NOMINAL PVALUE 0.05 AND QVALUE 0.05
ego2 <- enrichGO(gene         = FullGene_list2$entrezgene,
                 OrgDb         = org.Hs.eg.db,
                 ont           = c("BP", "CC", "MF"),
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.05,
                 qvalueCutoff  = 0.05,
                 readable      = TRUE)
head(ego2)

write.csv(ego2,"clusterProfiler_module.csv")

