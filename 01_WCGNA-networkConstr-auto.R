#=====================================================================================
#
#  Code chunk 1
#
#=====================================================================================
rm(list=ls())

# Display the current working directory
getwd();
# If necessary, change the path below to the directory where the data files are stored. 
# "." means current directory. On Windows use a forward slash / instead of the usual \.
workingDir = ".";
setwd(workingDir); 
# Load the WGCNA package
library(WGCNA)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
# Allow multi-threading within WGCNA. This helps speed up certain calculations.
# At present this call is necessary for the code to work.
# Any error here may be ignored but you may want to update WGCNA if you see one.
# Caution: skip this line if you run RStudio or other third-party R environments. 
# See note above.
allowWGCNAThreads() 
#enableWGCNAThreads()
# Load the data saved in the first part
lnames = load(file = "stat_fast_90percentobo2_log2_dataInput.RData");
#The variable lnames contains the names of loaded variables.
lnames


#=====================================================================================
#
#  Code chunk 2
#
#=====================================================================================


# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=30, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(ExprDataFinal, powerVector = powers, verbose = 5)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndies[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold ",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale"));
text(sft$fitIndices[,1], -sign(sft$fitIndies[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.8,col="red")

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold ",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))

text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")


#=====================================================================================
#
#  Code chunk 3
#
#=====================================================================================

##ks edits added "maxBlockSize = 16305," to make a single block (for better clustering)/. If you 
##have 8GB, you can try maxBlockSize about 12000.  You will get fewer 
##blocks and most likely a better clustering/modules. 
#while using avgFPKMobo1_ fed data power=6/ fast data power=7
#while using .9percentFPKMobo1_ fast data power=10 fed data power=9(minModuleSize=100(default is 30))
#while using .9percentFPKMobo1ABS_ fast data power=9 fed data power=8(minModuleSize=100(default is 30))
##ks edits ends==============================
net = blockwiseModules(ExprDataFinal, power = 10,maxBlockSize = 23000,
                       TOMType = "signed", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = FALSE,
                       saveTOMFileBase = "FED_2_DEseq2_VarianceStabilizingTransformation", 
                       verbose = 3)

#always chech if mergeCutHeight = 0.25, sometimes I change it to 0.1!!!!!!!!!!!!!
#=====================================================================================
#
#  Code chunk 4
#
#=====================================================================================


# open a graphics window

sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)


#=====================================================================================
#
#  Code chunk 5
#
#=====================================================================================


moduleLabels = net$colors
#moduleColors = labels2colors(net$colors)
moduleColors = labels2colors(moduleLabels)
##ksedits--------------------------------------------
#good to have the module numbers and the number of genes present in each module information store.
#note: columns number 0 refers to the genes that were unassigned to any module
moduleList=as.data.frame(table(moduleLabels))
colorList=as.data.frame(table(moduleColors))
ModuleListColor=cbind(moduleList,colorList)
#pdf(file="ModuleNumbers_moduleSize_moduleColors_FAST_1_r0.8_fpkmobo1.pdf", wi=8,he=12)
#library(gridExtra)
#grid.table(moduleList)
#dev.off()
moduleColorList=as.data.frame(table(moduleColors))
write.csv(ModuleListColor, file="ModuleNumbers_moduleSize_pc_nc_cutoffs_power10.csv")
#write.csv(moduleColorList, file="ModuleColors_moduleSize_FAST_r0.8_fpkmobo1ABS.csv")
#ks edits ends---------------------------------------------
#ks edit========
#moduleColors = mergedColors[net$blockGenes[[1]]]
#ks edit ends=============
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
save(MEs, moduleLabels, moduleColors, geneTree, 
     file = "Network_fast_90percentobo2_log2_dataInput.RData.RData")


