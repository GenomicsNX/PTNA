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
# Load the expression and trait data saved in the first part
lnames = load(file = "Input.RData");
#The variable lnames contains the names of loaded variables.
lnames
# Load network data saved in the second part.
lnames = load(file = "Network.RData");
lnames

#===================================================================================================
#
#   use this part if you already have TOM file and would like to get top200 genes in each module
#
#===================================================================================================
#ks edits 1---------------------------------------------------
#generate annotation file
annot = read.csv(file = "stat_Fats_log2.csv");
annot = annot[,c(2:4)]
head(annot)
#modGenes = annot$gene_name[match(modProbes, annot$gene_id)];
#ks edits 2---------------------------------------------
#if you have already saved the TOM file load here:
gnames = load(file = "TOM.RData");
gnames
# Recalculate topological overlap if needed
TOM = TOMsimilarityFromExpr(ExprDataFinal, power = 10);
#save TOM file
save(TOM, file = "TOMmax_Fast.RData")

probes = names(ExprDataFinal)
#modules = "blue";
# Select module probes
intModules = c("turquoise", "blue", "brown","yellow","green","red","pink","black","grey")
intModules = c("turquoise", "blue", "brown","yellow","green","red","black","pink","magenta","purple","greenyellow",
               "lightcyan","tan","salmon","cyan","grey60","lightgreen","lightyellow","darkred","royalblue")
#32 modules
intModules = c("turquoise", "blue", "brown","yellow","green","red","black","pink","magenta","purple","greenyellow",
               "lightcyan","tan","salmon","cyan","grey60","lightgreen","lightyellow","darkred","royalblue",
               "white", "darkorange", "skyblue", "orange", "midnightblue", "steelblue", "saddlebrown", 
               "darkgreen", "darkgrey", "darkturquoise", "paleturquoise","violet","grey")
#35 modules
intModules = c("turquoise", "blue", "brown","yellow","green","red","black","pink","magenta","purple","greenyellow",
               "lightcyan","tan","salmon","cyan","grey60","lightgreen","lightyellow","darkred","royalblue",
               "white", "darkorange", "skyblue", "orange", "midnightblue", "steelblue","saddlebrown", 
               "darkgreen", "darkgrey", "darkturquoise", "paleturquoise","darkolivegreen","darkmagenta","violet")
                
#40 modules
intModules = c("turquoise", "blue", "brown","yellow","green","red","black","pink","magenta","purple","greenyellow",
               "lightcyan","tan","salmon","cyan","grey60","lightgreen","lightyellow","darkred","royalblue",
               "white", "darkorange", "skyblue", "orange", "midnightblue", "steelblue", "saddlebrown", "violet", 
               "darkgreen", "darkgrey", "darkturquoise", "paleturquoise", "darkolivegreen", "darkmagenta", 
               "plum1", "orangered4", "sienna3", "skyblue3", "yellowgreen")
#43 modules
intModules = c("turquoise", "blue", "brown","yellow","green","red","black","pink","magenta","purple","greenyellow",
               "lightcyan","tan","salmon","cyan","grey60","lightgreen","lightyellow","darkred","royalblue",
               "white", "darkorange", "skyblue", "orange", "midnightblue", "steelblue", "saddlebrown", "violet", 
               "darkgreen", "darkgrey", "darkturquoise", "paleturquoise", "darkolivegreen", "darkmagenta", "mediumpurple3",
               "plum1", "orangered4", "sienna3", "skyblue3", "yellowgreen", "lightcyan1", "lightsteelblue1", "plum1")
#59 modules
intModules = c("turquoise", "blue", "brown","yellow","green","red","black","pink","magenta","purple","greenyellow",
               "lightcyan","tan","salmon","cyan","grey60","lightgreen","lightyellow","darkred","royalblue",
               "white", "darkorange", "skyblue", "orange", "midnightblue", "steelblue", "saddlebrown", "violet", 
               "darkgreen", "darkgrey", "darkturquoise", "paleturquoise", "darkolivegreen", "darkmagenta", "mediumpurple3",
               "plum1", "orangered4", "sienna3", "skyblue3", "yellowgreen", "plum1", "plum2", "thistle1", "thistle2", "salmon4",
               "maroon", "lightpink4","lightcyan1","lightsteelblue1", "mediumpurple3", "palevioletred3", "navajowhite2", "bisque4",
               "brown4", "darkslateblue","darkorange2", "floralwhite", "ivory", "honeydew1")
#49 modules
intModules = c("turquoise", "blue", "brown","yellow","green","red","black","pink","magenta","purple","greenyellow",
               "lightcyan","tan","salmon","cyan","grey60","lightgreen","lightyellow","darkred","royalblue",
               "white", "darkorange", "skyblue", "orange", "midnightblue", "steelblue", "saddlebrown", "violet", 
               "darkgreen", "darkgrey", "darkturquoise", "paleturquoise", "darkolivegreen", "darkmagenta", "mediumpurple3",
               "plum1", "orangered4", "sienna3", "skyblue3", "yellowgreen",  
               "lightcyan1","lightsteelblue1", "mediumpurple3", "bisque4",
               "brown4", "darkslateblue","darkorange2", "floralwhite", "ivory",
               "thistle1", "thistle2","palevioletred3","salmon4","plum2")
# (OR) 49+ "thistle1", "thistle2","palevioletred3","salmon4","plum2"
#32 modules
intModules = c("turquoise", "blue", "brown","yellow","green","red","black","pink","magenta","purple","greenyellow",
               "lightcyan","tan","salmon","cyan","grey60","lightgreen","lightyellow","darkred","royalblue",
               "white", "darkorange", "skyblue", "orange", "midnightblue", "steelblue", "saddlebrown","paleturquoise", 
               "darkturquoise", "darkgreen", "darkgrey", "grey")

probes = names(ExprDataFinal)

#==========================export to cytoscape========================

# (OR) # Select modules
modules = "black";
# Select module probes
probes = names(ExprDataFinal)
inModule = is.finite(match(moduleColors, modules));
modProbes = probes[inModule];
#ks edits---------------------------------------------------
annotMod=annot[inModule,]
modGenes=annotMod[,3]
#modGenes = annot$gene_name[match(modProbes, annot$gene_id)];
#ks edits ends---------------------------------------------
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule];
#modTOMtop=modTOM[top, top];
dimnames(modTOM) = list(modProbes, modProbes)
# Export the network into edge and node list files Cytoscape can read
cyt = exportNetworkToCytoscape(modTOM,
  edgeFile = paste("CytoscapeInput-edges-OneCondition0", paste(modules, collapse="-"), ".9percentFED2ABStop.txt", sep=""),
  nodeFile = paste("CytoscapeInput-nodes-OneCondition0", paste(modules, collapse="-"), ".9percentFED2ABStop.txt", sep=""),
  weighted = TRUE,
  threshold = 0,
  nodeNames = modProbes,
  altNodeNames = modGenes,
  nodeAttr = moduleColors[inModule]);


