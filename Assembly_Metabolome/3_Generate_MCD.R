# Code adapted from Danczak et al. 2020, https://github.com/danczakre/Meta-Metabolome_Ecology


### Generating the molecular characteristics dendrogram
# RED 2020; robert.danczak@pnnl.gov; danczak.6@gmail.com

##pw##

options(digits = 10)

require(phangorn) # For tree based functions
require(ggtree) # For tree visualization
require(vegan) # For vegdist
library(tidyverse)

# ################## #
#### Load in data ####
# ################## #

# Set sample name
Sample_Name = "Abisko"

# Load in data

mol = read.csv("output/Processed_Abisko_Mol.csv", row.names = 1)

### Ensuring that isotopic peaks are removed
if(length(which(colnames(data) %in% "C13")) > 0){
  w = row.names(data)[which(data$C13 > 0)]
  
  if(length(w) > 0){
    data = data[-which(row.names(data) %in% w),]
  }
  
  rm("w")
}

# Removing peaks that have no formula assignments

mol = mol[-which(mol$MolForm %in% NA),]

# Setting objects for useful parameters
Mol.Info = mol[,c("C", "H", "O", "N", "S", "P", "DBE", "AI_Mod", "kdefect.CH2"), drop = F]
Mol.Ratio = mol[,c("OtoC_ratio", "HtoC_ratio", "NtoC_ratio", "PtoC_ratio", "NtoP_ratio")]


# ##################### #
#### Generating tree ####
# ##################### #

# Pairwise distance between peaks
Mol.Info = as.data.frame(apply(Mol.Info, 2, scale), row.names = row.names(Mol.Info)) # Generating a distance matrix based upon the provided parameters

# Create tree
##Adding na.rm = TRUE to avoid error (only for active analysis)##
tree = as.phylo(hclust(vegdist(Mol.Info, "euclidean", na.rm = TRUE), "average")) # Converting the distance matrix into a tree

# Writing tre
write.tree(tree, paste(Sample_Name, "_MCD_UPGMA.tre", sep = ""))




