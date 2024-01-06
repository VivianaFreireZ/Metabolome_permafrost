# Code adapted from Danczak et al. 2020, https://github.com/danczakre/Meta-Metabolome_Ecology

### Calculating null bMNTD variants
# RED 2020; robert.danczak@pnnl.gov; danczak.6@gmail.com

# Modified from https://github.com/stegen/Stegen_etal_ISME_2013/blob/master/bNTI_Local_Machine.r; Stegen et al, 2013
# Modified from https://github.com/stegen/Stegen_etal_ISME_2013/blob/master/Raup_Crick_Abundance.r; Stegen et al, 2013


range = 1:999
Sample_Name = "Abisko"
tree.type = "TWCD"

#-----------------#

library(vegan)
library(picante)

###################################
#### Data Loading and cleaning ####
###################################

data = read.csv("output/Processed_Abisko_Data.csv", row.names = 1) # Importing the organismal data  
tree = read.tree("output/Abisko_Weighted_All-Trans_UPGMA.tre") # Importing the tree

# Creating necessary directories

if(!dir.exists(paste0(tree.type, "_Null_Results"))){
  dir.create(paste0(tree.type, "_Null_Results"))
}


####################################
#### Beginning the bNTI Process ####
####################################

# Converting to presence/absence
data[data>1] = 1

# Matching the tree to the newly rarefied OTU dataset
phylo = match.phylo.data(tree, data)

# Running cophenetic outside of the for-loop
coph = cophenetic(phylo$phy)

# Calculating the bMNTD for 999 random distributions
print(paste(date(), " - Start for loop"))

for(i in range){
  bMNTD.rand = as.matrix(comdistnt(t(phylo$data), taxaShuffle(coph), abundance.weighted = F, exclude.conspecifics = F))
  write.csv(bMNTD.rand, paste(tree.type, "_Null_Results_", Sample_Name, "_", tree.type, "_bMNTD_rep", i, ".csv", sep = ""), quote = F)
  rm("bMNTD.rand")
  
  print(c(date(),i))
} # Performing the calculations on using the OTU table but with randomized taxonomic affiliations

print(paste(date(), " - End for loop"))



