###################################################
### Clade niche
###################################################

setwd(".//Acaena niche evolution")
source("06_Clade pairing.R")

### Unload libraries
loadedNamespaces() %>% unloadNamespace()


library(tidyverse)
library(grid)
library(gridExtra)
library(ggplot2)

###################################################
###  Climate data preparation
###################################################

da1 <- read.csv("Y:\\acaena_bioclim_landcover_history_inclNAonland.csv")
d <- da1[is.na(da1$landCoverChange) == F, ]

# sp names
sname <- colnames(d)[grepl("^Acaena", colnames(d))]
# Acaena emittens and A. minor have no occurrence in primary and/or secondary habitats.
s <- sname[ - c(6, 13)]

for(i in sname){
  d[is.na(d[,i]),i] <- 0
}

# get env. corrdinates (PCA axes)
pca <- prcomp(d[, paste("bioclim", c(1, 6, 12, 15), sep = "")],
              center = TRUE,
              scale. = TRUE)
scores <- data.frame(d[, c(colnames(d)[grep("^bioclim", colnames(d))], sname,
                           "x", "y" #, "preLandcover", "currentLandcover", "landCoverChange"
                           )], pca$x[, 1:2])

extent_x = c(min(scores$PC1), max(scores$PC1))
extent_y = c(min(scores$PC2), max(scores$PC2))



########################################################################################
### Calculate node niche
########################################################################################

# Modify species names in phylogentic distance file
a <- unlist(strsplit(rownames(nodes), "_Ac"))
a2 <- gsub("_EU352216", "", a) %>% gsub("_AY634821", "", .) %>% gsub("novae-", "novae.", .)
rownames(nodes) <- grepl("_", a2) %>% a2[.]

### Import species name codes
codes <- read.csv("Y:\\traits.csv")
codes$X <- gsub("novae_", "novae.",codes$X)
codes2 <- codes[codes$X %in% rownames(nodes), ]


########################################################################################
### Plot clade niche of internal nodes
########################################################################################

source(".//generateClimateDataOfClades.R")
source(".//plotClimateSpaceWithSpNameList.R")

### Run
# The following nodes must be eliminated. See "res" object.

number <- (1:max(acaena$edge))[-c(8,10,13,14,17,18,21:29,34,35,39,40,42,45)]

cladedata <- lapply(number, function(i){
  
  clades <- generateClimateDataOfClades(i, allnodesister, scores)
  
  # ploTwoGroupWithSpNames(background = scores,
  #                        axis1 = "PC1", axis2 = "PC2", # Names of coordinates
  #                        data1 = clades[[1]], data2 = clades[[3]], # Dataframes of two groups of points
  #                        col1 = "green", col2 = "purple",
  #                        nodeName = clades[[5]],
  #                        sisnodeName = clades[[4]],
  #                        nodeNumber = i,
  #                        extent_x, extent_y)
  
  return(clades)
  }
  
  )

save(cladedata, file = "./cladePairData.data")
