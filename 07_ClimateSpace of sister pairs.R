###################################################
### Clade niche
###################################################
setwd(".//Acaena niche evolution")
source("06_Clade pairing.R")

library(dplyr)
library(grid)
library(gridExtra)
library(ggplot2)

###################################################
###  Climate data preparation
###################################################

# da1 <- read.csv("Y:\\acaena_bioclim_landcover_history_inclNAonland.csv")
# d <- da1[is.na(da1$landCoverChange) == F, ]
# 
# # sp names
# sname <- colnames(d)[grepl("^Acaena", colnames(d))]
# # Acaena emittens and A. minor have no occurrence in primary and/or secondary habitats.
# s <- sname[ - c(6, 13)]
# 
# for(i in sname){
#   d[is.na(d[,i]),i] <- 0
# }
# 
# # get env. corrdinates (PCA axes)
# pca <- prcomp(d[, paste("bioclim", c(1, 6, 12, 15), sep = "")],
#               center = TRUE,
#               scale. = TRUE)
# scores <- data.frame(d[, c(colnames(d)[grep("^bioclim", colnames(d))], sname,
#                            "x", "y" #, "preLandcover", "currentLandcover", "landCoverChange"
#                            )], pca$x[, 1:2])
# save(scores, file = ".//Scores.data")

load("..//Scores.data")
extent_x = c(min(scores$PC1), max(scores$PC1))
extent_y = c(min(scores$PC2), max(scores$PC2))


########################################################################################
### Node niche
########################################################################################

# Modify species names in phylogentic distance file
a <- unlist(strsplit(rownames(nodes), "_Ac"))
a2 <- gsub("_EU352216", "", a) %>% gsub("_AY634821", "", .) %>% gsub("novae-", "novae.", .)
rownames(nodes) <- grepl("_", a2) %>% a2[.]

### Make species name codes
codes2 <- makeTag(colnames(scores), "Acaena")


########################################################################################
### Plot clade niche of internal nodes
########################################################################################

source(".//generateClimateDataOfClades.R")
source(".//plotClimateSpaceWithSpNameList.R")

### Run
# The following nodes must be eliminated. See "res" object.
# Check which species are not shared between phylogenetic tree and occurrence record data
rownames(nodes)[!(rownames(nodes) %in% colnames(scores))]
colnames(scores)[!(colnames(scores) %in% rownames(nodes))]

# No occurrence data for the following nodes
noOccSpNo <- which(!(rownames(nodes) %in% colnames(scores)))
# No occurrence data for their sister clades
noOccSisSpNo <- which(allnodesister %in% noOccSpNo)

## The following nodes must be eliminated. Because no occurrence data for them or their sister clades.
number <- (1:max(acaena$edge))[-c(noOccSpNo,
                                 noOccSisSpNo, 17, 45,
                                 24, # Oldest ancestor node
                                 25, 26, # Outgroup nodes
                                 27 # Sister nodes of Outgroup nodes
)]

cladedata <- lapply(number, function(i){
  
  clades <- generateClimateDataOfClades(i, acaena, allnodesister, scores,
                                        nodes = nodes, tips = tips, spnameCodes = codes2)
  
  return(clades)
  }
  
  )

save(cladedata, file = ".//cladePairData.data")

i=1
clades=cladedata[[1]]
ploTwoGroupWithSpNames(background = scores,
                       axis1 = "PC1", axis2 = "PC2", # Names of coordinates
                       data1 = clades[[1]], data2 = clades[[3]], # Dataframes of two groups of points
                       col1 = "green", col2 = "purple",
                       nodeName = clades[[2]],
                       sisnodeName = clades[[4]],
                       nodeNumber = i,
                       extent_x, extent_y,
                       save = TRUE
)

sisplots <- ploTwoGroupWithSpNames(background = scores,
                       axis1 = "PC1", axis2 = "PC2", # Names of coordinates
                       data1 = clades[[1]], data2 = clades[[3]], # Dataframes of two groups of points
                       col1 = "green", col2 = "purple",
                       nodeName = clades[[2]],
                       sisnodeName = clades[[4]],
                       nodeNumber = i,
                       extent_x, extent_y,
                       save = FALSE
)