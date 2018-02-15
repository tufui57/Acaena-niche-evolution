###################################################
### Clade niche overlap/volume
###################################################

load(".//Scores.data")
load(".//cladePairData.data")

setwd(".//Acaena niche evolution")

library(ecospat)
# tidyverse loads too much DLLs. Load tidyr instead.
library(tidyr)
library(nichePlot)

###################################################
# Calculate niche overlap between sister clades
###################################################

scho <- list()

for(i in 1:length(cladedata)){

  scho[[i]] <- SchoenerD_ecospat(scores, "PC1", "PC2",
                                 cladedata[[i]][[1]],  cladedata[[i]][[3]] # Data of target clade pair
                                 )

}


scholist <- lapply(1:length(scho), function(i){
  # Convert list to dataframe
  x <- data.frame(scho[[i]])

  # Add colnames
  colnames(x) <- c("ecospat.corrected", "ecospat.uncorrected")

  # Get Shoenner's D. I don't use I.
  return(x["D",])

  }
)

# Add column in schonner's D dataframe showing clade numbers
scholist <- do.call(rbind, scholist)
nodeNo <- sapply(cladedata, "[[", 5) %>% strsplit(., "\ ")

scholist$node1 <- as.numeric(do.call(rbind, lapply(nodeNo, "[[",1)))
scholist$node2 <- as.numeric(do.call(rbind, lapply(nodeNo, "[[",2)))

write.csv(scholist, ".//clade_schoennerD_acaena.csv")

###################################################
### Clade niche volume of clades
###################################################

nichevol <- list()
for(i in 1:length(cladedata)){
  
  nichevol[[i]] <- SchoenerD_ecospat(
    
    # I calculate just one clade of each sister clade pair.
    data1 = cladedata[[i]][[1]], # data of nodes including internal nodes
    background = scores, data2 = scores, "PC1", "PC2"
  )
  
}

source(".//SchonnerDdataframeFormat.r")
nichevoldata <- SchonnerDdataframeFormat(nichevol)


###### Find which nodes I haven't calculated niche volume
# I need another clade in each pair, but I need avoid duplicated one. 

source("06_Clade pairing.R")

# Total number of nodes in tree 
noTreeNode <- max(
  c(acaena$edge[,1], acaena$edge[,2])
                  )
no1 <- (1:noTreeNode)[!(1:noTreeNode %in% nichevoldata$node1)]


# Nodes with no occurrence data, ancestor nodes for all species, nodes with overlapping species
noNotClac <- no1[!(no1 %in% c(10,18,13,14,21:28,42,44,45))]

source("generateClimateDataOfClades.R")

missingcladedata <- list()

for(i in noNotClac){
  cladeScore <- generateClimateDataOfTargetNode(i, acaena,
                                allnodesister, scores, nodes, tips)
  cladedata <- cladeScore[cladeScore$targetClade == 1,]
  
  missingcladedata[[i]] <- cladedata
 
  }

nichevolmissinglist <- list()
for(i in noNotClac){
  
  nichevolmissinglist[[i]] <- SchoenerD_ecospat(background = scores, 
                                       data1 = missingcladedata[[i]], data2 = scores,
                                       "PC1", "PC2"
                                       )
}

nichevoldatamissing <- SchonnerDdataframeFormat(nichevolmissinglist)

# Omit missing data rows
missingvol <- (nichevoldatamissing$ecospat.corrected != 999) %>% nichevoldatamissing[., ]
nichevolume <- rbind(nichevoldata %>% select(-node2), missingvol)

library(dplyr)

write.csv(nichevolume, "Y://clade_nicheVolume_acaena.csv")


