###################################################
### Clade niche overlap/volume
###################################################

load(".//Scores.data")
load(".//cladePairData.data")

setwd(".//Acaena niche evolution")

library(ecospat)
library(dplyr)
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
### Calculate niche volume of clades
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

source("06_Clade_pairing.R")

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


#################################################################################
### Calculate node ages
#################################################################################

### Distance between all combinations of tips
distances <- dist.nodes(acaena)

# Calculate species ages
ages <- sapply(nodes[[1]], function(i){
  # Phylogenetic distance list has 0 (distance to themselves)
  (distances[,i] > 0) %>% distances[., i] %>% min
}
)

ages <- as_tibble(ages) %>% mutate(., spname = row.names(nodes))

### Calculate internal node ages
nodeage <- cbind(rep(NA, length(branching.times(acaena))), branching.times(acaena))
colnames(nodeage) <- c("spname", "value")


agesTipAndNode <- rbind(ages, nodeage)


###################################################################
###  Dataframe of Clade niche volume & age
###################################################################

volume <- read.csv(".//clade_nicheVolume_acaena.csv")
extractAges <- agesTipAndNode[rownames(agesTipAndNode) %in% volume$node1,]

ageVolData <- cbind(extractAges, volume[order(volume$node1), ])

colnames(ageVolData)[c(1,2,4)] <- c("speciesAge", "spname", "nicheVolume")

write.csv(ageVolData[, c("node1", "spname", "nicheVolume", "speciesAge")], "NicheVolume_age.csv")

###################################################
### Dataframe of Clade niche overlap & phylogenetic distances
###################################################

dis <- data.frame(distance2)
overlap <- read.csv(".//clade_schoennerD_acaena.csv")

overlapPdData <- cbind(overlap, dis[dis$node %in% overlap$node1, ])[ ,c("node1", "node2", 
                                                                        "ecospat.corrected", "distance")]

colnames(overlapPdData)[3:4] <- c("nicheOverlap", "phyloDistance")

# Species name of node
overlapspnames <- (nodes$nodelabel %in% overlapPdData$node1) %>% rownames(nodes)[.]
c(overlapspnames, rep(NA, nrow(overlapPdData) - length(overlapspnames)))

overlapPdData <- mutate(overlapPdData, node1name = 
                          c(overlapspnames, rep(NA, nrow(overlapPdData) - length(overlapspnames)))
)

# Show sister species list
for(i in overlapPdData$node1){
  print(i)
  print(rownames(nodes)[allnodesister[[i]]])
}


write.csv(overlapPdData, "Nicheovrlap_PD.csv")
