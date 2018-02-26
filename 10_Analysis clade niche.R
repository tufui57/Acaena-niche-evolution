###################################################
### Clade niche analysis
###################################################

library(phytools)
library(ape)

##############################################################################
### Data preparation
##############################################################################

genus_name <- "Acaena"

source(paste(".//", genus_name, " niche evolution//06_Clade_pairing.R", sep = ""))
source(paste(".//Acaena niche evolution//plotAnalysis_clade_niche.R", sep = ""))

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


########################################
###  Dataframe of Clade niche volume & age
########################################

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
#########################################################################
### Plots
#########################################################################
library(dplyr)

ageVolData <- read.csv("NicheVolume_age.csv")
overlapPdData <- read.csv("Nicheovrlap_PD.csv")

#########################################################################
### Sister species pairs' Phylogenetic distances ~ niche overlap of occurrence records
#########################################################################

### Node numbers of sister species pairs
sispairs <- c(1,20,3,5,16)
sisOverlapPd <- (overlapPdData$node1 %in% sispairs) %>% overlapPdData[., ]

m <- lm(phyloDistance ~ nicheOverlap, sisOverlapPd)

myplot <- plotAnalysis(data = sisOverlapPd, 
                       m = m, 
                       xv = "nicheOverlap", yv = "phyloDistance", 
                       nodeNumber = "node1", showStats = T,
                       xlabname = "Niche overlap of occurrence records", ylabname = "Phylogenetic distances between sister species pairs"
)

# save
ggsave(paste("Y:\\sister_pd_nicheoverlap_legend.png", sep = ""), plot = myplot,
       width = 300, height = 210, units = 'mm')

rm(myplot, m)

#########################################################################
### Sister species pair - their ancestor's sister node
### Phylogenetic distances ~ niche overlap of occurrence records
#########################################################################

ancSisNode <- sapply(sispairs, function(i){
  
  ancestor <- acaena$edge[which(i == acaena$edge[, 2])]
  ancestorSisNode <- distance2[distance2[,"node"] == ancestor, "sisterNode"]
  return(ancestorSisNode)
  
}
)


ancsisOverlapPd <- rbind(
  ((overlapPdData$node1 %in% ancSisNode) %>% overlapPdData[., ]),
  ((overlapPdData$node2 %in% ancSisNode) %>% overlapPdData[., ])
)

dup <- duplicated(ancsisOverlapPd$nicheOverlap) %>% which
ancsisOverlapPd <- ancsisOverlapPd[-dup, ]

m <- lm(phyloDistance ~ nicheOverlap, sisOverlapPd)

myplot <- plotAnalysis(data = sisOverlapPd, 
                       m = m, 
                       xv = "nicheOverlap", yv = "phyloDistance", 
                       nodeNumber = "node1", showStats = T,
                       xlabname = "Niche overlap of occurrence records", ylabname = "Phylogenetic distances between sister species pairs"
)

# save
ggsave(paste("Y:\\sisterAunt_pd_nicheoverlap_legend.png", sep = ""), plot = myplot,
       width = 300, height = 210, units = 'mm')

rm(myplot, m)

#########################################################################
### Clade Phylogenetic distances ~ niche overlap of occurrence records
#########################################################################

### Eliminate outlier, SAC, of clade age

overlapPd <- overlapPdData[which(overlapPdData$node1 != 19), ]
m <- lm(phyloDistance ~ nicheOverlap, overlapPd)

myplot <- plotAnalysis(data = overlapPd, 
                       m = m, 
                       xv = "nicheOverlap", yv = "phyloDistance", 
                       nodeNumber = "node1", showStats = T,
                       xlabname = "Niche overlap of occurrence records", ylabname = "Phylogenetic distances between clades"
)

# save
ggsave(paste("Y:\\clade_pd_nicheoverlap_legend.png", sep = ""), plot = myplot,
       width = 300, height = 210, units = 'mm')

rm(myplot, m)


### Leave outlier of clade age in data
m <- lm(phyloDistance ~ nicheOverlap, overlapPdData)

myplot <- plotAnalysis(data = overlapPdData, m = m, xv = "nicheOverlap", yv = "phyloDistance",
                       nodeNumber = "node1", showStats = T,
                       xlabname = "Niche overlap of occurrence records", ylabname = "Phylogenetic distances between clades")

# save
ggsave(paste("Y:\\clade_pd_nicheoverlap_outlier.png", sep = ""), plot = myplot,
       width = 300, height = 210, units = 'mm')

rm(myplot, m)

#########################################################################
### Clade ages ~ niche volume of occurrence records
#########################################################################

### Eliminate outlier
outlier <- which(max(ageVolData$speciesAge) == ageVolData$speciesAge)
ageVol <- ageVolData[-outlier,]

m <- lm(speciesAge ~ nicheVolume, ageVol)

myplot <- plotAnalysis(data=ageVolData[-outlier,], m=m, xv = "nicheVolume", yv = "speciesAge", 
                       nodeNumber = "node1", showStats = T,
                       xlabname = "Niche volume of occurrence records", ylabname = "Clade age")

# save
ggsave(paste("Y:\\clade_age_nicheVolume_legend_acaena.png", sep = ""), plot = myplot,
       width = 300, height = 210, units = 'mm')

rm(myplot, m)


### Leave outlier
m <- lm(speciesAge ~ nicheVolume, ageVolData)

myplot <- plotAnalysis(data=ageVolData, m=m, xv = "nicheVolume", yv = "speciesAge", 
                       nodeNumber = "node1", showStats = T,
                       xlabname = "Niche volume of occurrence records", ylabname = "Clade age")

# save
ggsave(paste("Y:\\clade_age_nicheVolume_outlier.png", sep = ""), plot = myplot,
       width = 300, height = 210, units = 'mm')

rm(myplot, m)
