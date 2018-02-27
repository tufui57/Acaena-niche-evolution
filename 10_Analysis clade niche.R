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

myplot <- plotAnalysis(data = sisOverlapPd, 
                       yv = "nicheOverlap", xv = "phyloDistance", 
                       nodeNumbercol = "node1", showStats = T,
                       ylabname = "Niche overlap of occurrence records", 
                       xlabname = "Phylogenetic distances between sister species pairs",
                       label.point = TRUE
) +
  ylim(0, 1)

# save
ggsave(paste("Y:\\sister_pd_nicheoverlap_legend.png", sep = ""), plot = myplot,
       width = 300, height = 210, units = 'mm')

rm(myplot)

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

myplot <- plotAnalysis(data = ancsisOverlapPd,
                       yv = "nicheOverlap", xv = "phyloDistance", 
                       nodeNumbercol = "node1", showStats = T,
                       ylabname = "Niche overlap of occurrence records", 
                       xlabname = "Phylogenetic distances between sister species pairs",
                       label.point = TRUE
) +
  ylim(0, 1)

# save
ggsave(paste("Y:\\sisterAunt_pd_nicheoverlap_legend.png", sep = ""), plot = myplot,
       width = 300, height = 210, units = 'mm')

rm(myplot)

#########################################################################
### Clade Phylogenetic distances ~ niche overlap of occurrence records
#########################################################################

### Eliminate outlier, SAC, of clade age

overlapPd <- overlapPdData[which(overlapPdData$node1 != 19), ]

myplot <- plotAnalysis(data = overlapPd,
                       yv = "nicheOverlap", xv = "phyloDistance", 
                       nodeNumbercol = "node1", showStats = T,
                       ylabname = "Niche overlap of occurrence records",
                       xlabname = "Phylogenetic distances between clades",
                       label.point = TRUE
) +
  ylim(0, 1)

# save
ggsave(paste("Y:\\clade_pd_nicheoverlap_legend.png", sep = ""), plot = myplot,
       width = 300, height = 210, units = 'mm')

rm(myplot)


### Leave outlier of clade age in data

myplot <- plotAnalysis(data = overlapPdData, 
                       yv = "nicheOverlap", xv = "phyloDistance",
                       nodeNumbercol = "node1", showStats = T,
                       ylabname = "Niche overlap of occurrence records", 
                       xlabname = "Phylogenetic distances between clades",
                       label.point = TRUE
) +
  ylim(0, 1)

# save
ggsave(paste("Y:\\clade_pd_nicheoverlap_outlier.png", sep = ""), plot = myplot,
       width = 300, height = 210, units = 'mm')

rm(myplot)

#########################################################################
### Clade ages ~ niche volume of occurrence records
#########################################################################

### Eliminate outlier
outlier <- which(max(ageVolData$speciesAge) == ageVolData$speciesAge)

myplot <- plotAnalysis(data = ageVolData[-outlier,],
                       yv = "nicheVolume", xv = "speciesAge", 
                       nodeNumber = "node1", showStats = T,
                       ylabname = "Niche volume of occurrence records",
                       xlabname = "Clade age",
                       label.point = TRUE
) +
  ylim(0, 1)

# save
ggsave(paste("Y:\\clade_age_nicheVolume_legend_acaena.png", sep = ""), plot = myplot,
       width = 300, height = 210, units = 'mm')

rm(myplot)


### Leave outlier
myplot <- plotAnalysis(data = ageVolData,
                       yv = "nicheVolume", xv = "speciesAge", 
                       nodeNumber = "node1", showStats = T,
                       ylabname = "Niche volume of occurrence records",
                       xlabname = "Clade age",
                       label.point = TRUE
) +
  ylim(0, 1)
# save
ggsave(paste("Y:\\clade_age_nicheVolume_outlier.png", sep = ""), plot = myplot,
       width = 300, height = 210, units = 'mm')

rm(myplot)
