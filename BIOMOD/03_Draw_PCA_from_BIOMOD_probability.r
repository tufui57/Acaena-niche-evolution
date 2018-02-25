############################################################################################################
############################   Get Schoenner's D from two groups of probability
############################################################################################################

library(dplyr)
library(dismo)
library(raster)
source(".\\Chionochloa niche evolution\\makeTag.R")
source(".//Acaena niche evolution//plotAnalysis_clade_niche.R")

# Load prediction of ensemble model
load("Y://ensemblePrediction.data")

source(".\\Acaena niche evolution\\06_Clade pairing.R")

spname <- gsub("\\.", "\\_", names(pred)) %>% 
  gsub("var_", "var.", .) %>% 
  gsub("subsp_", "subsp.", .) %>% 
  gsub("novae_", "novae.", .)

code <- makeTag_separate(spname, "Acaena", "_")

# First clade in sister pairs
sispairs <- c(1,20,3,5,16)

probD <- list()
for(i in sispairs){
  # Species name codes
  nodeName = pull(code[code$X %in% rownames(nodes)[i], ], X)
  sisnodeName = pull(code[code$X %in% rownames(nodes)[allnodesister[[i]]], ], X)
  
  probANS <- (spname == nodeName) %>% pred[.]
  probDUM <- (spname == sisnodeName) %>% pred[.]
  ### Use dismo::nicheOverlap
  probD[[i]] <- nicheOverlap(probANS[[1]], probDUM[[1]], stat = 'D', mask = TRUE, checkNegatives = TRUE) 
  
}



############################################################################################################
##### Compare Schoenner's D from two groups of probability and the one from occurrence records
############################################################################################################

# Import data
overlapPdData <- read.csv("Nicheovrlap_PD.csv")

### Node numbers of sister species pairs
sispairs <- c(1,20,3,5,16)
sisOverlapPd <- (overlapPdData$node1 %in% sispairs) %>% overlapPdData[., ]

overlaps <- cbind(sisOverlapPd, unlist(probD))
colnames(overlaps)[ncol(overlaps)] <- "probD"


m <- lm(probD ~ nicheOverlap, overlaps)
myplot <- plotAnalysis(data = overlaps, 
                       m = m, 
                       xv = "nicheOverlap", yv = "probD", 
                       nodeNumber = "node1", showStats = T,
                       xlabname = "Niche overlap of occurrence records", ylabname = "Niche overlap of model prediction"
)

# save
ggsave(paste("Y:\\sister_nicheoverlap_acaena.png", sep = ""), plot = myplot,
       width = 300, height = 210, units = 'mm')

rm(myplot, m)


#########################################################################
### Sister species pairs' Phylogenetic distances ~ niche overlap of predictions
#########################################################################

m <- lm(probD ~ phyloDistance, overlaps)
myplot <- plotAnalysis(data = overlaps, 
                       m = m, 
                       xv = "phyloDistance", yv = "probD", 
                       nodeNumber = "node1", showStats = T,
                       xlabname = "Phylogenetic distances", ylabname = "Niche overlap of model prediction"
)

# save
ggsave(paste("Y:\\sister_predNicheoverlap_phyloDistance_acaena.png", sep = ""), plot = myplot,
       width = 300, height = 210, units = 'mm')

rm(myplot, m)

