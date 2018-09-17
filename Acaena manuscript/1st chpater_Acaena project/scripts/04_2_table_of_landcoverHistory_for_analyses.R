############################################################################################################
#########                   Analyses of geographical land cover change                        ###############
############################################################################################################

library(dplyr)
source("Y://R scripts//Scripts for Practice and ones written by someone else//cbind_fill.r")

########################################################
##### Make summary of land cover change history
########################################################

### Summary of pre-human land cover for each sp
dat <- read.csv(paste("Y://", genus_name, "_bioclim_landcover_history_worldclim1_1km17sep.csv", sep=""))
dat <- dat[is.na(dat$landCoverChange) == F, ]

# sp names
spname <- grepl(paste("^", genus_name, sep=""), colnames(dat)) %>% colnames(dat)[.]
# # C. flavicans f. temata, C. nivifera have no occurrence in primary and/or secondary habitats.
# spname <- spname[ - c(17, 21)]

### Create tables of species current/pre-human landcover

r <- list()

for(i in spname){
  # Extract species occurrence column
  spcol <- colnames(dat)[grep("landcover", colnames(dat), ignore.case = TRUE)]
  
  sd <- dat[is.na(dat[,i]) == F, spcol]
  # Table of c(number of occurrences in pre-human landcover, current landcover, unfeasible land cover change)
  prel <- c(nrow(sd), sum(sd[, "preLandcover"] == 1 & sd[, "currentLandcover"] == 4), table(sd[, "preLandcover"]))
  
  if(length(prel) < 4){
    # The species has no records in pre-human forest
    if(names(prel)[3] == "1"){
      prel["2"] <- 0
    }else{
      # The species has no records in pre-human non-forest
      prel["1"] <- 0
    }
  }
  
  names(prel) <- c("number of 1km occurrence cells", "NF-nonPotentialHabitat", "pre-NF", "pre-nonF")
  a <- c(prel, table(sd[, "currentLandcover"]))
  r[[i]] <- a
    }

# Fill NA columns
spsum <- cbind.fill(r)
colnames(spsum) <- spname
rownames(spsum)[5:8] <- c("current-NF", "current-exoticForest","current-nonF", "current-nonPotential")

### Create table of land cover change
landcover <- lapply(spname, function(i){
  spcol <- colnames(dat)[grep("landcover", colnames(dat), ignore.case=TRUE)]
  sd <- dat[is.na(dat[,i]) == F, spcol]
  return(c(nrow(sd), table(sd$landCoverChange)))
}
)

landcover2 <- cbind.fill(landcover)
colnames(landcover2) <- spname

write.csv(t(rbind(spsum, landcover2)), paste(".//", genus_name, "_landscapeChangeHistory.csv", sep=""))

