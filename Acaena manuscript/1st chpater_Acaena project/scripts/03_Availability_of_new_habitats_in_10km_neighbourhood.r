###############################################################################################
### Calculate average availability per primary occurrence in X km neighbourhood
###############################################################################################
# Local availability of secondary open habitat = average size of secondary open area within 10km square of occurrence records in primary open area and forests
# Quantification steps:
# First, for each species, we counted number of 1 km grid cells in secondary open area within 10 km-radius square of neighbourhood area from an occurrence in forest or primary open area. 
# The number of cells is size of available habitat within 10 km-radius neighbourhood of the occurrence.
# We repeated counting size of available habitat within 10 km-radius neighbourhood for all the other occurrence records in forest and primary open area.
# When neighbourhood area of an occurrence record was overlapped with the neighbourhood of another occurrence record, we did not count cells in the overlapped area.

library(dplyr)

genus_name = "Acaena"
dat <- read.csv(paste("Y://", genus_name, "_bioclim_landcover_history_worldclim1_1km13sep.csv", sep=""
))

# Omit NA
dat <- mutate(dat, landCoverChange = as.character(dat$landCoverChange))
dat2 <- dat[is.na(dat$landCoverChange) == FALSE, ]

# ID number of cells
dat$id <- 1:nrow(dat)

# sp names
spname <- colnames(dat)[grepl("^Acaena", colnames(dat))]
# A. emittens and minor have no occurrence in primary or secondary open habitats.

#####################################################################################################################################
### Availability of secondary open habitat within a m radius squire from occurrences in current forest and primary open area
#####################################################################################################################################

Available_SecondaryOpenHabitat <- function(data, # data of points
                                           i, # species number in "spname"
                                           a # half length (unit of the distance is the unit of coordinate) of the squire which you want to count the number of avairable secondary open cells.
){
  
  if(is.factor(data$landCoverChange) == TRUE){
      data <- mutate(data, landCoverChange = as.character(dat$landCoverChange))
      data2 <- data[is.na(data$landCoverChange) == FALSE, ]
  }else{
    data2 <- data[is.na(data$landCoverChange) == FALSE, ]
  }

  
  # Cells having occurrences of the target
  spdat<- data2[!is.na(data2[, spname[i]]), c("x", "y", "landCoverChange")]
  
  # Pre-human occurrence cells; cells in primary open habitat or native forest
  # Occurrence in current exotic forests and native forests that changed from pre-human open habitat (nonF-NF) aren't considered as pre-human occurrences.
  old <- spdat[spdat[, "landCoverChange"] == "nonF-nonF" | spdat[, "landCoverChange"] == "NF-NF", ]
  old$id <- 1:nrow(old)
  
  # Cells in secondary open habitat
  secondary <- data2[data2[, "landCoverChange"] == "NF-nonF", colnames(spdat)]
  
  # Count the number of secondary open cells within squire of the distance "a(m)"
  availSecondary <- lapply(1:nrow(old), function(x){
    secondarylat <- secondary[ (secondary[, "y"] <= old[x, "y"] + a & secondary[, "y"] >= old[x, "y"] - a), ]
    secondarylatlon <- secondarylat[ (secondarylat[, "x"] <= old[x, "x"] + a & secondarylat[, "x"] >= old[x, "x"] - a), ]
    secondarylatlon$oldCellID <- rep(x, nrow(secondarylatlon))
    return(secondarylatlon)
  }
  )
  
  # Merge all list into one data frame
  secondaryInclDuplicated <- do.call(rbind, availSecondary)
  # Remove duplicated cells of secodanry open area
  res <- unique(secondaryInclDuplicated[, c("x","y")])
  
  return(res)
  
  rm(res)
}

### Run for local(10km) and regional(100km) availability

result <- list()
for(i in (1:length(spname))){
  
    result[[i]] <- sapply(c(10000, 100000), function(j){
      # If Available_SecondaryOpenHabitat(dat, i, 10000) doesn't stop within half hour,
      # check if "dat[, c(x,y,landCoverChange)]" has NA.
      Available_SecondaryOpenHabitat(dat, i, j)
  })
}

# Count the number of available secondary open cells for each pre-human occurrence cells
res <- data.frame(cbind(sapply(result, function(x){length(x[[1]])}),
                        sapply(result, function(x){length(x[[3]])})
                        )
                  )
rownames(res) <- spname

# Total number of cells in secondary open habitat
secondary <- dat2[dat2[, "landCoverChange"] == "NF-nonF", ]
res2 <- cbind(res, res[,1]/nrow(secondary), res[,2]/nrow(secondary))
colnames(res2) <- c("10km", "100km",  "10km %", "100km %")
write.csv(res2, "Y://availability_from_forest_primary_10_100km13sep.csv")

