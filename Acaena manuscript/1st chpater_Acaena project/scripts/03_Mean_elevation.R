############################################################################################################
# Elevation data
############################################################################################################

# Import data
dat <- read.csv(paste("Y://", genus_name, "_bioclim_landcover_history_worldclim1_1km14sep.csv", sep="")
         )

spname <- colnames(dat)[grepl("^Acaena", colnames(dat))]

ele <- sapply(spname, function(x){
  sprow <- which(dat[, x] == 1)
  alt <- mean(dat[sprow, "value"])
  return(alt)
  }
)

write.csv(ele, ".//elevation.csv")
