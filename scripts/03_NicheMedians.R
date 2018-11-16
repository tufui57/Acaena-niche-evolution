
library(ggplot2)
library(devtools)

### Import PCA scores
load(paste(".\\Scores_", genus_name, "_landcover18sep.data", sep = ""))


########################################################
### Medians of climate niche across all habitats
########################################################

med1 <- lapply(spname, function(x){
  spscores <- scores[scores[, x] == 1, c("PC1", "PC2")]
  sphist <- hist(spscores$PC1)
  return(median(spscores$PC1))
}
)

med2 <- lapply(spname, function(x){
  spscores <- scores[scores[, x] == 1, c("PC1", "PC2")]
  sphist <- hist(spscores$PC2)
  return(median(spscores$PC2))
}
)

# Combine PC medians
med.data <- cbind(spname, unlist(med1), unlist(med2))
colnames(med.data) <- c("spname", "median.of.temp", "median.of.prec")
write.csv(med.data, file = paste(".//", genus_name, "_nicheBreadths.csv", sep = ""))


########################################################################
### Medians of climate niche in primary vs. secondary occurrences
########################################################################

# Primary occurrences
primaryOcc <- scores[scores$landCoverChange == "NF-NF" | scores$landCoverChange == "nonF-nonF",]
# Secondary occurrences
secondaryOcc <- scores[scores$landCoverChange == "NF-nonF",]

# Medians for primary occurrences
med1.primary <- lapply(spname, function(x){
  spscores <- primaryOcc[primaryOcc[, x] == 1, c("PC1", "PC2")]
  sphist <- hist(spscores$PC1)
  return(median(spscores$PC1))
}
)

med2.primary <- lapply(spname, function(x){
  spscores <- primaryOcc[primaryOcc[, x] == 1, c("PC1", "PC2")]
  sphist <- hist(spscores$PC2)
  return(median(spscores$PC2))
}
)

# Medians for 2ndary occurrences
med1.2ndary <- lapply(spname, function(x){
  spscores <- secondaryOcc[secondaryOcc[, x] == 1, c("PC1", "PC2")]
  sphist <- hist(spscores$PC1)
  return(median(spscores$PC1))
}
)

med2.2ndary<- lapply(spname, function(x){
  spscores <- secondaryOcc[secondaryOcc[, x] == 1, c("PC1", "PC2")]
  sphist <- hist(spscores$PC2)
  return(median(spscores$PC2))
}
)

# Combine PC medians
med.data <- cbind(spname, unlist(med1.primary), unlist(med2.primary), unlist(med1.2ndary), unlist(med2.2ndary))
colnames(med.data) <- c("spname", "primary.temp", "primary.prec", "secondary.temp", "secondary.prec")
write.csv(med.data, file = paste(".//", genus_name, "_nicheMedians.csv", sep = ""))

