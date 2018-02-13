########################################################################################
### Quantification of niche volume
########################################################################################
### Quantify niche volume as niche overlap with imaginary species having occurrence records in all grid cells across NZ

########################################################################################
### Acaena
########################################################################################

# Data import
da1 <- read.csv("Y:\\acaena_bioclim_landcover_history_inclNAonland.csv")
d <- da1[is.na(da1$landCoverChange) == F, ]

# Create species distriobution having full occurrences across NZ
d$Acaena_allNZ <- 1

# sp names
sname <- colnames(d)[grepl("^Acaena", colnames(d))]
# Acaena emittens and A. minor have no occurrence in primary and/or secondary habitats.
s <- sname[ - c(6, 13)]

for(i in sname){
  d[is.na(d[,i]),i] <- 0
}

# get env. corrdinates (PCA axes)
pca <- prcomp(d[, paste("bioclim", c(1, 6, 12, 15), sep = "")],
              center = TRUE,
              scale. = TRUE)
scores <- data.frame(d[, c(colnames(d)[grep("^bioclim", colnames(d))], sname, "Acaena_allNZ", "x", "y")], pca$x[, 1:2])

### Original method in ecospat to calculate Schonner's D
library(ecospat)

SchoenerD <- function(scores, gen1, gen2) {
  
  # Extract data of two target species
  scores.gen1 <-  scores[scores[, gen1] == 1,c("PC1", "PC2")]
  scores.gen2 <-  scores[scores[, gen2] == 1,c("PC1", "PC2")]
  scores.clim <- scores[, c("PC1", "PC2")]
  
  # calculation of occurence density and test of niche equivalency and similarity
  z1 <- ecospat.grid.clim.dyn(scores.clim, scores.clim, scores.gen1, R = 100)
  z2 <- ecospat.grid.clim.dyn(scores.clim, scores.clim, scores.gen2, R = 100)
  
  res <- list()
  ## Schoener D
  res[[1]] <- unlist(ecospat.niche.overlap(z1, z2, cor = T))
  res[[2]] <- unlist(ecospat.niche.overlap(z1, z2, cor = F))
  # Name
  name_genera <- paste(strsplit(gen1, ".csv"),strsplit(gen2, ".csv"), sep="_")
  names(res) <- c(paste(name_genera, "corrected"), paste(name_genera, "not corrected"))
  
  return(res)
}

D <- lapply(sname[-19], SchoenerD, scores=scores, gen2="Acaena_allNZ")
D2 <- unlist(D)

write.csv(D2, "Y:\\acaena_niche_volume.csv")

# Collate data
d <- read.csv("Y:\\acaena_niche_volume.csv")
d2 <- d[grepl("allNZ corrected.D$", d$X),"x"]
d3 <- d[grepl("allNZ not corrected.D$", d$X),"x"]
d4 <- d[grepl("allNZ corrected.I$", d$X),"x"]
d5 <- d[grepl("allNZ not corrected.I$", d$X),"x"]
d6 <- data.frame(cbind(d2,d3,d4,d5))
colnames(d6) <- c("corrected.D","not corrected.D","corrected.I","not corrected.I")
rownames(d6) <- gsub("_Acaena_allNZ.*","", d[grepl("allNZ corrected.D$", d$X),"X"])
write.csv(d6, "Y:\\Niche change of lineages\\nicheVolume.csv")

########################################################################################
### Quantification of niche breadths
########################################################################################

scores.r <- data.frame(d[, c(colnames(d)[grep("^bioclim", colnames(d))], sname, "landCoverChange", "x", "y")], pca$x[, 1:2])

niche.breadth <- function(gen1) {
  
  # Current niche breadths
  scores.gen1 <-  scores.r[scores.r[, gen1] == 1, c("landCoverChange","PC1", "PC2")]
  range.pc1 <- (max(scores.gen1$PC1) - min(scores.gen1$PC1))
  range.pc2 <- (max(scores.gen1$PC2) - min(scores.gen1$PC2))
  current <- c(range.pc1, range.pc2)
  
  # Niche breadths within Primary open habitat
  scores.1 <-  scores.gen1[scores.gen1[, "landCoverChange"] == "nonF-nonF", c("PC1", "PC2")]
  p.pc1 <- (max(scores.1$PC1) - min(scores.1$PC1))
  p.pc2 <- (max(scores.1$PC2) - min(scores.1$PC2))
  primary <- c(p.pc1, p.pc2)
  
  # Niche breadths within Secondary open habitat
  scores.2 <-  scores.gen1[scores.gen1[, "landCoverChange"] == "NF-nonF", c("PC1", "PC2")]
  s.pc1 <- (max(scores.2$PC1) - min(scores.2$PC1))
  s.pc2 <- (max(scores.2$PC2) - min(scores.2$PC2))
  secondary <- c(s.pc1, s.pc2)
  
  res <- rbind(current, primary, secondary)
  res<-data.frame(res)
  colnames(res) <- c("PC1", "PC2")
  return(res)
}

res <- lapply(sname[-19], niche.breadth)
names(res) <- sname[-19]
res2 <- do.call(cbind, res)

write.csv(res2, "Y:\\acaena_niche_breadth.csv")
