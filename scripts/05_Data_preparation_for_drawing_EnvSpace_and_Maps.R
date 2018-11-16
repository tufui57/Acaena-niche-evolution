###############################################################
## Data preparation for drawing maps and PCA
###############################################################

library(dplyr)
library(ggplot2)
library(gridExtra)
library(extrafont)
library(grid)
library(raster)
library(rgdal)

source(".\\functions\\F_plot_map_and_PCA.r")

if(genus_name == "Chionochloa"){
  genus_tag <- "chion"
  
  # data import
  chdata <- read.csv(paste("Y://1st chpater_Acaena project//meta data//", genus_name, "_bioclim_landcover_history_inclNAonland.csv", sep = ""))
}

if(genus_name == "Acaena"){
  genus_tag <- "acaena"
  # data import
  chdata <- read.csv(paste("Y://1st chpater_Acaena project//meta data//", genus_name, "_bioclim_landcover_history_worldclim1_1km18sep.csv", sep = ""))
  
}

chdata <- chdata[is.na(chdata$landCoverChange) == F, ]

# species name
spname <- grepl(genus_name, colnames(chdata)) %>% colnames(chdata)[.]

for(i in spname){
  chdata[is.na(chdata[,i]),i] <- 0
}


# Eliminate species without occurrence records
spname <- (colSums(chdata[, spname]) > 0) %>% spname[.]


# Get PCA axes
pca <- prcomp(chdata[, paste("bioclim", c(1, 6, 12, 15), sep = "")],
              center = TRUE,
              scale. = TRUE)
scores <- data.frame(chdata[, c(colnames(chdata)[grep("^bio", colnames(chdata))], spname,
                           "x", "y", 
                           "landCoverChange", "currentLandcover", "preLandcover")],
                     pca$x[, 1:2])
scores$landCoverChange <- factor(scores$landCoverChange)

# Convert landcover ID to names
scores$pre <- factor(ifelse(scores$preLandcover == 1, "NF", "nonF"))
scores$current <- factor(ifelse(scores$currentLandcover == 1, "NF",
                             ifelse(scores$currentLandcover == 3, "nonF", "non potential habitat" # EF(2) is also non potential habitat
                             ))
)


# reference rasters
ref <- raster("Y://GIS map and Climate data//pre-human_landcover1km.bil")

# Outline of NZ
path = "Y:\\GIS map and Climate data\\lds-nz-coastlines-and-islands-polygons-topo-150k-SHP\\nz-coastlines-and-islands-polygons-topo-150k.shp"
LAYERS <- ogrListLayers(path)
nzland <- readOGR(path, LAYERS)
# Crop extent of polygon
nzland <- crop(nzland, ref)


extent_x = c(min(scores$PC1), max(scores$PC1))
extent_y = c(min(scores$PC2), max(scores$PC2))

