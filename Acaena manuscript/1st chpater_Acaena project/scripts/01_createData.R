############################################################################################################
# Create land cover change history and species occurrence data
############################################################################################################

############################################################################################################
### NOTE; it's impossible to create the exact same landcover data using R as those using ArcGIS,
### because R doesn't have the same fucntion to assign raster values as Arc.
### Current manuscript (17Sep2018) is based on data using Arc.
############################################################################################################

############################################################################################################
# Preparation
############################################################################################################

library(raster)
library(rgdal)
library(SDMTools)
library(dplyr)

# Give genus name
genus_name <- "Acaena"

# Give raster resolution (km)
reso <- 1

Worldclim <- 1
# Give path of WORLDCLIM raster
if(Worldclim == 1){

  # Worldclim ver.1.4
  path <-"Y:\\GIS map and Climate data\\worldclim\\bio_411"
}
if(Worldclim == 2){
 # Worldclim ver.2
  path <- "Y:\\GIS map and Climate data\\worldclim\\wc2.0_2.5m_bio"
}

setwd("C:\\Users\\nomur\\Documents")

############################################################################################################
# Resample pre-human land cover raster from resolution 100m to x km
############################################################################################################
### Create a reference raster with 1km2 resolution and NZTM Geographic Coordinate System (GCS)

if(file.exists(paste("Y:\\GIS map and Climate data\\pre-human_landcover", reso, "km.bil", sep="")) == FALSE){
  
  # Import (100m)^2 raster data
  pre <- raster("Y:\\GIS map and Climate data\\lris-new-zealand-potential-vegetation-grid-version-GTiff\\new-zealand-potential-vegetation-grid-version.tif")
  
  # Resample raster to x km
  
  # You need assign the raster value occupying the maximum area within a cell.
  # Use "fun = modal" which gives "mode" (the most common value).
  ref.raster <- aggregate(pre, fun = modal, na.rm = T, fact = reso*10)
  writeRaster(ref.raster, paste("Y:\\GIS map and Climate data\\pre-human_landcover", reso, "km.bil", sep=""), format = "EHdr")
  
}

# WARNING; the following reference raster doesn't have enough raster values when you run this code in computers with insufficient disc space.
# But it works as a reference raster to give reference extent and coordination system.

ref.raster <- raster(
  # paste("Y:\\GIS map and Climate data\\pre-human_landcover", reso,"km.bil", sep=""
  #       )
  
  # This raster was converted from pre-human raster with ArcGIS using "majority" method for raster value assignment
  "Y://GIS map and Climate data//newzealandpotentialvegetatio1.bil"
)

############################################################################################################
# Resample current land cover polygon to raster at resolution 5km
############################################################################################################

pre <- ref.raster

current_ras <- raster(
  # paste("Y:\\GIS map and Climate data\\current_landcover", reso, "km.bil", sep=""
  # )
  
  # This raster was converted from current land cover polygon with ArcGIS using "maximum area" method for raster value assignment
  "Y://GIS map and Climate data//lcdbv41landcoverdata2.bil"
)

# Resample to adjust extent and number of cells
current_ras2  <- resample(current_ras, pre, method = "ngb")

############################################################################################################
# Convert WORLDCLIM tiff to raster
# Project coordinate system of bioclim rasters from WGS84 to NZTM
############################################################################################################

proj4stringNZTM <- proj4string(current_ras)

source(".\\functions\\F01_project_resample_WORLDCLIM.R")

############################################################################################################
# Resample species occurrence raster on 1km to 5km resolution 
############################################################################################################

# Import species occurrence data
source(".\\functions\\F02_clean_up_species_records13Sep.R")

### Get rid of species whose occurren resords < 5
dat2 <- dat[sapply(dat, nrow) >= 5]

# Import reference raster of WGS
refWGS <- raster("Y:\\GIS map and Climate data\\worldclim\\bio_411\\bio1_411.bil")

# Import function
source(".//functions//F04_convert_occurrencePoints_to_raster.R")

spRaster <- lapply(dat2, project_and_convert_occurrencePoints_to_raster, refWGS = pre, val = "occurrence")


############################################################################################################
# Convert raster to dataframe
############################################################################################################

# Make raster stack
bio_land <- stack(c(bioNZ, spRaster, pre, current_ras2))

source(".//Acaena niche evolution//Acaena manuscript//1st chpater_Acaena project//scripts//02_elevation.R")

res <- data.frame(cbind(coordinates(bio_land2[[1]]), values(bio_land2)))

############################################################################################################
# Create land cover history
############################################################################################################

source(".//functions//F03_convert_landcoverChange.R")

d <- landCoverChange(res, prehuman_landcover="layer.1", current_landcover="layer.2")

# Delete NA rows of BIOCLIM data
d2 <- d[!is.na(d$bioclim1), ]

write.csv(d2, file = paste("Y://", genus_name, "_bioclim_landcover_history_worldclim",
                           Worldclim, "_", reso, "km17sep.csv", sep=""
                           )
          )
