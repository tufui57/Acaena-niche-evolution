##############################################################################################
###     Draw maps for each sp coloured by land use change     ###
##############################################################################################

library(raster)
library(rgdal)
library(maptools)

# Data import
d1 <- read.csv(paste("Y://1st chpater_Acaena project//meta data//", genus_name, "_bioclim_landcover_history_worldclim1_1km18sep.csv", sep = ""))

# sp name
spname <- colnames(d1)[grep("^Acaena", colnames(d1))]

# Reference rasters
ref <- raster("Y://GIS map and Climate data//current_landcover1km.bil")

# Outline of NZ
path = "Y:\\GIS map and Climate data\\lds-nz-coastlines-and-islands-polygons-topo-150k-SHP\\nz-coastlines-and-islands-polygons-topo-150k.shp"
LAYERS <- ogrListLayers(path)
nzland <- readOGR(path, LAYERS)

# ###############################################################
# ## Point map for land use change of sp 
# ###############################################################
# 
# pointPlot_sp <- function(d, # subset data for a species including land use change column 
#                          title) {
#   
#   d <- d[!is.na(d[, "landCoverChange"]),]
#   d <- d[ - which(d$landCoverChange == "NF-nonPotentialHabitat"| d$landCoverChange == "nonF-nonPotentialHabitat"| d$landCoverChange == "NF-EF" | d$landCoverChange == "NF-NF" | d$landCoverChange == "nonF-EF" | d$landCoverChange == "nonF-NF"),]
#   
#   # Convert land use change column to numeric
#   d$changeNo <- NA
#   d[d$landCoverChange == "nonF-nonF", "changeNo"] <- 1
#   d[d$landCoverChange == "NF-nonF", "changeNo"] <- 2
#   
#   # create point
#   pts <- d[, c("x", "y")]
#   
#   # point coordinate system setting
#   coordinates(pts) <- d[, c("x", "y")]
#   proj4pts <- proj4string(ref)
#   proj4string(pts) <- CRS(proj4pts)
#   # land use change column
#   pts$changeNo <- d$changeNo
#   
#   #####################
#   # Plot
#   #####################
#   png(filename = paste("Y:\\landCoverChange_", title, ".png", sep = ""),
#       width = 500, height = 710)
#   par(cex = 0.8)
#   plot(pts,
#        # colour by group and add alpha
#        bg = c(rgb(0, 0, 1, 0.5), rgb(1, 0, 0, 0.5))[pts$changeNo],
#        pch = 21, axes = T,
#        # no outline
#        col = NA,
#        main = paste(title, "\n Number of occurrence cells = ", nrow(d)),
#        cex.main = 1.8,
#        # set extent
#        xlim = extent(ref)[1:2], ylim = extent(ref)[3:4]
#   )
#   # add legend manually
#   legend(x = 1100000, y = 6200000,
#          pch = 21, col = c("red", "blue"), legend = c("Secondary\n open area", "Primary\n open area"),
#          cex = 1.5)
#   
#   # Outline of NZ
#   plot(nzland, add = TRUE)
#   
#   dev.off()
# }
# 
# ## all sp
# # extract sp data
# dall <- d1[rowSums(d1[, spname]) > -18,]
# pointPlot_sp(dall, "All studied sepcies of Acaena")
# 
# ###############################################################
# ## Point map for land use change of sp 
# ###############################################################
# 
# d2 <- lapply(spname, function(i){d1[d1[,i]==1,]})
# lapply(1:length(spname), function(i){pointPlot_sp(d2[[i]], title=paste(spname[i], "2ndCategories", sep="_"))})
#   
#   

#########################################################################
## Raster plotting for Pre/post-human forest/non-forest diatribution
#########################################################################

# create new raster having the same dimentions as reference raster (ex. pre-human map)
rast <- raster(ncol = ncol(ref), nrow = nrow(ref))
extent(rast) <- extent(ref)

prepost <- function(d, title, legendLabels, saveImage = T ) {
    # create point
    pts <- d[, c("x", "y")]

    # point coordinate system setting
    coordinates(pts) <- d[, c("x", "y")]
    proj4pts <- proj4string(ref)
    proj4string(pts) <- CRS(proj4pts)
    # land use change column
    pts$changeNo <- d$changeNo

    # rasterize
    prast <- rasterize(pts, ref, field = pts$changeNo, fun = mean)

    # plot
    if(saveImage==T){
      
      png(filename = paste("Y:\\", title, ".png", sep = ""),
          width = 550, height = 710)
      
      plot(prast, col = c("green", "brown"),
           main = paste(title),
           axes = F, box=FALSE, cex.main = 4,
           legend = FALSE
      )
      # add legend manually
      legend(x = 1100000, y = 6200000, cex = 3.25,
             fill = c("green", "brown", "white"), legend = legendLabels)
      
      # Outline of NZ
      plot(nzland, add = TRUE)
      
      dev.off()
    }else{
      plot(prast, col = c("green", "brown"),
           main = paste(title),
           axes = F, box=FALSE, cex.main = 4,
           legend = FALSE
      )
      # add legend manually
      legend(x = 1100000, y = 6200000, cex = 3.25,
             fill = c("green", "brown", "white"), legend = legendLabels)
      
      # Outline of NZ
      plot(nzland, add = TRUE)
    }

    rm(prast, pts, d)
}

##############################
## Pre-human
##############################
d <- d1[!is.na(d1[, "preLandcover"]),]
d$changeNo <- d$preLandcover

prepost(d, "Pre-human land cover", c("Native Forest", "Non-Forest", "Others"))

##############################
## Current
##############################
d2 <- d1[!is.na(d1[, "currentLandcover"]),]
d3 <- d2[ - which(d2$currentLandcover == 2 | d2$currentLandcover == 4),]
d3$changeNo <- d3$currentLandcover

prepost(d3, "Current land cover",  c("Native Forest", "Non-Forest", "Others"))

##############################
## All sp land use change
##############################
dall <- d1[rowSums(d1[, spname], na.rm = T) > 0,]
dall <- dall[!is.na(dall[, "landCoverChange"]),]
dall <- dall[ - which(dall$landCoverChange == "NF-nonPotentialHabitat"| dall$landCoverChange == "nonF-nonPotentialHabitat"| dall$landCoverChange == "NF-EF" | dall$landCoverChange == "NF-NF" | dall$landCoverChange == "nonF-EF" | dall$landCoverChange == "nonF-NF"),]

# Convert land use change column to numeric
dall$changeNo <- NA
dall[dall$landCoverChange == "NF-nonF", "changeNo"] <- 1
dall[dall$landCoverChange == "nonF-nonF", "changeNo"] <- 2

landuseRasterMap <- function(dall, title, saveImage = T) {
# create point
  pts <- dall[, c("x", "y")]
  
  # point coordinate system setting
  coordinates(pts) <- dall[, c("x", "y")]
  proj4pts <- proj4string(ref)
  proj4string(pts) <- CRS(proj4pts)
  # land use change column
  pts$changeNo <- dall$changeNo
  
  # rasterize
  prast <- rasterize(pts, rast, field = pts$changeNo, fun = mean)
  
  # plot
  if(saveImage==T){
    
  png(filename = paste("Y:\\", title, ".png", sep = ""),
      width = 550, height = 710)
  
  plot(prast, col = c("red", "blue"),
       axes = F, box = FALSE, cex.main = 4,
       main = paste(title),
       legend = FALSE
  )
  # add legend manually
  legend(x = 1100000, y = 6200000, cex = 3.25,
         fill = c("red", "blue", "white"), legend = c("Secondary\n open area", "Primary\n open area", "Others"))
  
  # Outline of NZ
  plot(nzland, add = TRUE)
  
  dev.off()
  }else{
    plot(prast, col = c("red", "blue"),
         axes = F, box = FALSE, cex.main = 4,
         main = paste(title),
         legend = FALSE
    )
    # add legend manually
    legend(x = 1100000, y = 6200000, cex = 3.25,
           fill = c("red", "blue", "white"), legend = c("Secondary\n open area", "Primary\n open area", "Others"))
    
    # Outline of NZ
    plot(nzland, add = TRUE)
  }
}

# landuseRasterMap(dall, "Land use change of Acaena")

##############################
## Land use change of NZ
##############################
dland <- d1[!is.na(d1[, "landCoverChange"]),]
dland <- dland[ - which(dland$landCoverChange == "NF-nonPotentialHabitat"| dland$landCoverChange == "nonF-nonPotentialHabitat"| dland$landCoverChange == "NF-EF" | dland$landCoverChange == "NF-NF" | dland$landCoverChange == "nonF-EF" | dland$landCoverChange == "nonF-NF"),]

# Convert land use change column to numeric
dland$changeNo <- NA
dland[dland$landCoverChange == "NF-nonF", "changeNo"] <- 1
dland[dland$landCoverChange == "nonF-nonF", "changeNo"] <- 2

# landuseRasterMap(dland, "Land cover change")


##########################################################################################
## Three Maps in one panel; pre-/post-land use and Land use change of NZ
##########################################################################################
png("Y://three_maps.png", width = 1600, height = 800)

par(mfrow=c(1,3))
prepost(d, "Pre-human land cover", c("Native Forest", "Non-Forest", "Others"), saveImage = F)
prepost(d3, "Current land cover",  c("Native Forest", "Non-Forest", "Others"), saveImage = F)
landuseRasterMap(dland, "Land cover change", saveImage = F)

dev.off()

####################################
## All sp land use change (point)
####################################
pointPlot_sp <- function(d, title) {
  
  # create point
  pts <- d[, c("x", "y")]
  
  # point coordinate system setting
  coordinates(pts) <- d[, c("x", "y")]
  proj4pts <- proj4string(ref)
  proj4string(pts) <- CRS(proj4pts)
  # land use change column
  pts$changeNo <- d$changeNo
  
  #####################
  # Plot
  #####################
  png(filename = paste("Y:\\landuseChange_", title, ".png", sep = ""),
      width = 550, height = 710)
  par(cex = 0.8)
  plot(pts,
       # colour by group and add alpha
       bg = c(rgb(0, 0, 1, 0.5), rgb(1, 0, 0, 0.5))[pts$changeNo],
       pch = 21, axes = TRUE,
       # no outline
       col = NA,
       main = paste(title, "\n total N = ", nrow(d)),
       # set extent
       xlim = extent(ref)[1:2], ylim = extent(ref)[3:4]
  )
  # add legend manually
  legend(x = 1200000, y = 6200000, pch = 21, col = c("red", "blue"), legend = c("NF-nonF", "nonF-nonF"), cex = 0.8)
  
  # Outline of NZ
  plot(nzland, add = TRUE)
  
  dev.off()
}

pointPlot_sp(dall, "Land cover change of all studied Acaena")
