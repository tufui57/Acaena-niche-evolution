##############################################################################################
###     Draw maps for species pairs
##############################################################################################

library(ggplot2)
library(raster)
library(rgdal)
library(maptools)
library(gridExtra)

# Data import
d1 <- read.csv("Y://acaena_bioclim_landcover_history_inclNAonland.csv")

### reference rasters
ref <- raster("Y://GIS map and Climate data//current_landcover1km.bil")

extent_x = extent(ref)[1:2]
extent_y = extent(ref)[3:4]

# sp name
spname <- colnames(d1)[grep("^Acaena", colnames(d1))]

# Outline of NZ
path = "Y:\\GIS map and Climate data\\lds-nz-coastlines-and-islands-polygons-topo-150k-SHP\\nz-coastlines-and-islands-polygons-topo-150k.shp"
LAYERS <- ogrListLayers(path)
nzland <- readOGR(path, LAYERS)

# Next the shapefile has to be converted to a dataframe for use in ggplot2
shapefile_df <- fortify(nzland)

###############################################################
## Point map for land use change of sp
###############################################################

spPairMap <- function(data, # occurrence record data
                         spnames) {
  
  spdata1 <- data[which(data[, spnames[1]] == 1),]
  spdata2 <- data[which(data[, spnames[2]] == 1),]
  
  pMain <- ggplot() +
    # plot all NZ data points
    geom_path(data = shapefile_df, 
              aes(x = long, y = lat, group = group),
              color = 'gray', fill = 'white', size = .2) +
    # point of each sp
    geom_point(data = spdata1, aes(x = x, y = y), color = "green", alpha = 0.3) +
    # point of each sp
    geom_point(data = spdata2, aes(x = x, y = y), color = "purple", alpha = 0.3) +
    # extent
    xlim(extent_x) +
    ylim(extent_y) +
    guides(colour = guide_legend(override.aes = list(size = 5, shape=16, alpha=0.7))) +
    # legend position inside plot
    theme(legend.justification = c(1, 1), legend.position = c(1, 1),
          panel.background = element_rect(fill = 'gray96')
    )
  
  pRightTop <- ggplot() +
    # plot all NZ data points
    geom_path(data = shapefile_df, 
              aes(x = long, y = lat, group = group),
              color = 'gray', fill = 'white', size = .2) +
    # extent
    xlim(extent_x) +
    ylim(extent_y) +
    # point of each sp
    geom_point(data = spdata1, aes(x = x, y = y), color = "green", alpha = 0.3) +
    ggtitle(spnames[1]) +
    theme(panel.background = element_rect(fill = 'gray96'))
  
  
  pRightBottom <- ggplot() +
    # plot all NZ data points
    geom_path(data = shapefile_df, 
              aes(x = long, y = lat, group = group),
              color = 'gray', fill = 'white', size = .2) +
    # extent
    xlim(extent_x) +
    ylim(extent_y) +
    # point of each sp
    geom_point(data = spdata2, aes(x = x, y = y), color = "purple", alpha = 0.3)+
    ggtitle(spnames[2]) +
    theme(panel.background = element_rect(fill = 'gray96'))
  
  ### Plot tow species niche on one figure 
  png(filename = paste("Y:\\map_", spnames[1], spnames[2], ".png"), width = 900, height = 1300)
  plot(pMain)
  dev.off()
  
  ### Plot two species niche separately
  png(filename = paste("Y:\\map_", spnames[1], spnames[2], "_separate.png"), width = 450, height = 1300)
  # change font size
  theme_set(theme_gray(base_size = 18))
  # Plot in multiple panels
  grid.arrange(pRightTop, pRightBottom, ncol = 1, nrow = 2)
  dev.off()
  
}

spPairMap(d1,  c("Acaena_tesca", "Acaena_buchananii"))
