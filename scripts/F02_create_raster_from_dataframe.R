
library(raster)
library(rgdal)
library(maptools)


#########################################################################
## Convert dataframe to raster
#########################################################################

convert_dataframe_to_raster <- function(ref, # reference raster
         dat, # dataframe
         coordinateCols, # colnames of coordinates in "dat"
         cellvalueCol # colname of raster cell values in "dat"
){
  # create points
    pts <- dat[, coordinateCols]

    # point coordinate system setting
    coordinates(pts) <- dat[, coordinateCols]
    proj4pts <- proj4string(ref)
    proj4string(pts) <- CRS(proj4pts)
    # Cell value column
    pts$cellVlaue <- dat[, cellvalueCol]

    # rasterize
    prast <- rasterize(pts, ref, field = pts$cellVlaue, fun = mean)

    return(prast)
}


