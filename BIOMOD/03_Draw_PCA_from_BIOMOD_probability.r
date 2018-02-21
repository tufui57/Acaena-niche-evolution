SchoenerD_ecospat.prob <-
  function(background,
           axis1, 
           axis2,
           data1, 
           data2,
           R = 100 # Resolution of background
  ) {
    
    background.clim <- background[, c(axis1, axis2)]
    
    # # calculation of occurence density and test of niche equivalency and similarity
    # z1 <- ecospat::ecospat.grid.clim.dyn(background.clim, background.clim, data1[ ,c(axis1, axis2)], R = 100)
    # z2 <- ecospat::ecospat.grid.clim.dyn(background.clim, background.clim, data2[ ,c(axis1, axis2)], R = 100)
    # 
    res <- list()
    ## Schoener D
    res[[1]] <- unlist(ecospat::ecospat.niche.overlap(z1, z2, cor = T))
    res[[2]] <- unlist(ecospat::ecospat.niche.overlap(z1, z2, cor = F))
    ## Name
    return(res)
  }






### Test data for probability from Acaena anserinifolia
setwd("Y:\\BIOMOD for Grid\\medianRaster\\medianRaster_1km")

list.files(".")
test <- raster(list.files(".")[1])
plot(test)