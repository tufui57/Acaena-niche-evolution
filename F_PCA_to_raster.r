####################################################################
### Function to Rasterize PCA scores of probability
####################################################################

raster_from_dataframe <- function(probdata, # PCA score with probability data
                                  spname, #
                                  rasterSize
){
  # Get probability of the species
  prob1 <- pred[names(pred) == gsub("_",".", i)]
  
  data <- probdata[, c("PC1", "PC2", colnames(probdata)[grep("^prob", colnames(probdata))])]
  colnames(data)[1:2] <- c("x","y")
  
  # Extract target species occurrence records
  points <- data[,c("x", "y")]
  # Set coordinates
  coordinates(points) <- data[,c("x", "y")]
  # Put values in point object
  points$sp <- data[ ,colnames(probdata)[grep("^prob", colnames(probdata))]]
  
  # Prepare empty raster
  r <- raster(ncol = round((max(prob1$PC1) - min(prob1$PC1)) / rasterSize),
              nrow = round((max(prob1$PC2) - min(prob1$PC2)) / rasterSize)
  )
  extent(r) <- extent(points)
  
  # Rasterize
  sp_raster <- rasterize(points, r, points$sp, fun=mean)
  
  return(sp_raster)
}

