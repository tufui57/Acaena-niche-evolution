############################################################################################################
# Convert WORLDCLIM tiff to raster
############################################################################################################

### Extract bioclim data
# Import bioclim 2.5 min (5km) raster data
files <- paste(path, 
               list.files(path),
               sep="\\"
)

### Worldclim 1.4
if(files[grep("bil$", files)] %>% length != 0){
  # Import worldclim
  bio <- lapply(files[grep("bil$", files)], raster)
  
  # Change names to match occurrence data columns
nam <- sapply(sapply(files[grep("bil$", files)], strsplit, "\\\\|.bil"), "[[", 5)
nam2 <- gsub("_411", "", gsub("bio", "bioclim", nam))
names(bio) <- nam2

}else{
  # Import worldclim
  bio <- lapply(files[grep("tif$", files)], raster)
### Worldclim 2
  fil <- sapply(sapply(files[grep("tif$", files)], strsplit, "\\\\|.tif"), "[[", 5)
  nam <- sapply(strsplit(fil, "_"), "[[", 4) %>% gsub("^0", "", .)
  nam2 <- paste("bioclim", nam, sep="")
  names(bio) <- nam2
  
  
}




############################################################################################################
# Project coordinate system of bioclim rasters from WGS84 to NZTM
############################################################################################################

### Change extent to smaller one
# you can use crop() or trim()
e <- extent(165, 180, -50, -32)
bio_cropped <- lapply(bio, crop, e)

### Project raster from WGS84 to NZTM
bioNZ <- lapply(bio_cropped, function(x){
  
  # Projection to new GCS
  dat <- projectRaster(x, crs = proj4stringNZTM)
  
  # Resample raster from 30sec to 1km.
  # 30sec in NZ is NOT 1km.
  dat2 <- resample(dat, ref.raster, method = "ngb")
  return(dat2)
}
)
