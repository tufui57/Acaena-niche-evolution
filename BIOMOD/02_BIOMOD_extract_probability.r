############################################################################################################
############################         BIOMOD RESULT FORMATTING               ################################
############################################################################################################
##
## Input of this script:
##  "Acaena project\\allNZ_bioclim_sp_landuseAllExtent_modified_test.csv" or "pcaScores_test.data"
##
## Input of this script:
##  ASCII_files\\lcdv1000mcombined.txt
##  BIOMOD\\SPNAME_probability_median.tif
##
## Output of this script:
##  BIOMOD\\Probability_median_over_7types.csv
##  "Acaena project\\allNZ_bioclim_sp_landuseChange_BIOMODprediction.csv"
############################################################################################################

library(raster)
library(rgdal)
library(SDMTools)
library(reshape2)

setwd("Y:\\")

##########################################################################################
####################      BIOMOD RESULT MEDIAN                ############################
##########################################################################################

# prepare reference raster, 1km cell based
post <- raster("ASCII_files\\lcdv1000mcombined.txt")


# projection between WGS84 and NZTM
projectRaster_resample <- function(data, gcs, ref.raster){
  
  # projection to new Geographic Coordinate System (GCS)
  dataP <- projectRaster(data, crs = gcs)
  proj4string(dataP) <- gcs
  
  # resample
  dataPR <- resample(dataP, ref.raster, method = "ngb")
  
  return(dataPR)
}

### prepare extent to change it smaller
e <- extent(165, 180, -50, -32)

# get path of raster data
files <- list.files("BIOMOD")
files2 <- files[grep("tif$", files)]

# make data frame with the same row number as reference raster values
prob <- data.frame(matrix(NA, nrow = 1477950, ncol = 18))
colnames(prob) <- sapply(strsplit(files[grep("tif$", files)], "\\_projection"), "[[", 1)

i <- files2[15]
### raster projection and formatting
for(i in files2[15:16]){
  
  path <- paste("BIOMOD\\", i, sep = "")
  # import raster
  ras <- raster(path)
  
  # crop raster
  ras2 <- crop(ras, e)
  
  # project raster
  projNZ <- proj4string(post)
  ras3 <- projectRaster_resample(ras2, projNZ, post)
  
  coln <- strsplit(i, "\\_projection")[[1]][1]
  # get raster values
  prob[, coln] <- values(ras3)
  
}

# change colnames to match with another data
#da3 <- read.csv("Acaena project\\allNZ_bioclim_sp_landuseAllExtent_modified_test.csv")
colnames(prob) <- colnames(da3[,26:43])

colnames(prob) <- paste(colnames(prob), "_probMedian", sep = "")
  
write.csv(prob,"BIOMOD\\Probability_median_over_7types.csv")

##########################################################################################
#################           BIOMOD RESULT FORMATTING              ########################
##########################################################################################

da3 <- read.csv("Acaena project\\allNZ_bioclim_sp_landuseAllExtent_modified_test.csv")
prob <- read.csv("BIOMOD\\Probability_median_over_7types.csv")

dat <- cbind(da3, prob)


##########################################################################################
######   Check agreement of row order of median probabilities     ########################
##########################################################################################

## If the following shows proper NZ maps, the data row order is right.
plot(raster(matrix(dat$Acaena_anserinifolia_probMedian, nrow = 1475, ncol = 1002)))

## Row order of probMedian columns is wrong. Transpose by the former script for landuse data formatting.

post <- raster("ASCII_files\\lcdv1000mcombined.txt")

trasposeOccurrence <- function(data, spname, ref.raster) {

    if (file.exists(paste("", spname, ".asc", sep = "")) == FALSE) {
        spr <- raster(t(matrix(data[, spname], nrow = 1002, ncol = 1475)))
        extent(spr) <- extent(ref.raster)
        writeRaster(spr, file = paste("", spname, ".asc", sep = ""))
    }

    x <- read.table(paste("", spname, ".asc", sep = ""), skip = 6)
    spcol <- data.frame(melt(x)[, 2])
    colnames(spcol) <- spname
    spcol[spcol[, 1] == -9999,] <- NA
    return(spcol)
}

MedianCols <- colnames(dat)[grep("_probMedian$", colnames(dat))]

dat2 <- lapply(MedianCols,
                  trasposeOccurrence, data = dat, ref.raster = post)

##   Check agreement
plot(raster(matrix(dat2[[1]][,1], nrow = 1475, ncol = 1002)))

dat3 <- do.call(cbind, dat2)
dat[, grep("_probMedian$", colnames(dat))] <- dat3

write.csv(dat, "Acaena project\\allNZ_bioclim_sp_landuseChange_BIOMODprediction.csv")






### Error for Acaena_pallida
### Replace its column with modified column.
dat22 <- read.csv("Acaena project\\allNZ_bioclim_sp_landuseChange_BIOMODprediction.csv")

dat2 <- trasposeOccurrence(spname = MedianCols[15], data = dat, ref.raster = post)
dat22$Acaena_pallida_probMedian <- dat2[,1]

write.csv(dat22, "Acaena project\\allNZ_bioclim_sp_landuseChange_BIOMODprediction3.csv")





