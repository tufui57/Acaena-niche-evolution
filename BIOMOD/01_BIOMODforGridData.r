############################################################################################################
############################                    BIOMOD                     #################################
############################################################################################################
## 
## This script's references are biomod2 documentation and Matthew's script (00_SDM_read in data_05.r)
## Do NOT run this script on laptops or any slow computers. Processing would be too slow to finish in 24h.
##
## Output of this script:
##  functons in biomod2 package create folders and store all results as default option.
##  See biomod2 documentation for folder path
## 
############################################################################################################


##################################################################################################################
########################        Load libraries
##################################################################################################################

library(rgdal)
library(raster)
library(dplyr)

######################################################################################################
## NOTE! MAKE SURE package "biomod2" is the latest version! IF NOT, INSTALL AGAIN!
## Otherwise, you get errors on model developing!
######################################################################################################
library(biomod2)

# Load functions to develop BIOMOD
source(".\\Acaena niche evolution\\BIOMOD\\Create_Package_BIOMOD.r")

############ NOTE ########################################################################
# BIOMOD_Projection() & BIOMOD_Modeling() write raster files in tempdir()
# MAKE SURE TEMPORARY RASTER FOLDER IS NOT FULL! 
##########################################################################################

# clear temporary raster folder 
tempFilePath <- list.files(rasterOptions()$tmpdir, full.names = T)
file.remove(tempFilePath)


##################################################################################################################
########################        1km grid data import
##################################################################################################################
## Keep working directory same, because all data needed for restoration is saved there.
if(dir.exists("Y://BIOMOD for Grid2") == FALSE){
  dir.create("Y://BIOMOD for Grid2")
}
setwd("Y://BIOMOD for Grid2")

### Set arguments
# data frame of occurrence data and climate data
datapath <- "Y://Acaena project//acaena_bioclim_landcover_1km.csv"
# character string of target genus name
genus_name <- "Acaena"
# Raster of climate data
rasterpath <- "Y:\\GIS map and Climate data\\Acaena_Bioclim_1km_NZTM.data"



# Import data frame of bioclim and occurrence records
climate.occ <- read.csv(datapath)
names(climate.occ)[names(climate.occ) %in% c("x","y")] <- c("NZTMlon", "NZTMlat")

# Raster of bioclim and occurrence records
data.ras <- (load(rasterpath) %>% get)

# Species name
spname <- colnames(climate.occ)[grepl(genus_name, colnames(climate.occ))]

# Remove species with < 1 presence records.
spname <- spname[(climate.occ3[, spname] %>% colSums) > 1]

# Remove NA
climate.occ2 <- (is.na(climate.occ$bioclim1) == F) %>% climate.occ[.,]
climate.occ3 <- (rowSums(climate.occ2[, spname], na.rm = T) > 0) %>% climate.occ2[.,]
for(i in spname){climate.occ3[is.na(climate.occ3[, i]), i] <- 0}

# Environmental variables extracted from BIOCLIM and converted into NZTM.
myExpl <- stack(data.ras[[c("bioclim1", "bioclim6", "bioclim12", "bioclim15")]])

#########################################################
## Takes a while to run... Don't run on laptop! Slow!
#########################################################

try(lapply(spname, runBiomod, data = climate.occ3, myExpl = myExpl, folder.name = "23Feb18"),
    silent = FALSE)

##########################################################
## NOTE! MAKE SURE TEMPORARY RASTER FOLDER IS NOT FULL! ##
##########################################################


#################################################################################################################
########################        PROJECT THE MODELS ONTO CURRENT CLIMATE SPACE
#################################################################################################################

## keep working directory same, because all data needed for restoration is saved there.
setwd("Y://BIOMOD for Grid2")

# get folder names of target genus
folders <- list.dirs(getwd(), full.names = FALSE, recursive = F) %>% grepl(genus_name, .) %>% list.dirs(getwd(), full.names = FALSE, recursive = F)[.]

#########################################################
## Takes a while to run... slow but wait!
#########################################################

lapply(folders, biomodProjection_fromSavedBiomodModels,
       modelname = "23Feb18",
       proj.name = "23Feb18"
       )

#################################################################################################################
########################   Model evaluations     
#################################################################################################################

# myBiomodModelOut <- biomodProjection_fromSavedBiomodModels(
#   folders[1],
#   modelname = "23Feb18",
#   proj.name = "23Feb18"
# )
# get_evaluations(myBiomodModelOut)
# get_calib_lines(myBiomodModelOut)
# get_variables_importance(myBiomodModelOut)


#################################################################################################################
########################        Plot BIOMOD projections
#################################################################################################################

## PLOTS THE PROJECTIONS
setwd("Y:\\BIOMOD for Grid2")
# get folder names
folders <- list.dirs(getwd(), full.names = FALSE, recursive = F) %>% grepl(genus_name, .) %>% list.dirs(getwd(), full.names = FALSE, recursive = F)[.]

projectionPlot <- function(spname, # species name
                           proj.name
){
  # load projection data
  modelname <- load(paste(".\\", spname, "\\proj_", proj.name, "\\", spname, ".", proj.name,  ".projection.out",
                  sep=""))
  model <- get(modelname)
  
  ## plot each projection separately 
  proj_val <- get_predictions(model)
  
  for (i in 1:length(proj_val@layers)) {
    
    png(filename = paste(".\\", spname, "\\proj_", proj.name, "\\", names(subset(proj_val, i)), ".png", sep=""), 
        height = 900, width = 750, units = "px")
    plot(proj_val[[i]])
    title(names(subset(proj_val, i)))
    dev.off()
  }
  
}

lapply(folders, projectionPlot, proj.name = "23Feb18")

#################################################################################################################
########################        BIOMOD ensamble models
#################################################################################################################

# Get folder names
folders <- list.dirs(getwd(), full.names = FALSE, recursive = F) %>% grepl(genus_name, .) %>% list.dirs(getwd(), full.names = FALSE, recursive = F)[.]

ensembleModelling_projection <- function(spname, # species name
                                folder.name,
                                BIOMODproj.name, ensambleProj.name
) {
  # load projection data
  files <- list.files(paste(".//", spname, "//proj_", ensambleProj.name, sep = ""), full.names = T)
  ensamblefile <- (files  %>% grepl("out$", .) %>% files[.])
  
  if(length(ensamblefile) > 0){
    
    print("The ensamble prediction result file already exists.")
  
    }else{
    
    # Load BIOMOD.model.out
    mod <- load(paste(".\\", spname, "\\", spname, ".", folder.name, ".models.out",
                      sep = "")
                )
    model <- get(mod)
    
    ## Ensemble modelling
    myBiomodEM <- BIOMOD_EnsembleModeling(modeling.output = model,
                                          chosen.models = 'all',
                                          em.by = 'all',
                                          eval.metric = c('TSS'),
                                          eval.metric.quality.threshold = c(0.7),
                                          # Models.eval.meth must be one or two, because temporary raster folder can't store data of models for more than 3 evaluation metrics.
                                          # If you run this on computer with bigger storage for the folder,  it may run without the error. (In writeBin(as.vector(v[start:end, ]), x@file@con, size = x@file@dsize) :problem writing to connection)
                                          models.eval.meth = c('TSS'), #, 'ROC', 'ACCURACY'
                                          prob.mean = TRUE,
                                          prob.cv = FALSE,
                                          prob.ci = FALSE,
                                          prob.ci.alpha = 0.05,
                                          prob.median = FALSE,
                                          committee.averaging = FALSE,
                                          prob.mean.weight = TRUE,
                                          prob.mean.weight.decay = 'proportional'
    )
    
    # Load BIOMOD.projection.output
    modelname <- 
      load(paste(".\\", spname, "\\proj_", BIOMODproj.name, "\\", spname, ".", BIOMODproj.name, ".projection.out",
                 sep = "")
      )
    projModel <- get(modelname)
    
    # Creating the ensemble projections
    # Nothing is returned by this function, but specific projection files () are saved on the hard drive projection folder. 
    BIOMOD_EnsembleForecasting(projection.output = projModel,
                                        EM.output = myBiomodEM,
                                        proj.name = ensambleProj.name
    )
    
    
  }
}


lapply(folders, ensembleModelling_projection, 
       folder.name = "23Feb18", BIOMODproj.name = "23Feb18", ensambleProj.name = "24Feb18_ensamble")


## Plot ensemble projection

EMprojectionPlot <- function(spname, # species name
                             proj.name # file location of ensamble model projection.out
) {
  # load projection data
  files <- list.files(paste(".//", spname, "//proj_", proj.name, sep = ""), full.names = T)
  proj <- (files  %>% grepl("grd$", .) %>% files[.] %>% raster)
  
  png(filename = paste(".//ensemble_projection//", names(proj), ".png", sep = ""),
      height = 900, width = 750, units = "px"
  )
  plot(proj)
  title(names(proj))
  dev.off()

}

lapply(folders, EMprojectionPlot, proj.name = "24Feb18_ensamble")
