############################################################################################################
############################   Get probability from BIOMOD ensemble models
############################################################################################################

setwd("Y://BIOMOD for Grid2")

## Plot ensemble projection
folders <- list.dirs("Y://BIOMOD for Grid2//", full.names = FALSE, recursive = F)
folders <- (grepl("Acaena", folders) %>% folders[.])


get_EMprojection <- function(spname, # species name
                             proj.name # file location of ensamble model projection.out
) {
  # load projection data
  files <- list.files(paste(".//", spname, "//proj_", proj.name, sep = ""), full.names = T)
  proj <- (files  %>% grepl("grd$", .) %>% files[.] %>% raster)
  
  return(proj)
}

pred <- lapply(folders, get_EMprojection, proj.name = "22Feb18_ensamble")
names(pred) <- folders

save(pred, file = "Y://ensemblePrediction.data")
