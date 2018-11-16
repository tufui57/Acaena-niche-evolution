##########################################################################################
####################            DATA CLEANING           ##################################
##########################################################################################

### libraries
library(raster)
library(dismo)
library(dplyr)

##########################################################################################
####################           IMPORT DATA           #####################################
##########################################################################################

sp <- read.csv(
  paste(".\\Acaena niche evolution\\Acaena manuscript\\1st chpater_Acaena project\\meta data\\rdata_", genus_name, "_all-presence.csv", sep = "")
)

##########################################################################################
####################   POINT VALIDATION  #################################################
##########################################################################################

sp2<-sp
# ### remove rows having NAs
# sp2 <- sp %>%
#   filter(!is.na( bio1_411 ))

### make sure all NA has been removed
apply(sp2, 2, function(x){sum(is.na(x))})


### in bioclim data, temperature in * 10C and precipitation in mm
### check -> http://www.worldclim.org/formats


### Change data type of taxa
sp2 <- sp2 %>%
  mutate(taxa = as.character(taxa))


### Delete all data out of NZ by coodinate
sp3 <- sp2 %>%
  filter(lat >= -48 & lat <= -33 & lon <= 179 & lon >= 165)

##########################################################################################
####################  NAME CLEANING    ###################################################
##########################################################################################

########################################
# Correct wrong names
########################################

### Check all names in http://www.nzpcn.org.nz/
levels(as.factor(sp3$taxa))

chnames <- read.csv(
  paste("Y:\\1st chpater_Acaena project\\meta data\\NZPCN_Flora_",
                    genus_name, ".csv", sep = "")
  )

chnames2 <- as.character(chnames$Species)
chnames3 <- sort(gsub(" ", "_", chnames2))

# Names of Chionochloa in original data
unique(sp3$taxa)
# Accepted scientific names of Chionochloa among them
sum(unique(sp3$taxa) %in% chnames3 == F)
unique(sp3$taxa)[unique(sp3$taxa) %in% chnames3]

##########################################################################################
####################   Point cleaning     ################################################
##########################################################################################

########################################
# Correct wrong occurrence points
########################################

###### Acaena ######

# Acaena_agnipila
## The data has one point only for Acaena_agnipila_var._aequispina and no points for other 2 var.
## Acaena_agnipila_var._aequispina should be discarded, looking wrong (according to KL).

# Acaena_minor  
## Acaena_minor has 2 var. But these var. have only 6 points (Acaena_minor_var._antarctica has 4 ocean points) in total.
## Acaena_minor_var._minor should be discarded, looking wrong (according to KL).


# Acaena_microphylla (all subsp. was in this sp in the past)
## Acaena_microphylla_var._microphylla distributes in volcano in north island (according to KL).
## Acaena_microphylla_var._pauciglochidiata distributes along coastline in Otago region (according to KL).
## Acaena_rorida is synomym of Acaena microphylla var. microphylla, can be integrated into Acaena_microphylla_var._microphylla.
## discard Acaena_microphylla, but each subsp. is analysed separately.

sp4 <- sp3 %>%
  mutate(taxa = ifelse(grepl("Acaena_rorida",  as.character(taxa)), "Acaena_microphylla_var._microphylla", taxa)) %>%
  filter(taxa != "Acaena_microphylla" )


sp5<-sp4
# ##### remove rows with 0
# sp5 <- sp4 %>% 
#   filter(bio12_411 != 0)

########################################
# Detect wrong occurrences manually 
########################################

### Point detection for removal was done in GIS. Altitude & precipitation can help to find possible mistake.
# Wrong points detected in GIS

# sp5 <- sp4[-which(sp4$FID %in% c(250, 4268, 4398, 4402, 6672, 10108, 7108)), ]

dat <- split(sp5, list(sp5$taxa))

# ########################################
# # Omit outliers in precipitation
# ########################################
# 
# 
# plot.prec <- function(data){lapply(data, function(x){plot(x$bio13_411, x$NZL1_alt, 
#                                                           main = x$taxa[1], 
#                                                           xlab = "Precipitation of Wettest Month",
#                                                           ylab = "Alt"
# )
# }
# )
# }
# 
# # remove outliers
# aca2 <- aca
# 
# plot.prec(aca["Acaena_buchananii"])
# aca2[["Acaena_buchananii"]] <- aca[["Acaena_buchananii"]] %>%
#   filter(bio13_411 < 300)
# 
# plot.prec(aca["Acaena_glabra"])
# aca2[["Acaena_glabra"]] <- aca[["Acaena_glabra"]] %>%
#   filter(NZL1_alt > 400)
# 
# plot.prec(aca["Acaena_profundeincisa"])
# aca2[["Acaena_profundeincisa"]] <- aca[["Acaena_profundeincisa"]] %>%
#   filter((bio13_411 < 500) & (NZL1_alt > 100))                            # In dplyr(C++ format), && doesn't make sense.



# ##### Pass the resultant data
# dat <- aca2

