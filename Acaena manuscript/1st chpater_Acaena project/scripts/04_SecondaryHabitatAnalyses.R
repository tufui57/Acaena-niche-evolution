########################################################
### Regression figures for ratio of 2ndary open habitat 
########################################################

###################################################################################
### Reproducible data set for analyses (on construction at 3.July.2018) 
###################################################################################

genus_name = "Acaena"
genus_tag = "acaena"

# Import data
if(
  file.exists(paste(".//Acaena niche evolution//Acaena manuscript//1st chpater_Acaena project//meta data//", genus_tag, "_data_analyses.csv", sep = ""))
){
  d <- read.csv(paste(".//Acaena niche evolution//Acaena manuscript//1st chpater_Acaena project//meta data//", genus_tag, "_data_analyses.csv", sep = ""))
}else{
  source(".//Acaena niche evolution//Acaena manuscript//1st chpater_Acaena project//scripts//04_2_table_of_landcoverHistory_for_analyses.R")
  source(".//Acaena niche evolution//Acaena manuscript//1st chpater_Acaena project//scripts//04_2_calculate_indices_for_analyses.R")
  d <- read.csv(paste(".//", genus_name, "_data_analyses.csv", sep = ""))
}

library(dplyr)
source(".//functions//F_speciseNameCleaning_spnameFromPhylogenyTree.r")
d$tag <- makeTag_separate(d$spname, genus_name, "_")
d$tag <- d$tag$tag %>% toupper

library(ggplot2)
source(".//functions//F_plotAnalysis_clade_niche.R")

### Merge indices
breath <- read.csv(paste(".//", genus_name, "_nicheBreadths.csv", sep = ""))
trait <- read.csv(".//Acaena niche evolution//Acaena manuscript//1st chpater_Acaena project//raw data//traits.csv")

# combination of dispersal type and life form
trait$type <- paste(trait$barbSection, trait$lifeForm, sep="_")
trait$type[trait$type != "Ancistrum_Stoloniferous" & trait$type != "Microphyllae_Rhizomatous"] <- "other"
trait$type <- as.factor(trait$type)

schoenerD <- read.csv(".//shoenerD17sep_r100.csv")
### Does NA in Schoener's D mean no overlap? ###
schoenerD$corrected.D[is.na(schoenerD$corrected.D)] <- 0

avail <- read.csv(".//availability_from_forest_primary_10_100km17sep.csv")

nicheVol <- read.csv(".//Acaena_nicheVolume.csv")
colnames(nicheVol)[2] <- "niche_volume"

meanEle <- read.csv(".//elevation.csv")
colnames(meanEle)[2] <- "alt"


dat <- merge(d, breath) %>% 
  merge(., trait, by.x = "spname", by.y = "X") %>% 
  merge(., schoenerD, by.x = "spname", by.y = "X") %>%
  merge(., avail, by.x = "spname", by.y = "X") %>% 
  merge(., meanEle, by.x = "spname", by.y = "X") %>%
  merge(., nicheVol, by.x = "spname", by.y = "X")

write.csv(dat, file = ".//Acaena niche evolution//Acaena manuscript//1st chpater_Acaena project//meta data//Acaena_data_analyses17sep.csv")

#################################################################################
### Use original data set (but unreproducible) for regression figures
#################################################################################
# dat <- read.csv("Y://1st chpater_Acaena project//meta data//Acaena_elevation.csv")

dat <- read.csv(".//Acaena niche evolution//Acaena manuscript//1st chpater_Acaena project//meta data//Acaena_data_analyses14sep.csv")


summary(glm(proportionSecondaryHabitat ~ log10.total + 
              alt + X10km.. + PreferenceOpen +
              niche_volume + corrected.D + 
              type +
              median.of.prec + median.of.temp,
            data = dat)
)


##########################################################################
###  Proportion of Secondary Open Habitat - current range size
##########################################################################

myplot <- plotAnalysis(dat, genus_name = genus_name,
                       xv = "log10.current_range", yv = "proportionSecondaryHabitat", showStats = T,
                       xlabname = "log10(Current range size)", ylabname = "Proportion of secondary open habitat",
                       nodeNumbercol = "tag", label.point = T
                       )
  

# save
ggsave(paste("Y:\\log10CurrentRangeSize_proportionSecondary.png", sep = ""), plot = myplot,
       width = 200, height = 140, units = 'mm')

rm(myplot)


#########################################################################
### Proportion of Secondary Open Habitat ~ niche volume
#########################################################################

myplot <- plotAnalysis(data = dat, genus_name = genus_name,
                       xv = "niche_volume", yv = "proportionSecondaryHabitat",  showStats = T,
                       xlabname = "Species niche volume", ylabname = "Proportion of secondary open habitat",
                       nodeNumbercol = "tag", label.point = T
                       )
# Save
ggsave(paste("Y:\\proportionSecondary_NicheVolume.png", sep = ""), plot = myplot,
       width = 200, height = 140, units = 'mm')

rm(myplot)


#############################################################################################
###   Proportion of secondary open habitat - availability of secondary open habitat   ###
#############################################################################################

### Local(10km) availability

myplot <- plotAnalysis(dat, genus_name,
                       xv = "X10km_average_abailability", yv = "proportionSecondaryHabitat", showStats = F,
                       xlabname = "Availability of secondary open habitat within 10 km neighbourhood", ylabname = "Proportion of secondary open habitat",
                       nodeNumbercol = "tag", label.point = T
                       )

# save
ggsave(paste("Y:\\ProportionSecondary_10kmAvailability.png", sep = ""), plot = myplot,
       width = 200, height = 140, units = 'mm')

rm(myplot)

