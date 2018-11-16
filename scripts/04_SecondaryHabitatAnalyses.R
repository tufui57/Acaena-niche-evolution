########################################################
### Regression figures for ratio of 2ndary open habitat 
########################################################

###################################################################################
### Reproducible data set for analyses (on construction at 3.July.2018) 
###################################################################################

genus_name = "Acaena"
genus_tag = "acaena"

### Generate data, if no data file for analysis existed.
source(".//Acaena niche evolution//Acaena manuscript//1st chpater_Acaena project//scripts//04_2_table_of_landcoverHistory_for_analyses.R")
source(".//Acaena niche evolution//Acaena manuscript//1st chpater_Acaena project//scripts//04_2_calculate_indices_for_analyses.R")

### Import data
d <- read.csv(paste(".//", genus_name, "_data_analyses.csv", sep = ""))

### Library & source
library(dplyr)
source(".//functions//F_speciseNameCleaning_spnameFromPhylogenyTree.r")

### Name tag
d$tag <- makeTag_separate(d$spname, genus_name, "_")
d$tag <- d$tag$tag %>% toupper

### Merge indices
breath <- read.csv(paste(".//", genus_name, "_nicheBreadths.csv", sep = ""))
trait <- read.csv("Y://1st chpater_Acaena project//raw data//traits.csv")

# Name cleaning
name_clean <- function(x){gsub("novae_zelandiae", "novae.zelandiae", x)}
trait$X <- sapply(trait$X, name_clean)

# combination of dispersal type and life form
trait$type <- paste(trait$barbSection, trait$lifeForm, sep="_")
trait$type[trait$type != "Ancistrum_Stoloniferous" & trait$type != "Microphyllae_Rhizomatous"] <- "other"
trait$type <- as.factor(trait$type)

schoenerD <- read.csv("Y://shoenerD18sep_r100.csv")
### Does NA in Schoener's D mean no overlap? ###
schoenerD$corrected.D[is.na(schoenerD$corrected.D)] <- 0

avail <- read.csv("Y://availability_from_forest_primary_10_100km18sep.csv")

nicheVol <- read.csv(".//Acaena_nicheVolume.csv")
colnames(nicheVol)[2] <- "niche_volume"

meanEle <- read.csv(".//elevation18sep.csv")
colnames(meanEle)[2] <- "alt"



dat <- merge(d, breath) %>% 
  merge(., trait, by.x = "spname", by.y = "X") %>% 
  merge(., schoenerD, by.x = "spname", by.y = "X") %>%
  merge(., avail, by.x = "spname", by.y = "X") %>% 
  merge(., meanEle, by.x = "spname", by.y = "X") %>%
  merge(., nicheVol, by.x = "spname", by.y = "X")

write.csv(dat, file = "Y://1st chpater_Acaena project//meta data//Acaena_data_analyses18sep.csv")

#################################################################################
### Use original data set (but unreproducible) for regression figures
#################################################################################

# Previous data derived from wrong land cover data processing
# dat <- read.csv("Y://1st chapter_Acaena project//meta data//Acaena_elevation.csv")

### Library & source
library(ggplot2)
source(".//functions//F_plotAnalysis_clade_niche.R")

dat <- read.csv("Y://1st chapter_Acaena project//Acaena manuscript//meta data//Acaena_data_analyses18sep.csv")

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

myplot <- plotAnalysis(dat, genus_name = "",
                       xv = "log10.total", yv = "proportionSecondaryHabitat", showStats = F,
                       xlabname = "log10(Current range size)", ylabname = "Proportion of occurrence records \n in secondary open habitat",
                       nodeNumbercol = "tag.y", label.point = T,
                       cex=22)

# Save
ggsave(paste("Y:\\log10CurrentRangeSize_proportionSecondary.png", sep = ""), plot = myplot,
       width = 220, height = 150, units = 'mm')

rm(myplot)


#########################################################################
### Proportion of Secondary Open Habitat ~ niche volume
#########################################################################

myplot <- plotAnalysis(data = dat, genus_name = "",
                       xv = "niche_volume", yv = "proportionSecondaryHabitat",  showStats = F,
                       xlabname = "Species niche volume", ylabname = "Proportion of occurrence records \n in secondary open habitat",
                       nodeNumbercol = "tag.y", label.point = T,
                       cex=22)
# Save
ggsave(paste("Y:\\proportionSecondary_NicheVolume.png", sep = ""), plot = myplot,
       width = 230, height = 150, units = 'mm')

rm(myplot)


#############################################################################################
###   Proportion of secondary open habitat - availability of secondary open habitat   ###
#############################################################################################

### Local(10km) availability

myplot <- plotAnalysis(dat, "",
                       xv = "X10km..", yv = "proportionSecondaryHabitat", showStats = F,
                       xlabname = "Availability of secondary open habitat", ylabname = "Proportion of occurrence records \n in secondary open habitat",
                       nodeNumbercol = "tag.y", label.point = T,
                       cex=22)

# save
ggsave(paste("Y:\\ProportionSecondary_10kmAvailability.png", sep = ""), plot = myplot,
       width = 230, height = 150, units = 'mm')

rm(myplot)

