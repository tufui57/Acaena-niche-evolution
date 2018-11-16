##########################################################################
### Stacked bar plot to show that Acaena occurrs in shlublabds
##########################################################################

library(reshape2)
library(ggplot2)

##########################################################################
### Data preparation
##########################################################################

# Import data of current land cover
d <- read.csv(paste("Y://1st chpater_Acaena project//meta data//", genus_name, "_bioclim_landcover_history_worldclim1_1km18sep.csv", sep = ""))

# Current land cover
# Reference of LCDB land cover; LCDB v4.1 - Land Cover Database version 4.1, Mainland New Zealand
d$currentLCDBlandcover <- sapply(d$layer.2, # d$current_landcover1km, 
                               function(x){
                                 # native forest
                                 ifelse(x == 54, "Broadleaved Indigenous Hardwoods",
                                        ifelse(x == 69, "Indigenous Forest",
                                               
                                 ifelse(x %in% c(64,68,71), "exoticForest",
                                        
                             #"nonForest" incl. shrubland and gravely habitat         
                                 ifelse(x == 15, "Alpine Grass/Herbfield", 
                                        ifelse(x == 40, "High Producing Exotic Grassland",
                                               ifelse(x == 41, "Low Producing Grassland",
                                                      ifelse(x == 43, "Tall Tussock Grassland",
                                                             ifelse(x == 44, "Depleted Grassland",
                                                                    ifelse(x == 10, "Sand or Gravel",
                                                                           ifelse(x == 12, "Landslide",
                                                                                  ifelse(x ==  16, "Gravel or Rock",
                                                                                         ifelse(x == 47, "Flaxland",
                                                                                                ifelse(x == 50, "Fernland",
                                                                                                       ifelse(x == 51, "Gorse and/or Broom",
                                                                                                              ifelse(x == 52, "Manuka and/or Kanuka",
                                                                                                                     ifelse(x == 55, "Sub Alpine Shrubland",
                                                                                                                            ifelse(x == 56, "Mixed Exotic Shrubland",
                                                                                                                                   ifelse(x == 58, "Matagouri or Grey Scrub",
                                 ifelse(x %in% c(0,1,2,5,6,14,20,21,22,30,33,45,46), "nonPotentialAcaenaHabitat", NA
                                                      )))))))))))))))))))
                               }
)

## If you need non-potential habitat from data frame.
d$calc <- ifelse(d$currentLCDBlandcover == "exoticForest" | d$currentLCDBlandcover == "nonPotentialAcaenaHabitat", NA, 1)
d2 <- d[!is.na(d$calc), ]

# # otherwise. just 
# d2<-d

# Replace NA with 0
d3 <- data.frame(sapply(d2, function(x){
  ifelse(is.na(x), 0, x)
})
)

# Count the number of occurrences by land cover type
landcoverTable <- lapply(colnames(d3)[grepl("^Acaena", colnames(d3))], function(x){
  dat <- d3[d3[, x] == 1, ]
  numberOcc <- nrow(dat)
  b <- table(dat$currentLCDBlandcover)
  c <- b / numberOcc
  return(cbind(b, c))
}
)

# Calculate ratio of each land cover type
landcoverRatioTable <- do.call(cbind, lapply(landcoverTable, function(d){d[,2]}))
colnames(landcoverRatioTable) <- colnames(d3)[grepl("^Acaena", colnames(d3))]

ratio <- data.frame(t(landcoverRatioTable))

#############################################################################################
### Stacked bar plot showing how much ratio each habitat type in forest and non-forest have.
#############################################################################################
ratio$spname <- gsub("novae.", "novae_",rownames(ratio))

# Load data of preference for open habiatat and proportion of secondary open
propSecondayOpen <- read.csv("Y://1st chpater_Acaena project//meta data//Acaena_elevation.csv")

# # Proportion of secondary open habitat = occ in secondary habitat / occ in primary and secondary habitat
# propSecondayOpen$proportionSecondaryHabitat <- (d$NF.nonF / (d$NF.nonF + d$nonF.nonF))
# 
# # total = no. of all occurrence cells
# propSecondayOpen$total <- d$Occurrence1kmCells
# propSecondayOpen$log10.total <- log10(d$total)
# 
# # Initial range size of all sp = occurrences in old habitat
# propSecondayOpen$initialRangeSize <- d$nonF.nonF
# propSecondayOpen$log10.initialRange <- log10(d$initialRangeSize)
# 
# Preference for open habitat = sum of occ in open (NF.nonF + nonF.nonF) / (forest + open)
propSecondayOpen$PreferenceOpen <- 
  sum(d$currentLandcover == "openHabitat") / sum(d$currentLandcover == "openHabitat"|d$currentLandcover == "nativeForest")
point <- propSecondayOpen[, c("spname", "tag", "proportionSecondaryHabitat", "Preference_open_habitat")]

### Convert data for drowning stacked bar plot
orderP <- order(point$Preference_open_habitat, decreasing = T)

# Change order to descending order of proportion of Secondary open Habitat
ratio2 <- merge(ratio, point, by= "spname")
ratio2$spname <- factor(ratio2$spname)
ratio3 <- ratio2[orderP, !(colnames(ratio2) %in% c("spname", "proportionSecondaryHabitat", "Preference_open_habitat"))]

# Species name shown in the figure should be name tag.
pl.m <- melt(ratio3, id.vars = "tag")

# Change order of species name tag in order to plot bars in the descending order of proportion of Secondary open Habitat
orderP <- order(point$Preference_open_habitat, decreasing = T)
pl.m$tag <- factor(pl.m$tag, levels = ratio2[orderP,"tag"])


################################################################################
# Change order of LCDB land cover types
################################################################################
pl.m$variable <- factor(gsub("\\.", "\\ ", pl.m$variable))

landcoverType <- as.character(levels(pl.m$variable))

# Grassland
grass <- grep("Grass|Gravel|Landslide", landcoverType)

# Forest
forest <- grep("woods|Forest", landcoverType)

# Shrubland
shrub <- which(!(1:length(landcoverType) %in% c(forest,grass)))

# Change order of factors 
pl.m$variable <- factor(pl.m$variable, levels = landcoverType[c(forest, grass, shrub)])

# Create colour gradient ranging from color 1 and color 2
goldfunc <- colorRampPalette(c("gold", "gold4"))
bluefunc <- colorRampPalette(c("dodgerblue", "dodgerblue4"))
blackfunc <- colorRampPalette(c("black", "white"))

# colour list for plot
cols <- c("olivedrab2", "olivedrab4", 
       goldfunc(8),
       bluefunc(7)#, blackfunc(3)
       )

plotS <- ggplot(pl.m) +
  geom_bar(aes(x = tag, y = value, fill = variable), stat = "identity") +
  scale_fill_manual(name = "LCDB land cover", values = cols) +
  ylab("Proportion") +
  xlab("Acaena species") +
  geom_point(data = point, aes_string(x = "tag", y = "proportionSecondaryHabitat"), size = 2, colour = "black") +
  # geom_point(data = point, aes_string(x = "tag", y = "Preference_open_habitat"), size = 2, colour = "red") +
  theme(panel.background = element_rect(fill = 'white'))

ggsave("Y://LCDBhabitatTypes_preferenceOrder.png", plotS , width = 300, height = 200, units = "mm")
