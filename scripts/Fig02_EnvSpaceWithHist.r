###################################################
### Schoener's D for multiple variables
###################################################
## Basic  Schoener's D is calculated for single variable. 
## But here, I want to calculate Shoener's D for multiple variables.


library(ggplot2)
source("Y:\\1st chapter_Acaena project\\Acaena manuscript\\scripts\\03_2_function_NicheOverlap_EnvSpace_Map.R")
source("Y:\\1st chapter_Acaena project\\Acaena manuscript\\scripts\\F05_EnvSpaceWithHist.r")

genus_name = "Acaena"

### Data import

if(genus_name == "Acaena"){
  da1 <- read.csv("Y:\\1st chapter_Acaena project\\Acaena manuscript\\meta data\\Acaena_bioclim_landcover_history_worldclim1_1km18sep.csv")
  
}else{
  da1 <- read.csv("Y:\\2nd chapter_phylogentic niche conservation\\meta data\\Chionochloa_bioclim_landcover_history_worldclim1_1km.csv"
  )
  }

d <- da1[is.na(da1$landCoverChange) == F, ]

# Species names
sname <- colnames(d)[grepl(paste("^", genus_tag, sep=""), colnames(d))]

for(i in sname){
  d[is.na(d[,i]),i] <- 0
}

# get env. corrdinates (PCA axes)
pca <- prcomp(d[, paste("bioclim", c(1, 6, 12, 15), sep = "")],
              center = TRUE,
              scale. = TRUE)
scores <- data.frame(d[, c(colnames(d)[grep("^bioclim", colnames(d))], sname,
                           "x", "y", "preLandcover", "currentLandcover", "landCoverChange")], pca$x[, 1:2])
scores$landCoverChange <- factor(scores$landCoverChange)
scores$pre <- factor(ifelse(scores$preLandcover == 1, "NF", "nonF"))
scores$post <- factor(ifelse(scores$currentLandcover == 1, "NF",
                             ifelse(scores$currentLandcover == 3, "nonF", "non potential habitat" # EF(2) is also non potential habitat
                             ))
)
extent_x = c(min(scores$PC1), max(scores$PC1))
extent_y = c(min(scores$PC2), max(scores$PC2))


#########################################################################################
## Plot environmental space with schoener D value and histgrams of each axis
#########################################################################################

#############################
## Plot by gglot multipanel
#############################

# traits data
t <- read.csv("Y:\\1st chapter_Acaena project\\Acaena manuscript\\raw data\\traits.csv")
# schoener D
sch <- read.csv("Y:\\1st chapter_Acaena project\\Acaena manuscript\\meta data\\Acaena_data_analyses18sep.csv")
sch$spname <- as.character(sch$spname)

###########################################
##   Environmental space of each species
###########################################
lapply(sname, function(i){
  envPlot_sp(scores, i, sch[sch[,"spname"]==i, "corrected.D"], #schoener D
             save=T
             )
  })

## Reset par()
#resetPar <- function() {
#dev.new()
#op <- par(no.readonly = TRUE)
#dev.off()
#op
#}
#par(resetPar())

############################################################################################################
##   Environmental space of Acaena occurrences in priamry and secondary open area
############################################################################################################

dnz <- scores[rowSums(scores[, sname]) > 0,]

# primary open
ol <- dnz[dnz[, "landCoverChange"] == "nonF-nonF", c("PC1", "PC2", "landCoverChange")]
# secondary open
ne <- dnz[dnz[, "landCoverChange"] == "NF-nonF", c("PC1", "PC2", "landCoverChange")]
# current open, primary and secondary open
al <- rbind(ol, ne)

# Calculate Schonner's D between primary and secondary open of Acaena
source("Y:\\1st chapter_Acaena project\\Acaena manuscript\\scripts\\03_2_function_NicheOverlap_EnvSpace_Map.R")
sch <- try(SchoenerD(scores, ne, ol, "Acaena"), silent=T)

# Plots
allspPlot(al, D = round(sch[[1]][1], 2),
          title = "Climate niche of Acaena in primary and secondary open habitats", 
          coln = "landCoverChange",
          extent_x = extent_x, extent_y = extent_y)
SeparatedPrimarySecondary(ne, ol, "all studied Acaena species", extent_x, extent_y, 
                          col1 = "blue", col2 = "red",
                          save = TRUE)

rm(al,ol,ne)

###########################
# Env space change for NZ
###########################

# primary open
ol <- scores[scores[, "landCoverChange"] == "nonF-nonF", c("PC1", "PC2", "landCoverChange")]
# secondary open
ne <- scores[scores[, "landCoverChange"] == "NF-nonF", c("PC1", "PC2", "landCoverChange")]
# primary and secondary open
al <- rbind(ol, ne)

# Calculate Schonner's D between primary and secondary open of Acaena
# source(".\\Acaena niche evolution\\Acaena manuscript\\1st chpater_Acaena project\\scripts\\03_2_function_NicheOverlap_EnvSpace_Map.R")
sch <- try(SchoenerD(scores, ne, ol, "Acaena"), silent=T)

allspPlot(al, D = round(sch[[1]][1], 2),
          title = "New Zealand", 
          coln = "landCoverChange",
          extent_x = extent_x, extent_y = extent_y)
SeparatedPrimarySecondary(ne, ol, "New Zealand", extent_x, extent_y, 
                          col1 = "blue", col2 = "red",
                          save = TRUE)

rm(al,ol,ne)

########################################################################################
# Env space of forest and open habitat for pre-human and current times
########################################################################################

#############################
# Env space for pre-human
#############################

forest <-  scores[scores[, "pre"] == "NF", ]
nonf <-  scores[scores[, "pre"] == "nonF", ]

# Calculate Schonner's D between forest and non-forest
sch <- try(SchoenerD(scores, forest, nonf, "pre-human"), silent=T)

pre_post_envSpace(scores, coln = "pre", D = round(sch[[1]][1], 2), title = "Environmental space of Pre-human")
SeparatedPrimarySecondary(nonf, forest, "Environmental space of Pre-human", extent_x, extent_y,
                         col1 = "green", col2 = "brown",
                         save = TRUE)
rm(forest,nonf)

#############################
# Env space for current
#############################

forest <-  scores[scores[, "post"] == "NF", ]
nonf <-  scores[scores[, "post"] == "nonF", ]

# Calculate Schonner's D between forest and non-forest
sch <- try(SchoenerD(scores, forest, nonf, "current"), silent=T)

pre_post_envSpace(scores, "post", D = round(sch[[1]][1], 2), "Environmental space of current")
SeparatedPrimarySecondary(nonf, forest, "Environmental space of current", extent_x, extent_y,
                          col1 = "green", col2 = "brown",
                          save = TRUE)
rm(forest,nonf)


############################################################################################################
##   Environmental space of forest vs. secondary open area
############################################################################################################

# current forest
fo <- scores[scores[, "landCoverChange"] == "NF-NF", c("PC1", "PC2", "landCoverChange")]
# secondary open
ne <- scores[scores[, "landCoverChange"] == "NF-nonF", c("PC1", "PC2", "landCoverChange")]
# current forest and secondary open
al <- rbind(fo, ne)

forestSecondaryPlot(al, "landCoverChange", "Native forest and Secondary open")
SeparatedPrimarySecondary(ne, fo, "Forest and Secondary open", extent_x, extent_y, 
                          col1 = "green", col2 = "blue",
                          save = TRUE)
