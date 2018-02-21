############################################################################################################
############################         BIOMOD RESULT VISUALIZING               ################################
############################################################################################################
##
## Input of this script:
##  Acaena project\\allNZ_bioclim_sp_landuseChange_BIOMODprediction.csv
##
############################################################################################################

library(png)
library(dplyr)
library(gridExtra)
library(ggplot2)
source(".//Chionochloa niche evolution//makeTag.R")


setwd("Y:\\")

##########################################
######   Import data     
##########################################

dat <- read.csv(".\\Acaena project\\each sp prob csv\\allNZ_bioclim_sp_landuseChange_BIOMODprediction.csv")

### Add species name tag
acaenaName <- colnames(dat) %>% grepl("Acaena", .) %>% colnames(dat)[.]
x <- acaenaName[!(grepl("_probMedian$", acaenaName))] %>% 
  makeTag(., "Acaena")

colnames(x) <- c("sp", "tag")

##########################################################################################
####################       PCA gradated by sp probability         ########################
##########################################################################################

# NA omit
dat2 <- dat[dat$bioclim1 != -9999,]

#####################################################################################################
### NOTE; PCA for model prediction need use same BIOCLIM variables as those used to make the models
#####################################################################################################
pca <- prcomp(dat2[, paste("bioclim", c(1, 6, 12, 15), sep = "")],
              center = TRUE,
              scale. = TRUE
              )

# Extract columns PC 1&2 and probability
scores <- data.frame(
  grep("^bioclim", colnames(dat2)) %>% colnames(dat)[.] %>% dat2[, .],
  pca$x[, 1:2]
  )

# prepare pca data for plots
scores2 <- data.frame(
  grep("_probMedian$", colnames(dat2)) %>% colnames(dat2)[.] %>% dat2[, .],
  scores
  )

# PCA plot
spname <- x$sp

for (i in spname) {

    # you can use string type as colums name with aes_string() instead of aes()
    myplot <- ggplot(scores2, aes_string(x = "PC1", y = "PC2", colour = paste(i, "_probMedian", sep = ""))) +
      geom_point(aes_string(colour = paste(i, "_probMedian", sep = "")), size = 0.1) +
      ggtitle(i) +
      theme(legend.title = element_text(size = 5), plot.title = element_text(size = 15),
            # change panel background colour
            panel.background = element_rect(fill = 'gray70', colour = 'white')
            ) +
      scale_colour_gradientn(
        # same colour as maps, default colour for raster is terrain.color.
        colours = rev(terrain.colors(1000)) # gradient colour can be assigned as a vector. c("green","yellow","salmon","gray90")
        ) 

    ima <- readPNG(
      paste("Acaena project\\figures\\PCA\\landuse PCA\\selected vars\\", i, "_landuse_var1_6_12_15.png", sep = ""),
      FALSE
    )
    myplot2 <- myplot + annotation_raster(ima, ymin = -4.6, ymax= -2.5, xmin = 1, xmax = Inf #, interpolate=T
                               )

    ggsave(paste("", i, "_BIOMOD_probPCAwithOccurrence.png", sep = ""), plot = myplot2,
         width = 300, height = 210, units = 'mm')

}




##########################################################################################
############       PCA gradated with 4 colours by sp probability         #################
##########################################################################################

setwd("Y:\\")

# PCA
dat <- read.csv("Acaena project\\allNZ_bioclim_sp_landuseChange_BIOMODprediction.csv")

# NA omit
dat2 <- dat[dat$bioclim1 != -9999,]


pca <- prcomp(dat2[, paste("bioclim", c(1, 6, 12, 15), sep = "")],
              center = TRUE,
              scale. = TRUE)

# data formatting for plotting
scores <- data.frame(dat2[, c(colnames(dat)[grep("^bioclim", colnames(dat2))], "postlandUse", "prelandUse", "change")], pca$x[, 1:2])

scores$prelandUse <- factor(scores$prelandUse)
scores$postlandUse <- factor(scores$postlandUse)
scores$change <- factor(scores$change)

# prepare pca data for plots
scores2 <- data.frame(dat2[, colnames(dat2)[grep("_probMedian$", colnames(dat2))]], scores)

# convert probability values into 4 classes
acaN <- grepl("_probMedian$", colnames(dat2))




