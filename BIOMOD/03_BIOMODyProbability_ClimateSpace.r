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

setwd("Y:\\")

##########################################################################################
######   PROBABILITY BOXPLOTS CATEGOLIZED BY LANDSCAPE CHANGE     ########################
##########################################################################################

dat <- read.csv(".\\Acaena project\\each sp prob csv\\allNZ_bioclim_sp_landuseChange_BIOMODprediction.csv")

### Add species name tag
acaenaName <- colnames(dat) %>% grepl("Acaena", .) %>% colnames(dat)[.]
acaenaName[!(grepl("_probMedian$", acaenaName))] %>% 
  makeTag(., "Acaena")


# x <- colnames(dat)[grep("_probMedian$", colnames(dat))] %>% 
#   data.frame(
#     cbind(
#       as.character(acaenaName$tag), gsub("_probMedian", "", .)
#       )
#     )
colnames(x) <- c("tag", "sp")

###################################
### 1 panel of boxplots per sp
###################################

spPlots <- function(data, number, spnumbers){
  
  # prepare .png
  png(filename = paste("speciesBoxplots", number, ".png", sep = ""), width = 1880, height = 970, units = "px")
  
  # Setup the panels
  layout(t(1:4))
  par(oma = c(5, 4, 1, 0), mar = rep(0.5, 4), cex = 1.5)
  # `mar` controls the space around each boxplot group
  
  lapply(colnames(data)[grep("_probMedian$", colnames(data))][spnumbers],
         function(i) {
           form <- as.formula(paste(i, "~ change"))
           box <- boxplot(form, data = data, ylim = c(0, 1000), axes = F, col = c("salmon", "tan", "darkolivegreen","salmon", "tan", "salmon"))
            
           axis(1, at = 1:6, labels = levels(data$change), las = 2)
           mtext(x[x$sp==gsub("_probMedian", "", i),"tag"], 3, 0, cex = 1.5)
           
           text(c(1:6), 1000, paste("n = ", table(data[data[, gsub("_probMedian", "", i)] == 1, "change"]), sep = ""))
           
           if (i == colnames(data)[grep("_probMedian$", colnames(data))][spnumbers[1]]) {
             axis(2, las = 1)
             mtext("Median Occurrence Probability", 2, 3, cex = 1.5)
           }
         })
  
  dev.off()
}

spPlots(dat, 1, 1:4)
spPlots(dat, 2, 5:8)
spPlots(dat, 3, 9:12)
spPlots(dat, 4, 13:16)
spPlots(dat, 5, 17:18)

####################################################
### 1 panel of boxplots per landscape change type
####################################################

b$sect <- paste(b$barbSection, b$lifeForm)

# change levels of landscpae change
lev <- with(dat, levels(change))
lev[lev == "f-ef"] <- "NF-EF"
lev[lev == "f-n"] <- "NF-nonF"
lev[lev == "f-nf"] <- "NF-NF"
lev[lev == "n-ef"] <- "nonF-EF"
lev[lev == "n-n"] <- "nonF-nonF"
lev[lev == "n-nf"] <- "nonF-NF"
dat <- within(dat, levels(change) <- lev)

### Plot
lapply(levels(dat$change)[-5],
       function(i) {
         # Extract data for the landscape cahnge type
         ldat <- dat[dat$change == i,]
         
         # prepare .png
         png(filename = paste("", i, "landscapeBoxplots.png", sep = ""), width = 1880, height = 970, units = "px")
         
         # Setup the panels
         layout(t(1:18))
         par(oma = c(2, 4, 4, 0), mar = rep(0, 4), cex = 1.5)
         # `mar` controls the space around each boxplot group
         
         # Plot all the boxes
         for (j in colnames(dat)[grep("_probMedian$", colnames(dat))]) {
           
           if (j == "Acaena_novae.zelandiae_probMedian") {
             bcol="lightblue"
           } else {
             bcol <- ifelse(b[b$X == gsub("_probMedian", "", j), "sect"] == "Ancistrum Stoloniferous", "lightblue",
                            ifelse(b[b$X == gsub("_probMedian", "", j), "sect"] == "Microphyllae Rhizomatous", "salmon", "brown"))
           }
           
           box <- boxplot(ldat[, j], ylim = c(0, 1000), axes = FALSE, col = bcol)
           mtext(x[x$sp == gsub("_probMedian", "", j), "tag"], 1, 0, col=bcol, cex = 1.5)
           
           text(1, 1000, #box$stats[nrow(box$stats), 1] + 0.5,
                paste("n = ", sum(ldat[, gsub("_probMedian", "", j)] == 1, na.rm = T), sep = ""), cex = 1.5)
           
           if (j == colnames(dat)[grep("_probMedian$", colnames(dat))][1]) {
             axis(2, las = 1)
             mtext("Median Occurrence Probability", 2, 3, cex = 1.5)
           }
         }
         title(i, outer = TRUE)
         dev.off()
       })


##########################################################################################
####################       PCA gradated by sp probability         ########################
##########################################################################################

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

# PCA plot
spn <- gsub("_probMedian", "", colnames(scores2)[grep("_probMedian$", colnames(scores2))])

# library for adding another PCA figure on the followings.
library(png)

for (i in spn) {

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




