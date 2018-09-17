###################################################
### Schoener's D for multiple variables
###################################################
## Basic  Schoener's D is calculated for single variable. 
## But here, I want to calculate Shoener's D for multiple variables.


library(ggplot2)
source(".\\03_2_function_NicheOverlap_EnvSpace_Map.R")

genus_name = "Acaena"

### Data import
# da1 <- read.csv("..\\meta data\\Acaena_bioclim_landcover_history_inclNAonland.csv")

da1 <- read.csv(paste("Y://", genus_name, "_bioclim_landcover_history_worldclim1_1km14sep.csv", sep=""
)
)

d <- da1[is.na(da1$landCoverChange) == F, ]

# sp names
sname <- colnames(d)[grepl("^Acaena", colnames(d))]
# Acaena emittens and A. minor have no occurrence in primary and/or secondary habitats.
s <- sname[ - c(6, 13)]

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
t <- read.csv("..\\raw data\\traits.csv")
# schoener D
sch <- read.csv("..\\meta data\\shoenerD.csv")

###########################################
##   Environmental space of each species
###########################################
lapply(sname, function(i){
  envPlot_sp(scores, i, sch[sch[,"X"]==i, "corrected.D"], #schoener D
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


allspPlot <- function(al, spname, D, extent_x, extent_y, save = TRUE) {
  
  pMain <- ggplot() +
    # plot all NZ data points
    geom_point(data = scores, aes(PC1, PC2), color = 'gray90', alpha = 0.25) +
    # point of each sp
    geom_point(data = al, aes(PC1, PC2, colour = landCoverChange), alpha = 0.1) +
    # extent
    xlim(extent_x) +
    ylim(extent_y) +
    # change point colour and legend title and texts
    scale_colour_manual(
      # title
      name = paste("Schoener's D =", D, "\nN =", nrow(al)),
      # Name of each legend factor. 
      # This must be same factors as factors in "colname" of ggplot(aes(colour = colname)), otherwise no legend will be drawn.
      breaks = c("nonF-nonF", "NF-nonF"),
      # Change name of points
      labels = c("Primary open area", "Secondary open area"),
      # colours
      values = c("red", "blue")
    ) +
    guides(colour = guide_legend(override.aes = list(size = 5, shape=16, alpha=0.7))) +
    # legend position inside plot
    theme(legend.justification = c(1, 1), legend.position = c(1, 1),
          panel.background = element_rect(fill = 'gray96')
    )
  
  pTop <- ggplot(al, aes(x = PC1)) +
    geom_histogram(data = subset(al, landCoverChange == 'NF-nonF'), fill = "red", alpha = 0.35) +
    geom_histogram(data = subset(al, landCoverChange == 'nonF-nonF'), fill = "blue", alpha = 0.35) +
    xlim(extent_x) +
    xlab(expression(hot %<->% cold)) +
    theme(axis.text.x = element_blank(),
          axis.title.y = element_text(""),
          axis.ticks.x = element_blank()
    ) +
    ggtitle(paste("Environmental space of", spname)) +
    theme(panel.background = element_rect(fill = 'gray96'))
  
  
  pRight <- ggplot(al, aes(x = PC2)) +
    geom_histogram(data = subset(al, landCoverChange == 'NF-nonF'), fill = "red", alpha = 0.35) +
    geom_histogram(data = subset(al, landCoverChange == 'nonF-nonF'), fill = "blue", alpha = 0.35) +
    xlim(extent_y) +
    xlab(expression(dry %<->% wet)) +
    theme(axis.text.y = element_blank(),
          axis.title.y = element_text(angle = 270),
          axis.ticks.y = element_blank()
    ) +
    coord_flip() +
    theme(panel.background = element_rect(fill = 'gray96'))
  
  pEmpty <- ggplot(scores, aes(PC1, PC2)) +
    geom_blank() +
    theme(axis.text = element_blank(),
          axis.title = element_blank(),
          line = element_blank(),
          panel.background = element_blank())
  
  if (save == TRUE) {
    png(filename = paste("Y:\\dist_", spname, "_EnvSpace.png"), width = 900, height = 630)
    # change font size
    theme_set(theme_gray(base_size = 18))
    # Plot in multiple panels
    grid.arrange(pTop, pEmpty, pMain, pRight,
                 ncol = 2, nrow = 2, widths = c(3, 1), heights = c(1, 3))
    dev.off()
  } else {
    res <- list(pMain, pTop, pRight, pEmpty)
    return(res)
  }
  
}

# Plots
allspPlot(al, "Climate niche of Acaena in primary and secondary open habitats", 0.22, extent_x, extent_y)

SeparatedPrimarySecondary <- function(ne, ol, spname, extent_x, extent_y, save = TRUE) {
  
  pPrimary <- ggplot() +
    # plot all NZ data points
    geom_point(data = scores, aes(PC1, PC2), color = 'gray90', alpha = 0.25) +
    # point of each sp
    geom_point(data = ol, aes(PC1, PC2), colour = "blue", size=0.5, alpha = 0.1) +
    # extent
    xlim(extent_x) +
    ylim(extent_y) +
    # legend position inside plot
    theme(legend.position = "none",
          panel.background = element_rect(fill = 'gray96')
    )
  
  pSecondary <- ggplot() +
    # plot all NZ data points
    geom_point(data = scores, aes(PC1, PC2), color = 'gray90', alpha = 0.25) +
    # point of each sp
    geom_point(data = ne, aes(PC1, PC2), colour = "red", size=0.5, alpha = 0.1) +
    # extent
    xlim(extent_x) +
    ylim(extent_y) +
    # legend position inside plot
    theme(legend.position = "none",
          panel.background = element_rect(fill = 'gray96')
    )
  
  if (save == TRUE) {
    png(filename = paste("Y:\\", spname, "_separated_EnvSpace.png", sep=""), width = 400, height = 600)
    # change font size
    theme_set(theme_gray(base_size = 18))
    # Plot in multiple panels
    grid.arrange(pPrimary, pSecondary,
                 ncol = 1, nrow = 2)
    dev.off()
  } else {
    res <- list(pPrimary, pSecondary)
    return(res)
  }
}

SeparatedPrimarySecondary(ne, ol, "all studied Acaena species", extent_x, extent_y, save = TRUE)

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

allspPlot(al,"New Zealand", 0.18, extent_x, extent_y, save = TRUE)
SeparatedPrimarySecondary(ne, ol, "New Zealand", extent_x, extent_y, save = TRUE)

rm(al,ol,ne)

########################################################################################
# Env space of forest and open habitat for pre-human and current times
########################################################################################

pre_post_envSpace <- function(scores, coln, title) {
  
  if (length(levels(scores[, coln])) > 2) {
    scores <-  scores[scores[, coln] == "NF" | scores[, coln] == "nonF", ]
  }
  
  pMain <- ggplot(scores, aes_string("PC1", "PC2", colour = coln)) +
    xlim(extent_x) +
    ylim(extent_y) +
    # alpha
    geom_point(alpha = 0.1) +
    # change point colour and legend title and texts
    scale_colour_manual(
      # title
      name = NULL,
      # Name of each legend factor. 
      # This must be same factors as factors in "colname" of ggplot(aes(colour = colname)), otherwise no legend will be drawn.
      breaks = c("NF", "nonF"),
      label = c("Native forest", "Non-Forest"),
      # Colours
      values = c("green", "brown")
    ) +
    guides(colour = guide_legend(override.aes = list(size = 5, shape=16, alpha=0.7))) +
    
    # Legend position inside plot
    theme(legend.justification = c(1, 1), legend.position = c(1, 1),
          panel.background = element_rect(fill = 'gray96')
    )
  
  pTop <- ggplot(scores, aes(x = PC1)) +
    geom_histogram(data = scores[scores[, coln] == 'NF',], fill = "green", alpha = 0.35) +
    geom_histogram(data = scores[scores[, coln] == 'nonF',], fill = "brown", alpha = 0.35) +
    xlim(extent_x) +
    xlab(expression(hot %<->% cold)) +
    theme(axis.text.x = element_blank(),
          axis.title.y = element_text(""),
          axis.ticks.x = element_blank()
    ) +
    ggtitle(paste(title, "forest and non-forest")) +
    theme(panel.background = element_rect(fill = 'gray96'))
  
  
  pRight <- ggplot(scores, aes(x = PC2)) +
    geom_histogram(data = scores[scores[, coln] == 'NF',], fill = "green", alpha = 0.35) +
    geom_histogram(data = scores[scores[, coln] == 'nonF',], fill = "brown", alpha = 0.35) +
    xlim(extent_y) +
    xlab(expression(dry %<->% wet)) +
    theme(axis.text.y = element_blank(),
          axis.text.x = element_text(angle = 270, vjust = 0.25),
          axis.title.y = element_text(angle = 270),
          axis.ticks.y = element_blank()
    ) +
    coord_flip() +
    theme(panel.background = element_rect(fill = 'gray96'))
  
  pEmpty <- ggplot(scores, aes(PC1, PC2)) +
    geom_blank() +
    theme(axis.text = element_blank(),
          axis.title = element_blank(),
          line = element_blank(),
          panel.background = element_blank())
  
  png(filename = paste("Y:\\dist_", title, "_EnvSpace.png", sep = ""), width = 900, height = 630)
  # Change font size
  theme_set(theme_gray(base_size = 18))
  # Plot in multiple panels
  grid.arrange(pTop, pEmpty, pMain, pRight,
               ncol = 2, nrow = 2, widths = c(3, 1), heights = c(1, 3))
  dev.off()
  
}

SeparatedForestNonforest <- function(forest, nonf, spname, extent_x, extent_y, save = TRUE) {
  
  pForest <- ggplot() +
    # Plot all NZ data points
    geom_point(data = scores, aes(PC1, PC2), color = 'gray90', alpha = 0.25) +
    # Point of each sp
    geom_point(data = forest, aes(PC1, PC2), colour = "green", size = 0.5, alpha = 0.1) +
    # Extent
    xlim(extent_x) +
    ylim(extent_y) +
    # No legend
    theme(legend.position = "none",
          panel.background = element_rect(fill = 'gray96')
    )
  
  pNonforest <- ggplot() +
    # plot all NZ data points
    geom_point(data = scores, aes(PC1, PC2), color = 'gray90', alpha = 0.25) +
    # point of each sp
    geom_point(data = nonf, aes(PC1, PC2), colour = "brown", size = 0.5, alpha = 0.1) +
    # extent
    xlim(extent_x) +
    ylim(extent_y) +
    # No legend
    theme(legend.position = "none",
          panel.background = element_rect(fill = 'gray96')
    )
  
  if (save == TRUE) {
    png(filename = paste("Y:\\", spname, "separated_EnvSpace.png", sep=""), width = 400, height = 600)
    # change font size
    theme_set(theme_gray(base_size = 18))
    # Plot in multiple panels
    grid.arrange(pForest, pNonforest,
                 ncol = 1, nrow = 2)
    dev.off()
  } else {
    res <- list(pForest, pNonforest)
    return(res)
  }
}


#############################
# Env space for pre-human
#############################

forest <-  scores[scores[, "pre"] == "NF", ]
nonf <-  scores[scores[, "pre"] == "nonF", ]

pre_post_envSpace(scores, "pre", "Environmental space of Pre-human")
SeparatedForestNonforest(forest, nonf, "Environmental space of Pre-human", extent_x, extent_y, save = TRUE)

rm(forest,nonf)

#############################
# Env space for current
#############################

forest <-  scores[scores[, "post"] == "NF", ]
nonf <-  scores[scores[, "post"] == "nonF", ]

pre_post_envSpace(scores, "post", "Environmental space of current")
SeparatedForestNonforest(forest, nonf, "Environmental space of current", extent_x, extent_y, save = TRUE)

rm(forest,nonf)
