#######################################################################################
## Functions to draw environmental space with histograms and point maps
#######################################################################################

###############################################################
## Plot niche space with histograms
###############################################################

get_occurrences_in_openHabitat <- function(species, # character string of species name
                                           scores){
  # subset data for PCA axes of species
  scores.s <- scores[scores[, species] == 1, ]
  
  # Find rows of Non open habitat
  nonOpenHabitat <- (scores.s$landCoverChange == "NF-nonPotentialHabitat"| scores.s$landCoverChange == "nonF-nonPotentialHabitat"| scores.s$landCoverChange == "NF-EF" | scores.s$landCoverChange == "NF-NF" | scores.s$landCoverChange == "nonF-EF" | scores.s$landCoverChange == "nonF-NF")
  
  if(sum(nonOpenHabitat) > 0){
    scores.s2 <- scores.s[ - which(nonOpenHabitat),]
  }else{
    scores.s2 <- scores.s
  }
  
  return(scores.s2)
}

###############################################################
## Plot niche space with histograms
###############################################################

# Plot just occurrence records in open habitat and colour by land cover
niche_plot_colourByLandcover <- function(species, # character string of species name
                       scores # data of background
                       ) {
  
  scores.sp <- get_occurrences_in_openHabitat(species, scores)
  
  extent_x = c(min(scores$PC1), max(scores$PC1))
  extent_y = c(min(scores$PC2), max(scores$PC2))
  
  
  # Plot niche space
  pMain <- ggplot() +
    # plot all NZ data points
    geom_point(data = scores, aes(PC1, PC2), color = 'gray90', alpha = 0.25) +
    # point of each sp
    geom_point(data = scores.sp, aes(PC1, PC2, colour = landCoverChange), alpha = 0.1) +
    # extent
    xlim(extent_x) +
    ylim(extent_y) +
    # change point colour and legend title and texts
    scale_colour_manual(
      name="",
      # Name of each legend factor. 
      # This must be same factors as factors in "colname" of ggplot(aes(colour = colname)), otherwise no legend will be drawn.
      breaks = c("NF-nonF", "nonF-nonF"),
      label = c("Secondary open area", "Primary open area"),
      # colours
      values = c("red", "blue")
    ) +
    # legend position inside plot
    theme(axis.title = element_text(size = 15),
          legend.position = "none",
          panel.background = element_rect(fill = 'gray96')
    )
  
  # Histogram of PC1 (x axis)
  pTop <- ggplot(scores.sp, aes(x = PC1)) +
    geom_histogram(data = subset(scores.sp, landCoverChange == 'NF-nonF'), fill = "red", alpha = 0.35) +
    geom_histogram(data = subset(scores.sp, landCoverChange == 'nonF-nonF'), fill = "blue", alpha = 0.35) +
    xlim(extent_x) +
    xlab(expression(hot %<->% cold)) +
    theme(axis.text.x = element_blank(),
          axis.title.x = element_text(size = 15),
          axis.title.y = element_text(""),
          axis.ticks.x = element_blank()
    ) +
    theme(axis.title = element_text(size = 15),
          panel.background = element_rect(fill = 'gray96'))
  
  # Histogram of PC2 (y axis)
  pRight <- ggplot(scores.sp, aes(x = PC2)) +
    geom_histogram(data = subset(scores.sp, landCoverChange == 'NF-nonF'), fill = "red", alpha = 0.35) +
    geom_histogram(data = subset(scores.sp, landCoverChange == 'nonF-nonF'), fill = "blue", alpha = 0.35) +
    xlim(extent_y) +
    xlab(expression(dry %<->% wet)) +
    theme(
      axis.text.y = element_blank(),
      axis.text.x = element_text(angle = 270, vjust = 0.25),
      axis.title.x = element_text(size = 15),
      axis.title.y = element_text(angle = 270, size = 15),
      axis.ticks.y = element_blank()
    ) +
    coord_flip() +
    theme(panel.background = element_rect(fill = 'gray96'))
  
  # Empty panel
  pEmpty <- ggplot(scores, aes(PC1, PC2)) +
    geom_blank() +
    theme(axis.title = element_blank(),
          axis.text = element_blank(),
          line = element_blank(),
          panel.background = element_blank()
    )
  
  # Plot in multiple panels
  combinedP <- grid.arrange(pTop, pEmpty, pMain, pRight,
                            ncol = 2, nrow = 2, widths = c(3, 1), heights = c(1, 3))
  return(combinedP)
}

# Plot all occurrence records regardless of habitat
niche_plot_monoColour <- function(species, # character string of species name
                                         scores # data of background
) {
  
  scores.sp <- scores[scores[, species] == 1, ]
  extent_x = c(min(scores$PC1), max(scores$PC1))
  extent_y = c(min(scores$PC2), max(scores$PC2))
  
  # Plot niche space
  pMain <- ggplot() +
    # plot all NZ data points
    geom_point(data = scores, aes(PC1, PC2), color = 'gray90', alpha = 0.25) +
    # point of each sp
    geom_point(data = scores.sp, aes(PC1, PC2), colour = "blue", alpha = 0.1) +
    # extent
    xlim(extent_x) +
    ylim(extent_y) +
    # legend position inside plot
    theme(axis.title = element_text(size = 15),
          legend.position = "none",
          panel.background = element_rect(fill = 'gray96')
    )
  
  # Histogram of PC1 (x axis)
  pTop <- ggplot(scores.sp, aes(x = PC1)) +
    geom_histogram(data = scores.sp, fill = "blue", alpha = 0.35) +
    xlim(extent_x) +
    xlab(expression(hot %<->% cold)) +
    theme(axis.text.x = element_blank(),
          axis.title.x = element_text(size = 15),
          axis.title.y = element_text(""),
          axis.ticks.x = element_blank()
    ) +
    theme(axis.title = element_text(size = 15),
          panel.background = element_rect(fill = 'gray96'))
  
  # Histogram of PC2 (y axis)
  pRight <- ggplot(scores.sp, aes(x = PC2)) +
    geom_histogram(data = scores.sp, fill = "blue", alpha = 0.35) +
    xlim(extent_y) +
    xlab(expression(dry %<->% wet)) +
    theme(
      axis.text.y = element_blank(),
      axis.text.x = element_text(angle = 270, vjust = 0.25),
      axis.title.x = element_text(size = 15),
      axis.title.y = element_text(angle = 270, size = 15),
      axis.ticks.y = element_blank()
    ) +
    coord_flip() +
    theme(panel.background = element_rect(fill = 'gray96'))
  
  # Empty panel
  pEmpty <- ggplot(scores, aes(PC1, PC2)) +
    geom_blank() +
    theme(axis.title = element_blank(),
          axis.text = element_blank(),
          line = element_blank(),
          panel.background = element_blank()
    )
  
  # Plot in multiple panels
  combinedP <- grid.arrange(pTop, pEmpty, pMain, pRight,
                            ncol = 2, nrow = 2, widths = c(3, 1), heights = c(1, 3))
  return(combinedP)
}

###############################################################
# Plot point map
###############################################################

# Plot just occurrence records in open habitat and colour by land cover
map_plot_colourByLandcover <- function(species, # character string of species name
                     data # data for map
){
  
  # subset data for a species including land use change column
  d.s <- get_occurrences_in_openHabitat(species, data)
  
  # Convert land use change column to numeric
  d.s$changeNo <- NA
  d.s[d.s$landCoverChange == "nonF-nonF", "changeNo"] <- 1
  d.s[d.s$landCoverChange == "NF-nonF", "changeNo"] <- 2
  
  # create point
  pts <- d.s[, c("x", "y")]
  
  # point coordinate system setting
  coordinates(pts) <- d.s[, c("x", "y")]
  proj4pts <- proj4string(ref)
  proj4string(pts) <- CRS(proj4pts)
  # land use change column
  pts$changeNo <- d.s$changeNo
  
  # Plot map
  pMap <- ggplot() +
    geom_polygon(data = nzland, aes(x = long, y = lat, group=group), colour = "gray50", fill='gray90') +
    geom_point(data=d.s, aes(x = x, y = y, color = landCoverChange), alpha = 0.1) +
    
    scale_colour_manual(
      # title
      name = paste("N =", nrow(d.s)),
      breaks = c("nonF-nonF", "NF-nonF"),
      label = c("Primary \n open area","Secondary \n open area"),
      # colours
      values = c("red", "blue")
    ) +
    guides(colour = guide_legend(override.aes = list(size = 5, shape=16, alpha=0.7))) +
    
    labs(x="", y="") +
    #ggtitle(gsub("_", " ", spname[i])) +
    # legend position inside plot at topleft
    theme(legend.text = element_text(size=15),
          legend.title = element_text(size=15),
          plot.title = element_text(family = "Times New Roman", size = 20),
          legend.justification = c(1, 1), legend.position = c(0.3, 1),
          panel.background =  element_blank(), #element_rect(fill = 'gray96'),
          #axis.text.y = element_text(angle = 90, hjust = 0.5)
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()
    )
  
  return(pMap)
  
}

# Plot all occurrence records regardless of habitat
map_plot_monoColour <- function(species, # character string of species name
                                data # data for map
){
  
  # subset data for a species
  d.s <- data[data[, species] == 1, ]
  
  # create point
  pts <- d.s[, c("x", "y")]
  
  # point coordinate system setting
  coordinates(pts) <- d.s[, c("x", "y")]
  proj4pts <- proj4string(ref)
  proj4string(pts) <- CRS(proj4pts)

  # Plot map
  pMap <- ggplot() +
    geom_polygon(data = nzland, aes(x = long, y = lat, group = group), colour = "gray50", fill = 'gray90') +
    geom_point(data = d.s, aes(x = x, y = y), color = "blue", alpha = 0.1) +
    ggtitle(paste("N =", nrow(d.s))) +
    guides(colour = guide_legend(override.aes = list(size = 5, shape = 16, alpha = 0.7))) +
    labs(x = "", y = "") +
    # legend position inside plot at topleft
    theme(legend.text = element_text(size=15),
          legend.title = element_text(size=15),
          plot.title = element_text(family = "Times New Roman", size = 20),
          legend.justification = c(1, 1), legend.position = c(0.3, 1),
          panel.background =  element_blank(), #element_rect(fill = 'gray96'),
          #axis.text.y = element_text(angle = 90, hjust = 0.5)
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()
    )
  
  return(pMap)
  
}


