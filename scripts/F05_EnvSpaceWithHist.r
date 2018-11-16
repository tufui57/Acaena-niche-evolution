
### Function to draw forest vs. non-forest climate space
pre_post_envSpace <- function(scores, # scores for climate points
                              coln, # Column names of data to colour climate points
                              D, # Schoener's D
                              title # Plot title
                              ){
  
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
    name = paste("Schoener's D =", D),
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
  
  envSpace_histograms(pMain, # ggplot object
                                scores1 = scores[scores[, coln] == 'NF',], # scores for climate points
                                scores2 = scores[scores[, coln] == 'nonF',], # scores for climate points
                                coln = coln, # Column names of data to colour climate points
                                title = title, # Plot title
                                fill.cols = c("green", "red")# vector of two colour names of points
                                )
}


### Function to draw histograms
envSpace_histograms <- function(pMain, # ggplot object
  scores1, # scores for climate points
                                scores2, # scores for climate points
                              coln, # Column names of data to colour climate points
                              title, # Plot title
                              fill.cols # vector of two colour names of points
) {
  
  pTop <- ggplot(scores, aes(x = PC1)) +
    geom_histogram(data = scores1, fill = fill.cols[[1]], alpha = 0.35) +
    geom_histogram(data = scores2, fill = fill.cols[[2]], alpha = 0.35) +
    xlim(extent_x) +
    xlab(expression(hot %<->% cold)) +
    theme(axis.text.x = element_blank(),
          axis.title.y = element_text(""),
          axis.ticks.x = element_blank()
    ) +
    ggtitle(title) +
    theme(panel.background = element_rect(fill = 'gray96'))
  
  
  pRight <- ggplot(scores, aes(x = PC2)) +
    geom_histogram(data = scores1, fill = fill.cols[[1]], alpha = 0.35) +
    geom_histogram(data = scores2, fill = fill.cols[[2]], alpha = 0.35) +
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
  
  png(filename = paste("Y:\\", title, ".png", sep = ""), width = 900, height = 630)
  # Change font size
  theme_set(theme_gray(base_size = 18))
  # Plot in multiple panels
  grid.arrange(pTop, pEmpty, pMain, pRight,
               ncol = 2, nrow = 2, widths = c(3, 1), heights = c(1, 3))
  dev.off()
  
}

### Function to draw primary vs. secondary open habitats' climate spaces
allspPlot <- function(al, # scores for climate points
                      D, # Schoenner's D
                      coln, # Column names of data to colour climate points
                      title, # file title to save the plot
                      extent_x, extent_y, save = TRUE) {
  
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
    guides(colour = guide_legend(override.aes = list(size = 5, shape = 16, alpha = 0.7))) +
    # legend position inside plot
    theme(legend.justification = c(1, 1), legend.position = c(1, 1),
          panel.background = element_rect(fill = 'gray96')
    )
  envSpace_histograms(pMain, coln,
                      scores1 = al[al[, coln] == 'NF-nonF',],
                      scores2 = al[al[, coln] == 'nonF-nonF',],
                      title = title,
                      fill.cols = c("red", "blue")
  )
  
}


### Function to plot two climate spaces separately
SeparatedPrimarySecondary <- function(ne, ol, 
                                      spname, # file name to save the resultant plot
                                      extent_x, extent_y, 
                                      col1, # character;colour of "ne" points 
                                      col2, # character;colour of "ol" points 
                                      save = TRUE) {
  
  pPrimary <- ggplot() +
    # plot all NZ data points
    geom_point(data = scores, aes(PC1, PC2), color = 'gray90', alpha = 0.25) +
    # point of "ol"
    geom_point(data = ol, aes(PC1, PC2), colour = col1, size=0.5, alpha = 0.1) +
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
    # point of "ne"
    geom_point(data = ne, aes(PC1, PC2), colour = col2, size=0.5, alpha = 0.1) +
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


### Function to draw current forest vs. secondary open habitats' climate spaces
forestSecondaryPlot <- function(scores, 
                                coln, # Column names of data to colour climate points
                                title # file title to save the plot
                                ) {

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
      breaks = c("NF-NF", "NF-nonF"),
      label = c("Native forest", "Secondary open"),
      # Colours
      values = c("green", "blue")
    ) +
    guides(colour = guide_legend(override.aes = list(size = 5, shape=16, alpha=0.7))) +

    # Legend position inside plot
    theme(legend.justification = c(1, 1), legend.position = c(1, 1),
          panel.background = element_rect(fill = 'gray96')
    )

  envSpace_histograms(pMain, coln,
                      scores1 = scores[scores[, coln] == 'NF-NF',],
                      scores2 = scores[scores[, coln] == 'NF-nonF',],
                      title = title,
                      fill.cols = c("green", "blue")
  )

}