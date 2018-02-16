####################################################################################################
## Plots for two groups of points with background and lists of species names or node number on phylogenetic tree
####################################################################################################

################################################################
## Function; Plots for two groups of points
################################################################

plotTwoGroups <- function(background,
                     axis1,axis2,
                     data1, data2,
                     col1, col2, # colours of points
                     extent_x, extent_y
) {
                              
  
  pMain <- ggplot() +
    # plot all NZ data points
    geom_point(data = background, aes_string(axis1, axis2), color = 'gray90', alpha = 0.25) +
    # point of group 1
    geom_point(data = data1, aes_string(axis1, axis2), color = col1, alpha = 0.3) +
    # point of group 2
    geom_point(data = data2, aes_string(axis1, axis2), color = "purple", alpha = 0.3) +
    # extent
    xlim(extent_x) +
    ylim(extent_y) +
    guides(colour = guide_legend(override.aes = list(size = 5, shape=16, alpha=0.7))) +
    # legend position inside plot
    theme(legend.justification = c(1, 1), legend.position = c(1, 1),
          panel.background = element_rect(fill = 'gray96')
    )
  
  return(pMain)
}

################################################################
## Function; Plots for one group of points
################################################################

plotOneGroup <- function(background,
                     axis1, axis2,
                     dat, # Dataframe of a group of points
                     col,
                     extent_x, extent_y
) {
  p <- ggplot() +
    # plot all NZ data points
    geom_point(data = background, aes_string(axis1, axis2), color = 'gray90', alpha = 0.25) +
    # point of each sp
    geom_point(data = dat, aes_string(axis1, axis2), color = col, alpha = 0.3) +
    theme(panel.background = element_rect(fill = 'gray96'))
  
  return(p)
  }



################################################################
## Function;  Create legend to show a list of species in the clade
################################################################

plotSpNameList <- function(background,
                         axis1, axis2,
                         col,
                         nodeName, # list of taxon names in the point group
                         extent_x, extent_y
) {
  p <- ggplot(background, aes_string(axis1, axis2)) +
    geom_blank() +
    theme(axis.text = element_blank(),
          axis.title = element_blank(),
          line = element_blank(),
          panel.background = element_blank()) +
    # Species list
    annotation_custom(grob = textGrob(paste(nodeName, collapse  = "\n")))
  
  return(p)
  
  }
  

################################################################
## Final function; Plot all
################################################################

ploTwoGroupWithSpNames <- function(background,
                   axis1, axis2, # Names of coordinates
                   data1, data2, # Dataframes of two groups of points
                   col1, col2,
                   nodeNumber, # Character of target node and sister node numbers
                   nodeName, sisnodeName,
                   extent_x, extent_y,
                   save = TRUE
) {
  
  pMain <- plotTwoGroups(background, 
                    axis1, axis2,
                    data1, data2,
                    col1, col2,
                    extent_x, extent_y
                    )
  
  pRightTop <- plotOneGroup(background,  
                            axis1, axis2,
                            data1,
                            col = col1,
                            extent_x, extent_y)
  
  ## Create legend to show a list of species in the clade
  # Species list of clade 1
  pclade1 <- plotSpNameList(background,
                          axis1, axis2,
                          col,
                          nodeName,
                          extent_x, extent_y
                          )
    
    
  pRightBottom <- plotOneGroup(background,  
                               axis1, axis2,
                               data2,
                               col = col2,
                               extent_x, extent_y)
  
  ## Create legend to show a list of species in the clade
  # Species list of clade 2
  pclade2 <- plotSpNameList(background,
                          axis1, axis2,
                          col,
                          sisnodeName,
                          extent_x, extent_y
  )

  if(save == TRUE){
    ### Plot tow species niche on one figure 
    png(filename = paste("Y:\\niche_", nodeNumber, ".png"), width = 900, height = 630)
    plot(pMain)
    dev.off()
    
    ### Plot two species niche separately
    png(filename = paste("Y:\\niche_", nodeNumber, "_separate.png"), width = 450, height = 630)
    # change font size
    theme_set(theme_gray(base_size = 18))
    # Plot in multiple panels
    grid.arrange(pRightTop, pclade1, pRightBottom, pclade2, ncol = 2, nrow = 2, widths = c(3,1), heights = c(1,1))
    dev.off()
    
  }else{
    return(pMain)
  }
  
}
