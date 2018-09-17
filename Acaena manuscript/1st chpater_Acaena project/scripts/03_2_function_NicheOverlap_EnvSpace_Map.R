#####################################################################################
### 1. Schoenner D for multiple variables
### 2. Draw env space colouerd by land use change (new/old) with histgrams of axes
### 3. Draw env space for pre/post-human coloured by land use (NF, nonF, others)
#####################################################################################

library(ecospat)
library(ggplot2)
library(gridExtra)
library(raster)
library(rgdal)
library(maptools)


#######################
## new/old habitats
#######################

newold <- function(spname) {

    dsp <- scores[scores[, spname] == 1,]

    # old
    if (sum(dsp[, "landCoverChange"] == "nonF-nonF") != 0) {
        dsp_o <- dsp[dsp[, "landCoverChange"] == "nonF-nonF",]
    } else {
        dsp_o <- NA
    }

    # new
    if (sum(dsp[, "landCoverChange"] == "NF-nonF") != 0) {
        dsp_n <- dsp[dsp[, "landCoverChange"] == "NF-nonF",]
    } else {
        dsp_n <- NA
    }

    res <- list()
    res[[1]] <- dsp
    res[[2]] <- dsp_o # old
    res[[3]] <- dsp_n # new

    return(res)
}


####################################################################
## calc Schoener D by ecospat.niche.overlap()
####################################################################

################
## Schoener D
################
SchoenerD <- function(scores, ne, ol, sp) {
    scores.clim <- scores[, c("PC1", "PC2")]
    scores.sp1 <- ne[, c("PC1", "PC2")]
    scores.sp2 <- ol[, c("PC1", "PC2")]

    # calculation of occurence density and test of niche equivalency and similarity
    z1 <- ecospat.grid.clim.dyn(scores.clim, scores.clim, scores.sp1, R = 100)
    z2 <- ecospat.grid.clim.dyn(scores.clim, scores.clim, scores.sp2, R = 100)

    res <- list()
    ## Schoener D
    res[[1]] <- unlist(ecospat.niche.overlap(z1, z2, cor = T))
    res[[2]] <- unlist(ecospat.niche.overlap(z1, z2, cor = F))
    names(res) <- c(paste(sp, "corrected"), paste(sp, "not corrected"))

    return(res)
}




####################################################################
###       Draw new/old env space with histgrams of axes
####################################################################

## ggplot in multipanel

spPlot <- function(al, spname, D, extent_x, extent_y, save = TRUE) {

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
    name = paste("Schoener's D =", D),
    # Name of each legend factor. 
    # This must be same factors as factors in "colname" of ggplot(aes(colour = colname)), otherwise no legend will be drawn.
    breaks = c("NF-nonF", "nonF-nonF"),
    label = c("Secondary open area", "Primary open area"),
    # colours
    values = c("red", "blue")
    ) +
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
    ggtitle(paste("Land use change of", spname)) +
    theme(panel.background = element_rect(fill = 'gray96'))


    pRight <- ggplot(al, aes(x = PC2)) +
    geom_histogram(data = subset(al, landCoverChange == 'NF-nonF'), fill = "red", alpha = 0.35) +
    geom_histogram(data = subset(al, landCoverChange == 'nonF-nonF'), fill = "blue", alpha = 0.35) +
    xlim(extent_y) +
    xlab(expression(wet %<->% dry)) +
    theme(axis.text.y = element_blank(),
        axis.title.x = element_text(""),
        axis.ticks.y = element_blank()
        ) +
    coord_flip() +
    theme(panel.background = element_rect(fill = 'gray96'))

    pEmpty <- ggplot(scores, aes(PC1, PC2)) +
    geom_blank() +
    ggtitle(paste("N =", nrow(al))) +
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


################################################
## Env space change for each sp
################################################

 envPlot_sp <- function(scores, # data
                        i, #species name
                        sch, #schoener D
                        save # if you want to save the plot or just show it in plot window
                        ) {
    
    ### data preparation
    dsp <- scores[scores[, i] == 1,]

    # old
    if (sum(dsp[, "landCoverChange"] == "nonF-nonF") != 0) {
        ol <- dsp[dsp[, "landCoverChange"] == "nonF-nonF", c("PC1", "PC2", "landCoverChange")]
    } else {
        ol <- NA
    }

    # new
    if (sum(dsp[, "landCoverChange"] == "NF-nonF") != 0) {
        ne <- dsp[dsp[, "landCoverChange"] == "NF-nonF", c("PC1", "PC2", "landCoverChange")]
    } else {
        ne <- NA
    }
    # new & old
    if (sum(scores[scores[, i] == 1, "landCoverChange"] == "nonF-nonF") == 0) {
        al <- ne
    } else {
        if (sum(scores[scores[, i] == 1, "landCoverChange"] == "NF-nonF") == 0) {
            al <- ol
        }
        else {
            # for histgrams
            al <- rbind(ol, ne)
        }
    }

  spPlot(
    al = al,
    extent_x = c(min(scores$PC1), max(scores$PC1)),
    extent_y = c(min(scores$PC2), max(scores$PC2)),
    D = round(sch, digit = 2),
    spname = i, save = save
  )

 }
 
 
 ################################################
 # Env space for pre/post-human
 ################################################
 
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
       # colours
       values = c("green", "brown")
     ) +
     guides(colour = guide_legend(override.aes = list(size = 5, shape=16, alpha=0.7))) +
     
     # legend position inside plot
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
   # change font size
   theme_set(theme_gray(base_size = 18))
   # Plot in multiple panels
   grid.arrange(pTop, pEmpty, pMain, pRight,
                ncol = 2, nrow = 2, widths = c(3, 1), heights = c(1, 3))
   dev.off()
   
 }
 
 