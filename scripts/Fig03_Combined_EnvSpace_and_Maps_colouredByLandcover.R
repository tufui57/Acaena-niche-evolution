###############################################################
## Data preparation 
###############################################################

source(".//Acaena niche evolution//Acaena manuscript//1st chpater_Acaena project//scripts//05_Data_preparation_for_drawing_EnvSpace_and_Maps.r")

###############################################################
## Plot map and niche space 
###############################################################


lapply(seq(1,length(spname),3), function(i){
  
  png(paste("Y://niche_map_3sp", spname[i], ".png", sep=""), width = 800, height = 1300)
  
  pMain1 <- niche_plot_colourByLandcover(spname[i], scores)
  pMap1 <- map_plot_colourByLandcover(spname[i], chdata)
  title1 <- textGrob(gsub("_", " ", spname[i]), gp = gpar(fontface = "bold", cex = 1.75))
  spi <- grid.arrange(pMap1, pMain1, top = title1, ncol = 2, widths = c(3,5))
  
  pMain2 <- niche_plot_colourByLandcover(spname[i+1], scores)
  pMap2 <- map_plot_colourByLandcover(spname[i+1], chdata)
  title2 <- textGrob(gsub("_", " ", spname[i+1]), gp=gpar(fontface="bold", cex=1.75))
  spi2 <- grid.arrange(pMap2, pMain2, top = title2, ncol = 2, widths = c(3,5))
  
  pMain3 <- niche_plot_colourByLandcover(spname[i+2], scores)
  pMap3 <- map_plot_colourByLandcover(spname[i+2], chdata)
  title3 <- textGrob(gsub("_", " ", spname[i+2]), gp=gpar(fontface="bold", cex=1.75))
  spi3 <- grid.arrange(pMap3, pMain3, top = title3, ncol = 2, widths = c(3,5))
  
  # Plot in multiple panels
  #grid.arrange(pMap1, pMain1, pMap2, pMain2, pMap3, pMain3, ncol = 2, nrow = 3, widths = c(3,5))
  grid.arrange(spi, spi2, spi3, nrow = 3)
  
  dev.off()
  
})
