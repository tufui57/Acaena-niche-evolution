############################################################################################################
##########       Climate space with predicted probability                     ##############################
############################################################################################################
library(RColorBrewer)
library(ggplot2)

# Data import
alld <- read.csv("Y://Plateau model//alldata_Acaena5kmGrid.csv")
d <- alld[is.na(alld$bio1) == F, ]

# sp names
sname <- gsub("_pa$", "", colnames(d)[grepl("^Acaena.*pa$", colnames(d))])

# get env. corrdinates (PCA axes)
pca <- prcomp(d[, paste("bio", c(1, 6, 12, 15), sep = "")],
              center = TRUE,
              scale. = TRUE)

myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))

for(i in sname[-2]){
  
  prob <- read.csv(paste("Y://Plateau model//plateau_prob_current_", i, ".csv", sep = ""))
  prob$prob.abs <- abs(prob$prob - 1000)
  scores <- data.frame(d[, c(paste("bio", c(1, 6, 12, 15), sep = ""), "NZTMlon", "NZTMlat")], pca$x[, 1:2], prob$prob.abs)
  
  ### Plot PCA
  myplot <- ggplot(scores, aes(x=PC1, y=PC2, colour = prob.prob.abs)) +
    geom_point(aes(colour = prob.prob.abs), size=0.3) +
    ggtitle(i) +
    theme(legend.title=element_text(size=5), plot.title = element_text(size=15)) +
    scale_colour_gradientn(name="Environemntal suitability",
                           limits = c(0,1000), # Data value range of scales should be shared among all species
                           colours = rev(terrain.colors(1000)),
                           breaks=c(0,500,1000), labels=format(c(0,500,1000))
                           )
  
  # Save the plot
  ggsave(paste("Y://Plateau model//", i, "_PlateauPCA.png", sep = ""), plot = myplot,
         width=300, height=210, units='mm')
  
}
  