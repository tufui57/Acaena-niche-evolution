###################################################
### Clade niche analysis
###################################################

library(phytools)
library(ape)

##############################################################################
### Data preparation
##############################################################################

###########################
### Internal node ages
###########################

# Import phylogenetic tree data
acaena <- read.nexus("Y:\\Niche change of lineages\\Niche evolution of open habitat species in islands\\Phylogenetic data and trees\\From Angela\\NZ_Acaena_BEAST_output_6gene.tree")

#################################################################################
### Extract the set of terminal edge lengths associated with these tips.
#################################################################################
# Extract names of edges (i.e. tips and taxa)
tips <- acaena$tip.label
## first get the node numbers of the tips
nodes <- sapply(tips, function(x,y) which(y==x), y=acaena$tip.label)

## Distance between all combinations of tips
distances <- dist.nodes(acaena)

## Species age
ages <- sapply(nodes, function(i){
  min(distances[distances[,i]>0,i])
}
)

ages <- read.csv("Y:\\Niche change of lineages\\Niche evolution of open habitat species in islands\\acaena_species_ages.csv")

### Calculate internal node ages
nodeage <- cbind(rep(NA, length(branching.times(acaena))), branching.times(acaena))
colnames(nodeage) <- colnames(ages)
ages2 <- rbind(ages, nodeage)


########################################
###  Clade niche volume - Clade age
########################################

volume <- read.csv("Y://Niche change of lineages//Niche evolution of open habitat species in islands//clade_nicheVolume_acaena.csv")
ages3 <- ages2[rownames(ages2) %in% volume$node1,]

data.age.vol <- cbind(ages3, volume[order(volume$node1), ])


###################################################
### Clade niche overlap - Clade distances
###################################################

source("Y:\\R scripts\\3 Acaena niche evolution\\06_Clade pairing.R")

dis <- data.frame(distance2)
overlap <- read.csv("Y://Niche change of lineages//Niche evolution of open habitat species in islands//clade_schoennerD_acaena.csv")

data <- cbind(overlap, dis[dis$node %in% overlap$node1, ])[ ,c("node1", "node2", 
                                                      "ecospat.corrected", "ecospat.uncorrected", "distance")]

data[order(data$distance),]

########################################################
### Function for plotting
########################################################
library(ggplot2)


plotAnalysis <- function(data, 
                         m, # linear model object
                         xv, yv, # column names of responce and explanatory variable
                         xlabname, ylabname, # axes names for plot
                         showStats = T # TRUE; Show p value and slope of linear model and colour points, FALSE; No stat values and black points
){
  
  if(showStats == T){
    myplot <- ggplot(data, aes_string(x = xv, y = yv, label = "nodeNumber", colour = "sprichness")) +
      geom_point() +
      # text label for points
      geom_text(size=5) +
      # change xy labels
      labs(x = xlabname, y = ylabname) +
      # change text size
      theme(text = element_text(size = 20),
            axis.text.x = element_text(size = 20)) +
      # drow LM line & confident intervals 
      stat_smooth(method = "lm", col = "red") +
      # show stats result as title
      labs(title = paste("Adj R2 =", signif(summary(m)$adj.r.squared, digits = 2),
                         "Intercept =", signif(m$coef[[1]], 2),
                         " Slope =", signif(m$coef[[2]], 2),
                         " P =", signif(summary(m)$coef[2, 4], 2))) +
      theme(panel.background = element_rect(fill = "gray95"))
  } else {
    myplot <- ggplot(data, aes_string(x = xv, y = yv)) +
      geom_point() +
      # change xy labels
      labs(x = xlabname, y = ylabname) +
      # change text size
      theme(text = element_text(size = 20),
            axis.text.x = element_text(size = 20)) +
      # drow LM line & confident intervals 
      stat_smooth(method = "lm", col = "red") +
      theme(panel.background = element_rect(fill = "gray95"), legend.position="none")
  }
  return(myplot)
}

#########################################################################
### Clade Phylogenetic distances ~ niche overlap of occurrence records
#########################################################################

# Eliminate outlier of clade age
m <- lm(distance~ecospat.corrected, data[which(data$node1 != 19),])

myplot <- plotAnalysis(data=data[which(data$node1 != 19),], m=m, xv = "ecospat.corrected", yv = "distance", showStats = T,
                       xlabname = "Niche overlap of occurrence records", ylabname = "Phylogenetic distances between clades")

# save
ggsave(paste("Y:\\clade_pd_nicheoverlap_legend.png", sep = ""), plot = myplot,
       width = 300, height = 210, units = 'mm')

rm(myplot, m)


# Leave outlier of clade age in data
m <- lm(distance~ecospat.corrected, data)

myplot <- plotAnalysis(data=data, m=m, xv = "ecospat.corrected", yv = "distance", showStats = T,
                       xlabname = "Niche overlap of occurrence records", ylabname = "Phylogenetic distances between clades")

# save
ggsave(paste("Y:\\clade_pd_nicheoverlap_outlier.png", sep = ""), plot = myplot,
       width = 300, height = 210, units = 'mm')

rm(myplot, m)

#########################################################################
### Clade ages ~ niche volume of occurrence records
#########################################################################

outlier <- which(max(data.age.vol$x) == data.age.vol$x)

# Eliminate outlier of phylogeny
m <- lm(data.age.vol$x[-outlier] ~ data.age.vol$corrected.D[-outlier])

myplot <- plotAnalysis(data=data.age.vol[-outlier,], m=m, xv = "corrected.D", yv = "x", showStats = T,
                       xlabname = "Niche volume of occurrence records", ylabname = "Clade age")

# save
ggsave(paste("Y:\\clade_age_nicheVolume_legend.png", sep = ""), plot = myplot,
       width = 300, height = 210, units = 'mm')

rm(myplot, m)


# Eliminate outlier of phylogeny
m <- lm(data.age.vol$x ~ data.age.vol$corrected.D)

myplot <- plotAnalysis(data=data.age.vol, m=m, xv = "corrected.D", yv = "x", showStats = T,
                       xlabname = "Niche volume of occurrence records", ylabname = "Clade age")

# save
ggsave(paste("Y:\\clade_age_nicheVolume_outlier.png", sep = ""), plot = myplot,
       width = 300, height = 210, units = 'mm')

rm(myplot, m)
