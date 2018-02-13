#########################################################################################################
### Tests for phylogenetic distances, species ages, niche volume and Schoener's D
#########################################################################################################


###################################################################
### Linear model test for phylogenetic distances and Schoener's D
###################################################################

# ## Import data
# phy <- read.csv("Y:\\Niche change of lineages\\Niche evolution of open habitat species in islands\\phylogenetic_distances.csv")
# d <- read.csv("Y:\\Niche change of lineages\\Niche evolution of open habitat species in islands\\schoennerD_Acaena.csv")
# # Wrong name because of mistake of name conversion in niche overlap calculation
# d$X <- gsub("Acaenallida", "Acaena_pallida", d$X)
# d$X <- gsub("Acaena_microphylla_var.uciglochidiata", "Acaena_microphylla_var._pauciglochidiata", d$X)
# 
# # Split species names of combinations for niche overlap calculation
# comname <- sapply(as.character(d$X), function(x){
#   a <- strsplit(x, "_Acaena")
#   c(a[[1]][1], paste("Acaena", a[[1]][2], sep = ""))
# }
# )
# 
# dmat <- data.frame(cbind(t(comname), d$ecospat.corrected, d$corrected))
# colnames(dmat) <- c("spname1", "spname2", "ecospatD", "PlateauD")
# 
# 
# # Modify species names in phylogentic distance file
# a <- sapply(strsplit(colnames(phy), "_Ac"), "[[", 1)
# a2 <- gsub("_AY634821", "", gsub("_EU352216", "", a))
# 
# colnames(phy) <- a2
# rownames(phy) <- a2[-1]
# 
# # Extract shared species list
# dsp <- unique(as.character(rbind(t(comname[1,]), t(comname[2,]))))
# splist <- a2[a2 %in% dsp]
# 
# phy2 <- phy[rownames(phy)[rownames(phy) %in% splist],colnames(phy)[colnames(phy) %in% splist]]
# dmat2 <- dmat[dmat$spname1 %in% splist,]
# dmat3 <- dmat2[dmat2$spname2 %in% splist,]
# 
# 
# # Put phylogenetic distance and niche overlap into one dataframe
# dmat3$pd <- NA
# for(i in 1:nrow(dmat3)){
#   te<-dmat3$spname2[i]
#   dmat3$pd[i] <- phy2[rownames(phy2) == dmat3$spname1[i], as.character(te)]
#   
# }
# 
# write.csv(dmat3, "Y:\\Niche change of lineages\\ShonnerD_phylogenetic_distances.csv")

########################################################
### Function for plotting
########################################################
library(ggplot2)

dmat <- read.csv("Y:\\Niche change of lineages\\Niche evolution of open habitat species in islands\\ShonnerD_phylogenetic_distances.csv")

plotAnalysis <- function(data, 
                         m, # linear model object
                         xv, yv, # column names of responce and explanatory variable
                         xlabname, ylabname, # axes names for plot
                         showStats = T # TRUE; Show p value and slope of linear model and colour points, FALSE; No stat values and black points
){
  
  if(showStats == T){
    myplot <- ggplot(data, aes_string(x = xv, y = yv, colour = "spname1")) +
      geom_point(aes(colour = spname1)) +
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
      geom_point(aes(colour = spname1)) +
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
### Phylogenetic distances ~ niche overlap of occurrence records
#########################################################################

m <- lm(pd~ecospatD, dmat)

myplot <- plotAnalysis(data=dmat, m=m, xv = "ecospatD", yv = "pd", showStats = F,
                       xlabname = "Schonner's D of occurrence records", ylabname = "Phylogenetic distances")

# save
ggsave(paste("Y:\\pd_occD_nolegend.png", sep = ""), plot = myplot,
       width = 300, height = 210, units = 'mm')

rm(myplot, m)


#########################################################################
### Phylogenetic distances ~ niche overlap of Plateau prediction
#########################################################################

m <- lm(pd~PlateauD, dmat)

myplot <- plotAnalysis(data=dmat, m=m, xv = "PlateauD", yv = "pd", showStats = F,
                       xlabname = "Schonner's D of habitat suitability estimated by Plateau", ylabname = "Phylogenetic distances")

# save
ggsave(paste("Y:\\pd_PlateauD_nolegend.png", sep = ""), plot = myplot,
       width = 300, height = 210, units = 'mm')

rm(myplot, m)


########################################################################################
### Node ages
########################################################################################

# library(ape)
# # Import phylogenetic tree data
# acaena <- read.nexus("Y:\\Niche change of lineages\\Niche evolution of open habitat species in islands\\Phylogenetic data and trees\\From Angela\\NZ_Acaena_BEAST_output_6gene.tree")
# 
# ### Extract the set of terminal edge lengths associated with these tips.
# # Extract names of edges (i.e. tips and taxa)
# tips <- acaena$tip.label
# ## first get the node numbers of the tips
# nodes <- sapply(tips, function(x,y) which(y==x), y=acaena$tip.label)
# 
# ## Distance between all combinations of tips
# distances <- dist.nodes(acaena)
# 
# ## Speices age
# ages <- sapply(nodes, function(i){
#   min(distances[distances[,i]>0,i])
# }
# )
# 
# write.csv(ages, "Y:\\Niche change of lineages\\Niche evolution of open habitat species in islands\\acaena_species_ages.csv")

ages <- read.csv("Y:\\Niche change of lineages\\Niche evolution of open habitat species in islands\\acaena_species_ages.csv")
d <- read.csv("Y:\\Niche change of lineages\\Niche evolution of open habitat species in islands\\Acaena_nicheVolume.csv")

# Clean up species names
a <- sapply(strsplit(as.character(ages$X), "_Ac"), "[[", 1)
a2 <- gsub("_AY634821", "", gsub("_EU352216", "", a))
ages$X <- a2

dat <- cbind(d[as.character(d$X) %in% ages$X, c("X","corrected.D")], 
      ages[ages$X %in% as.character(d$X),"x"]
      )

colnames(dat) <- c("spname1", "occ.corrected.D", "sp.age")

### Import proportion of secondary habitat
second <- read.csv("Y://Acaena project//data_analyses.csv")

dat2 <- merge(second, dat, by.x = "spname", by.y = "spname1")
colnames(dat2)[1] <- "spname1"

#########################################################################
### Node ages ~ niche volume
#########################################################################

m <- lm(occ.corrected.D ~ sp.age, dat2)

myplot <- plotAnalysis(data=dat2, m=m, xv = "occ.corrected.D", yv = "sp.age", showStats = F,
                       xlabname = "Species niche volume", ylabname = "Species age")

# save
ggsave(paste("Y:\\nicheVolume_spAge_noLagend.png", sep = ""), plot = myplot,
       width = 300, height = 210, units = 'mm')

rm(myplot, m)

summary(lm(occ.corrected.D ~ sp.age, dat))



#########################################################################
### Proportion of Secondary Open Habitat ~ niche volume
#########################################################################

m <- lm(proportionSecondaryHabitat ~ occ.corrected.D, dat2)

myplot <- plotAnalysis(data=dat2, m=m, xv = "proportionSecondaryHabitat", yv = "occ.corrected.D", showStats = T,
                       xlabname = "Proportion of secondary open habitat", ylabname = "Species age")

# save
ggsave(paste("Y:\\acaena_proportionSecondary_spNicheVolume.png", sep = ""), plot = myplot,
       width = 300, height = 210, units = 'mm')

rm(myplot, m)


#########################################################################
### Proportion of Secondary Open Habitat ~ species age
#########################################################################

m <- lm(proportionSecondaryHabitat ~ sp.age, dat2)

myplot <- plotAnalysis(data=dat2, m=m, xv = "proportionSecondaryHabitat", yv = "sp.age", showStats = T,
                       xlabname = "Proportion of secondary open habitat", ylabname = "Species age")

# save
ggsave(paste("Y:\\acaena_proportionSecondary_spAge.png", sep = ""), plot = myplot,
       width = 300, height = 210, units = 'mm')

rm(myplot, m)