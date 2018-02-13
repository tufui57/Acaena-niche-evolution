###################################################
### Clade niche
###################################################

source("Y:\\R scripts\\3 Acaena niche evolution\\06_Clade pairing.R")

# Data import
da1 <- read.csv("Y:\\acaena_bioclim_landcover_history_inclNAonland.csv")
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

extent_x = c(min(scores$PC1), max(scores$PC1))
extent_y = c(min(scores$PC2), max(scores$PC2))



########################################################################################
### Calculate node niche
########################################################################################

library(grid)
library(gridExtra)
library(ggplot2)
# Will nest this package into another, phytools?
library(getDescendants)


# Import phylogeny tree data
acaena <- read.nexus("Y:\\Niche change of lineages\\Niche evolution of open habitat species in islands\\Phylogenetic data and trees\\From Angela\\NZ_Acaena_BEAST_output_6gene.tree")


## first get the node numbers of the tips
nodes <- data.frame(sapply(acaena$tip.label, function(x,y) which(y==x), y=acaena$tip.label))
colnames(nodes)<-"nodelabel"

# Modify species names in phylogentic distance file
a <- sapply(strsplit(rownames(nodes), "_Ac"), "[[", 1)
a2 <- gsub("_AY634821", "", gsub("_EU352216", "", a))
a2 <- gsub("novae-", "novae.", a2)
rownames(nodes) <- a2

### Import species name codes
codes <- read.csv("Y:\\traits.csv")
codes$X <- gsub("novae_", "novae.",codes$X)
codes2 <- codes[codes$X %in% rownames(nodes), ]


########################################################################################
### Plot clade niche of internal nodes
########################################################################################

generate_clades <- function(i, pl = F){
  ### Print node status
  print(paste("Target node is", i, rownames(nodes)[i]))
  print(paste("Sister node is", distance2[i,"sisterNode"]))
  print("Sister node contains")
  cat(paste(rownames(nodes)[res[[i]]], collapse  = "\n"))
  
  
  ## Find columns of species in the target node
  
  # If the target node has no descendant species, i.e. the node is a terminal tip
  if(sum(getDescendants(acaena, node = i)) == 0){
    # If the target node has no occurrence records
    if(sum(colnames(scores) %in% rownames(nodes)[i]) == 0){
      
      stop("The target node has no occurrence records")
      
      }else{
    cladesp <- (colnames(scores) == rownames(nodes)[i])
    scores$targetClade <- scores[,rownames(nodes)[i]]
    
      }
  }else{
    
    cladesp <- colnames(scores) %in% rownames(nodes)[getDescendants(acaena, node = i)]
    
    if(sum(cladesp) > 1){
    
    # if the target node has multiple descendant species
    # Create a column showing clade occurrence records
    scores$targetClade <- ifelse(rowSums(scores[,cladesp]) >= 1, 1, 0)
    }else{
      
      # if the target node has just one descendant
      scores$targetClade <- scores[,cladesp]
      }
    

  }
  
  cladeScore <- scores[,c("PC1", "PC2", "targetClade")]
  
  
  ## Find columns of species in the sister node
  sisCladeSp <- colnames(scores) %in% rownames(nodes)[res[[i]]]
  ## Create a column showing clade occurrence records
  if(sum(sisCladeSp) > 1){
    scores$sisClade <- ifelse(rowSums(scores[,sisCladeSp]) >= 1, 1, 0)
  }else{
    if(sum(sisCladeSp) == 0){
      cat(paste("Sister species", rownames(nodes)[res[[i]]], "has no occurrence records."))
      ### stop() doesn't skip the failed for loop, so you have to re-run the loop manually.
      stop("Script stops!")
      
      }else{
      scores$sisClade <- scores[,sisCladeSp]
    }
  }
  sisCladeScore <- scores[,c("PC1", "PC2", "sisClade")]
  
  ## Data
  nodeNumber = i
  # Species name codes
  nodeName = as.character(codes2[codes2$X %in% colnames(scores)[cladesp], "tag"])
  sisnodeNumber = distance2[i,"sisterNode"]
  sisnodeName = as.character(codes2[codes2$X %in% rownames(nodes)[res[[i]]], "tag"])
  cladedata1 <- cladeScore[cladeScore$targetClade == 1,] # Clade 1 PCA data
  cladedata2 <- sisCladeScore[sisCladeScore$sisClade == 1, ]  # Clade 2 PCA data
  
  # Output results
  clades <- list()
  clades[[1]] <- cladedata1
  clades[[2]] <- nodeName
  
  clades[[3]] <- cladedata2
  clades[[4]] <- sisnodeName
  
  clades[[5]] <- paste(nodeNumber, sisnodeNumber)
  
  if(pl == T){
    ## Plot
    SpPair(scores, # Backgroud PCA data
           cladedata1, # Target clade PCA data
           cladedata2, # Sister clade PCA data
           nodeNumber, nodeName, sisnodeName, extent_x, extent_y
    )
  }
  
  return(clades)
  
  
}

SpPair <- function(scores, # Backgroud PCA data
                   cladedata1, # Target clade PCA data
                   cladedata2, # Sister clade PCA data
                   nodeNumber, nodeName, sisnodeName,
                   extent_x, extent_y
                   ) {
  
  pMain <- ggplot() +
    # plot all NZ data points
    geom_point(data = scores, aes(PC1, PC2), color = 'gray90', alpha = 0.25) +
    # point of each sp
    geom_point(data = cladedata1, aes(PC1, PC2), color = "green", alpha = 0.3) +
    # point of each sp
    geom_point(data = cladedata2, aes(PC1, PC2), color = "purple", alpha = 0.3) +
    # extent
    xlim(extent_x) +
    ylim(extent_y) +
    guides(colour = guide_legend(override.aes = list(size = 5, shape=16, alpha=0.7))) +
    # legend position inside plot
    theme(legend.justification = c(1, 1), legend.position = c(1, 1),
          panel.background = element_rect(fill = 'gray96')
    )
  
  pRightTop <- ggplot() +
    # plot all NZ data points
    geom_point(data = scores, aes(PC1, PC2), color = 'gray90', alpha = 0.25) +
    # point of each sp
    geom_point(data = cladedata1, aes(PC1, PC2), color = "green", alpha = 0.3) +
    theme(panel.background = element_rect(fill = 'gray96'))
  
  ## Create legend to show a list of species in the clade
  # Species list of clade 1
  pclade1 <- ggplot(scores, aes(PC1, PC2)) +
    geom_blank() +
    theme(axis.text = element_blank(),
          axis.title = element_blank(),
          line = element_blank(),
          panel.background = element_blank()) +
    # Species list
    annotation_custom(grob = textGrob(paste(nodeName, collapse  = "\n")))
  
  pRightBottom <- ggplot() +
    # plot all NZ data points
    geom_point(data = scores, aes(PC1, PC2), color = 'gray90', alpha = 0.25) +
    # point of each sp
    geom_point(data = cladedata2, aes(PC1, PC2), color = "purple", alpha = 0.3) +
    theme(panel.background = element_rect(fill = 'gray96'))
  
  ## Create legend to show a list of species in the clade
  # Species list of clade 2
  pclade2 <- ggplot(scores, aes(PC1, PC2)) +
    geom_blank() +
    theme(axis.text = element_blank(),
          axis.title = element_blank(),
          line = element_blank(),
          panel.background = element_blank()) +
    # Species list
    annotation_custom(grob = textGrob(paste(sisnodeName, collapse  = "\n")))
  
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
  
}



### Run
# The following nodes must be eliminated. See "res" object.

number <- (1:max(acaena$edge))[-c(8,10,13,14,17,18,21:29,34,35,39,40,42,45)]

result <- lapply(number, generate_clades, pl = T)

