###################################################
### Clade niche overlap/volume
###################################################

source("..\\Acaena niche evolution\\07_ClimateSpace of sister pairs.R")

library(ecospat)
library(adehabitatMA)
library(adehabitatHR)
library(raster)

##########################################
# Calculate Schonner's D
##########################################

scho <- list()

for(i in 1:length(cladedata)){

  scho[[i]] <- SchoenerD_ecospat(scores, "PC1", "PC2", cladedata[[i]][[1]],  cladedata[[i]][[3]])
  
}

dat <-list()
for(i in 1:length(scho)){
  x <- data.frame(scho[[i]])
  colnames(x) <- c("ecospat.corrected", "ecospat.uncorrected")
  dat[[i]] <- x[1,]
}
# Name result of schonner's D
names(dat) <- sapply(cladedata, "[[", 5)
# make colums of node numbers
nodeNo <- strsplit(names(dat), "\ ")
dat$node1 <- as.numeric(do.call(rbind, lapply(nodeNo, "[[",1)))
dat$node2 <- as.numeric(do.call(rbind, lapply(nodeNo, "[[",2)))

write.csv(do.call(rbind, dat), "Y://clade_schoennerD_acaena.csv")

###################################################
### Clade niche volume
###################################################

D <- lapply(lapply(cladedata, "[[", 1), SchoenerD, scores = scores, cladedata2 = scores)

# Collate data
d2 <- do.call(rbind, lapply(D, "[[",1))
d3 <- do.call(rbind, lapply(D, "[[",2))
d4 <- data.frame(cbind(d2,d3))
colnames(d4) <- c("corrected.D","corrected.I","not corrected.D", "not corrected.I")
rownames(d4) <- lapply(cladedata, "[[", 5)

# make colums of node numbers
nodeNo <- strsplit(rownames(d4), "\ ")
d4$node1 <- as.numeric(do.call(rbind, lapply(nodeNo, "[[",1)))


###### Find which nodes I haven't calculated niche volume
no1 <- (1:45)[!(1:45 %in% d4$node1)]

# Nodes with no occurrence data, ancestor nodes for all species, nodes with overlapping species
no2 <- no1[!(no1 %in% c(10,18,13,14,21:28,42,45))]

D2 <- list() 
for(i in no2){
  
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
  
  cladedata1 <- cladeScore[cladeScore$targetClade == 1,]
  
  D <- SchoenerD(scores = scores, cladedata1 = cladedata1, cladedata2 = scores)
  
  D2[[i]]<-D
}


# Collate data
dd2 <- do.call(rbind, lapply(D2, "[[",1))
dd3 <- do.call(rbind, lapply(D2, "[[",2))
dd4 <- data.frame(cbind(dd2,dd3))
colnames(dd4) <- c("corrected.D","corrected.I","not corrected.D", "not corrected.I")
dd4$node1 <- no2

rbind(d4, dd4)

write.csv(rbind(d4, dd4), "Y://clade_nicheVolume_acaena.csv")


########################################################################################
### Count descendant species richness of nodes
########################################################################################

########################################################################################
### Plot niche of sister species pair and its sister nodes
########################################################################################


sispairs <- c(1,20,3,5,16)

plot_sister_ancestor <- function(i){
  
  ### Print node status
  print(paste("Target node is", i, rownames(nodes)[i]))
  print(paste("Sister node is", distance2[i,"sisterNode"]))
  print("Sister node contains")
  cat(paste(rownames(nodes)[res[[i]]], collapse  = "\n"))
  
  
  ## Find parent node of target sister species pair
  ancestor <- acaena$edge[which(i == acaena$edge[,2])]
  ancestorSisNode <- distance2[distance2[,"node"] == ancestor,"sisterNode"]
  
  print(paste("Sister node of parent node of target sister species pair is", ancestorSisNode))
  
  ## Niche of sister node of parent node
  
  # If the target node has no descendant species, i.e. the node is a terminal tip
  if(sum(getDescendants(acaena, node = ancestorSisNode)) == 0){
    
    ancSisClade <- (colnames(scores) == rownames(nodes)[ancestorSisNode])
    scores$ancSisClade <- scores[,rownames(nodes)[ancestorSisNode]]
      
    }else{
    
    ancSisClade <- colnames(scores) %in% rownames(nodes)[getDescendants(acaena, node = ancestorSisNode)]
    
    if(sum(ancSisClade) > 1){
      
      # if the target node has multiple descendant species
      # Create a column showing clade occurrence records
      scores$targetClade <- ifelse(rowSums(scores[,ancSisClade]) >= 1, 1, 0)
    }else{
      
      # if the target node has just one descendant
      scores$ancSisClade <- scores[,ancSisClade]
    }
    
    
  }
  
  ancSisScore <- scores[scores[,"ancSisClade"] == 1, c("PC1", "PC2", "ancSisClade")]
  
  
  ## Target species niche
  cladesp <- colnames(scores) %in% rownames(nodes)[i]
  scores$targetClade <- scores[,cladesp]
  cladeScore <- scores[scores[,"targetClade"] == 1, c("PC1", "PC2", "targetClade")]
  
  
  ## Niche of sister species
  sisCladeSp <- colnames(scores) %in% rownames(nodes)[res[[i]]]
  scores$sisClade <- scores[,sisCladeSp]
  sisCladeScore <- scores[scores[,"sisClade"] == 1, c("PC1", "PC2", "sisClade")]
  
  
  ## Data
  nodeNumber = i
  # Species name codes
  nodeName = as.character(codes2[codes2$X %in% colnames(scores)[cladesp], "tag"])
  sisnodeNumber = distance2[i,"sisterNode"]
  sisnodeName = as.character(codes2[codes2$X %in% rownames(nodes)[res[[i]]], "tag"])
  
  ancestorSisName = ancestorSisNode
  
  ### Plot niche
  
  pMain <- ggplot() +
    # plot all NZ data points
    geom_point(data = scores, aes(PC1, PC2), color = 'gray90', alpha = 0.25) +
    # point of sister node of ancestor
    geom_point(data = ancSisScore, aes(PC1, PC2), color = "blue", alpha = 0.3) +
    
    # point of each sp
    geom_point(data = cladeScore, aes(PC1, PC2), color = "green", alpha = 0.3) +
    # point of each sp
    geom_point(data = sisCladeScore, aes(PC1, PC2), color = "purple", alpha = 0.3) +
    # extent
    xlim(extent_x) +
    ylim(extent_y) +
    # title
    ggtitle(paste(ancestorSisName, nodeName, sisnodeName, sep = " ")) +
    guides(colour = guide_legend(override.aes = list(size = 5, shape=16, alpha=0.7))) +
    # legend position inside plot
    theme(legend.justification = c(1, 1), legend.position = c(1, 1),
          panel.background = element_rect(fill = 'gray96')
    )
  
  ### Plot tow species niche on one figure 
  png(filename = paste("Y:\\niche_", ancestorSisName, nodeName, sisnodeName, ".png"), width = 900, height = 630)
  plot(pMain)
  dev.off()
  
}


lapply(sispairs, plot_sister_ancestor)
