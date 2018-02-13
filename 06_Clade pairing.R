########################################################################################
### Phylogenetic distances between nodes
########################################################################################


library(adegenet)
library(geiger)
library(phytools)


### Extract the set of descendant terminal edge (tips) of an internal node.
getDescendants<-function(tree,node,curr=NULL){
  if(is.null(curr)) curr<-vector()
  daughters<-tree$edge[which(tree$edge[,1]==node),2]
  curr<-c(curr,daughters)
  w<-which(daughters>=length(tree$tip))
  if(length(w)>0) for(i in 1:length(w)) 
    curr<-getDescendants(tree,daughters[w[i]],curr)
  return(curr)
}


# Import phylogeny tree data
acaena <- read.nexus("Y:\\Niche change of lineages\\Niche evolution of open habitat species in islands\\Phylogenetic data and trees\\From Angela\\NZ_Acaena_BEAST_output_6gene.tree")

########################################################################################
### Get sister nodes
########################################################################################

### Get sister species (tip) or sister group of a target species
# Extract names of edges (i.e. tips and taxa)
tips <- acaena$tip.label

## first get the node numbers of the tips
nodes <- data.frame(sapply(tips, function(x,y) which(y==x), y=acaena$tip.label))
colnames(nodes)<-"nodelabel"

## Second, find sister node of a target species
res<-list()
for(i in nodes$nodelabel){
  sis <- getSisters(acaena,i)
  
  desnode <- getDescendants(acaena, node = sis) 
  
  if(sum(desnode) == 0){
    res[[i]] <- sis
  }else{
    destip <- desnode[desnode %in% nodes$nodelabel]
    
    res[[i]] <- destip 
  }
  
}

# Name a result list
names(res) <- rownames(nodes)

### Get a sister group of a internal node

# Get node numbers of internal nodes
innode <- (length(acaena$tip.label) + 1):max(acaena$edge)
for(i in innode){
  sis <- getSisters(acaena, i)
  
  if(sum(sis) == 0){
    res[[i]] <- "the farmost ancestor"
  }else{
    desnode <- getDescendants(acaena, node = sis)
    destip <- desnode[desnode %in% nodes$nodelabel]
    res[[i]] <- destip 
  }
}



########################################################################################
### Plot edge length on trees
########################################################################################

### Phylogenetic distance between nodes
# The order of the node pair doesn't matter. dist.nodes(acaena)[i, getSisters(acaena, i)] == dist.nodes(acaena)[getSisters(acaena, i), i]
distance <- sapply(1:max(acaena$edge), function(i){
  c(i, getSisters(acaena, i), dist.nodes(acaena)[i, getSisters(acaena, i)])
})

distance2 <- do.call(rbind, distance)
colnames(distance2) <- c("node", "sisterNode","distance")


########################################################################################
### Calculate node niche
########################################################################################

# source("Y:\\R scripts\\3 Acaena niche evolution\\06_ClimateSpace of sister pairs.R")

# ########################################################################################
# ### Plot edge length on trees
# ########################################################################################
# 
# png("Y:\\Niche change of lineages\\Acaena phylogeny tree.png",width = 2000,height = 750)
# plot(acaena)
# edgelabels(round(acaena$edge.length,5),cex=0.75)
# nodelabels()
# tiplabels()
# dev.off()




