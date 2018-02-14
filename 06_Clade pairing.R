########################################################################################
### Phylogenetic distances between nodes
########################################################################################


# library(adegenet)
# library(geiger)
library(phytools)
library(nichePlot)

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

tipssister <- findSisterNode(acaena)

# Name a result list
names(tipssister) <- rownames(nodes)

### Get a sister group of internal nodes as well as terminal nodes

# Get node numbers of internal nodes
allnodesister <- GetInternalNodeNumber(acaena)

allnodesister[1:length(tipssister)] <- tipssister

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





