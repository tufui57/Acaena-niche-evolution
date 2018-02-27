########################################################################################
### Phylogeny tree drawing
########################################################################################

library(phytools)
library(dplyr)

source(".//Acaena niche evolution//BIOMOD//Create_Package_speciseNameCleaning.r")
# Import tree data
acaena <- read.nexus("Y:\\Niche change of lineages\\Niche evolution of open habitat species in islands\\Phylogenetic data and trees\\From Angela\\NZ_Acaena_BEAST_output_6gene.tree")

########################################################################################
### Plot edge length on trees
########################################################################################

acaena$tip.label <- clean_speciesname(acaena$tip.label)

extracted.tree <- extract.clade(acaena, 28)

png("Y:\\Acaena phylogeny tree.png", width = 1800, height = 1250)
plotTree(acaena, adj=0, label.offset=0.75, no.margin=T,cex = 1.5)
edgelabels(round(acaena$edge.length,5),
           adj = 2,
           cex = 1.5
           )
nodelabels(cex=1.5, adj = 2)
tiplabels(cex=1.5,adj = 1)
dev.off()


########################################################################################
### Calculate phylogenetic distances from DNA alignments
########################################################################################
# Read fasta file
matk2 <- read.dna(file = "Y:\\Acaena project\\Data from Greg\\Chionochloa_genetic_data\\Chiono_matK.fasta",
                  format = "fasta")
as.character(matk2)[1:5, 1:10]

# Calculate phylogenetic distances
dm.m <- dist.dna(matk2, model="raw", pairwise.deletion=TRUE)
temp <- as.data.frame(as.matrix(dm.m))
temp[is.na(temp)] <- 0
# Plot pairwise phylogenetic distances
table.paint(temp, cleg = 0, clabel.row = 0.5, clabel.col = 0.5)

### Develop phylogenetic tree
dm.m[!rowSums(!is.finite(dm.m)),]
njs(dm.m)


