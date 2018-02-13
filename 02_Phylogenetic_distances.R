library(adegenet)
library(geiger)
library(phytools)

### Calculate phylogenetic distances from phylogenetic trees
# Output of BEAST
# BEAST output has all of tree candidates estimated on the process of MCMC sampling.
#bea <- read.nexus("Y:\\Chiono_matK.trees")

# Output of TreeAnnotator
# TreeAnnotator output is the one chosen among tree candidates estimated on the process of MCMC sampling.
# Species names must not have space. Replace space with _ (underscore).
acaena <- read.nexus("Y:\\Niche change of lineages\\Phylogenetic data and trees\\From Angela\\NZ_Acaena_BEAST_output_6gene.tree")

########################################################################################
### Phylogenetic distances
########################################################################################

### Extract the set of terminal edge lengths associated with these tips.
# Extract names of edges (i.e. tips and taxa)
tips <- acaena$tip.label
## first get the node numbers of the tips
nodes <- sapply(tips, function(x,y) which(y==x), y=acaena$tip.label)

## Distance between all combinations of tips
distances <- dist.nodes(acaena)

dis <- distances[nodes,nodes]

colnames(dis) <- names(sort(nodes))
rownames(dis) <- names(sort(nodes))

write.csv(dis, "Y:\\Niche change of lineages\\phylogenetic_distances.csv")

########################################################################################
### Plot edge length on trees
########################################################################################

png("Y:\\Niche change of lineages\\Acaena phylogeny tree.png",width = 2000,height = 750)
plot(acaena)
edgelabels(round(acaena$edge.length,5),cex=0.75)
nodelabels()
tiplabels()
dev.off()

########################################################################################
### Visualizing phylogenetic distance matrix
########################################################################################

phy <- read.csv("Y:\\Niche change of lineages\\phylogenetic_distances.csv")
dis <- data.matrix(phy[,-1])

dim <- ncol(dis)

# Prepare species names
tag <- read.csv("Y:\\numbers.csv")
nam <- c(as.character(tag$tag), as.character(phy[21:23,1]))
nam2 <- c(nam[2:10],"A.magellanica",nam[11:12],"MIN_a","MIN_m",nam[14:16],"ROR",nam[17:21])

png("Y:\\Niche change of lineages\\Acaena phylogenetic distance matrix.png",width = 1500,height = 1000)
image(1:dim, 1:dim, dis, axes = FALSE, xlab="", ylab="", col = terrain.colors(100))

axis(1, 1:dim, nam2, cex.axis = 1, las=3)
axis(2, 1:dim, nam2, cex.axis = 1, las=1)

text(expand.grid(1:dim, 1:dim), sprintf("%0.4f", dis), cex=0.6)
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


########################################################################################
### Rates of niche evolution in Cooney2016
########################################################################################
tree<-acaena
trait<-rnorm(length(nodes), mean=0, sd=1)

names(trait)<-tree$tip.label
wn <- fitContinuous(tree, trait, model="white")

ou<-fitContinuous(tree, trait, model="OU") # The fitConinuous function gives an error message signaling that things are screwed up. (But maybe as an empiricist one says "Hey, whatever? things are always giving error messages.?")

bm<-fitContinuous(tree, trait, model="BM")

## Extract likelihoods from all models

lou<-ou$opt$lnL

lwn<-wn$opt$lnL

lbm<-bm$opt$lnL


pchisq(2*(lou-lbm), df=1, lower.tail=F)


