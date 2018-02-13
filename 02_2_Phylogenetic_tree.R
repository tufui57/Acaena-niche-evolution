########################################################################################
### Phylogeny tree drawing
########################################################################################

library(adegenet)
library(geiger)
library(phytools)

# Import tree data
acaena <- read.nexus("Y:\\Niche change of lineages\\Niche evolution of open habitat species in islands\\Phylogenetic data and trees\\From Angela\\NZ_Acaena_BEAST_output_6gene.tree")

########################################################################################
### Plot edge length on trees
########################################################################################

# Modify species names in phylogentic distance file
a <- sapply(strsplit(acaena$tip.label, "_Ac"), "[[", 1)
a2 <- gsub("_AY634821", "", gsub("_EU352216", "", a))
a2 <- gsub("novae-", "novae.", a2)

acaena$tip.label <- a2

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


