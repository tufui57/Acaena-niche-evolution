########################################################################################
### Rates of niche evolution in Cooney2016
########################################################################################

library(phytools)
library(dplyr)

source(".//Acaena niche evolution//06_Clade pairing.R")
source(".//Acaena niche evolution//BIOMOD//Create_Package_speciseNameCleaning.r")

# Import tree data
acaena <- read.nexus("Y:\\Niche change of lineages\\Niche evolution of open habitat species in islands\\Phylogenetic data and trees\\From Angela\\NZ_Acaena_BEAST_output_6gene.tree")

########################################################################################
### Plot edge length on trees
########################################################################################

acaena$tip.label <- clean_speciesname(acaena$tip.label)

tree <- acaena
trait <- rnorm(length(nodes), mean=0, sd=1)

names(trait) <- tree$tip.label
wn <- fitContinuous(tree, trait, model="white")

ou<-fitContinuous(tree, trait, model="OU") # The fitConinuous function gives an error message signaling that things are screwed up. (But maybe as an empiricist one says "Hey, whatever? things are always giving error messages.?")

bm<-fitContinuous(tree, trait, model="BM")

## Extract likelihoods from all models

lou<-ou$opt$lnL

lwn<-wn$opt$lnL

lbm<-bm$opt$lnL


pchisq(2*(lou-lbm), df=1, lower.tail=F)


