


### Test data for probability from Acaena anserinifolia
dat <- read.csv(".\\Acaena project\\each sp prob csv\\allNZ_bioclim_sp_landuseChange_BIOMODprediction.csv")

probANS <- scores2[, c("Acaena_anserinifolia_probMedian","PC1", "PC2")]
probDUM <- scores2[, c("Acaena_dumicola_probMedian","PC1", "PC2")]


library(dismo)
### Use dismo::nicheOverlap
nicheOverlap(probANS, probDUM, stat='D', mask=TRUE, checkNegatives=TRUE) 
