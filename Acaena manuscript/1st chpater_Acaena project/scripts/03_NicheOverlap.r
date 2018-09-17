###################################################
### Schoener's D for multiple variables
###################################################

source(".\\Acaena niche evolution\\Acaena manuscript\\1st chpater_Acaena project\\scripts\\03_2_function_NicheOverlap_EnvSpace_Map.R")

genus_name = "Acaena"

# Data import
# da1 <- read.csv("..\\meta data\\Acaena_bioclim_landcover_history_inclNAonland.csv")

da1 <- read.csv(paste("Y://", genus_name, "_bioclim_landcover_history_worldclim1_1km17sep.csv", sep=""
                      )
)

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
scores$landCoverChange <- factor(scores$landCoverChange)
scores$pre <- factor(ifelse(scores$preLandcover == 1, "NF", "nonF"))
scores$post <- factor(ifelse(scores$currentLandcover == 1, "NF",
        ifelse(scores$currentLandcover == 3, "nonF", "non potential habitat" # EF(2) is also non potential habitat
               ))
        )
extent_x = c(min(scores$PC1), max(scores$PC1))
extent_y = c(min(scores$PC2), max(scores$PC2))

### Import PCA scores
save(scores, file = paste(".\\Scores_", genus_name, "_landcover13sep.data", sep = ""))


############################################################################################################
## Calculation of niche overlap between environmental spaces of priamry and 2ndary open habitat 
############################################################################################################
## Schoener's D is calculated as niche overlap by ecospat.niche.overlap().


# The following function can't handle species having no occurrence in new/old habitats
D <- lapply(s, function(x) {

    if ((sum(scores[scores[, x] == 1, "landCoverChange"] == "nonF-nonF") != 0) & (sum(scores[scores[, x] == 1, "landCoverChange"] == "NF-nonF") != 0)) {
        ol <- newold(x)[[2]]
        ne <- newold(x)[[3]]
        
        re <- try(SchoenerD(scores, ne, ol, x), silent=T)
        return(re)
        
    } else {
        res <- list()
        res[[1]] <- NA
        res[[2]] <- NA
        return(res)
    }
})


m <- matrix(NA, length(sname), 4)
colnames(m) <- c("corrected D", "corrected I", "not corrected D", "not corrected I")
rownames(m) <- sname
for (i in 1:length(D)) {
    if (is.character(D[[i]])) {
        m[i,] <- rep(NA, 4)
    } else {
        m[i, 1:2] <- D[[i]][[1]]
        m[i, 3:4] <- D[[i]][[2]]
    }
}

write.csv(m, "Y://shoenerD17sep_r100.csv")
