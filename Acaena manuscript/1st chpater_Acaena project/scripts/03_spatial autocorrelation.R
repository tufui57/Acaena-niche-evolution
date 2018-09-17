#############################################
### Get Moran's I for raster data
#############################################

library(spdep)

setwd("Y:\\")
alld <- read.csv("Y://acaena_bioclim_landcover_history_inclNAonland.csv")

pa.colnames <- colnames(alld)[grepl("^Acaena.*", colnames(alld))]
# Replace NA with 0
for(i in  pa.colnames){
  alld[is.na(alld[,i]),i] <- 0
}

# Work out the adjacencies (and number of them) for the grid; d2 needs to be large enough to include cells which share an edge but not a point.
# Thus, if the width of each grid cell is x, the value for d2 needs to be somewhere between (but not including) x and sqrt(2)*x. Using
# d2=1.2x should be fine. d1 just needs to be a number less than x; I've suggested x/10 here.
# Vector of lists of adjacient cells
neighs <- dnearneigh(as.matrix(alld[,c("x","y")]), d1=4900, d2=5200)
# weight list
lw <- nb2listw(neighs, style='B', zero.policy=T)
print(lw, zero.policy=TRUE) # Add the zero.policy=TRUE argument to each subsequent command too.

###################################
### Moran's I for 1km cell data
###################################

resM <- lapply(
  pa.colnames,
  function(x){
    # Moran's I test
    resMt <- moran.test(alld[, x], lw, randomisation=FALSE, zero.policy=T)
    
    print(gsub("_pa", "", x))
    print(resMt)
    return(resMt)
  }
)

resM2 <- sapply(1:length(resM), function(x){resM[[x]]$estimate[1]})
names(resM2) <- pa.colnames

### Moran's I of 2ndary occurrences
alldc <- alld[!is.na(alld$landCoverChange),] 
neighs2nd <- dnearneigh(as.matrix(alldc[alldc$landCoverChange=="NF-nonF",c("x","y")]), d1=4900, d2=5200)
# weight list
lw2nd <- nb2listw(neighs2nd, style='B', zero.policy=T)
print(lw2nd, zero.policy=TRUE) # Add the zero.policy=TRUE argument to each subsequent command too.

resM2nd <- lapply(
  pa.colnames,
  function(x){
    # Moran's I test
    resMt <- moran.test(alldc[alldc$landCoverChange=="NF-nonF", x], lw2nd, randomisation=FALSE, zero.policy=T)
    
    print(x)
    print(resMt)
    return(resMt)
  }
)

### Moran's I of primary occurrences
neighs1 <- dnearneigh(as.matrix(alldc[alldc$landCoverChange=="nonF-nonF",c("x","y")]), d1=4900, d2=5200)
# weight list
lw1 <- nb2listw(neighs1, style='B', zero.policy=T)
print(lw1, zero.policy=TRUE) # Add the zero.policy=TRUE argument to each subsequent command too.

resM1 <- lapply(
  pa.colnames,
  function(x){
    # Moran's I test
    resMt <- moran.test(alldc[alldc$landCoverChange=="nonF-nonF", x], lw1, randomisation=FALSE, zero.policy=T)
    
    print(x)
    print(resMt)
    return(resMt)
  }
)

resM3 <- 
  cbind(resM2, 
        sapply(1:length(resM2nd), function(x){resM2nd[[x]]$estimate[1]}),
        sapply(1:length(resM1), function(x){resM1[[x]]$estimate[1]})
)

colnames(resM3)<-c("all", "second", "primary")
write.csv(resM3, "Y://MoranI1km_131117.csv")


#############################################################
###    Weighted matrix for spatial autocorrelation
#############################################################
wm <- nb2mat(neighs, style='B', zero.policy=T)

n <- length(data[, pa.colnames[1]])
y <- data[, pa.colnames[1]]
ybar <- mean(y)
# Compute (yi-ybar)(yj-ybar) for all pairs.
dy <- y - ybar
g <- expand.grid(dy, dy)
yiyj <- g[,1] * g[,2]

# Make a matrix of the multiplied pairs
pm <- matrix(yiyj, ncol=n)
# And multiply this matrix with the weights to set to zero the value for the pairs that are not adjacent.
pmw <- pm * wm
