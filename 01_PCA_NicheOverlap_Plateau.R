
library(adehabitatMA)
library(adehabitatHR)
library(raster)

# Data import
alld <- read.csv("Y://Plateau model//alldata_Acaena5kmGrid.csv")
d <- alld[is.na(alld$bio1) == F, ]

# sp names
sname <- colnames(d)[grepl("^Acaena.*pa$", colnames(d))]

# get env. corrdinates (PCA axes)
pca <- prcomp(d[, paste("bio", c(1, 6, 12, 15), sep = "")],
              center = TRUE,
              scale. = TRUE)

prob <- read.csv("Y://Plateau model//plateau_prob_current_Acaena_agnipila.csv")
prob$prob.abs <- abs(prob$prob - 1000)
scores <- data.frame(d[, c(paste("bio", c(1, 6, 12, 15), sep = ""), "NZTMlon", "NZTMlat")], pca$x[, 1:2], prob$prob.abs)

glob <- scores[, c("PC1", "PC2")]
glob1<- scores[, c("PC1", "PC2")]
sp.scores<-scores


grid.clim.dyn_from_model.prediction <-
function (glob, glob1, sp.scores, R, th.sp = 0.50, th.env = 0.50){

  # Data preparation
  glob <- as.matrix(glob)
  glob1 <- as.matrix(glob1)
  sp <- as.matrix(sp.scores[, c("PC1", "PC2")])
  l <- list()

  if (ncol(glob) != 2) {
    cat("Arguments \"glob\" and \"glob1\" should have 2 columns of coordinates")
  }
  
  if (ncol(glob) == 2) {
    
    ## Create cell density raster
    xmin <- min(glob[, 1])
    xmax <- max(glob[, 1])
    ymin <- min(glob[, 2])
    ymax <- max(glob[, 2])
    glob1r <- data.frame(cbind((glob1[, 1] - xmin)/abs(xmax - xmin), (glob1[, 2] - ymin)/abs(ymax - ymin)))
    spr <- data.frame(cbind((sp[, 1] - xmin)/abs(xmax - xmin), 
                            (sp[, 2] - ymin)/abs(ymax - ymin)))
    mask <- ascgen(SpatialPoints(cbind((1:R)/R, (1:R)/R)), 
                   nrcol = R - 2, count = FALSE)
    
    glob1.dens <- kernelUD(SpatialPoints(glob1r[, 1:2]), 
                           grid = mask, kern = "bivnorm")
    glob1.dens <- raster(xmn = xmin, xmx = xmax, ymn = ymin, 
                         ymx = ymax, matrix(glob1.dens$ud, nrow = R))
    ## Resample probability raster on the R x R resolution and extent of PC1 and 2
    # Create raster in 5km resolution first.
    ref.raster <- raster("Y:\\GIS map and Climate data\\worldclim\\bio_5km_grid\\bioclim1NZTM5km.bil")
    
    points <- sp.scores[, c("PC1", "PC2")]
    # Set coordinates
    coordinates(points) <- sp.scores[,c("PC1", "PC2")]
    # Put values in point object
    points$prob <- sp.scores$prob.prob.abs
    sp.dens <- rasterize(points, glob1.dens, field = points$prob)
    
    ## Resample raster to (R x R) resolution
    x <- seq(from = min(glob[, 1]), to = max(glob[, 1]), 
             length.out = R)
    y <- seq(from = min(glob[, 2]), to = max(glob[, 2]), 
             length.out = R)
    glob1r <- extract(glob1.dens, glob1)
    # "th.env" % of quantile value of background cell density
    Z.th <- quantile(glob1r, th.env)
    # Calculation of raster values
    glob1.dens[glob1.dens < Z.th] <- 0
    Z <- glob1.dens * nrow(glob1)/cellStats(glob1.dens, "sum")
    
    ## Repeat the same for sp occ density raster
    spr <- extract(sp.dens, sp)
    z.th <- quantile(spr, th.sp)
    sp.dens[Z == 0] <- 0
    sp.dens[sp.dens < z.th] <- 0
    z <- sp.dens * nrow(sp)/cellStats(sp.dens, "sum")
    
    # Calculate corrected/uncorrected Z
    z.uncor <- z/cellStats(z, "max")
    w <- z.uncor
    w[w > 0] <- 1
    z.cor <- z/Z
    z.cor[is.na(z.cor)] <- 0
    z.cor <- z.cor/cellStats(z.cor, "max")
    
    # Put all results in list l
    l$x <- x
    l$y <- y
    l$z <- z
    l$z.uncor <- z.uncor
    l$z.cor <- z.cor
    l$Z <- Z
    l$w <- w
    l$sp.dens <- sp.dens
  }
  return(l)
}

### Original method in ecospat to calculate Schonner's D
library(ecospat)

SchoenerD <- function(scores, gen1, gen2) {
  
  # Extract data of two target species
  scores.gen1 <-  scores[scores[, gen1] == 1,c("PC1", "PC2")]
  scores.gen2 <-  scores[scores[, gen2] == 1,c("PC1", "PC2")]
  scores.clim <- scores[, c("PC1", "PC2")]
  
  # calculation of occurence density and test of niche equivalency and similarity
  z1 <- ecospat.grid.clim.dyn(scores.clim, scores.clim, scores.gen1, R = 100)
  z2 <- ecospat.grid.clim.dyn(scores.clim, scores.clim, scores.gen2, R = 100)
  
  res <- list()
  ## Schoener D
  res[[1]] <- unlist(ecospat.niche.overlap(z1, z2, cor = T))
  res[[2]] <- unlist(ecospat.niche.overlap(z1, z2, cor = F))
  # Name
  name_genera <- paste(strsplit(gen1, ".csv"),strsplit(gen2, ".csv"), sep="_")
  names(res) <- c(paste(name_genera, "corrected"), paste(name_genera, "not corrected"))
  
  return(res)
}

com <- combn(sname[-2],2)

result <- list()

for(i in 1:ncol(com)){
  
  scores <- data.frame(d[, c(paste("bio", c(1, 6, 12, 15), sep = ""),sname, "NZTMlon", "NZTMlat")], pca$x[, 1:2], prob$prob.abs)
  gen1=com[,i][1]
  gen2=com[,i][2]
  gen1.name = gsub("_pa$", "", gen1)
  gen2.name = gsub("_pa$", "", gen2)
  
  ### Niche overlap for occurrence records
  r <- SchoenerD(scores, gen1, gen2)
  
  ### Niche overlap for Plateau model estimates
  # agnipila
  prob <- read.csv(paste("Y://Plateau model//plateau_prob_current_", gen1.name, ".csv", sep=""))
  prob$prob.abs <- abs(prob$prob - 1000)
  scores.gen1 <- data.frame(d[, c(paste("bio", c(1, 6, 12, 15), sep = ""), "NZTMlon", "NZTMlat")], pca$x[, 1:2], prob$prob.abs)
  
  gen1.dens <- grid.clim.dyn_from_model.prediction(scores.gen1[, c("PC1", "PC2")], scores.gen1[, c("PC1", "PC2")], 
                                             scores.gen1, R=100, th.sp = 0.50, th.env = 0.50)
  
  # juvenca
  prob <- read.csv(paste("Y://Plateau model//plateau_prob_current_", gen2.name, ".csv", sep=""))
  prob$prob.abs <- abs(prob$prob - 1000)
  scores.gen2 <- data.frame(d[, c(paste("bio", c(1, 6, 12, 15), sep = ""), "NZTMlon", "NZTMlat")], pca$x[, 1:2], prob$prob.abs)
  
  gen2.dens <- grid.clim.dyn_from_model.prediction(scores.gen2[, c("PC1", "PC2")], scores.gen2[, c("PC1", "PC2")], 
                                             scores.gen2, R=100, th.sp = 0.50, th.env = 0.50)
  res <- list()
  ## Schoener D
  res[[1]] <- unlist(ecospat.niche.overlap(gen1.dens, gen2.dens, cor = T))
  res[[2]] <- unlist(ecospat.niche.overlap(gen1.dens, gen2.dens, cor = F))
  # Name
  name_genera <- paste(strsplit(gen1, ".csv"),strsplit(gen2, ".csv"), sep="_")
  names(res) <- c(paste(name_genera, "corrected"), paste(name_genera, "not corrected"))
  
  result[[i]] <- data.frame(c(r,res))
  
}

dat <-list()
for(i in 1:length(result)){
  x <- data.frame(result[[i]])
  colnames(x) <- c("ecospat.corrected", "ecospat.uncorrected", "corrected", "uncorrected")
  dat[[i]] <- x[1,]
}
names(dat)<-sapply(result, function(x){
  a <- gsub("_pa", "", colnames(x)[1])
  gsub(".corrected","", a)
  }
  )


write.csv(do.call(rbind, dat), "Y://schoennerD_Acaena.csv")

