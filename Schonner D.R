#########################################################################
### How Shoenner's D works?
#########################################################################

# Create 2 pairs of imaginary species which have same size of species niche space and overlapped area,
# but have very different number of records.

# Scarce species pair
ax <- rep(c(100,150,200),3)
ay <- c(100,100,100,150,150,150,200,200,200)
a <- data.frame(ax,ay)

bx <- rep(c(200,250,300),3)
b <- data.frame(bx,ay)

# Dense species pair
cx <- rep(seq(100,200,10),length(seq(100,200,10)))
cy <- unlist(lapply(seq(100,200,10), function(i) rep(i,length(seq(100,200,10)))))
c <- data.frame(cx,cy)

dx <- rep(seq(200,300,10),length(seq(200,300,10)))
d <- data.frame(dx,cy)

bacx <- rep(seq(90,310,10), length(seq(90,210,10)))
bacy <- unlist(lapply(seq(90,210,10), function(i) rep(i,length(seq(90,310,10)))))
background <- data.frame(bacx,bacy)

plot(a, xlim = c(0,400), ylim = c(0,400))
points(b, col="red")
points(c, col="blue")
points(d, col="green")
points(background, col="yellow")
### Schonner's D
### Original method in ecospat to calculate Schonner's D
library(ecospat)
library(adehabitatMA)
SchoenerD <- function(scores.clim, scores.gen1, scores.gen2) {

  # calculation of occurence density and test of niche equivalency and similarity
  z1 <- ecospat.grid.clim.dyn(scores.clim, scores.clim, scores.gen1, R = 50)
  z2 <- ecospat.grid.clim.dyn(scores.clim, scores.clim, scores.gen2, R = 50)
  plot(z1$z, main="Density raster of gen1 occurrences")
  plot(z1$Z, main="Density raster of background")
  
  plot(z2$z, main="Density raster of gen2 occurrences")
  plot(z2$Z, main="Density raster of background")
  
  res <- list()
  ## Schoener D
  res[[1]] <- unlist(ecospat.niche.overlap(z1, z2, cor = T))
  res[[2]] <- unlist(ecospat.niche.overlap(z1, z2, cor = F))
  
  return(res)
}

ab <- SchoenerD(background, a, b)
cd <- SchoenerD(background, c, d)

z3 <- ecospat.grid.clim.dyn(scores.clim, scores.clim, c, R = 50)
z4 <- ecospat.grid.clim.dyn(scores.clim, scores.clim, d, R = 50)
plot(z1$z, main="Density raster of a occurrences")
plot(z3$z, main="Density raster of c")

plot(z2$z, main="Density raster of b occurrences")
plot(z4$z, main="Density raster of d")

library(raster)
ecospat.niche.overlap<-function (z1, z2, cor) 
{
  l <- list()
  if (cor == FALSE) {
    p1 <- as.matrix(z1$z.uncor)/sum(as.matrix(z1$z.uncor))
    p2 <- as.matrix(z2$z.uncor)/sum(as.matrix(z2$z.uncor))
  }
  if (cor == TRUE) {
    p1 <- as.matrix(z1$z.cor)/sum(as.matrix(z1$z.cor))
    p2 <- as.matrix(z2$z.cor)/sum(as.matrix(z2$z.cor))
  }
  D <- 1 - (0.5 * (sum(abs(p1 - p2))))
  H <- sqrt(sum((sqrt(p1) - sqrt(p2))^2))
  I <- 1 - (H^2)/2
  l$D <- D
  l$I <- I
  return(l)
}
