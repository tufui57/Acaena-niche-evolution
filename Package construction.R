
#######################################################
### Package construction
#######################################################

### Install libraries to package construction first.
### Libraries
library(roxygen2)
library(devtools)



# Remove package whose name is same as the one you're going to create to avoid overwrite or duplication.
# Replace "packagename" and "functionname" your object names. 

packagename <- "nichePlot"
system("rm -Rf nichePlot")
remove.packages(packagename)

# Create package.
create(packagename)
# Show files in the package
list.files(path = packagename, recursive=TRUE, include.dirs=TRUE)

## Add your function

#######################################################
### Function 1
#######################################################

functionname <- "SchoenerD_ecospat"
### Calculation of niche overlap 
# Original method in ecospat to calculate Schonner's D

SchoenerD_ecospat <- function(background,
                              axis1, 
                              axis2,
                              data1, 
                              data2,
                              R = 100 # Resolution of background
) {
  
  background.clim <- background[, c(axis1, axis2)]
  
  # calculation of occurence density and test of niche equivalency and similarity
  z1 <- ecospat.grid.clim.dyn(background.clim, background.clim, data1[ ,c(axis1, axis2)], R = 100)
  z2 <- ecospat.grid.clim.dyn(background.clim, background.clim, data2[ ,c(axis1, axis2)], R = 100)
  
  res <- list()
  ## Schoener D
  res[[1]] <- unlist(ecospat.niche.overlap(z1, z2, cor = T))
  res[[2]] <- unlist(ecospat.niche.overlap(z1, z2, cor = F))
  ## Name
  return(res)
}

# Create documentation of the package.
# Then put yout function and the documentation in the package.

x <- c("#' A SchoenerD_ecospat Function",
     "#'",
     "#' This function calculates Schoener's D, index of niche overlap between two groups of points",
     "#' @param background Data frame containing background area",
     "#' @param axis1 Colum name of coordinate 1", 
     "#' @param axis2 Colum name of coordinate 2",
     "#' @param data1 Data frame of group 1", 
     "#' @param data2 Data frame of group 2",
     "#' @param R Resolution of grid cells. Default is 100",
     "#' @keywords niche overlap",
     "#' @export",
     "#' @examples",
     "#' SchoenerD()",
     "")

write(x, file =  
paste(packagename, "/R/", functionname, ".R", sep = "")
)

dump(functionname, file = paste(packagename, '/R/SchoenerD_ecospat.R', sep = ""), append = TRUE)

#######################################################
### Function 2
#######################################################

### Extract the set of descendant terminal edge (tips) of an internal node.

getDescendants <- function(
  tree,
  node,
  curr=NULL){
  
  if(is.null(curr)) curr <- vector()
  daughters <- tree$edge[which(tree$edge[,1]==node),2]
  curr <- c(curr,daughters)
  w <- which(daughters>=length(tree$tip))
  
  if(length(w) > 0) for(i in 1:length(w)) 
    curr <- getDescendants(tree,daughters[w[i]], curr)
  
  return(curr)
}

x<-c("#' getDescendants",
     "#'",
     "#' This function allows you to get descendant node id of the target node.",
     "#' @param tree phylo object.",
     "#' @param node node id number",
     "#' @param curr Default to NULL",
     "#' @keywords phytools",
     "#' @export",
     "#' @examples",
     "#' acaena <- read.nexus(NZ_Acaena_BEAST_output_6gene.tree)",
     "#' getDescendants(acaena, node = 20)",
     "")

write(x, file =  
        paste(packagename, "/R/getDescendants.R", sep = "")
)

dump(getDescendants, file = paste(packagename, '/R/getDescendants.R', sep = ""), append = TRUE)

#######################################################
### Function 3
#######################################################

# Get node numbers of all internal nodes
GetInternalNodeNumber <- function(tree){
  
  innode <- (length(tree$tip.label) + 1):max(tree$edge)
  
  for(i in innode){
    sis <- getSisters(tree, i)
    
    if(sum(sis) == 0){
      res[[i]] <- "the farmost ancestor"
    }else{
      desnode <- getDescendants(tree, node = sis)
      destip <- desnode[desnode %in% nodes$nodelabel]
      res[[i]] <- destip 
    }
  }
  
  return(res)
  
}

x <- c("#' GetInternalNodeNumber",
     "#'",
     "#' This function gets node numbers of all internal nodes in the tree.",
     "#' @param tree phylo object.",
     "#' @keywords phytools",
     "#' @export",
     "#' @examples",
     "#' acaena <- read.nexus(NZ_Acaena_BEAST_output_6gene.tree)",
     "#' GetInternalNodeNumber(acaena)",
     "")


write(x, file =  
        paste(packagename, "/R/GetInternalNodeNumber.R", sep = "")
)

dump(GetInternalNodeNumber, file = paste(packagename, '/R/GetInternalNodeNumber.R', sep = ""), append = TRUE)

#######################################################
### Function 4
#######################################################
## Find sister node of a target species

findSisterNode <- function(tree){
  
  # Extract names of edges (i.e. tips and taxa)
  tips <- tree$tip.label
  
  ## first get the node numbers of the tips
  nodes <- data.frame(sapply(tips, function(x,y) which(y == x), y = tree$tip.label))
  colnames(nodes) <- "nodelabel"
  
  ## Finally, find sister node of a target species
  res <- list()
  
  for(i in nodes$nodelabel){
    sis <- getSisters(tree,i)
    
    desnode <- getDescendants(tree, node = sis) 
    
    if(sum(desnode) == 0){
      res[[i]] <- sis
    }else{
      destip <- desnode[desnode %in% nodes$nodelabel]
      
      res[[i]] <- destip 
    }
    
  }
  
  return(res)
}

x<-c("#' findSisterNode",
     "#'",
     "#' This function gets sister node id of all nodes.",
     "#' @param tree phylo object.",
     "#' @keywords phytools",
     "#' @export",
     "#' acaena <- read.nexus(NZ_Acaena_BEAST_output_6gene.tree)",
     "#' findSisterNode(acaena)",
     "")

write(x, file =  
        paste(packagename, "/R/findSisterNode.R", sep = "")
)

dump(findSisterNode, file = paste(packagename, '/R/findSisterNode.R', sep = ""), append = TRUE)









# Remove current function
setwd(paste("./", packagename, sep=""))
      
document()

rm(functionname)
detach(paste("package:", packagename, sep = ""), unload=TRUE)


# Install created package

setwd("..")
install(packagename)

