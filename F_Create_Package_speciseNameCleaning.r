
###############################################################################################
### Get a list of species names from data frame with species names in its column names
###############################################################################################

list_spname <- function(dat, # data frame with species names in its column names
                          genus_name, # genus name
                          separation # syntax separating genus and species name that you want to have in the returned object
                          ){
  spnames <- (colnames(dat) %>% grepl(genus_name, .)) %>% colnames(dat)[.]
  
  return(spnames)
}

#############################################################################
### Find species which has no occurrence records
#############################################################################

list_spWithNoOcc <- function(dat, # data frame of occurrence records
                        genus_name # genus name
                        ){
  spnames <- (colnames(dat) %>% grepl(genus_name, .)) %>% colnames(dat)[.]
  
  noocc <- (colSums(data[,spnames]) == 0) %>% spnames[.]
  
  return(noocc)
}

#############################################################################
### Clean species names
#############################################################################

# Syntax separating between subspecies, variant and form will be changed to "_", e.g., "Acaena_tesca_subsp._tesca"  
# Gene tags must be removed from sepcies names
clean_speciesname <- function(spnames # vector of character
){
  ### Acaena
  # Modify species names in phylogentic distance file
  a <- unlist(strsplit(spnames, "_Ac"))
  a2 <- gsub("_EU352216", "", a) %>% gsub("_AY634821", "", .) %>% gsub("novae-", "novae.", .)
  spnames <- grepl("_", a2) %>% a2[.]
  
  ### Chionochloa
  # Modify names
  spnames <- gsub("subsp", "subsp.", spnames) %>% 
    gsub("var_", "var._", .) %>% 
    gsub("Chionochloa_flavicans", "Chionochloa_flavicans_f._flavicans", .) %>% 
    gsub("Chionochloa_rubra_subsp._rubra", "Chionochloa_rubra_var._rubra", .)
  
  return(spnames)
}

#############################################################################
### Create species name tag from species name
#############################################################################

makeTag_separate <- function(data, # vector of full species names 
                             genus_name, # genus name should be got rid of
                             separate # syntax which separates genus and species names e.g. "_" between "Acaena_agnipila"
){
  ### Import species name codes
  spname <- grepl(genus_name, data) %>% data[.]
  codes <- gsub(paste(genus_name, separate, sep = ""), "", spname) %>% 
    gsub(paste("subsp.", separate, sep = ""), "", .) %>% 
    gsub(paste("var.", separate, sep = ""), "", .) %>%     
    gsub(paste("subsp", separate, sep = ""), "", .) %>% 
    gsub(paste("var", separate, sep = ""), "", .)
  
  spname <- (codes %>% substring(., 1, last = 3) %>% mutate(as_tibble(spname), tag = .))
  colnames(spname)[1] <- "X"
  subsp <- codes %>% 
    strsplit(., separate) %>% 
    lapply(., function(x){
      ifelse(is.na(x[2]), "", x[2])
    }) %>% 
    substring(., 1, last = 3)
  
  spname <- lapply(1:length(subsp), function(i){
    paste(spname[i,"tag"], subsp[i], sep = "_")
  }
  ) %>% unlist %>% 
    gsub(paste(separate, "$", sep=""), "", .) %>% 
    mutate(spname, tag = .) 
  
  return(spname)
  
}

#############################################################################
### Get the species name from node number
#############################################################################

get_spname_from_nodeID <- function(node, # node ID number
                                   tree
){
  tips <- tree$tip.label
  
  ## first get the node numbers of the tips
  nodes <- data.frame(sapply(tips, function(x,y) which(y == x), y = tree$tip.label))
  colnames(nodes) <- "nodelabel"
  
  rownames(nodes)<- clean_speciesname(rownames(nodes))
  
  nodeName <- (nodes == node) %>% rownames(nodes)[.]
  
  return(nodeName)
}

#############################################################################
### Get the node number from species name
#############################################################################

get_nodeID_from_spname <- function(spname, # species name
                       tree
                       ){
  # Get a list of tip (species) names from tree
  tips <- tree$tip.label
  
  ## Get the node numbers of the tips
  nodes <- data.frame(sapply(tips, function(x,y) which(y == x), y = tree$tip.label))
  colnames(nodes) <- "nodelabel"
  
  rownames(nodes)<- clean_speciesname(rownames(nodes))
  
  nodeID <- (rownames(nodes) == spname) %>% nodes[.,]

  return(nodeID)
}

#############################################################################
### Get the node number of the closest ancestor (parent) of the node
#############################################################################

library(phylobase)
get_parentNodeID <- function(node, # Node number or species name. Both work.
                                    tree
){
  
  # Get "phylo4" object from "phylo" object
  tree2 <- extractTree(tree)
  return(ancestor(tree2, node))
}


########################################################
### Get a list of sister species pairs in the tree
########################################################

list_sisterSpPairs <- function(tree # Phylogeny tree object
                                   ){
  # find sister node of a target species
  tipssister <- findSisterNode(tree)

  sistersp <- sapply(tipssister, function(x){
    if(length(x) == 1){
      return(x)
    }else{
      return("NA")
    }
  }
  )
  
  sistersp2 <- cbind(1:length(tipssister), as.numeric(sistersp))
  colnames(sistersp2) <- c("nodes", "sisterNodes")
  return(sistersp2)

  }

########################################################
### Get a name of sister species of the node
########################################################

get_sisterSpNames <- function(node, # Node number or species name
                              tree
                              ){
  if(is.character(node)){
    node <- get_nodeID_from_spname(node, tree)
  }
  
  sislist <- findSisterNode(tree)
  sis <- c(sislist[node], get_spname_from_nodeID(sislist[node], tree))
  return(sis)
  
}

########################################################
### Count the number of species within the clade
########################################################

count_spRichness <- function(i, # Node number
                          tree
){
  # If the target node has no descendant species, i.e. the node is a terminal tip/node
  if( getDescendants(tree, node = i) %>% length <= 1 ){
    
    print("The target node is terminal node")
    return(1)
    
  }else{
    
    sprich <- getDescendants(tree, node = i) %in% 1:length(tree$tip.label) %>% sum
    return(sprich)
  }
}

########################################################
### Get a node ID of the closest clade (aunt)
########################################################

get_auntNodeID <- function(node, # Node number or species name. Both work.
                             tree
){
  parentnode <- get_parentNodeID(node, tree)
  aunt <- getSisters(tree, parentnode)
  return(aunt)
}
