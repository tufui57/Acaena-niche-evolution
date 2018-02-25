


list_spname <- function(dat, # data frame with species names in its column names
                          genus_name, # genus name
                          separation # syntax separating genus and species name that you want to have in the returned object
                          ){
  spnames <- (colnames(dat) %>% grepl(genus_name, .)) %>% colnames(dat)[.]
  
  return(spnames)
}

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

makeTag_separate <- function(data, # vector of full species names 
                             genus_name, # genus name should be got rid of
                             separate # syntax which separates genus and species names e.g. "_" between "Acaena_agnipila"
){
  ### Import species name codes
  spname <- grepl(genus_name, data) %>% data[.]
  codes <- gsub(paste(genus_name, separate, sep=""), "", spname) %>% 
    gsub(paste("subsp.", separate, sep = ""), "", .) %>% 
    gsub(paste("var.", separate, sep = ""), "", .)
  
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

get_nodeID_from_spname <- function(spname, # species name
                       tree
                       ){
  tips <- tree$tip.label
  
  ## first get the node numbers of the tips
  nodes <- data.frame(sapply(tips, function(x,y) which(y == x), y = tree$tip.label))
  colnames(nodes) <- "nodelabel"
  
  rownames(nodes)<- clean_speciesname(rownames(nodes))
  
  nodeID <- (rownames(nodes) == spname) %>% nodes[.,]

  return(nodeID)
  }


list_sisterSpPairs <- function(tree # Phylogeny tree
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


get_sisterSpNames <- function(node, # Node number or species name
                              tree
                              ){
  if(is.character(node)){
    node <- get_nodeID(node, tree)
  }
  
  sislist <- findSisterNode(tree)
  sis <- c(sislist[node], get_spname_from_nodeID(sislist[node], tree))
  return(sis)
}


