###########################################################################
### Get climate data of single species from node ID number
###########################################################################


generateClimateDataOfTargetNode <- function(i, # Node number
                                        tree, # tree object
                                        allnodesister, # List of discendant nodes of its sister node
                                        scores, # background data containing occurrence data of target taxa
                                        nodes,
                                        tips
                                        ){
  
  # If the target node has no descendant species, i.e. the node is a terminal tip/node
  if( getDescendants(tree, node = i) %>% length <= 1 ){
    
    # If the target node has no occurrence records
    if(colnames(scores) %in% rownames(nodes)[i] %>% sum == 0){
      
      stop("The target node has no occurrence records")
      
    }else{
      descendantColumn <- (colnames(scores) == rownames(nodes)[i])
      scores$targetClade <- rownames(nodes)[i] %>% scores[, .]
      
    }
  }else{
    
    descendants <- getDescendants(tree, node = i) %>% rownames(nodes)[.]
    descendantColumn <- colnames(scores) %in% descendants
    
    if(sum(descendantColumn) > 1){
      
      # if the target node has multiple descendant species
      # Create a column showing clade occurrence records
      scores$targetClade <- ifelse(scores[, descendantColumn] %>% rowSums >= 1, 1, 0)
      
    }else{
      
      # if the target node has just one descendant
      scores$targetClade <- scores[, descendantColumn]
    }
    
  }
  
  cladeScore <- scores[, c("PC1", "PC2", "targetClade")]
  
  return(cladeScore)
}




###########################################################################
### Get climate data of all species within the node from node ID number
###########################################################################

generateClimateDataOfClades <- function(i, # node number
                                        tree, # tree object
                                        allnodesister, # List of discendant nodes of its sister node
                                        scores, # background data containing occurrence data of target taxa
                                        nodes,
                                        tips,
                                        spnameCodes
                                        ){
  
  if(i %>% is.numeric == FALSE){
    stop("Use number for i, not character")
  }
  
  ### Print node status
  print(paste("Target node is", i, rownames(nodes)[i]))
  print(paste("Sister node is", distance2[i,"sisterNode"]))
  print("Sister node contains")
  cat(paste(rownames(nodes)[allnodesister[[i]]], collapse  = "\n"))
  
  ################################################
  ## Find columns of species in the target node
  ################################################
  
  cladeScore <- generateClimateDataOfTargetNode(i, tree, allnodesister, scores, nodes, tips)
  
  ################################################
  ## Find columns of species in the sister node
  ################################################
  
  sisdescendants <- allnodesister[[i]] %>% rownames(nodes)[.]
  sisdescendantColumn <- colnames(scores) %in% sisdescendants
  
  ## Create a column showing clade occurrence records
  if(sum(sisdescendantColumn) > 1){
    
    scores$sisClade <- ifelse(scores[, sisdescendantColumn] %>% rowSums >= 1,
                              1, 0)
  
    }else{
          if(sum(sisdescendantColumn) == 0){
            
            cat(paste("Sister species", allnodesister[[i]] %>% rownames(nodes)[.], "has no occurrence records."))
            stop("Script stops!")
      
    }else{
      scores$sisClade <- scores[, sisdescendantColumn]
    }
  }
  sisCladeScore <- scores[,c("PC1", "PC2", "sisClade")]
  
  ################################################
  ## Put all results in list object
  ################################################
  
  nodeNumber = i
  # Species name codes
  # Get descendant species name of the node
  if( getDescendants(tree, node = i) %>% length <= 1 ){
    
    descendantColumn <- (colnames(scores) == rownames(nodes)[i])
    
    }else{
    
    descendants <- getDescendants(tree, node = i) %>% rownames(nodes)[.]
    descendantColumn <- colnames(scores) %in% descendants
    
  }
  nodeName = pull(spnameCodes[spnameCodes$X %in% colnames(scores)[descendantColumn], ], tag)
  
  sisnodeNumber = distance2[i, "sisterNode"]
  sisnodeName = pull(spnameCodes[spnameCodes$X %in% rownames(nodes)[allnodesister[[i]]], ], tag)
  cladedata1 <- (cladeScore$targetClade == 1) %>% cladeScore[.,] # Clade 1 PCA data
  cladedata2 <- (sisCladeScore$sisClade == 1) %>% sisCladeScore[., ]  # Clade 2 PCA data
  
  # Output results
  clades <- list()
  clades[[1]] <- cladedata1
  clades[[2]] <- nodeName
  
  clades[[3]] <- cladedata2
  clades[[4]] <- sisnodeName
  
  clades[[5]] <- paste(nodeNumber, sisnodeNumber)
  
  return(clades)
  
  
}

##########################################################################################################################
### Get climate data of aunt node of sister species from the node ID number of target sister species
##########################################################################################################################

get_climatedata_of_auntNode <- function(i, # node ID of target sister species
                                        tree, # tree object
                                        scores # background data containing occurrence data of target taxa
                                        ){
  
  ## Find parent node of target sister species pair
  ancestor <- tree$edge[which(i == tree$edge[, 2])]
  ancestorSisNode <- distance2[distance2[, "node"] == ancestor, "sisterNode"]
  
  print(paste("Sister node of parent node of target sister species pair is", ancestorSisNode))
  
  ## Niche of sister node of ancestor node
  ansSisScore <- generateClimateDataOfTargetNode(
    ancestorSisNode, tree, allnodesister, scores, nodes, tips)
  ansSisScore2 <- ansSisScore[ansSisScore[, "targetClade"] == 1, c("PC1", "PC2", "targetClade")]
  
  return(ansSisScore2)
}
