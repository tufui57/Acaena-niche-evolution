
source("06_Clade pairing.R")

generateClimateDataOftargetNode <- function(i, # internal node number
                                        acaena, # tree object
                                        allnodesister, # List of discendant nodes of its sister node
                                        scores, # background data containing occurrence data of target taxa
                                        nodes = nodes,
                                        tips = tips
                                        ){

  
  # If the target node has no descendant species, i.e. the node is a terminal tip/node
  if( getDescendants(acaena, node = i) %>% sum == 0 ){
    
    # If the target node has no occurrence records
    if(colnames(scores) %in% rownames(nodes)[i] %>% sum == 0){
      
      stop("The target node has no occurrence records")
      
    }else{
      descendantColumn <- (colnames(scores) == rownames(nodes)[i])
      scores$targetClade <- rownames(nodes)[i] %>% scores[, .]
      
    }
  }else{
    
    descendants <- getDescendants(acaena, node = i) %>% rownames(nodes)[.]
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




generateClimateDataOfClades <- function(i, # internal node number
                                        acaena, # tree object
                                        allnodesister, # List of discendant nodes of its sister node
                                        scores # background data containing occurrence data of target taxa
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
  
  # If the target node has no descendant species, i.e. the node is a terminal tip/node
  if( getDescendants(acaena, node = i) %>% sum == 0 ){
    
    # If the target node has no occurrence records
    if(colnames(scores) %in% rownames(nodes)[i] %>% sum == 0){
      
      stop("The target node has no occurrence records")
      
    }else{
      descendantColumn <- (colnames(scores) == rownames(nodes)[i])
      scores$targetClade <- rownames(nodes)[i] %>% scores[, .]
      
    }
  }else{
    
    descendants <- getDescendants(acaena, node = i) %>% rownames(nodes)[.]
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
  
  ################################################
  ## Find columns of species in the sister node
  ################################################
  
  sisdescendants <- allnodesister[[i]] %>% rownames(nodes)[.]
  sisdescendantColumn <- colnames(scores) %in% sisdescendants
  
  ## Create a column showing clade occurrence records
  if(sum(sisdescendantColumn) > 1){
    
    scores$sisClade <- ifelse(scores[,sisdescendantColumn] %>% rowSums >= 1,
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
  nodeName = as.character(codes2[codes2$X %in% colnames(scores)[descendantColumn], "tag"])
  sisnodeNumber = distance2[i, "sisterNode"]
  sisnodeName = as.character(codes2[codes2$X %in% rownames(nodes)[allnodesister[[i]]], "tag"])
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
