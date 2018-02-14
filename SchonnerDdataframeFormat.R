
##########################################
# Schonner's D dataframe formatting
##########################################

SchonnerDdataframeFormat <- function(scho){
  
  scholist <- lapply(1:length(scho), function(i){
    # Convert list to dataframe
    x <- data.frame(scho[[i]])
    
    # Add colnames
    colnames(x) <- c("ecospat.corrected", "ecospat.uncorrected")
    
    # Get Shoenner's D. I don't use I.
    return(x["D",])
    
  }
  )
  
  # Add column in schonner's D dataframe showing clade numbers
  scholist <- do.call(rbind, scholist)
  nodeNo <- sapply(cladedata, "[[", 5) %>% strsplit(., "\ ")
  
  scholist$node1 <- as.numeric(do.call(rbind, lapply(nodeNo, "[[",1)))
  scholist$node2 <- as.numeric(do.call(rbind, lapply(nodeNo, "[[",2)))
  
  return(scholist)
  
  
}

