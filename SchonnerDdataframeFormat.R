
##########################################
# Schonner's D dataframe formatting
##########################################

SchonnerDdataframeFormat <- function(scho, # List of correlated/uncorrelated Schoenner's D and I
                                     colname = c("ecospat.corrected.D", "ecospat.uncorrected.D"), # names of columns
                                     chooseIndex = "D" # choose overlap index from D and I
                                     ){

  scholist <- lapply(1:length(scho), function(i){
    
    ### If list has null objects
    
    if(scho[[i]] %>% length == 0){
      
      return(c(999, 999))
      
    }else{
      
      # Convert list to dataframe
      x <- data.frame(scho[[i]])
      
      # Add colnames
      colnames(x) <- colname
      
      # Get the index you need
      return(x[chooseIndex, ])
      
    }
    }
  )
  
  # Add column in schonner's D dataframe showing clade numbers
  scholist <- do.call(rbind, scholist)
  nodeNo <- sapply(cladedata, "[[", 5) %>% strsplit(., "\ ")
  
  if(((scholist[colname[1]] == 999) %>% sum) == 0){
    
    # If there is a NA row, the shoenner's D list is niche volume list.  
    scholist$node1 <- as.numeric(do.call(rbind, lapply(nodeNo, "[[",1)))
    scholist$node2 <- as.numeric(do.call(rbind, lapply(nodeNo, "[[",2)))
    
  }else{
    scholist$node1 <- 1:length(scho)
  }

  
  return(scholist)
  
}
