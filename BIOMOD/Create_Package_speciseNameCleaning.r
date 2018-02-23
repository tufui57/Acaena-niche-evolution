


getSpnameList <- function(dat, # data frame with species names in its column names
                          genus_name, # genus name
                          separation # syntax separating genus and species name that you want to have in the returned object
                          ){
  spnames <- (colnames(dat) %>% grepl(genus_name, .)) %>% colnames(dat)[.]
  
}


clean_speciesname <- function(spname # vector of character
                              ){
  
}


get_nodeID <- function(){
  
  
}


get_sisternodeID <- function(){
  
  
}


get_sisterSpPairs <- function(){
  
  
}

#chion
sispairs <- c(9,29,20,33,15,12,4,6,31,30)

overlapPdData[overlapPdData$node1 %in% sispairs,]




# Show sister species list
for(i in overlapPdData$node1){
  print(i)
  print(rownames(nodes)[allnodesister[[i]]])
}
