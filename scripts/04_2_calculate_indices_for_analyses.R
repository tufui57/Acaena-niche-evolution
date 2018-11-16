########################################################
### Calculate indices for analysis
########################################################

# Land cover change data
gd <- read.csv(paste(".//", genus_name, "_landscapeChangeHistory.csv", sep=""))

# Add land cover type "others"
gd$nonF.others <- (gd$nonF.EF + gd$nonF.nonPotentialHabitat)
gd$NF.others <- (gd$NF.EF + gd$NF.nonPotentialHabitat)

# Replace NA with 0
gd[is.na(gd$NF.NF)] <- 0

geoIncTable <- function(d, filepath){
  ### Calculate habitats increase
  
  # Proportion of secondary open habitat = occ in secondary habitat / occ in primary and secondary habitat
  d$proportionSecondaryHabitat <- (d$NF.nonF / (d$NF.nonF + d$nonF.nonF))
  
  # total = no. of all occurrence cells
  d$total <- d$number.of.1km.occurrence.cells
  d$log10.total <- log10(d$total)
  
  # Preference for open habitat = sum of occ in open (NF.nonF + nonF.nonF) / (forest + open)
  d$PreferenceOpen <- rowSums(d[, c("NF.nonF", "nonF.nonF")]) / rowSums(d[,c("NF.nonF", "NF.NF", "nonF.nonF", "nonF.NF")])
  
  write.csv(d[, c("X", # species names
                   "total", "log10.total",
                   "proportionSecondaryHabitat", 
                   "PreferenceOpen"
  )],
  file = filepath)
  
  return(d)
}

geoIncTable(gd, paste(".//", genus_name, "_data.csv", sep=""))

########################################################
### Add species age, niche volume/overlap
########################################################

gd <- read.csv(paste(".//", genus_name, "_data.csv", sep = ""))

# species age and niche volume
if(file.exists(paste("Y:\\NicheVolume_age_", genus_name,".csv", sep="")) == FALSE){
  source(".//Chionochloa niche evolution//09_2_DataPreparation_for_Analysis.R")
}
ages <- read.csv(paste("Y:\\NicheVolume_age_", genus_name,".csv", sep=""))
colnames(gd)[colnames(gd) == "X"] <- "spname"

# Modify species name of species age dataframe
a <- gsub("subsp", "subsp.", ages$X)
# Chionochloa flaviscans in phylogeny tree and Chionochloa flaviscans f. flaviscans in occurrence records are treated as the same species?
ages$X <- gsub("Chionochloa_flavicans", "Chionochloa_flavicans_f._flavicans", gsub("var", "var.", a))


# Extract just species age. Because "ages" object doesn't have niche volume data for species which doesn't appear in phylogeny trees
gd$spname <- as.character(gd$spname)
ages$spname <- as.character(ages$spname)
ages2 <- ages[!is.na(ages$spname),]

gd$speciesAge <- sapply(gd$spname, function(x){
  if(((x == ages2$spname) %>% sum) == 0){
    return(NA)
  }else{  
    a <- ages2[x == ages2$spname, "speciesAge"]
    return(a)
  }

  }
  )
            
write.csv(gd, paste(".//", genus_name, "_data_analyses.csv", sep = ""))
