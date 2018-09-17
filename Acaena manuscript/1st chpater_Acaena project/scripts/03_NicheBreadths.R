
library(ggplot2)
library(devtools)

### Import PCA scores
load(paste(".\\Scores_", genus_name, "_landcover13sep.data", sep = ""))

########################################################
### Quantify species climate niche breadths
########################################################

show_PCaxes_range <- function(spname,
         scores){
  spscores <- scores[scores[, spname] == 1, c("PC1", "PC2")]
  
  print(spname)
  
  sphist <- hist(spscores$PC1)
  
  maxPC1 <- which(sphist$counts == max(sphist$counts))
  
  print(paste("Maximum counts in 7 ranges of PC1 \ between",
              sphist$breaks[maxPC1], "and", sphist$breaks[maxPC1 + 1], " \ is",
              sphist$counts[maxPC1])
  )
  
  print(
    paste("Median of PC1 is", median(spscores$PC1)
    )
  )
  
  
  sphist2 <- hist(spscores$PC2)
  
  maxPC2 <- which(sphist2$counts == max(sphist2$counts))
  
  print(paste("Maximum counts in 7 ranges of PC2 \ between",
              sphist$breaks[maxPC2], "and", sphist$breaks[maxPC2 + 1], " \ is",
              sphist$counts[maxPC2])
  )
  
  print(
    paste("Median of PC2 is", median(spscores$PC2)
          )
    )
  
  counts <- cbind(sphist$breaks, sphist$counts, sphist2$breaks, sphist2$counts)
  colnames(counts) <- c("PC1_breaks", "PC1_counts", "PC2_breaks", "PC2_counts")
  return(counts)
  
}
spname <- colnames(scores)[grepl("^Acaena", colnames(scores))]
lapply(spname, show_PCaxes_range, scores)

med1 <- lapply(spname, function(x){
  spscores <- scores[scores[, x] == 1, c("PC1", "PC2")]
  sphist <- hist(spscores$PC1)
  return(median(spscores$PC1))
}
)

med2 <- lapply(spname, function(x){
  spscores <- scores[scores[, x] == 1, c("PC1", "PC2")]
  sphist <- hist(spscores$PC2)
  return(median(spscores$PC2))
}
)

# Combine PC medians
med.data <- cbind(spname, unlist(med1), unlist(med2))
colnames(med.data) <- c("spname", "median.of.temp", "median.of.prec")
write.csv(med.data, file = paste(".//", genus_name, "_nicheBreadths.csv", sep = ""))


