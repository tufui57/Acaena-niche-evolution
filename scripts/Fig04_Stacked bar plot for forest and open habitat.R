##########################################################################
### Stacked bar plot to show that Acaena is open habitat plants
##########################################################################

genus_name="Acaena"

library(reshape2)
library(ggplot2)
library(dplyr)

source(".//functions//F_speciseNameCleaning_spnameFromPhylogenyTree.r")

##########################################################################
### Data preparation
##########################################################################

# geographical landscape change
landcover <- read.csv(paste("Y://1st chpater_Acaena project//meta data//", genus_name, "_bioclim_landcover_history_worldclim1_1km.csv", sep = ""))


# Replace NA with 0
landcover <- mutate(landcover, landCoverChange = as.character(landcover$landCoverChange))
landcover2 <- landcover[!(as.character(landcover$landCoverChange) %>% is.na), ]

dat <- matrix(nrow = 4, ncol = sum(grepl("^Acaena", colnames(landcover2)))
              )
dat <- data.frame(dat)

for(i in colnames(landcover2)[grep("^Acaena", colnames(landcover2))]){
 dat[,i] <- c(
   ((landcover2[landcover2$landCoverChange == "NF-NF", i] == 1) %>% sum(., na.rm = T)),
   ((landcover2[landcover2$landCoverChange == "NF-nonF", i] == 1) %>% sum(., na.rm = T)),
   ((landcover2[landcover2$landCoverChange == "nonF-nonF", i] == 1) %>% sum(., na.rm = T)),
   ((landcover2[landcover2$landCoverChange == "nonF-NF", i] == 1) %>% sum(., na.rm = T))
 )
}

rownames(dat) <- c("NF.NF", "NF.nonF", "nonF.nonF", "nonF.NF")

d <- rbind(dat,
      # Current open habitat
      colSums(dat[c("NF.nonF", "nonF.nonF"), ]),
      # Current forest excluding exotic forest ("nonF.EF", "NF.EF")
      colSums(dat[c("NF.NF", "nonF.NF"), ]),
      # Total occurrences for indice calculation
      colSums(dat)
      )
d2 <- d[c("5", "6","nonF.nonF","NF.nonF","7"),19:ncol(d)]
rownames(d2) <- c("open", "forest", "primary", "secondary", "total")

ratio <- rbind(d2["open",]/d2["total",], d2["forest",]/d2["total",], d2["secondary",]/(d2["primary",] + d2["secondary",]))
rownames(ratio) <- c("open/total", "forest/total", "ProportionOfSecondaryOpen")

data <- t(rbind(d2,ratio)) %>% data.frame

# Subset data by species names
trait <- read.csv("Y://1st chpater_Acaena project//raw data//traits.csv")

data2 <- data[rownames(data) %in% trait$X, ]

write.csv(data2, file = "Y://ProportionOfSecondary11Oct.csv")

##########################################################################
### Prepare data for stacked bar plot
##########################################################################

data <- read.csv("Y://ProportionOfSecondary11Oct.csv")

### Show how much ratio forest and open habitat have
pl <- data[order(data$open.total, decreasing = T),]

# # Put NA in proportion of secondary open habitat which the species has no occurences in primary open habitat
# pl[pl$primary == 0, "ProportionOfSecondaryOpen"] <- NA

# Species name shown in the figure should be name tag.
tag <- makeTag_separate(pl$X, genus_name, "_") %>% data.frame
pl$spname <- toupper(tag$tag)

# Convert data for drowning stacked bar plot
pl.m <- melt(pl[, c("spname", "open.total", "forest.total")], id.vars="spname")
# Data to add as points on the bar plot
point <- pl[, c("spname","ProportionOfSecondaryOpen")]


# Change order of blocks in bars
pl.m$variable <- relevel(pl.m$variable, "forest.total")

pl.m$name <- factor(pl.m$spname, levels = pl.m$spname[order(pl$forest.total)])
plotStacked <- ggplot(pl.m) +
  geom_bar(aes(x = name, y = value, fill = variable), stat = "identity") +
  scale_fill_manual(values = c("gray60","gray90"),
                    name = "",
                    # Name of each legend factor. 
                    # This must be same factors as factors in "colname" of ggplot(aes(colour = colname)), otherwise no legend will be drawn.
                    breaks = c("forest.total", "open.total"),
                    label = c("Forest", "Open")
                    ) +
  ylab("Proportion") +
  xlab("Acaena species") +
  geom_point(data = point, aes_string(x = "spname", y = "ProportionOfSecondaryOpen"), size = 2, colour = "black") +
  theme(panel.background = element_rect(fill='white'),
        text = element_text(size = 20),
        axis.text = element_text(size = 20),
        axis.text.x = element_text(angle = 45, hjust = 1)
        )

ggsave("Y://stacked_bar.png", plotStacked, width = 300, height = 200, units = "mm")



