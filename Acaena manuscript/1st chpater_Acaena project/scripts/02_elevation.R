########################################################################################
### Mean elevation of occurrence records
########################################################################################

genus_name="Acaena"

### Import data
alt <- read.csv(".\\Acaena niche evolution\\Acaena manuscript\\1st chpater_Acaena project\\meta data\\Acaena_bio_alt.csv")

# sp names
spname <- grepl(paste("^", genus_name, sep=""), colnames(alt)) %>% colnames(alt)[.]

ras <- project_and_convert_occurrencePoints_to_raster(alt, refWGS = pre, val = "NZL1_alt")

# Make raster stack
bio_land2 <- stack(c(bio_land, ras[[4]]))

