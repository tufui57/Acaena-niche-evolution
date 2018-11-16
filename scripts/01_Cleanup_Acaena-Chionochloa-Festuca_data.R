## code created by: Gregory T. Nelson
## last edit: 11.05.2016
## code created and modified by: Miki Nomura since 12.05.2016
## accessing distributional data from Barb Anderson

## GTN notes: This code imports, checks, and formats a distributional data set for Acaena, Chionochloa, and Festuca species obtained from herbaria data, National Vegetation survey data, personal observations from Kelvin Lloyd, data from Barbara Anderson, and from the Global Biodiversity Information Facility.

# files imported: "cleandata_Acaena-Chionochloa-Festuca_BARB.csv"
#                 "cleandata_Acaena-Chionochloa-Festuca_GBIF.csv"
#                 "cleandata_Acaena-Chionochloa-Festuca_herbaria.csv"
#                 "cleandata_Acaena-Chionochloa-Festuca_NVS.csv"
# files exported: "rdata_Acaena_all-presence.csv"
#                 "rdata_Chinoochloa_all-occurances.csv"
#                 "rdata_Festuca_all-presence.csv"

# code edited by: NA
# last edit: NA
# editor notes: NA




##########################################################################################
####################   SET UP THE STUFF   ################################################
##########################################################################################


# clear the work space
rm(list = ls())


# set directory for MN student desktop
setwd(".//Acaena niche evolution//Acaena manuscript//1st chpater_Acaena project//raw data//Data from Greg//Greg distribution data")


# install some flash packages
#install.packages("rgbif")
#install.packages("dplyr")


# load all the cool libraries
# the command in these packages are written in C++ format
# library("plyr")
library("dplyr")


##########################################################################################
####################   IMPORT THE DATA   #################################################
##########################################################################################


### get the distributional data from BARB
data.BARB1 <- read.csv("cleandata_Acaena-Chionochloa-Festuca_BARB.csv")

### get the distributional data from herbaria
data.herb1 <- read.csv("cleandata_Acaena-Chionochloa-Festuca_herbaria.csv")

### get the distributional data from NVS
data.NVS1 <- read.csv("cleandata_Acaena-Chionochloa-Festuca_NVS.csv")



##########################################################################################
####################   FORMAT THE DATA   #################################################
##########################################################################################


###### BARB's data ######

# make a taxa column that matches the phylogenetic tree object taxa labels
# use nested if.else() statements to create the taxa name depending on what is in the subspecies, variety, and forma columns
# # make a species column that combines genus and species information for simple analysis 
data.BARB2 <- data.BARB1 %>%
  mutate(taxa = ifelse(is.na(subspecies) & is.na(variety) & is.na(forma), paste(genus, species, sep = "_"), ifelse(is.na(subspecies) & is.na(variety), paste(genus, species, "f.", forma, sep = "_"), ifelse(!is.na(variety), paste(genus, species, "var.", variety, sep = "_"), paste(genus, species, "subsp.", subspecies, sep = "_"))))) %>%
  mutate(species = paste(genus, species, sep = "_"))



# rename lon and lat
# select only columns of interest for modeling (no metadata)
# remove redundant data rows per taxa (plot data are assigned coordinates by plot, so co-occuring species will have the same latitude and longitudes)
# remove any records that are missing the desired data types
data.BARB3 <- data.BARB2 %>%
  mutate(lat = latitude) %>%
  mutate(lon = longitude) %>%
  distinct(taxa, lon, lat) %>%
  select(species, taxa, lon, lat) %>%
  na.omit()



# some taxa labels contain question marks, indicating that the identification has low confidence
# any taxa label with a question mark (?) is removed from the data set
data.BARB4 <- data.BARB3 %>%
  filter(!grepl("\\?", taxa))



# change the incorrect spelling some taxa labels
# "Chionochloa_australi" is not a real species and assumed to be a typing mistake for "Chinochloa_australis"
data.BARB5 <- data.BARB4 %>%
  mutate(taxa = ifelse(grepl("Chionochloa_australi", taxa), "Chionochloa_australis", taxa)) %>%
  mutate(species = ifelse(grepl("Chionochloa_australi", species), "Chionochloa_australis", species))



# seperate into three data files, one for each genus

data.BARB.Acaena <- data.BARB5 %>%
  filter(grepl("Acaena_*", species))

data.BARB.Chion <- data.BARB5 %>%
  filter(grepl("Chionochloa_*", species))

data.BARB.Festu <- data.BARB5 %>%
  filter(grepl("Festuca_*", species))







###### herbaria data ######

# make a taxa column that matches the phylogenetic tree object taxa labels
# use nested if.else statements to create the taxa name depending on what is in the subspecies, variety, and forma columns
# make a species column that combines genus and species information for simple analysis
data.herb2 <- data.herb1 %>%
  mutate(species = as.character(species)) %>%
  mutate(taxa = ifelse(is.na(subspecies) & is.na(variety) & is.na(forma), paste(genus, species, sep = "_"), ifelse(is.na(subspecies) & is.na(variety), paste(genus, species, "f.", forma, sep = "_"), ifelse(!is.na(variety), paste(genus, species, "var.", variety, sep = "_"), paste(genus, species, "subsp.", subspecies, sep = "_"))))) %>%
  mutate(species = paste(genus, species, sep = "_"))




# select only columns of interest for modeling (no metadata)
# remove redundant data rows per taxa (plot data are assigned coordinates by plot, so co-occuring species will have the same latitude and longitudes)
# remove any records that are missing the desired data types
data.herb3 <- data.herb2 %>%
  select(species, taxa, lon, lat) %>%
  distinct(taxa, lon, lat) %>%
  na.omit()



# change incorrect taxa labels
# "Chionochloa_crassuiscula" is not a valid name and is assumed to be a typing mistake for "Chionochloa_crassiuscula"
# "Acaena_anserinifolia_var._sericeinitens" is not a current name and has been revised to be "Acaena_anserinifolia"
data.herb4 <- data.herb3 %>%
  mutate(taxa = ifelse(grepl("Chionochloa_crassuiscula", taxa), "Chionochloa_crassiuscula", taxa)) %>%
  mutate(species = ifelse(grepl("Chionochloa_crassuiscula", species), "Chionochloa_crassiuscula", species)) %>%
  mutate(taxa = ifelse(grepl("Chionochloa_crassuiscula_subsp._crassuiscula", taxa), "Chionochloa_crassiuscula_subsp._crassiuscula", taxa)) %>%
  mutate(taxa = ifelse(grepl("Chionochloa_crassuiscula_subsp._directa", taxa), "Chionochloa_crassiuscula_subsp._directa", taxa)) %>%
  mutate(taxa = ifelse(grepl("Chionochloa_crassuiscula_subsp._torta", taxa), "Chionochloa_crassiuscula_subsp._torta", taxa)) %>%
  mutate(taxa = ifelse(grepl("Acaena_anserinifolia_var._sericeinitens", taxa), "Acaena_anserinifolia", taxa))




# seperate into three data files, one for each genus

data.herb.Acaena <- data.herb4 %>%
  filter(grepl("Acaena_", species))

data.herb.Chion <- data.herb4 %>%
  filter(grepl("Chionochloa_", species))

data.herb.Festu <- data.herb4 %>%
  filter(grepl("Festuca_", species))








###### NVS data ######

# make a taxa column that matches the phylogenetic tree object taxa labels
# use nested if.else() statements to create the taxa name depending on what is in the subspecies, variety, and forma columns
# make a species column that combines genus and species information for simple analysis
data.NVS2 <- data.NVS1 %>%
  mutate(taxa = ifelse(is.na(subspecies) & is.na(variety) & is.na(forma), paste(genus, species, sep = "_"), ifelse(is.na(subspecies) & is.na(variety), paste(genus, species, "f.", forma, sep = "_"), ifelse(!is.na(variety), paste(genus, species, "var.", variety, sep = "_"), paste(genus, species, "subsp.", subspecies, sep = "_"))))) %>%
  mutate(species = paste(genus, species, sep = "_"))



# some of the records have outdated species names, they will be updated
# also want to remove records that do not contain species level information (genus only)
# "Acaena_pusilla" and "Acaena_viridor" are old names and have been updated to "Acaena_ansernifolia"
# "Aceana_hirstula" is an old name and is now "Acaena_profundeincisa"
# "Aceana_sericeinitens" is an old name and probobly refers to "Acaena_profundeincisa"
# "Acaena_ovina" is an old name and is now "Acaena_agnipila"
# "Acaena_caesiiglauca_subsp._pilosa" is an old name and is now "Acaena_caesiiglauca"
# "Acaena_minor_subsp._minor" is incorrect as there are two varieties, not subspecies of A. minor
data.NVS3 <- data.NVS2 %>%
  mutate(taxa = ifelse(grepl("Acaena_hirsutula", taxa), "Acaena_profundeincisa", taxa)) %>% 
  mutate(taxa = ifelse(grepl("Acaena_sericeinitens", taxa), "Acaena_profundeincisa", taxa)) %>%
  mutate(taxa = ifelse(grepl("Acaena_pusilla", taxa), "Acaena_anserinifolia", taxa)) %>%
  mutate(taxa = ifelse(grepl("Acaena_viridior", taxa), "Acaena_anserinifolia", taxa)) %>%
  mutate(taxa = ifelse(grepl("Acaena_ovina", taxa), "Acaena_agnipila", taxa)) %>%
  mutate(species = ifelse(grepl("Acaena_hirsutula", species), "Acaena_profundeincisa", species)) %>%
  mutate(species = ifelse(grepl("Acaena_sericeinitens", species), "Acaena_profundeincisa", species)) %>%
  mutate(species = ifelse(grepl("Acaena_pusilla", species), "Acaena_anserinifolia", species)) %>%
  mutate(species = ifelse(grepl("Acaena_viridior", species), "Acaena_anserinifolia", species)) %>%
  mutate(species = ifelse(grepl("Acaena_ovina", species), "Acaena_agnipila", species)) %>%
  mutate(taxa = ifelse(grepl("Acaena_caesiiglauca_subsp._pilosa", taxa), "Acaena_caesiiglauca", taxa)) %>%
  mutate(taxa = ifelse(grepl("Acaena_minor_subsp._minor", taxa), "Acaena_minor_var._minor", taxa))



# select only columns of interest for modeling (no metadata)
# remove redundant data rows per taxa (plot data are assigned coordinates by plot, so co-occuring species will have the same latitude and longitudes)
# remove any records that are missing the desired data types
data.NVS4 <- data.NVS3 %>%
  select(species, taxa, lon, lat) %>%
  distinct(taxa, lon, lat) %>%
  na.omit()




# seperate into three data files, one for each genus

data.NVS.Acaena <- data.NVS4 %>%
  filter(grepl("Acaena_*", species)) %>%
  filter(!grepl("^Acaena_$", species))

data.NVS.Chion <- data.NVS4 %>%
  filter(grepl("Chionochloa_*", species)) %>%
  filter(!grepl("^Chionochloa_$", species))

data.NVS.Festu <- data.NVS4 %>%
  filter(grepl("Festuca_*", species)) %>%
  filter(!grepl("^Festuca_$", species))







##########################################
# combine for each genera


###### combine for Acaena ######

# combine the herbaria, NVS, and Barbara's data frames
data.Acaena1 <- do.call("rbind", list(data.herb.Acaena, data.NVS.Acaena, data.BARB.Acaena))


# remove redundant values per taxa
# arrange the dataset alphabetially by taxa
data.Acaena2 <- data.Acaena1 %>%
  distinct(taxa, lon, lat) %>%
  arrange(taxa)




###### combine for Chionochoa ######

# combine the herbaria, NVS, and Barbara's data frames
data.Chion1 <- do.call("rbind", list(data.herb.Chion, data.NVS.Chion, data.BARB.Chion))


# remove redundant values per taxa
# arrange the dataset alphabetially by taxa
data.Chion2 <- data.Chion1 %>%
  distinct(taxa, lon, lat) %>%
  arrange(taxa)



# edit occurance data following H.E. Conner (1991)
# remove data points that depart greatly from the described distributions (change taxa to NA)
# assign subspecies/variety levels when they are missing following described distributions
# arrange the dataset alphabetially by taxa
# remove redundant values per taxa
# remove NA values
data.Chion3 <- data.Chion2 %>%
  mutate(taxa = as.character(taxa)) %>%
  mutate(taxa = ifelse(taxa == "Chionochloa_acicularis" & lat > -43, "NA", taxa)) %>%
  mutate(taxa = ifelse(taxa == "Chionochloa_australis" & lat < -43.5, "NA", taxa)) %>%
  mutate(taxa = ifelse(taxa == "Chionochloa_australis" & lon > 173.5, "NA", taxa)) %>%
  mutate(taxa = ifelse(taxa == "Chionochloa_beddiei" & lat > -39, "NA", taxa)) %>%
  mutate(taxa = ifelse(taxa == "Chionochloa_bromoides" & lat < -40, "NA", taxa)) %>%
  mutate(taxa = ifelse(taxa == "Chionochloa_conspicua" & lat < -40.22116 & lon < 174.57457, "Chionochloa_conspicua_subsp._conspicua", taxa)) %>%
  mutate(taxa = ifelse(taxa == "Chionochloa_conspicua", "Chionochloa_conspicua_subsp._cunninghamii", taxa)) %>%
  mutate(taxa = ifelse(taxa == "Chionochloa_conspicua_subsp._cunninghamii" & lat < -40.22116 & lon < 174.57457, "Chionochloa_conspicua_subsp._conspicua", taxa)) %>%
  mutate(taxa = ifelse(taxa == "Chionochloa_conspicua_subsp._conspicua" & lon > 174.5614, "Chionochloa_conspicua_subsp._cunninghamii", taxa)) %>%
  mutate(taxa = ifelse(taxa == "Chionochloa_crassiuscula" & lat > -42, "NA", taxa)) %>%
  mutate(taxa = ifelse(taxa == "Chionochloa_crassiuscula" & lat < -46.7, "Chionochloa_crassiuscula_subsp._crassiuscula", taxa)) %>%
  mutate(taxa = ifelse(taxa == "Chionochloa_crassiuscula_subsp._crassiuscula" & lat > -46.7, "Chionochloa_crassiuscula", taxa)) %>%
  mutate(taxa = ifelse(taxa == "Chionochloa_crassiuscula", "Chionochloa_crassiuscula_subsp._torta", taxa)) %>%
  mutate(taxa = ifelse(taxa == "Chionochloa_crassiuscula_subsp._torta" & lon > 167.30 & lat < -45.63, "Chionochloa_crassiuscula_subsp._directa", taxa)) %>%
  mutate(taxa = ifelse(taxa == "Chionochloa_crassiuscula_subsp._torta" & lon > 167.94 & lat < -44.9, "Chionochloa_crassiuscula_subsp._directa", taxa)) %>%
  mutate(taxa = ifelse(taxa == "Chionochloa_defracta" & lon < 172, "NA", taxa)) %>%
  mutate(taxa = ifelse(taxa == "Chionochloa_flavescens" & lat < -45, "NA", taxa)) %>%
  mutate(taxa = ifelse(taxa == "Chionochloa_flavescens" & lon > 174.5614, "Chionochloa_flavescens_subsp._flavescens", taxa)) %>%
  mutate(taxa = ifelse(taxa == "Chionochloa_flavescens" & lat > -42.2657, "Chionochloa_flavescens_subsp._lupeola", taxa)) %>%
  mutate(taxa = ifelse(taxa == "Chionochloa_flavescens" & lon > 170, "Chionochloa_flavescens_subsp._brevis", taxa)) %>%
  mutate(taxa = ifelse(taxa == "Chionochloa_flavescens" & lon < 170, "Chionochloa_flavescens_subsp._hirta", taxa)) %>%
  mutate(taxa = ifelse(taxa == "Chionochloa_flavescens_subsp._brevis" & lat > -41.5, "NA", taxa)) %>%
  mutate(taxa = ifelse(taxa == "Chionochloa_flavescens_subsp._brevis" & lat < -45, "NA", taxa)) %>%
  mutate(taxa = ifelse(taxa == "Chionochloa_flavescens_subsp._flavescens" & lat < -42, "Chionochloa_flavescens_subsp._hirta", taxa)) %>%
  mutate(taxa = ifelse(taxa == "Chionochloa_flavescens_subsp._lupeola" & lon > 173, "Chionochloa_flavescens_subsp._flavescens", taxa)) %>%
  mutate(taxa = ifelse(taxa == "Chionochloa_flavescens_subsp._flavescens" & lat > -38, "NA", taxa)) %>%
  mutate(taxa = ifelse(taxa == "Chionochloa_flavescens_subsp._lupeola" & lat < -42.5, "NA", taxa)) %>%
  mutate(taxa = ifelse(taxa == "Chionochloa_flavicans" & lat < -40.5, "NA", taxa)) %>%
  mutate(taxa = ifelse(taxa == "Chionochloa_flavicans_f._flavicans" & lon > 176.9 & lat < -39.5, "Chionochloa_flavicans_f._temata", taxa)) %>%
  mutate(taxa = ifelse(taxa == "Chionochloa_flavicans", "Chionochloa_flavicans_f._flavicans", taxa)) %>%
  mutate(taxa = ifelse(taxa == "Chionochloa_juncea" & lon > 177, "NA", taxa)) %>%
  mutate(taxa = ifelse(taxa == "Chionochloa_lanea" & lat > -46, "NA", taxa)) %>%
  mutate(taxa = ifelse(taxa == "Chionochloa_macra" & lon > 175, "NA", taxa)) %>%
  mutate(taxa = ifelse(taxa == "Chionochloa_nivifera" & lat > -45, "NA", taxa)) %>%
  mutate(taxa = ifelse(taxa == "Chionochloa_oreophila" & lon > 168.9 & lat < -45, "NA", taxa)) %>%
  mutate(taxa = ifelse(taxa == "Chionochloa_pallens" & lon > 175, "Chionochloa_pallens_subsp._pallens", taxa)) %>%
  mutate(taxa = ifelse(taxa == "Chionochloa_pallens" & lon < 173.113173 & lat > -41.723047, "Chionochloa_pallens_subsp._pallens", taxa)) %>%
  mutate(taxa = ifelse(taxa == "Chionochloa_pallens" & lat > -41.586232, "Chionochloa_pallens_subsp._pallens", taxa)) %>%
  mutate(taxa = ifelse(taxa == "Chionochloa_pallens" & lat < -44, "Chionochloa_pallens_subsp._cadens", taxa)) %>%
  mutate(taxa = ifelse(taxa == "Chionochloa_pallens" & lat > -43.06667, "Chionochloa_pallens_subsp._pilosa", taxa)) %>%
  mutate(taxa = ifelse(taxa == "Chionochloa_pallens" & lon > 171.396206, "Chionochloa_pallens_subsp._pilosa", taxa)) %>%
  mutate(taxa = ifelse(taxa == "Chionochloa_pallens" & lon < 168.574162, "Chionochloa_pallens_subsp._cadens", taxa)) %>%
  mutate(taxa = ifelse(taxa == "Chionochloa_pallens" & lon > 170.285680, "Chionochloa_pallens_subsp._pilosa", taxa)) %>%
  mutate(taxa = ifelse(taxa == "Chionochloa_pallens" & lon < 170.285680, "Chionochloa_pallens_subsp._cadens", taxa)) %>%
  mutate(taxa = ifelse(taxa == "Chionochloa_pallens_subsp._cadens" & lat > -43, "NA", taxa)) %>%
  mutate(taxa = ifelse(taxa == "Chionochloa_pallens_subsp._cadens" & lat < -46, "NA", taxa)) %>%
  mutate(taxa = ifelse(taxa == "Chionochloa_pallens_subsp._cadens" & lat < -45.083224 & lon > 168.51370, "NA", taxa)) %>%
  mutate(taxa = ifelse(taxa == "Chionochloa_pallens_subsp._pallens" & lat < -42, "NA", taxa)) %>%
  mutate(taxa = ifelse(taxa == "Chionochloa_pallens_subsp._pilosa" & lat < -44.5, "NA", taxa)) %>%
  mutate(taxa = ifelse(taxa == "Chionochloa_pallens_subsp._pilosa" & lat > -41.524147 & lon < 173.146919, "NA", taxa)) %>%
  mutate(taxa = ifelse(taxa == "Chionochloa_rigida" & lat > -42, "NA", taxa)) %>%
  mutate(taxa = ifelse(taxa == "Chionochloa_rigida" & lat < -43.428703 & lon > 171.208116, "Chionochloa_rigida_subsp._rigida", taxa)) %>%
  mutate(taxa = ifelse(taxa == "Chionochloa_rigida" & lat < -45.305527 & lon > 168.747179, "Chionochloa_rigida_subsp._rigida", taxa)) %>%
  mutate(taxa = ifelse(taxa == "Chionochloa_rigida" & lat < -44.403834 & lon > 169.669159, "Chionochloa_rigida_subsp._rigida", taxa)) %>%
  mutate(taxa = ifelse(taxa == "Chionochloa_rigida" & lon < 167.865, "Chionochloa_rigida_subsp._amara", taxa)) %>%
  mutate(taxa = ifelse(taxa == "Chionochloa_rigida" & lat > -44.5535 & lon < 168.3097, "Chionochloa_rigida_subsp._amara", taxa)) %>%
  mutate(taxa = ifelse(taxa == "Chionochloa_rigida" & lat > -43.7341 & lon < 170.1431, "Chionochloa_rigida_subsp._amara", taxa)) %>%
  mutate(taxa = ifelse(taxa == "Chionochloa_rigida" & lat > -44.4039 & lon < 168.6515, "Chionochloa_rigida_subsp._amara", taxa)) %>%
  mutate(taxa = ifelse(taxa == "Chionochloa_rigida" & lat > -43.9246 & lon < 169.6770, "Chionochloa_rigida_subsp._amara", taxa)) %>%
  mutate(taxa = ifelse(taxa == "Chionochloa_rigida" & lat > -44.7260 & lon < 168.3097, "Chionochloa_rigida_subsp._amara", taxa)) %>%
  mutate(taxa = ifelse(taxa == "Chionochloa_rigida" & lat > -44.5 & lon < 169.100, "Chionochloa_rigida_subsp._amara", taxa)) %>%
  mutate(taxa = ifelse(taxa == "Chionochloa_rigida" & lat < -46.7, "Chionochloa_rigida_subsp._amara", taxa)) %>%
  mutate(taxa = ifelse(taxa == "Chionochloa_rigida_subsp._amara" & lat < -45.5 & lon > 168.8, "NA", taxa)) %>%
  mutate(taxa = ifelse(taxa == "Chionochloa_rigida_subsp._rigida" & lat < -45.65 & lon < 168.2, "NA", taxa)) %>%
  mutate(taxa = ifelse(taxa == "Chionochloa_rubra_subsp._rubra", "Chionochloa_rubra_var._rubra", taxa)) %>%
  mutate(taxa = ifelse(taxa == "Chionochloa_rigida", "NA", taxa)) %>%
  mutate(taxa = ifelse(taxa == "Chionochloa_rubra" & lat > -42 & lon > 173.8, "Chionochloa_rubra_var._rubra", taxa)) %>%
  mutate(taxa = ifelse(taxa == "Chionochloa_rubra" & lat < -45.1, "Chionochloa_rubra_subsp._cuprea", taxa)) %>%
  mutate(taxa = ifelse(taxa == "Chionochloa_rubra" & lat > -42.0 & lon < 173.0, "Chionochloa_rubra_subsp._occulta", taxa)) %>%
  mutate(taxa = ifelse(taxa == "Chionochloa_rubra" & lon > 173.15, "Chionochloa_rubra_var._rubra", taxa)) %>%
  mutate(taxa = ifelse(taxa == "Chionochloa_rubra_var._rubra" & lat > -39.4 & lon < 174.5, "Chionochloa_rubra_var._inermis", taxa)) %>%
  mutate(taxa = ifelse(taxa == "Chionochloa_rubra_var._rubra" & lat < -43.4, "NA", taxa)) %>%
  mutate(taxa = ifelse(taxa == "Chionochloa_rubra_subsp._cuprea" & lat > -41, "NA", taxa)) %>%
  mutate(taxa = ifelse(taxa == "Chionochloa_rubra", "NA", taxa)) %>%
  mutate(taxa = ifelse(taxa == "Chionochloa_spiralis" & lon > 167.8, "NA", taxa)) %>%
  mutate(taxa = ifelse(taxa == "Chionochloa_teretifolia" & lat < -47, "NA", taxa)) %>%
  mutate(taxa = ifelse(taxa == "Chionochloa_vireta" & lon > 172, "NA", taxa)) %>%
  filter(!grepl("NA", taxa)) %>%
  mutate(taxa = as.factor(taxa)) %>%
  arrange(taxa, lat, lon) %>%
  distinct(taxa, lat, lon)





###### combine for Festuca ######

# combine the herbaria, NVS, and Barbara's data frames
data.Festu1 <- do.call("rbind", list(data.herb.Festu, data.NVS.Festu, data.BARB.Festu))


# remove redundant values per taxa
# arrange the dataset alphabetially by taxa
data.Festu2 <- data.Festu1 %>%
  distinct(taxa, lon, lat) %>%
  arrange(taxa)






##########################################################################################
####################   CHECK THE DATA   ##################################################
##########################################################################################

# what does the formatted data look like?
head(data.Acaena2); str(data.Acaena2)
head(data.Chion3); str(data.Chion3)
head(data.Festu2); str(data.Festu2)



# what are all of our taxa?
levels(as.factor(data.Acaena2$taxa))
levels(as.factor(data.Chion3$taxa))
levels(as.factor(data.Festu2$taxa))



# what are all of our species?
levels(as.factor(data.Acaena2$species))
levels(as.factor(data.Chion3$species))
levels(as.factor(data.Festu2$species))




# looks good






##########################################################################################
####################   PLOT THE DATA   ###################################################
##########################################################################################

# load all the cool libraries
library("sp")
library("raster")
library("maps")
library("maptools")
library("mapdata")



##### make a  map for each genus #####



### create .pdf geographic map for each Acaena taxa
pdf("plot_Acaena_all-presence-data.pdf")

# set the first page to be all the Acaena occurances
# plot the NZ map
map(database = "nzHires", xlim = c( 165,180), ylim = c(-47.5,-33.6), fill = T, col = "grey", border = "black")

# plot all the Acaena points
points(data.Acaena2$lon, data.Acaena2$lat, cex = c(0.25), pch = 21, lwd = 0.40, col = "red")

# make the summary legend for first page
legend("topleft", c("data from NVS, NZVH, and K.L. obsv", paste("taxa = ", length(levels(unique(as.factor(data.Acaena2$taxa))))), paste( "n = ", length(data.Acaena2$taxa), sep = "")), bty = "n", pch = ("*"), col = 1, cex = 1)


# do a run for each taxa
for(i in 1:length(unique(data.Acaena2$taxa))){


# choose the taxa of interest through each run
todo <- which(data.Acaena2$taxa == unique(data.Acaena2$taxa)[i])
data.todo <- data.Acaena2[todo, ]

# plot all the Acaena points in balck for each run 
map(database = "nzHires", xlim = c( 165,180), ylim = c(-47.5,-33.6), fill = T, col = rgb(.7, 0.7, 0.7, 0.9), border = "black")
points(data.Acaena2$lon, data.Acaena2$lat, cex = 0.7, pch = 19, col = rgb(.2, 0.2, 0.2, 0.1))

# plot the points for only the taxa of interest in red for each run
points(data.todo$lon, data.todo$lat, cex = 0.7, pch = 20, col = "red")

# add in the legend for each run
legend("topleft", c(paste(data.todo$taxa[1]), paste( "n = ", length(data.todo$taxa), sep = "")), bty = "n", pch = ("*"), col = 1, cex = 1)

}

dev.off()






### create .pdf geographic map for each Chionochloa taxa
pdf("plot_Chionochloa_all-presence-data.pdf")

# plot the points of the entire data set for first page
# plot the map
map(database = "nzHires", xlim = c( 165,180), ylim = c(-47.5,-33.6), fill = T, col = "grey", border = "black")

# plot all the points
points(data.Chion3$lon, data.Chion3$lat, cex = c(0.25), pch = 21, lwd = 0.40, col = "red")

# make the summary legend
legend("topleft", c("data from NVS, NZVH, and K.L. obsv", paste("taxa = ", length(levels(unique(data.Chion3$taxa)))), paste( "n = ", length(data.Chion3$taxa), sep = "")), bty = "n", pch = ("*"), col = 1, cex = 1)


# do a run for each taxa
for(i in 1:length(unique(data.Chion3$taxa))){


# choose the taxa of interest through each run
todo <- which(data.Chion3$taxa == unique(data.Chion3$taxa)[i])
data.todo <- data.Chion3[todo, ]

# plot the points for all Chionochloa records in black
map(database = "nzHires", xlim = c( 165,180), ylim = c(-47.5,-33.6), fill = T, col = rgb(.7, 0.7, 0.7, 0.9), border = "black")
points(data.Chion3$lon, data.Chion3$lat, cex = 0.7, pch = 19, col = rgb(.2, 0.2, 0.2, 0.1))

# plot the points for only the taxa of interest in red
points(data.todo$lon, data.todo$lat, cex = 0.7, pch = 20, col = "red")

# add in the legend
legend("topleft", c(paste(data.todo$taxa[1]), paste( "n = ", length(data.todo$taxa), sep = "")), bty = "n", pch = ("*"), col = 1, cex = 1)

}

dev.off()






### create .pdf geographic map for each Festuca taxa
pdf("plot_Festuca_all-presence-data.pdf")

# set the first page to be all the Festuca occurances
# plot the NZ map
map(database = "nzHires", xlim = c( 165,180), ylim = c(-47.5,-33.6), fill = T, col = "grey", border = "black")

# plot all the Festuca points
points(data.Festu2$lon, data.Festu2$lat, cex = c(0.25), pch = 21, lwd = 0.40, col = "red")

# make the summary legend for first page
legend("topleft", c("data from NVS, NZVH, and K.L. obsv", paste("taxa = ", length(levels(unique(as.factor(data.Festu2$taxa))))), paste( "n = ", length(data.Festu2$taxa), sep = "")), bty = "n", pch = ("*"), col = 1, cex = 1)


# do a run for each taxa
for(i in 1:length(unique(data.Festu2$taxa))){


# choose the taxa of interest through each run
todo <- which(data.Festu2$taxa == unique(data.Festu2$taxa)[i])
data.todo <- data.Festu2[todo, ]

# plot the points for all Chionochloa records in black
map(database = "nzHires", xlim = c( 165,180), ylim = c(-47.5,-33.6), fill = T, col = rgb(.7, 0.7, 0.7, 0.9), border = "black")
points(data.Acaena2$lon, data.Acaena2$lat, cex = 0.7, pch = 19, col = rgb(.2, 0.2, 0.2, 0.1))

# plot the points for only the taxa of interest in red
points(data.todo$lon, data.todo$lat, cex = 0.7, pch = 20, col = "red")

# add in the legend
legend("topleft", c(paste(data.todo$taxa[1]), paste( "n = ", length(data.todo$taxa), sep = "")), bty = "n", pch = ("*"), col = 1, cex = 1)

}

dev.off()









##########################################################################################
####################   EXPORT THE DATA   #################################################
##########################################################################################


##### export a data file for each genus #####

# Acanea
write.csv(data.Acaena2, file = "rdata_Acaena_all-presence.csv", row.names = FALSE)



# Chionochloa
write.csv(data.Chion3, file = "rdata_Chionochloa_all-occurances.csv", row.names = FALSE)



# Festuca
write.csv(data.Festu2, file = "rdata_Festuca_all-presence.csv", row.names = FALSE)







##########################################################################################
####################   OTHER STUFF   #####################################################
##########################################################################################




