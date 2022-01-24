#### R Script AL544 zooplankton distribution plot
# Author: H.Hauss (hhauss@geomar.de)
# Plot integrated zooplankton/Copepoda biomass from Ecotaxa tsv export
#AL544 was a cruise in September 2020. Mesozooplankton was collected using a 
#WP-2 plankton net (200µm mesh size), vertical hauls from a couple of meters over
#bottom to surface.
#samples were scanned using an Epson V750 pro scanner and imported into Ecotaxa

#load libraries (install first)
library(ggplot2)
library(maps)
library(mapdata)
library(scatterpie)

## set working directory -- change this to where you saved your .tsv data
setwd("V:/Daten/Cruises/AL544")  

## list all ecotaxa in a folder
file_list <- list.files(pattern="*.tsv") # create list of all .tsv files in that folder
print(file_list)
##read in all ecotaxa tsvs into one dataframe and select variables to keep
data_raw <- do.call(rbind,lapply(file_list,read.csv, header=TRUE, quote = '"', sep = "\t"))[ ,c('sample_id','object_annotation_category','object_annotation_hierarchy', 'object_date', 'object_time','object_lat', 'object_lon','object_depth_min', 'object_depth_max', 'object_area', 'sample_volconc')]

## remove detritus and other non-living organisms from dataframe
data <-data_raw[-grep("not-living", data_raw$object_annotation_hierarchy),] 
rm(data_raw)

## check format of variables (numeric, character, factor, etc.)
str(data)

## rename some variables for easier/shorter names
names(data) = gsub(pattern = "object_", replacement = "", x = names(data))
##date format
data$date <- as.Date(as.character(data$date), format =  "%Y%m%d")
##time format
data$time <- sprintf("%06d", data$time)
data$time <- format(strptime(data$time, format="%H%M%S"), format = "%H:%M")
##calculate image area from pixels (pixel size = 10.6µm)
data$area_mm2 <- data$area*0.00011236


##convert area to biomass for different taxa -- these relationships are from Lehette & Hernandez-leon and
##may be not too accurate for Baltic zooplankton
data$biomass_ug <- with(data, 
                      ifelse(grepl("Calanoida", annotation_hierarchy),45.25*data$area_mm2^1.59,
                          ifelse(grepl("Annelida", annotation_hierarchy),43.38*data $area_mm2^1.54,
                              ifelse(grepl("Chaetognatha", annotation_hierarchy),23.45*data $area_mm2^1.19,
                                  ifelse(grepl("Gastropoda", annotation_hierarchy),43.38*data $area_mm2^1.54,
                                      44.78*data$area_mm2^1.56)))))

##build summary dataframe with total abundance per volume
summary <- aggregate(annotation_category ~ (sample_id+date+lat+lon+depth_min+depth_max+sample_volconc), data, length)
names(summary)[names(summary) == 'annotation_category'] <- 'total_count'
summary$total_abundance_m3 <- summary$total_count/summary$sample_volconc
#integrate (total abundance under a square meter of the sea)
summary$total_abundance_m2 <- summary$total_abundance_m3*summary$depth_max

##same for biomass
total_bm <- aggregate(biomass_ug ~ (sample_id+sample_volconc), data, sum)
total_bm$biomass_um_m3 <- total_bm$biomass/total_bm$sample_volconc
total_bm = subset(total_bm, select = -c(biomass_ug,sample_volconc))
###aggregate abundance for copepods only
## filter copepods only
copepods <- dplyr::filter(data, grepl('Copepoda', annotation_hierarchy))
##plot a histogram of the weight frequency distribution
hist(copepods$biomass_ug, main="Copepod Biomass Distribution",
     xlab="Biomass in µg DW",
     xlim=c(0,20),breaks=500)
##aggregate copepod biomass to samples
total_cop_bm <- aggregate(biomass_ug ~ (sample_id+sample_volconc), copepods, sum)
total_cop_bm$copepod_biomass_um_m3 <- total_cop_bm$biomass_ug/total_cop_bm$sample_volconc
total_cop_bm = subset(total_cop_bm, select = -c(biomass_ug,sample_volconc) )

## filter cladocerans only
cladocera <- dplyr::filter(data, grepl('Cladocera', annotation_hierarchy))
hist(cladocera$biomass_ug, main="Cladocera Biomass Distribution",
     xlab="Biomass in µg DW",
     xlim=c(0,20),breaks=100)
total_clado_bm <- aggregate(biomass_ug ~ (sample_id+sample_volconc), cladocera, sum)
total_clado_bm$cladocera_biomass_um_m3 <- total_clado_bm$biomass_ug/total_clado_bm$sample_volconc
total_clado_bm = subset(total_clado_bm, select = -c(biomass_ug,sample_volconc) )

##join copepods and cladocerans on summary dataframe
# merge the biomass data frames into summary dataframe by sample ID
summary <- merge(summary,total_bm, by="sample_id")
summary <- merge(summary,total_cop_bm, by="sample_id")
summary <- merge(summary,total_clado_bm, by="sample_id")

#calculate integrated biomass in mg per m2
summary$total_biomass_mg_m2 <- summary$biomass_um_m3*summary$depth_max/1000
summary$copepod_biomass_mg_m2 <- summary$copepod_biomass_um_m3*summary$depth_max/1000
summary$cladocera_biomass_mg_m2 <- summary$cladocera_biomass_um_m3*summary$depth_max/1000

##calculate fraction from total for copepods, cladocerans and others
summary$copepod_fraction = summary$copepod_biomass_mg_m2/summary$total_biomass_mg_m2
summary$cladocera_fraction = summary$cladocera_biomass_mg_m2/summary$total_biomass_mg_m2
summary$other_fraction = 1 - summary$copepod_fraction - summary$cladocera_fraction
##export summary data in case this shall be published (e.g. on Pangaea)
write.table(summary, file = "AL544_biomass.txt")

################plot on map
###load coastline data for Europe, using the map_data function
Area <- map_data("world2Hires", region =c("Germany", "Denmark", "Poland", "Sweden"))


#plot map of total biomass as bubble plot
ggplot() +
  theme_bw()+
  geom_polygon(fill="darkgoldenrod1", color = "black", data = Area, aes(x=long, y = lat, group = group)) + ##this adds the coastlines
  geom_point(data = summary, aes(x = lon, y = lat, size = total_biomass_mg_m2), color = "red", alpha=0.7) +
  coord_quickmap(xlim=c(10,15), ylim=c(53.5,55.5))  + ##map boundaries
  labs(title = "AL544",x = "Longitude °W", y = "Latitude °N", size = "mg / m²")                                         ##axis labels

#plot map with piecharts to show fraction of different groups
ggplot() +
  theme_bw()+
  geom_polygon(fill="darkgoldenrod1", color = "black", data = Area, aes(x=long, y = lat, group = group)) +  ## coastlines                                                       ##map boundaries
  labs(title = "AL544",x = "Longitude °W", y = "Latitude °N", size = "Biomass (mg / m²)")  +                          ##axis labels
  geom_point(data = summary, aes(x = lon, y = lat, size = total_biomass_mg_m2), pch=21, fill = "transparent", color= "black") +
  geom_scatterpie(aes(x=lon, y=lat, group = sample_id, r = ((total_biomass_mg_m2/3.14)^0.5)/100), 
                  data = summary, cols = colnames(summary[,c(17:19)])) +
  scale_fill_manual(name="Group",
                    labels=c("Copepoda","Cladocera","other"), values=c("green1","coral1","darkslateblue")) +
  coord_fixed(xlim=c(10,15), ylim=c(53.5,55.5))


ggsave("V:/Daten/Cruises/AL544/biomass_pie.pdf", width = 7, height = 5, bg = "transparent")