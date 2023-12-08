# install required packages and define universal parameters
if (!requireNamespace("rinat", quietly = TRUE)) install.packages("rinat")
require(rinat)
if (!requireNamespace("taxize", quietly = TRUE)) install.packages("taxize")
require(taxize)
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
require(dplyr)
if (!requireNamespace("plyr", quietly = TRUE)) install.packages("plyr")
require(plyr)
if (!requireNamespace("lubridate", quietly = TRUE)) install.packages("lubridate")
require(lubridate)
if (!requireNamespace("tidyr", quietly = TRUE)) install.packages("tidyr")
require(tidyr)
if (!requireNamespace("sf", quietly = TRUE)) install.packages("sf")
require(sf)
if (!requireNamespace("arcgisbinding", quietly = TRUE)) install.packages("arcgisbinding")
require(arcgisbinding)
if (!requireNamespace("rgbif", quietly = TRUE)) install.packages("rgbif")
require(rgbif)
if (!requireNamespace("purrr", quietly = TRUE)) install.packages("purrr")
require(purrr)
if (!requireNamespace("data.table", quietly = TRUE)) install.packages("data.table")
require(data.table)

# load the arcgis license
arc.check_product()

# set output geodatabase
output_gdb <- "H:/Projects/RKM_Invasives/Tiers/Data/Data.gdb"

#get current year and number of years of data
years <- 2000:year(now())

################################################################################
#### prepare PA species list with ITIS TNS IDs and iNat IDs - 
#this doesn't have to be done each time if we already have a prepared species list with all necessary IDs. if so, skip to load step
################################################################################

#paths to the NY and PA ultimate species lists
ny_species <- "C:/Users/mmoore/OneDrive - Western Pennsylvania Conservancy/Projects/RKM Invasives/NY_ultimate_species_list_2023.csv"
pa_species <- "C:/Users/mmoore/OneDrive - Western Pennsylvania Conservancy/Projects/RKM Invasives/PA_ultimate_species_list_2023.csv"

#read in species lists
ny_species_df <- read.csv(ny_species)
pa_species_df <- read.csv(pa_species)

#use taxize package to get ITIS_TSNs for matching SNAMEs
pa_species_df$ITIS_TSN <- get_tsn(pa_species_df$PA_sciName,"scientific")

#join itis TSNs and iNat IDs from NY to PA species list - join by iMap Network Species ID
pa_species_join <- pa_species_df %>%
                    left_join(select(ny_species_df,iMap_Network_SpeciesID, itis_tsn, iNat.Taxon_ID_NS), by=c("Network_Species_List_Id"="iMap_Network_SpeciesID"), na_matches="never") 


#cleanup - change "NULL" characters to NA and change itis id columns to integer
pa_species_join[pa_species_join == "NULL"] <- NA
pa_species_join$itis_tsn <- as.integer(as.character(pa_species_join$itis_tsn))
pa_species_join$ITIS_TSN <- as.integer(as.character(pa_species_join$ITIS_TSN))


pa_species_join <- pa_species_join %>%
  mutate(ITIS_TSN = coalesce(ITIS_TSN,itis_tsn))

pa_species_join$gbif_key <- get_gbifid(sci=pa_species_join$PA_sciName)

write.csv(pa_species_join,"C:/Users/mmoore/OneDrive - Western Pennsylvania Conservancy/Projects/RKM Invasives/PA_ultimate_species_list_2023_ids.csv")

################################################################################
#### load species list if already prepared
################################################################################

pa_species_list <- "C:/Users/mmoore/OneDrive - Western Pennsylvania Conservancy/Projects/RKM Invasives/PA_ultimate_species_list_2023_ids.csv"
pa_species_join <- read.csv(pa_species_list)

################################################################################
##### DOWNLOAD DATA FROM INATURALIST
################################################################################

species_list <- as.vector(pa_species_join$iNat.Taxon_ID_NS)

a <- list()
k <- NULL
for(x in 1:length(species_list)){
  #get metadata on the number of occurrences
  print(paste("getting metadata from iNaturalist for ",species_list[x],".", sep="") )
  try(k <- get_inat_obs(taxon_id=species_list[x], bounds=c(38.2963065, -82.5470831, 43.7707612,-72.1215712), geo=TRUE, meta=TRUE, quality="research")) # this step first queries iNat to see if there are any records present, if there are it actually downloads them.
  Sys.sleep(10) # this is too throttle our requests so we don't overload their servers
  if(is.list(k)){
    print(paste("There are ", k$meta$found, " records on iNaturalist", sep=""))
    if(k$meta$found>0 && k$meta$found<=10000){
      a[[x]] <- get_inat_obs(taxon_id=species_list[x], bounds=c(38.2963065, -82.5470831, 43.7707612,-72.1215712), geo=TRUE, quality="research", maxresults = k$meta$found)
      k <- NULL
    }
    else if(k$meta$found>10000){
      p <- NULL
      for(y in years){
        #get metadata on the number of occurrences
        print(paste("getting metadata from iNaturalist for year: ",y,".", sep="") )
        try(p <- get_inat_obs(taxon_id=species_list[x], bounds=c(38.2963065, -82.5470831, 43.7707612,-72.1215712), geo=TRUE, meta=TRUE, quality="research", year=y) ) # this step first queries iNat to see if there are any records present, if there are it actually downloads them.
        Sys.sleep(10) # this is too throttle our requests so we don't overload their servers
        if(is.list(p)){
          print(paste("There are ", p$meta$found, " records on iNaturalist", sep=""))
          if(p$meta$found>0 && p$meta$found<=10000){
            a[[x]] <- get_inat_obs(taxon_id=species_list[x], bounds=c(38.2963065, -82.5470831, 43.7707612,-72.1215712), geo=TRUE, quality="research", year=y, maxresults = p$meta$found)
            k <- NULL
          }}}
    } else {}
  } else {
    print("No records found")
  }
}

# convert to a data frame
inatrecs <- ldply(a)

###### DO SOME STUFF TO DROP THOSE WITHOUT CREATIVE COMMONS LICENSING AND ANYTHING ELSE?
# drop those without creative comments license and those that are captive/cultivated and those with uncertainty >300m
inatrecs[inatrecs==""]<-NA
inat_df <- inatrecs %>% 
  drop_na(license) %>%
  filter(captive_cultivated!='true') %>%
  filter(is.na(public_positional_accuracy)|public_positional_accuracy<=300)

inat_sf <- st_as_sf(inat_df, coords=c("longitude","latitude"), crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
inat_sf <- st_transform(inat_sf, crs=st_crs("ESRI:102008")) # reproject to the north america albers equal area
arc.write(path=paste(output_gdb,"inat_export",sep="/"), inat_sf, overwrite=TRUE) # write a feature class into the geodatabase

################################################################################
######## DOWNLOAD DATA FROM GBIF
################################################################################

keys <- as.list(pa_species_join$gbif_key)
keys <- keys[!is.na(keys)]

datalist = vector("list", length = 0)

# break up 100mi buffer on PA into six sections - otherwise, too many records to download
# bounding_geoms <- c('POLYGON ((-76.055890042999977 41.027789386000052, -75.949174252999967 43.809236436000049, -79.587646130999985 43.876029018000054, -79.571056897999938 41.091774731000044, -76.055890042999977 41.027789386000052))',
#                     'POLYGON ((-76.152113308999958 38.285799703000066, -76.055890042999977 41.027789386000052, -79.571056897999938 41.091774731000044, -79.556126897999945 38.352101811000068, -76.152113308999958 38.285799703000066))',
#                     'POLYGON ((-79.571056897999938 41.091774731000044, -79.587646130999985 43.876029018000054, -82.651499611999952 43.862992074000033, -82.55647294399995 41.077566560000037, -79.571056897999938 41.091774731000044))',
#                     'POLYGON ((-79.556126897999945 38.352101811000068, -79.571056897999938 41.091774731000044, -82.55647294399995 41.077566560000037, -82.471075827999982 38.340795451000076, -79.556126897999945 38.352101811000068))',
#                     'POLYGON ((-76.055890042999977 41.027789386000052, -72.621582254999964 40.880540860000053, -72.262611202999949 43.649947729000075, -75.949174252999967 43.809236436000049, -76.055890042999977 41.027789386000052))',
#                     'POLYGON ((-76.152113308999958 38.285799703000066, -72.945985547999953 38.144619584000054, -72.621582254999964 40.880540860000053, -76.055890042999977 41.027789386000052, -76.152113308999958 38.285799703000066))')
# 
# for(bounds in bounding_geoms){
#   print(bounds)
  dat <- occ_data(
    taxonKey=keys, 
    limit=10000, # modify if needed, fewer will make testing go faster 
    hasCoordinate=TRUE,
    geometry='POLYGON ((-80.577111647999971 42.018959019000079, -80.583025511999949 39.690462536000041, -77.681987232999973 39.68735201800007, -75.761816590999956 39.690666106000037, -75.678308913999956 39.790810226000076, -75.53064649099997 39.815101786000071, -75.411566911999955 39.776679135000052, -75.101245089999964 39.880029385000057, -75.09383042199994 39.944216030000064, -74.690932882999959 40.133570156000076, -74.690425973999936 40.17528313400004, -74.893196517999968 40.350896889000069, -74.914505704999954 40.415842984000051, -75.012247039999977 40.448477402000037, -75.004556583999943 40.522413349000033, -75.134560399999941 40.623471625000036, -75.136516799999981 40.723392383000032, -75.002409694999983 40.867515299000047, -75.082051382999964 40.971575944000051, -74.830463730999952 41.152763058000062, -74.768212647999974 41.271891205000031, -74.640518995999969 41.358839422000074, -74.709416559999966 41.454495330000043, -74.826329023999961 41.475865789000068, -74.936988959999951 41.521739840000066, -75.018029425999941 41.617276498000081, -75.012709979999954 41.733926517000043, -75.061642930999938 41.85481505100006, -75.218658916999971 41.904656042000056, -75.336705265999967 42.017618624000079, -77.511689405999959 42.017704281000078, -79.721693517999938 42.024739989000068, -79.715980736999938 42.353623043000027, -80.577111647999971 42.018959019000079))', # simplified boundary of Pennsylvania. # simplified boundary of Pennsylvania.
    year='2000,2023'
  )
  
  dat <- dat[dat!="no data found, try a different search"]
  dd <- rbindlist(lapply(dat, function(x) x$data), fill = TRUE, use.names = TRUE)
  print(nrow(dd))
  datalist[[bounds]] <- dd
}

big_data <- do.call(rbind, datalist, fill=TRUE)

#exclude records with >300m uncertainty
gbif_df <- dd[!(dd$coordinateUncertaintyInMeters >300),]
#remove records with null lat/long
gbif_df <- gbif_df %>%
  drop_na(decimalLatitude) %>%
  drop_na(decimalLongitude)

# create a spatial layer
gbif_sf <- st_as_sf(gbif_df, coords=c("decimalLongitude","decimalLatitude"), crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
gbif_sf <- st_transform(gbif_sf, crs=st_crs("ESRI:102008")) # reproject to the north america albers equal area
arc.write(path=paste(output_gdb,"gbif_export",sep="/"), gbif_sf, overwrite=TRUE) # write a feature class into the geodatabase

################################################################################
######## DOWNLOAD DATA FROM NAS
################################################################################














#gbif geometry for PA
geometry='POLYGON ((-80.577111647999971 42.018959019000079, -80.583025511999949 39.690462536000041, -77.681987232999973 39.68735201800007, -75.761816590999956 39.690666106000037, -75.678308913999956 39.790810226000076, -75.53064649099997 39.815101786000071, -75.411566911999955 39.776679135000052, -75.101245089999964 39.880029385000057, -75.09383042199994 39.944216030000064, -74.690932882999959 40.133570156000076, -74.690425973999936 40.17528313400004, -74.893196517999968 40.350896889000069, -74.914505704999954 40.415842984000051, -75.012247039999977 40.448477402000037, -75.004556583999943 40.522413349000033, -75.134560399999941 40.623471625000036, -75.136516799999981 40.723392383000032, -75.002409694999983 40.867515299000047, -75.082051382999964 40.971575944000051, -74.830463730999952 41.152763058000062, -74.768212647999974 41.271891205000031, -74.640518995999969 41.358839422000074, -74.709416559999966 41.454495330000043, -74.826329023999961 41.475865789000068, -74.936988959999951 41.521739840000066, -75.018029425999941 41.617276498000081, -75.012709979999954 41.733926517000043, -75.061642930999938 41.85481505100006, -75.218658916999971 41.904656042000056, -75.336705265999967 42.017618624000079, -77.511689405999959 42.017704281000078, -79.721693517999938 42.024739989000068, -79.715980736999938 42.353623043000027, -80.577111647999971 42.018959019000079))', # simplified boundary of Pennsylvania.

#inat geometry for PA
bounds=c(39.7198, -80.519891, 42.26986,	-74.689516)