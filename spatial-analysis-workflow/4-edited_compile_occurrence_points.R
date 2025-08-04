################################################################################
# Load libraries
################################################################################

my.packages <- c('tidyverse','textclean','CoordinateCleaner','terra','countrycode'
                 # if you'd like to flag by elevation, you'll need:
                 ,'elevatr',"sf"
)
# versions I used (in the order listed above): 2.0.0, 0.9.3, 2.0-20, 1.7-29, 1.5.0, 0.4.5, 1.0-13
#install.packages (my.packages) #Turn on to install current versions
lapply(my.packages, require, character.only=TRUE)
rm(my.packages)

# be extra sure we are using dplyr when we mean to
select <- dplyr::select
filter <- dplyr::filter
mutate <- dplyr::mutate
rename <- dplyr::rename
count <- dplyr::count


################################################################################
# Set working directory
################################################################################

getwd()

# create folder for output data
data_out <- "taxon_points_raw"
if(!dir.exists(file.path(standard_occ,data_out)))
  dir.create(file.path(standard_occ,data_out), 
             recursive=T)

# assign folder where you have input data (saved in 3-get_occurrence_data.R)
data_in <- "input_datasets"

################################################################################
# Read in and compile occurrence point data
################################################################################

# read in raw occurrence data from 3-get_occurrence_data.R
file_list <- list.files(file.path(standard_occ,data_in), 
                        pattern = ".csv", full.names = T)
file_dfs <- lapply(file_list, read.csv, header = T, na.strings = c("","NA"),
                   colClasses = "character")
length(file_dfs)

# stack all datasets using bind_rows, which keeps non-matching columns
#   and fills with NA; 'Reduce' iterates through list and merges with previous;
#   this may take a few minutes if you have lots of data
all_data <- Reduce(bind_rows, file_dfs)
rm(file_dfs, file_list)
nrow(all_data)
table(all_data$database)

##### trims white space but doesn't seem to be the issue #####
# all_data$taxon_name <- trimws(all_data$taxon_name)
# taxon_list$taxon_name <- trimws(taxon_list$taxon_name)

################################################################################
# Filter by target taxa
################################################################################

# read in target taxa list
taxon_list <- read.csv(file.path(taxa_dir,"target_taxa_with_synonyms.csv"),
                       header=T, colClasses="character",na.strings=c("","NA"))
target_taxa <- unique(taxon_list$taxon_name)
# if needed, add columns that separate out taxon name
taxon_list <- taxon_list %>%
  separate("taxon_name",c ("genus","species","infra_rank","infra_name"),
           sep=" ", remove=F, fill="right")
# create species name column
taxon_list$species_name <- paste(taxon_list$genus,taxon_list$species)


## join data to taxon list
all_data <- all_data %>% dplyr::select(-genus,-specificEpithet)
all_data <- left_join(all_data, taxon_list, by = c("taxon_name"))


## join again just by species name if no taxon match
need_match <- all_data[which(is.na(all_data$taxon_name_status)),]
nrow(need_match)
# remove columns from first taxon name match
need_match <- need_match[,1:(ncol(all_data)-ncol(taxon_list)+1)]
# rename column for matching and get just the genus and specific epithet
need_match <- need_match %>% rename(taxon_name_full = taxon_name)
need_match$taxon_name <- str_split_i(need_match$taxon_name_full, " var.| f.",i = 1) 
  # removed | subsp. from the str split chain to keep in nesophila 
# new join
need_match <- left_join(need_match, taxon_list, by = c("taxon_name"))
need_match$taxon_name <- need_match$taxon_name_full
need_match <- need_match %>% select(-taxon_name_full)
# bind together new matches and previously matched
matched <- all_data[which(!is.na(all_data$taxon_name_status)),]
all_data <- bind_rows(matched,need_match)
table(all_data$taxon_name_status)

# setdiff(colnames(all_data), colnames(taxon_list))

# check names that got excluded.....
still_no_match <- all_data[which(is.na(all_data$taxon_name_status)),]
nrow(still_no_match)
sort(unique(still_no_match$taxon_name))
table(still_no_match$database)
# ^ if you want to keep any of these, add them to your target taxa list
rm(matched,need_match,still_no_match)

############ had to change all the left_join commands above to join just by taxon_name
# as still_no_match was coming up with taxa in my list 

all_data$taxon_name_accepted <- all_data$taxon_name_accepted.y
sort(unique(all_data$taxon_name_accepted))
# keep only rows for target taxa
all_data <- all_data %>% filter(!is.na(taxon_name_status))
nrow(all_data)

### see which target taxa have no occurrence data:
unique(taxon_list$taxon_name_accepted)[
  !(unique(taxon_list$taxon_name_accepted) %in% (unique(all_data$taxon_name_accepted)))]


################################################################################
# Add UID (universal id) column 
################################################################################

# now we add a UID that can be used to manually flag points in the maps created 
#   later in script 6-visualize_occurrence_points.R

# !!!
## IT IS VERY IMPORTANT to be careful here IF you are adding additional 
#   occurrence data after already compiling and manually vetting (selecting
#   good/bad points) - you need to be sure any new data are sorted to the *end*
#   of your dataset before adding the UIDs, otherwise the IDs you've already
#   used for vetting will change (not good!). If this is your first time running
#   this script, you're all good.
## if needed, sort data so any new points are at the end; you can use whatever
#   field you'd like, but here is an example using the database column:
all_data$database <- factor(all_data$database,
                            levels = c("BIEN","Ex_situ","FIA","GBIF","iDigBio",
                                       "IUCN_RedList","NorthAm_herbaria",
                                       "<NEW_DATA>")) #change to your dataset name
all_data <- all_data %>% arrange(database)

# add UID
nms <- names(all_data)
digits <- nchar(as.character(nrow(all_data)))
all_data <- all_data %>% 
  mutate(UID=paste0('id', sprintf(paste0("%0",digits,"d"),
                                  1:nrow(all_data)))) %>% 
  select(c('UID', all_of(nms)))


################################################################################
# Standardize/check some key columns
################################################################################

# create localityDescription column
#   "NAs introduced by coercion" warning ok
all_data <- all_data %>%
  unite("localityDescription",
        c(locality,municipality,county,stateProvince,country,
          locationNotes,verbatimLocality), remove = F, 
        sep = " | ") %>%
  mutate(decimalLatitude=as.numeric(decimalLatitude),
         decimalLongitude=as.numeric(decimalLongitude))
# get rid of NAs but keep pipes, so you can split back into parts if desired
all_data$localityDescription <- mgsub(all_data$localityDescription,
                                      c("NA "," NA"), "")
# if no locality info at all, make it NA
all_data$localityDescription <- gsub("| | | | | | | |", NA,
                                     all_data$localityDescription, fixed = T)
# check it
head(unique(all_data$localityDescription))

# check year column
#   if ex situ years were not standardized, most of those get removed here :(
#   "NAs introduced by coercion" warning ok
all_data$year <- as.numeric(all_data$year)
# remove values less than 1500 or greater than current year
all_data$year[which(all_data$year < 1500 |
                      all_data$year > as.numeric(format(Sys.time(),"%Y")))] <- NA
sort(unique(all_data$year))

# check basis of record column
unique(all_data$basisOfRecord)
all_data$basisOfRecord[which(is.na(all_data$basisOfRecord))] <- "UNKNOWN"

# check establishment means
unique(all_data$establishmentMeans)
# use this example if you need to recode any mistakes:
#all_data <- all_data %>%
#  mutate(establishmentMeans = recode(establishmentMeans,
#    "Introduced" = "INTRODUCED",
#    "Uncertain" = "UNCERTAIN",
#    "Native" = "NATIVE"))
all_data$establishmentMeans[which(is.na(all_data$establishmentMeans))] <-
  "UNCERTAIN"

## check validity of lat and long; remove invalid coordinates
# if coords are both 0, set to NA
zero <- which(all_data$decimalLatitude == 0 & all_data$decimalLongitude == 0)
all_data$decimalLatitude[zero] <- NA; all_data$decimalLongitude[zero] <- NA; rm(zero)
# flag non-numeric and not available coordinates, and when not valid, i.e.:
#   lat > 90, lat < -90, lon > 180, and lon < -180
coord_test <- cc_val(all_data, lon = "decimalLongitude", lat = "decimalLatitude",
                     value = "flagged", verbose = TRUE)
# try switching lat and long for invalid points and check validity again
all_data[!coord_test,c("decimalLatitude","decimalLongitude")] <-
  all_data[!coord_test,c("decimalLongitude","decimalLatitude")]
coord_test <- cc_val(all_data, lon = "decimalLongitude",lat = "decimalLatitude",
                     value = "flagged", verbose = TRUE)
# mark these as flagged
all_data$flag <- NA
all_data[!coord_test,]$flag <- "Coordinates invalid"
rm(coord_test)

################################################################################
# Add lat-long country code column & flag points that are not on land
################################################################################

# we clip the occurrence points to the countries layer since some analyses rely
#   on knowing which country (or sometimes state) the point falls within
# IMPROVEMENT OPPORTUNITY: a buffer could be added to the countries layer, to 
#   catch points that are just outside the feature boundary but still could
#   be on land, since the features aren't so exact (would make them too large)

# read in world countries layer created in 1-prep_gis_layers.R
world_ctry <- vect(file.path(gis_dir,"world_countries_10m",
                             "world_countries_10m.shp"))
# select just the 2-digit country code column we need
world_ctry <- world_ctry[,"iso_a2"]

# add country polygon data to each point based on lat-long location
# turn occurrence point data into a spatial object
have_coord <- all_data %>% filter(!is.na(decimalLatitude) & !is.na(decimalLongitude))
pts_spatial <- vect(cbind(have_coord$decimalLongitude,have_coord$decimalLatitude),
                    atts=have_coord, crs="+proj=longlat +datum=WGS84")
no_coord <- anti_join(all_data,have_coord)
# intersect with countries layer
# if this takes too long, you may need to download the 50m layer in script 
#   1-prep_gis_layers.R instead of the 10m layer
pts_spatial <- terra::intersect(pts_spatial,world_ctry)
# if the following two lines are not zero, you should look at script
#  1-prep_gis_layers.R and fill in an additional missing code(s)
nrow(pts_spatial[which(is.na(pts_spatial$iso_a2)),])
nrow(pts_spatial[which(pts_spatial$iso_a2=="-99"),])
# bring all the data back together and flag water points
on_land <- as.data.frame(pts_spatial)
on_land <- on_land %>% rename(latlong_countryCode = iso_a2)
land_id <- unique(on_land$UID)
in_water <- have_coord %>% filter(!(UID %in% land_id)); nrow(in_water)
in_water$flag <- "Not on land"
have_coord <- full_join(on_land,in_water)
all_data2 <- full_join(have_coord,no_coord)
rm(have_coord,pts_spatial,no_coord,on_land,land_id,in_water)
nrow(all_data2)


################################################################################
# Select records that have valid coordinates and are not in the water
################################################################################
# edited section to include Araucaria luxurians as it grows on shore
################################################################################

# first, if you're working at the taxon level, add infrataxon records to their 
#   parent species too
add_again <- all_data2 %>% filter(grepl("var\\.|subsp\\.", taxon_name_accepted))
unique(add_again$taxon_name_accepted)
add_again$taxon_name_accepted <- gsub(" var\\.*\\s.+", "",add_again$taxon_name_accepted)
unique(add_again$taxon_name_accepted)
all_data2 <- rbind(all_data2,add_again)
table(all_data2$database)

# write file of raw data before selecting only geolocated records;
#   this will be used for the GapAnalysis package's summary of occurrences
# if you don't plan to run the GapAnalysis package, you can skip this step
#   since it saves a pretty large file if lots of points for your taxa!
write.csv(all_data2, file.path(standard_occ,
                               paste0("all_occurrence_data_raw_", Sys.Date(), 
                                      ".csv")), row.names = F)

# separate out points with locality description but no valid lat-long
 table(all_data2$flag)
 locality_pts <- all_data2 %>% filter(!is.na(localityDescription) & 
                                      flag == "Coordinates invalid")
 nrow(locality_pts)
#  save the file, for reference; you can geolocate these records if you want, 
#   but usually only necessary for rare taxa without enough lat-long records
write.csv(locality_pts, file.path(standard_occ,
                               paste0("could_attempt_geolocation_", 
                                         Sys.Date(), ".csv")), row.names = F)

############################################################################
# have # this out as some points are along rivers - manually check instead 
# e.g if coordinates from RBGE then probably reliable etc
# and e.g if there's lots of points probably look into keeping unless wrong etc
############################################################################
# separate out points that will be removed because they're in the water; can
#   reference this, and potentially move points to the land if desired
# water_pts <- all_data2 %>% filter(flag == "Not on land" & 
                                   # taxon_name != "Araucaria luxurians" &
                                   # taxon_name != "Araucaria scopulorum")
# nrow(water_pts)
# save the file, for reference; you can geolocate these records if you want, 
#   but usually only necessary for rare taxa without enough lat-long records
# write.csv(water_pts, file.path(standard_occ,
                        #     paste0("removed_water_points_", Sys.Date(), ".csv")),
#          row.names = F)

### EDITED VERSION create final subset that is only points with valid lat-long on land
 geo_pts <- all_data2 %>% filter(is.na(flag) | 
  (flag == "Not on land" & taxon_name == "Araucaria luxurians"))
 nrow(geo_pts)
# see how many points are from each database
 table(geo_pts$database)

### create final subset that is only points with valid lat-long on land
# geo_pts <- all_data2 %>% filter(is.na(flag))
# nrow(geo_pts)
# see how many points are from each database
# table(geo_pts$database)


################################################################################
# Standardize country code column for checking against lat-long later
################################################################################

# country name to 3 letter ISO code
# fix some issues first (can add anything that is not matched unambiguously in
# the 'find codes for names' step directly below; then rerun from here)
geo_pts$country <- mgsub(geo_pts$country,
                         c("AU","AUS","Australia - Norfolk Island","Bresil",
                           "Chili","CHL","NC","New South Wales","BR","CL","NZL","Australie"
                           ,"NZ","Hawaii","North America"),#matched unambiguously
                         c("Australia","Australia","Australia","Brazil",
                           "Chile","Chile","New Caledonia","Australia","Brazil","Chile",
                           "New Zealand","Australia","New Zealand","USA","USA")) #country name to use
# find codes for names
country_codes1 <- as.data.frame(sort(unique(geo_pts$country))) %>%
  add_column(iso3c = countrycode(sort(unique(geo_pts$country)),
                                 origin="country.name", destination="iso2c")) %>%
  rename("country" = "sort(unique(geo_pts$country))",
         "countryCode_standard" = "iso3c")
# add to data
geo_pts <- left_join(geo_pts,country_codes1)

# country code to standard 2 letter ISO code
# fix some issues (can add anything that is not matched unambiguously in the
# 'check codes' step directly below; then rerun from here)
geo_pts$countryCode <- mgsub(geo_pts$countryCode,
                             c("AUT","CAN","CZE","DEU","MEX","NOR","SWE","USA","ZZ"), #matched unambiguously
                             c("AU","CA","CZ","DE","MX","NO","SE","US",NA)) #2-letter code to use
geo_pts$countryCode <- str_to_upper(geo_pts$countryCode)
# check codes 
country_codes2 <- as.data.frame(sort(unique(geo_pts$countryCode))) %>%
  add_column(iso3c = countrycode(sort(unique(geo_pts$countryCode)),
                                 origin="iso2c", destination="iso2c")) %>%
  rename("countryCode" = "sort(unique(geo_pts$countryCode))",
         "countryCode_standard2" = "iso3c")
# add to data
geo_pts <- left_join(geo_pts,country_codes2)
geo_pts <- geo_pts %>% 
  unite("countryCode_standard",
        c("countryCode_standard","countryCode_standard2"),
        sep=";",remove=T,na.rm=T) %>%
  separate("countryCode_standard","countryCode_standard",sep=";",extra="drop")
geo_pts$countryCode_standard[which(geo_pts$countryCode_standard == "")] <- NA
# see how many records have each country code:
sort(table(geo_pts$countryCode_standard))


################################################################################
# (Optionally) Add elevation data to points
################################################################################
# install.packages("elevatr")

# create an sf object of your coordinates
sf_points <- st_as_sf(geo_pts, 
                      coords = c("decimalLongitude","decimalLatitude"),
                      crs = 4326)
# add elevation column to the points
#   this can take a couple minutes or longer if you have many points; if it 
#   takes too long, you can move this step to script #5, inside the 
#   taxon-by-taxon loop
add_elev <- get_elev_point(locations = sf_points,
                           proj = "+proj=longlat +datum=WGS84",
                           src = "aws")
# add elevation column to our main data frame
geo_pts$elevationInMeters <- add_elev$elevation


################################################################################
# Remove spatial duplicates
################################################################################

## NOTE THAT IF YOUR TARGET TAXA ARE NOT WIDESPREAD, IT MAY BE BETTER TO SKIP
## THIS STEP, ESPECIALLY IF YOU'RE USING RESTRICTING FILTERS -- E.G. ONLY USING
## POINTS WITH A SMALL COORDINATE UNCERTAINTY -- SINCE YOU MAY NEED ALL THE 
## DATA YOU HAVE! This step was designed for taxa with many coordinates.


################################################################################
# Select final set of columns
################################################################################

## set header/column name order and keep only necessary columns
keep_col <- c( 
  #universal ID and source database
  "UID","database",
  #taxon
  "taxon_name","scientificName","taxon_name_accepted",
  "taxonRank","infraspecificEpithet","taxonIdentificationNotes",
  #event
  "year","basisOfRecord",
  #record-level
  "nativeDatabaseID","datasetName","publisher",
  "rightsHolder","references",
  "issue","recordedBy",
  #occurrence
  "establishmentMeans","individualCount",
  #location
  "decimalLatitude","decimalLongitude","coordinateUncertaintyInMeters",
  "elevationInMeters", # if you added elevation column
  "geolocationNotes","localityDescription","locality","verbatimLocality",
  "locationNotes","municipality","county","stateProvince",
  "country","countryCode_standard","latlong_countryCode",
  #additional optional taxa metadata
  "rl_category","all_native_dist_iso2","match_info"
)
geo_pts <- geo_pts[,keep_col]

colnames(geo_pts)


# take a look
head(as.data.frame(geo_pts))
nrow(geo_pts)
table(geo_pts$database)

# write file if you'd like; not necessary since we write taxon-level files
write.csv(geo_pts, file.path(standard_occ,
 paste0("all_occurrence_points_compiled_", Sys.Date(), ".csv")), row.names = F)



################################################################################
# Create taxon-level summary table
################################################################################

# summarize results for each target taxon
# count lat-long records
count_geo <- geo_pts %>% count(taxon_name_accepted)
names(count_geo)[2] <- "num_latlong_records"
# count records with invalid lat-long but have locality description
count_locality <- locality_pts %>% count(taxon_name_accepted)
names(count_locality)[2] <- "num_locality_records"
# count records that fell outside country feature boundaries
# count_water<- water_pts %>% count(taxon_name_accepted)
# names(count_water)[2] <- "num_water_records"
# make table
pieces <- list(count_geo,count_locality,
               data.frame(taxon_name_accepted=unique
                          (taxon_list$taxon_name_accepted)))
summary <- Reduce(full_join, pieces)
summary <- summary %>% 
  mutate_all(~replace(., is.na(.), 0)) %>%
  arrange(!is.na(num_latlong_records), num_latlong_records)
head(as.data.frame(summary))
# write file
write.csv(summary, file.path(standard_occ,
                             paste0("summary_of_occurrences_", 
                                    Sys.Date(), ".csv")),row.names = F)



################################################################################
# Split by taxa to save
################################################################################

# split records to create one CSV for each target taxon
sp_split <- split(geo_pts, as.factor(geo_pts$taxon_name_accepted))
names(sp_split) <- gsub(" ","_",names(sp_split))
names(sp_split) <- gsub("\\.","",names(sp_split))

# write files
lapply(seq_along(sp_split), function(i) write.csv(sp_split[[i]],
                                                  file.path(standard_occ,data_out,
                                                            paste0(names
                                                                   (sp_split)[[i]],
                                                                   ".csv")),
                                                  row.names = F))
