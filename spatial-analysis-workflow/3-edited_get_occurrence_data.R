################################################################################
# Load libraries
################################################################################

# install.packages("ridigbio")
# install.packages('BIEN')
my.packages <- c('tidyverse','textclean','data.table','rgbif','ridigbio','BIEN')
# versions I used (in the order listed above): 2.0.0, 0.9.3, 1.14.8, 3.7.7, 0.3.6, 1.2.6
# install.packages(my.packages) #Turn on to install current versions
lapply(my.packages, require, character.only=TRUE)
rm(my.packages)

################################################################################
# Set working directory
################################################################################

getwd()

# create folder for output data
if(!dir.exists(file.path(standard_occ,"input_datasets")))
  dir.create(file.path(standard_occ,"input_datasets"), 
             recursive=T)
data_out <- "input_datasets"

################################################################################
# Load or create target taxa list
################################################################################

# read in taxa list
taxon_list <- read.csv(file.path(taxa_dir,"target_taxa_with_synonyms.csv"),
                       header=T, colClasses="character",na.strings=c("","NA"))
head(taxon_list); nrow(taxon_list)

######################## only need extra GBIF data for taxa with no IUCN data 
##################### (excluding Wollemia as not mapping anyway)
species_of_interest <- c(
  "Agathis australis",
  "Agathis borneensis",
  "Agathis dammara",
  "Agathis flavescens",
  "Agathis kinabaluensis",
  "Agathis lanceolata",
  "Agathis lenticula",
  "Agathis macrophylla",
  "Agathis montana",
  "Agathis moorei",
  "Agathis orbicula",
  "Agathis ovata",
  "Agathis robusta subsp. nesophila",
  "Araucaria angustifolia",
  "Araucaria araucana",
  "Araucaria goroensis",
  "Araucaria heterophylla",
  "Araucaria humboldtensis",
  "Araucaria luxurians",
  "Araucaria montana",
  "Araucaria muelleri",
  "Araucaria nemorosa",
  "Araucaria rulei",
  "Araucaria schmidii",
  "Araucaria scopulorum",
  "Wollemia nobilis"
)

taxon_list <- taxon_list[taxon_list$taxon_name_accepted %in% species_of_interest, ]

# list of target taxon names
taxon_names <- sort(taxon_list$taxon_name)



# create list of species names as well, to do initial removal of non-target 
#   taxa when downloading at the genus level; full name match in next script
# target_sp_names <- unique(sapply(taxon_list$taxon_name, function(x) {
#  if (grepl("^Agathis robusta", x)) {
 #   x  # keep full name
 # } else {
 #   unlist(strsplit(x, " var. | subsp. | f. "))[1]  # extract base name
 # }
# }))
# view(target_sp_names)

### you can also create a target taxa list by hand instead...
# include synonyms if you want to find them:
#taxon_names <- c("Asimina incana","Asimina longifolia","Asimina triloba")
# if you have any infrataxa, just keep the genus and specific epithet here:
#target_sp_names <- c("Asimina incana","Asimina longifolia","Asimina triloba")

################################################################################
# Download & do basic standardization of occurrence data from each database
################################################################################

## You can pick and choose any/all sections below (A-F), depending on
#  which databases you'd like to use
## Some sections have two options for downloading data: manually via
#  the website, or automatically using the API; choose whichever works for you
## Note that if you have a long taxa list, the downloaded data may be very
#  large; if this makes it unworkable for your computer setup, you may need to 
#  skip writing the files and move directly to 4-compile_occurrence_data.R to
#  compile and filter (edits will be needed in this script and script 4 to make 
#  that work)

###
###############
###############################################
### A) Global Biodiversity Information Facility (GBIF)
###    https://www.gbif.org
###############################################
###############
###

# create new folder if not already present
if(!dir.exists(file.path(raw_occ,"GBIF")))
  dir.create(file.path(raw_occ,"GBIF"), recursive=T)

###
### OPTION 1: automatic download via API
### (can go down to option 2 -manual download- if this isn't working)
###

# load GBIF account user information
# if you don't have account yet, go to https://www.gbif.org then click
# "Login" in top right corner, then click "Register"
# either read in a text file with username, password, and email (one on each
#   line) or manually fill in below (if you're not saving this script publicly):
# login <- read_lines(log_loc)
user  <- "hbrown15" #username
pwd   <- "Gate.GBIF2003" #password
email <- "brownhelena15@gmail.com" #email
# rm(login)

# get GBIF taxon keys for all taxa in target list
keys <- sapply(taxon_names,function(x) name_backbone(name=x)$speciesKey,
               simplify = "array"); keys
# remove duplicate and NULL keys
keys_nodup <- keys[
 # !duplicated(keys) & 
    keys != "NULL"]
# create data frame of keys and matching taxon_name
gbif_codes <- map_df(keys_nodup,~as.data.frame(.x),.id="taxon_name")
names(gbif_codes)[2] <- "speciesKey"
# create vector of keys as input into gbif download
gbif_taxon_keys <- vector(mode="numeric")
for(i in 1:length(keys_nodup)){
  gbif_taxon_keys <- c(gbif_taxon_keys,keys_nodup[[i]][1])
}; sort(gbif_taxon_keys)
rm(i)
# download GBIF data (Darwin Core Archive format)

gbif_download <- occ_download(
  pred_in("taxonKey", gbif_taxon_keys),
  ## you can also add additional filters:
  pred_in("basisOfRecord", c("PRESERVED_SPECIMEN",
  #    "HUMAN_OBSERVATION","OBSERVATION",
  #    "UNKNOWN","MACHINE_OBSERVATION",
  "MATERIAL_SAMPLE")),
  #    "LITERATURE")),
  pred("hasCoordinate", TRUE),
  pred("hasGeospatialIssue", FALSE),
  ## "SIMPLE_CSV" format is another option
  format = "DWCA",
  user=user,pwd=pwd,
  email=email)
rm(user, pwd, email)
# load gbif data just downloaded
# download and unzip before reading in
download_key <- gbif_download
# must wait for download to complete before continuing;
# it may take a while (up to 3 hours) if you have a large taxa list;
# function below will pause script until the download is ready, though it
#   seems like its not working quite right for me right now?
# you can also log in on the GBIF website and go to your profile to see the
#   progress of your download 
occ_download_wait(download_key, status_ping=10, quiet=TRUE)
# get download when its ready then unzip and read in
occ_download_get(key=download_key[1],
                 path=file.path(raw_occ,"GBIF"))
unzip(zipfile=paste0(
  file.path(raw_occ,"GBIF",download_key[1]),".zip"),
  files="occurrence.txt",
  exdir=file.path(raw_occ,"GBIF"))
# create file list for next step (standardizing)
file_list <- list.files(path = file.path(raw_occ,"GBIF"), 
                        pattern = "occurrence", full.names = T)
length(file_list) #1

# !!! NOW SCROLL DOWN TO "STANDARDIZE THE DATA" SECTION !!!


###
### STANDARDIZE THE DATA
### Run this section no matter which option your chose for getting the data !!
###

# loop through file(s) [you may have more than one if you downloaded manually]
 for(i in 1:length(file_list)){
  # read in data
  gbif_raw <- fread(file_list, quote="", na.strings="")
  print(paste0("Total number of records: ",nrow(gbif_raw)))
  # remove any genus-level records
  gbif_raw <- gbif_raw %>% filter(taxonRank != "GENUS")
  # keep only necessary columns
  gbif_raw <- gbif_raw %>% select(
    # Taxon
    "scientificName","family","genus","specificEpithet","taxonRank",
    "infraspecificEpithet",
    # *concatenate into taxonIdentificationNotes
    "identificationRemarks","identificationVerificationStatus","identifiedBy",
    "taxonRemarks",
    # Event
    "year",
    # Record-level
    "basisOfRecord","gbifID","datasetName","publisher","institutionCode",
    "rightsHolder","license","references","informationWithheld","issue",
    # Occurrence
    "recordedBy","establishmentMeans",
    # Location
    "decimalLatitude","decimalLongitude","coordinateUncertaintyInMeters",
    "county","municipality","stateProvince","higherGeography","countryCode",
    # *concatenate into golocationNotes
    "georeferencedBy","georeferencedDate","georeferenceProtocol",
    "georeferenceRemarks","georeferenceSources",
    "georeferenceVerificationStatus",
    # *concatenate into locationNotes
    "locality","verbatimLocality","associatedTaxa","eventRemarks","fieldNotes",
    "habitat","locationRemarks","occurrenceRemarks","occurrenceStatus"
  ) %>% rename(nativeDatabaseID = gbifID)
  # add database column
  gbif_raw$database <- "GBIF"
  # combine a few similar columns
  gbif_raw <- gbif_raw %>% unite("taxonIdentificationNotes",
                                 identificationRemarks:taxonRemarks,na.rm=T,remove=T,sep=" | ")
  gbif_raw$taxonIdentificationNotes <-
    gsub("^$",NA,gbif_raw$taxonIdentificationNotes)
  gbif_raw <- gbif_raw %>% unite("locationNotes",
                                 associatedTaxa:occurrenceStatus,na.rm=T,remove=T,sep=" | ")
  gbif_raw$locationNotes <- gsub("^$",NA,gbif_raw$locationNotes)
  gbif_raw <- gbif_raw %>% unite("geolocationNotes",
                                 georeferencedBy:georeferenceVerificationStatus,na.rm=T,remove=T,sep=" | ")
  gbif_raw$geolocationNotes <- gsub("^$",NA,gbif_raw$geolocationNotes)
  # create taxon_name column
  gbif_raw$taxon_name <- NA
 # subsp <- gbif_raw %>% filter(taxonRank == "SUBSPECIES")
 # subsp$taxon_name <- paste(subsp$genus,subsp$specificEpithet,"subsp.",
      #                      subsp$infraspecificEpithet)
#  var <- gbif_raw %>% filter(taxonRank == "VARIETY")
#  var$taxon_name <- paste(var$genus,var$specificEpithet,"var.",
     #                     var$infraspecificEpithet)
#  form <- gbif_raw %>% filter(taxonRank == "FORM")
#  form$taxon_name <- paste(form$genus,form$specificEpithet,"f.",
     #                      form$infraspecificEpithet)
#  spp <- gbif_raw %>% filter(taxonRank == "SPECIES")
#  spp$taxon_name <- paste(spp$genus,spp$specificEpithet)
  
  # For each taxon rank, filter and build taxon_name only if rows exist
  subsp <- gbif_raw %>% filter(taxonRank == "SUBSPECIES")
  if(nrow(subsp) > 0){
    subsp$taxon_name <- paste(subsp$genus, subsp$specificEpithet, "subsp.",
                              subsp$infraspecificEpithet)
  }
  
  var <- gbif_raw %>% filter(taxonRank == "VARIETY")
  if(nrow(var) > 0){
    var$taxon_name <- paste(var$genus, var$specificEpithet, "var.",
                            var$infraspecificEpithet)
  }
  
  form <- gbif_raw %>% filter(taxonRank == "FORM")
  if(nrow(form) > 0){
    form$taxon_name <- paste(form$genus, form$specificEpithet, "f.",
                             form$infraspecificEpithet)
  }
  
  spp <- gbif_raw %>% filter(taxonRank == "SPECIES")
  if(nrow(spp) > 0){
    spp$taxon_name <- paste(spp$genus, spp$specificEpithet)
  }
  
  
  gbif_raw <- Reduce(rbind,list(subsp,var,form,spp))
  rm(subsp,var,form,spp)
  # fix hybrid names as needed (I don't think any were found)
  print(sort(unique(gbif_raw$taxon_name)))
 # gbif_raw$taxon_name <- mgsub(gbif_raw$taxon_name,
                   #            c("Asimina xnashii"),
                   #            c("Asimina x nashii"))
  # create species_name column
  gbif_raw$species_name <- NA
  gbif_raw$species_name <- sapply(gbif_raw$taxon_name, function(x)
    unlist(strsplit(x," var. | subsp. | f. "))[1])
  # keep only target species
  gbif_raw2 <- gbif_raw %>%
    filter(species_name %in% target_sp_names)
  print(paste0("Number of records for target taxa: ",nrow(gbif_raw)))
  print("Species removed; if you want to keep any of these, add them to your target taxa list...")
  print(unique(setdiff(gbif_raw,gbif_raw2)[,"taxon_name"]))
  # check a few standards and recode if needed
  # establishmentMeans
  # check and add below as needed
  print(sort(unique(gbif_raw2$establishmentMeans))) 
  gbif_raw2 <- gbif_raw2 %>%
    mutate(establishmentMeans = recode(establishmentMeans,
                                       "Uncertain" = "UNCERTAIN",
                                       "Native" = "NATIVE",
                                       "Introduced" = "INTRODUCED"))
  # write file
  write.csv(gbif_raw2, file.path(standard_occ,data_out,
                                 paste0("just_gbif",i,".csv")),row.names=FALSE)
  rm(gbif_raw); rm(gbif_raw2)
  
 }

###
###############
###############################################
# C) IUCN Red List of Threatened Species
#    https://www.iucnredlist.org
###############################################
###############
###

# create new folder if not already present
if(!dir.exists(file.path(raw_occ,"IUCN_RedList")))
  dir.create(file.path(raw_occ,"IUCN_RedList"), recursive=T)

## First, download raw data:
# Go to https://www.iucnredlist.org/search
# Click "Login/Register" in top bar; create an account if you don't have one,
#   then log in to your account
# Open the "Taxonomy" tab in the left bar, type in your target genus name
#   and check the box next to the genus name when it comes up - note that
#   sometimes there is a glitch where the box is grayed out; you can navigate
#   to the genus you want by going through the tree, starting with Plantae. 
#   Remember to search for any synonym genera as well, as applicable.
# Or, alternatively, if you are just looking for a few
#   taxa you can search for them individually.
# You should be able to add each genus/taxon to your search so only
#   one file needs to be exported.
# In the far-left bar, scroll down and, if desired, check
#   "Subspecies and varieties" (if you want infrataxa).
# Click the gray "Download" button and select "Range data - Points (CSV)"
#   then fill in the prompts in the popup window.
# Go to https://www.iucnredlist.org/account to find your query
# Click "Download" next to your query
# Move the folder you downloaded into the "IUCN_RedList" folder in
#   occurrence_data > raw_occurrence_data
# Open the folder you just added and pull the points_data.csv file out into 
#   the IUCN_RedList folder


# read in data
redlist_raw <- read.csv(file.path(raw_occ,
                                  "IUCN_RedList","points_data.csv"), colClasses = "character",
                        na.strings=c("", "NA"), strip.white=T, fileEncoding="UTF-8")
nrow(redlist_raw)
# create taxon_name column
spp <- redlist_raw %>% filter(is.na(subspecies))
spp$taxon_name <- spp$sci_name
subsp <- redlist_raw %>% filter(!is.na(subspecies))
if(nrow(subsp!=0)){
  subsp$taxon_name <- paste(subsp$sci_name,"subsp.",subsp$subspecies) }
redlist_raw <- rbind(subsp,spp)
sort(unique(redlist_raw$taxon_name))
# create species_name column
redlist_raw$species_name <- NA
# redlist_raw$species_name <- sapply(redlist_raw$taxon_name, function(x)
  # unlist(strsplit(x," subsp. "))[1])
redlist_raw$species_name <- redlist_raw$taxon_name
# keep only target species
redlist_raw <- redlist_raw %>%
  filter(species_name %in% target_sp_names)
# combine a few similar columns
redlist_raw$rightsHolder <- paste(redlist_raw$citation,redlist_raw$year)
# keep only necessary columns
redlist_raw <- redlist_raw %>%
  select(taxon_name,tax_comm,event_year,
         basisofrec,origin,dec_lat,dec_long,dist_comm,source,
         compiler,rightsHolder,presence,subspecies,id_no)
# check a few standards and recode if needed
# establishmentMeans
redlist_raw <- redlist_raw %>%
  mutate(origin = recode(origin,
                         "1" = "NATIVE",
                         "2" = "REINTRODUCED",
                         "3" = "INTRODUCED",    "4" = "VAGRANT",
                         "5" = "UNCERTAIN",
                         "6" = "ASSISTED_COLONISATION"))
# issue
redlist_raw <- redlist_raw %>%
  mutate(presence = recode(presence,
                           "1" = "EXTANT",
                           "2" = "EXTANT",
                           "3" = "POSSIBLY_EXTANT",
                           "4" = "POSSIBLY_EXTINCT",
                           "5" = "EXTINCT",
                           "6" = "PRESENCE_UNCERTAIN")) %>%
  # remove extinct rows
  filter(presence != "EXTINCT")
# basisOfRecord
redlist_raw$basisofrec <- str_to_lower(redlist_raw$basisofrec)
redlist_raw$basisofrec <- gsub(" |_","",redlist_raw$basisofrec)
unique(redlist_raw$basisofrec) # check and add below as needed
redlist_raw <- redlist_raw %>%
  mutate(basisofrec = recode(basisofrec,
                             "humanobservation" = "HUMAN_OBSERVATION",
                             "preservedspecimen" = "PRESERVED_SPECIMEN",
                             "literature" = "LITERATURE",
                             "expert" = "HUMAN_OBSERVATION",
                             "fossilspecimen" = "FOSSIL_SPECIMEN",
                             "livingspecimen" = "LIVING_SPECIMEN",
                             "unknown" = "UNKNOWN",
                             "liturature" = "LITERATURE",
                             "materialsample" = "MATERIAL_SAMPLE",
                             "observation" = "OBSERVATION",
                             "specimen" = "PRESERVED_SPECIMEN",
                             .missing = "UNKNOWN"))
# rename to fit standard
redlist_raw <- redlist_raw %>%
  rename(taxonIdentificationNotes = tax_comm,
         year = event_year,
         basisOfRecord = basisofrec,
         decimalLatitude = dec_lat,
         decimalLongitude = dec_long,
         locality = dist_comm,
         datasetName = source,
         recordedBy = compiler,
         issue = presence,
         establishmentMeans = origin,
         taxonRank = subspecies,
         nativeDatabaseID = id_no,
         references = compiler)
# add database column & publisher columns
redlist_raw$database <- "IUCN_RedList"
redlist_raw$publisher <- "IUCN Red List of Threatened Species"
# write file
write.csv(redlist_raw, file.path(standard_occ,data_out,
                                 "redlist.csv"), row.names=FALSE)
rm(redlist_raw)


###
###############
###############################################
# G) Ex situ accession-level data: wild collection locations
# These data are compiled in 2-compile_exsitu_data.R and here we are simply 
# formatting them for use in the occurrence point dataset
###############################################
###############
###

# read in ex situ data we saved in 2-compile_exsitu_data.R (edit to match your 
#   file names)
exsitu_raw1 <- read.csv(file.path(raw_occ,"Ex_situ",
                                  "ExSitu_Compiled_Post-Geolocation_2025-06-19.csv"), colClasses = "character",
                        na.strings=c("", "NA"), strip.white=T, fileEncoding="UTF-8")
exsitu_raw2 <- read.csv(file.path(raw_occ,"Ex_situ",
                                  "ExSitu_Dead_2025-06-16.csv"), colClasses = "character",
                        na.strings=c("", "NA"), strip.white=T, fileEncoding="UTF-8")
exsitu_raw <- bind_rows(exsitu_raw1,exsitu_raw2)
nrow(exsitu_raw)
rm(exsitu_raw1,exsitu_raw2)
# rename columns to fit standard
exsitu_raw <- exsitu_raw %>%
  rename("taxon_name" = "taxon_name",
         "scientificName" = "taxon_full_name_orig",
         "specificEpithet" = "species",
         "taxonRank" = "infra_rank",
         "infraspecificEpithet" = "infra_name",
         "taxonIdentificationNotes" = "taxon_verif",
         "year" = "coll_year",
         "nativeDatabaseID" = "UID",
         "publisher" = "data_source",
         "rightsHolder" = "inst_short",
         "references" = "acc_num",
         "issue" = "latlong_flag",
         "individualCount" = "num_indiv",
         "decimalLatitude" = "lat_dd",
         "decimalLongitude" = "long_dd",
         "coordinateUncertaintyInMeters" = "latlong_uncertainty",
         "verbatimLocality" = "all_locality",
         "stateProvince" = "state",
         "establishmentMeans" = "prov_type") %>%
  unite("recordedBy", c("coll_name","coll_num"), 
        remove=T, sep="; ", na.rm=T) %>%
  unite("geolocationNotes", c("latlong_det","geolocated_by","latlong_notes"), 
        remove=T, sep="; ", na.rm=T) %>%
  unite("locationNotes", c("germ_type"), 
        remove=T, sep="; ", na.rm=T)
# establishmentMeans recode
sort(unique(exsitu_raw$establishmentMeans))
exsitu_raw <- exsitu_raw %>%
  mutate(establishmentMeans = recode(establishmentMeans,
                                     "H" = "CULTIVATED",
                                     "H?" = "CULTIVATED",
                                     "N" = "NATIVE",
                                     "NG" = "UNCERTAIN",
                                     "U" = "UNCERTAIN",
                                     "W" = "NATIVE",
                                     "Z" = "NATIVE"
  ))
# add a few informational columns
exsitu_raw$database <- "Ex_situ"
exsitu_raw$basisOfRecord <- "HUMAN_OBSERVATION"
exsitu_raw$datasetName <- exsitu_raw$rightsHolder
# keep only necessary columns
exsitu_raw <- exsitu_raw %>%
  select(database,taxon_name,taxon_name_accepted,scientificName,genus,specificEpithet,taxonRank,
         infraspecificEpithet,taxonIdentificationNotes,year,basisOfRecord,
         nativeDatabaseID,datasetName,publisher,rightsHolder,references,
         issue,recordedBy,establishmentMeans,individualCount,decimalLatitude,
         decimalLongitude,coordinateUncertaintyInMeters,geolocationNotes,
         locality,verbatimLocality,locationNotes,municipality,county,
         stateProvince,country,match_info)
head(exsitu_raw)
# write file
write.csv(exsitu_raw, file.path(standard_occ,data_out,
                                "exsitu.csv"), row.names=FALSE)
rm(exsitu_raw)

