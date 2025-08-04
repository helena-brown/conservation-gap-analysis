################################################################################
# Load libraries
################################################################################

my.packages <- c('tidyverse','textclean','CoordinateCleaner','tools','terra')
# versions I used (in the order listed above): 2.0.0, 0.9.3, 2.0-20, 4.3.0, 1.7-29
#install.packages (my.packages) #Turn on to install current versions
lapply(my.packages, require, character.only=TRUE)
rm(my.packages)

################################################################################
# Set working directory
################################################################################

getwd()

# create folder for output data
data_out <- "taxon_points_ready-to-vet"
if(!dir.exists(file.path(standard_occ,data_out)))
  dir.create(file.path(standard_occ,data_out), 
             recursive=T)

# assign folder where you have input data (saved in 4-compile_occurrence_points.R)
data_in <- "taxon_points_raw"


################################################################################
# Read in data
################################################################################

# read in target taxa list
taxon_list <- read.csv(file.path(taxa_dir,"target_taxa_with_synonyms.csv"),
                       header=T, colClasses="character",na.strings=c("","NA"))

# read in world countries layer created in 1-prep_gis_layers.R
world_polygons <- vect(file.path(gis_dir,"world_countries_10m",
                                 "world_countries_10m.shp"))

# read in urban areas layer created in 1-prep_gis_layers.R
urban_areas <- vect(file.path(gis_dir,"urban_areas_50m",
                              "urban_areas_50m.shp"))

################################################################################
# Iterate through taxon files and flag potentially suspect points
################################################################################

## set threshold for coordinate uncertainty flag; anything greater than this 
##  value (meters) will be flagged
max_uncertainty <- 1000

# list of taxon files to iterate through
taxon_files <- list.files(path=file.path(standard_occ,data_in), 
                          ignore.case=FALSE, full.names=FALSE, recursive=TRUE)
target_taxa <- file_path_sans_ext(taxon_files)

# start a table to add summary of results for each species
summary_tbl <- data.frame(
  taxon_name_accepted = "start", 
  unflagged_pts = "start", 
  selected_pts = "start", 
  .cen = "start", 
  .urb = "start",
  .inst = "start",
  .con = "start", 
  .outl = "start", 
  .nativectry = "start", 
#  .yr1950 = "start", 
 # .yr1980 = "start",
 # .yrna = "start",
 # .unc = "start",
  .elev = "start",
  .rec = "start",
  stringsAsFactors=F)


## iterate through each species file to flag suspect points...

for (i in 1:length(target_taxa)){
  
  taxon_file <- target_taxa[i]
  taxon_nm <- mgsub(taxon_file, c("_","_var_","_subsp_"), 
                    c(" "," var. "," subsp. "))
  
  # bring in records
  taxon_now <- read.csv(file.path(standard_occ,data_in,
                                  paste0(taxon_file, ".csv")))
  # print the taxon name we're working with
  cat("-----\n","Starting ", taxon_nm, ", taxon ", i, " of ", length(target_taxa), ".\n", sep="")
  cat("Number of records: ",nrow(taxon_now),"\n",sep="")
  
  
  # now we will go through a set of tests to flag potentially suspect records...
  
  
  ### FLAG RECORDS THAT HAVE COORDINATES NEAR COUNTRY AND STATE/PROVINCE CENTROIDS
  flag_cen <- CoordinateCleaner::cc_cen(
    taxon_now,
    lon = "decimalLongitude", lat = "decimalLatitude",
    # buffer = radius around country/state centroids (meters); default=1000
    buffer = 1000, value = "flagged")
  taxon_now$.cen <- flag_cen
  
  
  ## FLAG RECORDS THAT HAVE COORDINATES IN URBAN AREAS
   if(nrow(taxon_now)<2){
    taxon_now$.urb <- NA
    print("Taxa with fewer than 2 records will not be tested.")
  } else {
    taxon_now <- as.data.frame(taxon_now)
    flag_urb <- CoordinateCleaner::cc_urb(
      taxon_now,
      lon = "decimalLongitude",lat = "decimalLatitude",
      ref = urban_areas, value = "flagged")
    taxon_now$.urb <- flag_urb
  }
  
  
  ### FLAG RECORDS THAT HAVE COORDINATES NEAR BIODIVERSITY INSTITUTIONS
  flag_inst <- CoordinateCleaner::cc_inst(
    taxon_now,
    lon = "decimalLongitude", lat = "decimalLatitude",
    # buffer = radius around biodiversity institutions (meters); default=100
    buffer = 100, value = "flagged")
  taxon_now$.inst <- flag_inst
  
  
  ### COMPARE THE COUNTRY LISTED IN THE RECORD VS. THE LAT-LONG COUNTRY; flag
  #   when there is a mismatch; CoordinateCleaner package has something like 
  #   this but also flags when the record doesn't have a country...
  #   you could use that function if you want to flag NAs, or edit below:
  taxon_now <- taxon_now %>% mutate(.con=(ifelse(
    (as.character(latlong_countryCode) == as.character(countryCode_standard) &
       !is.na(latlong_countryCode) & !is.na(countryCode_standard)) |
      is.na(latlong_countryCode) | is.na(countryCode_standard), TRUE, FALSE)))
  cat("Testing country listed\n",sep="")
  cat("Flagged ",length(which(!taxon_now$.con))," records.\n",sep="")
  
  
  ### FLAG SPATIAL OUTLIERS
  taxon_now <- as.data.frame(taxon_now)
  flag_outl <- CoordinateCleaner::cc_outl(
    taxon_now,
    lon = "decimalLongitude",lat = "decimalLatitude",
    species = "taxon_name_accepted", 
    # read more about options for the method and the multiplier:
    #   https://www.rdocumentation.org/packages/CoordinateCleaner/versions/2.0-20/topics/cc_outl
    # if you make the multiplier larger, it will flag less points.
    # the default is 5; you may need to experiment a little to see what works
    #   best for most of your target taxa (script #6 helps you view flagged pts)
    method = "quantile", mltpl = 7, 
    value = "flagged")
  taxon_now$.outl <- flag_outl
  
  
  ### CHECK LAT-LONG COUNTRY AGAINST "ACCEPTED" NATIVE COUNTRY DISTRUBUTION; 
  #   flag when the lat-long country is not in the list of native countries;
  #   we use the native countries compiled in 1-get_taxa_metadata.R, which
  #   combines the IUCN Red List, BGCI GlobalTreeSearch, and manually-added data
  native_ctrys <- unique(unlist(strsplit(as.character(taxon_now$all_native_dist_iso2), "; ")))
  if(!is.na(native_ctrys[1])){
    # flag records where native country doesn't match record's coordinate location
    taxon_now <- taxon_now %>% 
      mutate(.nativectry=(ifelse(latlong_countryCode %in% native_ctrys, 
                                 TRUE, FALSE)))
  } else {
    # if no country data for the taxon, leave unflagged
    taxon_now$.nativectry <- TRUE
  }
  cat("Testing native countries\n",sep="")
  cat("Flagged ",length(which(!taxon_now$.nativectry))," records.\n",sep="")
  
  
  ### FLAG OLDER RECORDS, based on two different year cutoffs (1950 & 1980)
  # 1950 cutoff
 # taxon_now <- taxon_now %>% mutate(.yr1950=(ifelse(
  #  (as.numeric(year)>1950 | is.na(year)), TRUE, FALSE)))
#  cat("Testing year < 1950\n",sep="")
 # cat("Flagged ",length(which(!taxon_now$.yr1950))," records.\n",sep="")
  # 1980 cutoff
#  taxon_now <- taxon_now %>% mutate(.yr1980=(ifelse(
 #   (as.numeric(year)>1980 | is.na(year)), TRUE, FALSE)))
 # cat("Testing year < 1980\n",sep="")
 # cat("Flagged ",length(which(!taxon_now$.yr1980))," records.\n",sep="")
  
  
  ### FLAG RECORDS THAT DON'T HAVE A YEAR PROVIDED
 # taxon_now <- taxon_now %>% mutate(.yrna=(ifelse(
  #  !is.na(year), TRUE, FALSE)))
 # cat("Testing year NA\n",sep="")
 # cat("Flagged ",length(which(!taxon_now$.yrna))," records.\n",sep="")
  
  
  ### FLAG RECORDS WITH COORDINATE UNCERTAINTY ABOVE YOUR THRESHOLD (FLAGS NA!) 
 # taxon_now <- taxon_now %>% mutate(.unc=(ifelse(
   # as.numeric(coordinateUncertaintyInMeters)<max_uncertainty & 
   #   !is.na(coordinateUncertaintyInMeters), TRUE, FALSE)))
#  cat("Testing coordinate uncertainty\n",sep="")
#  cat("Flagged ",length(which(!taxon_now$.unc))," records.\n",sep="")
  
  
  ### FLAG RECORDS OUTSIDE YOUR ELEVATION RANGE, BY TAXON
  ###   for this you need an 'elevation_range' column in your 
  ###   target_taxa_with_synonyms.csv file
  if("elevation_range" %in% colnames(taxon_now)){
    elev_range_raw <- str_squish(as.character(taxon_now$elevation_range))[1] 
    # added as.character to avoid errors
    if(!is.na(elev_range_raw)){
      elev_range <- unique(unlist(strsplit(elev_range_raw, "-")))
      elev_min <- as.numeric(elev_range[1])
      elev_max <- as.numeric(elev_range[2])
      taxon_now <- taxon_now %>% mutate(.elev=(ifelse(
        (as.numeric(elevationInMeters) >= elev_min &
           as.numeric(elevationInMeters) <= elev_max) |
          is.na(elevationInMeters), TRUE, FALSE)))
    } else {
      taxon_now <- taxon_now %>% mutate(.elev=TRUE)
    }
  } else {
    taxon_now <- taxon_now %>% mutate(.elev=TRUE)  # <<--- Ensure .elev is always added
  }
    cat("Testing elevation\n",sep="")
    cat("Flagged ",length(which(!taxon_now$.elev))," records.\n",sep="")
  
  
  
  ### FLAG RECORDS THAT ARE NOT CURRENT/NATIVE BASED ON TWO COLUMNS:
  ###   basisOfRecord AND/OR establishmentMeans
  taxon_now <- taxon_now %>% mutate(.rec=(ifelse(
    basisOfRecord != "FOSSIL_SPECIMEN" & 
      basisOfRecord != "LIVING_SPECIMEN" &
      establishmentMeans != "INTRODUCED" & 
      establishmentMeans != "MANAGED" &
      establishmentMeans != "CULTIVATED", 
    TRUE, FALSE)))
  cat("Testing basisOfRecord & establishmentMeans\n",sep="")
  cat("Flagged ",length(which(!taxon_now$.rec))," records.\n\n",sep="")
  
  
  
  # create some subsets to count how many records are in each, for summary table...
  
  # count of completely unflagged points
  total_unflagged <- taxon_now %>%
    filter(.cen & 
              .urb &
              .inst & .con & .outl & .nativectry 
           # & .yr1950 & .yr1980  & .yrna 
           &
             # .unc 
              .elev & .rec)
  
  # OPTIONAL count of unflagged points based on selected filters you'd like 
  #   to use, to see how many points there are; change as desired; commented- 
  #   out lines are the filters you don't want to use
  select_unflagged <- taxon_now %>%
    filter(
      .cen & 
        .inst & 
        .outl &
        .con & 
        .urb & 
        .nativectry &
      #  .yr1950 & 
      #  .yr1980 & 
      #  .yrna &
      #  .unc &
        .elev &
        .rec
    )
  
  # add data to summary table
  summary_add <- data.frame(
    taxon_name_accepted = taxon_nm,
    unflagged_pts = nrow(total_unflagged),
    selected_pts = nrow(select_unflagged),
    .cen = sum(!taxon_now$.cen),
    .urb = sum(!taxon_now$.urb),
    .inst = sum(!taxon_now$.inst),
    .con = sum(!taxon_now$.con),
    .outl = sum(!taxon_now$.outl),
    .nativectry = sum(!taxon_now$.nativectry),
  #  .yr1950 = sum(!taxon_now$.yr1950),
  #  .yr1980 = sum(!taxon_now$.yr1980),
  #  .yrna = sum(!taxon_now$.yrna),
  #  .unc = sum(!taxon_now$.unc),
    .elev = sum(!taxon_now$.elev),
    .rec = sum(!taxon_now$.rec),
    stringsAsFactors=F)
  summary_tbl[i,] <- summary_add
  
  # WRITE NEW FILE
  write.csv(taxon_now, file.path(standard_occ,data_out,
                                 paste0(taxon_file, "gbif.csv")), row.names=FALSE)
  
}


write.csv(summary_tbl, file.path(standard_occ,
                               paste0("justgbif_Summary_Table_5Flag_", Sys.Date(),
                                      ".csv")), row.names=FALSE)



# add summary of points to summary we created in 4-compile_occurrence_points.R
file_nm <- list.files(path = file.path(standard_occ),
                      pattern = "justgbif_summary_of_occurrences_2025-07-11.csv", 
                      full.names = T)
orig_summary <- read.csv(file_nm, colClasses = "character")
# keep just the first four columns, in case you're running this script a second time
orig_summary <- orig_summary %>% select(taxon_name_accepted)
summary_tbl2 <- full_join(orig_summary,summary_tbl,by="taxon_name_accepted")
summary_tbl2

# write summary table
write.csv(summary_tbl2, file.path(standard_occ,
                                  paste0("justgbif_summary_of_occurrences_", Sys.Date(),
                                         ".csv")),row.names = F)

