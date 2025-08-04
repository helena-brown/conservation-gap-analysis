# load packages
my.packages <- c('tidyverse','textclean','CoordinateCleaner','terra')
# versions I used (in the order listed above): 2.0.0, 0.9.3, 2.0-20, 1.7-29
# install.packages(my.packages) #Turn on to install current versions
lapply(my.packages, require, character.only=TRUE)
rm(my.packages)

# be sure we're using dplyr when we want to
select <- dplyr::select
rename <- dplyr::rename
filter <- dplyr::filter
group_by <- dplyr::group_by
mutate <- dplyr::mutate
distinct <- dplyr::distinct



################################################################################
# 1. Read in and stack all ex situ accessions data from survey of collections
################################################################################

# read in data from multiple surveys and stack, or just read in from one folder.
# this function also adds columns for 1) the file name [often equivalent to the
# "inst_short" institution nickname] 2) a submission year, 3) an accession
# number if one isn't given
### CHANGE BASED ON FOLDER(S) AND YEAR(S) YOU HAVE...
all_data <- read.csv(file.path(exsitu_dir, 
                                      "input_cleaning_data.csv"), header = TRUE)




# to check the column headings 
sort(colnames(all_data))


################################################################################
# 1. Read in and stack all ex situ accessions data from survey of collections
################################################################################


# remove leading, trailing, and middle (e.g., double space) white space
all_data <- as.data.frame(lapply(all_data, function(x) str_squish(x)), 
                          stringsAsFactors=F)
# replace "" cells with NA in whole dataset
all_data[all_data == ""] <- NA
# add a general data source column
all_data$data_source <- "ex_situ_BG_survey"
# fix some common lat/long character issues before replacing non-ascii characters
all_data$orig_lat <- mgsub(all_data$orig_lat,
                           c("\"","ç","d","'","°","º","Â","¼","،","¡","_","´",
                             "*","À","?","`")," ")
all_data$orig_long <- mgsub(all_data$orig_long,
                            c("\"","ç","d","'","°","º","Â","¼","،","¡","_","´",
                              "*","À","?","`")," ")
# now replace all non-ascii characters
all_data <- as.data.frame(lapply(all_data,replace_non_ascii), stringsAsFactors=F)

# add additional necessary columns, some of which are in the optional data we 
#   add from public databases like Genesys
all_data$latlong_flag <- NA
all_data$latlong_det <- NA
all_data$latlong_uncertainty <- NA
all_data$geolocated_by <- NA
all_data$latlong_notes <- NA



################################################################################
# 4. Save raw output of compiled data for target genera
#    (can be used to look for hybrids/cultivars, which are removed in next step)
################################################################################

# create new version of data before big edits, so we can more easily go back; 
# we do this throughout but I'll just mention it here; this does take more 
# computational power, so if you have a lot of data you may want to skip this 
# and edit the script to remove
all_data2 <- all_data

# read in target taxa list
taxon_list <- read.csv(file.path(taxa_dir,"target_taxa_with_synonyms.csv"),
                       header=T, colClasses="character",na.strings=c("","NA"))
str(taxon_list)

# preserve original taxon name before we edit
all_data2$taxon_full_name_orig <- all_data2$taxon_full_name

### CHECK OUT THE HYBRID COLUMN...
sort(unique(all_data2$hybrid))
## NOTE that this is only really necessary if you have a hybrid in your target
#   taxa list. If you don't, all you need to do is make sure there is nothing
#   in the hybrid column that does not mark the record as a hybrid (e.g. 
#   "species") and remove that text; the following will need edits !!!...

# create concatenated taxon_full_name column
all_data2 <- tidyr::unite(all_data2, "taxon_full_name_concat",
                          c(genus,hybrid,species,infra_rank,infra_name,cultivar),
                          sep=" ", remove=F,
                          na.rm=T)

# make sure there is a space b/w the author and taxon name in taxon_full_name
all_data2$taxon_full_name <- gsub("\\("," (",all_data2$taxon_full_name)
all_data2$taxon_full_name <- str_squish(all_data2$taxon_full_name)
unique(all_data2$taxon_full_name)

# add space after periods in taxon_full_name
# all_data2$taxon_full_name <- gsub(".",". ",all_data2$taxon_full_name,fixed=T)
# all_data2$taxon_full_name <- str_squish(all_data2$taxon_full_name)
# didn't do this 

################################################################################
# 5. Further standardize taxon name, then keep data for target taxa only
################################################################################

### manually fix taxon name issues you notice...
all_data2$taxon_full_name <- mgsub(all_data2$taxon_full_name,
 c("Araucaria angustifolia (Bertol. )Kunze","Araucaria heterophylla (Salisb. )
   Franco",
    "Agathis macrophylla (Lindl. )Mast. (syn. A. brownie (Lem. )L. H. Bailey", 
    "Araucaria angustifolia (bertol. ) Kuntze (syn. A. brasiliensis A. Rich)",
    "Agathis spec.","Agathis celebica Sulawesi'","Agathis Australis","Agathis 
   spec.",
    "Agathis robusta subsp. Nesophila"),
 c("Araucaria angustifolia (Bertol. ) Kunze","Araucaria heterophylla (Salisb. ) 
   Franco",
    "Agathis macrophylla (Lindl. ) Mast. ","Araucaria angustifolia (bertol. ) 
   Kuntze","Agathis sp.",
    "Agathis celebica 'Sulawesi'","Agathis australis","Agathis sp.","Agathis 
   robusta subsp. nesophila"))

all_data2$taxon_full_name_concat <- mgsub(all_data2$taxon_full_name_concat,
  c("Agathis celebica Sulawesi'"), c("Agathis celebica 'Sulawesi'"))

all_data2$cultivar <- mgsub(all_data2$cultivar, c("Glauca'"), c("'Glauca'"))

all_data2$genus <- mgsub(all_data2$genus, c("AGATHIS", "ARAUCARIA","WOLLEMIA"), 
                             c("Agathis","Araucaria","Wollemia"))

# have also done more manual editing within original google sheet and reuploaded 
# when easier 

# first change hybrid notation temporarily: remove space so it stays together
all_data2$taxon_full_name <- gsub(" x "," +",all_data2$taxon_full_name)

# separate out taxon full name and trim white space again
# this warning is ok: "Expected 9 pieces. Additional pieces discarded..."
all_data2 <- all_data2 %>% separate("taxon_full_name",
                                    c("genus_new","species_new","extra1","extra2",
                                      "extra3","extra4","extra5","extra6","extra7"), 
                                    sep=" ", extra="warn",
                                    remove=F, fill="right")
all_data2 <- as.data.frame(lapply(all_data2, str_squish), stringsAsFactors=F)
# replace genus_new with genus, since we fixed that up in the previous section
all_data2$genus_new <- all_data2$genus


## REMOVE RECORDS WITHOUT SPECIFIC EPITHET

# remove records with no/non-standard specific epithet (mostly cultivars or
#   'sp.' or questionable '?') by looking in species name column;
# if you WANT cultivars, you may need to edit this!
all_data3 <- all_data2 %>%
  filter(!grepl("\"",species_new) &
           !grepl("\'",species_new) &
           !grepl("\\[",species_new) &
           !grepl("\\(",species_new) &
           !grepl("\\.",species_new) &
           !grepl("[A-Z]",species_new) &
           !grepl("[0-9]",species_new) &
           !grepl("\\?",species_new) &
           !is.na(species_new))
nrow(all_data3)
# see names for records removed; can add anything you want to fix to the
#   "manually fix taxon name issues you notice" section above and rerun 
#   starting where you first created all_data2
sort(unique(suppressMessages(anti_join(all_data2,all_data3))$taxon_full_name))

# to view which rows have been removed
# removed_rows_2_3 <- setdiff(all_data2, all_data3)
# print(removed_rows_2_3)

## FIND INFRATAXA

# look for infrataxa key words
# make data in all "extra" columns lower case
sp_col <- grep("^species_new$", colnames(all_data3))
all_data3[,sp_col:(sp_col+7)] <- as.data.frame(sapply(
  all_data3[,sp_col:(sp_col+7)], tolower), stringsAsFactors=F)
# create matrix of all "extra" species name columns, to search for
#   infraspecific key words
search.col <- matrix(cbind(all_data3$extra1,all_data3$extra2,all_data3$extra3,
                           all_data3$extra4,all_data3$extra5,all_data3$extra6,
                           all_data3$extra7),
                     nrow=nrow(all_data3))
# search the "extra" column matrix for matches to infraspecific key words
matches_i <- which(search.col=="variety"|search.col=="var"|search.col=="var."|
                     search.col=="v"|search.col=="v."|
                     search.col=="subspecies"|search.col=="subsp"|
                     search.col=="subsp."|search.col=="ssp"|search.col=="ssp."|
                     search.col=="subs."|search.col=="spp."|search.col=="sub."|
                     search.col=="infra"|
                     search.col=="forma"|search.col=="form"|search.col=="fma"|
                     search.col=="fo"|search.col=="fo."|search.col=="f"|
                     search.col=="f.",arr.ind=T)
matches_i[,2] <- matches_i[,2]+sp_col
# create new infra_rank column and fill with "extra" contents that matched
#   infraspecific key words
all_data3$infra_rank_new <- NA
all_data3[matches_i[,1],"infra_rank_new"] <- all_data3[matches_i]
unique(all_data3$infra_rank_new) # check results

# create new infra_name column and fill with next column over from "extra"
#   contents that matched infraspecific key word
all_data3$infra_name_new <- NA
matches_i[,2] <- matches_i[,2]+1
all_data3[matches_i[,1],"infra_name_new"] <- all_data3[matches_i]

# standardize infraspecific rank names
all_data3$infra_rank_new <- replace(all_data3$infra_rank_new,
                                    grep("^variety$|^var$|^v$|^v.$",
                                         all_data3$infra_rank_new), "var.")
all_data3$infra_rank_new <- replace(all_data3$infra_rank_new,
                                    grep("^subspecies$|^subsp$|^ssp$|
                                         ^ssp.$|^subs.$|^spp.$|^sub.$|^infra$",
                                         all_data3$infra_rank_new), "subsp.")
all_data3$infra_rank_new <- replace(all_data3$infra_rank_new,
                                    grep("^forma$|^form$|^fma$|^fo$|^fo.$|^f$",
                                         all_data3$infra_rank_new), "f.")

# names we identified as infraspecific 
unique(all_data3[which(!is.na(all_data3$infra_rank_new)),"taxon_full_name"])


## CREATE FINAL TAXON FULL NAME FOR FILTERING

# create new taxon full name column
all_data3$taxon_full_name <- NA
# select rows with infraspecific name and concatenate
yes_infra <- which(!is.na(all_data3$infra_rank_new) &
                     !is.na(all_data3$infra_name_new))
all_data3$taxon_full_name[yes_infra] <- paste(all_data3$genus_new[yes_infra],
                                              all_data3$species_new[yes_infra], 
                                              all_data3$infra_rank_new[yes_infra],
                                              all_data3$infra_name_new[yes_infra],
                                              sep=" ")
# select rows without infraspecific name and concatenate
all_data3$taxon_full_name[-yes_infra] <- paste(all_data3$genus_new[-yes_infra],
                                               all_data3$species_new[-yes_infra],
                                               sep=" ")
# switch hybrid symbol back to " x "
all_data3$taxon_full_name <- gsub(" \\+"," x ",all_data3$taxon_full_name)
sort(unique(all_data3$taxon_full_name))

# code to find a specific text piece
# which(grepl("NA nobilis", all_data3$taxon_full_name))

## FILTER OUT NON-TARGET TAXA

# rename some taxon name columns to preserve originals
all_data3 <- all_data3 %>%
  rename(taxon_name = taxon_full_name,
         genus_orig = genus,
         species_orig = species,
         infra_rank_orig = infra_rank,
         infra_name_orig = infra_name)
all_data3 <- all_data3 %>%
  rename(genus = genus_new,
         species = species_new,
         infra_rank = infra_rank_new,
         infra_name = infra_name_new)

# join dataset to taxa list
### ROUND 1: match full taxon name
# add genus_species column to taxon list
taxon_list$genus_species <- NA
taxon_list$genus_species <- sapply(taxon_list$taxon_name, function(x)
  unlist(strsplit(x," var. | subsp. | f. "))[1])
# join by taxon name
all_data4 <- full_join(all_data3,taxon_list)
table(all_data4$taxon_name_status) 
### ROUND 2: match just by species name
need_match2 <- all_data4[which(is.na(all_data4$taxon_name_status)),]
nrow(need_match2)
# remove columns from first taxon name match
need_match2 <- need_match2[,1:(ncol(all_data4)-ncol(taxon_list)+1)]
# remove taxon_name col from taxon data so it doesn't match
taxon_list_sp <- taxon_list %>% 
 arrange(taxon_name_status) %>%
 distinct(genus_species,.keep_all=T) %>%
 select(-taxon_name)
# create genus_species column
need_match2$genus_species <- paste(need_match2$genus,need_match2$species)
# new join by genus_species
need_match2 <- left_join(need_match2,taxon_list_sp); nrow(need_match2)
# bind together new matches and previously matched
matched1 <- all_data4[which(!is.na(all_data4$taxon_name_status)),]
all_data4 <- rbind(matched1,need_match2)
table(all_data4$taxon_name_status)
head(all_data4)

### CHECK UNMATCHED SPECIES; ADD TO SYNONYM LIST AS NECESSARY ###
check <- all_data4 %>% filter(is.na(taxon_name_status))
check <- data.frame(taxon_name = sort(unique(check$taxon_name)))
nrow(check); check

# all unmatched are just my not target taxa

# write file for checking, as desired
# IF YOU FIND MISSPELLINGS AND/OR ADDITIONAL SYNONYMS, YOU CAN ADD THEM TO
#   YOUR TARGET TAXA LIST AND GO BACK TO THE START OF SECTION 4 AND RUN AGAIN 
#   FROM THERE
write.csv(check, file.path(standard_exsitu,
                           paste0("ExSitu_Unmatched_Species_", Sys.Date(), 
                                  ".csv")),row.names = F)


# keep only matched names
all_data5 <- all_data4 %>% 
  filter(!is.na(taxon_name_status) & !is.na(inst_short))
nrow(all_data5)
# see target taxa with no data
unique(taxon_list$taxon_name_accepted)[
  !(unique(taxon_list$taxon_name_accepted) %in% (unique(all_data5$
                                                          taxon_name_accepted)))]
# character (0)

# see any names with an x in the full name but not the accepted name
unique(all_data5[which(grepl(" x ",all_data5$taxon_full_name_orig) &
                         !grepl(" x ",all_data5$taxon_name_accepted)),
                 c("taxon_full_name_orig","taxon_name_accepted")])
# there is 1 Araucaria araucana x angustifolia F2
# these are complicated hybrids we don't want; remove
all_data5 <- all_data5 %>%
  filter(!(grepl(" x ",all_data5$taxon_full_name_orig) &
             !grepl(" x ",all_data5$taxon_name_accepted)))


# final part for removing duplicate data - from previous years or networks
all_data6 <- all_data5
# have left out the rest of this code as don't have date duplicates
  # have potentially network duplicates e.g ICCP

# summary of genera for each institution
gen_summary <- all_data6 %>%
  arrange(genus) %>%
  rename(genera = genus) %>%
  group_by(inst_short) %>%
  mutate(
    genera = paste(unique(genera), collapse = '; ')) %>%
  ungroup() %>%
  distinct(inst_short,genera)
gen_summary

# write file
write.csv(gen_summary, 
          file.path(standard_exsitu,
                    paste0("Genera_Institutions_Summary_", Sys.Date(), ".csv")),
          row.names = F)




################################################################################
# 6. Standardize important columns
################################################################################

########################### taken from 1. of this 2-compile.R script 

# read in institution metadata file, for comparing list to data read in
inst_data <- read.csv(file.path(exsitu_dir,
                                "respondent_institution_data_table.csv"), 
                      stringsAsFactors = F) 
## CHECK ALL INSTITUTIONS ARE PRESENT...
### INSTITUTIONS IN THE DATA YOU READ IN BUT NOT IN METADATA TABLE
# this should be "character(0)"
unique(all_data$inst_short)[!(unique(all_data$inst_short) %in% 
                                unique(inst_data$inst_short))]
### INSTITUTIONS IN METADATA TABLE BUT NOT IN THE DATA YOU READ IN:
# this should just be parent R files (have multiple institutions)
unique(inst_data$inst_short)[!(unique(inst_data$inst_short) %in% 
                                 unique(all_data$inst_short))]

# list.files("exsitu_data") - to check what files are uploaded correctly 

##########################

# add institution metadata
str(inst_data)
# SELECT COLUMNS YOU WANT ADDED TO ALL DATA:
inst_data <- inst_data %>%
  select(inst_short,inst_country,inst_lat,inst_long,inst_type) %>%
  arrange(inst_lat) %>%
  distinct(inst_short,.keep_all=T) %>%
  arrange(inst_short); inst_data
all_data7 <- left_join(all_data6,inst_data)
str(all_data7)


######################
## A) Provenance type
######################

# save original version of column for reference
all_data7$orig_prov_type <- all_data7$prov_type
all_data7$prov_type <- str_to_lower(all_data7$prov_type)

### look at column contents and CHANGE THE SECTIONS BELOW AS NEEDED...
sort(unique(all_data7$prov_type))
# standardize the column by searching for keywords & replacing with standard value
# first remove any confusing words/phrases
all_data7$prov_type <- mgsub(all_data7$prov_type,
                             c("not of known wild origin","could be cultivated",
                               "not of wild source",
                               "wildsourceunsure","cultivated plant of known ",
                               "natural source","unknown","g/u")," ")
# ex wild (Z)
all_data7$prov_type <- ifelse(grepl(paste(
  c("indirect","ex wild","^z$","cultivated from wild","c ex w","cw",
    "g from w plant","g ex w","cultivatedpropagatedfromwildsource",
    "(indirect) wild origin"),
  collapse = "|"), all_data7$prov_type),"Z",all_data7$prov_type)
# wild (W)
all_data7$prov_type <- ifelse(grepl(paste(
  c("wild","wld","collect","^w$","^w\\*","^\\(w\\)$","wd","w\\?","genetic",
    "100","110","130", "sao paulo w", "w s queensland", "w n queensland"),
  collapse = "|"), all_data7$prov_type),"W",all_data7$prov_type)
  # native to site (N)
all_data7$prov_type <- ifelse(grepl(paste(
  c("original to site","spontaneous","^n$","existing on site","native"),
  collapse = "|"), all_data7$prov_type),"N",all_data7$prov_type)
  # unknown (U)
all_data7$prov_type <- ifelse(grepl(paste(
  c("^\\(u\\)$","^u$","^u ","unsure","insufficient data","unknown","\\?","un",
    "need to confirm","breeding","400","410","416","donated","999","purchased","ul",
    "a"),
  collapse = "|"), all_data7$prov_type),"U",all_data7$prov_type)
  # cultivated (H)
all_data7$prov_type <- ifelse(grepl(paste(
  c("cultiva","garden","^c$","^g$","^g ","^h$","horticult","landrace","clone",
    "300","500","^g\\.", "cultivated", "cultivation", "new zealand, dow","G",
    "first release of wollemi","wollemi np","e"),
  collapse = "|"), all_data7$prov_type),"H",all_data7$prov_type)
## check one last time to be sure you got everything
sort(unique(all_data7$prov_type))
  # not given (NG) ; everything else
all_data7$prov_type <- ifelse(all_data7$prov_type!= "W" &
  all_data7$prov_type != "Z" & all_data7$prov_type != "H" &
  all_data7$prov_type != "N" & all_data7$prov_type != "U",
  "NG",all_data7$prov_type)
all_data7$prov_type[which(is.na(all_data7$prov_type))] <- "NG"

# which(grepl(" ", all_data7$prov_type))


# check results
table(all_data7$prov_type)


############################
## B) Number of Individuals
############################

## IF NEEDED: change # of individuals for rows that say dead/removed, etc., if
# you have a 'condition' column
sort(unique(all_data7$current_cond))
## CHANGE BASED ON PHRASE(S) IN YOUR DATA...
# all_data7[which(all_data7$condition=="Dead"),]$num_indiv <- "0"
# all_data7[which(all_data7$num_indiv=="all dead"),]$num_indiv <- "0"
# my current_cond is full of other information

# make everything lowercase; see what you have
all_data7$num_indiv <- str_to_lower(all_data7$num_indiv)
sort(unique(all_data7$num_indiv))

## IF NEEDED: replace unwanted characters
all_data7$num_indiv <- mgsub(all_data7$num_indiv,
                             c("\\?",","," plants"," pieces","ca\\.","alive",
                               "\\+",
                               " in terra",";[0-9][0-9]",";1",";2",";3",";4",
                               ";5",";7",";0",";"),
                             c(""), fixed=F)
all_data7$num_indiv <- str_squish(all_data7$num_indiv)
sort(unique(all_data7$num_indiv))

all_data7[which(all_data7$num_indiv=="large number (hundreds)"),]$num_indiv <- NA


# save version where we identify which didn't have # of individuals provided
all_data7$num_indiv[which(!grepl("^[0-9]+$",all_data7$num_indiv))] <- "Unknown"
all_data7$orig_num_indiv <- all_data7$num_indiv

# now we change the type to numeric and replace NA with 1, so we can use this
#   in calculations and when combining duplicate records later
# "NAs introduced by coercion" warning ok
all_data7$num_indiv <- as.numeric(all_data7$num_indiv)
all_data7$num_indiv[which(is.na(all_data7$num_indiv))] <- 1

# check results
sort(unique(all_data7$num_indiv))

# remove records with no individuals (first save as separate file)
no_indiv <- all_data7[which(all_data7$num_indiv == 0),]
nrow(no_indiv)
write.csv(no_indiv, file.path(standard_exsitu,
                              paste0("ExSitu_Dead_", Sys.Date(), ".csv")),
          row.names = F)
# save to in situ data folder as well
if(!dir.exists(file.path(raw_occ,"Ex_situ")))
  dir.create(file.path(raw_occ,"Ex_situ"),
             recursive=T)
write.csv(no_indiv, file.path(raw_occ,"Ex_situ",
                              paste0("ExSitu_Dead_", Sys.Date(), ".csv")),
          row.names = F)
# remove records with no individuals
all_data7 <- all_data7[which(all_data7$num_indiv > 0),]
nrow(all_data7)

write.csv(all_data7, file.path(raw_occ,"Ex_situ",
                              paste0("All_Data7_Export_Check_", Sys.Date(), 
                                     ".csv")),
          row.names = F)

#############################
## C) Latitude and Longitude
#############################

all_data8 <- all_data7

# create temporary ID col for use in this section
all_data8$temp_id <- seq.int(nrow(all_data8))

# preserve original lat and long columns
all_data8$lat_dd <- all_data8$orig_lat
all_data8$long_dd <- all_data8$orig_long

# replace comma with decimal (European notation)
all_data8$lat_dd <- mgsub(all_data8$lat_dd, c(","), ".")
all_data8$long_dd <- mgsub(all_data8$long_dd, c(","), ".")
# separate values if lat and long both ended up in the lat column
all_data8[which(grepl("\\. ",all_data8$lat_dd)),] <-
  separate(all_data8[which(grepl("\\. ",all_data8$lat_dd)),], col = lat_dd,
           into = c("lat_dd","long_dd"), sep = "\\. ", remove = FALSE)

# replace unwanted characters
## latitude
# replace unnecessary characters so we just have numbers
all_data8$lat_dd <- mgsub(all_data8$lat_dd,
                          c("N","\\","/","M","A",": ","E","AZ","R","d","a"," \\.",
                            " W")," ")
# remove leading zero
all_data8$lat_dd[which(grepl("^ *[0][1-9]+",all_data8$lat_dd))] <- gsub(
  "^ *[0]","",all_data8$lat_dd[which(grepl("^ *[0][1-9]+",all_data8$lat_dd))])
all_data8$lat_dd[which(grepl("^S *[0][1-9]+",all_data8$lat_dd))] <- gsub(
  "^S *[0]","-",all_data8$lat_dd[which(grepl("^S *[0][1-9]+",all_data8$lat_dd))])
# add negative sign if south and remove "S"
all_data8$lat_dd[grep("S",all_data8$lat_dd,ignore.case=T)] <-
  paste("-",all_data8$lat_dd[grep("S",all_data8$lat_dd,ignore.case=T)],sep="")
all_data8$lat_dd <- gsub("S","",all_data8$lat_dd)
all_data8$lat_dd <- gsub("--","-",all_data8$lat_dd)
# remove double spaces or leading/trailing whitespace
all_data8$lat_dd <- str_squish(all_data8$lat_dd)
all_data8$lat_dd[all_data8$lat_dd==""] <- NA
sort(unique(all_data8$lat_dd))
# can check source of specific values that aren't formatted correctly

which(grepl("#", all_data8$lat_dd))
all_data8$lat_dd <- mgsub(all_data8$lat_dd, c("#"), c(NA))

#all_data8[which(all_data8$lat_dd == "422538"),]
## longitude
# replace unnecessary characters so we just have numbers
all_data8$long_dd <- mgsub(all_data8$long_dd,
                           c("E","\\","/","NR","d","A","a"," .","o","O","#")," ")
# remove leading zero
all_data8$long_dd[which(grepl("^ *[0][1-9]+",all_data8$long_dd))] <- gsub(
  "^ *[0]","",all_data8$long_dd[which(grepl("^ *[0][1-9]+",all_data8$long_dd))])
# add negative sign if west and remove "W"
all_data8$long_dd[which(grepl("^W *[0][1-9]+",all_data8$long_dd))] <- gsub(
  "^W *[0]","-",all_data8$long_dd[which(grepl("^W *[0][1-9]+",
                                              all_data8$long_dd))])
all_data8$long_dd[grep("W",all_data8$long_dd,ignore.case=T)] <-
  paste("-",all_data8$long_dd[grep("W",all_data8$long_dd,ignore.case=T)],sep="")
all_data8$long_dd <- gsub("W","",all_data8$long_dd)
all_data8$long_dd <- mgsub(all_data8$long_dd,c("--","- "),"-")
# remove double spaces or leading/trailing whitespace
all_data8$long_dd <- str_squish(all_data8$long_dd)
all_data8$long_dd[all_data8$long_dd==""] <- NA
sort(unique(all_data8$long_dd))

all_data8$long_dd <- mgsub(all_data8$long_dd, c("-0.2","N"), c(NA, NA))


# convert decimal-minutes-seconds (dms) to decimal degrees (dd)
#   [d, m, and s must be in the same cell, with 1 space between each value]
#   format = ## ## ## (DMS) OR ## ##.### (DM)
# mark rows that need to be converted
convert <- all_data8 %>% filter(grepl(" ",lat_dd) | grepl(" ",long_dd))
nrow(convert)
no_convert <- anti_join(all_data8, convert)
# separate by dec_min_sec and deg_dec_min then convert to decimal degrees
# latitude
unique(convert$lat_dd)
dms <- convert[which(str_count(convert$lat_dd," ") == 2),]; nrow(dms)
ddm <- convert[which(str_count(convert$lat_dd," ") == 1),]; nrow(ddm)
other <- convert[which((str_count(convert$lat_dd," ") != 1 &
                          str_count(convert$lat_dd," ") != 2) | is.na(str_count(convert$lat_dd," "))),]
nrow(other)
dms$lat_dd = measurements::conv_unit(dms$lat_dd, from = 'deg_min_sec', to = 'dec_deg')
ddm$lat_dd = measurements::conv_unit(ddm$lat_dd, from = 'deg_dec_min', to = 'dec_deg')
convert <- rbind(dms,ddm,other); nrow(convert)
# longitude
unique(convert$long_dd)
dms <- convert[which(str_count(convert$long_dd," ") == 2),]; nrow(dms)
ddm <- convert[which(str_count(convert$long_dd," ") == 1),]; nrow(ddm)
other <- convert[which((str_count(convert$long_dd," ") != 1 &
                          str_count(convert$long_dd," ") != 2) | is.na(str_count(convert$long_dd," "))),]
nrow(other)
dms$long_dd = measurements::conv_unit(dms$long_dd, from = 'deg_min_sec', to = 'dec_deg')
ddm$long_dd = measurements::conv_unit(ddm$long_dd, from = 'deg_dec_min', to = 'dec_deg')
convert <- rbind(dms,ddm,other); nrow(convert)
# join everything back together
all_data8 <- rbind(no_convert,convert); nrow(all_data8)


# check validity of lat and long
all_data8$lat_dd <- as.numeric(all_data8$lat_dd)
all_data8$long_dd <- as.numeric(all_data8$long_dd)
# if coords are both 0, set to NA
zero <- which(all_data8$lat_dd == 0 & all_data8$long_dd == 0)
all_data8$lat_dd[zero] <- NA; all_data8$long_dd[zero] <- NA
# flag non-numeric, not available, and invalid (lat > 90 | lat < -90 |
#   lon > 180 | lon < -180) coordinates 
coord_test <- cc_val(all_data8, lon = "long_dd",lat = "lat_dd",
                     value = "flagged", verbose = TRUE)
# try switching lat and long for invalid points and check validity again
all_data8[!coord_test,c("lat_dd","long_dd")] <-
  all_data8[!coord_test,c("long_dd","lat_dd")]
coord_test <- cc_val(all_data8,lon = "long_dd",lat = "lat_dd",
                     value = "flagged",verbose = TRUE)

sort(unique(all_data8$lat_dd))
all_data8$lat_dd[all_data8$lat_dd == -730059.00000] <- -73.0059
all_data8$lat_dd[all_data8$lat_dd == -730048.00000] <- -73.0048
all_data8$lat_dd[all_data8$lat_dd == -714233.00000] <- -71.4233


sort(unique(all_data8$long_dd))
all_data8$long_dd[all_data8$long_dd == -390809.000000] <- -39.0809
all_data8$long_dd[all_data8$long_dd == -374929.000000] <- -37.4929
all_data8$long_dd[all_data8$long_dd == -374905.000000] <- -37.4905

format(sort(unique(all_data8$long_dd)), scientific = FALSE)

sort(unique(all_data8$country))
# IF the taxon is native to N/S America, make longitude value negative (you need
#   a taxon_region column for this, or could update code to work with countries!)...
# edited this just to multiply with -1 as issues with numeric vs characters and NA
# coercion
all_data8 <- all_data8 %>% 
  mutate(long_dd = if_else(long_dd > 0 & country %in% c(
    "ARGENTINA: Provincia de Neuquen", "Bermuda", "Brazil/Chile", "Chile South",
    "CHL", "Hawaii", "Mexico", "North America", "Brazil", "Bresil", "Chili",
    "CL", "United States", "Argentina", "Brazil South", "Chile", "28B", "28E", 
    "27D" ),
  -1 * long_dd,
  long_dd ))


# make coords NA if they are still flagged
coord_test <- cc_val(all_data8, lon = "long_dd", lat = "lat_dd",
                     value = "flagged", verbose = TRUE)
all_data8[!coord_test,"lat_dd"] <- NA
all_data8[!coord_test,"long_dd"] <- NA
nrow(all_data8)


### now we'll work just with the geolocated points for a bit
all_data8$latlong_flag <- ""
have_coord <- all_data8 %>% filter(!is.na(lat_dd) & !is.na(long_dd))
nrow(have_coord)
no_coord <- anti_join(all_data8,have_coord)
no_coord$latlong_country <- ""
# add country-level information to fix potentially-switched lat/longs
# turn occurrence point data into a spatial object
geo_pts_spatial <- vect(cbind(have_coord$long_dd, have_coord$lat_dd),
                        atts=have_coord, crs="+proj=longlat +datum=WGS84")
# read in world countries layer created in 1-prep_gis_layers.R
world_polygons <- vect(file.path(gis_dir,"world_countries_10m",
                                 "world_countries_10m.shp"))

# select just the country name column we need
world_polygons <- world_polygons[,"admin"]
# add country polygon data to each point based on lat-long location
geo_pts <- terra::intersect(geo_pts_spatial,world_polygons)
# try switching lat and long for points in Antarctica
on_land <- as.data.frame(geo_pts); nrow(on_land)
on_land[which(on_land$admin == "Antarctica"),c("lat_dd","long_dd")] <-
  on_land[which(on_land$admin == "Antarctica"),c("long_dd","lat_dd")]
head(on_land)
# remove columns from first join and join to countries again
on_land <- on_land %>% select(-admin)
geo_pts_spatial <- vect(cbind(on_land$long_dd,on_land$lat_dd),
                        atts=on_land, crs="+proj=longlat +datum=WGS84")
on_land <- as.data.frame(terra::intersect(geo_pts_spatial,world_polygons))
on_land <- on_land %>% rename(latlong_country = admin)
nrow(on_land)
# check if points are in water and mark
# get points that have coords but didn't fall in a country;
# these are in the water
land_id <- unique(on_land$temp_id)
in_water <- have_coord %>% filter(!(temp_id %in% land_id))
in_water <- as.data.frame(in_water)
if(nrow(in_water >0)){
  in_water$latlong_flag <- "Given lat-long is in the water" 
  in_water$latlong_country <- "" }
nrow(in_water)
# bind all the points back together
all_data9 <- rbind(on_land,in_water,no_coord)
nrow(all_data9)


# mark lat-long for records with same inst lat-long and wild lat-long
all_data9$lat_round <- round(all_data9$lat_dd,digits=1)
all_data9$long_round <- round(all_data9$long_dd,digits=1)
all_data9$inst_lat_round <- round(all_data9$inst_lat,digits=1)
all_data9$inst_long_round <- round(all_data9$inst_long,digits=1)
garden_latlong <- all_data9 %>% filter(lat_round == inst_lat_round &
                                         long_round == inst_long_round & prov_type != "N")
unique(garden_latlong$inst_short)
nrow(garden_latlong)
all_data9[which(all_data9$temp_id %in% garden_latlong$temp_id),]$latlong_flag <-
  "Given lat-long is at institution, use only if native to grounds"

# there are 4 rows in garden_latlong so have changed in all_data9 manually below

all_data9$lat_round[all_data9$inst_short == "NaplesBG" & 
                      all_data9$acc_num == "2014-1353*B"] <- NA

all_data9$long_round[all_data9$inst_short == "NaplesBG" & 
                      all_data9$acc_num == "2014-1353*B"] <- NA

all_data9$lat_round[all_data9$inst_short == "	
FranciscoJavierClavijeroBG" & 
                      all_data9$acc_num == "1982-519"] <- NA

all_data9$long_round[all_data9$inst_short == "	
FranciscoJavierClavijeroBG" & 
                       all_data9$acc_num == "1982-519"] <- NA

all_data9$lat_round[all_data9$inst_short == "BermudaBG"] <- NA

all_data9$long_round[all_data9$inst_short == "BermudaBG"] <- NA


# you can make lat-long NA for these without checking them first, if desired:
# manually did above
# all_data9[all_data9$UID %in% garden_latlong$UID,]$lat_dd <- NA
# all_data9[all_data9$UID %in% garden_latlong$UID,]$long_dd <- NA
table(all_data9$latlong_flag) 

sum(all_data9$UID %in% garden_latlong$UID)

# add latlong_det (latlong determination) column
all_data9$latlong_det[which(all_data9$prov_type == "H")] <- "N/A (horticultural)"
all_data9$latlong_det[which(!is.na(all_data9$lat_dd) &
                              !is.na(all_data9$long_dd))] <- "Given in original record"
all_data9$latlong_det[which(all_data9$latlong_det == "")] <- NA
table(all_data9$latlong_det)


# where prov_type is "N/A (horticultural)" but lat-long is given, change to "H?"
# create new prov type column
all_data9$prov_type[which(all_data9$latlong_det == "Given in original record" &
                            all_data9$prov_type == "H")] <- "H?"
table(all_data9$prov_type)


#######################
## D) Collection year
#######################

# this is not usually vital and takes some effort depending on your dataset;
# if you decide it's not necessary to standardize collection year, just comment 
# out this section; if you do want to standardize, you'll need to play with 
# updating the following section for your data


#####################
## E) Lineage number
#####################

# remove lin_num when same as acc_num
all_data9[which(all_data9$acc_num == all_data9$lin_num),]$lin_num <- NA


###########################
## F) Locality description
###########################

# create all_locality column (concatenating all locality data into one column 
#   for easy reference, especially when geolocating)
all_data9$latitude <- round(all_data9$lat_dd,digits=4)
all_data9$longitude <- round(all_data9$long_dd,digits=4)
all_data9 <- unite(all_data9, "all_locality",
                   c(locality,municipality,county,state,country,orig_source,
                     lin_num,coll_num,coll_name,coll_year,
                     latitude,longitude,notes),sep = " | ",remove = F)
# remove NA in concatenated locality column
all_data9$all_locality <- gsub("NA","",all_data9$all_locality)
# if no locality info at all, make it NA
all_data9$all_locality[which(all_data9$all_locality ==
                               " |  |  |  |  |  |  |  |  |  |  |  | ")] <- NA
head(all_data9$all_locality)


##################################################################
## H) Combine individuals (same institution and accession number)
##      so that everything is (hopefully) at the accession-level
##################################################################

## some institutions provided data at the accession level and some provided it
## at the individual level; we want everything at the accession level for
## later analyses, so we will try to combine individuals from the same 
## accession; this can be tricky...

# preserve original acc_num before removing individual signifiers in a minute
all_data9$orig_acc_num <- all_data9$acc_num

all_data9 <- all_data9 %>%
  mutate(UID = paste(inst_short, acc_num, prov_type, taxon_name_accepted, 
                     sep = "~"))


# have made a new all_data10 instead of keeping it on 9 just in case 
sort(names(all_data9))
## KEEP ONLY NECESSARY COLUMNS
#   change this as needed based on the columns you want;
keep_col <- c(
  "UID","inst_short","taxon_name_accepted","acc_num","prov_type",
  "lat_dd","long_dd","latlong_flag","latlong_det", "latlong_uncertainty",
  "geolocated_by","latlong_notes",
  "all_locality","locality","municipality","county","state","country",
  "latlong_country","assoc_sp",
  "orig_source","lin_num","coll_name","coll_num","coll_year",
  "num_indiv","germ_type",
  "notes","data_source",
  "taxon_name","taxon_full_name_orig","taxon_full_name_concat",
  "genus","species","infra_rank","infra_name","hybrid","cultivar",
  "taxon_name_status","taxon_verif",
  "inst_country","inst_lat","inst_long","inst_type",
  "orig_prov_type","orig_acc_num","orig_num_indiv","orig_lat","orig_long",
  "rl_category","all_native_dist_iso2"
)
all_data10 <- all_data9[,keep_col]

# save version without duplicates combined, in cases needed for reference
write.csv(all_data10, file.path(standard_exsitu,
                               paste0("ExSitu_Compiled_No-Dups-Combined_", 
                                      Sys.Date(), ".csv")),row.names = F)

# fill any spaces in acc_num with dashes
all_data10$acc_num <- gsub(" ","-",all_data10$acc_num)


# creating a standardised 8 number accession for duplicate code below
all_data10$acc_num_clean <- sapply(all_data10$orig_acc_num, function(x) {
  # Check for NA or empty input
  if (is.na(x) || x == "") {
    return(NA)
  }
  
  # Extract only digits (keep as character string)
  digits_only <- gsub("\\D", "", x)
  
  # If result is empty, return NA
  if (nchar(digits_only) == 0) {
    return(NA)
  }
  
  # Pad with leading zeros to ensure 13-digit string
  sprintf("%013s", digits_only)
})

sum(is.na(all_data10$acc_num_clean))


# combine duplicates (same acc_num) - this shouldn't happen but can
# all_data11 <- all_data10 %>%
 # group_by(inst_short,acc_num,taxon_name_accepted) %>%
#  mutate(num_indiv = sum(as.numeric(num_indiv)),
       #  orig_num_indiv = paste(unique(orig_num_indiv),collapse="; "),
       #  germ_type = paste(unique(germ_type),collapse="; "),
       #  inst_country = paste(unique(inst_country),collapse="; "),
       #  orig_acc_num = paste(unique(orig_acc_num),collapse="; "),
       #  all_locality = paste(unique(all_locality),collapse="; ")) %>%
 # ungroup() %>%
 # distinct(inst_short,acc_num,taxon_name_accepted,.keep_all=T)
#nrow(all_data11)

# removed_rows_10_11 <- setdiff(all_data10, all_data11)
# print(removed_rows_10_11)

##############################
## I) Add Universal ID column
##############################

# create UID with institution name, accession number, provenance type, and 
# taxon name; also remove duplicates based on new UID and sums individuals

# create UID
nms <- c("inst_short", "acc_num_clean", "prov_type", "taxon_name_accepted", 
         "num_indiv", "orig_lat", "Locality")

all_data11 <- all_data10 %>%
  mutate(UID = paste(inst_short, acc_num_clean, prov_type, 
                     taxon_name_accepted, sep = "~")) %>%
  group_by(UID) %>%
  summarise(
    inst_short = first(inst_short),
    acc_num_clean = first(acc_num_clean),
    prov_type = first(prov_type),
    taxon_name_accepted = first(taxon_name_accepted),
    orig_lat = first(orig_lat),
    Locality = first(locality),
    num_indiv = sum(as.numeric(num_indiv), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  select(c("UID", all_of(nms)))

nrow(all_data11)
# all_data11 is a summary of each target taxa and how many at the garden ######
# IE IT DOESN"T SORT DUPLICATES 


##############################
# Duplicates code specific to our dataset - made by Isabelle Grovenstein (MSc coursemate 24/25)
##############################

#--------------Duplicates-----------------
data_sel <- cbind(
  all_data10[, keep_col, drop = FALSE],
  acc_num_clean = all_data10$acc_num_clean
)

data_all <- data_sel

library(purrr)
library(dplyr)
library(stringr)
library(purrr)

# -----------------------------------------------------------------------------
# 1. Find all matching pairs (same taxon_name, no identical acc_num,
#    and an 8+ character substring of one acc_num appears in the other's all_locality)
# -----------------------------------------------------------------------------
# helper: does any substring of acc (length>=min_len) appear in loc?
acc_in_locality <- function(acc, loc, min_len = 13) {
  if (is.na(acc) || is.na(loc) || nchar(acc) < min_len) return(FALSE)
  subs <- map_chr(1:(nchar(acc) - min_len + 1),
                  ~ substr(acc, ., . + min_len - 1))
  any(str_detect(loc, fixed(subs, ignore_case = TRUE)))
}
# helper: extract that first matching substring
first_acc_match <- function(acc, loc, min_len = 13) {
  subs <- map_chr(1:(nchar(acc) - min_len + 1),
                  ~ substr(acc, ., . + min_len - 1))
  matched <- subs[str_detect(loc, fixed(subs, ignore_case = TRUE))]
  matched[1]
}

# build all unordered index pairs
pairs <- combn(nrow(data_all), 2, simplify = FALSE)

# scan for valid matches
match_pairs <- map_dfr(pairs, function(idx) {
  i <- idx[1]; j <- idx[2]
  a_i <- as.character(data_all$acc_num_clean[i])
  a_j <- as.character(data_all$acc_num_clean[j])
  loc_i <- data_all$all_locality[i]
  loc_j <- data_all$all_locality[j]
  tx_i  <- data_all$taxon_name[i]
  tx_j  <- data_all$taxon_name[j]
  
  # must share the same non-NA taxon_name
  if (is.na(tx_i) || is.na(tx_j) || tx_i != tx_j) return(NULL)
  # no exact acc_num match
  if (!is.na(a_i) && !is.na(a_j) && a_i == a_j) return(NULL)
  
  # test substring both ways
  i_in_j <- acc_in_locality(a_i, loc_j, 8)
  j_in_i <- acc_in_locality(a_j, loc_i, 8)
  if (!(i_in_j || j_in_i)) return(NULL)
  
  # record whichever substring triggered
  ms <- if (i_in_j) first_acc_match(a_i, loc_j, 8)
  else            first_acc_match(a_j, loc_i, 8)
  
  tibble(row1 = i, row2 = j, match_str = ms)
})

# -----------------------------------------------------------------------------
# 2. Prepare data_all for merging
# -----------------------------------------------------------------------------
data_work <- data_all %>%
  mutate(
    row_id         = row_number(),
    orig_num_indiv = num_indiv  # save original
  )

# partner lookup for later join
partner_lookup <- data_work %>%
  select(partner_id = row_id,
         partner_UID = UID,
         partner_orig = orig_num_indiv)

# -----------------------------------------------------------------------------
# 3a. Expand into primary/secondary roles
# -----------------------------------------------------------------------------
matches_long <- match_pairs %>%
  transmute(
    row_id      = row1,
    partner_id  = row2,
    match_str,
    role        = "primary"
  ) %>%
  bind_rows(
    match_pairs %>%
      transmute(
        row_id      = row2,
        partner_id  = row1,
        match_str,
        role        = "secondary"
      )
  ) %>%
  left_join(partner_lookup, by = "partner_id")

# -----------------------------------------------------------------------------
# 3b. **De‐duplicate** so each row_id appears only once**  
#     (here: take the first match per row_id)
# -----------------------------------------------------------------------------
matches_unique <- matches_long %>%
  group_by(row_id) %>%
  slice(1) %>%
  ungroup()

# -----------------------------------------------------------------------------
# 4. Merge back and compute new columns
# -----------------------------------------------------------------------------
result <- data_work %>%
  left_join(matches_unique, by = "row_id") %>%
  mutate(
    match_info = if_else(
      is.na(match_str),
      NA_character_,
      paste(UID,
            partner_UID,
            match_str,
            orig_num_indiv,
            partner_orig,
            sep = ";")
    ),
    num_indiv = case_when(
      is.na(role)         ~ as.character(orig_num_indiv),
      role == "primary"   ~ as.character(orig_num_indiv + partner_orig),
      role == "secondary" ~ "dup"
    )
  ) %>%
  select(-row_id, -partner_id, -role, -match_str,
         -partner_UID, -partner_orig)

# `result` now has each original row exactly once,
# with:
#  • orig_num_indiv (your preserved original count)
#  • match_info (only for the first matched partner, or NA)
#  • num_indiv updated (sum for primary, “dup” for secondary, else original)
result
#---------------------------END

matched_rows_dups <- result %>%
  filter(!is.na(match_info))

 write.csv(matched_rows_dups, file.path(standard_exsitu,
                                paste0("13_acc_dups_ouput_test", 
                                      Sys.Date(), ".csv")),row.names = F)
 
 

 ##############################################################################
 # sorting out ROBUSTA ISSUE - as changed all original names that were just robusta
 # to robusta subsp. nesophila
 # and fine to run as no real nesophila in this ex situ data
 ##############################################################################

 
result_no_robusta <- result[result$taxon_name_accepted 
                             != "Agathis robusta subsp. nesophila", ]
 
write.csv(result_no_robusta, file.path(standard_exsitu,
                                        paste0("All_ExSitu_Compiled_withdups_", 
                                               Sys.Date(), ".csv")),row.names = F)


################################################################################
# 8. [OPTIONAL] Explore georeferencing needs and save file for geolocation
# original input was data_sel but changed to result_no_robusta so had match info column
################################################################################

# if desired, you can try to manually find the latitude and longitude for 
# records with wild collection locality descriptions; the remaining sections
# help with this process

### explore georeferencing needs
# table with...
#   species name
#   num wild acc
#   num non-H acc w/ coords
#   num non-H acc with no coords & yes locality info
geo_needs <- result_no_robusta %>%
  group_by(taxon_name_accepted) %>%
  summarize(
    num_acc = sum(!is.na(taxon_name_accepted)),
    num_wild = sum(prov_type == "W"),
    NotH_YesCoords = sum(!is.na(lat_dd) & prov_type != "H"),
    NotH_NoCoords_YesLocality = sum(is.na(lat_dd) & 
                                      !is.na(all_locality) & 
                                      prov_type != "H"),
    Percent_NonH_NeedGeo = (sum(is.na(lat_dd) & 
                                  !is.na(all_locality) & 
                                  prov_type != "H") 
                            / sum(prov_type != "H")*100)
  )
head(geo_needs,n=20)
# write file
write.csv(geo_needs, file.path(standard_exsitu,
                               paste0("ExSitu_Geolocation_Needs_Summary_",
                                      Sys.Date(), ".csv")),
          row.names = F)

# select records that may need geolocation
#   (no lat-long, yes locality, prov type not H)
#   (also add flagged records: water or at institution)
need_geo <- result_no_robusta %>%
  filter((is.na(lat_dd) & prov_type != "H" &
            !is.na(all_locality) & all_locality != "NA") |
           latlong_flag!="")
nrow(need_geo)


# optionally, flag records for taxa that are highest priority for geolocating;
#   e.g., threatened and/or have less than 15 wild accessions
#   (you can choose whatever threshold(s) you want)
# thresholds
rl_threat <- c("Critically Endangered","Endangered","Vulnerable")
ns_threat <- c("G1 (Critically Imperiled)","G2 (Imperiled)","G3 (Vulerable)")
few_wild <- geo_needs[geo_needs$num_wild<15,]$taxon_name_accepted
# get list of priority taxa
priority_taxa <- taxon_list %>%
  filter(rl_category %in% rl_threat
           %in% ns_threat |
           taxon_name %in% few_wild) %>%
  distinct(taxon_name_accepted)
priority_taxa <- priority_taxa[,1]
priority_taxa <- paste(priority_taxa, collapse="|")
# flag priority taxa
need_geo$priority <- NA
need_geo[which(grepl(priority_taxa,need_geo$taxon_name_accepted)),]$priority <- "Priority"
table(need_geo$priority)

# replace NA with "" for easier viewing when geolocating
need_geo$lat_dd <- as.character(need_geo$lat_dd)
need_geo$long_dd <- as.character(need_geo$long_dd)
need_geo[is.na(need_geo)] <- ""

# write file
write.csv(need_geo, file.path(standard_exsitu,
                              paste0("ExSitu_Need_Geolocation_", Sys.Date(),
                                     ".csv")),row.names = F)






##############################################################################



### NOW MANUALLY GEOLOCATE !
### INSTRUCTIONS FOR GEOLOCATING:
### https://docs.google.com/document/d/1RBUD6-ogLc7PRVkDJSIEKzkgvC6xekj9Q_kl1vzxhCs/edit?usp=share_link

### When you're done geolocating, save your file with the same name but add 
### "_Geolocated" to the end
### We will read this in next!

################################################################################
# 9. (If you geolocated) Add geolocated data, after manual geolocation;
#    save final output file
################################################################################

# read in all compiled ex situ data (exported above)
exsitu <- read.csv(file.path(standard_exsitu,
                             "All_ExSitu_Compiled_withdups_2025-06-19.csv"), 
                   header = T, colClasses="character", na.strings = c("NA",""))

##############################
# remove all dup rows and use that as exsitu object , as the geo_raw below has 
# all dups removed too 
# this way dups are saved in google drive but removed here for further occurrence analysis 
######################################

exsitu <- exsitu[exsitu$num_indiv != "dup", ]


# read in geolocated dataset
geo_raw <- read.csv(file.path(standard_exsitu,
                              "ExSitu_Need_Geolocation_2025-06-19_Geolocated.csv"), 
                    header = T, colClasses="character", na.strings = c("NA",""))

# check this is just NA, meaning no "priority" records are not geolocated
unique(geo_raw[which(is.na(geo_raw$latlong_det)),"priority"])

############# to add another UID number version so below the geolocated data doesn't get 
# put into matching UID rows that have no corresponding geo data 
# geo_raw <- geo_raw %>%
 #  mutate(UID_row_number = row_number())

# exsitu <- exsitu %>%
 #  mutate(UID_row_number = row_number())



# add geolocated coordinates to ex situ data
# separate UID row
geolocated <- separate_rows(geo_raw, UID, sep=" \\| ")
# keep only edited columns and records that have latlong_det filled in
geolocated <- geolocated %>%
  select(UID,prov_type,lat_dd,long_dd,latlong_uncertainty,latlong_det,
         geolocated_by) %>%
  filter(!is.na(latlong_det))
head(geolocated)
table(geolocated$latlong_det)
# select geolocated rows in full dataset and remove cols we want to add
exsitu_geo <- exsitu %>%
  filter(UID %in% geolocated$UID) %>%
  select(-prov_type,-lat_dd,-long_dd,-latlong_uncertainty,-latlong_det,
         -geolocated_by)
# these two values should be the same:
nrow(exsitu_geo)
nrow(geolocated)


# add geolocation data; should just print "Joining with `by = join_by(UID)`"
exsitu_geo <- full_join(exsitu_geo,geolocated)
# join geolocated rows with rest of ex situ rows
exsitu_no_geo <- exsitu %>%
  filter(!(UID %in% exsitu_geo$UID))
nrow(exsitu_no_geo)
exsitu_all <- bind_rows(exsitu_no_geo,exsitu_geo)
nrow(exsitu_all)
table(exsitu_all$latlong_det)



#########################
# edited code above and created code below with chat - not useful now
#########################

# separate UID row
# geolocated <- separate_rows(geo_raw, UID, sep=" \\| ")

# add UID_row_number to geolocated
# geolocated <- geolocated %>%
#  group_by(UID) %>%
#  mutate(UID_row_number = row_number()) %>%
#  ungroup()

# add UID_row_number to exsitu
# exsitu <- exsitu %>%
#  group_by(UID) %>%
 # mutate(UID_row_number = row_number()) %>%
  #ungroup()

# select geolocated rows in full dataset and remove cols we want to add
# exsitu_geo <- exsitu %>%
 # filter(paste0(UID, "_", UID_row_number) %in% paste0(geolocated$UID, "_", 
# geolocated$UID_row_number)) %>%
 # select(-prov_type,-lat_dd,-long_dd,-latlong_uncertainty,-latlong_det,
  # -geolocated_by)

# these two values should be the same:
# nrow(exsitu_geo)
# nrow(geolocated)


# diff_in_data <- anti_join(geolocated, exsitu_geo, by = c("UID", "UID_row_number"))

# add geolocation data; join by both UID and UID_row_number
# exsitu_geo <- full_join(exsitu_geo, geolocated, by = c("UID", "UID_row_number"))

# join geolocated rows with rest of ex situ rows
# exsitu_no_geo <- exsitu %>%
 # filter(!(paste0(UID, "_", UID_row_number) %in% paste0(exsitu_geo$UID,
# "_", exsitu_geo$UID_row_number)))

# nrow(exsitu_no_geo)

# exsitu_all <- bind_rows(exsitu_no_geo, exsitu_geo)
# nrow(exsitu_all)
# table(exsitu_all$latlong_det)


###########################################
# end of edited section
##########################################


## FINAL VERSION: SELECT ONLY THE COLUMNS YOU WANT
# this is the same as the version created above in step 6H, just copying here
#   for easy access if not running the whole script at once
keep_col <- c(
  # key data
  "UID","inst_short","taxon_name_accepted","acc_num","prov_type",
  # locality
  "lat_dd","long_dd","latlong_flag","latlong_det","latlong_uncertainty",
  "geolocated_by",
  "all_locality","locality","municipality","county","state","country",
  "latlong_country","assoc_sp",
  # source
  "orig_source","lin_num","coll_name","coll_num","coll_year",
  # material info
  "num_indiv","germ_type",
  # other metadata
  "notes","data_source",
  # taxon name details
  "taxon_name","taxon_full_name_orig","taxon_full_name_concat",
  "genus","species","infra_rank","infra_name","hybrid","cultivar",
  "taxon_name_status","taxon_verif",
  # institution metadata
  "inst_country","inst_lat","inst_long","inst_type",
  # original versions of columns, for reference
  "orig_prov_type","orig_acc_num","orig_num_indiv","orig_lat","orig_long",
  # OPTIONAL additional taxon metadata
  "rl_category","all_native_dist_iso2","match_info"
)
exsitu_all <- exsitu_all[,keep_col]

# write final file
write.csv(exsitu_all, file.path(standard_exsitu,
                                paste0("FINAL_ExSitu_Compiled_Post-Geolocation_", Sys.Date(), ".csv")), 
          row.names = F)
# write to in situ folder also
write.csv(exsitu_all, file.path(raw_occ, "Ex_situ",
                                paste0("ExSitu_Compiled_Post-Geolocation_", Sys.Date(), ".csv")), 
          row.names = F)
