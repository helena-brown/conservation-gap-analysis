library(knitr)
library(kableExtra)
library(RColorBrewer)

################################################################################
# IUCN INTRODUCTION FIGURE
################################################################################

IUCN <- read.csv(file.path(analysis_dir,"IUCN_intro_figure2.csv"), 
                        header=T)

IUCN <- IUCN[, -1]

# gradient color 
colors <- brewer.pal(9, "OrRd")[4:6]

# making the table 
kable(IUCN, booktabs = TRUE,
      col.names = c("IUCN Red List Status",
                    "Number of target taxa",
                    "% of total Araucariaceae taxa",
                    "% of total conifer taxa")) %>%
  kable_styling(font_size = 30,
                latex_options = "scale_down") %>%
  row_spec(5, background = colors[1]) %>%
  row_spec(4, background = colors[2]) %>%
  row_spec(3, background = colors[3])

## other colours 
 # row_spec(3, color = 'white', background = '#FF4500') %>%
 # row_spec(4, color = 'white', background = '#FF8C00') %>%
 # row_spec(5, color = 'white', background = '#FFA500')



################################################################################
# SUMMARY STATS
################################################################################

library(dplyr)
library(readr)
library(purrr)
library(stringr)

data_in <- "taxon_points_final"

csv_files <- list.files(path = "~/Documents/taxon_points_final", pattern = ".csv", full.names = TRUE)


# read and combine all CSV files
exsitu_summary <- lapply(csv_files, function(file) {
  read.csv(file, stringsAsFactors = FALSE, colClasses = "character")
})

exsitu_summary <- bind_rows(exsitu_summary) %>%
  filter(database == "Ex_situ")

##################################
# ex situ lat long stats
##################################

# unique coordinates in lat or long for each taxa ##############################

unique_lat_summary <- exsitu_summary %>%
  group_by(taxon_name_accepted) %>%
  summarise(
    unique_decimalLatitude = list(unique(decimalLatitude)), # list of unique lat_dd
    num_unique_lat_dd = n_distinct(decimalLatitude),  # count of unique lat_dd
    .groups = "drop"
  )

print(unique_lat_summary)



# sum number of individuals for each taxa with lat and long ####################

exsitu_summary$individualCount <- as.numeric(exsitu_summary$individualCount)

summary_counts_individuals <- exsitu_summary %>%
  group_by(taxon_name_accepted) %>%
  summarise(total_individualCount = sum(individualCount, na.rm = TRUE))


# sum number of individuals from all data (GBIF etc) for each taxa with lat and long ####################

data_in <- "taxon_points_final"

csv_files <- list.files(path = "~/Documents/taxon_points_final", pattern = ".csv", full.names = TRUE)


# read and combine all CSV files
individuals_summary <- lapply(csv_files, function(file) {
  read.csv(file, stringsAsFactors = FALSE, colClasses = "character")
})

individuals_summary <- bind_rows(individuals_summary)

individuals_summary$individualCount <- as.numeric(individuals_summary$individualCount)

individuals_summary$individualCount[is.na(individuals_summary$individualCount)] <- 1

summary_counts_individuals_all <- individuals_summary %>%
  group_by(taxon_name_accepted) %>%
  summarise(total_individualCount = sum(individualCount, na.rm = TRUE))

##################################
# ex situ individual vs geographic coverage table 
##################################

coverage <- read.csv(file.path(analysis_dir,"individuals_compared_geographic.csv"), 
                           header=T)


# making the table 
coverage %>%
  kbl(col.names = c("Target taxon",
                    "Total number of all known occurrence points mapped",
                    "Number of ex situ occurence points mapped",
                    "Percentage individual occurrence coverage",
                    "Percentage geographic coverage")) %>%
  kable_classic(full_width = F) %>%
  kable_styling(font_size = 35) %>%
  row_spec(1, background = "#F0F0F0") %>%
  row_spec(3, background = "#F0F0F0") %>%
  row_spec(5, background = "#F0F0F0") %>%
  row_spec(7, background = "#F0F0F0") %>%
  row_spec(9, background = "#F0F0F0") %>%
  row_spec(11, background = "#F0F0F0") %>%
  row_spec(13, background = "#F0F0F0") %>%
  row_spec(15, background = "#F0F0F0") %>%
  row_spec(17, background = "#F0F0F0") %>%
  row_spec(19, background = "#F0F0F0") %>%
  row_spec(21, background = "#F0F0F0") %>%
  row_spec(23, background = "#F0F0F0") %>%
  row_spec(25, background = "#F0F0F0") %>%
  column_spec(1, italic = T)



# unique coordinates in lat or long for each taxa ##############################

exsitu_summary$coordinateUncertaintyInMeters <- as.numeric(
  exsitu_summary$coordinateUncertaintyInMeters)

# Group by taxon_name_accepted and calculate sum and mean
uncertainty_summary <- exsitu_summary %>%
  group_by(taxon_name_accepted) %>%
  summarise(
    total_uncertainty = sum(coordinateUncertaintyInMeters, na.rm = TRUE),
    average_uncertainty = mean(coordinateUncertaintyInMeters, na.rm = TRUE)
  )


##################################
# ex situ all stats 
##################################


# sum number of individuals for each taxa ######################################

exsitu_all_summary <- read.csv(file.path(
  raw_occ,"ExSitu_Compiled_Post-Geolocation_2025-06-19.csv"), 
                 header=T)

summary_counts_all_individuals <- exsitu_all_summary %>%
  group_by(taxon_name_accepted) %>%
  summarise(total_num_indiv = sum(num_indiv, na.rm = TRUE))

sum(!is.na(exsitu_all_summary$lat_dd))
# 1714 with NA
# 2415 with lat and long 

# sum number of institutions for each taxa #############

exsitu_compiled <- read.csv(file.path(
  standard_exsitu,"All_ExSitu_Compiled_withdups_2025-06-19.csv"), 
                 header=T)

unique_inst_exsitu <- exsitu_compiled %>%
  group_by(taxon_name_accepted) %>%
  summarise(
    unique_inst_short = list(unique(inst_short)),
    num_unique_inst_short = n_distinct(inst_short),  
    .groups = "drop"
  )







################################################################################
# SUMMARY STATS TABLE
################################################################################

##################################
# New Caledonia
##################################


stats_table_NC <- read.csv(file.path(analysis_dir,"stats_table_NC.csv"), 
                 header=T)


# making the table 
stats_table_NC %>%
  kbl(col.names = c("Target taxon",
                    "Number of ex situ individuals mapped",
                    "Number of ex situ wild collected locations mapped",
                    "Average uncertainty buffer radius for geolocated occurrences",
                    "Total number of individuals reported in ex situ wild collections",
                    "Percentage of total number of ex situ individuals mapped",
                    "Number of institutions with ex situ collections of the taxon")) %>%
  kable_classic(full_width = F) %>%
  kable_styling(font_size = 35) %>%
  row_spec(1, background = "#F0F0F0") %>%
  row_spec(3, background = "#F0F0F0") %>%
  row_spec(5, background = "#F0F0F0") %>%
  row_spec(7, background = "#F0F0F0") %>%
  row_spec(9, background = "#F0F0F0") %>%
  row_spec(11, background = "#F0F0F0") %>%
  row_spec(13, background = "#F0F0F0") %>%
  column_spec(1, italic = T)


##################################
# SE Asia
##################################


stats_table_SEA <- read.csv(file.path(analysis_dir,"stats_table_SEA.csv"), 
                           header=T)


# making the table 
stats_table_SEA %>%
  kbl(col.names = c("Target taxon",
                    "Number of ex situ individuals mapped",
                    "Number of ex situ wild collected locations mapped",
                    "Average uncertainty buffer radius for geolocated occurrences",
                    "Total number of individuals reported in ex situ wild collections",
                    "Percentage of total number of ex situ individuals mapped",
                    "Number of institutions with ex situ collections of the taxon")) %>%
  kable_classic(full_width = F) %>%
  kable_styling(font_size = 35) %>%
  row_spec(1, background = "#F0F0F0") %>%
  row_spec(3, background = "#F0F0F0") %>%
  row_spec(5, background = "#F0F0F0") %>%
  row_spec(7, background = "#F0F0F0") %>%
  column_spec(1, italic = T) 


##################################
# South America
##################################

stats_table_SA <- read.csv(file.path(analysis_dir,"stats_table_SA.csv"), 
                            header=T)


stats_table_SA %>%
  kbl(col.names = c("Target taxon",
                    "Number of ex situ individuals mapped",
                    "Number of ex situ wild collected locations mapped",
                    "Average uncertainty buffer radius for geolocated occurrences",
                    "Total number of individuals reported in ex situ wild collections",
                    "Percentage of total number of ex situ individuals mapped",
                    "Number of institutions with ex situ collections of the taxon")) %>%
  kable_classic(full_width = F) %>%
  kable_styling(font_size = 35) %>%
  row_spec(1, background = "#F0F0F0") %>%
  column_spec(1, italic = T)



##################################
# Pacific Islands and Australia
##################################


stats_table_PIAA <- read.csv(file.path(analysis_dir,"stats_table_PIAA.csv"), 
                            header=T)


# making the table 
stats_table_PIAA %>%
  kbl(col.names = c("Target taxon",
                    "Number of ex situ individuals mapped",
                    "Number of ex situ wild collected locations mapped",
                    "Average uncertainty buffer radius for geolocated occurrences",
                    "Total number of individuals reported in ex situ wild collections",
                    "Percentage of total number of ex situ individuals mapped",
                    "Number of institutions with ex situ collections of the taxon")) %>%
  kable_classic(full_width = F) %>%
  kable_styling(font_size = 35) %>%
  row_spec(1, background = "#F0F0F0") %>%
  row_spec(3, background = "#F0F0F0") %>%
  column_spec(1, italic = T)


#######################
# combined table
#####################



stats_table_NC <- read.csv(file.path(analysis_dir,"stats_table_NC.csv"), header=TRUE) %>%
  mutate(Region = "New Caledonia")

stats_table_SEA <- read.csv(file.path(analysis_dir,"stats_table_SEA.csv"), header=TRUE) %>%
  mutate(Region = "South-East Asia")

stats_table_SA <- read.csv(file.path(analysis_dir,"stats_table_SA.csv"), header=TRUE) %>%
  mutate(Region = "South America")

stats_table_PIAA <- read.csv(file.path(analysis_dir,"stats_table_PIAA.csv"), header=TRUE) %>%
  mutate(Region = "Pacific Islands, Australia, and New Zealand")

stats_table_NC <- stats_table_NC %>%
  mutate(across(everything(), as.character))

stats_table_SEA <- stats_table_SEA %>%
  mutate(across(everything(), as.character))

stats_table_SA <- stats_table_SA %>%
  mutate(across(everything(), as.character))

stats_table_PIAA <- stats_table_PIAA %>%
  mutate(across(everything(), as.character))

# Combine all rows into one dataframe
combined_stats <- bind_rows(stats_table_NC, stats_table_SEA, stats_table_SA, stats_table_PIAA)

# Reorder columns to put Region first (optional)
combined_stats <- combined_stats %>%
  select(Region, everything())


# change colour opacity 

hex_to_rgba <- function(hex, alpha = 0.5) {
  # Remove the # if present
  hex <- gsub("#", "", hex)
  
  # Convert hex to RGB
  r <- strtoi(substr(hex, 1, 2), 16L)
  g <- strtoi(substr(hex, 3, 4), 16L)
  b <- strtoi(substr(hex, 5, 6), 16L)
  
  # Return rgba string with opacity
  sprintf("rgba(%d, %d, %d, %.2f)", r, g, b, alpha)
}


# Render combined table
combined_stats %>%
  kbl(col.names = c("Region",
                    "Target taxon",
                    "Number of ex situ individuals mapped",
                    "Number of ex situ wild collected locations mapped",
                    "Average uncertainty buffer radius for geolocated occurrences",
                    "Total number of individuals reported in ex situ wild collections",
                    "Percentage of total number of ex situ individuals mapped",
                    "Number of institutions with ex situ collections of the taxon")) %>%
  kable_classic(full_width = F) %>%
  kable_styling(font_size = 30) %>%
  # Optional: add alternating row colors
  row_spec(1:13, background = hex_to_rgba("#F8766D", 0.6)) %>%
  row_spec(14:20, background = hex_to_rgba("#C77CFF", 0.6)) %>%
  row_spec(21:22, background = hex_to_rgba("#00BFC4", 0.6)) %>%
  row_spec(23:26, background = hex_to_rgba("#7CAE00", 0.6)) %>%
  column_spec(2, italic = TRUE) 



################################################################################
# CONSERVATION PRIORITISATION MATRIX 
################################################################################

consv_matrix_PA_exsitu <- read.csv(file.path(analysis_dir,"conservation_matrix_PA_exsitu.csv"), 
                        header=T)

consv_matrix_PA_exsitu <- consv_matrix_PA_exsitu[, -13]
consv_matrix_PA_exsitu <- consv_matrix_PA_exsitu[, -11]

consv_matrix_PA_exsitu$
  geographic.area[consv_matrix_PA_exsitu$geographic.area == 
                    "Pacific Islands and Australia"] <- "Pacific Islands, Australia, and New Zealand"

# colour coding score vulnerability - Wong colour palette - same as survey pie charts 
get_score_color <- function(score) {
  if (score <= 5) {
    "#009E73"  
  } else if (score <= 10) {
    "#F0E442"
  } else if (score <= 20) {
    "#E69F00"
  } else if (score <= 40) {
    "#D55E00"
  } else {
    "#CC79A7"  
  }
}

# apply colours 
consv_matrix_PA_exsitu <- consv_matrix_PA_exsitu %>%
  mutate(vulnerability.score = cell_spec(vulnerability.score,
                                                 background = sapply(vulnerability.score, get_score_color),
                                                 color = "white",
                                                 bold = TRUE,
                                                extra_css = "display: block; padding: 4px;"))

# gradient color 
color_5 <- brewer.pal(9, "OrRd")[6]
color_4 <- brewer.pal(9, "OrRd")[5]
color_3 <- brewer.pal(9, "OrRd")[4]
color_0 <- brewer.pal(9, "OrRd")[2]

# colour coding score vulnerability - OrRd colour palette 
get_score_color_IUCN <- function(score) {
  if (score == 5) {
    color_5  
  } else if (score == 4) {
    color_4
  } else if (score == 3) {
    color_3
  } else if (score == 0) {
    color_0
  } else {
    "#CC79A7"  
  }
}

# apply colours 
consv_matrix_PA_exsitu <- consv_matrix_PA_exsitu %>%
  mutate(IUCN.scoring = cell_spec(IUCN.scoring,
                                         background = sapply(IUCN.scoring, get_score_color_IUCN),
                                         color = "white",
                                         bold = TRUE,
                                        extra_css = "display: block; padding: 4px;"))


# making the table 
consv_matrix_PA_exsitu %>%
  kbl(col.names = c("Target taxon",
                    "Case study geographic area",
                    "IUCN scoring",
                    "Total number of living individuals in ex situ wild collections",
                    "Vulnerability score",
                    "Number of living individuals mapped for ex situ gap analysis",
                    "Percentage geographic coverage of ex situ collections",
                    "Percentage of all known occurrence points within protected areas",
                    "Number of institutions with ex situ collections of the taxon",
                    "Number of institutions with a conservation understanding of the taxon",
                    "IUCN category and criteria"),
                    escape = FALSE) %>%
  kable_classic(full_width = F) %>%
  kable_styling(font_size = 20) %>%
  row_spec(1, background = "#F0F0F0") %>%
  row_spec(3, background = "#F0F0F0") %>%
  row_spec(5, background = "#F0F0F0") %>%
  row_spec(7, background = "#F0F0F0") %>%
  row_spec(9, background = "#F0F0F0") %>%
  row_spec(11, background = "#F0F0F0") %>%
  row_spec(13, background = "#F0F0F0") %>%
  row_spec(15, background = "#F0F0F0") %>%
  row_spec(17, background = "#F0F0F0") %>%
  row_spec(19, background = "#F0F0F0") %>%
  row_spec(21, background = "#F0F0F0") %>%
  row_spec(23, background = "#F0F0F0") %>%
  row_spec(25, background = "#F0F0F0") %>%
 # row_spec(1:4, background = colors[3]) %>%
 # row_spec(5:13, background = colors[2]) %>%
 # row_spec(14:25, background = colors[1]) %>%
 # row_spec(26, background = color_2[1]) %>%
  column_spec(1, italic = T) %>%
  row_spec(0, bold = T)


################################################################################
# DIFFERENCE IN EOO TABLE
################################################################################

diff_in_EOO <- read.csv(file.path(analysis_dir,"difference_in_EOO.csv"), 
                         header=T)

diff_in_EOO %>%
  kbl(col.names = c("Target taxon",
                    "Extent of Occurrence (EOO) with no GBIF data",
                    "EOO with GBIF data",
                    "Percentage difference")) %>%
  kable_classic(full_width = F) %>%
  kable_styling(font_size = 10) %>%
  row_spec(1, background = "#F0F0F0") %>%
  row_spec(3, background = "#F0F0F0") %>%
  row_spec(5, background = "#F0F0F0") %>%
  row_spec(7, background = "#F0F0F0") %>%
  row_spec(9, background = "#F0F0F0") %>%
  row_spec(11, background = "#F0F0F0") %>%
  row_spec(13, background = "#F0F0F0") %>%
  row_spec(15, background = "#F0F0F0") %>%
  column_spec(1, italic = T)




################################################################################
# ECO AND GEO COVERAGE
################################################################################

eco_geo <- read.csv(file.path(analysis_dir,"eco_geo_coverage.csv"), 
                        header=T)

for (i in 1:nrow(eco_geo)) {
  
  # split csv file
  row_df <- eco_geo[i, , drop = FALSE]
  row_name <- make.names(eco_geo[i, 1])
  assign(row_name, row_df)
}
  

  # make table
  Araucaria.montana %>%
    kbl(col.names = c("Target taxon",
                      "Geographic coverage with small (20km) buffer",
                      "Geographic coverage with buffer size of taxon",
                      "Ecological coverage with small (20km) buffer",
                      "Ecological coverage with buffer size of taxon",
                      "Extent of Occurrence (EOO) (m\u00B2)")) %>%
    kable_classic(full_width = FALSE) %>%
    kable_styling(font_size = 10) %>%
    row_spec(1, background = "#F0F0F0") %>%
    column_spec(2, italic = T)
  
# copy and repeat for all other 15 taxa
  
  
  
################################################################################
# VULNERABILITY TABLES 
################################################################################
  
  library(dplyr)
  library(kableExtra)
  library(scales) 
  
  # colour coding score - Wong colour palette - same as survey pie charts 
  get_score_color <- function(score) {
    if (score <= 5) {
      "#009E73"  
    } else if (score <= 10) {
      "#F0E442"
    } else if (score <= 20) {
      "#E69F00"
    } else if (score <= 40) {
      "#D55E00"
    } else {
      "#CC79A7"  
    }
  }
  
  
  ##################################
  # New Caledonia
  ##################################
  
  
  vulnerability_NC <- read.csv(file.path(analysis_dir,"vulnerability_NC.csv"), 
                             header=T)
  
  
  # try replacing < with html format of < 
  vulnerability_NC <- vulnerability_NC %>%
    mutate(across(where(is.character), ~ gsub("<", "&lt;", .)))
  
  # apply colours 
  vulnerability_NC <- vulnerability_NC %>%
    mutate(Average.vulnerability.score = cell_spec(Average.vulnerability.score,
                                                   background = sapply(Average.vulnerability.score, get_score_color),
                                                   color = "white",
                                                   bold = TRUE,
                                                   extra_css = "display: block; padding: 4px;"))
  
  
  # make the table
  vulnerability_NC %>%
    kbl(col.names = c("Target taxon",
                      "Population size",
                      "Range/Endemism",
                      "Population decline",
                      "Fragmentation",
                      "Regeneration/Recruitment",
                      "Genetic variation/Integrity",
                      "Average vulnerability score"),
        escape = FALSE) %>%
    kable_classic(full_width = F) %>%
    kable_styling(font_size = 30) %>%
    row_spec(1, background = "#F0F0F0") %>%
    row_spec(3, background = "#F0F0F0") %>%
    row_spec(5, background = "#F0F0F0") %>%
    row_spec(7, background = "#F0F0F0") %>%
    row_spec(9, background = "#F0F0F0") %>%
    row_spec(11, background = "#F0F0F0") %>%
    row_spec(13, background = "#F0F0F0") %>%
    column_spec(1, italic = T)
  
  
  ##################################
  # SE Asia
  ##################################
  
  vulnerability_SEA <- read.csv(file.path(analysis_dir,"vulnerability_SEA.csv"), 
                               header=T)
  
  # try replacing < with html format of < 
  vulnerability_SEA <- vulnerability_SEA %>%
    mutate(across(where(is.character), ~ gsub("<", "&lt;", .)))
  
  
  # apply colours 
  vulnerability_SEA <- vulnerability_SEA %>%
    mutate(Average.vulnerability.score = cell_spec(Average.vulnerability.score,
                                                   background = sapply(Average.vulnerability.score, get_score_color),
                                                   color = "white",
                                                   bold = TRUE,
                                                   extra_css = "display: block; padding: 4px;"))
  
  
  # make the table
  vulnerability_SEA %>%
    kbl(col.names = c("Target taxon",
                      "Population size",
                      "Range/Endemism",
                      "Population decline",
                      "Fragmentation",
                      "Regeneration/Recruitment",
                      "Genetic variation/Integrity",
                      "Average vulnerability score"),
        escape = FALSE) %>%
    kable_classic(full_width = FALSE) %>%
    kable_styling(font_size = 30) %>%
    row_spec(1, background = "#F0F0F0") %>%
    row_spec(3, background = "#F0F0F0") %>%
    row_spec(5, background = "#F0F0F0") %>%
    row_spec(7, background = "#F0F0F0") %>%
    column_spec(1, italic = T) 
  
  
  ##################################
  # South America
  ##################################
  

  vulnerability_SA <- read.csv(file.path(analysis_dir, "vulnerability_SA.csv"), header = TRUE)
  
  # try replacing < with html format of < 
  vulnerability_SA <- vulnerability_SA %>%
    mutate(across(where(is.character), ~ gsub("<", "&lt;", .)))
  
  # apply colours 
  vulnerability_SA <- vulnerability_SA %>%
    mutate(Average.vulnerability.score = cell_spec(Average.vulnerability.score,
                                                   background = sapply(Average.vulnerability.score, get_score_color),
                                                   color = "white",
                                                   bold = TRUE,
                                                   extra_css = "display: block; padding: 4px;"))
  
  
  # make the table
  vulnerability_SA %>%
    kbl(col.names = c("Target taxon",
                      "Population size",
                      "Range/Endemism",
                      "Population decline",
                      "Fragmentation",
                      "Regeneration/Recruitment",
                      "Genetic variation/Integrity",
                      "Average vulnerability score"),
        escape = FALSE) %>%
    kable_classic(full_width = FALSE) %>%
    kable_styling(font_size = 30) %>%
    row_spec(1, background = "#F0F0F0") %>%
    column_spec(1, italic = TRUE)
  
  ##################################
  # Pacific Islands and Australia
  ##################################
  
  
  vulnerability_PIAA <- read.csv(file.path(analysis_dir,"vulnerability_PIAA.csv"), 
                               header=T)
  
  # try replacing < with html format of < 
  vulnerability_PIAA <- vulnerability_PIAA %>%
    mutate(across(where(is.character), ~ gsub("<", "&lt;", .)))
  
  # apply colours 
  vulnerability_PIAA <- vulnerability_PIAA %>%
    mutate(Average.vulnerability.score = cell_spec(Average.vulnerability.score,
                                                   background = sapply(Average.vulnerability.score, get_score_color),
                                                   color = "white",
                                                   bold = TRUE,
                                                   extra_css = "display: block; padding: 4px;"))
  
  
  # make the table
  vulnerability_PIAA %>%
    kbl(col.names = c("Target taxon",
                      "Population size",
                      "Range/Endemism",
                      "Population decline",
                      "Fragmentation",
                      "Regeneration/Recruitment",
                      "Genetic variation/Integrity",
                      "Average vulnerability score"),
        escape = FALSE) %>%
    kable_classic(full_width = F) %>%
    kable_styling(font_size = 30) %>%
    row_spec(1, background = "#F0F0F0") %>%
    row_spec(3, background = "#F0F0F0") %>%
    column_spec(1, italic = T)
  
  
################################################################################
# VULNERABILITY SCORING TABLE
################################################################################

  vulnerability_scoring <- read.csv(file.path(analysis_dir,"vulnerability_scoring.csv"), 
                                 header=T)
  
  vulnerability_scoring <- vulnerability_scoring %>%
    mutate(across(where(is.character), ~ gsub("<", "&lt;", .)))
  
  
  
  # make the table
  vulnerability_scoring %>%
    kbl(col.names = c("Category",
                      "Scoring",
                      "Population size",
                      "Range/Endemism",
                      "Population decline",
                      "Fragmentation",
                      "Regeneration/Recruitment",
                      "Genetic variation/Integrity"),
        escape = FALSE) %>%
    kable_classic(full_width = F) %>%
    kable_styling(font_size = 28) %>%
    row_spec(1, background = "#D55E00") %>%
    row_spec(2, background = "#E69F00") %>%
    row_spec(3, background = "#F0E442") %>%
    row_spec(4, background = "#009E73") %>%
    row_spec(5, background = "#56B4E9") %>%
    row_spec(6, background = "#56B4E9") %>%
    row_spec(7, background = "#999999") %>%
    column_spec(1, italic = T)
  
  
  ###########################
  # convert to rgba for less harsh colours
  ###########################
  
  vulnerability_scoring <- vulnerability_scoring %>%
    mutate(across(where(is.character), ~ gsub("<", "&lt;", .)))
  
  vulnerability_scoring %>%
    kbl(col.names = c("Category",
                      "Scoring",
                      "Population size",
                      "Range/Endemism",
                      "Population decline",
                      "Fragmentation",
                      "Regeneration/Recruitment",
                      "Genetic variation/Integrity"),
        escape = FALSE) %>%
    kable_classic(full_width = FALSE) %>%
    kable_styling(font_size = 28) %>%
    row_spec(1, background = "rgba(213, 94, 0, 0.6)") %>%    # softened orange-red
    row_spec(2, background = "rgba(230, 159, 0, 0.6)") %>%   # softened orange
    row_spec(3, background = "rgba(240, 228, 66, 0.6)") %>%  # softened yellow
    row_spec(4, background = "rgba(0, 158, 115, 0.6)") %>%   # softened green
    row_spec(5, background = "rgba(86, 180, 233, 0.6)") %>%  # softened blue
    row_spec(6, background = "rgba(86, 180, 233, 0.6)") %>%  # repeated blue
    row_spec(7, background = "rgba(153, 153, 153, 0.6)") %>% # softened gray
    column_spec(1, italic = TRUE)
  
  
  
  
###############################################################################
# THREAT TABLES
###############################################################################
 
 
 
 ##########################################
 # LOOP THROUGH ALL THREAT TABLES AND STATS 
 ##########################################
 

 
 
 threat_names <- c("NC", "SEA", "PIAA", "SA")
 results_list <- list()
 
 # Define your color palette function once
 color_5 <- brewer.pal(9, "OrRd")[6]
 color_4 <- brewer.pal(9, "OrRd")[5]
 color_3 <- brewer.pal(9, "OrRd")[4]
 color_0 <- brewer.pal(9, "OrRd")[2]
 
 get_score_color_threat <- function(score) {
   if (score == "high impact") {
     color_5  
   } else if (score == "moderate impact") {
     color_4
   } else if (score == "low impact") {
     color_3
   } else if (score == "not a significant threat") {
     color_0
   } else if (score == "cause of past decline") {
     "lightpink"  
   }
 }
 
 
 
 for (name in threat_names) {
   # Read data
   threat_df <- read.csv(file.path(analysis_dir, paste0("threat_", name, ".csv")), 
                         header = TRUE)
   
   # Clean whitespace
   threat_df <- threat_df %>%
     mutate(across(where(is.character), ~ gsub("\\s+", " ", trimws(.))))
   
   ############### Stats ###########################
   score_data <- threat_df %>% select(-Species.of.concern)
   
   top_cols_high <- colSums(score_data == "high impact", na.rm = TRUE) %>%
     sort(decreasing = TRUE)
   
   # create stats tables
   assign(paste0("top_cols_high_", name), 
          data.frame(
            Threat_Category = names(top_cols_high),
            High_Impact_Count = as.vector(top_cols_high)
          ),
          envir = .GlobalEnv)
   
   top_taxa_high_mod <- threat_df %>%
     mutate(high_mod_count = rowSums(score_data == "high impact" | 
                                     score_data == "moderate impact", 
                                     na.rm = TRUE)) %>%
     arrange(desc(high_mod_count)) %>%
     slice_head(n = 3) %>%
     select(Species.of.concern, high_mod_count)
   
   # create stats tables
   assign(paste0("top_taxa_high_mod_", name), top_taxa_high_mod, envir = .GlobalEnv)
   
   # apply colour coding
   threat_df_colored <- threat_df %>%
     mutate(across(
       .cols = -Species.of.concern,
       .fns = ~ cell_spec(" ",
                          background = sapply(.x, get_score_color_threat),
                          color = sapply(.x, get_score_color_threat),
                          bold = TRUE,
                          extra_css = "display: block; width: 100%; 
                          height: 100%; padding: 20px 0; margin: 0;")
     ))
   
   n_rows <- nrow(threat_df_colored)
   
   rows_to_color <- c(1,3,5,7,9,11)
   rows_to_color <- rows_to_color[rows_to_color <= n_rows]
   
   ############ Create and render the table #####################
   threat_table <- threat_df_colored %>%
     kbl(col.names = c("Target taxon",
                       "Human use of species",
                       "Agriculture/ Silviculture/ Ranching/ Grazing",
                       "Residential/ Development/ Mining/Roads",
                       "Habitat degradation",
                       "Habitat conversion",
                       "Tourism/ Recreation",
                       "Fire/ Disturbance regime modification",
                       "Invasive species impact",
                       "Climate change",
                       "Pests/ Pathogens",
                       "Extremely small restricted populations"),
         escape = FALSE) %>%
     kable_classic(full_width = FALSE) %>%
     kable_styling(font_size = 20) %>%
     column_spec(1, italic = TRUE) %>%
     row_spec(0, bold = TRUE)
   
   # loop as diff number of rows for each 
   for (r in rows_to_color) {
     threat_table <- row_spec(threat_table, r, background = "#F0F0F0")
   }
   
   print(threat_table)  # render in viewer
   
   # save
   assign(paste0("threat_", name, "_colored"), threat_df_colored, 
          envir = .GlobalEnv)
 }
 
 
 
 
 
 ##################
 # threat legend
 ###################
 
 threat_legend <- read.csv(file.path(analysis_dir,"threat_legend.csv"), 
                           header = TRUE)
 
 threat_legend <- threat_legend %>%
   mutate(across(where(is.character), ~ gsub("\\s+", " ", trimws(.))))
 
 # Define your color palette function once
 color_5 <- brewer.pal(9, "OrRd")[6]
 color_4 <- brewer.pal(9, "OrRd")[5]
 color_3 <- brewer.pal(9, "OrRd")[4]
 color_0 <- brewer.pal(9, "OrRd")[2]
 
 get_legend_colour <- function(score) {
   if (score == "high impact") {
     color_5  
   } else if (score == "moderate impact") {
     color_4
   } else if (score == "low impact") {
     color_3
   } else if (score == "not a significant threat") {
     color_0
   } else if (score == "cause of past decline") {
     "lightpink"  
   }
 }
 
 
 threat_legend <- threat_legend %>%
   mutate(colour = cell_spec(colour,
                             background = sapply(colour, get_legend_colour),
                             color = sapply(colour, get_legend_colour),
                             bold = TRUE,
                             extra_css = "display: block; width: 100%; 
                             height: 100%; padding: 20px 0; margin: 0;"))
 
 
 # make the table
 threat_legend %>%
   kbl(col.names = c("Key to threat impact matrices",
                     " "),
       escape = FALSE) %>%
   kable_classic(full_width = F) %>%
   kable_styling(font_size = 30) %>%
   row_spec(0, bold = TRUE)
 
 
 
 
 ##################
 # threat legend transposed
 ###################
 
 
 
 threat_legend <- read.csv(file.path(analysis_dir,"threat_legend.csv"), header = TRUE) %>%
   mutate(across(where(is.character), ~ gsub("\\s+", " ", trimws(.))))
 
 # Define colors
 color_5 <- brewer.pal(9, "OrRd")[6]
 color_4 <- brewer.pal(9, "OrRd")[5]
 color_3 <- brewer.pal(9, "OrRd")[4]
 color_0 <- brewer.pal(9, "OrRd")[2]
 
 get_legend_colour <- function(score) {
   if (score == "high impact") {
     color_5  
   } else if (score == "moderate impact") {
     color_4
   } else if (score == "low impact") {
     color_3
   } else if (score == "not a significant threat") {
     color_0
   } else if (score == "cause of past decline") {
     "lightpink"  
   }
 }
 
 # Apply colors before transpose
 threat_legend <- threat_legend %>%
   mutate(colour = cell_spec(colour,
                             background = sapply(colour, get_legend_colour),
                             color = sapply(colour, get_legend_colour),
                             bold = TRUE,
                             extra_css = "display: block; width: 100%; height: 100%; padding: 1px 0; margin: 10;"))
 
 # Transpose the data frame
 threat_legend_t <- as.data.frame(t(threat_legend))
 
 # Use the 2nd row as the column names (threat row)
 colnames(threat_legend_t) <- threat_legend_t[2, ]
 
 # Remove only the 2nd row (used as header), but keep the 1st row (colour) as a data row
 threat_legend_t <- threat_legend_t[-2, , drop = FALSE]
 
 
 # Print table with kable
 threat_legend_t %>%
   kbl(escape = FALSE) %>%
   kable_classic(full_width = F) %>%
   kable_styling(font_size = 30) %>%
   kable_styling() %>%
  
   row_spec(0, bold = TRUE)
