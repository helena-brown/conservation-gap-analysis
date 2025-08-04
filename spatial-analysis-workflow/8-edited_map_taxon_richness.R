################################################################################
# Load libraries
################################################################################

# load packages
my.packages <- c('tidyverse','sf','RColorBrewer','leaflet','rnaturalearth',
                 'countrycode','webshot')
# versions I used (in the order listed above): 2.0.0, 1.0-13, 1.1-3, 2.1.2, 0.3.3, 1.5.0
#install.packages (my.packages) #Turn on to install current versions
lapply(my.packages, require, character.only=TRUE)
rm(my.packages)

################################################################################
# Set working directory
################################################################################

# create folder for output data
data_out <- "taxon_richness_maps"
if(!dir.exists(file.path(analysis_dir,data_out)))
  dir.create(file.path(analysis_dir,data_out), 
             recursive=T)


################################################################################
# Create functions
################################################################################

# function to create table of taxon richness and join to polygon data, for mapping
richness.countries <- function(df,polygons){
  # see max number of country codes for one taxon
  count_codes <- sapply(df$all_native_dist_iso2,function(x) str_count(x, pattern = "; "))
  # create an array of separated country codes
  iso_a2 <- str_split_fixed(df$all_native_dist_iso2, "; ", n = (max(count_codes)+1))
  # sum to calculate richness
  richness <- as.data.frame(table(iso_a2))
  # merge polygons with taxon richness data
  merged <- sp::merge(polygons,richness)
  # return shapefile with taxon richness added
  return(merged)
}

# create leaflet map for country-level taxon richness
map.countries <- function(countries,title,legend_text,legend_labels,pal,
                          my_lat,my_long
                          #,my_zoom
                          ){
  map <- leaflet() %>%
    # add background for the map; see additional background options here:
    #		https://leaflet-extras.github.io/leaflet-providers/preview/
    addProviderTiles("CartoDB.PositronNoLabels") %>%
    # add countries to map
    addPolygons(data = countries,
                color = "grey", weight = 2, opacity = 1,
                fillColor = ~pal(countries$Freq), fillOpacity = 1) %>%
    # add legend with number of taxa represented by each color
    addLegend(values = countries$Freq,
              pal = pal, opacity = 1,
              title = legend_text,
              labFormat = function(type, cuts, p) {paste0(legend_labels)},
              position = "bottomright") %>%
    addControl(
      html = "<img src='north_arrow.png' style='width:25px;'>",
      position = "topright"
    )
    # add title
   # addControl(title, position = "topright") %>%
    # set center coordinate and zoom level of the output map
   # setView(lat = my_lat, lng = my_long, zoom = my_zoom)
  
  return(map)
}

# create leaflet map for state/province-level taxon richness
map.states <- function(state_richness,states_all,title,legend_text,
                       legend_labels,pal,centroids){
  map <- leaflet() %>%
    # add background for the map; see additional background options here:
    #		https://leaflet-extras.github.io/leaflet-providers/preview/
    addProviderTiles("CartoDB.PositronNoLabels") %>%
    # add state/province taxon richness to map
    addPolygons(data = state_richness,
                stroke = FALSE,
                fillColor = ~pal(state_richness$Freq), fillOpacity = 1) %>%
    # add all state/province borders
    addPolygons(data = states_all,
                color = "grey", weight = 2, opacity = 1,
                fillOpacity = 0) %>%
    # add labels for number of taxa per state
    addLabelOnlyMarkers(data = state_richness,
                        lng = ~longitude, lat = ~latitude, 
                        label = ~Freq,
                        labelOptions = 
                          labelOptions(noHide = TRUE, textOnly = TRUE,
                                       # can change the text size as needed!
                                       textsize = "12px", 
                                       style = list("font-weight" = "bold"))) %>%
    # add legend with number of taxa represented by each color
    addLegend(values = state_richness$Freq,
              pal = pal, opacity = 1,
              title = legend_text,
              labFormat = function(type, cuts, p) {paste0(legend_labels)},
              position = "bottomright") %>%
    addControl(
      html = "<img src='north_arrow.png' style='width:25px;'>",
      position = "topright"
    ) 
    # add title
   # addControl(title, position = "topright")
  
  return(map)
}


################################################################################
# 1. Country-level taxon richness
################################################################################

## note that these calculations are based on native countries for each target 
##  taxon, based on the native countries of occurrence pulled in script
##  1-get_taxa_metadata.R, which gets countries listed in IUCN Red List
##  assessments, BGCI's GlobalTreeSearch database, and any manual edits you add. 
##  Native countries are NOT pulled from the taxon's occurrence points; if you'd
##  like to do that, reference section 2 and create a new section as needed

# read in world countries layer from rnaturalearth package
world_polygons <- ne_countries(scale = 50, type = "countries", 
                               returnclass = "sf")

# read in target taxa list with native countries added in 1-get_taxa_metadata.R
taxon_list <- read.csv(file.path(taxa_dir,"target_taxa_with_synonyms.csv"), 
                       header=T, colClasses="character",na.strings=c("","NA"))
# keep just our accepted taxa
taxon_list <- taxon_list %>% filter(taxon_name_status == "accepted")

taxon_list$all_native_dist_iso2[taxon_list$taxon_name 
                                == "Agathis robusta subsp. nesophila"] <- "PG"
### calculate and map country-level taxon richness using functions we created above
# note you can do this for all taxa, or for subsets, as desired; below we create
#   maps for taxa assesses as Threatened on the IUCN Red List and endemic taxa;
#   change as you like!

####
## country-level richness for ALL taxa
####

## calculate richness
ctry_richness_all <- richness.countries(taxon_list,world_polygons)
## title text for map
# my_title <- "Example taxon richness map"
## text for legend title
my_legend_title <- paste0("Number of native","<br/>"," target taxa")
## color bins and labels
# can look at distribution of data to inform bins
hist(ctry_richness_all$Freq,breaks=20,xlim=c(0,50),ylim=c(0,25))
# BASED ON YOUR DATA, assign bin breaks and labels 
#   most color palettes work will with max of 9 bins (not counting Inf)
# you'll need to decide if you want to use natural breaks, even, exponential, etc.
bins <- c(1,2,3,4,5,6,7,8,10,Inf)
my_legend_labels <- c("1","2","3","4","5","6","7","8-9","10+")
# create color palette
#   see palette options by running display.brewer.all()
my_palette <- colorBin(palette = "OrRd", bins = bins,
                       domain = ctry_richness_all$Freq, 
                       reverse = F, na.color = "white")
## create map
map_richness_all <- map.countries(ctry_richness_all,
                                 # my_title,
                                 "",
                                  my_legend_title,
                                  my_legend_labels,my_palette,
                                 "",
                                  # center lat and long of map, plus zoom level.
                                  # these are important if you want to save
                                  #   automatically as a PNG (line 219-222), but 
                                  # YOU WILL NEED TO PLAY WITH THESE TO GET RIGHT
                                 # 48, -1, 2
                                 )
## view map
# map_richness_all
## save map as HTML (interactive)
htmlwidgets::saveWidget(map_richness_all,
                        file.path(analysis_dir,data_out,
                                  "country-level_taxon_richness_ALL.html"))
## optionally, save screenshot of map (PNG)
# webshot(file.path(analysis_dir,data_out,
                #  "country-level_taxon_richness_ALL.html"), 
      #  file.path(analysis_dir,data_out,
              #    "country-level_taxon_richness_ALL.png"))

