################################################################################
# Load libraries
################################################################################
# install.packages("spThin")
# load packages
# install.packages('dplyr')
my.packages <- c('tidyverse','textclean','terra','sf','spThin','rnaturalearth',
                 'leaflet','knitr','dplyr')
# versions I used (in the order listed above): 2.0.0, 0.9.3, 1.7-29, 1.0-13, 0.2.0, 0.3.3, 2.1.2
#install.packages (my.packages) #Turn on to install current versions
lapply(my.packages, require, character.only=TRUE)
rm(my.packages)

addLegend <- leaflet::addLegend

list.files("taxon_points_final")
################################################################################
# Create functions
################################################################################

# function to create ex situ coverage map, with ecoregions, buffers, and points
map.pa <- function(taxon,ctry_boundaries,protected_areas,pts,pts_in_pa){
  map <- leaflet() %>%
    ## background
    addProviderTiles(providers$CartoDB.PositronNoLabels) %>%
    ## taxon name bolded at top right
    addControl(
      html = paste0('<div style="
      #background: rgba(255,255,255,0.8); 
                          #   padding: 4px 8px; 
                          #   border-radius: 4px; 
                             font-style: italic;
                             font-weight: bold;
                            # box-shadow: 0 0 5px rgba(0,0,0,0.3); 
                             max-width: 150px; 
                             word-wrap: break-word;">
                    ', taxon, '
                  </div>'),
      position = "topright"
    ) %>%
    ## protected areas
    addPolygons(
      data = protected_areas,                   
      smoothFactor = 0.5,	weight = 1, opacity = 1, color = "#4b965a",
      fillOpacity = 0.3) %>%
    ## country boundaries
    addPolygons(
      data = ctry_boundaries, fillColor = "transparent",
      weight = 1.5, opacity = 0.5, color = "#969696") %>%
    ## occurrence points
    addCircleMarkers(
      data = pts, lng = ~decimalLongitude, lat = ~decimalLatitude,
      color = "#d6569c", fillOpacity = 1, stroke = F,  
      # you may want to change the radius
      radius = 4) %>%
    ## occurrence points in protected areas
    addCircleMarkers(
      data = pts_in_pa,  lng = ~decimalLongitude, lat = ~decimalLatitude,
      color = "#231185", fillOpacity = 1, stroke = F,
      # you may want to change the radius
      radius = 4) %>%
    ## add scale bar
    addScaleBar(position = "bottomright",
                options = scaleBarOptions(maxWidth = 150)) %>%
    ## add legend
    addLegend(labels =
                c("Protected areas (WDPA, 2025)",
                  "Taxon occurrence <b>within</b> protected areas",
                  "Taxon occurrence <b>outside</b> protected areas"),
              colors = c("#4b965a","#231185","#d6569c"),
              position = "bottomright", opacity = 1) %>% 
    ## adding north arrow
    addControl(
      html = "<img src='north_arrow.png' style='width:25px;'>", data_in,
      position = "topright"
    ) %>%
    ## set view (long and lat) and zoom level, for when map initially opens
     setView(0, 0, zoom = 2) 
}

# calculate percent of occurrence points in global protected areas (WDPA)
#		first thin points, then see which are in PAs, then calculate percent
 pts.in.pa <- function(pts,lat_col,long_col,taxon_col,
                       # thin.km,
                       my.proj,pa.list){
  # first thin points; thin.km gives the distance between thinned pts
 # thin_pts <-
  #  spThin::thin(
    #  loc.data=pts,
    #  lat.col=lat_col,
     # long.col=long_col,
    #  spec.col=taxon_col,
    #  thin.par=thin.km, #kilometers b/w pts
   #   reps=100, #right now we are just doing once, but should prob do more?
     # locs.thinned.list.return=T,
    #  write.files=F,
    #  write.log.file=F)
#  thinned_pts <- thin_pts[[1]]
  # count number of points after thinning
 # num_thinned <- nrow(thinned_pts)
  # make thinned points a terra spatial object
   
  thinned_pts <- pts
   
  thinned_spatial <- vect(cbind(thinned_pts$decimalLongitude, thinned_pts$decimalLatitude), 
                          crs = my.proj)
  # see which points are in protected areas polygons; we iterate through the 
  #   list of protected areas polygons
  pt_in_pa_list <- lapply(pa.list, function(x) extract(x,thinned_spatial))
  # create one object with all pts found in PAs
  pt_in_pa <- do.call(rbind, pt_in_pa_list)
  pt_in_pa <- pt_in_pa[!is.na(pt_in_pa$WDPAID),]
  # join to pt data to get the lat and long attached again
  thinned_pts$id.y <- 1:nrow(thinned_pts)
  pt_in_pa <- left_join(pt_in_pa,thinned_pts) 
  pt_in_pa <- dplyr::distinct(pt_in_pa, id.y, .keep_all = TRUE)
  # count number of points in PAs
  in_pa <- nrow(pt_in_pa)
  # create list of values to return; dataframe of stats & actual points in PAs
  num_thinned <- nrow(thinned_pts)
  pt_in_pa_stats <- data.frame(
    taxon = pts$taxon_name_accepted[1],
    num_pt_unthinned = nrow(pts),
    km_thin_value = 0,
    num_pt_after_thin = 0,
    num_pt_in_pa = in_pa,
    percent_pt_in_pa = round(((in_pa/num_thinned)*100),digits=2))
  return_list <- list(pt_in_pa_stats,pt_in_pa)
  return(return_list)
 }




# create buffers around points, using specified projection
create.buffers <- function(df,radius,pt_proj,buff_proj,boundary){
  # turn occurrence point data into a SpatVector
  spat_pts <- vect(df, geom=c("decimalLongitude", "decimalLatitude"),
                   crs=pt_proj)
  # reproject to specified projection
  proj_df <- project(spat_pts,buff_proj)
  # place buffer around each point, then dissolve into one polygon
  buffers <- buffer(proj_df,width=radius)
  buffers <- aggregate(buffers,dissolve = TRUE)
  # clip by boundary so they don't extend into the water
  boundary <- project(boundary,buff_proj)
  buffers_clip <- crop(buffers,boundary)
  # return buffer polygons
  return(buffers_clip)
}


################################################################################
# Set working directory
################################################################################

# assign folder where you have input occurrence point data
data_in <- "taxon_points_final"

## DECIDE if you'd like to make maps, or run calculations only
#   no maps = FALSE
#   yes maps = TRUE
make_maps <- TRUE

if(make_maps){
  # create folder for output maps
  maps_out <- "protected_areas_maps"
  if(!dir.exists(file.path(analysis_dir,maps_out)))
    dir.create(file.path(analysis_dir,maps_out), 
               recursive=T)
}

# create folder for protected area layer downloads
pa_dir <- "protected_areas_downloads"
if(!dir.exists(file.path(gis_dir,pa_dir)))
  dir.create(file.path(gis_dir,pa_dir), 
             recursive=T)


################################################################################
# Download and read in protected areas layers
################################################################################

# you can download the global protected areas layer but it's very large
#   https://www.protectedplanet.net/en/thematic-areas/wdpa?tab=WDPA
## we are going to download PA layers by country instead (which
##  are still very large), so we can just get them for countries that have
##  target taxa...

# CHOOSE target countries (countries with target taxa), using the 
#   3-digit ISO code --> https://www.iban.com/country-codes
#   Note that Puerto Rico (PRI) is separate from the US
target_countries <- c("ARG","CHL","BRA","PRY","NCL","VUT","IDN","MYS","PHL")

###################################################
# NZL, SLB, BRN, NFK downloaded manually as wasn't working - maybe URL being built wrong 
manual_countries <- c("NZL","SLB","BRN","NFK","PNG")

# SET current month, using 3-letter abbreviation
current_month <- "Jul"
# SET current year
current_year <- "2025"

# set timeout value for when download.file will stop; default is 60 seconds
#   and this is too short for some of the larger protected areas files; we will
#   change to 600, but you can increase further if needed
getOption('timeout') # check what you have right now
options(timeout=600000) # set to new value

# download protected areas layers
pa_list <- list()
for(i in 1:length(target_countries)){
  
  # set up file paths
  shp_path <- paste0("https://d1gam3xoknrgr2.cloudfront.net/current/WDPA_WDOECM_",
                     current_month,current_year,"_Public_",target_countries[i],
                     "_shp.zip")
  dest_dir <- file.path(gis_dir,pa_dir)
  
  print(paste0("---Getting protected areas for ",target_countries[i],"---"))
  
  # download from web
  print("Downloading")
  if(grepl('http',shp_path)){
    dest <- paste0(dest_dir,"/",basename(shp_path))
    utils::download.file(shp_path, destfile = dest)
  }
  
  # unzip file downloaded
  print("Unzipping main download")
  if(grepl('.tgz$|.tar.gz$',dest)){
    utils::untar(dest, exdir = dest_dir)
  } else if(grepl('.zip$',dest)){
    utils::unzip(dest, exdir = dest_dir)
  } else{
    stop('Unsupported filetype')
  }
  
  # each protected areas download is split into three shapefiles;  
  #   we now unzip and read in each
  n <- 0
  shapes <- list()
  while(n < 3){
    path_now <- paste0(gsub(".zip","",dest),"_",n,".zip")
    dest_dir_now <- gsub(".zip","",path_now)
    dir.create(dest_dir_now)
    print(paste0("Unzipping part ",n))
    if(grepl('.tgz$|.tar.gz$',path_now)){
      utils::untar(path_now, exdir = dest_dir_now)
    } else if(grepl('.zip$',path_now)){
      utils::unzip(path_now, exdir = dest_dir_now)
    } else{
      stop('Unsupported filetype')
    }
    print(paste0("Reading in part ",n))
    shp_path <- paste0(gsub(".zip","",dest_dir_now),"/",
                       gsub(".zip","",basename(dest)),"-polygons.shp")
    temp_shape <- vect(shp_path)
    shapes[[n+1]] <- temp_shape
    n <- n + 1
  }
  
  # combine all shapefiles
  print("Combining into one shapefile")
  one_file <- Reduce(rbind,shapes)
  
  # add as item to list of country-level PAs
  pa_list[[i]] <- one_file
  
}

################################################### chat code edited from this code #####
# Now read in the manually downloaded shapefiles 

for(j in 1:length(manual_countries)){
  country_code <- manual_countries[j]
  print(paste0("---Reading manually downloaded protected areas for ",country_code,"---"))
  
  shapes <- list()
  n <- 0
  while(n < 3){
    # build path to folder: e.g. WDPA_WDOECM_Jul2025_Public_NZL_shp_0
    folder_now <- file.path(gis_dir,"protected_areas_downloads",
                            paste0("WDPA_WDOECM_",current_month,current_year,
                                            "_Public_",country_code,"_shp_",n))
    
    print(paste0("Reading in part ",n))
    # shapefile name: WDPA_WDOECM_Jul2025_Public_NZL_shp-polygons.shp
    shp_path_manual <- file.path(folder_now, paste0("WDPA_WDOECM_",current_month,current_year,
                                                    "_Public_",country_code,"_shp-polygons.shp"))
    
    temp_shape_manual <- vect(shp_path_manual)
    shapes[[n+1]] <- temp_shape_manual
    
    n <- n + 1
  }
  
  print("Combining into one shapefile")
  one_file_manual <- Reduce(rbind, shapes)
  
  pa_list[[length(pa_list)+1]] <- one_file_manual
}
################################################### chat code edited from this code #####

if(make_maps){
  # merge all PAs
  # we only need this for the map; for calculations we iterate instead
  # these functions can take a while...
  all_pa <- Reduce(rbind,pa_list)
}

rm(shapes,one_file,temp_shape,temp_shape_manual,i,j,n,
   country_code,shp_path_manual)

# if you'd like, you can remove local downloads:
#unlink(file.path(main_dir,gis_dir,pa_dir), recursive = TRUE)



################################################################################
# Set up standards and read in additional polygon data if mapping
################################################################################

# define projections
#	points will be WGS84
pt.proj <- "+proj=longlat +datum=WGS84"
# for calculations, we need something with units in meters and equal area
calc.proj <- "+proj=eqearth +datum=WGS84"

# read in target taxa list
# you can also simply create the "target_taxa" and "target_files" lists
#   manually, if desired
taxon_list <- read.csv(file.path(taxa_dir,"target_taxa_with_synonyms.csv"), 
                       header=T, colClasses="character",na.strings=c("","NA"))
# list of accepted taxa to cycle through
taxon_list <- read.csv(file.path(taxa_dir,"target_taxa_with_synonyms.csv"), 
                       header=TRUE, colClasses="character", na.strings=c("","NA"))

# get unique accepted taxa names as character vector
target_taxa <- unique(taxon_list$taxon_name_accepted)

# remove unwanted taxa
# excluded_taxa <- c("Agathis flavescens", "Agathis kinabaluensis",
             #      "Agathis lenticula","Agathis orbicula",
             #      "Agathis robusta subsp. nesophila",
             #      "Araucaria goroensis","Araucaria nemorosa","Wollemia nobilis")

# target_taxa <- target_taxa[!target_taxa %in% excluded_taxa]
# make taxon list with underscores added where spaces, to format for reading/
#   writing when we cycle through in our loop below
target_files <- unique(mgsub(target_taxa, 
                             c(" ","var.","subsp."), c("_","var","subsp")))

if(make_maps){
  
  # if you have very widespread taxa (e.g. span one third or more of the US)
  #   make a list of those widespread taxa to skip mapping (output too large)
  # taxa_no_map <- c("Asimina parviflora","Asimina triloba","Juglans cinerea",
                #   "Juglans nigra")
  
  taxa_no_map <- character(0)
  
  # read in world countries layer created in 1-prep_gis_layers.R
  # this will be used to clip buffers so they're not in the water, and
  #   also for adding to the map as a layer
  world_poly_clip <- vect(file.path(gis_dir,"world_countries_10m",
                                    "world_countries_10m.shp"))
  world_poly_sf <- st_as_sf(world_poly_clip)
  
}


################################################################################
## Calculate & map protected areas coverage
################################################################################

# start summary table for analysis results
# we add each target taxon as we go along
summary_tbl <- data.frame(
  taxon = "start",
  num_pt_unthinned = "start",
  km_thin_value = "start",
  num_pt_after_thin = "start",
  num_pt_in_pa = "start",
  percent_pt_in_pa = "start",
  stringsAsFactors=F)

### CYCLE THROUGH TARGET TAXA TO CALCULATE PROTECTED AREAS COVERAGE

for(i in 1:length(target_taxa)){
  
  ## can test with one taxon first if you'd like - skip loop line above and
  ##  uncomment next line
  # i <- 1
  
  # print progress
  cat("\nStarting", target_taxa[i], "\n")
  
  ### READ IN AND PREP POINT DATA
  
  ## read in occurrence points (includes ex situ)
  occ_pts <- read.csv(file.path(data_in,
                                paste0(target_files[i],".csv")), 
                      na.strings=c("","NA"), stringsAsFactors = F)
  nrow(occ_pts)
  
  ### CALC PROTECTED AREAS COVERAGE
  
  # if you have a lot of points (~10K+) you may get "Error: vector memory 
  #   exhausted"; I followed the instructions here to fix (MacOS)...
  #   https://stackoverflow.com/questions/51295402/r-on-macos-error-vector-memory-exhausted-limit-reached
 # pa_analysis_output <- pts.in.pa(occ_pts, "decimalLatitude", "decimalLongitude",
                              #    "taxon_name_accepted", 1, pt.proj, pa_list)
  
  pa_analysis_output <- pts.in.pa(
    occ_pts,
    lat_col = "decimalLatitude",
    long_col = "decimalLongitude",
    taxon_col = "taxon_name_accepted",
    my.proj = pt.proj,
    pa.list = pa_list
  )
  
  # add row to summary table
  summary_tbl[i,] <- pa_analysis_output[[1]]
  
  ### CREATE MAP
  
  if(make_maps){
    if(!(target_taxa[i] %in% taxa_no_map)){
      
      # skip Agathis australis
     # if (target_taxa[i] == "Agathis australis") {
      #  cat("\nSkipping", target_taxa[i], "\n")
      #  next
    #  }
      
      # points in protected areas df
      pts_in_pa_now <- pa_analysis_output[[2]]
      # add buffer around occurrence points and crop PAs to that buffer
      occ_buff <- create.buffers(occ_pts,100000,pt.proj,pt.proj,world_poly_clip)
      # First crop PA polygons to buffer area
      pa_crop <- crop(all_pa, occ_buff)
      
      # Then clip marine parts using country land polygons
      # First make sure both are in the same CRS
      pa_crop <- project(pa_crop, crs(world_poly_clip))
      land_only <- intersect(pa_crop, world_poly_clip)
      
      # if pa_crop is not empty, aggregate and convert; else create empty sf object
      if (nrow(pa_crop) > 0) {
        pa_crop <- aggregate(pa_crop, by = "ISO3")
        pa_crop <- st_as_sf(pa_crop)
      } else {
        cat("\nNo protected areas overlap buffer for ", target_taxa[i], 
            " — creating map without protected areas.\n")
        
        # create an empty sf object with same CRS as world_poly_sf
        pa_crop <- st_sf(geometry = st_sfc(), crs = st_crs(world_poly_sf))
      }
      
      # create and save the map anyway
      map <- map.pa(target_taxa[i], world_poly_sf, pa_crop, occ_pts, pts_in_pa_now); map
      
      # save map
      htmlwidgets::saveWidget(map,file.path(analysis_dir,maps_out,
                                            paste0(target_files[i],
                                                   "__mac_notthinned_protected_area_coverage_map",
                                                   ".html")))
      }
    }
  }


## write summary table
# summary_tbl
# write.csv(summary_tbl, file.path(analysis_dir,
                             #    paste0("_not_thinned_protected_area_coverage__",Sys.Date(),".csv")),
        #  row.names = F)

