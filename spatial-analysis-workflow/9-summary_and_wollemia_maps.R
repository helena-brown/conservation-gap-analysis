################################################################################
# WOLLEMIA INSTITUTION MAPPING
################################################################################
my.packages <- c('ggplot2','rnaturalearth','rnaturalearthdata','dismo',
                 'sf','geodata','raster',
                 'ggspatial','terra','dismo','sp','rJava','tidyr',"ggtext")
lapply(my.packages, require, character.only=TRUE)
rm(my.packages)

# making lats and longs into spatial objects
woll_inst_map <- read.csv(file.path(standard_exsitu,
                                    "Wollemia_institutions_mapping.csv"), 
                                    header=T)

# colour coding the points for genetic diversity pre post 2022
institution_categories <- woll_inst_map %>%
  group_by(inst_short) %>%
  summarize(
    has_pre = any(year_category == "pre 2022", na.rm = TRUE),
    has_post = any(year_category == "post 2022", na.rm = TRUE)
  ) %>%
  mutate(color_group = case_when(
    has_pre & has_post ~ "#009E73",      # both → green
    !has_pre & has_post ~ "#E69F00",     # only post → orange
    has_pre & !has_post ~ "#CC79A7",     # only pre → pink
    !has_pre & !has_post ~ "#CC79A7"     # NA → pink
  ))

# Join this back to the original data
woll_inst_map <- woll_inst_map %>%
  left_join(institution_categories %>% select(inst_short, color_group), by = "inst_short")


# st_as_sf = standard data as spatial object
#turn the data into spatial points by using the coordinates
woll_points <- st_as_sf(woll_inst_map, coords = c("inst_long", "inst_lat"), 
                        crs = 4326, remove=FALSE)

# for getting coordinates to zoom in for map
extent(woll_points)

####### results #########
# class      : Extent 
# xmin       : -123.2511 
# xmax       : 174.9058 
# ymin       : -43.38987 
# ymax       : 64.1403 
########################

# plotting the map
world <- ne_countries(scale = "medium", returnclass = "sf")

# map with view of all points no scale bar 
ggplot(data = world) + 
  geom_sf(fill = "white", color = "gray60") +
  geom_sf(data = woll_points, aes(color = color_group)) +            
  coord_sf(xlim = c(-140, 190), ylim = c(-60, 85), expand = FALSE)  +
  scale_color_manual(
    values = c(
      "#009E73" = "#009E73",  # both pre & post 2022
      "#E69F00" = "#E69F00",  # post 2022 only
      "#CC79A7" = "#CC79A7"   # pre 2022 or NA
    ),
    labels = c(
      "#009E73" = "Individuals from Wollemi meta-collection and 
      previous Wollemia nobilis distributions",
      "#E69F00" = "Individuals from Wollemi meta-collection",
      "#CC79A7" = "Individuals from previous collections or unknown origin"
    )) +
  labs(color = expression("Origin of " * italic("Wollemia nobilis") * " ex situ collection")) +
  scale_x_continuous(name = "Longitude") +
  scale_y_continuous(name = "Latitude") +
  annotation_north_arrow(location = "tr", 
                         style = north_arrow_fancy_orienteering,
                         height = unit(1, "cm"),               
                         width = unit(1, "cm")) +                
  theme(
    panel.background = element_rect(fill = "lightblue"),
    panel.grid = element_blank()
  )




# map with view of australia
ggplot(data = world) + 
  geom_sf(fill = "white", color = "gray60") +
  geom_sf(data = woll_points, aes(color = color_group), size = 3) +            
  coord_sf(xlim = c(110, 180), ylim = c(-50, -10), expand = FALSE)  +
  scale_color_manual(
    values = c(
      "#009E73" = "#009E73",  # both pre & post 2022
      "#E69F00" = "#E69F00",  # post 2022 only
      "#CC79A7" = "#CC79A7"   # pre 2022 or NA
    ),
    labels = c(
      "#009E73" = "Individuals from Wollemi meta-collection and 
      previous Wollemia nobilis distributions",
      "#E69F00" = "Individuals from Wollemi meta-collection",
      "#CC79A7" = "Individuals from previous collections or unknown origin"
    )) +
  labs(color = expression("Origins of " * italic("Wollemia nobilis") * " ex situ collection")) +
  annotation_scale() + 
  scale_x_continuous(name = "Longitude") +
  scale_y_continuous(name = "Latitude") +
  annotation_north_arrow(location = "tr", 
                         style = north_arrow_fancy_orienteering,
                         height = unit(1, "cm"),               
                         width = unit(1, "cm")) +                
  theme(
    panel.background = element_rect(fill = "lightblue"),
    panel.grid = element_blank()
  )


# map with view of europe
ggplot(data = world) + 
  geom_sf(fill = "white", color = "gray60") +
  geom_sf(data = woll_points, aes(color = color_group), size = 3) +            
  coord_sf(xlim = c(-30, 35), ylim = c(70, 35), expand = FALSE)  +
  scale_color_manual(
    values = c(
      "#009E73" = "#009E73",  # both pre & post 2022
      "#E69F00" = "#E69F00",  # post 2022 only
      "#CC79A7" = "#CC79A7"   # pre 2022 or NA
    ),
    labels = c(
      "#009E73" = "Individuals from Wollemi meta-collection and 
      previous Wollemia nobilis distributions",
      "#E69F00" = "Individuals from Wollemi meta-collection",
      "#CC79A7" = "Individuals from previous collections or unknown origin"
    )) +
  labs(color = expression("Origins of " * italic("Wollemia nobilis") * " ex situ collection")) +
  annotation_scale() + 
  scale_x_continuous(name = "Longitude") +
  scale_y_continuous(name = "Latitude") +
  annotation_north_arrow(location = "tr", 
                         style = north_arrow_fancy_orienteering,
                         height = unit(1, "cm"),               
                         width = unit(1, "cm")) +                
  theme(
    panel.background = element_rect(fill = "lightblue"),
    panel.grid = element_blank()
  )








################################################################################
# ALL INSTITUTION MAPPING - where any ex situ data is held, coordinates or none
# - splitting by case study
################################################################################


# making lats and longs into spatial objects
all_inst_map <- read.csv(file.path(standard_exsitu,
                                    "All_institutions_mapping.csv"), 
                          header=T)


############### South America ############### 

south_america_taxa_inst_map <- all_inst_map %>%
  filter(taxon_name_accepted %in% c("Araucaria araucana", 
                                    "Araucaria angustifolia"))

south_america_points <- st_as_sf(south_america_taxa_inst_map, 
                                 coords = c("inst_long", "inst_lat"), 
                   crs = 4326, remove=FALSE)

# for getting coordinates to zoom in for map
extent(south_america_points)





############### SE ASIA ############### 


south_east_asia_taxa_inst_map <- all_inst_map %>%
  filter(taxon_name_accepted %in% c("Agathis borneensis", 
                                    "Agathis dammara",
                                    "Agathis flavescens",
                                    "Agathis kinabaluensis",
                                    "Agathis lenticula",
                                    "Agathis orbicula",
                                    "Agathis robusta subsp. nesophila"))

south_east_asia_points <- st_as_sf(south_east_asia_taxa_inst_map, 
                                 coords = c("inst_long", "inst_lat"), 
                                 crs = 4326, remove=FALSE)
extent(south_east_asia_points)




############### New Caledonia ############### 


new_caledonia_taxa_inst_map <- all_inst_map %>%
  filter(taxon_name_accepted %in% c("Agathis lanceolata",
                                    "Agathis montana",
                                    "Agathis moorei",
                                    "Agathis ovata",
                                    "Araucaria goroensis",
                                    "Araucaria humboldtensis",
                                    "Araucaria luxurians",
                                    "Araucaria montana",
                                    "Araucaria muelleri",
                                    "Araucaria nemorosa",
                                    "Araucaria rulei",
                                    "Araucaria schmidii",
                                    "Araucaria scopulorum"))

new_caledonia_points <- st_as_sf(new_caledonia_taxa_inst_map, 
                                   coords = c("inst_long", "inst_lat"), 
                                   crs = 4326, remove=FALSE)
extent(new_caledonia_points)


############### Pacific ############### 


pacific_taxa_inst_map <- all_inst_map %>%
  filter(taxon_name_accepted %in% c("Agathis australis",
                                    "Agathis macrophylla",
                                    "Araucaria heterophylla",
                                    "Wollemia nobilis"))

pacific_points <- st_as_sf(pacific_taxa_inst_map, 
                                      coords = c("inst_long", "inst_lat"), 
                                      crs = 4326, remove=FALSE)
extent(pacific_points)




############### PLOTTING THE MAP ############### 


# add columns then combine into one data frame

south_america_points$region <- "South America"
south_east_asia_points$region <- "South-East Asia"
new_caledonia_points$region <- "New Caledonia"
pacific_points$region <- "Pacific Islands, Australia and New Zealand"

all_inst_points <- rbind(
  south_america_points,
  south_east_asia_points,
  new_caledonia_points,
  pacific_points
)

# plotting the map
world <- ne_countries(scale = "medium", returnclass = "sf")

# map with view of all points no scale bar 
ggplot(data = world) + 
  geom_sf() +
  geom_sf(data = all_inst_points, aes(color = region)) +            
  coord_sf(xlim = c(-130, 180), ylim = c(-50, 70), expand = FALSE)  + 
  scale_x_continuous(name = "Longitude") +
  scale_y_continuous(name = "Latitude")+
  # annotation_scale() + 
  annotation_north_arrow(location = "tr", 
                         style = north_arrow_fancy_orienteering,
                         height = unit(1, "cm"),               
                         width = unit(1, "cm")) +                
  theme(
    panel.background = element_rect(fill = "lightblue"),
    panel.grid = element_blank())


####### make pie chart map ##########
library(scatterpie)

pie_data <- all_inst_points %>%
  st_drop_geometry() %>%
  group_by(inst_long, inst_lat, region) %>%
  summarise(count = n(), .groups = "drop") %>%
  pivot_wider(names_from = region, values_from = count, values_fill = 0)


# plot pie charts 
ggplot(data = world) +
  geom_sf(fill = "white", color = "gray60") +
  geom_scatterpie(data = pie_data, 
                  aes(x = inst_long, y = inst_lat, group = paste(inst_long, inst_lat)),
                  cols = setdiff(names(pie_data), c("inst_long", "inst_lat")),
                  pie_scale = 0.4) + 
  labs(fill = "Geographic area representation") +
  coord_sf(xlim = c(-180, 180), ylim = c(-60, 90), expand = FALSE) +
  scale_x_continuous(name = "Longitude") +
  scale_y_continuous(name = "Latitude") +
  annotation_north_arrow(location = "tr", 
                         style = north_arrow_fancy_orienteering,
                         height = unit(1, "cm"),               
                         width = unit(1, "cm")) +
  theme(
    panel.background = element_rect(fill = "lightblue"),
    panel.grid = element_blank()
  )




####################################################################################
####### CHAT EDITED CODE ########################################################
###################################################################################



####################################################################################
####### plot for Australia ########################################################
###################################################################################


# Load world map
world <- ne_countries(scale = "medium", returnclass = "sf")

# Prepare pie data
pie_data <- all_inst_points %>%
  st_drop_geometry() %>%
  group_by(inst_long, inst_lat, region) %>%
  summarise(count = n(), .groups = "drop") %>%
  pivot_wider(names_from = region, values_from = count, values_fill = 0)

# -----------------------------------------------
# Step 1: Classify points into "offset" and "not"
# -----------------------------------------------

# You can use bounding boxes for rough filtering:
pie_data <- pie_data %>%
  mutate(offset_group = case_when(
    # North America
    inst_long > -150 & inst_long < -50 & inst_lat > 10 & inst_lat < 80 ~ TRUE,
    # Europe
    inst_long > -20 & inst_long < 40 & inst_lat > 35 & inst_lat < 70 ~ TRUE,
    # Australia
      inst_long > 110 & inst_long < 155 & inst_lat > -45 & inst_lat < -10 ~ TRUE,
    # New Zealand
      inst_long > 165 & inst_long < 180 & inst_lat > -50 & inst_lat < -30 ~ TRUE,
    TRUE ~ FALSE
  ))

# ---------------------------------------------------
# Step 2: Add offsets only to the selected regions
# ---------------------------------------------------

pie_data <- pie_data %>%
  mutate(
    x_plot = ifelse(offset_group, inst_long + runif(n(), -10, 10), inst_long),
    y_plot = ifelse(offset_group, inst_lat + runif(n(), -10, 10), inst_lat)
  )

# Step 3: Prepare line data (only for offset points)
line_data <- pie_data %>%
  filter(offset_group) %>%
  select(inst_long, inst_lat, x_plot, y_plot)


ggplot(data = world) +
  geom_sf(fill = "white", color = "gray60") +
  
  # Optional: original institution points for clarity
  geom_point(data = pie_data, aes(x = inst_long, y = inst_lat), size = 0.95, color = "black") +
  
  # Connecting lines (only for offset pies)
  geom_segment(data = line_data,
               aes(x = inst_long, y = inst_lat,
                   xend = x_plot, yend = y_plot),
               color = "gray40", linetype = "dashed", linewidth = 0.4) +
  
  # Pie charts (offset or not depending on location)
  geom_scatterpie(data = pie_data,
                  aes(x = x_plot, y = y_plot, group = paste(inst_long, inst_lat)),
                  cols = setdiff(names(pie_data), c("inst_long", "inst_lat", "x_plot", "y_plot", "offset_group")),
                  pie_scale = 0.25) +
  
  coord_sf(xlim = c(110, 190), ylim = c(-60, -10), expand = FALSE) +
  labs(fill = "Geographic area representation") +
  annotation_north_arrow(location = "tr",
                         style = north_arrow_fancy_orienteering,
                         height = unit(1, "cm"),
                         width = unit(1, "cm")) +
  scale_x_continuous(name = "Longitude") +
  scale_y_continuous(name = "Latitude") +
  annotation_north_arrow(location = "tr", 
                         style = north_arrow_fancy_orienteering,
                         height = unit(1, "cm"),               
                         width = unit(1, "cm")) + 
  annotation_scale() + 
  theme(
    panel.background = element_rect(fill = "lightblue"),
    panel.grid = element_blank()
  )





####################################################################################
####### plot for North America #####################################################
###################################################################################


# Load world map
world <- ne_countries(scale = "medium", returnclass = "sf")

# Prepare pie data
pie_data <- all_inst_points %>%
  st_drop_geometry() %>%
  group_by(inst_long, inst_lat, region) %>%
  summarise(count = n(), .groups = "drop") %>%
  pivot_wider(names_from = region, values_from = count, values_fill = 0)

# -----------------------------------------------
# Step 1: Classify points into "offset" and "not"
# -----------------------------------------------

# You can use bounding boxes for rough filtering:
pie_data <- pie_data %>%
  mutate(offset_group = case_when(
    # North America
    inst_long > -150 & inst_long < -65 & inst_lat > 25 & inst_lat < 55 ~ TRUE,
    # Europe
    inst_long > -20 & inst_long < 40 & inst_lat > 35 & inst_lat < 70 ~ TRUE,
    # Australia
    inst_long > 110 & inst_long < 155 & inst_lat > -45 & inst_lat < -10 ~ TRUE,
    # New Zealand
    inst_long > 165 & inst_long < 180 & inst_lat > -50 & inst_lat < -30 ~ TRUE,
    TRUE ~ FALSE
  ))

# ---------------------------------------------------
# Step 2: Add offsets only to the selected regions
# ---------------------------------------------------

pie_data <- pie_data %>%
  mutate(
    x_plot = ifelse(offset_group, inst_long + runif(n(), -10, 10), inst_long),
    y_plot = ifelse(offset_group, inst_lat + runif(n(), -20, 10), inst_lat)
  )

# Step 3: Prepare line data (only for offset points)
line_data <- pie_data %>%
  filter(offset_group) %>%
  select(inst_long, inst_lat, x_plot, y_plot)


ggplot(data = world) +
  geom_sf(fill = "white", color = "gray60") +
  
  # Optional: original institution points for clarity
  geom_point(data = pie_data, aes(x = inst_long, y = inst_lat), size = 0.95, color = "black") +
  
  # Connecting lines (only for offset pies)
  geom_segment(data = line_data,
               aes(x = inst_long, y = inst_lat,
                   xend = x_plot, yend = y_plot),
               color = "gray40", linetype = "dashed", linewidth = 0.4) +
  
  # Pie charts (offset or not depending on location)
  geom_scatterpie(data = pie_data,
                  aes(x = x_plot, y = y_plot, group = paste(inst_long, inst_lat)),
                  cols = setdiff(names(pie_data), c("inst_long", "inst_lat", "x_plot", "y_plot", "offset_group")),
                  pie_scale = 0.2) +
  
  coord_sf(xlim = c(-30, -140), ylim = c(60, -40), expand = FALSE) +
  labs(fill = "Geographic area representation") +
  annotation_north_arrow(location = "tr",
                         style = north_arrow_fancy_orienteering,
                         height = unit(1, "cm"),
                         width = unit(1, "cm")) +
  scale_x_continuous(name = "Longitude") +
  scale_y_continuous(name = "Latitude") +
  annotation_north_arrow(location = "tr", 
                         style = north_arrow_fancy_orienteering,
                         height = unit(1, "cm"),               
                         width = unit(1, "cm")) + 
  annotation_scale() + 
  theme(
    panel.background = element_rect(fill = "lightblue"),
    panel.grid = element_blank()
  )




####################################################################################
####### plot for EUROPE #####################################################
###################################################################################


# Load world map
world <- ne_countries(scale = "medium", returnclass = "sf")

# Prepare pie data
pie_data <- all_inst_points %>%
  st_drop_geometry() %>%
  group_by(inst_long, inst_lat, region) %>%
  summarise(count = n(), .groups = "drop") %>%
  pivot_wider(names_from = region, values_from = count, values_fill = 0)

# -----------------------------------------------
# Step 1: Classify points into "offset" and "not"
# -----------------------------------------------

# You can use bounding boxes for rough filtering:
pie_data <- pie_data %>%
  mutate(offset_group = case_when(
    # North America
    inst_long > -150 & inst_long < -65 & inst_lat > 25 & inst_lat < 55 ~ TRUE,
    # Europe
    inst_long > -3.5 & inst_long < 18 & inst_lat > 35 & inst_lat < 62 ~ TRUE,
    # Australia
    inst_long > 110 & inst_long < 155 & inst_lat > -45 & inst_lat < -10 ~ TRUE,
    # New Zealand
    inst_long > 165 & inst_long < 180 & inst_lat > -50 & inst_lat < -30 ~ TRUE,
    TRUE ~ FALSE
  ))

# ---------------------------------------------------
# Step 2: Add offsets only to the selected regions
# ---------------------------------------------------

pie_data <- pie_data %>%
  mutate(
    x_plot = ifelse(offset_group, inst_long + runif(n(), -10, 10), inst_long),
    y_plot = ifelse(offset_group, inst_lat + runif(n(), -7, 7), inst_lat)
  )

# Step 3: Prepare line data (only for offset points)
line_data <- pie_data %>%
  filter(offset_group) %>%
  select(inst_long, inst_lat, x_plot, y_plot)

n_lines <- nrow(line_data)
n_half <- ceiling(n_lines / 2)

light_greys <- colorRampPalette(c("grey40", "grey60"))(n_half)
dark_greys  <- colorRampPalette(c("grey10", "grey33"))(n_lines - n_half)

# Interleave dark and light greys for contrast
alternating_greys <- as.vector(rbind(sample(dark_greys), sample(light_greys)))
line_data$line_color <- alternating_greys[1:n_lines]


ggplot(data = world) +
  geom_sf(fill = "white", color = "gray60") +
  
  # Optional: original institution points for clarity
  geom_point(data = pie_data, aes(x = inst_long, y = inst_lat), size = 0.95, color = "black") +
  
  # Connecting lines (only for offset pies)
  geom_segment(data = line_data,
               aes(x = inst_long, y = inst_lat,
                   xend = x_plot, yend = y_plot),
               color = line_data$line_color,
               linetype = "dashed", linewidth = 0.4) +
  
  # Pie charts (offset or not depending on location)
  geom_scatterpie(data = pie_data,
                  aes(x = x_plot, y = y_plot, group = paste(inst_long, inst_lat)),
                  cols = setdiff(names(pie_data), c("inst_long", "inst_lat", "x_plot", "y_plot", "offset_group")),
                  pie_scale = 0.09) +
  
  coord_sf(xlim = c(-30, 45), ylim = c(70, 32), expand = FALSE) +
  labs(fill = "Geographic area representation") +
  annotation_north_arrow(location = "tr",
                         style = north_arrow_fancy_orienteering,
                         height = unit(1, "cm"),
                         width = unit(1, "cm")) +
  scale_x_continuous(name = "Longitude") +
  scale_y_continuous(name = "Latitude") +
  annotation_north_arrow(location = "tr", 
                         style = north_arrow_fancy_orienteering,
                         height = unit(1, "cm"),               
                         width = unit(1, "cm")) + 
  annotation_scale() + 
  theme(
    panel.background = element_rect(fill = "lightblue"),
    panel.grid = element_blank()
  )



#########################
# alternative plot for europe code
#########################


world <- ne_countries(scale = "medium", returnclass = "sf")

# Prepare pie data
pie_data <- all_inst_points %>%
  st_drop_geometry() %>%
  group_by(inst_long, inst_lat, region) %>%
  summarise(count = n(), .groups = "drop") %>%
  pivot_wider(names_from = region, values_from = count, values_fill = 0)

# -----------------------------------------------
# Step 1: Classify points into "offset" and "not"
# -----------------------------------------------

pie_data <- pie_data %>%
  mutate(offset_group = case_when(
    # North America
    inst_long > -150 & inst_long < -65 & inst_lat > 25 & inst_lat < 55 ~ TRUE,
    # Europe
    inst_long > -3.5 & inst_long < 18 & inst_lat > 35 & inst_lat < 62 ~ TRUE,
    # Australia
    inst_long > 110 & inst_long < 155 & inst_lat > -45 & inst_lat < -10 ~ TRUE,
    # New Zealand
    inst_long > 165 & inst_long < 180 & inst_lat > -50 & inst_lat < -30 ~ TRUE,
    TRUE ~ FALSE
  ))

# ---------------------------------------------------
# Step 2: Custom offsets
# ---------------------------------------------------

pie_data <- pie_data %>%
  mutate(
    is_uk = inst_long > -10 & inst_long < 2 & inst_lat > 49 & inst_lat < 61,
    is_east_of_belgium = inst_long > 4.5
  )

pie_data <- pie_data %>%
  mutate(
    x_plot = case_when(
      offset_group & is_uk ~ inst_long + runif(n(), -15, 0),               # UK: west only
      offset_group & is_east_of_belgium ~ inst_long + runif(n(), 0, 10),   # East of Belgium: east only
      offset_group ~ inst_long + runif(n(), -10, 10),                      # Others: random
      TRUE ~ inst_long
    ),
    y_plot = case_when(
      offset_group & is_uk ~ inst_lat + runif(n(), 0, 7),                  # UK: north only
      offset_group ~ inst_lat + runif(n(), -7, 7),                         # Others: random
      TRUE ~ inst_lat
    )
  )

# ---------------------------------------------------
# Step 3: Prepare connecting lines
# ---------------------------------------------------

line_data <- pie_data %>%
  filter(offset_group) %>%
  select(inst_long, inst_lat, x_plot, y_plot)

n_lines <- nrow(line_data)
n_half <- ceiling(n_lines / 2)

light_greys <- colorRampPalette(c("grey40", "grey60"))(n_half)
dark_greys  <- colorRampPalette(c("grey10", "grey33"))(n_lines - n_half)

alternating_greys <- as.vector(rbind(sample(dark_greys), sample(light_greys)))
line_data$line_color <- alternating_greys[1:n_lines]

# ---------------------------------------------------
# Plot
# ---------------------------------------------------

ggplot(data = world) +
  geom_sf(fill = "white", color = "gray60") +
  
  geom_point(data = pie_data, aes(x = inst_long, y = inst_lat), size = 0.95, color = "black") +
  
  geom_segment(data = line_data,
               aes(x = inst_long, y = inst_lat,
                   xend = x_plot, yend = y_plot),
               color = line_data$line_color,
               linetype = "dashed", linewidth = 0.4) +
  
  geom_scatterpie(data = pie_data,
                  aes(x = x_plot, y = y_plot, group = paste(inst_long, inst_lat)),
                  cols = setdiff(names(pie_data), c("inst_long", "inst_lat", "x_plot", "y_plot", "offset_group", "is_uk", "is_east_of_belgium")),
                  pie_scale = 0.09) +
  
  coord_sf(xlim = c(-30, 45), ylim = c(70, 32), expand = FALSE) +
  labs(fill = "Geographic area representation") +
  annotation_north_arrow(location = "tr",
                         style = north_arrow_fancy_orienteering,
                         height = unit(1, "cm"),
                         width = unit(1, "cm")) +
  scale_x_continuous(name = "Longitude") +
  scale_y_continuous(name = "Latitude") +
  annotation_scale() + 
  theme(
    panel.background = element_rect(fill = "lightblue"),
    panel.grid = element_blank()
  )




#####################
# institution map stats
########################


# how many institutions have each taxa 
taxa_inst_summary <- all_inst_map %>%
  filter(taxon_name_accepted %in% c(
    "Araucaria araucana", "Araucaria angustifolia",
    "Agathis borneensis", "Agathis dammara", "Agathis flavescens",
    "Agathis kinabaluensis", "Agathis lenticula", "Agathis orbicula",
    "Agathis robusta subsp. nesophila", "Agathis lanceolata", 
    "Agathis montana", "Agathis moorei", "Agathis ovata",
    "Araucaria goroensis", "Araucaria humboldtensis", "Araucaria luxurians",
    "Araucaria montana", "Araucaria muelleri", "Araucaria nemorosa",
    "Araucaria rulei", "Araucaria schmidii", "Araucaria scopulorum",
    "Agathis australis", "Agathis macrophylla", "Araucaria heterophylla",
    "Wollemia nobilis"
  )) %>%
  group_by(taxon_name_accepted) %>%
  summarise(num_institutions = n_distinct(inst_short))


# which geographic area has most institutions representing them
region_inst_summary <- all_inst_points %>%
  st_drop_geometry() %>%
  group_by(region) %>%
  summarise(num_institutions = n_distinct(inst_short))


# number of taxa per institutions
target_taxa <- c(
  "Araucaria araucana", "Araucaria angustifolia",
  "Agathis borneensis", "Agathis dammara", "Agathis flavescens",
  "Agathis kinabaluensis", "Agathis lenticula", "Agathis orbicula",
  "Agathis robusta subsp. nesophila", "Agathis lanceolata", 
  "Agathis montana", "Agathis moorei", "Agathis ovata",
  "Araucaria goroensis", "Araucaria humboldtensis", "Araucaria luxurians",
  "Araucaria montana", "Araucaria muelleri", "Araucaria nemorosa",
  "Araucaria rulei", "Araucaria schmidii", "Araucaria scopulorum",
  "Agathis australis", "Agathis macrophylla", "Araucaria heterophylla",
  "Wollemia nobilis"
)

# Now filter using this
taxa_per_institution <- all_inst_map %>%
  filter(taxon_name_accepted %in% target_taxa) %>%
  group_by(inst_short) %>%
  summarise(num_taxa = n_distinct(taxon_name_accepted)) %>%
  arrange(desc(num_taxa))



# taxa at only one garden

target_taxa <- c("Araucaria goroensis", "Agathis flavescens", "Agathis orbicula")

institutions_with_taxa <- all_inst_map %>%
  filter(taxon_name_accepted %in% target_taxa) %>%
  select(inst_short, taxon_name_accepted) %>%
  distinct() %>%
  arrange(taxon_name_accepted, inst_short)


taxa_UTurkuBG <- all_inst_map %>%
  filter(inst_short == "UTurkuBG") %>%
  select(taxon_name_accepted) %>%
  distinct() %>%
  arrange(taxon_name_accepted)
