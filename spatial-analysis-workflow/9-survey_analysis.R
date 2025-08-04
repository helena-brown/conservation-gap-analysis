################################################################################
# INSTITUTION QUESTION
################################################################################

library(dplyr)
library(tidyr)
library(ggplot2)

respondents <- read.csv(file.path(analysis_dir,"respondents_survey.csv"), 
                           header=T)

respondents <- respondents %>%
  filter(Responses > 0)

# create a color blind friendly palette 
cbPalette <- c("#0072B2","#E69F00", "#56B4E9", "#009E73", "#F0E442", "#D55E00"
               , "#CC79A7", "#999999")

# pie chart with legend 
ggplot(respondents, aes(x = "", y = Responses, fill = Institution_Type)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  scale_fill_manual(name = "Institution Type", values = cbPalette) +
  theme_void()


################################################################################
# TAXONOMY QUESTION
################################################################################


taxonomy <- read.csv(file.path(analysis_dir,"taxonomy_survey.csv"), 
                        header=T)

taxonomy <- taxonomy %>%
  filter(Responses > 0)

# create a color blind friendly palette 
cbPalette <- c("#0072B2","#E69F00", "#56B4E9", "#009E73", "#F0E442", "#D55E00"
               , "#CC79A7", "#999999")

# pie chart with legend 
ggplot(taxonomy, aes(x = "", y = Responses, fill = Taxonomy_Used)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  scale_fill_manual(name = "Taxonomic database used", values = cbPalette) +
  theme_void()



################################################################################
# SUMMARY OF TAXA CONSERVATION QUESTION
################################################################################

# install.packages('colorspace')
library(colorspace)

sum_consv <- read.csv(file.path(analysis_dir,
                                "taxa_conservation_summary_survey.csv"), 
                     header=T)

sum_consv <- sum_consv %>%
  filter(Responses > 0)

sum_consv$Target_taxa <- factor(sum_consv$Target_taxa)

# create named vector of labels with italic expressions
italic_labels <- setNames(
  paste0("italic('", levels(sum_consv$Target_taxa), "')"),
  levels(sum_consv$Target_taxa)
)

# pie chart with legend
ggplot(sum_consv, aes(x = "", y = Responses, fill = Target_taxa)) +
  geom_bar(stat = "identity", color = "white", size = 0.1, width = 1) +
  coord_polar("y", start = 0) +
  scale_fill_viridis_d(
    name = "Target taxa with conservation or threat \ninformation available from respondent institutions",
    option = "G",
    labels = parse(text = italic_labels)
  ) +
  theme_void()





################################################################################
# EACH TAXA CONSERVATION QUESTION
################################################################################

species_consv <- read.csv(file.path(analysis_dir,
                                    "species_conservation_actions_survey.csv"), 
                          header = TRUE)

# reformatting data
format_species_consv <- species_consv %>%
  pivot_longer(c("Collect.and.distribute.germplasm",
                 "Conservation.horticulture",
                 "Cryopreservation.and.or.micropropagation",
                 "Habitat.restoration",
                 "Implement.protection.policies.or.regulations",
                 "Occurrence.surveys.or.population.monitoring",
                 "Pollen.and.or.seed.banking",                  
                 "Population.reinforcement.or.introduction",
                 "Protect.and.or.manage.habitat",               
                 "Public.awareness.or.education",
                 "Research..Climate.Change",                    
                 "Research.Taxonomy",
                 "Research.Genetics",
                 "Research.Pests.and.Pathogens"),
               names_to = "Response",
               values_to = "Total")

# map original column names to clean labels and set factor levels
response_labels <- c(
  "Collect.and.distribute.germplasm" = "Collect and distribute germplasm",
  "Conservation.horticulture" = "Conservation horticulture",
  "Cryopreservation.and.or.micropropagation" = "Cryopreservation or micropropagation",
  "Habitat.restoration" = "Habitat restoration",
  "Implement.protection.policies.or.regulations" = "Implement protection policies or regulations",
  "Occurrence.surveys.or.population.monitoring" = "Occurrence surveys or population monitoring",
  "Pollen.and.or.seed.banking" = "Pollen or seed banking",
  "Population.reinforcement.or.introduction" = "Population reinforcement or introduction",
  "Protect.and.or.manage.habitat" = "Protect or manage habitat",
  "Public.awareness.or.education" = "Public awareness or education",
  "Research..Climate.Change" = "Research: Climate Change",
  "Research.Genetics" = "Research: Genetics",
  "Research.Pests.and.Pathogens" = "Research: Pests and Pathogens",
  "Research.Taxonomy" = "Research: Taxonomy"
)

format_species_consv$Response <- factor(response_labels[format_species_consv$Response],
                                        levels = response_labels)

# plot stacked bar with legend
ggplot(format_species_consv, aes(fill = Response, y = Total, x = Target_taxa)) + 
  geom_bar(position = "stack", stat = "identity", color = "white", size = 0.1,
           width = 0.8) +
  labs(y = "Number of Institutions", x = "Target taxon") +
  scale_fill_viridis_d(name = "Conservation activities carried out at respondent institutions", 
                       labels = levels(format_species_consv$Response),
                       option = "C") +
  theme(axis.line = element_line(colour = "black", size = 0.25),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = 'italic'))


################################################################################
# EACH TAXA CONSERVATION QUESTION
################################################################################

urgent_consv <- read.csv(file.path(analysis_dir,
                                   "most_urgent_conservation_action_survey.csv"), 
                         header = TRUE)

# reformatting data
format_urgent_consv <- urgent_consv %>%
  pivot_longer(c("Collect.and.distribute.germplasm",
                 "Conservation.horticulture",
                 "Cryopreservation.and.or.micropropagation",
                 "Habitat.restoration",
                 "Implement.protection.policies.or.regulations",
                 "Occurrence.surveys.or.population.monitoring",
                 "Pollen.or.seed.banking",
                 "Population.reinforcement.or.introduction",
                 "Protect.or.manage.habitat",
                 "Public.awareness.or.education",
                 "Research..Genetics",
                 "Research..Taxonomy",
                 "Research..Climate.Change",
                 "Research..Pests...Pathogens",
                 "Unknown",
                 "None..no.conservation.actions.currently.needed."),
               names_to = "Response",
               values_to = "Total")

# map original column names to cleaned labels and set factor levels
response_labels <- c(
  "Collect.and.distribute.germplasm" = "Collect and distribute germplasm",
  "Conservation.horticulture" = "Conservation horticulture",
  "Cryopreservation.and.or.micropropagation" = "Cryopreservation or micropropagation",
  "Habitat.restoration" = "Habitat restoration",
  "Implement.protection.policies.or.regulations" = "Implement protection policies or regulations",
  "Occurrence.surveys.or.population.monitoring" = "Occurrence surveys or population monitoring",
  "Pollen.or.seed.banking" = "Pollen or seed banking",
  "Population.reinforcement.or.introduction" = "Population reinforcement or introduction",
  "Protect.or.manage.habitat" = "Protect or manage habitat",
  "Public.awareness.or.education" = "Public awareness or education",
  "Research..Genetics" = "Research: Genetics",
  "Research..Taxonomy" = "Research: Taxonomy",
  "Research..Climate.Change" = "Research: Climate Change",
  "Research..Pests...Pathogens" = "Research: Pests and Pathogens",
  "Unknown" = "Unknown",
  "None..no.conservation.actions.currently.needed." = "None: no conservation actions currently needed."
)

format_urgent_consv$Response <- factor(response_labels[format_urgent_consv$Response],
                                       levels = response_labels)

# plot stacked bar with legend
ggplot(format_urgent_consv, aes(fill = Response, y = Total, x = Target_taxa)) + 
  geom_bar(position = "stack", stat = "identity", color = "white", size = 0.1,
           width = 0.8) +
  labs(y = "Number of Institutions", x = "Target taxon") +
  scale_fill_viridis_d(name = "Most urgent conservation action", 
                       labels = levels(format_urgent_consv$Response),
                       option = "C") + 
  theme(axis.line = element_line(colour = "black", size = 0.25),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = 'italic'))


################################################################################
# MAIN THREATS TO TAXA QUESTION
################################################################################

main_threats <- read.csv(file.path(analysis_dir,
                                   "main_threats_survey.csv"), 
                         header = TRUE)

library(dplyr)
library(tidyr)
library(ggplot2)


# reformatting data
format_main_threats <- main_threats %>%
  pivot_longer(c("Agriculture..silviculture..and.or.ranching",
                 "Climate.change",
                 "Development..mining..and.or.roads",
                 "Disturbance.regime.modification",
                 "Inbreeding.or.introgression",
                 "Invasive.species.competition",
                 "Pests.or.pathogens",
                 "Tourism.or.recreation",
                 "Wild.harvesting",
                 "Unknown"),
               names_to = "Response",
               values_to = "Total")

# map original column names to cleaned labels and set factor levels
response_labels <- c(
  "Agriculture..silviculture..and.or.ranching" = "Agriculture, silviculture and/or ranching",
  "Climate.change" = "Climate change",
  "Development..mining..and.or.roads" = "Development, mining and/or roads",
  "Disturbance.regime.modification" = "Disturbance regime modification",
  "Inbreeding.or.introgression" = "Inbreeding or introgression",
  "Invasive.species.competition" = "Invasive species competition",
  "Pests.or.pathogens" = "Pests or pathogens",
  "Tourism.or.recreation" = "Tourism or recreation",
  "Wild.harvesting" = "Wild harvesting",
  "Unknown" = "Unknown"
)

format_main_threats$Response <- factor(response_labels[format_main_threats$Response],
                                       levels = response_labels)


# plot stacked bar with legend
ggplot(format_main_threats, aes(fill = Response, y = Total, x = X)) + 
  geom_bar(position = "stack", stat = "identity", color = "white", size = 0.1,
           width = 0.8) +
  labs(y = "Number of Institutions", x = "Target taxon") +
  scale_fill_viridis_d(name = "Most significant threats to wild populations", 
                       labels = levels(format_main_threats$Response),
                       option = "D") + 
  theme(axis.line = element_line(colour = "black", size = 0.25),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = 'italic'))






###########################
# finding out which taxa have no data for all questions 
###########################

get_zero_taxa <- function(df) {
  df %>%
    group_by(across(1)) %>%  # group by the first column (taxon)
    summarise(total_sum = sum(Total, na.rm = TRUE)) %>%
    filter(total_sum == 0) %>%
    pull(1)  # return just the taxon names
}

# Apply to your datasets
zero_taxa_urgent <- get_zero_taxa(format_urgent_consv)
zero_taxa_species <- get_zero_taxa(format_species_consv)
zero_taxa_threats <- get_zero_taxa(format_main_threats)

# Get taxa that have zero across all three
taxa_all_zero <- Reduce(intersect, list(zero_taxa_urgent, zero_taxa_species, zero_taxa_threats))

print(taxa_all_zero)


###########################
# finding out which taxa have highest ans for all questions 
###########################

get_top_5_taxa_by_total <- function(df) {
  df %>%
    group_by(across(1)) %>%  # Group by first column (taxon)
    summarise(TotalSum = sum(Total, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(TotalSum)) %>%
    slice_head(n = 5)
}

top5_urgent   <- get_top_5_taxa_by_total(format_urgent_consv)
top5_species  <- get_top_5_taxa_by_total(format_species_consv)
top5_threats  <- get_top_5_taxa_by_total(format_main_threats)

cat(" Top 5 taxa in urgent_consv:\n"); print(top5_urgent)
cat("\n Top 5 taxa in species_consv:\n"); print(top5_species)
cat("\n Top 5 taxa in main_threats:\n"); print(top5_threats)



###########################
# finding out most common categories 
###########################


# Function to count the number of non-zero cells per column and return top N
get_top_responses_with_nonzero_counts <- function(df, top_n = 5) {
  df %>%
    filter(!is.na(Total) & Total != 0) %>%  # Keep only non-zero responses
    group_by(Response) %>%
    summarise(NonZeroCount = n(), .groups = "drop") %>%
    arrange(desc(NonZeroCount)) %>%
    slice_head(n = top_n)
}


top_n <- 5

cat(" Top", top_n, "most common non-zero responses:\n")

cat("\nurgent_consv:\n")
print(get_top_responses_with_nonzero_counts(format_urgent_consv, top_n))

cat("\nspecies_consv:\n")
print(get_top_responses_with_nonzero_counts(format_species_consv, top_n))

cat("\nmain_threats:\n")
print(get_top_responses_with_nonzero_counts(format_main_threats, top_n))


