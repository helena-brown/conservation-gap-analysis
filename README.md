# conservation-gap-analysis

In the spatial-analysis-workflow folder I have uploaded my edited files that I used to complete this analysis, and kept in the original files from the forked methodology for comparison

Below is the READ ME forked from the original methodology 

### BACKGROUND

This repository supports the Conservation Gap Analysis Methodology User Guide, which is an internal document of The Morton Arboretum, available to collaborators upon request; please email treeconservation@mortonarb.org to request access. Virtual training sessions for the methodology were held in June 2023, led by Emily Bruns and funded by The Morton Arboretum.

Developments to this methodology have been ongoing since 2016. Project partners have included Botanic Gardens Conservation International (BGCI), BGCI-US, the USDA Forest Service, the United States Botanic Garden, San Diego Botanic Garden, Missouri Botanic Garden, Fauna & Flora International, and Atlanta Botanic Garden. Thank you to all those who have contributed to developing the methodology, including but not limited to: Abby Meyer, Murphy Westwood, Sean Hoban, Silvia Alvarez-Clare, Amy Byrne, Daniel Carver, Dan Crowley, Audrey Denvir, David Gill, Kate Good, Colin Khoury, Jean Linsky, Gary Man, Nan McCarry, Ray Mims, David Pivorunas, Kirsty Shaw, Shannon Still, and Emily Warschefsky.

This work was funded by grants from the Institute of Museum and Library Services (awards MA-30-18-0273-18, MG-245575-OMS-20, and MG-251613-OMS-22), Magnolia Society International, Stanley Smith Horticultural Trust, United States Botanic Garden, National Science Foundation (award 1759759), and the USDA Forest Service (Cooperative Agreement 16-CA-11132546-045).

### CITATION

If you use code from this repository, please cite it in relevant reports/publications. Click "Cite this repository" in the right-hand panel to retrieve the citation.

### REPO OVERVIEW

The spatial-analysis-workflow folder contains a cohesive step-by-step workflow for completing the computational part of a gap analysis. There are additional optional analysis components outside R. The bonus-scripts folder holds scripts from various previous analyses, which are optional and do not flow as part of the workflow but can be used for inspiration as desired. Details are available in the Conservation Gap Analysis Methodology User Guide. Note that although the workflow is cohesive, it is not plug-and-play; most of the scripts are meant to be worked through slowly and considered and edited as-needed for the specific taxa, region, stakeholders, skillset, and goals of the specific project.

Briefly, the components of the gap analysis that can be completed using the spatial-analysis-workflow scripts are:
Script                               | Description
:----------------------------------- | :----------------------------------------------------------------------------------------------------------------------------------------------------
0-set_working_directory.R            | Set up main folder structure
1-get_taxa_metadata.R                | Compile metadata on target taxa: IUCN Red List category and native countries of occurrence from the IUCN Red List and BGCI GlobalTreeSeach
1-prep_gis_layers.R                  | Download and prepare spatial layers used throughout the workflow: world country boundaries, urban areas, ecoregions
2-compile_exsitu_data.R              | Compile and refine accession-level *ex situ* collections data from Genesys, FAO WIEWS, and (if applicable) a manually-executed survey of collections
3-get_occurrence_data.R              | Download wild occurrence data from six publicly-available databases: GBIF, iDigBio, IUCN Red List, SEINet Portal Network, BIEN, USFS FIA
4-compile_occurrence_data.R          | Compile and refine wild occurrence data downloaded in previous step, plus any manually-added data sources
5-flag_occurrence_points.R           | Flag potentially-suspect wild occurrence points
6-visualize_occurrence_points.R      | Create interactive HTML maps of wild occurrence points, for reviewing and retrieving IDs of points to remove
7-filter_occurrence_points.R         | Filter wild occurrence points using flags (script 5) and any manually-selected IDs to remove (script 6), producing a final set of points for analyses
8-calc_map_exsitu_coverage.R         | Calculate and map the geographic and ecological representation of *ex situ* collections for each target taxon, using a buffer-intersection method
8-calc_map_protected_areas.R         | Calculate and map the protected areas coverage for each target taxon, using a point-in-polygon method
8-map_taxon_richness.R               | Create maps of taxon richness at the country level using native countries of occurrence and at the state level using wild occurrence points
<img width=500/>                     | 
