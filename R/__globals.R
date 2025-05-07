## Global variables and pacakages

library(tidyverse)
library(sf)
library(stars)
library(terra)
library(fmesher)
ont.proj <- 3161
g <- glue::glue
date_compiled <- "2025-02-28"
run_a2 <- FALSE

source("R/__paths.R")

 

xml <- readr::read_lines(onlc_xml_path)
OLC_label_df <- 
  str_extract_all(xml, "\\<Label\\>.+\\<\\/Label\\>") |> 
  str_remove_all("\\<[\\/]*Label\\>") |> 
  tibble(label = _) |> 
  filter(label!="character(0)") |> 
  separate(label, into = c("value", "label"), sep = "\\s-\\s", extra = 'merge') |> 
  bind_rows(tibble(value = c("157", "247"), label =c("Unusable", "Cloud/Shadow")))

colour_pallette <- 
  read_delim(olcc_pallette, delim = " ",
             col_names = c("value", "R", "G", "B")) |> 
  rowwise() |> 
  mutate(colour = rgb(R, G, B, maxColorValue = 255))


get_project_files <- function(type){
  list.files(g(wt_db_loc),g("({type})"),
           recursive = T, full.names = T, ignore.case = T) %>% 
  .[str_detect(., "((Between)|(Birds_of_James_Bay_Lowlands)|(Shield-Lowlands)|(Winisk)|Extraction|Pilot)")]
}


## Set Spatial Extent --------------------------------
ont_buff <- BASSr::ontario |> st_buffer(10000)
### locations from OntarioAtlas3_Data_prep------------
ra_area_official <- read_sf(
  ra_approx_area_loc
) |> st_transform(ont.proj) 


ra_area <- 
  ra_area_official |> 
  st_intersection(ont_buff)

ra_buffer <- ra_area|> 
  st_buffer(50*1000)|> 
  st_intersection(ont_buff)



a_b <- ra_buffer|> st_area() |> units::set_units('km2')
a_a <- ra_area |> st_area() |> units::set_units('km2')



locs <- g("{rds_data_loc}/locations.rds") |> 
  read_rds() |> 
  filter(collection %in% c("WildTrax", "ONATLAS3PC", "CWS-ON-TAU Boreal Pilot Project") &
         str_detect(project, "Nocturnal|Resample", negate=T)) |> 
  st_filter(ra_buffer) 



nn <- nngeo::st_nn(locs, locs, k = 200, maxdist = 150, returnDist = T)
locs_neighbours <- map(1:nrow(locs), ~tibble(
  project = locs$project[[.x]],
  location = locs$location[[.x]],
                         n_neighbours = length(nn$nn[[.x]])-1,
                         neighbours = list(locs$location[nn$nn[[.x]][-1]]),
                         neighbour_proj = list(locs$project[nn$nn[[.x]][-1]]),
                         dists = list(nn$dist[[.x]][-1]))) |> 
  list_rbind() |> 
  unnest(c(neighbours, neighbour_proj, dists)) 

locs_neighbours |> filter(n_neighbours>0) |> janitor::tabyl(project)
locs_neighbours$neighbour_proj[[1]]

ll_ag <- 
locs |>  
  mutate(proj_loc = glue::glue("{project}__{location}__{collection}")) |> 
  dplyr::select(proj_loc) %>%
  bind_cols(st_coordinates(.))

locs_agg <- aggregate(ll_ag,ll_ag,
                      do_union = T, FUN = function(x) if(is.numeric(x)){mean(x)
                        }else{list(sort(unique(x)))}, 
                      join = function(x,y) st_is_within_distance(x,
                                                                 y, dist =10 ))

aggregated_locs <- 
locs_agg |> st_drop_geometry() |> 
  mutate(site_id_agg = as.numeric(factor(glue::glue("{X}_{Y}")))) |>
  distinct() |> 
  unnest(proj_loc) |> 
  separate(proj_loc, into = c("project", "location", "collection"), sep = "__") |> 
  st_as_sf(coords = c("X", "Y"), crs= ont.proj)


individual_locs <- distinct(aggregated_locs, site_id_agg, geometry)

  
spp_codes_all <- naturecounts::meta_species_taxonomy() |> 
  left_join(naturecounts::meta_species_codes() |> 
              filter(authority=="BSCDATA"),
            by = join_by(species_id)) |> 
  slice_min(rank, by = english_name, with_ties = F)


bound_hull <- fm_extensions(
  ra_buffer,
  convex = c(10000, 75000),
  concave = c(10000, 75000)
)


bound_ <- st_intersection(ra_buffer, BASSr::ontario)
mesh_inla_full <- fm_mesh_2d_inla(
  loc =  individual_locs,
  boundary = bound_hull,
  crs = fm_crs(st_as_sf(individual_locs)),
  # boundary = bnd,
  max.edge=c(30000,50000),# 50000
  min.angle=c(15,21 ),
  max.n=c(5000), ## Safeguard against large meshes.
  max.n.strict=c(5000 ), ## Don't build a huge mesh!
  cutoff=10000,#.0005, ## Filter away adjacent points.600
  # offset = c(100, 300))
  offset=c(-0.05,-0.1)) ## Offset for extra boundaries, if needed.



test_training_data <- read_rds(g("output/{date_compiled}_test_training_data.rds"))

train_locs <-  test_training_data$train_recordings |> 
  left_join(aggregated_locs,by = join_by(project, location, collection)) |> 
  st_as_sf() |> 
  dplyr::distinct(site_id_agg, geometry) |> 
  filter(!st_is_empty(geometry))

mesh_inla <- fm_mesh_2d_inla(
  loc =  train_locs,
  boundary = bound_hull,
  crs = fm_crs(st_as_sf(train_locs)),
  # boundary = bnd,
  max.edge=c(30000,50000),# 50000
  min.angle=c(15,21 ),
  max.n=c(5000), ## Safeguard against large meshes.
  max.n.strict=c(5000 ), ## Don't build a huge mesh!
  cutoff=10000,#.0005, ## Filter away adjacent points.600
  # offset = c(100, 300))
  offset=c(-0.05,-0.1)) ## Offset for extra boundaries, if needed.

rm(test_training_data)
rm(train_locs)

hull_d <- fm_extensions(
  st_as_sf(individual_locs),
  convex = c(20000, 50000),
  concave = c(20000, 50000)
)
mesh_il <- fm_mesh_2d_inla(
  loc =  individual_locs,
  boundary = hull_d,
  crs = fm_crs(st_as_sf(individual_locs)),
  # boundary = bnd,
  max.edge=c(20000,50000),# 50000
  min.angle=c(15,21 ),
  max.n=c(1000), ## Safeguard against large meshes.
  max.n.strict=c(1000 ), ## Don't build a huge mesh!
  cutoff=5000,#.0005, ## Filter away adjacent points.600
  # offset = c(100, 300))
  offset=c(-0.05,-0.1))  ## Offset for extra boundaries, if needed. |> ## Offset for extra boundaries, if needed.

mesh_data <-  mesh_il |> 
  fmesher::fm_as_sfc( format='int') |> 
  st_as_sf() 
# sf::st_cast("POLYGON") |> st_zm()


