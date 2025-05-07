library(tidyverse)
library(terra)
library(sf)
source("R/__globals.R")
fnlc_org <- rast(g("{spatial_raster_location}/1_FNLC/fnlc_orig30m"))


pr <- as.polygons(fnlc_org > -Inf)
pr_sf <- st_as_sf(pr) |> 
  st_cast( "POLYGON") |>mutate(rn = st_area(geometry)) |> 
  slice_max(rn)

write_sf(pr_sf |> select(), (g("{spatial_raster_location}/1_FNLC/FNLC_boundary.shp")))





