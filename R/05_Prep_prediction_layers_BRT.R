library(terra)
library(tidyverse)
library(sf)
source("R/__globals.R")
ld <- list.dirs(
  spatial_raster_location, recursive = F
)

## List of all layers
ll <- map(ld, list.dirs, recursive=F) |>
  list_c()

l <- ll[str_detect(ll, 'info', negate = T)]

r <- rast(l[str_detect(l, "1km")])
r100 <- rast(l[str_detect(l, "200")])

name_replacements <- function(x){
  x |>
    str_remove("y$") |>
   
    str_replace( "nfi_", "NFIS_") |>
    str_replace( "NFIS_(d2s|sp)", "NFIS2_\\1") |>
    str_replace( "sfi_", "SCANFI_") |>
    str_replace( "clm", "Climate") |>
    str_replace( "ohn_([f|s])wat_", "\\1wat_1_") |>
    str_replace( "Climateb", "Climate2") |>
    # str_replace("(olcc)_(\\d)", "\\1_cls\\2") |>
    
    # str_replace("ohn_swat_", "swat_1_") |>  
    # str_replace("ohn_fwat_", "fwat_1_") |> 
    str_replace("ohn_", "OHN_") |>
    str_replace("fire_fire", "fire_rec30m")|> 
    str_replace("insct_", "insct_rec30m_") |> 
    str_replace("harv_", "harv_rec30m_")
}

n100 <- names(r100) |> 
  str_replace( "_200", "_100") |>
  name_replacements()



n500 <- names(r) |>
  str_replace( "_1km", "_500") |>
  name_replacements() 
  



r_pred <- c(r, r100) |>
  terra::crop(ra_buffer) |> terra::mask(ra_buffer)

terra::writeCDF(r_pred,glue::glue("{prediction_layer_loc}/{Sys.Date()}_prediction_rasters.nc"), overwrite=T)

xy <- terra::crds(
  r_pred$fnlc_1_1km,
  df=T, na.all=T, na.rm=T )
r_pred_df <- as_tibble(r_pred)
names(r_pred_df) <- c(n500, n100)

if(nrow(r_pred_df) == nrow(xy)){
preds <- bind_cols(r_pred_df, xy |> setNames(c("X", "Y"))) |>
  mutate(across(c(geob_100,
                  geob_500,
                  geoq_100,
                  geoq_500,
                  dem_asp8_100,
                  dem_asp8_500,
                  fire_rec30m_100,
                  fire_rec30m_500,
                  harv_rec30m_500,
                  harv_rec30m_100,
                  insct_rec30m_500,
                  insct_rec30m_100
                  ), factor),
    project = NA, location = NA, offset = 0,
    total_100=NA,
    total_500=NA,
    random_val=NA)
}
write_rds(preds, g("{prediction_layer_loc}/{Sys.Date()}_prediction_rasters_df.rds"))
write_rds(xy, g("{prediction_layer_loc}/{Sys.Date()}_prediction_rasters_xy.rds"))
