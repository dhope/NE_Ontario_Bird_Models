library(terra)
library(tidyverse)
library(sf)
source("R/__globals.R")
ld <- list.dirs(
  spatial_raster_location,
  recursive = F
)


## TODO NAflag() for layers where zero is actually NA -------------

fnlc_boundary <- read_sf(g(
  "{spatial_raster_location}/1_FNLC/FNLC_boundary.shp"
))

## List of all layers
ll <- map(ld, list.dirs, recursive = F) |>
  list_c()

l <- ll[str_detect(ll, 'info|Annual', negate = T)]

fnlc <- rast(l[str_detect(l, "fnlc_\\d")]) |>
  mask(fnlc_boundary)

na_flags_zero <-
  "(niced|pg1mm|sumn|age|nfi_cls|d2s|pcp|mint|maxt|grows|sprec|fire|geob|geoq|harv|insct)_\\d"


r <- rast(l[
  str_detect(l, "1km") &
    str_detect(l, "fnlc_\\d", negate = T) &
    str_detect(l, na_flags_zero, negate = T)
])
r100 <- rast(l[
  str_detect(l, "200") &
    str_detect(l, "fnlc_\\d", negate = T) &
    str_detect(l, na_flags_zero, negate = T)
])


r_naz <- rast(l[
  str_detect(l, "1km") &
    str_detect(l, "fnlc_\\d", negate = T) &
    str_detect(l, na_flags_zero, negate = F)
])
r100_naz <- rast(l[
  str_detect(l, "200") &
    str_detect(l, "fnlc_\\d", negate = T) &
    str_detect(l, na_flags_zero, negate = F)
])
NAflag(r_naz) <- NAflag(r100_naz) <- 0

name_replacements <- function(x) {
  x |>
    str_remove("y$") |>

    str_replace("nfi_", "NFIS_") |>
    str_replace("NFIS_(d2s|sp)", "NFIS2_\\1") |>
    str_replace("sfi_", "SCANFI_") |>
    str_replace("clm", "Climate") |>
    str_replace("ohn_([f|s])wat_", "\\1wat_1_") |>
    str_replace("Climateb", "Climate2") |>
    # str_replace("(olcc)_(\\d)", "\\1_cls\\2") |>

    # str_replace("ohn_swat_", "swat_1_") |>
    # str_replace("ohn_fwat_", "fwat_1_") |>
    str_replace("ohn_", "OHN_") |>
    str_replace("fire_fire", "fire_rec30m") |>
    str_replace("insct_", "insct_rec30m_") |>
    str_replace("harv_", "harv_rec30m_")
}

n100 <- c(names(r100), names(r100_naz)) |>
  str_replace("_200", "_100") |>
  name_replacements()


n500 <- c(names(r), names(r_naz)) |>
  str_replace("_1km", "_500") |>
  name_replacements() |>
  str_replace("rd_dens_500", "rd_dens_1km")
n_fnlc <- names(fnlc) |>
  str_replace("_200", "_100") |>
  str_replace("_1km", "_500") |>
  name_replacements()


names_replaced <- c("X", "Y", n500, n100, n_fnlc)

r_pred <- c(r, r_naz, r100, r100_naz, fnlc) |>
  terra::crop(ra_buffer) |>
  terra::mask(ra_buffer)

terra::writeCDF(
  r_pred,
  glue::glue("{prediction_layer_loc}/{Sys.Date()}_prediction_rasters.nc"),
  overwrite = T,
  split = T
)

# xy <- terra::crds(
#   r_pred,
#   df=T, na.all=T, na.rm=T )
r_pred_df <- terra::as.data.frame(r_pred, xy = TRUE, na.rm = NA)
if (length(names(r_pred_df)) == length(names_replaced)) {
  names(r_pred_df) <- names_replaced
} else {
  rlang::abort("Bad name replacements")
}
# missing_all <-
# r_pred_df |> dplyr::select(-x,-y) |>
#   as.matrix() |>
#   is.na() |>
#   matrixStats::rowSums2()
#
# which(missing_all == 221)
#

cov <- read_rds(g("output/rds/{date_compiled}_spatial_covariates_data.rds"))
factor_names <- cov |> select(where(is.factor)) |> names()

names(cov)[!names(cov) %in% names(r_pred_df)]


# if(nrow(r_pred_df) == nrow(xy)){
preds <- #bind_cols(
  r_pred_df |> #, xy |> setNames(c("X", "Y"))) |>
  mutate(
    across(any_of(factor_names), factor),
    project = NA,
    location = NA,
    QPAD_offset = 0,
    # total_100=NA,
    # total_500=NA,
    random_val = NA
  )

write_rds(
  preds,
  g("{prediction_layer_loc}/{Sys.Date()}_prediction_rasters_df.rds")
)
# write_rds(xy, g("{prediction_layer_loc}/{Sys.Date()}_prediction_rasters_xy.rds"))

## Add in distance to ocean
# https://www150.statcan.gc.ca/n1/pub/16-510-x/16-510-x2024006-eng.htm
ocean_points <- read_rds("output/oceans_sf_from_raster.rds")
d2_ocean <- dplyr::select(r_pred_df, X, Y) |>
  st_as_sf(coords = c("X", "Y"), crs = ont.proj) |>
  nngeo::st_nn(ocean_points, returnDist = T) |>
  pluck('dist') |>
  list_c()

write_rds(d2_ocean, g("{prediction_layer_loc}/distance_2_ocean_prediction.rds"))
