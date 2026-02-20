source("R/__globals.R")
full_locs <- locs
locs <- individual_locs
locs_atlas2 <-
  g("{rds_data_loc}/locations.rds") |>
  read_rds() |>
  filter(
    collection %in%
      c("OBBA2PC") &
      str_detect(project, "Nocturnal|Resample", negate = T)
  ) |>
  st_filter(ra_buffer)
if (run_a2) {
  locs <- locs_atlas2
}


ld <- list.dirs(
  spatial_raster_location,
  recursive = F
)

## TODO NAflag() for layers where zero is actually NA -------------

library(terra)
## List of all layers
ll <- map(ld, list.dirs, recursive = F) |>
  list_c()

length(ll[str_detect(ll, "info$", negate = T)])


ll_no_proc_cls <- str_subset(ll, "_200$|1km$|_cls\\d+$|info$", negate = T)

ll_no_proc_cls |>
  str_extract("(?<=/)\\d+(?=_)") |>
  as.numeric() |>
  janitor::tabyl() |>
  arrange(1)


# OLCC --------------------------------------------------------------------
# Ontario Land Cover Compilation Version 2.0 (OMNRF 2014b)
# NOT Ontario Land Cover Version 1.0 (OMNR 2025)

olcc <- terra::rast(ll[
  str_detect(ll, "OLCC", negate = F) &
    str_detect(ll, "_200|_1km|info", negate = T)
])
olcc_extract <- function(olcc, dist_) {
  exactextractr::exact_extract(
    olcc,
    st_transform(
      st_buffer(locs, dist = dist_),
      st_crs(olcc)
    ),
    'sum'
  ) |>
    rowwise() |>
    mutate(total = sum(c_across(starts_with('sum.')))) |>
    ungroup() |>
    mutate(across(
      starts_with('sum.'),
      ~ {
        .x / total * 100
      }
    )) |>
    rename_with(.fn = ~ str_remove(.x, "sum."), starts_with('sum.')) |>
    rename_with(
      .fn = ~ glue::glue("{str_remove(.x, 'cls')}_{dist_}"),
      everything()
    )
}

olcc_500 <- olcc_extract(olcc, 500)
olcc_100 <- olcc_extract(olcc, 100)


cat_layers <- terra::rast(
  str_subset(ll_no_proc_cls, "lc|geo|wat_", negate = F) |>
    str_subset("asp8|geob|geoq|Quat_", negate = T)
)


cat_layers_extract <- function(layers, dist_) {
  map(
    names(layers),
    ~ {
      exactextractr::exact_extract(
        layers[[.x]],
        st_transform(
          st_buffer(locs, dist = dist_),
          st_crs(layers)
        ),
        'frac',
        force_df = T,
        default_value = -69
      ) |>
        mutate(across(
          contains("frac_"),
          ~ {
            .x * 100
          }
        )) |>
        rename_with(
          str_replace,
          pattern = "frac_(\\d+)",
          replacement = glue::glue("{.x}_\\1_{dist_}")
        ) |>
        as_tibble()
    }
  ) |>
    list_cbind() |>
    rename_with(str_remove, pattern = "_orig30m")
}

cat_layers_100 <- cat_layers_extract(cat_layers, 100)
cat_layers_500 <- cat_layers_extract(cat_layers, 500)


mod_layers <- str_subset(ll_no_proc_cls, "asp8|geob|geoq|/sp_") |> rast()

mod_layers_extract <- function(layers, dist_, f = 'mode') {
  map(
    names(layers),
    ~ {
      n <- glue::glue("{.x}_{dist_}")
      tibble(
        {{ n }} := factor(exactextractr::exact_extract(
          layers[[.x]],
          st_transform(
            st_buffer(locs, dist = dist_),
            st_crs(layers)
          ),
          f
        ))
      )
    }
  ) |>
    list_cbind() |>
    rename_with(str_remove, pattern = "_orig30m") |>
    rename_with(~ str_replace(.x, "sp_", "NFIS2_sp_"))
}


mod_ex_500 <- mod_layers_extract(mod_layers, 500)
mod_ex_100 <- mod_layers_extract(mod_layers, 100)

fire_500 <-
  mod_layers_extract(
    str_subset(ll_no_proc_cls, "fire|insct|harv") |>
      str_subset("1km|200", negate = T) |>
      rast(),
    500,
    'max'
  )
fire_100 <- mod_layers_extract(
  str_subset(ll_no_proc_cls, "fire|insct|harv") |>
    str_subset("1km|200", negate = T) |>
    rast(),
  100,
  'max'
)


fire_100$fire_rec30m_100 |> janitor::tabyl()
fire_100$insct_rec30m_100 |> janitor::tabyl()
fire_100$harv_rec30m_100 |> janitor::tabyl()


non_lc <- ll_no_proc_cls |>
  str_subset("lc|geo|/sp_|wat_|fire|insct|harv|asp|Annual", negate = T) #|>
# str_subset("lc|info|orig|1km|_200|asp8|wat_|fire|insct|harv|dem", negate = T)

non_mean <- non_lc |>
  str_subset("watd")

# non_lc <- ll[(str_detect(ll, "orig", negate = F) &
#                 +                 str_detect(ll, "lc|geo|/sp_|wat_|fire|insct|harv", negate = T))|
#                +                (str_detect(ll, "lc|info|orig|1km|_200|asp8|wat_|fire|insct|harv", negate = T) )
#              +              ]

get_non_lc <- function(r, dist_, f = 'mean',ii=NULL) {
  if(is_null(ii)){ii <- 1:nrow(locs)}
  map(
    r,
    ~ {
      r <- rast(.x)
      if (str_detect(.x, "age|SCANFI_cls|d2s")) {
        NAflag(r) <- 0
      }
      print(r)
      if (str_detect(.x, "dem")) {
        n <- str_extract(.x, "(dem_\\w+)")
      } else {
        n <- str_extract(
          .x,
          "((?<=DataLayers/{1,2}\\d{1,2}_)[\\w,\\W]+/\\w+(?=_(orig|rec)30m))"
        ) |>
          str_replace("/", "_")
      }
      n <- paste0(c(n, dist_), collapse = "_")
      exactextractr::exact_extract(
        r,
        st_transform(
          st_buffer(locs[ii,], dist = dist_),
          st_crs(r)
        ),
        f
      ) %>%
        tibble::tibble(
          {{ n }} := .
        ) |>
        dplyr::select(2)
    }
  ) |>
    bind_cols()
}

non_lc_r_m <- get_non_lc(str_subset(non_lc, "watd", negate = T), 100)
non_lc_r_min <- get_non_lc(non_mean, 100, f = 'min')
# cls_500 <- get_non_lc(str_subset(non_lc, "cls"), f = 'mode', 500)
# non_lc_r_500 <- get_non_lc(non_lc, 500)
non_lc_r_m5 <- get_non_lc(str_subset(non_lc, "watd", negate = T), 500)
non_lc_r_min_500 <- get_non_lc(non_mean, 500, f = 'min')
# cls_comp <- tibble(mean =non_lc_r_m5$NFIS_cls_500, mode = cls_500$NFIS_cls_500) |> mutate( rn=row_number())

# filter(cls_comp, mode ==0) |> slice_max(mean)

# To fix ---------
# dem_pslp # perspective (direction)

road_comb <- read_sf(
 road_layr_loc
) |>
  st_zm()
int_100 <- st_intersects(st_buffer(locs, 100), road_comb, sparse = F) |>
  matrixStats::rowSums2()
get_dens_road <- function(.x) {
  i <- st_intersection(.x, road_comb)

  if (length(i) == 0) {
    return(units::as_units(0, "m/ha"))
  }
  sum(st_length(i)) / units::set_units(st_area(.x), 'ha')
}

den_roads_500 <- st_buffer(locs, 500) |>
  rowwise() |>
  mutate(road_dens = get_dens_road(geometry)) |>
  ungroup()


cov <- bind_cols(
  cat_layers_100,
  olcc_100,
  non_lc_r_m,
  non_lc_r_min,
  cat_layers_500,
  olcc_500,
  non_lc_r_m5,
  non_lc_r_min_500,
  mod_ex_100,
  mod_ex_500,
  fire_100,
  fire_500
) |>
  dplyr::select(-starts_with('frac_'), -starts_with("total_\\d"))
cov$rd_exist_100 <- ifelse(int_100 == 0, 0, 1)
cov$rd_dens_1km <- den_roads_500$road_dens |> units::drop_units()

which_zeros <- function(string) {
  dplyr::select(cov, starts_with(string)) |>
    as.matrix() |>
    matrixStats::rowSums2() |>
    (\(x) x == 0)() |>
    which()
}
sw_strings <- c("fnlc", "clc20", "olcc", "olcb")
sw_missing <- map(sw_strings, which_zeros) |> setNames(nm = sw_strings)

for (ii in length(sw_strings)) {
  cov[
    sw_missing[[sw_strings[[ii]]]],
    str_subset(names(cov), sw_strings[[ii]])
  ] <- NA
}


if (run_a2) {
  write_rds(cov, g("output/rds/{Sys.Date()}_spatial_covariates_Atlas2.rds"))
} else {
  write_rds(cov, g("output/rds/{Sys.Date()}_spatial_covariates_data.rds"))
}


ocean_sf <- "output/oceans_sf_from_raster.rds"
if (!file.exists(ocean_sf)) {
  oce_eol <- rast(g("{OCE_EOL_LOC}/OCE_EOL_2023.tif")) %>%
    terra::crop(st_transform(ra_buffer, st_crs(.)))
  oce_eol_df <- terra::as.data.frame(oce_eol, xy = TRUE)
  oe_sf <- oce_eol_df |>
    filter(OCE_EOL_2023 == 20) |>
    st_as_sf(coords = c("x", "y"), crs = st_crs(oce_eol)) |>
    st_transform(ont.proj)
  write_rds(oe_sf, ocean_sf)
} else {
  oe_sf <- read_rds(ocean_sf)
}

dist_locs <- nngeo::st_nn(locs, oe_sf, returnDist = T)
dist_locs_a2 <- nngeo::st_nn(locs_atlas2, oe_sf, returnDist = T)
if (!run_a2) {
  dist_locs$dist |>
    list_c() |>
    write_rds(g("output/rds/{Sys.Date()}_dist2ocean.rds"))
} else {
  dist_locs_a2$dist |>
    list_c() |>
    write_rds(g("output/rds/{Sys.Date()}_dist2ocean_a2.rds"))
}

