library(tidyverse)
library(sf)
library(fmesher)
library(patchwork)
source("R/06_BRT_time.R")

# QPAD::load_BAM_QPAD(4)
library(tidymodels)

spatial_cov <-
  "output/rds/2025-10-23_spatial_covariates_data.rds" |>
  read_rds() |>
  mutate(d2O = read_rds("output/rds/2025-10-23_dist2ocean.rds")) |>
  bind_cols(individual_locs) |> #Spatial_covariates_data_14March2024.rds") |>
  distinct() |>
  st_as_sf() %>%
  bind_cols(as_tibble(st_coordinates(.))) |>
  dplyr::select(where(~ sum(is.na(.x)) != length(.x))) |>
  left_join(st_drop_geometry(aggregated_locs), by = join_by(site_id_agg))


train_locs <- read_rds(g(
  "output/{date_compiled}_test_training_data.rds"
))$train_recordings |>
  dplyr::distinct(project, location, recording_id)

pca_cov <- distinct(spatial_cov) |>
  filter(location %in% train_locs$location) |>
  dplyr::select(
    location,
    where(
      ~ {
        !is.factor(.x) & !is.character(.x)
      }
    )
  ) |>
  st_drop_geometry() |>
  dplyr::select(-X, -Y, -site_id_agg)


preds_r <- read_rds(g(
  "{prediction_layer_loc}/2025-10-22_prediction_rasters_df.rds"
)) |> #2025-02-11_prediction_rasters_df.rds") |>
  mutate(
    rn = row_number(),
    d2O = read_rds(g(
      "{prediction_layer_loc}/distance_2_ocean_prediction.rds"
    )),
    event_gr = "Dawn",
    source = NA,
    collection = NA,
    type = NA,
    Rec_length = factor(5),
    SiteN = NA,
    longitude = NA,
    latitude = NA,
    event_id = NA,
    QPAD_offset = NA,
    olcb_99_100 = NA,
    total_100 = NA,
    olcb_99_500 = NA,
    total_500 = NA
  )
# xy <- read_rds(g("{prediction_layer_loc}/2025-02-28_prediction_rasters_xy.rds") )

pca_pred_data <- preds_r |>
  mutate(location = glue::glue("loc-{1:nrow(preds_r)}")) |>
  dplyr::select(any_of(names(pca_cov)))


cov_scaled <-
  recipe(data = pca_cov, formula = ~.) %>%
  update_role(location, new_role = "id") %>%
  step_rm(
    contains("project"),
    contains("total"),
    contains("event_id"),
    contains("collection"),
    where(is.factor)
  ) |> #,X,Y ) |>
  step_novel() |>
  step_unknown(all_factor_predictors(), -starts_with('QPAD_offset')) |>
  step_impute_median(all_numeric_predictors()) |>
  step_mutate(across(
    contains("rec30") & where(is.factor),
    ~ case_when(
      as.numeric(as.character(.x)) == 0 ~ NA_integer_,
      as.numeric(as.character(year)) < as.numeric(as.character(.x)) ~
        NA_integer_,
      TRUE ~ as.numeric(as.character(year)) - as.numeric(as.character(.x))
    )
  )) |>
  step_zv(all_predictors(), -starts_with('QPAD_offset')) |>
  step_center(all_numeric_predictors(), -starts_with('QPAD_offset')) |>
  step_scale(
    all_numeric_predictors(),
    sds = 2,
    -starts_with('QPAD_offset')
  ) |>
  step_dummy(
    all_factor_predictors(),
    -starts_with('QPAD_offset'),
    one_hot = TRUE
  ) |>
  step_nzv(
    all_predictors(),
    ,
    -starts_with('QPAD_offset'),
    freq_cut = 95 / 5,
    unique_cut = 2
  ) |>
  step_corr(
    all_predictors(),
    -starts_with('QPAD_offset'),
    threshold = 0.95
  )


cov_scaled_prep <- prep(cov_scaled)
cov_scaled_out_wl <- bake(cov_scaled_prep, new_data = NULL) #|>
cov_scaled_out <- cov_scaled_out_wl |> dplyr::select(-location)

cov_preds <- bake(cov_scaled_prep, new_data = pca_pred_data) |>
  dplyr::select(-location)


p_zeros_preds <-
  cov_preds |>
  summarize(across(where(is.numeric), ~ sum(.x == 0) / nrow(cov_preds))) |>
  pivot_longer(cols = everything()) |>
  arrange(desc(value))

# sum(!names(cov_scaled_out) %in% names(cov_preds))

# calc_nt2 <- function(train, pred, removed_vars, pthresh = 0.95) {
#   train_df <- dplyr::select(train, -any_of(removed_vars))
#   pred_df <- dplyr::select(pred, -any_of(removed_vars))

#   trainMeans <- colMeans(train_df, na.rm = TRUE)
#   trainVar <- var(train_df, na.rm = TRUE)
#   trainMah <- mahalanobis(train_df, trainMeans, trainVar, na.rm = TRUE)
#   thresh <- quantile(trainMah, probs = pthresh, na.rm = TRUE)
#   trainMah[which(trainMah > thresh)] <- NA
#   # browser()
#   maxMah <- max(trainMah, na.rm = TRUE)
#   # trainMah <- trainMah/maxMah
#   naMah <- mahalanobis(pred_df, trainMeans, trainVar, na.rm = TRUE)
#   naMah / maxMah
# }

# reframe(pca_out, across(where(is.numeric), range)) |>
#   glimpse()

# nt2_full <- calc_nt2(cov_scaled_out, cov_preds, NULL, pthresh = 0.95)

# # This needs fixing still ----------------
# drop_var <- map(
#   names(cov_scaled_out),
#   ~ calc_nt2(cov_scaled_out, cov_preds, .x, pthresh = 0.95),
#   .progress = T
# )
# # calc_nt2 = calc_nt2, cov_scaled_out =  cov_scaled_out,
# # cov_preds = cov_preds))

# matrixStats::row()
# dv_df <- map(
#   1:length(names(cov_scaled_out)),
#   ~ {
#     you_tibble <- names(cov_scaled_out)[.x]
#     tibble({{ you_tibble }} := drop_var[[.x]])
#   }
# ) |>
#   list_cbind()

# # maxMah <-
# # dv_df |> #slice_head(n=100)|>
# #   mutate(rn = row_number()) |>
# #   rowwise() |>
# #   mutate(across(-rn, ~{.x-nt2_full[rn]})) |>
# #   dplyr::select(-rn) |>
# #   transmute(maxmah = max(c_across(everything())),
# #             Maxcol = colnames(dv_df)[which.max(c_across(everything()))]) |>
# #   ungroup()

# diff_mat <- map(
#   1:ncol(dv_df),
#   ~ {
#     dv_df[, .x] - nt2_full
#   },
#   .progress = T
# )
# diffs <- list_cbind(diff_mat)

# str(diffs)
# max_diffs <- as.matrix(diffs) |> matrixStats::rowMaxs()

# which_maxes <- max.col(diffs, ties.method = c("random"))
# max_diff_cols <- colnames(dv_df)[which_maxes] |> factor()

# r_pred <- rast(glue::glue(
#   "{prediction_layer_loc}/2025-10-21_prediction_rasters.nc"
# ))
# dim_r <- dim(r_pred$olcc_1_1km)
# res_r <- res(r_pred$olcc_1_1km)
# e_r <- ext(r_pred$olcc_1_1km)
# nt2 <- rast(r_pred$olcc_1_1km)
# miv <- rast(
#   e_r,
#   nrows = dim_r[1],
#   ncol = dim_r[2],
#   crs = crs(r_pred$olcc_1_1km)
# )

dsm.extrapolation <- dsmextra::compute_extrapolation(
  samples = cov_scaled_out_wl |>
    left_join(
      select(spatial_cov, location, x = X, y = Y),
      by = join_by(location)
    ) |>
    select(-location),
  covariate.names = names(cov_scaled_out),
  prediction.grid = bind_cols(cov_preds, preds_r[, c("X", "Y")]) |>
    rename(x = X, y = Y),
  # st_as_sf(coords = c("X", "Y"), crs = ont.proj),
  coordinate.system = st_crs(spatial_cov) #as(spatial_cov, "Spatial")@proj4string
)

write_rds(
  dsm.extrapolation,
  glue::glue("{BRT_output_loc}/DSM_Extrapolation_results.rds")
)

# dsm.extrapolation <- read_rds(
#   glue::glue("{BRT_output_loc}/DSM_Extrapolation_results.rds")
# )

cls <- tibble(
  id = 1:length(dsm.extrapolation$covariate.names),
  mic = dsm.extrapolation$covariate.names
)

dsm.extrapolation$data$univariate |>
  janitor::tabyl(mic_univariate) |>
  left_join(cls, by = join_by(mic_univariate == id)) |>
  arrange(desc(n)) |>
  relocate(mic, .before = 1) |>
  write_csv(g("{prediction_layer_loc}/{Sys.Date()}_MIC_univariate.csv"))


dsm.extrapolation$data$combinatorial |>
  janitor::tabyl(mic_combinatorial) |>
  left_join(cls, by = join_by(mic_combinatorial == id)) |>
  arrange(desc(n)) |>
  relocate(mic, .before = 1) |>
  write_csv(g("{prediction_layer_loc}/{Sys.Date()}_MIC_combinatorial.csv"))


# summary(dsm.extrapolation)
# plot(r_pred$sfi_cls_1km == 0)
# swatd_500 <- rast(g("{spatial_raster_location}/8_OHN/ohn_swatd_1km"))
# cls <- rast(g("{spatial_raster_location}/4_SCANFI/cls_orig30m"))
# NFIS_age_100 <- rast(g("{spatial_raster_location}/5_NFIS/nfi_age_200"))
# NFIS_d2s_500 <- rast(g("{spatial_raster_location}/13_NFIS2/nfi_d2s_1km"))
# nicedays <- rast(g("{spatial_raster_location}/7_Climate/clm_niced_1km"))
# plot(swatd_500 == 0)
# plot(NFIS_age_100 == 0)
# plot(nicedays == 0)

# ggplot(spatial_cov, aes(fnlc_1_500)) +
#   geom_density() +
#   geom_density(fill = 'red', alpha = 0.4, data = pca_pred_data) # +
# scale_x_continuous(transform = scales::log1p_trans())

# ggplot(spatial_cov, aes(NFIS_age_100)) +
#   geom_density() +
#   geom_density(fill = 'red', alpha = 0.4, data = pca_pred_data) +
#   scale_x_continuous(transform = scales::log1p_trans())

# dsm.extrapolation$rasters$ExDet$combinatorial |> plot()
mic_r <- dsm.extrapolation$rasters$mic$univariate
mic_r_c <- dsm.extrapolation$rasters$mic$combinatorial
levels(mic_r_c) <- levels(mic_r) <- cls
# plot(mic_r)
# dsm.extrapolation$covariate.names |> str()

terra::writeRaster(
  mic_r,
  g("{prediction_layer_loc}/{Sys.Date()}_MIC_univariate.tif")
)
foreign::write.dbf(
  levels(mic_r)[[1]],
  file = g("{prediction_layer_loc}/{Sys.Date()}_MIC_univariate.tif")
)

terra::writeRaster(
  mic_r_c,
  g("{prediction_layer_loc}/{Sys.Date()}_MIC_combinatorial.tif")
)
foreign::write.dbf(
  levels(mic_r_c)[[1]],
  file = g("{prediction_layer_loc}/{Sys.Date()}_MIC_combinatorial.tif")
)

walk(
  names(dsm.extrapolation$rasters$ExDet),
  ~ {
    terra::writeRaster(
      dsm.extrapolation$rasters$ExDet[[.x]],
      g("{prediction_layer_loc}/{Sys.Date()}_ExDet_{.x}.tif")
    )
  }
)

# dsmcc <-
#   dsmextra::compare_covariates(
#     extrapolation.type = "both",
#     extrapolation.object = dsm.extrapolation,
#     n.covariates = NULL,
#     create.plots = TRUE,
#     display.percent = TRUE
#   )

# obs.nearby <- dsmextra::compute_nearby(
#   samples = cov_scaled_out_wl |>
#     left_join(
#       select(spatial_cov, location, x = X, y = Y),
#       by = join_by(location)
#     ) |>
#     select(-location),
#   prediction.grid = bind_cols(cov_preds, preds_r[, c("X", "Y")]) |>
#     rename(x = X, y = Y),
#   coordinate.system = as(spatial_cov, "Spatial")@proj4string,
#   covariate.names = names(cov_scaled_out),
#   nearby = 1
# )

# in_mesh <- fmesher::fm_is_within(
#   st_as_sf(preds_r[, c("X", "Y")], coords = c("X", "Y"), crs = 3161),
#   mesh_inla
# )
# r2 <- terra::cellFromXY(r_pred[[1]], preds_r[in_mesh, c("X", "Y")])

# table_values <- tibble(
#   id = 1:length(levels(max_diff_cols)),
#   miv = levels(max_diff_cols)
# )

# nt2[r2] <- nt2_full[in_mesh] #preds_sc$nt2[in_mesh]
# miv[r2] <- as.numeric(max_diff_cols)[in_mesh]
# names(miv) <- "miv"
# levels(miv) <- table_values
# is.factor(miv)

# janitor::tabyl(max_diff_cols)

# COLS <- c("green", "lightgreen", "aquamarine", "blue", "orange", "brown", "red")

# BREAKS = c(0, 0.25, 0.5, 1, 2, 4, 8, 16)

# ggplot() +
#   tidyterra::geom_spatraster(data = nt2) +
#   labs(fill = "Nt2") +
#   tidyterra::scale_fill_whitebox_c(
#     palette = "muted",
#     breaks = BREAKS,
#     # labels = scales::label_number(suffix = "indiv/ha"),
#     n.breaks = 8,
#     limits = c(1, max(nt2_full)),
#     guide = guide_legend(reverse = TRUE)
#   )

# miv_ex <- miv
# miv_ex[which.lyr(nt2 <= 1)] <- NA
# miv_ex <- terra::droplevels(miv_ex)

# janitor::tabyl(values(miv_ex)[, 1]) |>
#   as_tibble() |>
#   rename(id = `values(miv_ex)[, 1]`) |>
#   left_join(table_values) |>
#   filter(!is.na(miv)) |>
#   arrange(desc(n))

# ggplot() +
#   tidyterra::geom_spatraster(data = miv_ex) +
#   labs(fill = "Nt2") +
#   tidyterra::scale_fill_whitebox_d(
#     palette = "muted", #breaks = c(0.05, 0.1, 0.2, 0.3, 0.8),
#     # labels = scales::label_number(suffix = "indiv/ha"),
#     # n.breaks = 8,
#     # limits = c(0, 1),
#     guide = guide_legend(reverse = TRUE)
#   )
