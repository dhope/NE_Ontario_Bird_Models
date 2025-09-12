library(tidyverse)
library(sf)
library(fmesher)
library(patchwork)
source("R/06_BRT_time.R")

# QPAD::load_BAM_QPAD(4)
library(tidymodels)

wt_locations <- read_rds(g("output/rds/{date_compiled}_wt_locations.rds")) |> 
  distinct()
spatial_cov <- read_rds("output/rds/2024-08_30_Baldwin.rds") |>bind_cols(wt_locations) |>  #Spatial_covariates_data_14March2024.rds") |> 
  distinct() |> st_as_sf() %>%
  bind_cols(as_tibble(st_coordinates(.))) |>
  dplyr::select(where(~sum(is.na(.x))!=length(.x)))


train_locs <-  read_rds(g("output/{date_compiled}_test_training_data.rds"))$train_recordings |> 
  dplyr::distinct(project, location, recording_id)

pca_cov <- distinct(spatial_cov) |> 
  filter(location %in% train_locs$location) |> 
  dplyr::select(location,where(~{!is.factor(.x) & !is.character(.x)})) |> 
  st_drop_geometry() |> dplyr::select(-X, -Y, -site_id_agg)
  

preds_r <- read_rds(g("{prediction_layer_loc}/2025-02-28_prediction_rasters_df.rds") ) |> 
  mutate(rn = row_number()) 
xy <- read_rds(g("{prediction_layer_loc}/2025-02-28_prediction_rasters_xy.rds") )

pca_pred_data <- preds_r |> 
  mutate(location = glue::glue("loc-{1:nrow(preds_r)}")) |>
  dplyr::select( any_of(names(pca_cov)) )



cov_scaled <-  
  recipe(data = pca_cov, formula = ~ .) %>%
  update_role(location, new_role = "id") %>%
  # step_zv(all_predictors()) |> 
  step_nzv(all_predictors(), freq_cut = 95/5, unique_cut = 10) |>
  step_impute_median(all_numeric_predictors()) |> 
  step_corr(all_numeric_predictors(), threshold = 0.85) 


cov_scaled_prep <- prep(cov_scaled)
cov_scaled_out <- bake(cov_scaled_prep, new_data = NULL) |>   
  dplyr::select(-location)

cov_preds <- bake(cov_scaled_prep, new_data = pca_pred_data) |> 
  dplyr::select(-location)


p_zeros_preds <- 
  cov_preds |> 
  summarize(across(where(is.numeric), 
                   ~sum(.x==0)/nrow(cov_preds))) |> 
  pivot_longer(cols = everything()) |> 
  arrange(desc(value))

sum(!names(cov_scaled_out) %in% names(cov_preds))


calc_nt2 <- function(train, pred, removed_vars, pthresh = 0.95){
  train_df <- dplyr::select(train, -any_of(removed_vars))
  pred_df <- dplyr::select(pred, -any_of(removed_vars))
  
  trainMeans <- colMeans(train_df, na.rm = TRUE)
  trainVar <- var(train_df, na.rm = TRUE)
  trainMah <- mahalanobis(train_df, trainMeans,
                        trainVar, na.rm = TRUE)
  thresh <- quantile(trainMah, probs =pthresh, na.rm = TRUE)
  trainMah[which(trainMah > thresh)] <- NA
  # browser()
  maxMah <- max(trainMah, na.rm = TRUE)
  # trainMah <- trainMah/maxMah
  naMah <- mahalanobis(pred_df, 
                       trainMeans,
                       trainVar,
                       na.rm = TRUE)
  naMah/maxMah
}

reframe(pca_out, across(where(is.numeric), range)) |> 
  glimpse()


nt2_full <- calc_nt2(cov_scaled_out, cov_preds,NULL, pthresh = 0.95)



# This needs fixing still ----------------
drop_var <- map(names(cov_scaled_out),
                ~calc_nt2(cov_scaled_out, cov_preds, .x, pthresh=0.95), .progress = T)
                                                   # calc_nt2 = calc_nt2, cov_scaled_out =  cov_scaled_out,
                                                   # cov_preds = cov_preds))

matrixStats::row()
dv_df <- map(1:length(names(cov_scaled_out)),
             ~{you_tibble <- names(cov_scaled_out)[.x]
               tibble({{you_tibble}}:=drop_var[[.x]])
               }) |> list_cbind() 
  
# maxMah <- 
# dv_df |> #slice_head(n=100)|>
#   mutate(rn = row_number()) |>
#   rowwise() |> 
#   mutate(across(-rn, ~{.x-nt2_full[rn]})) |> 
#   dplyr::select(-rn) |> 
#   transmute(maxmah = max(c_across(everything())),
#             Maxcol = colnames(dv_df)[which.max(c_across(everything()))]) |> 
#   ungroup()


diff_mat <- map(1:ncol(dv_df), ~{dv_df[,.x] - nt2_full}, .progress = T)
diffs <- list_cbind(diff_mat)

str(diffs)
max_diffs <- as.matrix(diffs) |> matrixStats::rowMaxs()

      
which_maxes <- max.col(diffs, ties.method = c("random"))
max_diff_cols <- colnames(dv_df)[which_maxes] |> factor()    


r_pred <- rast(glue::glue("{prediction_layer_loc}/2025-02-28_prediction_rasters.nc"))
dim_r <- dim(r_pred[[1]])
res_r <- res(r_pred[[1]])
e_r <- ext(r_pred[[1]])
nt2 <- rast(r_pred[[1]]) 
miv <- rast(e_r,nrows= dim_r[1], ncol = dim_r[2],crs= crs(r_pred[[1]]) )

in_mesh <- fmesher::fm_is_within(st_as_sf(xy, coords = c("x", "y"), crs= 3161), mesh_inla)
r2 <- terra::cellFromXY(r_pred[[1]],xy[in_mesh,])



table_values <- tibble(id = 1:length(levels(max_diff_cols)), 
                       miv = levels(max_diff_cols))


nt2[r2] <- nt2_full[in_mesh] #preds_sc$nt2[in_mesh]
miv[r2] <- as.numeric(max_diff_cols)[in_mesh]
names(miv) <- "miv"
levels(miv) <- table_values
is.factor(miv)

janitor::tabyl(max_diff_cols)

COLS <- c("green", "lightgreen", "aquamarine", "blue", 
          "orange", "brown", "red")

BREAKS = c(0, 0.25, 0.5, 1, 2, 4, 8, 16)

ggplot() +
  tidyterra::geom_spatraster(data =nt2) +
  labs(fill = "Nt2") +
  tidyterra::scale_fill_whitebox_c(
    palette = "muted",breaks = BREAKS,
    # labels = scales::label_number(suffix = "indiv/ha"),
    n.breaks = 8, limits = c(1, max(nt2_full)),
    guide = guide_legend(reverse = TRUE)
  )


miv_ex <- miv
miv_ex[which.lyr(nt2<=1)] <- NA
miv_ex <- terra::droplevels(miv_ex)

janitor::tabyl(values(miv_ex)[,1]) |> 
  as_tibble() |> 
  rename(id = `values(miv_ex)[, 1]`) |> 
  left_join(table_values) |> 
  filter(!is.na(miv)) |> 
  arrange(desc(n))

ggplot() +
  tidyterra::geom_spatraster(data =miv_ex) +
  labs(fill = "Nt2") +
  tidyterra::scale_fill_whitebox_d(
    palette = "muted",#breaks = c(0.05, 0.1, 0.2, 0.3, 0.8),
    # labels = scales::label_number(suffix = "indiv/ha"),
    # n.breaks = 8, 
    # limits = c(0, 1),
    guide = guide_legend(reverse = TRUE)
  ) 



