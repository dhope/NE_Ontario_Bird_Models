library(tidyverse)
library(sf)
library(tidymodels)
library(offsetreg)
library(patchwork)
library(tictoc)

library(DALEXtra)
run_a2 <- FALSE
g <- glue::glue

if(str_detect(osVersion, "Windows")){
  source("R/__paths.R")
} else{
  source("R/__paths_linux.R")
}

# # Spatial covariates compiled in )3_Site_Data.R
# spatial_cov <-
#   "output/rds/2025-10-28_spatial_covariates_data.rds" |>
#   read_rds() |>
#   mutate(d2O = read_rds("output/rds/2025-10-28_dist2ocean.rds")) |>
#   bind_cols(individual_locs) |> #Spatial_covariates_data_14March2024.rds") |>
#   distinct() |>
#   st_as_sf() %>%
#   bind_cols(as_tibble(st_coordinates(.))) |>
#   dplyr::select(where(~ sum(is.na(.x)) != length(.x))) |>
#   left_join(st_drop_geometry(aggregated_locs), by = join_by(site_id_agg))
# 
# # Testing and training data compiled in 01_Withold_test_sites.R
# test_training_data <- read_rds(g(
#   "output/{date_compiled}_test_training_data.rds"
# ))

# # Pull out training locations
# train_locs <- test_training_data$train_recordings |>
#   # left_join(aggregated_locs, by = join_by(project, location, collection)) |>
#   st_as_sf() |>
#   dplyr::distinct(site_id_agg, geometry) |>
#   filter(!st_is_empty(geometry))

# Pull out training locations from train locations
# recordings <- test_training_data$train_recordings |>
#   # left_join(aggregated_locs, by = join_by(project, location, collection)) |>
#   filter(str_detect(
#     project,
#     "(Extraction)|(Nocturnal)|(Resample)",
#     negate = T
#   )) |>
#   # Assign recordings to day period
#   mutate(
#     Time_period = dplyr::case_when(
#       t2ss >= -60 & t2ss <= 150 ~ "Dusk",
#       t2sr >= -70 & t2sr <= 220 ~ "Dawn",
#       t2ss > 150 & t2sr < -70 ~ "Night",
#       abs(t2sr) > 220 ~ "Day",
#       is.na(t2sr) ~ "Missing time or location data",
#       TRUE ~ "Unk"
#     ),
#     dur = round(clip_length_min, 2),
#     # t2se is the time to either the sunrise or sunset depending on time period
#     t2se = dplyr::case_when(
#       Time_period == "Dusk" ~ t2ss,
#       Time_period == "Dawn" ~ t2sr,
#       TRUE ~ pmin(abs(t2ss), abs(t2sr))
#     ),
#     doy = yday(date),
#     recording_id = as.numeric(recording_id)
#   ) #|>

# rm(test_training_data, train_locs)

# # Load in prepared count data
# counts_full <- read_rds(g("{rds_data_loc}/counts.rds"))
# spp_list <- distinct(counts_full, species_name_clean, species_code) |>
#   filter(!is.na(species_code))

# Clean count data for the analysis and drop testing data
# counts <- counts_full |>
#   replace_na(list(collection = "WildTrax")) |>
#   filter(event_id %in% recordings$event_id) |>
#   filter(str_detect(
#     project,
#     "(Extraction)|(Nocturnal)|(Resample)",
#     negate = T
#   )) |>
#   dplyr::select(
#     event_id,
#     location,
#     project,
#     species_name_clean,
#     total_count,
#     species_code,
#     total_count_with_tmtt
#   ) |>
#   mutate(
#     y = ifelse(is.na(total_count_with_tmtt), total_count, total_count_with_tmtt)
#   ) |>
#   filter(!is.na(y))

# rm(counts_full)
# 
# qpad_offsets <- read_rds("output/QPAD_global_offsets.rds") |>
#   rename(max_dist = r, time_minutes = t)
# na_pops_offsets <- read_rds("output/na_pops_offsets.rds") |>
#   rename(species = spp)



run_brt <- function(spp) {
  spp_dir <- str_replace_all(spp, ' ', '_')
  out_dir_spatial <- g(
    "{BRT_output_loc_spatial}/{spp_dir}"
  )
  out_dir <- g(
    "{BRT_output_loc}/{spp_dir}"
  )
  gc()
df_std <- read_rds(g("{brt_spp_dat_loc}/{spp_dir}.rds"))

  cat(glue::glue("Starting Model runs for {spp}\n"))
  no_off <- all(is.na(df_std$QPAD_offset))
  xgb_spec <- function(no_off, n_trees) {
    {
      if (!no_off) {
        boost_tree_offset(
          trees = tune(), # tune(),#
          tree_depth = tune(),
          min_n = tune(),
          loss_reduction = tune(), ## first three: model complexity
          sample_size = tune(),
          mtry = NULL, ## randomness
          learn_rate = tune(),
          stop_iter = Inf ## step size
        ) |>
          set_engine(
            "xgboost_offset",
            # params = list(
            # device = 'cuda',
            tree_method = "hist",
            # sampling_method = 'gradient_based'
            # ),
            objective = 'count:poisson',
            offset_col = 'QPAD_offset'
          )
      } else {
        boost_tree(
          trees = tune(), #1500,#tune(),#
          tree_depth = tune(),
          min_n = tune(),
          loss_reduction = tune(), ## first three: model complexity
          sample_size = tune(),
          mtry = NULL, ## randomness
          learn_rate = tune(),
          stop_iter = Inf ## step size
        ) |>
          set_engine(
            "xgboost",
            # tree_method = "hist", # nthread=32,
            objective = 'count:poisson'
          )
      }
    } %>%
      set_mode("regression")
  }

  xgb_rec <-
    df_std %>%
    recipe() %>%
    {
      if (no_off) {
        step_rm(., QPAD_offset)
      } else {
        .
      }
    } |>
    update_role(y, new_role ='outcome') |> 
    update_role(-y, new_role= 'predictor') |> 
    update_role(event_id, new_role = "id") %>%
    # update_role(QPAD_offset, new_role = "offset") %>%

    step_rm(
      contains("location"),
      contains("project"),
      contains("total"),
      contains("d2s"),
      contains("insct"),
      contains("event_id"),
      contains("collection"),
      contains("OHN_swatd_500"),
      contains("dem_asp8")
    ) |> #,X,Y ) |>
    step_mutate(NFIS_is_forest_100 = factor(ifelse(is.na(NFIS_age_100), 1, 0)),
                NFIS_is_forest_500 = #across(contains("NFIS_age_500"),
                                            factor(ifelse(is.na(NFIS_age_100), 1, 0)),
    ) |> 
    step_novel() |>
    step_unknown(all_factor_predictors(), -starts_with('QPAD_offset')) |>
    step_zv(all_predictors(), -starts_with('QPAD_offset')) |>
    step_impute_linear(contains("Climate"),
                    impute_with = imp_vars(c(X,Y) ) ) |> 
    step_impute_median(all_numeric_predictors()) |>
    step_mutate(
      across(
        contains("rec30") & where(is.factor) & contains("fire"),
        ~ case_when(
          as.numeric(as.character(.x)) == 0 ~
            as.numeric(as.character(year)) - 1955,
          TRUE ~ as.numeric(as.character(year)) - as.numeric(as.character(.x))
        )
      ),
      across(
        contains("rec30") & where(is.factor) & contains("harvest"),
        ~ case_when(
          as.numeric(as.character(.x)) == 0 ~
            as.numeric(as.character(year)) - 1980,
          TRUE ~ as.numeric(as.character(year)) - as.numeric(as.character(.x))
        )
      )
    ) |>
    
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
      all_numeric_predictors(),
      -starts_with('QPAD_offset'),
      threshold = 0.95
    )

  # tic("Prep")
  # df_prepped <- prep(xgb_rec) |> bake(new_data = NULL) |> dplyr::select(-y)
  # 
  # toc()
  if (no_off) {
    xgb_wf <- xgb_rec |>
      workflow() %>%
      add_model(xgb_spec(no_off, 300))
  } else {
    xgb_wf <- xgb_rec |>
      workflow() %>%
      add_model(xgb_spec(no_off, 300))
  } #, formula = y ~  offset(QPAD_offset) + .) }

  xgb_grid <- grid_space_filling(
    type = "latin_hypercube", # grid_latin_hypercube(
    tree_depth(),
    trees(c(1L, 8000L)),
    min_n(),
    loss_reduction(),
    sample_size = sample_prop(),
    # finalize(mtry(), df_prepped),
    learn_rate = learn_rate(range = c(-10, -1)),
    # stop_iter(range = c(1, 500)),
    size = 30
  )

  withr::with_seed(123, {
    vb_folds <- vfold_cv(df_std)
  })
  xgb_param <-
    xgb_wf %>%
    extract_parameter_set_dials() %>%
    # update(mtry = finalize(mtry(), df_prepped)) |>
    update(sample_size = sample_prop())
  # library(future)
  # # options(future.globals.maxSize = 8000 * 1024^2)
  # plan(multicore, workers = 32)

  tic("Tune grid")
  xgb_res <- tune_grid(
    xgb_wf,
    resamples = vb_folds,
    # param_info = xgb_param,
    grid = xgb_grid,
    metrics = metric_set(
      poisson_log_loss,
      mae,
      rmse,
      ccc,
      poisson_log_loss,
      rsq
    ), #$,
    control = control_grid(
      verbose = F, #verbose_elim = T,
      parallel_over = 'everything'
    )
  )
  toc()

  # plan(sequential)
  # plan(multisession, workers = 10)
  tic("Tune bayes")
  xgb_bayes <- tune_bayes(
    xgb_wf,
    resamples = vb_folds,
    param_info = xgb_param,
    metrics = metric_set(poisson_log_loss, mae, rmse, ccc, rsq),
    initial = xgb_res,
    iter = 70,
    control = control_bayes(
      parallel_over = 'resamples',
      allow_par = T,
      seed = 3135,
      uncertain = 8,
      no_improve = 10
    )
  )
  toc()
  # plan(sequential)
  show_notes(.Last.tune.result)
  cat(glue::glue("Grid tune finished for {spp}\n"))
  metrics <- xgb_bayes %>%
    collect_metrics()

  write_csv(metrics, g("{out_dir}/{spp}_time_metrics.csv"))

  best_rmse <- select_best(xgb_bayes, metric = "poisson_log_loss")

  # m_d <- metrics[metrics$.metric == "poisson_log_loss", ] |>
  #   pivot_longer(cols = 1:6, names_to = "parameter", values_to = 'value') |>
  #   filter(mean<0.55)
  #   ggplot(m_d,
  #     aes(value, mean)) +
  #   geom_pointrange(aes(ymin = mean - std_err, ymax = mean + std_err)) +
  #   geom_point(
  #     data = best_rmse |>
  #       pivot_longer(cols = 1:6, names_to = "parameter", values_to = 'value') |>
  #       left_join(m_d),
  #     colour = 'red',
  #   ) +
  #   facet_wrap(~parameter, scales = 'free')

  tic("Finalize model")

  final_xgb <- finalize_workflow(
    xgb_wf,
    # update_model(xgb_wf, xgb_spec(no_off, 15000)),
    best_rmse
  )

  fit_workflow <- workflows::fit(final_xgb, df_std)
  toc()
  # extract_fit_engine(fit_workflow) |> collect_metrics()
  saveRDS(
    fit_workflow,
    file = glue::glue("{bundle_locs}/{spp}_time_bundle.rds")
  )
  # plan(sequential)
  cat(glue::glue("Final model fit for {spp}\n"))

 
    t_pred <-
      expand_grid(
        t2se = seq(-60, 240),
        event_gr = c("Dawn", "Dusk"),
        doy = 122:202
      ) |> bind_cols(df_std |> select(-t2se, -doy) |> 
                       summarize(across(where(is.numeric), \(x) median(x, na.rm=T)) )
                      ) |> 
      bind_rows(df_std |> select(-t2se, -doy) ) |>
      filter(!is.na(t2se))

    pred_t <- predict(fit_workflow, t_pred)
    t_pred$.pred <- pred_t$.pred

    time_plot <- ggplot(t_pred, aes(t2se, .pred, colour = event_gr)) +
      stat_summary(fun = 'mean', geom = 'line', aes(group = event_gr)) +
      rcartocolor::scale_colour_carto_d() +
      labs(x = "Time to Sun event", y = "Density", colour = NULL)

    doy_plot <- ggplot(
      t_pred,
      aes(ymd("2025-01-01") + doy, .pred, colour = event_gr)
    ) +
      stat_summary(fun = 'mean', geom = 'line', aes(group = event_gr)) +
      rcartocolor::scale_colour_carto_d() +
      labs(x = NULL, y = "Density", colour = NULL)

    ggsave(
      plot = time_plot / doy_plot,
      file = g("{out_dir}/{spp}_time_date_marginal.jpeg")
    )

    vi_tbl <- vip::vi(fit_workflow)
    write_csv(vi_tbl, g("{out_dir}/{spp}_time_vi.csv"))

    vip_plot <- vip::vip(fit_workflow, geom = "point", n = 15)

    ggsave(
      plot = vip_plot,
      filename = g("{out_dir}/{spp}_time_vi.jpeg"),
      width = 8,
      height = 6
    )
    try({
    explainer <- fit_workflow %>%
      DALEXtra::explain_tidymodels(
        data = df_std %>% select(-y),
        y = df_std$y,
      )
    # write_rds(explainer, glue::glue("{bundle_locs}/{spp}_time_explainer.rds"))
    tic()
    preds_org <- predict(fit_workflow, new_data = df_std) |> bind_cols(df_std)
    toc()
    # collect_metrics(fit_workflow)
    write_rds(preds_org, g("{out_dir}/{spp}_time_model_pred.rds"))

    vars_ <- vi_tbl |> slice_head(n = 15) |> pull(Variable)

    # DALEX::model_diagnostics(explainer, vars_[vars_ %in% names(explainer$data)])
    pd <- DALEX::model_profile(
      explainer = explainer,
      variables = vars_[vars_ %in% names(explainer$data)]
    )
    pdp <- plot(pd)
    ggsave(
      pdp,
      filename = g("{out_dir}/{spp}_time_pdp.jpeg"),
      width = 12,
      height = 12
    )

    resid_brt <- DALEX::model_performance(explainer)

    # create comparison plot of residuals for each model
    p1 <- plot(resid_brt)
    p2 <- plot(resid_brt, geom = "boxplot")

    resid_plt <- p1 + p2
    ggsave(
      plot = resid_plt,
      g("{out_dir}/{spp}_time_resid.jpeg"),
      width = 10,
      height = 6
    )
  })

  # run_brt_preds <- function(spp) {
  # Predict ------------------------------
  # out_dir_spatial <- g(
  #   "{BRT_output_loc_spatial}/{spp}"
  # )
  # out_dir <- g(
  #   "{BRT_output_loc}/{spp}"
  # )
  tic("Run predictions")
  
  period_to_use <- t_pred |> summarize(.pred = median(.pred), .by = c(event_gr)) |> 
    slice_max( order_by = .pred, 
                             n=1,with_ties = F) |> 
   pull(event_gr)
  
  doy_used <- t_pred |> filter( event_gr == period_to_use) |> 
    summarize(.pred = median(.pred), .by = c(doy)) |> 
    slice_max( order_by = .pred, 
               n=1,with_ties = F) |> 
    pull(doy)
  


  preds <- read_rds(g(
    "{prediction_layer_loc}/2025-10-22_prediction_rasters_df.rds"
  ))
  preds$d2O <- read_rds(g(
    "{prediction_layer_loc}/distance_2_ocean_prediction.rds"
  ))
  # preds_ez <- preds |>
  #   st_as_sf(coords = c("X", "Y"), crs = ont.proj) |>
  #   select() |>
  #   st_join(select(ontario_ez, on_er = ZONE_NAME))
  # preds$on_er <- preds_ez$on_er
  # rm(preds_ez)
  # xy <- read_rds("{prediction_layer_loc}/2025-02-28_prediction_rasters_xy.rds")

  # fit_workflow <- #bundle::unbundle(read_rds(
  #   readRDS(
  #     glue::glue("{bundle_locs}/{spp}_time_bundle.rds")
  #   )

  # 
  # library(future)
  # plan(multisession, workers = 32)
  p <- predict(
    fit_workflow,
    new_data = preds |>
      mutate(
        doy = doy_used,
        t2se = 0, #30,
        Time_period = period_to_use,
        year = factor(2026),
        event_gr = "Dawn",
        source = NA,
        collection = NA,
        type = NA,
        Rec_length = (5),
        SiteN = NA,
        site_id =NA,
        longitude = NA,
        latitude = NA,
        event_id = NA,
        QPAD_offset = NA,
        olcb_99_100 = NA,
        total_100 = NA,
        olcb_99_500 = NA,
        total_500 = NA
      )
  )
  toc()

  # plan(sequential)
  write_rds(p, glue::glue("{prediction_layer_loc}/{spp}_time_pred_brt.rds"))
  r_pred <- rast(g("{prediction_layer_loc}/2025-10-21_prediction_rasters.nc")) #2025-02-10_prediction_rasters.nc")
  pobs <- 1 - dpois(0, lambda = p$.pred)
  r_pobs <- predicted_raster <- rast(r_pred[[1]])
  r2 <- terra::cellFromXY(predicted_raster, preds[, c("X", "Y")])
  predicted_raster[r2] <- p$.pred
  r_pobs[r2] <- pobs
  predicted_raster <- mask(predicted_raster, aoi)
  r_pobs <- mask(r_pobs, aoi)
  names(r_pobs) <- names(predicted_raster) <- spp

  terra::writeRaster(
    predicted_raster,
    glue::glue("{out_dir_spatial}/{spp}_time_brt.tif"),
    overwrite = T
  )
  terra::writeRaster(
    r_pobs,
    glue::glue("{out_dir_spatial}/{spp}_pobs_time_brt.tif"),
    overwrite = T
  )
  
  min_q <- min(p, na.rm = T)

  mf_route <- read_sf(
    g("{mf_route_loc}/MFCAR_Route.shp")
  )

  mf_r <- crop(
    predicted_raster,
    st_bbox(st_transform(st_buffer(mf_route, 20000), st_crs(predicted_raster)))
  )
  # quat_loc <- st_crop(quat, st_bbox(st_transform(st_buffer(mf_route, 20000), st_crs(predicted_raster))))

  # cairo_pdf(g("{out_dir}/{spp}_time_map_MFCAR.pdf"), fallback_resolution = 300, height = 17, width = 11)
  # ggplot( ) + tidyterra::geom_spatraster(data = mf_r) +
  #   geom_sf(data = mf_route, alpha = 0.5, linetype =2) +
  #   tidyterra::scale_fill_whitebox_c(
  #     palette = "muted",#breaks = c(0.05, 0.1, 0.2, 0.3, 0.8),
  #     # labels = scales::label_number(suffix = "indiv/ha"),
  #     n.breaks = 8,#limits= c(min_q, quantile(p, 0.999, na.rm=T)),
  #     guide = guide_legend(reverse = TRUE))  +
  #   labs(fill = "individuals/ha", title = spp,
  #        colour = "#\nNon-zero\ncounts") +
  #   geom_sf(data = quat_loc, fill = NA, colour = 'white')
  # dev.off()
  blob <- read_sf(blob_loc) |>
    st_transform(ont.proj)
  # map_plot <- ggplot() +
  #   tidyterra::geom_spatraster(data = predicted_raster, maxcell = 2e6) +
  #   labs(fill = "individuals/ha", title = spp, colour = "#\nNon-zero\ncounts") +
  #   tidyterra::scale_fill_whitebox_c(
  #     palette = "muted", #breaks = c(0.05, 0.1, 0.2, 0.3, 0.8),
  #     # labels = scales::label_number(suffix = "indiv/ha"),
  #     n.breaks = 8,
  #     limits = quantile(p, c(0, 0.999), na.rm = T),
  #     guide = guide_legend(reverse = TRUE)
  #   ) +
  #   # geom_sf(data = napken_lake, shape =2 )+
  #   # geom_sf(data = locs_in) +
  #   geom_sf(data = ra_area, fill = NA, linetype = 2, colour = 'white') +
  #   geom_sf(data = mesh_data, fill = NA, linetype = 3, colour = 'black')

  can <- st_read(canada_shp_loc) %>%
    # filter(NAME %in% c("Manitoba", "Québec", )) |>
    st_transform(ont.proj, partial = F, check = T)

  inset <- ggplot() +
    geom_sf(data = can, fill = NA) +
    geom_sf(
      data = ra_area_official,
      fill = tidyterra::wiki.colors(1),
      alpha = 0.8
    ) +
    ggthemes::theme_map()

  map_plot <-
    ggplot() +
    tidyterra::geom_spatraster(
      data = mask(predicted_raster, ra_area),
      maxcell = 1e6
    ) +
    labs(
      fill = "individuals/ha",
      title = NULL,
      colour = "#\nNon-zero\ncounts"
    ) +
    tidyterra::scale_fill_wiki_c(
      # palette = "muted", #breaks = c(0.05, 0.1, 0.2, 0.3, 0.8),
      # labels = scales::label_number(suffix = "indiv/ha"),
      # n.breaks = 8,
      limits = quantile(p, c(0, 0.999), na.rm = T),
      guide = guide_legend(reverse = TRUE)
    ) +
    geom_sf(data = blob, fill = 'red', alpha = 0.4, colour = NA) +
    ggthemes::theme_map()

  # ggsave(
  #   "output/CONW_with_inset_rof.jpeg",
  #   test_plot +
  #     patchwork::inset_element(inset, 0.7, 0.6, 1.2, 1),
  #   width = 8,
  #   height = 6
  # )
  # ggsave("output/CONW_no_inset_rof.jpeg", test_plot, width = 8, height = 6)

  p_obs_plot <-
    ggplot() +
    tidyterra::geom_spatraster(
      data = mask(r_pobs, ra_area_official),
      maxcell = 2e6
    ) +
    labs(
      fill = "Probability\nof\nobservation",
      title = spp,
      colour = "#\nNon-zero\ncounts"
    ) +
    tidyterra::scale_fill_wiki_c(
      # palette = "muted", #breaks = c(0.05, 0.1, 0.2, 0.3, 0.8),
      # labels = scales::label_number(suffix = "indiv/ha"),
      n.breaks = 8,
      limits = quantile(pobs, c(0, 0.999), na.rm = T),
      guide = guide_legend(reverse = TRUE)
    ) +
    geom_sf(data = BASSr::ontario, fill = NA) +
    # geom_sf(data = napken_lake, shape =2 )+
    # geom_sf(data = locs_in) +
    geom_sf(data = ra_area, fill = NA, linetype = 2, colour = 'white') +
    # geom_sf(data = mesh_data, fill = NA, linetype = 3, colour = 'black') +
    coord_sf(
      xlim = st_bbox(ra_buffer)[c(1, 3)],
      ylim = st_bbox(ra_buffer)[c(2, 4)]
    )

  ggsave(
    plot = map_plot,
    g("{out_dir}/{spp}_time_map_{period_to_use}_{doy_used}.jpeg"),
    width = 6.73,
    height = 8.5
  )
  ggsave(
    plot = p_obs_plot,
    g("{out_dir}/{spp}_time_map_pobs_{period_to_use}_{doy_used}.jpeg"),
    width = 6.73,
    height = 8.5
  )
  # hist(p[[1]][p[[1]]<0.004], breaks = 100)
  # hasValues(predicted_raster)

  cat(glue::glue("Predictions complete for {spp}\n"))
}
