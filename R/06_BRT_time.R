library(tidyverse)
library(sf)
library(tidymodels)
library(offsetreg)
library(patchwork)
library(tictoc)

library(DALEXtra)
source("R/__globals.R")
ggplot2::theme_set(ggplot2::theme_minimal(base_size = 14, base_family = "Roboto Condensed") )
aoi <- read_rds( g("output/rds/{date_compiled}_Area_of_focus_mesh.rds"))
# Load data ---------

spatial_cov <- 
  "output/rds/2025-02-28_spatial_covariates_data.rds" |> 
  read_rds() |>bind_cols(individual_locs) |>  #Spatial_covariates_data_14March2024.rds") |> 
  distinct() |> st_as_sf() %>%
  bind_cols(as_tibble(st_coordinates(.))) |>
  dplyr::select(where(~sum(is.na(.x))!=length(.x))) |> 
  left_join(st_drop_geometry(aggregated_locs ),
            by = join_by(site_id_agg))

test_training_data <- read_rds(g("output/{date_compiled}_test_training_data.rds"))

train_locs <-  test_training_data$train_recordings |> 
  left_join(aggregated_locs,by = join_by(project, location, collection)) |> 
  st_as_sf() |> 
  dplyr::distinct(site_id_agg, geometry) |> 
  filter(!st_is_empty(geometry))


recordings <- test_training_data$train_recordings |> 
  left_join(aggregated_locs,by = join_by(project, location, collection)) |> 
  filter(str_detect(project, "(Extraction)|(Nocturnal)|(Resample)",negate=T)) |>
  mutate(
    Time_period =dplyr::case_when(
      t2ss >= -60 & t2ss <= 150~"Dusk",
      t2sr >= -70 & t2sr <= 220 ~"Dawn",
      t2ss > 150 & t2sr < -70 ~ "Night",
      abs(t2sr)>220 ~ "Day" ,
      is.na(t2sr) ~ "Missing time or location data",
      TRUE~"Unk" ),
    dur = round(clip_length_min, 2),
    t2se = dplyr::case_when(Time_period=="Dusk"~t2ss,
                            Time_period == "Dawn"~t2sr,
                            TRUE ~ pmin(abs(t2ss),
                                        abs(t2sr))),
    doy = yday(date),
    recording_id = as.numeric(recording_id)) #|> 

rm(test_training_data, train_locs)

counts <- read_rds(g("{rds_data_loc}/counts.rds")) |>
  replace_na(list(collection="WildTrax")) |> 
  filter(event_id %in% recordings$event_id ) |> 
  filter(str_detect(project, "(Extraction)|(Nocturnal)|(Resample)",negate=T)) |> 
  dplyr::select(event_id, location,  project, species_name_clean,total_count, 
                species_code,
                total_count_with_tmtt) |> 
  mutate(y = ifelse(is.na(total_count_with_tmtt), total_count,
                    total_count_with_tmtt)) |> 
  filter(!is.na(y))

qpad_offsets <- read_rds("output/QPAD_global_offsets.rds") |> 
  rename(max_dist =r, time_minutes = t)
na_pops_offsets <- read_rds("output/na_pops_offsets.rds") |> 
  rename(species =spp)  

prep_brt_data <- function(spp){
  flush.console()
  
  spp_name <- counts$species_name_clean[counts$species_code==spp & !is.na(counts$species_code)] |> 
    unique() 
  if(length(spp_name)!=1)rlang::abort("unable to identify species name")
  counts_spp <- filter(counts, species_name_clean == spp_name) |> 
    full_join(recordings |> dplyr::select(-geometry),
              by = join_by(event_id, location, project)) |> 
    replace_na(list(total_count_with_tmtt = 0, total_count = 0,
                    y=0))
  # offfiles <-g("output/species_rds/offsets_{spp}.rds")
  offsets_spp <- distinct(
    recordings, location, event_id,dur,species_code = spp) |> 
    mutate(o = 0, dur = round(dur)) 
  
  
  
  if(spp %in% na_pops_offsets$species[!is.na(na_pops_offsets$o)]) {
    offsets_spp <- left_join(
      offsets_spp |> dplyr::select(-o)
      ,
      filter(na_pops_offsets,species == spp&max_dist == Inf),
      by = join_by(species_code == species,
                   dur == time_minutes
      ) )
    
    
  }else{
    if(spp %in% qpad_offsets$species) {
      # offsets_spp <- filter(qpad_offsets, species == spp)
      offsets_spp <- left_join(
        offsets_spp|> dplyr::select(-o),
        filter(qpad_offsets,species == spp&max_dist == Inf),
        by = join_by(species_code == species,
                     dur == time_minutes
        ) )
    }
  }
  
  if(!"max_dist" %in% names(offsets_spp)){
    offsets_spp$max_dist <- Inf
  }
  
  out_dir_spatial <- g(
    "{BRT_output_loc_spatial}/{spp}"
  )
  out_dir <- g(
    "{BRT_output_loc}/{spp}"
  )
  dir.create(
    out_dir_spatial, recursive = T
  )
  dir.create(
    out_dir, recursive = T
  )
  
  setup_dat_0 <- counts_spp |> 
    left_join(offsets_spp |> 
                dplyr::select(event_id, offset =o),
              by = join_by( event_id)) |> 
    # filter(species_code==spp & str_detect(project, "Resample", negate=T) &
    #          recording_id %in% train_locs$recording_id) |>
    # 
    # left_join(recordings |> 
    #             dplyr::select(project:recording_id, date_time:t2se, dur),
    #           by = join_by(project, project_id, location, recording_id)
    # ) |>
    mutate(Rec_length = factor(as.numeric(round(dur,1)))) |> 
    mutate(
      SiteN = as.numeric(factor(location)),
      rec_id = as.numeric(as.factor(event_id)),
      event_gr = factor(case_when(
        is.na(date_time)~NA_character_,
        Time_period %in% c("Dawn", "Dusk", "Night")~Time_period,
        Time_period == "Day" & hour(date_time) <=10 ~"Dawn",
        Time_period == "Day" & hour(date_time) > 10 ~"Day",
        TRUE~"Other"
      ))
    ) |> 
    arrange(location,doy, Time_period,t2se) |> 
    mutate(    Row = row_number())
 
  
  included_times <- setup_dat_0 |> filter(y>0) |> 
    janitor::tabyl(event_gr) |> 
    filter(n>0)
  
  
  setup_dat_nested <- 
    setup_dat_0 |> 
    filter(event_gr %in% included_times$event_gr) |> 
    nest_by(event_gr) 
  
  
  rm(setup_dat_0)
  
  setup_dat <- 
    setup_dat_nested |> 
    filter(event_gr!="Day") |> unnest(data)
  
  df_std <- setup_dat |>st_drop_geometry() |> 
    ungroup() |> 
    dplyr::select(-c(project_id,  longitude, latitude,SiteN,
                     species_name_clean,total_count, total_count_with_tmtt,
                     species_code,  recording_id,#event_id,
                     tz,t2sr, t2ss,SamplingEventIdentifier,
                     Time_period, month, day, path, site_id, aru_id, date_time,
                     date,clip_length_min,dur,Row, rec_id
    )) %>% 
    mutate(random_val = rnorm(n=nrow(.))) |> 
    left_join(spatial_cov |> st_drop_geometry() 
              ,by = join_by(location,project, collection, site_id_agg) ) |> 
    dplyr::select(-site_id_agg)
  
  as.list.environment(environment())
  
}


run_brt <- function(spp){
  gc()
  
  list2env(prep_brt_data(spp), environment())
  
  
  cat(glue::glue("Starting Model runs for {spp}\n"))
  no_off <- all(is.na(df_std$offset))
  xgb_spec <- 
    {if(!no_off){
      boost_tree_offset(
        trees =1500,# tune(),#
        tree_depth = tune(),
        min_n = tune(),
        loss_reduction = tune(),                     ## first three: model complexity
        sample_size = tune(), mtry = tune(),         ## randomness
        learn_rate = tune() , stop_iter = tune()                       ## step size
      ) |> 
        set_engine( "xgboost_offset", params=list( nthread = 1,tree_method = "hist"),
                    objective = 'count:poisson', 
                    offset_col = 'offset')
    } else{
      boost_tree(
        trees = 1500,#tune(),#
        tree_depth = tune(),
        min_n = tune(),
        loss_reduction = tune(),                     ## first three: model complexity
        sample_size = tune(), 
        mtry = tune(),         ## randomness
        learn_rate = tune()   , stop_iter = tune()                       ## step size
      ) |> 
        set_engine("xgboost",tree_method = "hist",# nthread=32,
                   objective = 'count:poisson')
    }
    }%>%
    set_mode("regression")

  xgb_rec <- 
    df_std %>%
    recipe(y ~ .) %>%
    {
      if(no_off){
        step_rm(.,offset)
      } else{.}
    } |> 
    update_role(event_id, new_role = "id") %>%
    update_role(offset, new_role = "offset") %>%
    step_rm(contains("location"), contains("project")) |> 
    step_novel() |> 
    step_unknown(all_factor_predictors(), -starts_with('offset')) |> 
    step_dummy(all_factor_predictors(), -starts_with('offset'), one_hot=TRUE) |> 
    step_impute_median(all_numeric_predictors()) |> 
    step_zv(all_predictors(),-starts_with('offset')) |> 
    step_center(all_numeric_predictors(),-starts_with('offset')) |> 
    step_scale(all_numeric_predictors(), sds = 2,-starts_with('offset')) |> 
    step_corr(all_numeric_predictors(),,-starts_with('offset'), threshold = 0.85) 
  
  tic("Prep")
  df_prepped <- prep(xgb_rec) |> bake(new_data = NULL) |> dplyr::select(-y)
  toc()
  xgb_wf <-xgb_rec |>
    workflow() %>% 
    add_model(xgb_spec, formula = y ~  + .) 

    xgb_grid <- grid_space_filling(type = "latin_hypercube",# grid_latin_hypercube(
                                 tree_depth(),#trees(c(200L, 10000L)),
                                 min_n(),
                                 loss_reduction(),
                                 sample_size = sample_prop(),
                                  finalize(mtry(), df_prepped),
                                 learn_rate =learn_rate(range = c(-5, -0.5)),
                                           stop_iter(range = c(1, 500))
                                            ,
                                 size = 10 )
  
  
  
  withr::with_seed(123,{
    vb_folds <- vfold_cv(df_std)
  })
  xgb_param <- 
    xgb_wf %>% 
    extract_parameter_set_dials() %>% 
    update(mtry = finalize(mtry(), df_prepped)) |>
    update(sample_size = sample_prop())
  library(future)
  options(future.globals.maxSize = 8000 * 1024^2)
  plan(multisession, workers = 32)
  
  tic("Tune grid")
  xgb_res <- tune_grid(
    xgb_wf,
    resamples = vb_folds,
    # param_info = xgb_param,
    grid = xgb_grid,
    metrics = metric_set(poisson_log_loss,
                         mae, rmse, ccc, poisson_log_loss,rsq),#$,
    control = control_grid(verbose = F,#verbose_elim = T,
                           parallel_over = 'resamples')
  )
  toc()
  
    
  tic("Tune bayes")
  xgb_bayes <- tune_bayes(xgb_wf, 
                          resamples = vb_folds, 
                          param_info = xgb_param,
                          metrics =  metric_set(poisson_log_loss,mae, rmse, ccc, rsq),
                          initial = xgb_res,
                          iter = 70,
                          control = control_bayes(parallel_over = 'resamples', 
                                                  allow_par = T,seed = 3135,
                                                  uncertain = 8,
                                                  no_improve = 10))
  toc()
  plan(sequential)
  show_notes(.Last.tune.result)
  cat(glue::glue("Grid tune finished for {spp}\n"))
  metrics <- xgb_bayes %>%
    collect_metrics() 
  
  write_csv(metrics, g("{out_dir}/{spp}_time_metrics.csv"))
  

  
  best_rmse <- select_best(xgb_bayes,metric =  "poisson_log_loss") 
  
  tic("Finalize model")
  
  final_xgb <- finalize_workflow(
    xgb_wf,
    best_rmse
  )
  
  
  
  
  
  
  fit_workflow <- fit(final_xgb, df_std)
  toc()
  
  
  saveRDS(fit_workflow, file = glue::glue("{bundle_locs}/{spp}_time_bundle.rds"))

  
  cat(glue::glue("Final model fit for {spp}\n"))
  
  
  try({
    
    
    t_pred <- 
      expand_grid(t2se = seq(-60, 240), event_gr = c("Dawn", "Dusk"), doy = 122:202) |> 
      bind_rows(df_std) |> filter(is.na(X))
    
    pred_t <- predict(fit_workflow, t_pred)
    t_pred$.pred <- pred_t$.pred
    
    
    time_plot <- ggplot(t_pred, aes(t2se, .pred, colour = event_gr)) + 
       stat_summary(fun = 'mean', geom='line', aes(group = event_gr)) +
      rcartocolor::scale_colour_carto_d() +
      labs(x = "Time to Sun event", y = "Density", colour = NULL)
    
    doy_plot <- ggplot(t_pred, aes(ymd("2025-01-01") + doy, .pred, colour = event_gr)) + 
      stat_summary(fun = 'mean', geom='line', aes(group = event_gr)) +
      rcartocolor::scale_colour_carto_d() +
      labs(x = NULL, y = "Density", colour = NULL)
    
    ggsave(plot = time_plot/doy_plot,
           file = g("{out_dir}/{spp}_time_date_marginal.jpeg"))
    
    
    vi_tbl <- vip::vi(fit_workflow)
    write_csv(vi_tbl, g("{out_dir}/{spp}_time_vi.csv"))
    
    vip_plot <- vip::vip(fit_workflow,geom = "point", n= 15)
    
    ggsave(plot = vip_plot, 
           filename = g("{out_dir}/{spp}_time_vi.jpeg"),
           width = 8, height = 6)
    
    explainer <-  fit_workflow %>% 
      DALEXtra::explain_tidymodels(
        data = df_std %>% select(-y),
        y = df_std$y,
      )
    write_rds(explainer, glue::glue("{bundle_locs}/{spp}_time_explainer.rds") )
    
    preds_org <- predict(fit_workflow, new_data = df_std) |> bind_cols(df_std)
    write_rds(preds_org,  g("{out_dir}/{spp}_time_model_pred.rds"))
    
    vars_ <-  vi_tbl |> slice_head(n=15) |> pull(Variable)
    
    pd <-  DALEX::model_profile(explainer = explainer,
                                variables = vars_[vars_ %in% names(explainer$data)]
    )
    pdp <- plot(pd)
    ggsave(pdp, 
           filename= g("{out_dir}/{spp}_time_pdp.jpeg"),
           width = 12, height = 12)
    
    
    resid_brt <- DALEX::model_performance(explainer)
    
    # create comparison plot of residuals for each model
    p1 <- plot(resid_brt)
    p2 <- plot(resid_brt, geom = "boxplot")
    
    
    resid_plt <- p1 + p2
    ggsave(plot = resid_plt,
           g("{out_dir}/{spp}_time_resid.jpeg"),
           width = 10, height = 6)
    
  })
  
}



run_brt_preds <- function(spp){
  # Predict ------------------------------
  out_dir_spatial <- g(
    "{BRT_output_loc_spatial}/{spp}"
  )
  out_dir <- g(
    "{BRT_output_loc}/{spp}"
  )
  tic("Run predictions")
  
  preds <- read_rds("{prediction_layer_loc}/2025-02-28_prediction_rasters_df.rds")
  xy <- read_rds("{prediction_layer_loc}/2025-02-28_prediction_rasters_xy.rds")
  
  fit_workflow <- #bundle::unbundle(read_rds(
    readRDS(
    glue::glue("{bundle_locs}/{spp}_time_bundle.rds"))
  
    p <- predict(fit_workflow, new_data = preds |>
                   mutate(doy = NA,
                          t2se = NA,
                          Time_period = NA,
                          year = NA, 
                          event_gr = "Dawn",
                          source = NA,
                          collection = NA,
                          type = NA,
                          Rec_length = factor(5),
                          SiteN=NA,
                          longitude = NA,
                          latitude = NA,
                          event_id = NA
                   ))
  toc()
  
  write_rds(p, glue::glue("{prediction_layer_loc}/{spp}_time_pred_brt.rds"))
  r_pred <- rast("{prediction_layer_loc}/2025-02-28_prediction_rasters.nc")#2025-02-10_prediction_rasters.nc")
  
  
  
  predicted_raster <- rast(r_pred[[1]])
  r2 <- terra::cellFromXY(predicted_raster,xy)
  predicted_raster[r2] <- p
  predicted_raster <- mask(predicted_raster, aoi)
  names(predicted_raster) <- spp
  
  
  
  terra::writeRaster(predicted_raster, glue::glue("{out_dir_spatial}/{spp}_time_brt.tif"),overwrite=T)

  
  
  min_q <-  min( p,  na.rm=T)

    mf_route <- read_sf(
    g("{mf_route_loc}/MFCAR_Route.shp")
  )
  
  mf_r <- crop(predicted_raster, st_bbox(st_transform(st_buffer(mf_route, 20000), st_crs(predicted_raster))))
  
  cairo_pdf(g("{out_dir}/{spp}_time_map_MFCAR.pdf"), fallback_resolution = 300, height = 17, width = 11)
  ggplot( ) + tidyterra::geom_spatraster(data = mf_r) +
    geom_sf(data = mf_route, alpha = 0.5, linetype =2) +
    tidyterra::scale_fill_whitebox_c(
      palette = "muted",#breaks = c(0.05, 0.1, 0.2, 0.3, 0.8),
      # labels = scales::label_number(suffix = "indiv/ha"),
      n.breaks = 8,limits= c(min_q, quantile(p, 0.999, na.rm=T)),
      guide = guide_legend(reverse = TRUE))  +
    labs(fill = "individuals/ha", title = spp,
         colour = "#\nNon-zero\ncounts") 
  dev.off()
  
  map_plot <-  ggplot() +
    tidyterra::geom_spatraster(data = predicted_raster) +
    labs(fill = "individuals/ha", title = spp,
         colour = "#\nNon-zero\ncounts") +
    tidyterra::scale_fill_whitebox_c(
      palette = "muted",#breaks = c(0.05, 0.1, 0.2, 0.3, 0.8),
      # labels = scales::label_number(suffix = "indiv/ha"),
      n.breaks = 8,limits= c(min_q, quantile(p, 0.999, na.rm=T)),
      guide = guide_legend(reverse = TRUE)) +
    # geom_sf(data = napken_lake, shape =2 )+
    # geom_sf(data = locs_in) +
    geom_sf(data = ra_area, fill = NA, linetype =2, colour = 'white') +
    geom_sf(data = mesh_data, fill = NA, linetype = 3, colour = 'black')
  
  ggsave(plot = map_plot, g("{out_dir}/{spp}_time_map.jpeg"), width = 6.73, height = 8.5)
  # hist(p[[1]][p[[1]]<0.004], breaks = 100)
  # hasValues(predicted_raster)
  
  
  
  
  cat(glue::glue("Predictions complete for {spp}\n"))
}





