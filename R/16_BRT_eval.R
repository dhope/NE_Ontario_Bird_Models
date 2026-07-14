run_a2 <- FALSE
source("R/06_BRT_time.R")
source("R/__globals.R")
spp_cov_date <- "2026-02-16"
spp_cov_file <- g("output/rds/{spp_cov_date}_spatial_covariates_data.rds")
source("R/__data_prep.R")
print("Data prep complete--------------------------------")
use_training_data <- FALSE
if (!use_training_data) {

  
  extra_tests <- read_rds(g("{rds_data_loc}/../all_events.rds") )|> 
    filter(str_detect(project, "Crescent"))
  
  extra_counts <- read_rds(g("{rds_data_loc}/../counts.rds") )|> 
    filter(str_detect(project, "Crescent"))
  
  extra_locs <- locs_extra <- read_rds(g("{rds_data_loc}/../locations.rds") )|> 
    filter(str_detect(project, "Crescent"))
  
  
  
  extra_spat_var <- 
    read_rds("output/rds/2026-04-08_spatial_covariates_data_extra.rds") |> 
    mutate( d2O = read_rds(g("output/rds/2026-04-08_dist2ocean_extra.rds"))) |> 
    bind_cols(dplyr::select(
      extra_locs,
      site_id, geometry) ) |> 
    distinct() |>
    st_as_sf() %>%
    bind_cols(as_tibble(st_coordinates(.))) |>
    dplyr::select(where(~ sum(is.na(.x)) != length(.x))) 
  spatial_cov <- bind_rows(spatial_cov, extra_spat_var)
  
  
  
  test_training_data <- read_rds(g(
    "output/{date_compiled}_test_training_data.rds"
  ))
  test_locations <- test_training_data$test_locations
  test_recordings <- test_training_data$test_recordings
  train_ids <- test_training_data$train_recordings$event_id 
  rm(
    test_training_data,
    mesh_il,
    mesh_inla,
    mesh_inla_full,
    locs,
    ll_ag,
    locs_agg,
    locs_neighbours,
    run_brt,
    run_brt_preds
  )
  recordings <-
    g("{rds_data_loc}/all_events.rds") |>
    read_rds() |>
    filter(!event_id %in% train_ids) |>
    filter(collection %in% c("WildTrax", "ONATLAS3PC") &
             str_detect(location, "WHA-", negate=T) &
             str_detect(project, "Nocturnal", negate = T) ) |>
    
    # left_join(aggregated_locs,by = join_by(project, location, collection)) |>
    filter(str_detect(project, "(Extraction)|(Nocturnal)|(Resample)",negate=T) &
             (site_id %in% test_locations$site_id | event_id %in%
                test_recordings$event_id) ) |>
    bind_rows(extra_tests) |> 
    mutate(
      test_group = ifelse(site_id %in% test_locations$site_id, "Spatial", "Site"),
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
      Rec_length = factor(as.numeric(round(dur,1))),
      recording_id = as.numeric(recording_id)) %>%
  left_join(.,
  individual_locs |> filter(site_id %in% .$site_id) |>
    bind_rows(extra_locs |>st_as_sf() |>  distinct(site_id,geometry) ) |> 
    st_join(ontario_ez |> dplyr::select(on_er = ZONE_NAME)    ) ,
  by = join_by(site_id))
  
  # # rm(test_training_data, train_locs)
  # 
  counts <- read_rds(g("{rds_data_loc}/counts.rds")) |>
    filter(event_id %in% recordings$event_id ) |>
    filter(str_detect(project, "(Extraction)|(Nocturnal)|(Resample)",negate=T)) |>
    bind_rows(extra_counts) |> 
    dplyr::select(
      event_id,
      location,
      project,
      common_id,
      species_name_clean,
      total_count,
      species_scientific_name,
      species_code,
      total_count_with_tmtt
    ) |>
    replace_na(list(collection="WildTrax")) |>
   
    # dplyr::select(event_id, site_id,  project, 
    #               species_name_clean,total_count,
    #               common_id,
    #               total_count_with_tmtt) |>
    mutate(y = ifelse(is.na(total_count_with_tmtt), total_count,
                      total_count_with_tmtt)) |>
    filter(!is.na(y))
  # 
  # recordings <-
  #   g("{rds_data_loc}/all_events.rds") |>
  #   read_rds() |>
  #   # filter(location %in% aggregated_locs$location) |>
  #   filter(
  #     collection %in%
  #       c("WildTrax", "ONATLAS3PC") &
  #       str_detect(location, "WHA-", negate = T) &
  #       str_detect(project, "Nocturnal", negate = T)
  #   ) |>
  # 
  #   # left_join(aggregated_locs, by = join_by(project, location, collection)) |>
  #   filter(
  #     str_detect(project, "(Extraction)|(Nocturnal)|(Resample)", negate = T) &
  #       (site_id %in%
  #         test_locations$site_id |
  #         event_id %in%
  #           test_recordings$event_id)
  #   ) |>
  #   mutate(
  #     test_group = ifelse(
  #       site_id %in% test_locations$site_id,
  #       "Spatial",
  #       "Site"
  #     ),
  #     Time_period = dplyr::case_when(
  #       t2ss >= -60 & t2ss <= 150 ~ "Dusk",
  #       t2sr >= -70 & t2sr <= 220 ~ "Dawn",
  #       t2ss > 150 & t2sr < -70 ~ "Night",
  #       abs(t2sr) > 220 ~ "Day",
  #       is.na(t2sr) ~ "Missing time or location data",
  #       TRUE ~ "Unk"
  #     ),
  #     dur = round(clip_length_min, 2),
  #     t2se = dplyr::case_when(
  #       Time_period == "Dusk" ~ t2ss,
  #       Time_period == "Dawn" ~ t2sr,
  #       TRUE ~ pmin(abs(t2ss), abs(t2sr))
  #     ),
  #     doy = yday(date),
  #     recording_id = as.numeric(recording_id)
  #   ) #|>
  # 
  # # rm(test_training_data, train_locs)
  # 
  # counts <- read_rds(g("{rds_data_loc}/counts.rds")) |>
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
  #     common_id,
  #     total_count_with_tmtt
  #   ) |>
  #   mutate(
  #     y = ifelse(
  #       is.na(total_count_with_tmtt),
  #       total_count,
  #       total_count_with_tmtt
  #     )
  #   ) |>
  #   filter(!is.na(y))
  print("Loaded extra sites--------------------------------")
}
test_out_of_sample <- function(spp, write_to_file) {
  tic(g("{spp}:"))
  spp_dir <- str_replace_all(spp, " ", "_")
  out_dir_spatial <- g(
    "{BRT_output_loc_spatial}/{spp_dir}"
  )
  out_dir <- g(
    "{BRT_output_loc}/{spp_dir}"
  )
  # tic("Run predictions")

  list2env(prep_brt_data(spp), environment())
  if (use_training_data) {
    df_std <- read_rds(g("{brt_spp_dat_loc}/{spp_dir}.rds"))
    df_std$test_group <- "Train"
  }

  fit_workflow <- readRDS(glue::glue("{bundle_locs}/{spp_dir}_time_bundle.rds")) #)
  # bundle::unbundle(

  gen_tests <- function(test_group_) {
    if (!use_training_data) {
      test_dat <- df_std[
        df_std$location %in%
          (recordings$location[recordings$test_group == test_group_]),
      ]
    } else {
      test_dat <- df_std |> mutate(test_group = "Train")
    }
    test_predictions <- predict(fit_workflow, new_data = test_dat)

    test_res <- bind_cols(
      test_predictions,
      test_dat %>% select(y, QPAD_offset, test_group)
    ) |>
      mutate(
        expected_count = exp(log(.pred) - QPAD_offset),
        .resid = y - .pred,
        dpois = dpois(y, .pred),
        dpois_obs = 1 - dpois(0, .pred),
        obs_yn = factor(
          ifelse(y == 0, "Not Observed", "Observed"),
          levels = c("Observed", "Not Observed")
        )
      )
    # test_res

    test_lm <- lm(.pred ~ y, data = test_res, na.action = "na.exclude")

    test_metrics <- metric_set(rmse, rsq, mae, poisson_log_loss, ccc)
    get_eval_metrics <- function(df, spp) {
      tg <- df$test_group[[1]]
      df |>
        mutate(
          diff = abs(y - .pred)
        ) |>
        summarize(
          accuracy = mean(diff) / mean(y),
          precision = sd(y) / sd(.pred),
          cor_spearman = cor(.pred, y, method = "spearman"),
          cor_pearson = cor(.pred, y, , method = "pearson"),
          deviance = 2 *
            sum(
              ifelse(y == 0, 0, (y * log(y / .pred))) -
                (y - .pred)
            ),
          mean_residual = mean(abs(.resid)),
          discrim_intercept = test_lm$coefficients[1],
          discrim_slope = test_lm$coefficients[2],
        ) |>
        bind_cols(
          bind_rows(
            test_metrics(test_res, truth = y, estimate = .pred),
            roc_auc(test_res, obs_yn, dpois_obs)
          ) |>
            dplyr::select(-.estimator) |>
            pivot_wider(names_from = .metric, values_from = .estimate)
        ) |>
        mutate(species = spp, test = tg)
    }

    test_metrics <- get_eval_metrics(df = test_res, spp)

    if (sum(test_res$y) > 0) {
      test_roc_curve <-
        roc_curve(test_res, obs_yn, dpois_obs) |>
        ggplot(aes(x = 1 - specificity, y = sensitivity)) +
        geom_path() +
        geom_abline(lty = 3) +
        coord_equal() +
        theme_light() +
        labs(title = glue::glue("{spp} - {test_group_}"))

      # roc_auc(spatial_test_res, obs_yn, dpois_obs)
      test_predictions_plot <-
        ggplot(test_res, aes(x = y, y = .pred)) +
        # Create a diagonal line:
        geom_abline(lty = 2) +
        geom_point(alpha = 0.25) +
        geom_smooth(method = 'lm') +
        labs(
          y = "Predicted count",
          x = "Observed Count",
          title = glue::glue("{spp} - {test_group_}")
        ) +
        # Scale and size the x- and y-axis uniformly:
        coord_obs_pred()

      p_obs <-
        ggplot(test_res, aes(x = y, y = dpois_obs)) +
        ggdist::stat_dist_halfeye() +
        labs(
          y = "Probability of observing a bird",
          x = "Observed Count",
          title = glue::glue("{spp} - {test_group_}")
        )
    } else {
      p_obs <- NULL
      test_predictions_plot <- NULL
      test_roc_curve <- NULL
    }

    return(list(
      test_pred = test_predictions,
      test_metrics = test_metrics,
      plot_roc = test_roc_curve,
      plot_predictions = test_predictions_plot,
      plot_pobs = p_obs
    ))
  }

  output_test <- map(unique(df_std$test_group), gen_tests)
  names(output_test) <- unique(df_std$test_group)

  tr_output <- transpose(output_test)

  if (write_to_file) {
    write_csv(
      bind_rows(
        tr_output$test_metrics
      ),
      g("{out_dir}/{ifelse(use_training_data, 'train_', '')}test_metrics.csv")
    )

    write_csv(
      bind_rows(
        tr_output$test_pred,
      ),
      g(
        "{out_dir}/{ifelse(use_training_data, 'train_', '')}test_predictions.csv"
      )
    )

    withr::with_package('patchwork', {
      # browser()
      if (!all(list_c(map(tr_output$plot_roc, is.null)))) {
        try({
          ggsave(
            g(
              "{out_dir}/{ifelse(use_training_data, 'train_', '')}test_roc.jpeg"
            ),
            patchwork::wrap_plots(tr_output$plot_roc, nrow = 1),
            width = 12,
            height = 7
          )
        })

        try({
          ggsave(
            g(
              "{out_dir}/{ifelse(use_training_data, 'train_', '')}test_predictions.jpeg"
            ),
            patchwork::wrap_plots(tr_output$plot_predictions, nrow = 1),
            width = 12,
            height = 7
          )
        })

        try({
          ggsave(
            g(
              "{out_dir}/{ifelse(use_training_data, 'train_', '')}test_predictions_pobs.jpeg"
            ),
            patchwork::wrap_plots(tr_output$plot_pobs, nrow = 1),
            width = 12,
            height = 7
          )
        })
      }
    })
  }
  toc()
  if (!write_to_file) {
    return(output_test)
  }
}

rm_spp_pat <- c("Gull", "Loon", "Mallard", "Swan", "Harrier", "Eagle", "Teal", "Pintail", "Goldeneye", "Merganser", "Duck", "Merlin",
                "Scoter", "Tern", "Wigeon", "Gadwall") # Already ran these in testing
#
spp_to_run <- tibble(
                     species_name_clean=
                       bundle_locs |> 
                       list.files(full.names = F, recursive = T,
                                  "_time_bundle.rds") |>
                       str_remove("_time_bundle.rds") |> 
                       str_replace_all("_", " ") ) |> 
  # arrange(file_size) |>
  filter(str_detect(species_name_clean, 
                    glue::glue_collapse(rm_spp_pat, sep = "|"),
                    negate=T)
         & str_detect(species_name_clean, "^\\w{4}$", negate=T))

comp_metrics <- map(str_replace_all(spp_to_run$species_name_clean, "\\s", "_"), 
                    ~{
                      file.exists(g("{BRT_output_loc}/{.x}/{ifelse(use_training_data,'train_','')}test_metrics.csv"))
                    } ) |> list_c()
missing <- spp_to_run$species_name_clean[!comp_metrics]

list_modes <- list.files(bundle_locs, pattern = "\\w{4}_time_bundle.rds") |>
  str_extract("\\w{4}(?=_)")

evaluation_results <- map(
  missing,
  # spp_to_run$species_name_clean[-c(1,2)], #[(51+21):length(list_modes)],
  # str_subset(list_modes, "PAWA|ATSP", negate = T),
  safely(test_out_of_sample),
  write_to_file = T
) # ATSP needs to be run again
# e <- evaluation_results |> transpose() |> pluck( "error", .default = "No error")
# ii <- map_dbl(e, ~ifelse(is.null(.x), 0,1))
# spp_to_run$species_name_clean[ii==1]

metrics <-
  # map(
    list.files(
      BRT_output_loc,
      "test_metrics.csv",
      full.names = T,
      recursive = T
    ) |> stringr::str_subset("\\s", negate=T) |>
    readr::read_csv()
# 
readr::write_csv(metrics, "output/2026-04-30_BRT_metrics_combined.csv")



  