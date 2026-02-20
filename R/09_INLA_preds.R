library(tidyverse)
library(sf)
library(fmesher)
library(patchwork)
source("R/__globals.R")
# QPAD::load_BAM_QPAD(4)
library(INLA)
library(inlabru)
library(tictoc)
library(tidymodels)

locs_in <- st_filter(individual_locs, ra_area)


prep_predictions <- function(spp, save_objects, return_all = FALSE) {
  spp_dir <- str_replace(spp, " ", "_")
  if (run_a2) {
    out_dir_app <- "A2"
  } else {
    out_dir_app <- ""
  }
  out_dir_spatial <- g(
    "{INLA_output_loc_spatial}/{out_dir_app}/{spp_dir}"
  )
  out_dir <- g(
    "{INLA_output_loc}/{out_dir_app}/{spp_dir}"
  )
  out_dir_tmp <- g(
    "{INLA_output_loc_TMP}/{out_dir_app}/{spp_dir}"
  )

  res <- read_rds(g("{out_dir_tmp}/{spp}_inlabru_model.rds"))
  fam <- res$.args$family

  pca_cov <- read_rds(g("{out_dir_tmp}/pca_cov.rds"))
  pca_prep <- read_rds(g("{out_dir_tmp}/pca_bundle.rds")) |>
    bundle::unbundle() |>
    recipes::prep()
  basic_rec <- read_rds(g("{out_dir_tmp}/basic_bundle.rds")) |> 
    bundle::unbundle() |> recipes::prep()
  pls_prep <- read_rds(g("{out_dir_tmp}/pls_bundle.rds")) |>
    bundle::unbundle() |>
    recipes::prep()

  print(g("Preping predictions, {spp} --------------------------"))

  pred_date <- ifelse(run_a2, "2025-01-10", "2025-10-22") #"2025-02-28")

  preds_r <- read_rds(g(
    "{prediction_layer_loc}/{pred_date}_prediction_rasters_df.rds"
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
      # Rec_length = factor(5),
      SiteN = NA,
      longitude = NA,
      latitude = NA,
      event_id = NA,
      QPAD_offset = NA,
      olcb_99_100 = NA,
      total_100 = NA,
      olcb_99_500 = NA,
      total_500 = NA,
      site_id = NA
    )
  # xy <- read_rds(g(
  #   "{prediction_layer_loc}/{pred_date}_prediction_rasters_xy.rds"
  # ))

  # names(df_std)[names(df_std) %in% names(preds)]
  pca_pred_data <- preds_r |>
    mutate(site_id = glue::glue("loc-{rn}")) |>
    dplyr::select(any_of(names(pca_cov)))
  pca_preds <- recipes::bake(pca_prep, new_data = pca_pred_data) 

  pls_preds <- recipes::bake(pls_prep, new_data = pca_pred_data) 
  included_vars <- read_rds(g("{out_dir_tmp}/included_variables.rds"))
  
  basic_preds <- bake(basic_rec, new_data = pca_pred_data)|> 
    dplyr::select(site_id,all_of(list_c(included_vars))) |> 
    mutate(rn = preds_r$rn)

  pc_vars <- str_subset(names(pca_preds), "^PC\\d")
  pls_vars <- str_subset(names(pls_preds), "^PLS\\d")

  list2env(read_rds(g("{out_dir_tmp}/sd_means.rds")), environment())
  if (!run_a2) {
    t2se_mod <- read_rds(g("{out_dir_tmp}/t2se_mod.rds"))
  }
  offsets_spp <- read_rds(g("{out_dir_tmp}/offsets_spp.rds"))

  f <-
    res$bru_info$model$formula |> as.character() |> pluck(3)

  if (str_detect(f, "_PLS")) {
    p_vars <- pls_vars
  } else if (str_detect(f, "_PC")) {
    p_vars <- pc_vars
  } else if (str_detect(f, "_[51]00")){
    p_vars <- list_c(included_vars)
  }else {
    p_vars <- "project"
  }

  scaled_prediction <-
    {if(str_detect(f, "_[51]00")){
      basic_preds}else{preds_r}} |>
    bind_cols(pca_preds |> 
                dplyr::select( -contains("is_forest"))) %>%
    bind_cols(pls_preds|> 
                dplyr::select( -contains("is_forest"))) %>%
    dplyr::select(rn, which(names(.) %in% p_vars)) |>
    pivot_longer(cols = -rn, names_to = 'variable', values_to = "x") |>
    left_join(sd_mn, by = join_by(variable)) |>
    mutate(xx = (x - mn) / sd2) |>
    dplyr::select(rn, variable, xx) |>
    pivot_wider(names_from = variable, values_from = xx, id_cols = rn) |>
    arrange(rn)

  preds_sc <- preds_r |>
    mutate(X_sc = X / 10000, Y_sc = Y / 10000) %>%
    # bind_cols(pca_preds) %>%
    dplyr::select(X, Y, !which(names(.) %in% p_vars)) |> #sd_mn$variable)) |>
    bind_cols(scaled_prediction |> dplyr::select(-rn, -any_of("site_id"))) |>
    mutate(recording_id = factor("None")) |>
    st_as_sf(coords = c("X", "Y"), crs = ont.proj)

  preds_sc$RL <- factor(5)
  # preds_sc$t2se_scaled <- t2se_mod
  if (!run_a2) {
    time_period_prediction <- slice_max(
      t2se_mod,
      order_by = median,
      with_ties = F
    ) |>
      pull(event_gr)
  } else {
    time_period_prediction <- "Dawn"
  }
  period_to_use <-
    preds_sc$t2se_sc <- (0 -
      t2se_sd_mn$mean_t2se[
        t2se_sd_mn$event_gr == time_period_prediction &
          !is.na(t2se_sd_mn$event_gr)
      ]) /
      (2 *
        t2se_sd_mn$sd_t2se[
          t2se_sd_mn$event_gr == time_period_prediction &
            !is.na(t2se_sd_mn$event_gr)
        ])
  # preds_sc$doy_r <- doy_mod
  preds_sc$doy_r <- (161 - sd_mn$mn[sd_mn$variable == "doy"]) /
    (sd_mn$sd2[sd_mn$variable == "doy"])
  preds_sc$event_gr <- time_period_prediction
  preds_sc$X <- preds_r$X
  preds_sc$Y <- preds_r$Y
  preds_sc$site <- 437:(436 + nrow(preds_sc))
  preds_sc$QPAD_offset <- offsets_spp |>
    distinct(dur, max_dist, o) |>
    filter(dur == 5) |>
    pull(o)
  # pres_filt <- filter(preds_sc,!is.na(fnlc_9_100))
  tictoc::tic()
  nsamp <- 200 # number of posterior samples

  # samples from mgcv model
  # Generate estimates from SPDE estimate only
  # kappa_var <- 1/res$summary.hyperpar$mean[3]

  n_splits <- 10

  pred_test <-
    preds_sc |>
    filter(!is.na(clc20_1_100)) |>
    dplyr::select(
      rn,
      all_of(p_vars),
      X,
      Y,
      t2se_sc,
      event_gr,
      RL,
      QPAD_offset,
      doy_r
    ) %>%
    mutate(
      #chunk = cut_interval(rn, n=n_splits),
      rn_test = row_number()
    )

  # r <- tidyterra::as_spatraster(bind_cols(xy,preds_sc |>
  #                                           dplyr::select(rn,all_of(p_vars), X,Y, t2se_sc, event_gr,RL,offset,
  #                                                                      doy_r)), xycols = 1:2)
  #
  # rr <- preds_sc |>
  #   dplyr::select(rn,all_of(p_vars), X,Y, t2se_sc, event_gr,RL,offset,
  #                 doy_r) |> rasterize(y = r_pred)

  # if(p_vars == "project") p_vars <- n
  ff <- paste0(
    "~ Intercept + ",
    paste0("spat_cov_", p_vars, collapse = "+"),
    "+ alpha"
  ) |>
    as.formula()

  ff_sp <- ~ Intercept + alpha
  ff_g <- glue::glue(
    '{res$misc$configs$contents$tag[-c(1,2)] |> glue::glue_collapse(sep = " + ")}'
  )
  linear_or_spde <- ifelse(str_detect(ff_g, "linear"), 'linear', 'spat')

  if (run_a2 | str_detect(ff_g, "PC01|PLS01", negate = T)) {
    ff_counts <- paste0(
      "Intercept +   t2se + doy +",#o +
      paste0(linear_or_spde, "_cov_", p_vars, collapse = "+"),
      "+ alpha"
    ) #|>
  } else {
    ff_counts <- paste0(
      "Intercept + RecLength  + t2se + time_group + doy +", #+ o
      paste0(linear_or_spde, "_cov_", p_vars, collapse = "+"),
      "+ alpha"
    ) #|>
  }

  ff_complex <- paste0(
    "~{",
    "expect <- exp(",
    ff_counts,
    ")
                       expect_sp <- exp(Intercept +  alpha)
                      list(expect = expect,",
    "pobs = 1-dpois(0, expect),
                      expect_sp = expect_sp)}"
  ) |>
    as.formula()

  ff_sp <- paste0(
    "~{
                       expect_sp <- exp(Intercept +  alpha)
                      list(expect_sp = expect_sp,
                      pobs_sp = 1-dpois(0, expect_sp) )}"
  ) |>
    as.formula()

  ff_zi <- paste0(
    "~{
                 scaling_prob <- (1 - zero_probability_parameter_for_zero_inflated_poisson_1)
                lambda <- exp(",
    ff_counts,
    ")
                lambda_sp <- exp(Intercept +  alpha)
                expect <- scaling_prob * lambda
                expect_sp <- scaling_prob * lambda_sp
                list(
                expect = expect,
                pobs = 1 - dpois(0, expect),
                expect_sp = expect_sp)
                }"
  ) |>
    as.formula()

  ff_zi_sp <- paste0(
    "~{
                 scaling_prob <- (1 - zero_probability_parameter_for_zero_inflated_poisson_1)
               lambda_sp <- exp(Intercept +  alpha)
                expect_sp <- scaling_prob * lambda_sp
                list(
                pobs_sp = 1 - dpois(0, expect_sp),
                expect_sp = expect_sp)
                }"
  ) |>
    as.formula()

  out_list <- list(
    pred_test = pred_test,
    ff_complex = ff_complex,
    ff_zi = ff_zi,
    ff_zi_sp = ff_zi_sp,
    fam = fam,
    ff_g = ff_g,
    ff_sp = ff_sp
  )

  if (isTRUE(save_objects)) {
    write_rds(out_list, g("{out_dir_tmp}/temp_ob_for_pred_{spp}.rds"))
  } else {
    if (isTRUE(return_all)) {
      return(as.list(environment()))
    } else {
      return(out_list)
    }
  }
}

run_predictions <- function(spp, load_rds = FALSE, gen_map_outputs = TRUE) {
  spp_dir <- str_replace(spp, " ", "_")
  if (run_a2) {
    out_dir_app <- "A2"
  } else {
    out_dir_app <- ""
  }
  out_dir_tmp <- g(
    "{INLA_output_loc_TMP}/{out_dir_app}/{spp_dir}"
  )
  out_dir_tmp_preds <- g(
    "{INLA_preds_loc_TMP}/{out_dir_app}/{spp_dir}"
  )
  out_rds_file <- g("{out_dir_tmp_preds}/{spp}_predictions_full.rds")
  if (file.exists(out_rds_file)) {
    mod_time <- file.info(out_rds_file)$mtime
    if (ymd_hms(mod_time) > ymd_hms("2025-6-12 06:00:00")) {
      return(NULL)
    }
  }
  dir.create(out_dir_tmp_preds, recursive = T)
  if (isTRUE(load_rds)) {
    list2env(
      read_rds(g("{out_dir_tmp}/temp_ob_for_pred_{spp}.rds")),
      environment()
    )
    res <- read_rds(g("{out_dir_tmp}/{spp}_inlabru_model.rds"))
  } else {
    list2env(
      prep_predictions(spp, save_objects = FALSE, return_all = FALSE),
      environment()
    )
    res <- read_rds(g("{out_dir_tmp}/{spp}_inlabru_model.rds"))
  }
  gc()
  nsamp <- 300
  print(g("Starting predictions, {spp} --------------------------"))
  if (str_detect(ff_g, "PC01|PLS01|_[15]00")) {
    gc()
    tic(g("Running predictions, {spp} --------------------------"))

    preds <- inlabru::generate(
      res,
      pred_test,
      switch(fam, 'poisson' = ff_complex, "zeroinflatedpoisson1" = ff_zi),
      num.threads = 1,
      n.samples = nsamp
    ) # %>% exp()

    toc()
  } else {
    gc()
    tic(g("Running predictions sp only, {spp} --------------------------"))
    preds <- inlabru::generate(
      res,
      pred_test,
      switch(fam, 'poisson' = ff_sp, "zeroinflatedpoisson1" = ff_zi_sp),
      num.threads = 1,
      n.samples = nsamp
    ) # %>% exp()
    toc()
  }

  if (!isTRUE(gen_map_outputs)) {
    write_rds(preds, out_rds_file)
  } else {
    gen_maps(spp, preds = preds)
  }
  if (isTRUE(load_rds)) {
    file.remove(g("{out_dir_tmp}/temp_ob_for_pred_{spp}.rds"))
  }
}

gen_maps <- function(spp, preds = NULL) {
  spp_dir <- str_replace(spp, " ", "_")
  if (run_a2) {
    out_dir_app <- "A2"
  } else {
    out_dir_app <- ""
  }
  out_dir_tmp <- g(
    "{INLA_output_loc_TMP}/{out_dir_app}/{spp_dir}"
  )
  out_dir_tmp_preds <- g(
    "{INLA_preds_loc_TMP}/{out_dir_app}/{spp_dir}"
  )

  out_dir_spatial <- g(
    "{INLA_output_loc_spatial}/{out_dir_app}/{spp_dir}"
  )
  out_dir <- g(
    "{INLA_output_loc}/{out_dir_app}/{spp_dir}"
  )
  tictoc::tic()
  preds_file <- g("{out_dir_tmp_preds}/{spp}_predictions_full.rds")
  if (is.null(preds)) {
    preds_file <- g("{out_dir_tmp_preds}/{spp}_predictions_full.rds")

    preds <- read_rds(preds_file)
  }

  preds_t <- preds |>
    transpose()

  if ("expect" %in% names(preds_t)) {
    pres_counts <- preds_t |>
      pluck("expect") %>%
      do.call("cbind", .)

    p_obs_str <- "pobs"
    # }) %>% do.call("rbind", .)
  } else {
    p_obs_str <- "pobs_sp"
  }

  # kappa <- res$summary.hyperpar[str_detect(rownames(res$summary.hyperpar), "kappa"),]
  # preds_sd <-  1/rnorm(ncol(pres_counts), kappa$mean, kappa$sd)
  # preds_with_kappa <- apply(pres_counts,MARGIN = 1, function(x) log(x) +
  #         rnorm(length(x), 0,
  #               1/rnorm(1, kappa$mean, kappa$sd)) )
  #

  pres_counts_sp <- preds_t |>
    pluck("expect_sp") %>%
    do.call("cbind", .)

  cat(glue::glue("{spp}---------------------\n"))
  if (p_obs_str %in% names(preds_t)) {
    pres_obs <-
      preds_t |>
      pluck(p_obs_str) %>%
      do.call("cbind", .)
  } else {
    if (exists("pres_counts")) {
      pres_obs <- 1 - dpois(0, pres_counts)
    } else {
      pres_obs <- 1 - dpois(0, pres_counts_sp)
    }
  }
  tictoc::toc()

  rm(preds_t, preds)

  pred_date <- ifelse(run_a2, "2025-01-10", "2025-10-22") #"2025-02-28")

  # xy <- read_rds(g(
  #   "{prediction_layer_loc}/{pred_date}_prediction_rasters_xy.rds"
  # ))

  preds_r <- read_rds(g(
    "{prediction_layer_loc}/{pred_date}_prediction_rasters_df.rds"
  )) |> #2025-02-11_prediction_rasters_df.rds") |>
    mutate(rn = row_number()) |>
    dplyr::select(rn, clc20_1_100, X, Y)

  pred_test <-
    preds_r |>
    filter(!is.na(clc20_1_100))

  # prob obs -----
  pred_test$p_obs <- matrixStats::rowMedians(pres_obs, na.rm = T, digits = 4L) #apply((pres_obs), 1, mean)
  pred_test[, c("lci_p", "uci_p")] <- matrixStats::rowQuantiles(
    pres_obs,
    probs = c(.025, .975),
    na.rm = T
  ) #apply(modpred, 1, sd)

  # SP ------
  # pred_test$median_sp <- apply(pres_counts_sp, 1, median)
  pred_test$mean_sp <- matrixStats::rowMeans2(
    pres_counts_sp,
    na.rm = T,
    digits = 4L
  ) #apply(pres_counts_sp, 1, mean)
  pred_test$sd_sp <- matrixStats::rowSds(pres_counts_sp, na.rm = T, digits = 4L) # apply(pres_counts_sp, 1, sd)
  pred_test$mad_sp <- matrixStats::rowMads(
    pres_counts_sp,
    na.rm = T,
    digits = 4L
  ) #apply(pres_counts, 1, sd)
  pred_test[, c(
    "lci_count_sp",
    "median_sp",
    "uci_count_sp"
  )] <- matrixStats::rowQuantiles(
    pres_counts_sp,
    probs = c(.025, 0.5, .975),
    na.rm = T,
    digits = 4L
  ) #apply(modpred, 1, sd)
  # preds_sc$spatial_mean[!is.na(preds_sc$fnlc_10_500)] <- apply(modpred_xy, 1, mean)
  # preds_sc$spatial_sd[!is.na(preds_sc$fnlc_10_500)] <- apply(modpred_xy, 1, sd)
  # preds_sc$sd[!is.na(preds_sc$fnlc_10_500)] <-matrixStats::rowSds(x = M_predictions) #apply(modpred, 1, sd)
  # pred_test$mean_estimate[pred_test$mean_estimate>quantile(pred_test$mean_estimate, probs = 0.99)] <- NA

  # Counts ------
  if (exists("pres_counts")) {
    pred_test$mean_count <- matrixStats::rowMeans2(
      pres_counts,
      na.rm = T,
      digits = 4L
    ) #apply((pres_counts), 1, mean)
    # pred_test$mean_count_kappa <- matrixStats::rowMeans2(t(preds_with_kappa) , na.rm=T, digits = 4L) #apply((pres_counts), 1, mean)
    # pred_test$median_count <- apply(pres_counts, 1, median)
    pred_test$sd_count <- matrixStats::rowSds(
      pres_counts,
      na.rm = T,
      digits = 4L
    ) #apply(pres_counts, 1, sd)
    # pred_test$sd_count_kappa <- matrixStats::rowSds(t(preds_with_kappa), na.rm=T, digits = 4L)#apply(pres_counts, 1, sd)
    pred_test$mad_count <- matrixStats::rowMads(
      pres_counts,
      na.rm = T,
      digits = 4L
    ) #apply(pres_counts, 1, sd)

    pred_test[, c(
      "lci_count",
      "median_count",
      "uci_count"
    )] <- matrixStats::rowQuantiles(
      pres_counts,
      probs = c(.025, 0.5, .975),
      na.rm = T,
      digits = 4L
    ) #apply(modpred, 1, sd)

    # pred_test[, c("lci_sp", "uci_sp")] <-matrixStats::rowQuantiles(pred_sp, probs=c(.025, .975)) #apply(modpred, 1, sd)
    pred_test$ci_size_count <- pred_test$uci_count - pred_test$lci_count
  }

  # preds_sc$mean_estimate[!is.na(preds_sc$fnlc_10_500)] <-matrixStats::rowMeans2(M_predictions)# apply(modpred, 1, mean)
  # preds_sc$simulation[!is.na(preds_sc$fnlc_10_500)] <- apply(simulated_density, 1, median)
  # p_avg$sd <- p_sd

  # full_p <- left_join(preds_sc, p_avg)

  preds_sc <- left_join(
    preds_r |> select(rn),
    pred_test |>
      st_drop_geometry() |>
      dplyr::select(
        rn,
        tidyselect::any_of(c(
          "mean_count",
          "p_obs",
          "median_count",
          "median_sp",
          "mean_sp",
          "sd_sp",
          "sd_count",
          "mad_count",
          "mad_sp",
          "lci_count",
          "uci_count",
          "lci_p",
          "uci_p",
          "lci_count_sp",
          "uci_count_sp",
          "ci_size_count"
        ))
      ),
    by = join_by(rn)
  ) %>%
    {
      if (exists("pres_counts")) {
        mutate(
          .,
          cv_count = sd_count / mean_count,
          k_cv = sqrt(cv_count^2 / (1 + cv_count^2))
        )
      } else {
        .
      }
    } %>%
    {
      if (exists("pres_counts_sp")) {
        mutate(
          .,
          cv_spatial = mean_sp / sd_sp,
          k_cv_sp = sqrt(cv_spatial^2 / (1 + cv_spatial^2))
        )
      } else {
        .
      }
    }

  # k_cv_spatial = sqrt(cv_spatial^2/(1+cv_spatial^2)))#,
  # est_mod = ifelse(est>=log(0.25), log(0.25), est))

  rm(pres_counts, pres_counts_sp, pres_obs)
  # ggplot(p, aes(X_sc, Y_sc, colour =est )) + geom_point(alpha = 0.4)
  r_pred <- rast(glue::glue(
    "{prediction_layer_loc}/2025-10-21_prediction_rasters.nc"
  ))

  r_temp <- rast(r_pred[[1]]) #g("{pred_date}_prediction_rasters_1")]])

  in_mesh <- fmesher::fm_is_within(
    st_as_sf(preds_r[, c("X", "Y")], coords = c("X", "Y"), crs = 3161),
    mesh_inla
  )

  r2 <- terra::cellFromXY(r_temp, preds_r[, c("X", "Y")])
  # r2 <- terra::cellFromXY(predicted_raster, preds[, c("X", "Y")])
  # median_count <- #median_sp <-
  #   sd_r <- lci_r <- lci_p <- #lci_sp <- uci_sp <-
  #   uci_p<- uci_r <- mean_estimate <-
  #   p_obs <- cv_r <-  k_cv <-  ci_size_r <- r_temp
  names_outstack <- names(preds_sc)[-1]
  for (i in 1:length(names_outstack)) {
    j <- names_outstack[[i]]
    # assign(i, r_temp)
    print(j)
    tt <- r_temp
    tt[r2] <- preds_sc[[j]][in_mesh]

    if (i == 1) {
      outstack <- tt
    } else {
      add(outstack) <- tt
    }
    names(outstack)[[i]] <- j
    varnames(outstack)[[i]] <- j
    # assign(i, tt)

    rm(tt)
  }

  # sd_r[r2] <- preds_sc$sd_count[in_mesh]
  # lci_r[r2] <- preds_sc$lci_count[in_mesh]
  # lci_p[r2] <- preds_sc$lci_p[in_mesh]
  # # lci_sp[r2] <- preds_sc$lci_sp[in_mesh]
  # # uci_sp[r2] <- preds_sc$uci_sp[in_mesh]
  # uci_r[r2] <- preds_sc$uci_count[in_mesh]
  # uci_p[r2] <- preds_sc$uci_p[in_mesh]
  # mean_estimate[r2] <- preds_sc$mean_count[in_mesh]
  # median_count[r2] <- preds_sc$median_count[in_mesh]
  # # median_sp[r2] <- preds_sc$median_sp[in_mesh]
  # p_obs[r2] <- preds_sc$p_obs[in_mesh]
  # ci_size_r[r2] <- preds_sc$ci_size_count[in_mesh]
  # cv_r[r2] <- preds_sc$cv_count[in_mesh]
  # # cv_spatial[r2] <- preds_sc$cv_spatial
  # k_cv[r2] <- preds_sc$k_cv[in_mesh]
  # # k_cv_spatial[r2] <- preds_sc$k_cv_spatial

  #
  # outstack <- c(sd_r,
  #               lci_r,lci_p,
  #               # uci_sp,lci_sp,
  #               uci_r,uci_p,
  #               ci_size_r,
  #               mean_estimate,
  #               median_count,
  #               # median_sp,
  #               p_obs,
  #               cv_r,
  #               # cv_spatial,
  #               k_cv)
  # names(outstack) <- varnames(outstack) <- c(
  #   "sd_count", "lci_count", "lci_p",
  #   # "uci_sp", "lci_sp",
  #   "uci_count", "uci_p",
  #   "ci_size_count",
  #   "mean_count","median_count",
  #   # "median_sp",
  #   "p_obs",
  #   "cv_count", #"cv_spatial_only",
  #   "k_cv")

  # terra::writeCDF(outstack, filename = glue::glue("{out_dir_spatial}/Predictions_{spp}.nc"),overwrite=T, split=T)
browser()
  for (i in names_outstack) {
    terra::writeRaster(
      outstack[[i]],
      glue::glue("{out_dir_spatial}/{i}_{spp}.tif"),
      overwrite = T
    )
  }

  expectations <- str_subset(names_outstack, "mean|median|p_obs")
  errors <- str_subset(names_outstack, "mean|median|p_obs", negate = T)
  blob <- read_sf(blob_loc) |>
    st_transform(ont.proj)
  # locs_in <- locations[fmesher::fm_is_within(locations, mesh_inla),]

  plot_map <- function(i,  type, max_q = 0.999) {
    min_q <- 0
    if (str_detect(i, "_sp")) {
      min_q <- min(preds_sc[[i]], na.rm = T)
    }
    # print(
    #   ggplot() +
    #     tidyterra::geom_spatraster(data = outstack[[i]], maxcell = 1e6) +
    #     labs(fill = i, title = spp) +
    #     tidyterra::scale_fill_whitebox_c(
    #       palette = "muted", #breaks = c(0.05, 0.1, 0.2, 0.3, 0.8),
    #       # labels = scales::label_number(suffix = "indiv/ha"),
    #       n.breaks = 8,
    #       limits = c(min_q, quantile(preds_sc[[i]], 0.999, na.rm = T)),
    #       guide = guide_legend(reverse = TRUE)
    #     ) +
    #     geom_sf(data = napken_lake, shape = 2) +
    #     # geom_sf(data = locs_in) +
    #     geom_sf(data = ra_area, fill = NA, linetype = 2, colour = 'white') +
    #     geom_sf(data = mesh_data, fill = NA, linetype = 3, colour = 'black')
    pal <-   switch (type,
                     "expectation" = tidyterra::scale_fill_wiki_c(
                       limits =  c(min_q, quantile(preds_sc[[i]], max_q, na.rm = T)),
                       guide = guide_legend(reverse = TRUE)
                     ),
                     "uncertainty" =  tidyterra::scale_fill_princess_c(
                       palette = 'aura',
                       limits =  c(min_q, quantile(preds_sc[[i]], max_q, na.rm = T))) 
    )
    
      print(
        ggplot() +
        tidyterra::geom_spatraster(
          data = mask(outstack[[i]], ra_area),
          maxcell = 1e6
        ) +
        labs(fill = i, title = spp) + pal+
        geom_sf(data = blob, fill = 'red', alpha = 0.4, colour = NA) +
        ggthemes::theme_map()
      
    )
  }
  maps_ex <- g("{out_dir}/maps_expectations_{spp}.pdf")
  maps_uc <- g("{out_dir}/maps_uncertainty_{spp}.pdf")
  if (file.exists(maps_ex)) {
    file.remove(maps_ex)
  }
  if (file.exists(maps_uc)) {
    file.remove(maps_uc)
  }
  cairo_pdf(maps_ex)
  map(expectations, plot_map, type = 'expectation')
  dev.off()

  cairo_pdf(maps_uc)
  map(errors, plot_map, type = 'uncertainty')
  dev.off()
  toc()
  if (exists("preds_file")) file.remove(preds_file)
}
