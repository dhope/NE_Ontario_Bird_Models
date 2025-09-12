source("R/08_INLA_SPDE.R")


test_training_data <- read_rds(g("output/{date_compiled}_test_training_data.rds"))
test_locations <- test_training_data$test_locations
test_recordings <- test_training_data$test_recordings
rm(test_training_data, mesh_il, mesh_inla, mesh_inla_full, locs, ll_ag, locs_agg, locs_neighbours, run_brt, run_brt_preds)

recordings <- 
  g("{rds_data_loc}/all_events.rds") |> 
  read_rds() |> 
  filter(location %in% aggregated_locs$location) |> 
  filter(collection %in% c("WildTrax", "ONATLAS3PC") &
           str_detect(location, "WHA-", negate=T) & 
           str_detect(project, "Nocturnal", negate = T) ) |> 
  
  left_join(aggregated_locs,by = join_by(project, location, collection)) |> 
  filter(str_detect(project, "(Extraction)|(Nocturnal)|(Resample)",negate=T) &
           (site_id_agg %in% test_locations$site_id_agg | event_id %in%
              test_recordings$event_id) ) |>
  mutate(
    test_group = ifelse(site_id_agg %in% test_locations$site_id_agg, "Spatial", "Site"),
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
    recording_id = as.numeric(recording_id)) #|> 

# rm(test_training_data, train_locs)

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


inla_tests <- function(spp){
  tictoc::tic(g('{spp}'))
  out_dir_spatial <- g(
    "{INLA_output_loc_spatial}/{out_dir_app}/{spp}"
  )
  out_dir <- g(
    "{INLA_output_loc}/{out_dir_app}/{spp}"
  )
  out_dir_tmp <- g(
    "{INLA_output_loc_TMP}/{out_dir_app}/{spp}"
  )
  
  list2env( prep_inla_data(spp, run_a2 = run_a2), envir = environment())
  
  
  

  
  
  res <- read_rds(  g("{out_dir_tmp}/{spp}_inlabru_model.rds")) 
  df_org <- read_rds(g("{out_dir_tmp}/{spp}_inlabru_model_data.rds")) 
  
  fam <- res$.args$family
  
  pca_cov <- read_rds(g("{out_dir_tmp}/pca_cov.rds"))
  pca_prep <- read_rds(g("{out_dir_tmp}/pca_bundle.rds")) |> bundle::unbundle() |> 
    recipes::prep()
  pls_prep <- read_rds(g("{out_dir_tmp}/pls_bundle.rds")) |> bundle::unbundle() |> 
    recipes::prep()
  
  
  pca_cov_test <- (spatial_cov) |> 
    filter(location %in% setup_dat$location) |> 
    dplyr::select(location,where(~{!is.factor(.x) & !is.character(.x)}),
                  -site_id_agg, -X, -Y) |> 
    st_drop_geometry()
  
  # names(df_std)[names(df_std) %in% names(preds)]
  pca_pred_data <- pca_cov_test |> 
    # mutate(location = glue::glue("loc-{1:nrow(pca_cov_test)}")) |>
    dplyr::select( any_of(names(pca_cov)) )
  
  pca_preds <- recipes::bake(pca_prep, new_data = pca_pred_data) #|>
  pls_preds <- recipes::bake(pls_prep, new_data = pca_pred_data) |> 
    dplyr::select(-location)
  
  
    # dplyr::select(-location)
  pc_vars <- str_subset(names(pca_preds), "^PC\\d")
  pls_vars <- str_subset(names(pls_preds), "^PLS\\d")
  
  list2env(read_rds(g("{out_dir_tmp}/sd_means.rds")), environment())
  t2se_sd_mn <- filter(t2se_sd_mn,!is.na(event_gr))
  t2se_mod <- read_rds(g("{out_dir_tmp}/t2se_mod.rds"))
  offsets_spp <- read_rds(g("{out_dir_tmp}/offsets_spp.rds"))
  
  ff_g <- glue::glue('{res$misc$configs$contents$tag[-c(1,2)] |> glue::glue_collapse(sep = " + ")}')
  # ff_g <- glue::glue('{res$misc$configs$contents$tag[-c(1,2)] |>
  #  str_subset("kappa", negate = T) |> 
  #                    glue::glue_collapse(sep = " + ")}')
  
  
  # sum(!is.na(pca_preds))
  if(str_detect(ff_g,"_PLS")) {
    p_vars <- pls_vars
  }else if(str_detect(ff_g,"_PC")){
    p_vars <- pc_vars
  } else{p_vars <- "project"} 
  
  
  
  scaled_prediction <- 
    setup_dat |> 
    left_join(pca_preds |> bind_cols(pls_preds),
              by = join_by(location)) %>% 
    dplyr::select( location, doy,Rec_length,test_group,
                   event_gr,t2se,site_id_agg = site_id_agg.x,
                   X, Y,offset,
                  which(names(. )%in% p_vars)) %>%
    bind_cols(sd_mn |> 
                filter(variable %in% names(.) & variable !="t2se") |> 
                pivot_wider(names_from = variable,
                            values_from = c(sd2, mn))) |> 
    mutate(across(c(doy, starts_with("PC|PLS")),
           ~{(.x - get(glue::glue("mn_{cur_column()}")) )/
                get(glue::glue("sd2_{cur_column()}") )}
           ) ) |> 
   
    left_join(t2se_sd_mn,by = join_by(event_gr)) |> 
    mutate(t2se_sc = (t2se- mean_t2se)/(2*sd_t2se),
           doy_r = doy) |> 
    dplyr::select(-starts_with("sd"),
                  -starts_with("mn_"),
                  -starts_with("mean_")) |> 
    mutate(X_sc = X/10000,
           Y_sc = Y/10000,
           RL = Rec_length,
           y = setup_dat$y,
           site = factor(site_id_agg)) |> 
    filter(event_gr %in% t2se_sd_mn$event_gr) 
    
    
    
  

  # tictoc::tic()
  nsamp <- 200 # number of posterior samples
  
  # samples from mgcv model 
  # Generate estimates from SPDE estimate only
  # kappa_var <- 1/res$summary.hyperpar$mean[3]
  
  n_splits <- 10
  
  
  
  
  

  
  
 
  
  ff <- as.formula(
    paste0("~o + ", ff_g)
  )
  # nsamp <- 500
  gen_inla_test <- function(test_group_){
    tictoc:::tic(test_group_)
    if(test_group_=="Train"){
      d_test <- df_org |> 
        mutate(rn = as.character(row_number()),
               test_group = "Train")
    }else{
  d_test <- scaled_prediction[scaled_prediction$test_group==test_group_,] |> 
    mutate(rn = as.character(row_number()))
    }
    
    
  if(fam!="zeroinflatedpoisson1"){
  preds <- inlabru::generate(res, newdata=d_test,formula = ff,
                             num.threads =32,
                             n.samples = nsamp) %>% exp() 
  d_pois <- map(1:nrow(d_test),~{dpois(d_test$y[[.x]], (preds[.x,]) 
  )}) %>%
    do.call(rbind, .)
  
  test_log_scores <- -log(rowMeans(d_pois) )
  # p_zero <- dpois(0,exp(preds))
  
  preds_out <- 
    list(expect = (preds),
         prob_zero = dpois(0,(preds)),
         obs_prob = d_pois
         )
  
  } else{
    # if(str_detect(ff_g,"PC01|PLS01")){
    preds <- inlabru::generate(res, newdata=d_test,
                as.formula(paste0("~{
                  scaling_prob <- (1 - zero_probability_parameter_for_zero_inflated_poisson_1)
                  lambda <- exp(",
                    ff_g,"
                  ) 
                  expect_param <- lambda * exp(o)
                  expect <- scaling_prob * expect_param
                  variance <- scaling_prob * expect_param *
                    (1 + (1 - scaling_prob) * expect_param)
                  list(
                    lambda = lambda,
                    expect = expect,
                    variance = variance,
                    obs_prob = (1 - scaling_prob) * (y == 0) +
                      scaling_prob * dpois(y, expect_param),
                    prob_zero = (1 - scaling_prob) * (y == 0) +
                      scaling_prob * dpois(0, expect_param)
                  )
                  }")),
              num.threads =32,
              n.samples = nsamp)# %>% exp()
    # } else{
    #   preds <- inlabru::generate(res, newdata=d_test,
    #                              ~{
    #                                scaling_prob <- (1 - zero_probability_parameter_for_zero_inflated_poisson_1)
    #                                lambda <- exp(
    #                                  RecLength + t2se + time_group + kappa + doy + alpha +
    #                                   Intercept
    #                                ) 
    #                                expect_param <- lambda * exp(o)
    #                                expect <- scaling_prob * expect_param
    #                                variance <- scaling_prob * expect_param *
    #                                  (1 + (1 - scaling_prob) * expect_param)
    #                                list(
    #                                  lambda = lambda,
    #                                  expect = expect,
    #                                  variance = variance,
    #                                  obs_prob = (1 - scaling_prob) * (y == 0) +
    #                                    scaling_prob * dpois(y, expect_param),
    #                                  prob_zero = (1 - scaling_prob) * (y == 0) +
    #                                    scaling_prob * dpois(0, expect_param)
    #                                )
    #                              },
    #                              num.threads =32,
    #                              n.samples = nsamp)# %>% exp()
    # }
    
    pt <- preds |> transpose()
    preds_out <- map(pt, ~do.call(cbind, .x))
    test_log_scores <- -log(rowMeans(preds_out$obs_prob))
  }
  
    
    
 
    mean_est <- rowMeans(preds_out$expect)
  
  
  
  # values <- preds_out$expect |>   as_tibble() |>
  #   rowwise() |> 
  #   transmute(median = median(c_across(everything())),
  #             mean = mean(c_across(everything())),
  #             sd = sd(c_across(everything())),
  #             pred_var = mean + sd**2,
  #             var = var(c_across(everything())),
  #             lci = quantile(c_across(everything()),probs=c(.025)),
  #             uci = quantile(c_across(everything()),probs=c(.975))
  #            
  #             
  #             
  #             ) |> 
  #   ungroup() |> 
  #   mutate( log_score = test_log_scores)
  
  values <- tibble(
    median = matrixStats::rowMedians(preds_out$expect),
    mean =  matrixStats::rowMeans2(preds_out$expect),
    sd =  matrixStats::rowSds(preds_out$expect),
    var =  matrixStats::rowVars(preds_out$expect)) |> 
    mutate(pred_var = mean + sd **2,
           log_score = test_log_scores)
    values[,c("lci", "uci")] <- 
      matrixStats::rowQuantiles(preds_out$expect,
                                probs=c(.025, 0.975))
  
  
  # ecdf(exp(preds[1,]))(seq(0,10, by = 0.1)) |> plot()
  suppressMessages({
  pz_long <- 
    (preds_out$prob_zero) |> as_tibble( .name_repair = 'universal', rownames = "rn") |> 
    pivot_longer(, cols = -rn,
                 names_to = "iter", values_to = "p_zero") |> 
    mutate(iter = str_remove(iter,"...")) |> 
    left_join(
      d_test |> st_drop_geometry() |> 
        select(rn, y, test_group), by = join_by(rn)
    ) |> 
    mutate(  dpois_obs = 1-p_zero,
             obs_yn = factor(ifelse(y==0,"Not Observed", "Observed"),
                             levels = c("Observed", "Not Observed")))
  
  
  })
  
  
  if(sum(d_test$y) ==0) {any_observations <- FALSE}else{
    any_observations <- TRUE
  }
 
  if(any_observations){
  med_roc_curve <- 
    pz_long |> ungroup() |> 
    summarise(dpois_obs=median(dpois_obs), obs_yn=unique(obs_yn),
              .by = c(rn) ) |> 
    yardstick::roc_curve( obs_yn, dpois_obs)
  
  test_roc_curve <- 
    pz_long |> ungroup() |> 
    nest_by(iter) |> 
    ungroup() |> 
    slice_sample(n=20) |> 
    # unnest(cols = c(data)) |> 
    # group_by(iter) %>%
    rowwise() |> 
    mutate(roc = list(yardstick::roc_curve(data = data, obs_yn, dpois_obs) )  ) |>
    dplyr::select(-data) |> 
    unnest(cols = c(roc)) |> 
    ggplot(aes(x = 1 - specificity, y = sensitivity, group=iter)) +
    geom_path(alpha=0.1) +
    geom_path(data = med_roc_curve, group=1, colour = 'red')+
    geom_abline(lty = 3) +
    coord_equal() +
    theme_light() +
    labs(x = "False positive rate", y = "True positive rate",
         title = glue::glue("{spp} - {test_group_} "))} else{test_roc_curve <- NULL}
  

  test_preds <- bind_cols(d_test, values) |> 
    mutate(
      raw_E =y - median,
      AE = abs(raw_E),
      SE = (raw_E)^2,
      DS = (raw_E)^2 / var + log(var),
      PIT = 1-ppois(y, median) * c(NA_real_, 1)[1 + (d_test$y > 0)]
      # LG = log_score
    )
  
  p2 <- ggplot(test_preds) +
    stat_ecdf(aes(PIT), na.rm = TRUE) +
    scale_x_continuous(expand = c(0, 0)) +
    ggtitle(glue::glue("PIT: {spp} - {test_group_} "))
  
  raw_median_error <- 
  ggplot(test_preds, aes(median, raw_E) ) +
    geom_point() +
    ggtitle(glue::glue("Raw error: {spp} - {test_group_} ")) +
    labs(x = "Median estimate", y = "Median Residual")
  
 
  counts_vs_pred_density <- 
  ggplot() +
    geom_density( data=preds_out$expect[,sample(1:nsamp, 10,F)] |> 
                    as_tibble(rownames = "id") |> 
                    pivot_longer(cols = -id) %>% 
                    mutate(est_counts = rpois(nrow(.), (value))),
                  aes(est_counts,
                      group = name)) +
    geom_density(data = test_preds, aes(y), colour = 'red') +
    scale_x_continuous(trans = scales::log1p_trans()) +
    labs(title = g("counts vs estimates: {spp} - {test_group_}"))
  
    
  scores <- test_preds  |> 
    st_drop_geometry() |> 
    summarise(
      MAD = mean(AE),
      MAD_SD = sd(AE),
      MSE = sqrt(mean(SE)),
      MDS = mean(DS),
    ) 
  library(yardstick)
  test_metrics <- yardstick::metric_set(rmse, rsq, mae, poisson_log_loss)
  
  
  test_lm <- lm(median ~ y, data=test_preds, na.action = "na.exclude")
  
  output_median <- 
  test_preds |> 
    st_drop_geometry() |> 
    summarize(  accuracy = mean(AE)/mean(y) ,
                precision = sd(y)/sd(median),
                cor_spearman = cor(median, y, method="spearman"),
                cor_pearson = cor(median, y,, method="pearson"),
                deviance = 2* sum(ifelse( y == 0, 0, ( y * log( y/median))) -
                                    ( y - median)),
                mean_residual = mean(abs(median)),
                discrim_intercept = test_lm$coefficients[1],
                discrim_slope = test_lm$coefficients[2],
                
                
    ) |> bind_cols(
      bind_rows(test_metrics(test_preds, truth = y, estimate = median),
                roc_auc(  pz_long |> ungroup() |> 
                            summarise(dpois_obs=median(dpois_obs), obs_yn=unique(obs_yn),
                                      .by = c(rn) ), obs_yn, dpois_obs)) |>
        dplyr::select(-.estimator) |> 
        pivot_wider(names_from = .metric, values_from = .estimate)) |> 
    bind_cols(scores) |> 
    mutate(species = spp, test =test_group_) 
  suppressMessages({
  test_full_expect <-  preds_out$expect |> 
    as_tibble(, .name_repair = 'universal', rownames = "rn") |> 
    pivot_longer(, cols = -rn,
                 names_to = "iter", values_to = "expect") |> 
    mutate(iter = str_remove(iter,"...")) |> 
    left_join(
      d_test |> st_drop_geometry() |> 
        select(rn, y, test_group), by = join_by(rn)
    ) |> mutate(AE = abs(expect-y)) 
  })
  
  # test_glm <- lme4::lmer(expect~(y|iter), data = test_full_expect)
  # test_glm <- rlang::try_fetch({glmmTMB::glmmTMB(expect~(y|iter), data = test_full_expect)},
  #                              error = function(x) return(NULL))
  #   
  # if(!is.null(test_glm)){
  # 
  #   test_lmer_ref <- glmmTMB::ranef(test_glm) |> as_tibble() |>
  #   mutate(iter = as.numeric(as.character(grp))) |> arrange(iter) |>
  #   select(iter, term, condval) |> 
  #   pivot_wider(names_from = term, 
  #               values_from = condval, id_cols = iter)
  #   
  # } else{
  #   test_lmer_ref <- NULL
  # }
    
  
  metrics_iter <- 
  test_full_expect |> 
    nest_by(iter) |> 
    mutate(tm= list(test_metrics(data, truth = y, estimate = expect))) |> 
    dplyr::select(-data) |> unnest(tm) |> ungroup()  
  
  auc_iter <-  pz_long |> nest_by(iter) |> 
    mutate(auc = list(roc_auc(data, obs_yn, dpois_obs))) |> 
    dplyr::select(-data) |> unnest(auc) |> ungroup() 
    
  
  test_metrics_iter <- 
  test_full_expect |> 
    summarize(  accuracy = mean(AE)/mean(y) ,
                precision = sd(y)/sd(expect),
                cor_spearman = cor(expect, y, method="spearman"),
                cor_pearson = cor(expect, y,, method="pearson"),
                deviance = 2* sum(ifelse( y == 0, 0, ( y * log( y/expect))) -
                                    ( y - expect)),
                mean_residual = mean(abs(expect)), .by = iter) |> 
    # rowwise() |> 
    # mutate(
    #             discrim_intercept = test_lmer_ref$`(Intercept)`[test_lmer_ref$iter==iter],
    #             discrim_slope = test_lmer_ref$y[test_lmer_ref$iter==iter],
    #            
    #             
    #             
    # ) |> 
    left_join(
      bind_rows(metrics_iter, auc_iter) |> 
        pivot_wider(names_from = .metric, values_from = .estimate, id_cols = iter)
      , by = "iter") |>  
    mutate(species = spp, test =test_group_) 
  
  tictoc::toc()
  
  list(plot_roc_curve = test_roc_curve,
       plot_PIT = p2,
       median_metrics =output_median,
       iteration_metrics =test_metrics_iter,
       plot_counts_vs_pred_density=counts_vs_pred_density,
       plot_raw_median_error =raw_median_error,
       family = fam,
       test_group = test_group_
       )
  
  
}
  xxx <- c("Train","Spatial", "Site")
  # mirai::daemons(3)
  all_tests <- map(xxx,gen_inla_test)#,  .parallel = TRUE )
  names(all_tests) <- xxx
  
  t_test <- transpose(all_tests)
  
  median_metrics_all <- bind_rows(t_test$median_metrics)
  iter_metrics_all <- bind_rows(t_test$iteration_metrics)
  
  cairo_pdf(g("{out_dir}/roc_{spp}.pdf"))
  print(t_test$plot_roc_curve)
  dev.off()
  
  cairo_pdf(g("{out_dir}/PIT_{spp}.pdf"))
  print(t_test$plot_PIT)
  dev.off()
  
  cairo_pdf(g("{out_dir}/counts_density_{spp}.pdf"))
  print(t_test$plot_counts_vs_pred_density)
  dev.off()
  
  cairo_pdf(g("{out_dir}/median_error_{spp}.pdf"))
  print(t_test$plot_raw_median_error)
  dev.off()
  
  
  
   # write_rds(t_test, g("{out_dir_tmp}/PIT_{spp}.rds"))
   
   write_csv(median_metrics_all,g("{out_dir}/median_metrics_{spp}.csv") )
   write_rds(median_metrics_all,g("{out_dir_tmp}/median_metrics_{spp}.rds") )
   write_csv(iter_metrics_all,g("{out_dir}/iter_metrics_{spp}.csv") )
   write_rds(iter_metrics_all,g("{out_dir_tmp}/iter_metrics_{spp}.rds") )

  tictoc::toc()
}

# mirai::daemons(10)
# 
# crated_f <- carrier::crate(~inla_tests(.x))

spp_modeled <- 
list.files(INLA_output_loc_TMP, "_inlabru_model.rds", recursive = T) |> 
 stringr::str_extract("\\w{4}") |> unique() ## |> 
  # str_subset("ATSP", negate = T) |> 
  
spp_comp <- 
  list.files(INLA_output_loc, "iter_metrics_", recursive = T) |> 
  stringr::str_extract("\\w{4}") |> unique()


m <- walk(spp_modeled[!spp_modeled %in% spp_comp], inla_tests)


# library(tidyverse)
# 
# all_spp_mat <- 
# list.files(INLA_output_loc, pattern = "median_metrics", recursive = T, full.names = T) |> 
#   read_csv()
# 
# ggplot(all_spp_mat,
#        aes( factor(test, levels = c("Train", "Site", "Spatial")),roc_auc,
#            group = species)) +
#   geom_line() +
#   ggrepel::geom_label_repel(data = filter(all_spp_mat, test == "Spatial"),
#                             aes(label = species))
#   scale_colour_viridis_d()




