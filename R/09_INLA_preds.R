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



prep_predictions <- function(spp, save_objects){
  if(run_a2){
    out_dir_app <- "A2"
  } else{
    out_dir_app <- ""
  }
  out_dir_spatial <- g(
    "{INLA_output_loc_spatial}/{out_dir_app}/{spp}"
  )
  out_dir <- g(
    "{INLA_output_loc}/{out_dir_app}/{spp}"
  )
  out_dir_tmp <- g(
    "{INLA_output_loc_TMP}/{out_dir_app}/{spp}"
  )
  
  res <- read_rds(  g("{out_dir_tmp}/{spp}_inlabru_model.rds"))  
  fam <- res$.args$family
  pca_cov <- read_rds(g("{out_dir_tmp}/pca_cov.rds"))
  pca_prep <- read_rds(g("{out_dir_tmp}/pca_bundle.rds")) |> bundle::unbundle() |> 
    recipes::prep()
print(g("Preping predictions, {spp} --------------------------"))

pred_date <- ifelse(run_a2, "2025-01-10", "2025-02-28")

preds_r <- read_rds(g("{prediction_layer_loc}/{pred_date}_prediction_rasters_df.rds")) |> #2025-02-11_prediction_rasters_df.rds") |> 
  mutate(rn = row_number()) 
xy <- read_rds(g("{prediction_layer_loc}/{pred_date}_prediction_rasters_xy.rds"))

# names(df_std)[names(df_std) %in% names(preds)]
pca_pred_data <- preds_r |> 
  mutate(location = glue::glue("loc-{1:nrow(preds_r)}")) |>
  dplyr::select( any_of(names(pca_cov)) )
pca_preds <- recipes::bake(pca_prep, new_data = pca_pred_data) |> 
  dplyr::select(-location)
pc_vars <- str_subset(names(pca_preds), "^PC\\d")

list2env(read_rds(g("{out_dir_tmp}/sd_means.rds")), environment())
if(!run_a2) t2se_mod <- read_rds(g("{out_dir_tmp}/t2se_mod.rds"))
offsets_spp <- read_rds(g("{out_dir_tmp}/offsets_spp.rds"))




scaled_prediction <- 
  preds_r |> 
  bind_cols(pca_preds) %>% 
  dplyr::select(rn,which(names(. )%in% pc_vars)) |> 
  pivot_longer(cols = -rn,
               names_to = 'variable', values_to = "x") |> 
  left_join(sd_mn,by = join_by(variable)) |> 
  mutate(xx = (x-mn)/sd2) |> 
  dplyr::select(rn, variable, xx) |> 
  pivot_wider(names_from = variable, values_from = xx, id_cols = rn) |> 
  arrange(rn)



preds_sc <- preds_r |>  mutate(X_sc = X/10000,
                               Y_sc = Y/10000) %>% 
  # bind_cols(pca_preds) %>% 
  dplyr::select(X,Y,!which(names(. )%in% pc_vars)) |> #sd_mn$variable)) |> 
  bind_cols(scaled_prediction |> dplyr::select(-rn)
  ) |> mutate(recording_id = factor("None")) |> 
  st_as_sf(coords = c("X", "Y"), crs = ont.proj)

preds_sc$RL <- factor(5, levels = RL_Fac )
# preds_sc$t2se_scaled <- t2se_mod
if(!run_a2) {time_period_prediction <- slice_max(t2se_mod, order_by = median, with_ties = F) |> pull(event_gr)
} else{time_period_prediction <- "Dawn"}
period_to_use <- 
  preds_sc$t2se_sc <- (0-t2se_sd_mn$mean_t2se[t2se_sd_mn$event_gr==time_period_prediction&!is.na(t2se_sd_mn$event_gr)])/
  (2*t2se_sd_mn$sd_t2se[t2se_sd_mn$event_gr==time_period_prediction&!is.na(t2se_sd_mn$event_gr)]) 
# preds_sc$doy_r <- doy_mod
preds_sc$doy_r <- (161-sd_mn$mn[sd_mn$variable=="doy"])/
  (sd_mn$sd2[sd_mn$variable=="doy"])
preds_sc$event_gr <- time_period_prediction
preds_sc$X <- preds_r$X
preds_sc$Y <- preds_r$Y
preds_sc$site <- 437:(436+nrow(preds_sc))
preds_sc$offset <- offsets_spp |> distinct(dur, max_dist, o) |> 
  filter(dur==5) |> pull(o)
# pres_filt <- filter(preds_sc,!is.na(fnlc_9_100))
tictoc::tic()
nsamp <- 200 # number of posterior samples

# samples from mgcv model 
# Generate estimates from SPDE estimate only
# kappa_var <- 1/res$summary.hyperpar$mean[3]

n_splits <- 10

pred_test <- 
  preds_sc |> 
  filter( !is.na(clc20_1_100)) |> 
  dplyr::select(rn,all_of(pc_vars), X,Y, t2se_sc, event_gr,RL,offset,
                doy_r) %>%
  mutate(#chunk = cut_interval(rn, n=n_splits),
         rn_test = row_number())


ff <- paste0("~ Intercept + ",
             paste0("spat_cov_",pc_vars, collapse = "+"),
             "+ alpha") |> 
  as.formula()

ff_sp <- ~ Intercept +  alpha
ff_g <- glue::glue('{res$misc$configs$contents$tag[-c(1,2)] |> glue::glue_collapse(sep = " + ")}')

if(run_a2 | str_detect(ff_g,"PC01", negate = T)){
  ff_counts <- paste0("Intercept +  o + t2se + doy +",
                      paste0("spat_cov_",pc_vars, collapse = "+"),
                      "+ alpha") #|>
} else{
ff_counts <- paste0("Intercept + RecLength + o + t2se + time_group + doy +",
                    paste0("spat_cov_",pc_vars, collapse = "+"),
                    "+ alpha") #|>
}



ff_complex <-  paste0("~{",
                      "expect <- exp(",ff_counts, ")
                       expect_sp <- exp(Intercept +  alpha)
                      list(expect = expect,",
                      "pobs = 1-dpois(0, expect),
                      expect_sp = expect_sp)}") |> 
  as.formula()


ff_sp <-  paste0("~{
                       expect_sp <- exp(Intercept +  alpha)
                      list(expect_sp = expect_sp,
                      pobs_sp = 1-dpois(0, expect_sp) )}") |> 
  as.formula()

ff_zi <- paste0("~{
                 scaling_prob <- (1 - zero_probability_parameter_for_zero_inflated_poisson_1)
                lambda <- exp(", ff_counts,  ")
                lambda_sp <- exp(Intercept +  alpha)
                expect <- scaling_prob * lambda
                expect_sp <- scaling_prob * lambda_sp
                list(
                expect = expect,
                pobs = 1 - dpois(0, expect),
                expect_sp = expect_sp)
                }") |> as.formula()
                
ff_zi_sp <- paste0("~{
                 scaling_prob <- (1 - zero_probability_parameter_for_zero_inflated_poisson_1)
               lambda_sp <- exp(Intercept +  alpha)
                expect_sp <- scaling_prob * lambda_sp
                list(
                pobs_sp = 1 - dpois(0, expect_sp),
                expect_sp = expect_sp)
                }") |> as.formula()                







if(isTRUE(save_objects)){
  write_rds(list(
    pred_test = pred_test, ff_complex = ff_complex,
    ff_zi = ff_zi, ff_zi_sp =ff_zi_sp, fam = fam,ff_g = ff_g,
    ff_sp = ff_sp 
  ), g("{out_dir_tmp}/temp_ob_for_pred_{spp}.rds"))
} else{
  as.list(environment())
}
}

run_predictions <- function(spp){
  if(run_a2){
    out_dir_app <- "A2"
  } else{
    out_dir_app <- ""
  }
  out_dir_tmp <- g(
    "{INLA_output_loc_TMP}/{out_dir_app}/{spp}"
  )
  out_dir_tmp_preds <- g(
    "{INLA_preds_loc_TMP}/{out_dir_app}/{spp}"
  )
  out_rds_file <- g("{out_dir_tmp_preds}/{spp}_predictions_full.rds")
  if(file.exists(out_rds_file)) return(NULL)
  dir.create(out_dir_tmp_preds, recursive = T)
  list2env(prep_predictions(spp, save_objects = FALSE), environment())
nsamp <- 200
print(g("Running predictions, {spp} --------------------------"))
if(str_detect(ff_g, "PC01")){
  tic()
preds <- inlabru::generate(res, pred_test,
                  switch (fam,
                    'poisson' = ff_complex,
                    "zeroinflatedpoisson1" = ff_zi
                  ),num.threads =1,
                  n.samples = nsamp)# %>% exp() 

toc()

} else{
gc()
tic()
preds <- inlabru::generate(res, pred_test,
                           switch (fam,
                                   'poisson' = ff_sp,
                                   "zeroinflatedpoisson1" = ff_zi_sp
                           ),num.threads =1,
                           n.samples = nsamp)# %>% exp() 
toc()
}
write_rds(preds, out_rds_file)


}

gen_maps <- function(spp){
  if(run_a2){
    out_dir_app <- "A2"
  } else{
    out_dir_app <- ""
  }
  out_dir_tmp <- g(
    "{INLA_output_loc_TMP}/{out_dir_app}/{spp}"
  )
  out_dir_tmp_preds <- g(
    "{INLA_preds_loc_TMP}/{out_dir_app}/{spp}"
  )

  out_dir_spatial <- g(
    "{INLA_output_loc_spatial}/{out_dir_app}/{spp}"
  )
  out_dir <- g(
    "{INLA_output_loc}/{out_dir_app}/{spp}"
  )
tictoc::toc()

preds <- read_rds( g("{out_dir_tmp_preds}/{spp}_predictions_full.rds"))


pres_counts <-preds |> # map(1:n_splits,
                   # ~{pred_split[[.x]] |> 
  # list_flatten(pred_split) |>
  transpose() |>
  pluck(1) %>% do.call("cbind", .)
                   # }) %>% do.call("rbind", .)


cat(glue::glue("{spp}---------------------\n"))
pres_obs <-
preds |> 
                    transpose() |>
                    pluck(2) %>% do.call("cbind", .)
tictoc::toc()


pred_date <- ifelse(run_a2, "2025-01-10", "2025-02-28")

xy <- read_rds(g("{prediction_layer_loc}/{pred_date}_prediction_rasters_xy.rds"))

preds_r <- read_rds(g("{prediction_layer_loc}/{pred_date}_prediction_rasters_df.rds")) |> #2025-02-11_prediction_rasters_df.rds") |> 
  mutate(rn = row_number()) 

pred_test <- 
  preds_r |> 
  filter( !is.na(clc20_1_100))


pred_test$mean_count <- apply((pres_counts), 1, mean)
pred_test$p_obs <- apply((pres_obs), 1, mean)
# pred_test$mean_estimate[pred_test$mean_estimate>quantile(pred_test$mean_estimate, probs = 0.99)] <- NA
pred_test$median_count <- apply(pres_counts, 1, median)
# pred_test$median_sp <- apply(pred_sp, 1, median)
pred_test$sd_count <- apply(pres_counts, 1, sd)
# preds_sc$spatial_mean[!is.na(preds_sc$fnlc_10_500)] <- apply(modpred_xy, 1, mean)
# preds_sc$spatial_sd[!is.na(preds_sc$fnlc_10_500)] <- apply(modpred_xy, 1, sd)
# preds_sc$sd[!is.na(preds_sc$fnlc_10_500)] <-matrixStats::rowSds(x = M_predictions) #apply(modpred, 1, sd)
pred_test[, c("lci_count", "uci_count")] <-matrixStats::rowQuantiles(pres_counts, probs=c(.025, .975)) #apply(modpred, 1, sd)
pred_test[, c("lci_p", "uci_p")] <-matrixStats::rowQuantiles(pres_obs, probs=c(.025, .975)) #apply(modpred, 1, sd)
# pred_test[, c("lci_sp", "uci_sp")] <-matrixStats::rowQuantiles(pred_sp, probs=c(.025, .975)) #apply(modpred, 1, sd)
pred_test$ci_size_count <- pred_test$uci_count - pred_test$lci_count


# preds_sc$mean_estimate[!is.na(preds_sc$fnlc_10_500)] <-matrixStats::rowMeans2(M_predictions)# apply(modpred, 1, mean)
# preds_sc$simulation[!is.na(preds_sc$fnlc_10_500)] <- apply(simulated_density, 1, median)
# p_avg$sd <- p_sd


# full_p <- left_join(preds_sc, p_avg) 

preds_sc <- left_join(preds_r |> select(rn),
                      pred_test |> 
                        st_drop_geometry() |> 
                        dplyr::select(rn,
                                      mean_count,
                                      p_obs,
                                      median_count,
                                      # median_sp,
                                      sd_count,
                                      lci_count,
                                      uci_count,
                                      lci_p,
                                      uci_p,
                                      # lci_sp,
                                      # uci_sp,
                                      ci_size_count),
                      by = join_by(rn)
)

preds_sc <- preds_sc |> 
  mutate(
    cv_count = sd_count/mean_count,
    # cv_spatial = spatial_mean / spatial_sd,
    k_cv = sqrt(cv_count^2/(1+cv_count^2)) 
  )
# k_cv_spatial = sqrt(cv_spatial^2/(1+cv_spatial^2)))#,
# est_mod = ifelse(est>=log(0.25), log(0.25), est))


# ggplot(p, aes(X_sc, Y_sc, colour =est )) + geom_point(alpha = 0.4)
r_pred <- rast(glue::glue("{prediction_layer_loc}/{pred_date}_prediction_rasters.nc"))

r_temp <- rast(r_pred[[g("{pred_date}_prediction_rasters_1")]]) 

in_mesh <- fmesher::fm_is_within(st_as_sf(xy, coords = c("x", "y"), crs= 3161), mesh_inla)

r2 <- terra::cellFromXY(r_temp,xy[in_mesh,])
median_count <- #median_sp <- 
  sd_r <- lci_r <- lci_p <- #lci_sp <- uci_sp <- 
  uci_p<- uci_r <- mean_estimate <- 
  p_obs <- cv_r <-  k_cv <-  ci_size_r <- r_temp



sd_r[r2] <- preds_sc$sd_count[in_mesh]
lci_r[r2] <- preds_sc$lci_count[in_mesh]
lci_p[r2] <- preds_sc$lci_p[in_mesh]
# lci_sp[r2] <- preds_sc$lci_sp[in_mesh]
# uci_sp[r2] <- preds_sc$uci_sp[in_mesh]
uci_r[r2] <- preds_sc$uci_count[in_mesh]
uci_p[r2] <- preds_sc$uci_p[in_mesh]
mean_estimate[r2] <- preds_sc$mean_count[in_mesh]
median_count[r2] <- preds_sc$median_count[in_mesh]
# median_sp[r2] <- preds_sc$median_sp[in_mesh]
p_obs[r2] <- preds_sc$p_obs[in_mesh]
ci_size_r[r2] <- preds_sc$ci_size_count[in_mesh]
cv_r[r2] <- preds_sc$cv_count[in_mesh]
# cv_spatial[r2] <- preds_sc$cv_spatial
k_cv[r2] <- preds_sc$k_cv[in_mesh]
# k_cv_spatial[r2] <- preds_sc$k_cv_spatial



outstack <- c(sd_r,
              lci_r,lci_p,
              # uci_sp,lci_sp,
              uci_r,uci_p,
              ci_size_r,
              mean_estimate,
              median_count,
              # median_sp,
              p_obs,
              cv_r,
              # cv_spatial,
              k_cv)
names(outstack) <- varnames(outstack) <- c(
  "sd_count", "lci_count", "lci_p",
  # "uci_sp", "lci_sp",
  "uci_count", "uci_p",
  "ci_size_count",
  "mean_count","median_count",
  # "median_sp",
  "p_obs",
  "cv_count", #"cv_spatial_only", 
  "k_cv")


# terra::writeCDF(outstack, filename = glue::glue("{out_dir_spatial}/Predictions_{spp}.nc"),overwrite=T, split=T)


for (i in  names(outstack)){
  terra::writeRaster(outstack[[i]], glue::glue("{out_dir_spatial}/{i}_{spp}.tif"),overwrite=T)
}


expectations <- str_subset(names(outstack),"mean|median|p_obs")
errors <- str_subset(names(outstack),"mean|median|p_obs", negate = T)
napken_lake <- read_sf(
  nl_loc
) |> st_transform(ont.proj)
locations <- read_sf(
 place_names_loc
) |> st_transform(ont.proj) 

# locs_in <- locations[fmesher::fm_is_within(locations, mesh_inla),]

plot_map <- function(i){
  min_q <- 0
  if(str_detect(i, "_sp")) min_q <-  min( preds_sc[[i]],  na.rm=T)
  print(ggplot() +
          tidyterra::geom_spatraster(data = outstack[[i]]) +
          labs(fill = i, title = spp) +
          tidyterra::scale_fill_whitebox_c(
            palette = "muted",#breaks = c(0.05, 0.1, 0.2, 0.3, 0.8),
            # labels = scales::label_number(suffix = "indiv/ha"),
            n.breaks = 8,limits= c(min_q, quantile( preds_sc[[i]], 0.995, na.rm=T)),
            guide = guide_legend(reverse = TRUE)) +
          geom_sf(data = napken_lake, shape =2 )+
          # geom_sf(data = locs_in) +
          geom_sf(data = ra_area, fill = NA, linetype =2, colour = 'white') +
          geom_sf(data = mesh_data, fill = NA, linetype = 3, colour = 'black'))
  
}

cairo_pdf(g("{out_dir}/maps_expectations_{spp}.pdf") )
map(expectations, plot_map)
dev.off()

cairo_pdf(g("{out_dir}/maps_uncertainty_{spp}.pdf") )
map(errors, plot_map)
dev.off()


}
