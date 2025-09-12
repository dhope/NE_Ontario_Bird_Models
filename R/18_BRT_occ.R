df_std <- setup_dat |>#st_drop_geometry() |> 
  dplyr::select(Row, y,
                t2se_scaled,
                doy_r,RL,
                rec_id,
                SiteN,location,
                event_gr,
                offset,
                event_id,
                X_sc,
                Y_sc,X,Y,
                type,
                site_id_agg = site_id_agg.x,
                # sum_count_no_tmtt,
                geometry,
                # es_f,
                # all_of(inc_vi)
  ) |> 
  left_join(pca_out |> mutate(across(location, as.character)), by = join_by(location)) |> 
  left_join(pls_out  |> mutate(across(location, as.character)), by = join_by(location)) %>% {
    if(length(included_factors)>0){
      left_join(., distinct(spatial_cov) |> st_drop_geometry() |> 
                  select(location,{{included_factors}}),
                by = join_by(location))
    } else{
      .
    }
  } |> 
  dplyr::select(-location, -mean_count) |> 
  filter(!is.na(offset)) |>
  dplyr::select(where(~sum(is.na(.x))!=length(.x) ))  |> 
  # dplyr::select(-c(project_id, NA_codes, location_id,sum_count, 
  #                  species_code, species_common_name)) |> 
  mutate(event_id = factor(event_id)) |>
  mutate(across(where(is.numeric) & !c(X_sc, Y_sc, offset, y,event_gr,SiteN,
                                       geometry,X,Y,Row,site_id_agg,
                                       t2se_scaled, doy_r), arm::rescale),
         .keep = 'all') %>%  
  # rename(y = sum_count_no_tmtt) %>%
  mutate(random_val = rnorm(n=nrow(.)),
         pres = as.numeric(y>0),
         site = (factor(site_id_agg))) %>%  {
           if(run_a2){mutate(.,t2se_sc = t2se_scaled)} else{
             bind_cols(., model.matrix(~RL + t2se_scaled*event_gr, data = .) |> as_tibble() %>%
                         setNames(nm = c("Intercept","RL3", 
                                         "RL5", "RL10", "t2se_sc", "Dusk", 
                                         "t2se_dusk") ) )
           }}




det_var <- st_drop_geometry(df_std) |> 
  arrange(site) |> 
  # bind_cols(B_dawn |> as_tibble() |> rename_with(everything(),
  #                                           .fn = ~{paste0("t2se_Dawn_B_", .x)})) |> 
  # bind_cols(B_dusk |> as_tibble() |> rename_with(everything(),
  #                                           .fn = ~{paste0("t2se_Dusk_B_", .x)})) |> 
  dplyr::select(site,pres,  doy_r,RL,type, event_gr, 
                t2se_scaled,
                RL) |> 
                # RL3:t2se_dusk) |> 
  # t2se_sc,Dusk,t2se_dusk,
  # # starts_with("t2se_")
  #          ) |> 
  mutate(obs_n = row_number(), .by = c(site))  |> 
  arrange(site, obs_n)
det_var_names <- c('pres', 'doy_r', 'RL', 'type', "event_gr", "t2se_scaled")
det_vars <- map(det_var_names,
                ~{dplyr::select(det_var, site, obs_n, {{.x}}) |> 
                    pivot_wider(id_cols = site, values_from = {{.x}},
                                names_from = obs_n) |> 
                    arrange(site) |> 
                    column_to_rownames("site") |> 
                    as.matrix()
                } ) |> setNames(nm = occ_var_names)

occ_vars <- 
    distinct(
      select(st_drop_geometry(df_std),
             site,PLS01:PLS10, X,Y,
      )) |> arrange(site) |> column_to_rownames('site') |> 
  as.matrix()



XY_coords <-occ_vars[,c("X", "Y")] 


data.list <- list(y = det_vars$pres, 
                  occ.covs = occ_vars[,!colnames(occ_vars) %in% c("X", "Y")], 
                  det.covs = purrr::discard(det_vars, str_detect(names(det_vars),"pres")), 
                  coords = XY_coords)

# Number of batches
n.batch <- 400
# Batch length
batch.length <- 25
n.iter <- n.batch * batch.length
# Priors 
prior.list <- list(beta.normal = list(mean = 0, var = 2.72), 
                   alpha.normal = list(mean = 0, var = 2.72),
                   sigma.sq.ig = c(2, 2), 
                   phi.unif = c(3/1, 3/.1)) 
# Initial values
inits.list <- list(alpha = 0, beta = 0,
                   phi = 3 / .5, 
                   sigma.sq = 2,
                   w = rep(0, nrow(X)),
                   z = apply(y, 1, max, na.rm = TRUE))
# Tuning
tuning.list <- list(phi = 1) 

# n.samples <- 5000
n.burn <- 2000
n.thin <- 20
n.chains <- 3
n.report <- 1000

# Note that this is just a test case and more iterations/chains may need to 
# be run to ensure convergence.
library(spOccupancy)
out <- spPGOcc(det.formula = ~ doy_r + RL + type + t2se_scaled*event_gr, 
               occ.formula = ~ PLS01 + PLS02 + PLS03 + PLS04 + PLS05, 
               data = data.list, 
               # inits = inits.list, 
               n.batch = n.batch, 
               batch.length = batch.length, 
               priors = prior.list,
               cov.model = "exponential", 
               tuning = tuning.list, 
               NNGP = FALSE, 
               n.neighbors = 5, 
               search.type = 'cb', 
               n.report = n.report, 
               n.burn = n.burn, 
               n.chains = n.chains, 
               n.omp.threads = 10)

