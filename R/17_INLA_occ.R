library(tidyverse)
library(sf)
library(fmesher)
library(patchwork)
source("R/__globals.R")
library(INLA)
library(inlabru)
library(tictoc)
rn_vars <- c("X_")
test_training_data <- read_rds(g("output/{date_compiled}_test_training_data.rds"))

train_locs <-  test_training_data$train_recordings |> 
  left_join(aggregated_locs,by = join_by(project, location, collection)) |> 
  st_as_sf() |> 
  dplyr::distinct(site_id_agg, geometry) |> 
  filter(!st_is_empty(geometry))

run_a2 <- FALSE

if(run_a2){
  spp_cov_file <- "output/rds/2025-02-27_spatial_covariates_Atlas2.rds"
  individual_locs <- 
    g("{rds_data_loc}/locations.rds") |> 
    read_rds() |> 
    filter(collection %in% c( "OBBA2PC") &
             str_detect(project, "Nocturnal|Resample", negate=T)) |> 
    st_filter(ra_buffer) |> mutate(site_id_agg = location)
  aggregated_locs <- dplyr::select(individual_locs,site_id_agg)
  
  raw_recordings <- g("{rds_data_loc}/all_events.rds" ) |> 
    read_rds() |> 
    filter(collection %in% c( "OBBA2PC") & location %in% individual_locs$location) |> 
    mutate(date = ymd(glue::glue("{year}-{month}-{day}")),
           doy = yday(date))
  raw_counts <- read_rds(g("{rds_data_loc}/counts.rds")) |>
    replace_na(list(collection="WildTrax")) |> 
    filter(event_id %in% raw_recordings$event_id ) |> 
    filter(str_detect(project, "(Extraction)|(Nocturnal)|(Resample)",negate=T))
  out_dir_app <- "A2"
  
} else{
  spp_cov_file <- "output/rds/2025-02-28_spatial_covariates_data.rds"
  raw_recordings <- test_training_data$train_recordings |> 
    left_join(aggregated_locs,by = join_by(project, location, collection)) |> 
    filter(str_detect(project, "(Extraction)|(Nocturnal)|(Resample)",negate=T))
  
  raw_counts <- read_rds(g("{rds_data_loc}/counts.rds")) |>
    replace_na(list(collection="WildTrax")) |> 
    filter(event_id %in% raw_recordings$event_id ) |> 
    filter(str_detect(project, "(Extraction)|(Nocturnal)|(Resample)",negate=T))
  out_dir_app <- ""
}

spatial_cov <- 
  spp_cov_file |> 
  read_rds() |>bind_cols(individual_locs) |>  #Spatial_covariates_data_14March2024.rds") |> 
  distinct() |> st_as_sf() %>%
  bind_cols(as_tibble(st_coordinates(.))) |>
  dplyr::select(where(~sum(is.na(.x))!=length(.x))) |> 
  left_join(st_drop_geometry(aggregated_locs ),
            by = join_by(site_id_agg))





recordings <- raw_recordings |>
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
    recording_id = as.numeric(recording_id)) |> 
  mutate(Rec_length = factor(as.numeric(round(dur,1)))) |> 
  filter(Rec_length %in% c(1, 3, 5, 10))


rm(test_training_data, train_locs, raw_recordings)

counts <- raw_counts |> 
  dplyr::select(event_id, location,  project, species_name_clean,total_count, 
                species_code,
                total_count_with_tmtt) |> 
  mutate(y = ifelse(is.na(total_count_with_tmtt), total_count,
                    total_count_with_tmtt)) |> 
  filter(!is.na(y))

rm(raw_counts)

qpad_offsets <- read_rds("output/QPAD_global_offsets.rds") |> 
  rename(max_dist =r, time_minutes = t)
na_pops_offsets <- read_rds("output/na_pops_offsets.rds") |> 
  rename(species =spp)  


prep_inla_data <- function(spp, run_a2){
  if(!run_a2){
    spp_name <- spp_codes$species_name_clean[spp_codes$species_code==spp ]
    # counts$species_name_clean[counts$species_code==spp & !is.na(counts$species_code)] |> 
    # unique() 
    if(length(spp_name)!=1)rlang::abort("unable to identify species name")
  } else { spp_name <- "Palm Warbler"}
  counts_spp <- filter(counts, species_name_clean == spp_name | species_code == spp) |> 
    full_join(recordings %>%{
      if("geometry" %in% names(.)){
        dplyr::select(.,-geometry)} else{.}},
      by = join_by(event_id, location, project)) |> 
    replace_na(list(total_count_with_tmtt = 0, total_count = 0,
                    y=0))
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
  
  
  
  
  # tictoc::tic()
  
  
  
  
  setup_dat_0 <- counts_spp |> 
    left_join(offsets_spp |> 
                dplyr::select(event_id, offset =o),
              by = join_by( event_id)) |> 
    
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
    filter(event_gr!="Day") |> 
    rowwise() |>
    mutate(t2se_scaled = list(arm::rescale(data$t2se))) |> 
    unnest(c(data, t2se_scaled)) |> 
    mutate(doy_r = arm::rescale(doy)) |> 
    
    left_join(spatial_cov,by = join_by(location,project, collection) ) |>
    st_as_sf() |> 
    filter(!st_is_empty(geometry)) |>
    mutate(X_sc = X/10000,
           Y_sc = Y/10000) |> 
    
    ## !!! REMOVING REC_LENGTH 2.4 as there are only 8 points
    
    # filter(Rec_length != 2.4) |> 
    mutate(RL = droplevels(Rec_length),
           event_gr = droplevels(event_gr))
  
  
  as.list.environment(environment())
  
}


run_inlabru <- function(spp){
  
  print(glue::glue("Running {spp} --------------------------------------"))
  out_dir_spatial <- g(
    "{INLA_output_loc_spatial}/{out_dir_app}/{spp}"
  )
  out_dir <- g(
    "{INLA_output_loc}/{out_dir_app}/{spp}"
  )
  out_dir_tmp <- g(
    "{INLA_output_loc_TMP}/{out_dir_app}/{spp}"
  )
  dir.create(
    out_dir_spatial, recursive = T
  )
  dir.create(
    out_dir, recursive = T
  )
  dir.create(
    out_dir_tmp, recursive = T
  )
  
  list2env( prep_inla_data(spp, run_a2 = run_a2), envir = environment())
  
  t2se_sd_mn <- setup_dat_nested |>
    rowwise() |> 
    mutate(mean_t2se = mean(data$t2se),
           sd_t2se = sd(data$t2se),
    ) |> dplyr::select(-data)
  
  rm(setup_dat_nested)
  
  withNA <- 
    summarise(st_drop_geometry(setup_dat), across(everything(), 
                                                  ~sum(is.na(.x)))) %>% 
    pivot_longer(cols = everything()) |> 
    filter(value>0) |> 
    pull(name)
  
  library(tidymodels)
  pca_cov <- (spatial_cov) |> 
    filter(location %in% setup_dat$location) |> 
    dplyr::select(location,where(~{!is.factor(.x) & !is.character(.x)}),
                  -site_id_agg, -X, -Y) |> 
    left_join(
      summarise(counts_spp, mean_count = mean(y), .by = location),
      by = join_by(location)
    ) |> 
    st_drop_geometry()
  # library(usdm)
  # vifSel <- vifstep(data.frame(pca_cov[,-1]), th = 5,
  #                   size = 20000)
  
  factor_cov <- spatial_cov |> dplyr::select(where(~{is.factor(.x)| is.character(.x)}), -location,
                                             -geometry) |> st_drop_geometry() |> 
    names()
  
  brt_vi <- g("{BRT_output_loc}/{spp}/{spp}_time_vi.csv")
  if(file.exists(brt_vi)){
    top_vi <- read_csv(brt_vi) |> 
      filter(!Variable  %in% c("X" ,"Y", "doy", "t2se","year",
                               "Time_period_Dusk","total_100","event_gr_Dusk", "total_500"),
             str_detect(Variable,  "Rec_length",negate = T)) |> 
      arrange(desc(Importance)) |> 
      mutate(cVI = cumsum(Importance),
             dcVI = cVI- lag(cVI, default = 0))
    ranv <- top_vi$cVI[top_vi=="random_val"]
    
    inc_vi <- top_vi |> filter(cVI<min(0.8, ranv) & #!Variable %in% withNA &
                                 str_remove(Variable, "_X\\d+") %in% names(spatial_cov)) |>
      pull(Variable)
    
    
    included_factors <- inc_vi[inc_vi %in% factor_cov]
  } else {
    included_factors <- NULL
    file.create(g("{out_dir}/Factors_not_used.txt") )
  }
  
  
  
  
  dat_rec <-
    recipe(data = pca_cov, formula = ~ .) %>%
    update_role(location, new_role = "id") %>%
    update_role(mean_count, new_role = "outcome") %>%
    step_nzv(all_predictors(), freq_cut = 95/5, unique_cut = 10) |> 
    step_impute_median(all_numeric_predictors()) |> 
    step_corr(all_numeric_predictors(), threshold = 0.85) |> 
    step_center(all_predictors()) %>%
    step_scale(all_predictors()) 
  
  pca_rec <- dat_rec %>%
    step_rm(mean_count) |> 
    step_pca(all_numeric_predictors(), id = "pca", num_comp = 10)
  
  pls_rec <- dat_rec %>%
    step_pls(all_numeric_predictors(),outcome ='mean_count' , num_comp = 10)
  # embed::step_umap(all_numeric_predictors(), outcome = vars(avg_price_per_room))
  # step_spline_natural(arrival_date_num, deg_free = 10)
  
  pca_prep <- prep(pca_rec)
  pca_loading <- tidy(pca_prep, id="pca")
  pca_variances <- tidy(pca_prep, id = "pca", type = "variance")
  pca_out <- bake(pca_prep, new_data=NULL)
  
  pls_out <- prep(pls_rec) |> bake(new_data = NULL) 
  # pca <- prcomp(cov[,-1], retx=TRUE, center=TRUE, scale.=TRUE)
  # pca_pred <- predict(pca)
  
  pca_bun <- bundle::bundle(pca_prep)
  
  
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
                                         geometry,X,Y,Row,
                                         t2se_scaled, doy_r), arm::rescale),
           .keep = 'all') %>%  
    # rename(y = sum_count_no_tmtt) %>%
    mutate(random_val = rnorm(n=nrow(.)),
           pres = as.numeric(y>0),
           site = (factor(SiteN))) %>%  {
             if(run_a2){mutate(.,t2se_sc = t2se_scaled)} else{
               bind_cols(., model.matrix(~RL + t2se_scaled*event_gr, data = .) |> as_tibble() %>%
                           setNames(nm = c("Intercept","RL3", 
                                           "RL5", "RL10", "t2se_sc", "Dusk", 
                                           "t2se_dusk") ) )
             }}
  
  
  min_event <- df_std$t2se_scaled
  num_knots <- 8
  library(splines)
  knot_list <- quantile(min_event, probs = seq(0,1, length.out = num_knots))
  B <- (bs(min_event,
            knots=knot_list[-c(1,num_knots)],
            degree=3, intercept = FALSE)) # creating the B-splines
  dim(B)
  
  B_dawn <- B
  B_dawn[df_std$event_gr!="Dawn"] <- NA
  B_dusk <- B
  B_dusk[df_std$event_gr=="Dawn"] <- NA
  
  SSOM <- st_drop_geometry(df_std) |> 
    # bind_cols(B_dawn |> as_tibble() |> rename_with(everything(),
    #                                           .fn = ~{paste0("t2se_Dawn_B_", .x)})) |> 
    # bind_cols(B_dusk |> as_tibble() |> rename_with(everything(),
    #                                           .fn = ~{paste0("t2se_Dusk_B_", .x)})) |> 
    dplyr::select(SiteN,site,pres, Intercept, doy_r,#RL,type, #event_gr, 
                  RL3:t2se_dusk) |> 
                  # t2se_sc,Dusk,t2se_dusk,
                  # # starts_with("t2se_")
                  #          ) |> 
    mutate(obs_n = row_number(), .by = c(site, SiteN))  |> 
    arrange(site, obs_n) |> mutate(Intercept = 1) |> 
    pivot_wider(names_from = obs_n,values_from = c(
      pres, Intercept,
      doy_r, RL3, RL5, RL10, t2se_sc, Dusk, t2se_dusk), 
      id_cols = c(site,SiteN ),names_vary = 'slowest', names_sort = T) |> 
    arrange(site) |> 
    left_join(
      distinct(
        select(st_drop_geometry(df_std),
        site, SiteN, PC01:PLS10, X,Y,
      )) ,by = join_by(site, SiteN)
    )
  
  XY_coords <- SSOM %>% dplyr::select(c(X,Y)) %>% as.matrix()
  
  
  # write_rds(offsets_spp, g("{out_dir_tmp}/offsets_spp.rds"))
  # write_rds(pca_bun,g("{out_dir_tmp}/pca_bundle2.rds") )
  # write_rds(bundle::bundle(prep(pls_rec)),g("{out_dir_tmp}/pls_bundle2.rds") )
  # write_rds(pca_cov,g("{out_dir_tmp}/pca_cov2.rds") )
  
  sd_mn <- 
    setup_dat |>st_drop_geometry() |> 
    left_join(pca_out |> 
                mutate(across(location, as.character)),
              by = join_by(location)) |> 
    left_join(pls_out |> 
                mutate(across(location, as.character)),
              by = join_by(location)) |>
    filter(!is.na(offset)) |>
    mutate(X_sc = X/10000,
           Y_sc = Y/10000) |> 
    dplyr::select(where(~sum(is.na(.x))!=length(.x) ))  |> 
    dplyr::select(-any_of(c("project_id", #NA_codes, sum_count, 
                            "species_code") #species_common_name
    )
    ) |> 
    mutate(event_id = factor(event_id)) |> 
    summarize(across(where(is.numeric) & 
                       !c(X_sc, Y_sc, offset, y),
                     .fns = list(sd2 = ~{2*sd(.x, na.rm=T)},
                                 mn = ~mean(.x, na.rm=T)),
                     .names = "{.col}__{.fn}")
    ) |> 
    pivot_longer(cols = everything(),
                 values_to = 'value',
                 names_to = 'name') |> 
    separate(name, into=c("variable", "estimate_type"),
             sep = "__") |> 
    pivot_wider(names_from = estimate_type, values_from = value)
  
  # write_rds(list(sd_mn=sd_mn, 
  #                t2se_sd_mn=t2se_sd_mn,
  #                RL_Fac  = levels(df_std$RL)), g("{out_dir_tmp}/sd_means.rds"))
  
  sf_dat <- st_as_sf(setup_dat ) |> st_as_sfc()
  
  # make a set of distinct study sites for mapping
  site_map <- df_std %>%
    st_as_sf() |> 
    select(site, X, Y) %>%
    distinct() 
  # make a two extension hulls and mesh for spatial model
  
  
  
  
  
  # mesh_inla$n
  # make spde
  spde <- inla.spde2.pcmatern(
    mesh = mesh_inla,
    prior.range = c(50*1000, 0.1),
    prior.sigma = c(1, 0.5)
  )
  
  A_sp <- inla.spde.make.A(mesh = mesh_inla,loc = XY_coords)
  
  
  smoother_mesh <- seq(-3, 3, by = 0.1) 
  smoother_mesh2 <- seq(-3, 3, by = 0.25) 
  mesh1D <- fm_mesh_1d(smoother_mesh, boundary =  "free",degree = 2)
  mesh1D2 <- fm_mesh_1d(smoother_mesh2, boundary = "free", degree=2)
  
  
  the_spde_t2se <- inla.spde2.pcmatern(mesh1D,
                                       prior.range = c(0.05, 0.05),
                                       prior.sigma = c(0.1, 0.1) )
  the_spde_doy <- inla.spde2.pcmatern(mesh1D,
                                      prior.range = c(0.05, 0.05),
                                      prior.sigma = c(0.1, 0.1) )
  
  pc_vars <- c(str_subset(names(pca_out)[-1], "^PC\\d"),str_subset(names(pls_out)[-1], "^PLS\\d"))
  pca_vars <- c(str_subset(names(pca_out)[-1], "^PC\\d"))
  pls_vars <- str_subset(names(pls_out)[-1], "^PLS\\d")
  spde_inc <- vector('list', length = length(pc_vars))
  for(i in 1:length(pc_vars)) {
    r <- range(df_std[[pc_vars[[i]]]])
    assign(paste0("spde_", pc_vars[[i]]),inla.spde2.pcmatern(fm_mesh_1d(seq(r[[1]]-0.5,r[[2]]+0.5 ), boundary =  "free",degree = 2),
                                                             prior.range = c(0.25, 0.01), # 1% of range less than 0.01
                                                             prior.sigma = c(0.25, 0.01) ) )
    
    assign(paste0("A_", pc_vars[[i]]),inla.spde.make.A(mesh = get(paste0("spde_", pc_vars[[i]])), 
                                                       loc = SSOM[[pc_vars[[i]]]] ) )
    
    
  }

  # PC_sp <- inla.spde.make.A(mesh = spde_PC01,loc = SSOM$PC01)
  
  X_det <-  
    SSOM %>% dplyr::select(matches("(_\\d+$)")) |> 
    dplyr::select(!matches("pres_\\d")) |> 
    # dplyr::select(!matches("Intercept_\\d")) |> 
    # as.matrix()
    #starts_with("RL"),starts_with("type"), starts_with("doy_r"),starts_with("t2se_[Dawn|Dusk]")) |> 
    #   rename_with(.fn = ~str_replace(.x,"r_", "r."), everything()) %>% 
    as.matrix()
    # inla.Occupancy_detCov()
  
    
  
  X_occ2 <- SSOM %>% 
    dplyr::select(PC01:PC10) |> 
    mutate(Int_occ =1,
           across(.cols = starts_with("PC"), .fns = ~rep(NA, nrow(SSOM))),
           spatial_field = rep(NA,nrow(SSOM))) 
  
  
  
  data_list = as.list(X_occ2)
  data_list$Y = SSOM %>% dplyr::select(starts_with("pres_")) %>% as.matrix()
  data_list$X = X_det
  
  data_list_pls <- data_list
  names(data_list_pls) <- str_replace(names(data_list_pls), "PC", "PLS")
  
  
  formula_ssom <- as.formula(
    paste0("inla.mdata(Y,X) ~  -1 + Int_occ +",
           
           paste0("f(",pca_vars, ", model = ",
                  paste0("spde_", pca_vars),", A.local = ", 
                  paste0("A_", pca_vars), ")", collapse = "+"),
           " +
    f(spatial_field, model=spde,  A.local = A_sp)") )
  
  
  formula_ssom_pls <- as.formula(
    paste0("inla.mdata(Y,X) ~  -1 + Int_occ +",
           
           paste0("f(",pls_vars, ", model = ",
                  paste0("spde_", pls_vars),", A.local = ", 
                  paste0("A_", pls_vars), ")", collapse = "+"),
           " +
    f(spatial_field, model=spde,  A.local = A_sp)") )
  
  formula_ssom_sp <- inla.mdata(Y,X) ~  -1 + Int_occ +
    f(spatial_field, model=spde,  A.local = A_sp)
  
  
  
  model_ssom <- inla(formula_ssom,    # model formula
                     data= data_list,        # data 
                     family= 'occupancy', # model likelihood
                     
                     # priors
                     control.fixed =  list(prec = 1/2.72, prec.intercept = 1/2.72),
                     
                     # compute WAIC and DIC
                     control.compute = list(dic = TRUE, waic = TRUE, config = TRUE),
                     verbose = FALSE,
                     
                     # choose link functions for:
                     
                     # (i) the state process (control.link)
                     # (ii) the observation process (link.simple)
                     control.family = list(control.link = list(model = "logit"),
                                           link.simple = "logit",
                                           
                                           # priors for hyperparameters
                                           hyper = list(
                                             # Detection intercept prior
                                             beta1 = list(param = c(0,1), initial = -1),
                                             # Covariate 1 effect prior
                                             beta2 = list(param = c(0,1/2.72)),
                                             # Covariate 2 effect prior
                                             beta3 = list(param = c(0,1/2.72)),
                                             beta4 = list(param = c(0,1/2.72)),
                                             beta5 = list(param = c(0,1/2.72)),
                                             beta6 = list(param = c(0,1/2.72)),
                                             beta7 = list(param = c(0,1/2.72))
                                           )
                     ))
  
  
  model_ssom_pls <- inla(formula_ssom_pls,    # model formula
                     data= data_list_pls,        # data 
                     family= 'occupancy', # model likelihood
                     
                     # priors
                     control.fixed =  list(prec = 1/2.72, prec.intercept = 1/2.72),
                     
                     # compute WAIC and DIC
                     control.compute = list(dic = TRUE, waic = TRUE, config = TRUE),
                     verbose = FALSE,
                     
                     # choose link functions for:
                     
                     # (i) the state process (control.link)
                     # (ii) the observation process (link.simple)
                     control.family = list(control.link = list(model = "logit"),
                                           link.simple = "logit",
                                           
                                           # priors for hyperparameters
                                           hyper = list(
                                             # Detection intercept prior
                                             beta1 = list(param = c(0,1), initial = -1),
                                             # Covariate 1 effect prior
                                             beta2 = list(param = c(0,1/2.72)),
                                             # Covariate 2 effect prior
                                             beta3 = list(param = c(0,1/2.72)),
                                             beta4 = list(param = c(0,1/2.72)),
                                             beta5 = list(param = c(0,1/2.72)),
                                             beta6 = list(param = c(0,1/2.72)),
                                             beta7 = list(param = c(0,1/2.72))
                                           )
                     ))
   
  
 model_ssom_sp <- inla(formula_ssom_sp,    # model formula
                     data= purrr::discard(data_list, str_detect(names(data_list), "PC")),        # data 
                     family= 'occupancy', # model likelihood
                     
                     # priors
                     control.fixed =  list(prec = 1/2.72, prec.intercept = 1/2.72),
                     
                     # compute WAIC and DIC
                     control.compute = list(dic = TRUE, waic = TRUE, config = TRUE),
                     verbose = FALSE,
                     
                     # choose link functions for:
                     
                     # (i) the state process (control.link)
                     # (ii) the observation process (link.simple)
                     control.family = list(control.link = list(model = "logit"),
                                           link.simple = "logit",
                                           
                                           # priors for hyperparameters
                                           hyper = list(
                                             # Detection intercept prior
                                             beta1 = list(param = c(0,1), initial = -1),
                                             # Covariate 1 effect prior
                                             beta2 = list(param = c(0,1/2.72)),
                                             # Covariate 2 effect prior
                                             beta3 = list(param = c(0,1/2.72)),
                                             beta4 = list(param = c(0,1/2.72)),
                                             beta5 = list(param = c(0,1/2.72)),
                                             beta6 = list(param = c(0,1/2.72)),
                                             beta7 = list(param = c(0,1/2.72))
                                           )
                     ))
  
  summary(model_ssom)
  
  gen_ic <- function(mod, name){
    tribble(
      ~model, ~waic, ~dic,
      name, mod$waic$waic, mod$dic$dic)
  }
 
  
  ic_tab <- 
  tribble(
    ~mod, ~name,
    model_ssom, "PCA",
    model_ssom_pls, "PLS",
    model_ssom_sp, "SP"
  ) |> 
    pmap(gen_ic) |> 
    bind_rows() |> 
    mutate(deltawaic = waic-min(waic),
           deltadic = dic - min(dic)) |> 
    arrange(deltawaic)
  
  top_mod <- ic_tab$model[[1]]
  
  pca_cov <- read_rds(g("{out_dir_tmp}/pca_cov2.rds"))
  pca_prep <- read_rds(g("{out_dir_tmp}/pca_bundle2.rds")) |> bundle::unbundle() |> 
    recipes::prep()
  pls_prep <- read_rds(g("{out_dir_tmp}/pls_bundle2.rds")) |> bundle::unbundle() |> 
    recipes::prep()
  
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
  
  pls_preds <- recipes::bake(pls_prep, new_data = pca_pred_data) |> 
    dplyr::select(-location)
  
  pc_vars <- str_subset(names(pca_preds), "^PC\\d")
  pls_vars <- str_subset(names(pls_preds), "^PLS\\d")
  
  
  
  A_spat  <-  inla.spde.make.A(mesh= mesh_inla, loc = as.matrix(xy))
  # inla.spde.make.A(mesh = spde_PC01,loc = SSOM$PC01)
  samples_spat  <-  inla.posterior.sample(1000, model_ssom)
  
  
  if(top_mod=="PCA"){
    for(i in 1:length(pc_vars)) {
    assign(paste0("A_pred_", pc_vars[[i]]),inla.spde.make.A(mesh = get(paste0("spde_", pc_vars[[i]])), 
                                                       loc = pca_preds[[pc_vars[[i]]]] ) )
    
    
    }
  }
  
  if(top_mod=="PLS"){
    for(i in 1:length(pls_vars)) {
      assign(paste0("A_pred_", pls_vars[[i]]),inla.spde.make.A(mesh = get(paste0("spde_", pls_vars[[i]])), 
                                                               loc = pls_preds[[pls_vars[[i]]]] ) )
      
      
    }
  }
  
  
  
  
  func_spat <- function(...)
    
  {
    
    eta = paste0("(Int_occ + ",
            paste0("(get('A_pred_", switch (top_mod,
              "PCA" = pc_vars,
              "PLS" = pls_vars
            ), "') %*% get('", switch (top_mod,
                                       "PCA" = pc_vars,
                                       "PLS" = pls_vars
            ), "'))[,1]", collapse = " + ") , 
             " + (A_spat %*% spatial_field)[,1])")
    parse(text=eta) |> eval()
    
  }
  
  func_sp <-    function(...)
    
  {
    
    eta = Int_occ  + (A_spat %*% spatial_field)[,1]
    eta
    
  }
  
  if(top_mod=="SP"){
    eval_samples  <-  inla.posterior.sample.eval(func_sp, samples_spat)
  } else{
  eval_samples  <-  inla.posterior.sample.eval(func_spat, samples_spat, )
  }
  
  occ_probs  <-  inla.link.logit(eval_samples, inverse = T)
  
  pred_df_results <-  xy %>%
    mutate(eta_sd = matrixStats::rowSds(eval_samples),#apply(eval_samples,1,sd),
           eta_mu = matrixStats::rowMeans2(eval_samples),
           psi_mu = matrixStats::rowMeans2(occ_probs),
           quant_range = matrixStats::rowDiffs(matrixStats::rowQuantiles(eval_samples,probs = c(0.025,
                                                                               0.975)))[,1])
                               # function(x){quantile(x,0.975)-quantile(x,0.025)})) 
  
  # pred_df_results$quant_range <- pred_df_results$quant_range[,1]
  
  
  
  r_pred <- rast(glue::glue("{prediction_layer_loc}/{pred_date}_prediction_rasters.nc"))
  
  r_temp <- rast(r_pred[[g("{pred_date}_prediction_rasters_1")]]) 
  
  in_mesh <- fmesher::fm_is_within(st_as_sf(xy, coords = c("x", "y"), crs= 3161), mesh_inla)
  
  r2 <- terra::cellFromXY(r_temp,xy[in_mesh,])
  # median_count <- #median_sp <- 
  #   sd_r <- lci_r <- lci_p <- #lci_sp <- uci_sp <- 
  #   uci_p<- uci_r <- mean_estimate <- 
  #   p_obs <- cv_r <-  k_cv <-  ci_size_r <- r_temp
  names_outstack <- names(pred_df_results)[-c(1,2)]
  for(i in 1:length(names_outstack)){
    j <- names_outstack[[i]]
    # assign(i, r_temp)
    print(j)
    tt <- r_temp
    tt[r2] <- pred_df_results[[j]][in_mesh]
    
    if(i == 1){
      outstack <- tt
    } else{add(outstack) <- tt}
    names(outstack)[[i]] <- j
    varnames(outstack)[[i]] <- j
    # assign(i, tt)
    
    rm(tt)
  }
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  plot(outstack$psi_mu)
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
    smoother_mesh <- seq(-3, 3, by = 0.1) 
  smoother_mesh2 <- seq(-3, 3, by = 0.25) 
  mesh1D <- fm_mesh_1d(smoother_mesh, boundary =  "free",degree = 2)
  mesh1D2 <- fm_mesh_1d(smoother_mesh2, boundary = "free", degree=2)
  
  
  the_spde_t2se <- inla.spde2.pcmatern(mesh1D,
                                       prior.range = c(0.05, 0.05),
                                       prior.sigma = c(0.1, 0.1) )
  the_spde_doy <- inla.spde2.pcmatern(mesh1D,
                                      prior.range = c(0.05, 0.05),
                                      prior.sigma = c(0.1, 0.1) )
  
  pc_vars <- c(str_subset(names(pca_out)[-1], "^PC\\d"),str_subset(names(pls_out)[-1], "^PLS\\d"))
  pca_vars <- c(str_subset(names(pca_out)[-1], "^PC\\d"))
  pls_vars <- str_subset(names(pls_out)[-1], "^PLS\\d")
  spde_inc <- vector('list', length = length(pc_vars))
  for(i in 1:length(pc_vars)) {
    r <- range(df_std[[pc_vars[[i]]]])
    assign(paste0("spde_", pc_vars[[i]]),inla.spde2.pcmatern(fm_mesh_1d(seq(r[[1]]-0.5,r[[2]]+0.5 ), boundary =  "free",degree = 2),
                                                             prior.range = c(0.25, 0.01), # 1% of range less than 0.01
                                                             prior.sigma = c(0.25, 0.01) ) )
    
    
  }
  
  
  # the_spde2 <- the_spde
  # iid prior
  pc_prec <- list(prior = "pcprec", param = c(0.1, 0.1))
  
  # pcv <- pc_vars
  # pc_vars <- pcv[[1]]
  if(n_distinct(df_std$event_gr)==1){
    t2se_cr <- "t2se(t2se_sc, model=the_spde_t2se) +"
  } else{
    t2se_cr <- "t2se(t2se_sc, model=the_spde_t2se, replicate = event_gr, replicate_layer = 1) +
         time_group(event_gr, model = 'factor_contrast') +"
  }
  
  
  # pc_vars <- pcv[[1]]
  if(n_distinct(df_std$RL)==1){
    RL_cr <- ""
  } else{
    RL_cr <- "RecLength(RL, model = 'factor_contrast') +"
  }
  
  comp_str_base <- 
    paste0("~ o(offset, model = 'const') +
    Intercept(1) + ",
           RL_cr,
           t2se_cr, 
           "kappa(site, model = 'iid', constr = TRUE, hyper = list(prec = pc_prec)) +
    doy(doy_r, model = the_spde_doy) +
    alpha(geometry, model = spde) +
         ")
  
  comp_str_pca <- 
    paste0(comp_str_base,
           paste0("spat_cov_",pca_vars,"(",pca_vars, ", model = ", 
                  paste0("spde_", pca_vars),")", collapse = "+"))  
  # 'linear')", collapse = "+")) |> 
  # , paste0("spde_", inc_vi),")", collapse = "+")) |> 
  #    
  
  comp_str_pls <- 
    paste0(comp_str_base,
           paste0("spat_cov_",pls_vars,"(",pls_vars, ", model = ", 
                  paste0("spde_", pls_vars),")", collapse = "+")) 
  
  comp_str_pls_LINEAR <- 
    paste0(comp_str_base,
           paste0("linear_cov_",pls_vars,"(",pls_vars, ", model = 'linear'", 
                  ")", collapse = "+")) 
  comp_str_pca_LINEAR <- 
    paste0(comp_str_base,
           paste0("linear_cov_",pca_vars,"(",pca_vars, ", model = 'linear'", 
                  ")", collapse = "+")) 
  
  
  if(length(included_factors)>0){
    comp_str_pca <- paste0(comp_str_pca," + ",
                           paste0("factor_cov_",included_factors,"(",included_factors, ", model = 'factor_contrast'",")", collapse = "+") 
    )
    
    comp_str_pls <- paste0(comp_str_pls," + ",
                           paste0("factor_cov_",included_factors,"(",included_factors, ", model = 'factor_contrast'",")", collapse = "+") 
    )
    
    comp_str_pls_LINEAR <- paste0(comp_str_pls_LINEAR," + ",
                                  paste0("factor_cov_",included_factors,"(",included_factors, ", model = 'factor_contrast'",")", collapse = "+") 
    )
    
    comp_str_pca_LINEAR <- paste0(comp_str_pca_LINEAR," + ",
                                  paste0("factor_cov_",included_factors,"(",included_factors, ", model = 'factor_contrast'",")", collapse = "+") 
    )
  }
  
  comp_pca <- as.formula(comp_str_pca)
  comp_pls <- as.formula(comp_str_pls)
  comp_pls_linear <- as.formula(comp_str_pls_LINEAR)
  comp_pca_linear <- as.formula(comp_str_pca_LINEAR)
  
  
  comp_simple <-   as.formula(paste0("~ o(offset, model = 'const') +
    Intercept(1) + ",
                                     RL_cr,
                                     t2se_cr,
                                     "doy(doy_r, model = the_spde_doy) +
     kappa(site, model = 'iid', constr = TRUE, hyper = list(prec = pc_prec)) +
    alpha(geometry, model = spde)  "))
  
  # formula, with "." meaning "add all the model components":
  formula <-  y~.
  
  
  time_limit <- 0.5*60*60
  
  
  library(inlabru)
  flush.console()
  # bru_options_set(control.compute = list(cpo = TRUE))
  
  future::plan(future::multisession, workers = 32)
  inla.setOption(num.threads = "32:1")
  # job::job({
  
  inla.setOption(inla.timeout = time_limit)
  run_mod <- function(name, family_, comps_, ...){
    print(glue::glue("{spp} -- {name}", ))
    tictoc::tic()
    rlang::try_fetch({
      b <- bru(
        comps_,
        bru_obs(
          formula,
          family = family_,
          data = df_std
        ),
        options = list(
          control.compute = list(waic = TRUE, cpo = FALSE),
          control.inla = list(
            int.strategy = "eb",
            strategy = "adaptive"
          ),
          verbose = F
        ) )
      tictoc::toc()
      b
    }, 
    error = function(cnd){
      tictoc::toc()
      NULL
    })
    
  }
  
  
  res_tbl <- tribble(~name,~family_, ~comps_,
                     "poisson_pca", "poisson", list(comp_pca),
                     "poisson_pca_linear", "poisson", list(comp_pca_linear),
                     "poisson_pls", "poisson", list(comp_pls),
                     "poisson_pls_linear", "poisson", list(comp_pls_linear),
                     "zip_pca", "zeroinflatedpoisson1", list(comp_pca),
                     "zip_pca_linear", "zeroinflatedpoisson1", list(comp_pca_linear),
                     "zip_pls", "zeroinflatedpoisson1", list(comp_pls),
                     "zip_pls_linear", "zeroinflatedpoisson1", list(comp_pls_linear),
                     'poisson_sp',  "poisson", list(comp_simple),
                     "zip_sp", "zeroinflatedpoisson1", list(comp_simple)) #
  if(run_a2){
    # res_sp <- run_mod('poisson', comp_simple)
    # res <- run_mod('poisson', comp)
    res_tbl <- filter(res_tbl, family_ == "poisson") 
  }
  
  res_map <- pmap(res_tbl, run_mod) 
  names(res_map) <- res_tbl$name
  future::plan(future::sequential)
  
  
  get_waic <- function(mod_list, mod_names, criterion = "WAIC"){
    is_good <- map(mod_list, ~{!is.null(.x) & is.null(.x$error)}) |> 
      list_c()
    model <- mod_list[is_good]
    names <- mod_names[is_good]
    nmod <- sum(is_good)
    dic <- waic <- rep(NA, nmod)
    for (i in 1:nmod) {
      mod <- model[[i]]
      if (("DIC" %in% criterion) && is.null(mod$dic$dic)) {
        stop("Object ", i, " does not have DIC information.")
      }
      if (("WAIC" %in% criterion) && is.null(mod$waic$waic)) {
        stop("Object ", i, " does not have WAIC information.")
      }
      dic[i] <- mod$dic$dic
      waic[i] <- mod$waic$waic
    }
    if ("DIC" %in% criterion) {
      ord <- order(dic)
    }
    else if ("WAIC" %in% criterion) {
      ord <- order(waic)
    }
    dic <- dic[ord]
    waic <- waic[ord]
    names <- names[ord]
    ddic <- dic - min(dic)
    dwaic <- waic - min(waic)
    result <- data.frame(Model = names)
    if ("DIC" %in% criterion) {
      result <- cbind(result, data.frame(DIC = dic, Delta.DIC = ddic))
    }
    if ("WAIC" %in% criterion) {
      result <- cbind(result, data.frame(WAIC = waic, Delta.WAIC = dwaic))
    }
    return(result)
    
    
  }
  
  
  gc()
  
  # Model testing -----------------------
  waic_res <-  get_waic(res_map, res_tbl$name)
  write_csv(waic_res,g("{out_dir}/waic_{spp}_with_linear.csv") )
  res <- res_map[[waic_res$Model[[1]]]]
  
  if(!run_a2){
    ff_counts <- paste0("Intercept + RecLength + o + t2se + time_group + doy +",
                        paste0("spat_cov_",pc_vars, collapse = "+"),
                        "+ alpha") #|> 
    # as.formula()
    
    
    ff_complex <-  paste0("~{",
                          "expect <- exp(",ff_counts, ")\n",
                          "list(expect = expect,",
                          # "obs_prob = dpois(y, expect),",
                          "pobs = 1-dpois(0, expect))}") |> 
      as.formula()
    # 
    # if(str_detect(waic_res$Model[[1]], "pca")
    # f_pred <- 
    #   paste0("~{",
    #          "expect <- exp(",
    #          comp_str |> str_remove_all("\\([^\\)]+\\)")  |> str_remove_all("\\n") |> 
    #            str_remove_all("\\)") |> str_remove("~"), ")\n",
    #          "list(
    #     expect = expect,
    #     obs_prob = dpois(y, expect))}") |> 
    #   as.formula()
    
    
    x4pred_t2se <- expand_grid(      t2se_sc = seq(-1.,1., length.out = 100),
                                     event_gr = unique(df_std$event_gr),
                                     offset = offsets_spp |> 
                                       filter(dur==5 & max_dist ==Inf) |> 
                                       distinct(o) |> pull(o),
                                     RL = 5)
    pred_t2se <- predict(res,
                         x4pred_t2se,
                         ~{expect <-  Intercept + o + time_group + t2se + RecLength
                         1-dpois(0, exp(expect))
                         }
                         ,
                         n.samples = 200
    )  |> as_tibble() |> left_join(t2se_sd_mn,
                                   by = join_by(event_gr )) |> 
      mutate(t2se_modelled = t2se_sc,
             t2se_sc  = t2se_sc * (2*sd_t2se) + mean_t2se)
    
    
    time_plt <- 
      ggplot() +
      # gg(pred_t2se, aes(fill = event_gr, colour = event_gr,group = event_gr)) #+
      gg(pred_t2se, aes(colour  = event_gr)) +  
      facet_wrap(~event_gr, scales='free') +
      scale_fill_viridis_d(direction = -1) +
      scale_colour_viridis_d(direction = -1) +
      labs(x = "Time to sun event", y = "Prob of obs (5min) plus CI",
           colour = "",
           fill = "")
    
    x4pred_doy <-  expand_grid(      doy_r = seq(-1.,1., length.out = 100),
                                     offset = offsets_spp |> 
                                       filter(dur==5 & max_dist ==Inf) |> 
                                       distinct(o) |> pull(o),
                                     RL = 5)
    
    pred_doy <- predict(res,
                        x4pred_doy,
                        ~{expect <-  Intercept +o+  doy + RecLength
                        1-dpois(0, exp(expect))
                        },
                        
                        n.samples = 200
    ) |> as_tibble() |> bind_cols(sd_mn |> filter(variable =="doy")) |> 
      mutate(doy_r  = ymd("2025-01-01") -1 + doy_r * (sd2) + mn)
    
    doy_plt <- 
      ggplot() +
      # gg(pred_t2se, aes(fill = event_gr, colour = event_gr,group = event_gr)) #+
      gg(pred_doy) +  
      # facet_wrap(~event_gr, scales='free') +
      scale_fill_viridis_d(direction = -1) +
      scale_colour_viridis_d(direction = -1) +
      labs(x = "Time to sun event", 
           y = "Prob of obs (5min) plus CI",
           colour="",
           fill = "")
    
    
    
    x4 <- seq(-2,1.5, length.out = 100)
    
    
    
    
    
    
    t2se_mod <- pred_t2se |> 
      # filter( event_gr=="Dawn") |>
      summarize(median = median(median), .by = event_gr)
    write_rds(t2se_mod, g("{out_dir_tmp}/t2se_mod.rds"))
    # slice_max(by = event_gr,order_by = median, n=1, with_ties = F) 
    
    doy_mod <- pred_doy |> 
      slice_max(order_by = mean, n=1, with_ties = F) |> 
      pull(doy_r)
    
    ## REMOVING PCA PLTS FOR NOW ---------------------
    pred_pc <- function(mod,PC_var, f_var, nsamp=200){
      xx <- range(pull(df_std,{{PC_var}}) )
      xx4 <- seq(xx[1], xx[2], length.out = 100)
      pred_c1 <- predict(mod,tibble({{PC_var}}:= xx4),
                         as.formula(paste0(" ~ exp(Intercept +", f_var, ")")),
                         n.samples = nsamp)
      
      
      ggplot() +
        # gg(pred_t2se, aes(fill = event_gr, colour = event_gr,group = event_gr)) #+
        gg(pred_c1)
    }
    
    pcn <- str_pad(1:10, 2, 'left', 0)
    
    pcn_map <-
      tibble(PC_var = glue::glue("PC{pcn}"),
             f_var = glue::glue("spat_cov_PC{pcn}"))
    
    pca_plts <- pmap(pcn_map, pred_pc)
    
    ggsave(plot = time_plt,g(
      "{out_dir}/time_plot.jpeg"
    ))
    ggsave(plot = doy_plt,g(
      "{out_dir}/date_plot.jpeg"
    ))
  }
  
  # cairo_pdf(g("{out_dir}/PCA_patterns_{spp}.pdf") )
  # pca_plts
  # dev.off()
  
  
  #   # geom_point(data = cd, aes(x = x, y = count / exposure), cex = 2) +
  #   # geom_point(data = true.lambda, aes(x, y), pch = "_", cex = 9, col = "blue") +
  #   # coord_cartesian(xlim = c(0, 55), ylim = c(0, 6)) +
  #   xlab("x") +
  #   ylab("Intensity") +
  #   facet_wrap(~event_gr, scales = 'free')
  
  write_rds(res,  g("{out_dir_tmp}/{spp}_inlabru_model.rds"))  
  write_rds(df_std, g("{out_dir_tmp}/{spp}_inlabru_model_data.rds"))  
  gc()
  
}



