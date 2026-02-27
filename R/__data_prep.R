## To prepare data in a standard way between INLA and BRTs
## Used in 06_BRT_time.R & 08_INLA_SPDE.R 

# Spatial covariates compiled in )3_Site_Data.R
spatial_cov <-
  spp_cov_file |>
  read_rds() %>%  {
    if (!run_a2) {
      mutate(., d2O = read_rds(g("output/rds/{spp_cov_date}_dist2ocean.rds")))
    } else {
      mutate(., d2O = read_rds(g("output/rds/{spp_cov_date}_dist2ocean_a2.rds")))
    }
  } |>
  bind_cols(dplyr::select(
  individual_locs,
  site_id, geometry)) |> #Spatial_covariates_data_14March2024.rds") |>
  distinct() |>
  st_as_sf() %>%
  bind_cols(as_tibble(st_coordinates(.))) |>
  dplyr::select(where(~ sum(is.na(.x)) != length(.x))) #|>
  # left_join(st_drop_geometry(aggregated_locs), by = join_by(site_id_agg))

# Testing and training data compiled in 01_Withold_test_sites.R
test_training_data <- read_rds(g(
  "output/{date_compiled}_test_training_data.rds"
))


if(!run_a2){
raw_recordings <- test_training_data$train_recordings |>
  # left_join(aggregated_locs, by = join_by(project, location, collection)) |>
  filter(str_detect(
    project,
    "(Extraction)|(Nocturnal)|(Resample)",
    negate = T
  ))

# Load in prepared count data
counts_full <- read_rds(g("{rds_data_loc}/counts.rds")) |> mutate(event_id = glue::glue("wt-{event_id}"))
# spp_list <- distinct(counts_full, species_name_clean, species_code) |>
#   filter(!is.na(species_code))

raw_counts <- counts_full |>
  replace_na(list(collection = "WildTrax")) |>
  filter(event_id %in% raw_recordings$event_id) |>
  filter(str_detect(
    project,
    "(Extraction)|(Nocturnal)|(Resample)",
    negate = T
  ))
rm(counts_full)
}


recordings <- raw_recordings |>
  # Assign recordings to day period
  mutate(
    Time_period = dplyr::case_when(
      t2ss >= -60 & t2ss <= 150 ~ "Dusk",
      t2sr >= -70 & t2sr <= 220 ~ "Dawn",
      t2ss > 150 & t2sr < -70 ~ "Night",
      abs(t2sr) > 220 ~ "Day",
      is.na(t2sr) ~ "Missing time or location data",
      TRUE ~ "Unk"
    ),
    dur = round(clip_length_min, 2),
    t2se = dplyr::case_when(
      Time_period == "Dusk" ~ t2ss,
      Time_period == "Dawn" ~ t2sr,
      TRUE ~ pmin(abs(t2ss), abs(t2sr))
    ),
    doy = yday(date),
    recording_id = as.numeric(recording_id)
  ) 

counts <- raw_counts |>
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
  mutate(
    y = ifelse(is.na(total_count_with_tmtt), total_count, total_count_with_tmtt)
  ) |>
  filter(!is.na(y))

rm(test_training_data,  raw_recordings, raw_counts)


qpad_offsets <- read_rds("output/QPAD_global_offsets.rds") |>
  rename(max_dist = r, time_minutes = t)
na_pops_offsets <- read_rds("output/na_pops_offsets.rds") |>
  rename(species = spp)


prep_brt_data <- function(spp) {
  flush.console()
  
  # spp_name <- spp_list$species_name_clean[spp_list$species_code == spp]
  spp_name <- spp_codes |> filter(species_name_clean==spp & spp_group == "Bird species") |> 
    distinct(species_name_clean, common_id)
  all_relavent_nnames <- filter(spp_codes, species_name_clean == spp |
                                  common_id %in% spp_name$common_id ) |> 
    distinct(species_name_clean,species_common_name, common_id)
  # counts$species_name_clean[counts$species_code==spp & !is.na(counts$species_code)] |>
  # unique()
  if (nrow(spp_name) != 1) {
    rlang::abort("unable to identify species name")
  }
  # spp_name_sci <- naturecounts::search_species_code(
  #   spp_name$common_id,
  #   results = 'exact'
  # )$scientific_name[[1]]
  safe_dates_spp <- filter(
    safe_dates,
    (english_name %in% all_relavent_nnames$species_common_name | str_detect(species_id,spp_name$common_id)) & biol_region %in% c(4, 5) & level == 2
  ) |>
    mutate(
      on_er = case_when(
        biol_region == 4 ~ "Ontario Shield",
        biol_region == 5 ~ "Hudson Bay Lowlands",
        TRUE ~ NA_character_
      )
    ) |>
    summarize(
      start_doy = min(start_dt_julian),
      end_doy = max(end_dt_julian),
      .by = c(on_er, level)
    )
  if(nrow(safe_dates_spp)!=2){
    safe_dates_spp <-  expand_grid(on_er = c("Ontario Shield","Hudson Bay Lowlands"),
                                   safe_dates_spp |> dplyr::select(-on_er)    )
    
  }
  
  counts_spp <- filter(
    counts,
    species_name_clean == spp |
      common_id %in% spp_name$common_id 
  ) |>
    full_join(
      recordings |> dplyr::select(-any_of("geometry")),
      by = join_by(event_id, location, project)
    ) |>
    replace_na(list(total_count_with_tmtt = 0, total_count = 0, y = 0))
  # offfiles <-g("output/species_rds/offsets_{spp}.rds")
  offsets_spp <- distinct(
    recordings,
    location,
    event_id,
    dur,
    species_code = spp
  ) |>
    mutate(o = NA)#, dur = round(dur,2))
  
  # if (spp %in% na_pops_offsets$species[!is.na(na_pops_offsets$o)]) {
  #   
  #   offsets_spp <- left_join(
  #     offsets_spp |> dplyr::select(-o),
  #     filter(na_pops_offsets, species == spp & max_dist == Inf),
  #     by = join_by(species_code == species, dur == time_minutes)
  #   )
  # } else {
  #   if (spp %in% qpad_offsets$species) {
  #     # offsets_spp <- filter(qpad_offsets, species == spp)
  #     offsets_spp <- left_join(
  #       offsets_spp |> dplyr::select(-o),
  #       filter(qpad_offsets, species == spp & max_dist == Inf),
  #       by = join_by(species_code == species, dur == time_minutes)
  #     )
  #   }
  # }
  
  if (!"max_dist" %in% names(offsets_spp)) {
    offsets_spp$max_dist <- Inf
  }
  
  out_dir_spatial <- g(
    "{BRT_output_loc_spatial}/{str_replace(spp, '\\\\s', '_')}"
  )
  out_dir <- g(
    "{BRT_output_loc}/{str_replace(spp, '\\\\s', '_')}"
  )
  dir.create(
    out_dir_spatial,
    recursive = T
  )
  dir.create(
    out_dir,
    recursive = T
  )
  
  setup_dat_0 <- counts_spp |>
    left_join(
      offsets_spp |>
        dplyr::select(event_id, offset = o),
      by = join_by(event_id)
    ) |>
    # filter(species_code==spp & str_detect(project, "Resample", negate=T) &
    #          recording_id %in% train_locs$recording_id) |>
    #
    # left_join(recordings |>
    #             dplyr::select(project:recording_id, date_time:t2se, dur),
    #           by = join_by(project, project_id, location, recording_id)
    # ) |>
    mutate(Rec_length = dur) |> #factor(as.numeric(round(dur, 1)))) |>
    mutate(
      year = factor(year),
      SiteN = as.numeric(factor(location)),
      rec_id = as.numeric(as.factor(event_id)),
      event_gr = factor(case_when(
        is.na(date_time) ~ NA_character_,
        Time_period %in% c("Dawn", "Dusk", "Night") ~ Time_period,
        Time_period == "Day" & hour(date_time) <= 10 ~ "Dawn",
        Time_period == "Day" & hour(date_time) > 10 ~ "Day",
        TRUE ~ "Other"
      ))
    ) |>
    arrange(location, doy, Time_period, t2se) |>
    mutate(Row = row_number())
  
  included_times <- setup_dat_0 |>
    filter(y > 0) |>
    janitor::tabyl(event_gr) |>
    filter(n > 0)
  
  if (nrow(included_times) == 0) {
    included_times <- setup_dat_0 |> #filter(y>0) |>
      janitor::tabyl(event_gr) |>
      filter(n > 0 & !is.na(event_gr))
  }
  
  setup_dat_nested <-
    setup_dat_0 |>
    filter(event_gr %in% included_times$event_gr) |>
    nest_by(event_gr)
  
  rm(setup_dat_0)
  
  setup_dat <-
    setup_dat_nested |>
    filter(event_gr != "Day") |>
    unnest(data)
  
  df_std <- setup_dat |>
    st_drop_geometry() |>
    ungroup() |>
    dplyr::select(
      -c(
        # project_id,
        common_id,
        longitude,
        latitude,
        SiteN,
        # species_name_clean,
        species_scientific_name,
        total_count,
        total_count_with_tmtt,
        species_code,
        recording_id, #event_id,
        tz,
        t2sr,
        t2ss,
        SamplingEventIdentifier,
        Time_period,
        month,
        day,
        path,
        # site_id,
        aru_id,
        date_time,
        date,
        clip_length_min,
        dur,
        Row,
        rec_id,
        move_type,geo_id,cluster_id,n,source, collection,ProtocolType,
        n_obs,n_min, dist_meters, id, #QPAD_offset
      )
    ) %>%
    mutate(random_val = rnorm(n = nrow(.))) |>
    left_join(
      spatial_cov |> st_drop_geometry(),
      by = join_by( site_id)
    ) |>
    # dplyr::select(-site_id_agg) |>
    rename(QPAD_offset = offset) |>
    left_join(
      select(
        safe_dates_spp,
        on_er,
        start_doy,
        end_doy
      ),
      by = join_by(on_er)
    ) |>
    filter(doy >= start_doy & doy <= end_doy) |>
    select(-start_doy, -end_doy, -on_er)
  
  as.list.environment(environment())
}



