source("R/__globals.R")
library(napops)
library(tidyverse)

library(sf)

counts <- read_rds(g("{rds_data_loc}/counts.rds")) |>
  filter(location %in% aggregated_locs$location) |>
  filter(
    collection %in%
      c("WildTrax", "ONATLAS3PC") &
      str_detect(location, "WHA-", negate = T) &
      str_detect(project, "Nocturnal", negate = T)
  ) #|>
napops_spp <- list_species()


full_recordings <- g("{rds_data_loc}/all_events.rds") |>
  read_rds() |>
  filter(
    location %in% aggregated_locs$location,
    collection %in% c("ONATLAS3PC", "WildTrax")
  )


spp <- tibble(english_name = unique(counts$species_name_clean)) |>
  left_join(
    naturecounts::meta_species_taxonomy() |>
      filter(concept_source == "avibase"),
    by = join_by(english_name)
  ) |>
  left_join(
    naturecounts::meta_species_codes() |>
      filter(authority == "BSCDATA"),
    by = join_by(species_id)
  ) |>
  slice_min(rank, by = english_name, with_ties = F)

dur <- full_recordings$clip_length_min |> round(2) |> unique()


get_na_pops <- function(spp) {
  if (spp %in% napops_spp$Species) {
    na_info <- napops_spp |>
      filter(Species == spp) |>
      dplyr::select(Removal, Distance)
    if (na_info$Removal == 1) {
      cr <- cue_rate(spp, model = 1, od = 150, tssr = 0)$CR_est[[1]] #,
      #              samples = 1000,
      #              quantiles = c(0.025,0.15865,0.84135,  0.975)) |>
      # mutate(sd = CR_84.135 - CR_est)
    } else {
      cr <- NA
    }
    if (na_info$Distance == 1) {
      edr_sp <- edr(spp, model = 1, road = F, forest = 0.5) #,
      # samples = 1000,
      # quantiles = c(0.025,0.15865,0.84135,  0.975)) |>
      # mutate(sd = EDR_84.135 - EDR_est)
    } else {
      edr_sp <- tibble(EDR_est = NA)
    }

    # Calculate A and p, which jointly determine offset
    edr_est_ha <- edr_sp$EDR_est[[1]] / 100
    calc_values <- function(edr, cr, time_minutes, max_dist) {
      inf <- max_dist == Inf
      # edr_ha <- edr/100
      # browser()
      tibble(
        edr,
        cr,
        max_dist,
        time_minutes,
        A = if_else(inf, pi * edr**2, pi * max_dist**2), # c(pi*c(edr_est_ha, 1.0)^2)
        # ifelse(unlim, pi * cf0[2]^2, pi * MAXDIS[!OKq]^2)
        p = 1 - exp(-time_minutes * cr),
        q = if_else(
          inf,
          1,
          edr^2 * (1 - exp(-max_dist^2 / edr^2)) / max_dist^2
        ),
        o = log(A) + log(p) + log(q)
      )
      # q <- c(1,  edr_est_ha^2 * (1 - exp(-1^2/edr_est_ha^2))/1^2)
      # o <- log(A) + log(p) + log(q)
    }
    # browser()
    expand_grid(
      edr = edr_est_ha,
      cr = cr,
      time_minutes = dur,
      max_dist = c(0.5, 1, Inf)
    ) |>
      rowwise() |>
      mutate(
        qpad = list(calc_values(edr_est_ha, cr, time_minutes, max_dist))
      ) |>
      ungroup() |>
      dplyr::select(qpad) |>
      unnest(qpad) |>
      mutate(spp = spp)
  } else {
    tibble(spp)
  }
}

na_pops_offsets <- map(unique(spp$species_code), get_na_pops) |>
  list_rbind() |>
  left_join(
    spp |> distinct(english_name, species_code),
    by = join_by(spp == species_code)
  )


write_rds(na_pops_offsets, "output/na_pops_offsets.rds")
