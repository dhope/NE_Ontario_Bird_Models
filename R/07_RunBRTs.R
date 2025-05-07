source("R/06_BRT_time.R")
spp_codes <- 
summarize(counts, n=n(), n_sites = n_distinct(location), .by = c(species_name_clean)) |> 
  arrange(desc(n_sites)) |> left_join(
distinct(counts,species_name_clean, species_code) |> filter(!is.na(species_code)),
by = join_by(species_name_clean)) |> 
  filter(!is.na(species_code) & n_sites>10 )


srb <- purrr::safely(run_brt)
srbp <- purrr::safely(run_brt_preds)
x <- purrr::walk(spp_codes$species_code,srb )
x1 <- purrr::walk(spp_codes$species_code,srbp )