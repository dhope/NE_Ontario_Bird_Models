run_a2 <- FALSE
source("R/08_INLA_SPDE.R")

spp_comp <-
  list.files(
    BRT_output_loc,
    "time_model_pred.rds",
    recursive = T,
    full.names = F
  ) |>
  str_extract("\\w{4}") |>
  unique() |>
  sort() #|> str_subset("CONW", negate = T)

spp_more_10 <-
  summarize(
    counts,
    n = n(),
    n_sites = n_distinct(location),
    .by = c(species_name_clean, species_scientific_name)
  ) |>
  arrange(desc(n_sites)) |>
  left_join(
    distinct(counts, species_name_clean, species_code) |>
      filter(!is.na(species_code)),
    by = join_by(species_name_clean)
  ) |>
  filter(n_sites > 10)


if (run_a2) {
  spp_to_run <- spp_comp
} else {
  spp_to_run <- spp_more_10 |>
    filter(!is.na(species_code) & species_code %in% spp_comp) |>
    pull(species_code)
}


srs <- purrr::safely(run_inlabru)
x <- purrr::walk(spp_codes$species_code, srs)
