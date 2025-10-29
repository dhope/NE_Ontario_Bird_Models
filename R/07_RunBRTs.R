source("R/06_BRT_time.R")
spp_codes_full_run <-
  summarize(
    counts,
    n = n(),
    n_sites = n_distinct(location),
    .by = c(species_name_clean)
  ) |>
  arrange(desc(n_sites)) |>
  left_join(
    distinct(counts, species_name_clean, species_code) |>
      filter(!is.na(species_code)),
    by = join_by(species_name_clean)
  ) |>
  filter(!is.na(species_code) & n_sites > 10)

spp_completed <- list.files(
  path = BRT_output_loc,
  "\\w{4}_time_map",
  recursive = T
) |>
  str_extract("\\w{4}") |>
  unique()
spp_to_run <- spp_codes_full_run$species_code[
  !spp_codes_full_run$species_code %in% spp_completed
]

srb <- purrr::safely(run_brt)
x <- purrr::walk(spp_to_run, srb)
