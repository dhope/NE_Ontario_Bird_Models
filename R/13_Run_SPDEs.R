source("R/08_INLA_SPDE.R")

spp_comp <- 
  list.files(INLA_output_loc_TMP, "_inlabru_model.rds", recursive = T, full.names = F) |> 
  str_extract("\\w{4}")  |> unique() |> sort()#|> str_subset("CONW", negate = T)


spp_codes <- 
  summarize(counts, n=n(), n_sites = n_distinct(location), .by = c(species_name_clean)) |> 
  arrange(desc(n_sites)) |> left_join(
    distinct(counts,species_name_clean, species_code) |> filter(!is.na(species_code)),
    by = join_by(species_name_clean)) |> 
  filter(!is.na(species_code) & n_sites>10 & !species_code %in% spp_comp )


srs <- purrr::safely(run_inlabru)
x <- purrr::walk(spp_codes$species_code,
                 srs )
