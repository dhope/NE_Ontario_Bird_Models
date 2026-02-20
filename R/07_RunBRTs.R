source("R/06_BRT_time.R")
rm_spp_pat <- c("Gull", "Loon", "Mallard", "Swan", "Harrier", "Eagle", "Teal", "Pintail", "Goldeneye", "Merganser", "Duck", "Merlin",
                "Scoter", "Tern", "Wigeon")
spp_to_run <-
  summarize(
    counts,
    n = n(),
    n_sites = n_distinct(location),
    .by = c(species_name_clean, common_id)
  ) |>
  arrange(desc(n_sites)) |> 
  filter(str_detect(species_name_clean, glue::glue_collapse(rm_spp_pat, sep = "|"), negate=T) &n_sites>10)

# spp_completed <- list.files(
#   path = bundle_locs,
#   "\\w{4}_time",
#   recursive = T
# ) |>
#   str_extract("\\w{4}") |>
#   unique()
# spp_to_run <- spp_codes_full_run$species_code[
#   !spp_codes_full_run$species_code %in% spp_completed
# ]
library(future)
options(future.globals.maxSize = 16000 * 1024^2)
# plan(multisession, workers = 32)
# mirai::daemons(32)
plan(future.mirai::mirai_multisession, workers = 32, queue = TRUE)

# library(mirai)
# daemons(32)
srb <- purrr::safely(run_brt)
# , otherwise = {
#   plan(sequential)
#   plan(multico, workers = 32)}
x <- purrr::map(spp_to_run$species_name_clean[[1]], srb)
# plan(sequential)
# transpose(x) |> pluck('error')
