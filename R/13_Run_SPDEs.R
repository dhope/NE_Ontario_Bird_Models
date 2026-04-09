run_a2 <- FALSE
source("R/08_INLA_SPDE.R")

rm_spp_pat <- c("Gull", "Loon", "Mallard", "Swan", "Harrier", "Eagle", "Teal", "Pintail", "Goldeneye", "Merganser", "Duck", "Merlin",
                "Scoter", "Tern", "Wigeon", "Gadwall") # Already ran these in testing
# spp_to_run <-
#   summarize(
#     counts,
#     n = n(),
#     n_sites = n_distinct(location),
#     .by = c(species_name_clean, common_id)
#   ) |>
#   arrange(desc(n_sites)) |> 
#   filter(str_detect(species_name_clean, glue::glue_collapse(rm_spp_pat, sep = "|"), negate=T) &n_sites>10) |> 
#   filter(species_name_clean!="Song Sparrow")
spp_to_run <- tibble(file_size=brt_spp_dat_loc |> list.files(full.names = T) |> file.size(),
                     species_name_clean=
                       brt_spp_dat_loc |> list.files() |> str_remove(".rds") |> str_replace_all("_", " ") ) |> arrange(file_size) |>
  filter(str_detect(species_name_clean, glue::glue_collapse(rm_spp_pat, sep = "|"), negate=T))



spp_comp <-
  list.files(
    BRT_output_loc,
    "time_model_pred.rds",
    recursive = T,
    full.names = F
  ) |>
  str_subset("^\\w{4}/", negate = T) |>
  str_subset("\\s", negate = T) |>
  str_subset("CONI_prev", negate = T) |>
  str_extract("^\\w+(_)*\\w+(?=/)") |>
  unique() |>
  sort() |> #|> str_subset("CONW", negate = T)
  str_replace_all("_", " ")

# 
spp_comp_inla <-
  list.files(
    INLA_output_loc_TMP,
    "_inlabru_model.rds",
    recursive = !run_a2,
    full.names = F
  ) |>
  str_subset("^A2/", negate = !run_a2) |>
  str_subset("^A2/\\w{4}/", negate = T) |>
  str_subset("\\s", negate = T) |>
  str_subset("CONI_prev", negate = T) |>
  str_extract(glue::glue("{ifelse(run_a2, '(?<=A2/)','^')}\\w+(_)*\\w+(?=/)")) |>
  unique() |>
  sort() |> #|> str_subset("CONW", negate = T)
  str_replace_all("_", " ")

spp_brt_run <- spp_to_run$species_name_clean[str_remove_all(spp_to_run$species_name_clean, "\\'") %in% spp_comp] |> 
  sort()




missed <- spp_brt_run[!str_remove(spp_brt_run, "\\'") %in% spp_comp_inla]
# core <- read_csv(g("{core_path}/ECCC_Avian_Core_20241025.csv")) |>
#   janitor::clean_names() |>
#   filter(str_detect(technical_committees, "WF", negate = T ) & cdn_status!="DNO" )
#
# spp_comp[spp_comp %in% core$species_id]
# if (run_a2) {
#   spp_to_run <- spp_comp
# } else {
#   spp_to_run <- spp_comp#spp_more_10 |>
#     # filter(!is.na(species_code) & species_code %in% spp_comp) |>
#     # pull(species_code) |> unique()
# }

srs <- purrr::safely(run_inlabru)
time_limit <- 0.5 * 60 * 60
if (str_detect(osVersion, "Windows")) {
  future::plan(future::multisession, workers = 32)
} else {
  future::plan(future::multicore, workers = 32)
}
inla.setOption(inla.timeout = time_limit, num.threads = 32)
x <- purrr::map(spp_brt_run, srs, .progress = T)
future::plan(future::sequential)
sink("output.log")
print(x)
sink()
# a2comp <-
# list.files("D:/MODEL_OUTPUT/RoF/INLA3/A2", "model_parameter_summary", recursive = T) |>
#   str_extract("\\w{4}")
#
#
# x <- purrr::walk(spp_to_run[!(spp_to_run %in% a2comp)], srs)
