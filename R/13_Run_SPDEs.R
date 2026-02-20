run_a2 <- TRUE
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

core <- read_csv(g("{core_path}/ECCC_Avian_Core_20241025.csv")) |> 
  janitor::clean_names() |> 
  filter(str_detect(technical_committees, "WF", negate = T ) & cdn_status!="DNO" )

spp_comp[spp_comp %in% core$species_id]
# if (run_a2) {
#   spp_to_run <- spp_comp
# } else {
#   spp_to_run <- spp_comp#spp_more_10 |>
#     # filter(!is.na(species_code) & species_code %in% spp_comp) |>
#     # pull(species_code) |> unique()
# }


srs <- purrr::safely(run_inlabru)
time_limit <- 0.5 * 60 * 60
future::plan(future::mulitcore, workers = 32)
inla.setOption(inla.timeout = time_limit, num.threads = 32)
x <- purrr::walk(spp_comp, srs, .progress = T)
future::plan(future::sequential)
# a2comp <- 
# list.files("D:/MODEL_OUTPUT/RoF/INLA3/A2", "model_parameter_summary", recursive = T) |> 
#   str_extract("\\w{4}") 
# 
# 
# x <- purrr::walk(spp_to_run[!(spp_to_run %in% a2comp)], srs)
