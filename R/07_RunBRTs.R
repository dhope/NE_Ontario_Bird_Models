source("R/06_BRT_time.R")
rm_spp_pat <- c("Gull", "Loon", "Mallard", "Swan", "Harrier", "Eagle", "Teal", "Pintail", "Goldeneye", "Merganser", "Duck", "Merlin",
                "Scoter", "Tern", "Wigeon", 
                "Yellowlegs", "Sandhill", "Goldfinch") # Already ran these in testing
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

gen_spp_mod_list <- function(spp){
  spp_dir <- str_replace_all(spp, ' ', '_')
  df_std <- read_rds(g("{brt_spp_dat_loc}/{spp_dir}.rds"))
  n_sites <- df_std |> filter(y>0) |> pull(site_id) |> n_distinct()
  
  out_dir_spatial <- g(
    "{BRT_output_loc_spatial}/{spp_dir}"
  )
  out_dir <- g(
    "{BRT_output_loc}/{spp_dir}"
  )
  
  model_rds_file <- glue::glue("{bundle_locs}/{spp_dir}_time_bundle.rds")
  model_rds_file_old <- glue::glue("{bundle_locs}/{spp}_time_bundle.rds")
  
  if(n_sites<=20){
    # cat(glue::glue("Insufficient data, skipping {spp}\n\r"))
    # file.create(glue::glue("{out_dir}/skipped_due_insufficient_data.txt"))
    
    return(NULL)}
  if(file.exists(model_rds_file)||file.exists(model_rds_file_old)){
    # cat(glue::glue("Model rds already exists, skipping {spp}\n\r"))
    # file.create(glue::glue("{out_dir}/skipped_due_to_existing_rds.txt"))
    return(NULL)
  }
  
  spp
  
}

running_spp <- map(spp_to_run$species_name_clean, gen_spp_mod_list)
list_c(running_spp) |> write_lines("specieslist.txt")


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
# library(future)
# options(future.globals.maxSize = 8000 * 1024^2)
# if(str_detect(osVersion, "Windows")){
#   plan(multisession, workers = 32, gc=T)
# } else{
#   plan(multicore, workers = 32, gc=T)
# }

# mirai::daemons(32)
# plan(future.mirai::mirai_multisession, workers = 32, queue = TRUE)
# 
# # library(mirai)
# # daemons(32)
# srb <- purrr::safely(run_brt , otherwise = {
#   plan(sequential)})
# 
# 
# 
# 
# #   plan(multico, workers = 32)}
# # run_brt("Sandhill Crane")
# x <- purrr::map(spp_to_run$species_name_clean, srb)
# plan(sequential)
# transpose(x) |> pluck('error')
