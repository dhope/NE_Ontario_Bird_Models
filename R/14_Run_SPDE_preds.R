run_a2 <- FALSE
source("R/09_INLA_preds.R")




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
  # str_subset("American_Bittern", negate = T) |>
  # str_subset("Rusty_Blackbird", negate = T) |>
  str_extract(glue::glue("{ifelse(run_a2, '(?<=A2/)','^')}\\w+(_)*\\w+(?=/)")) |>
  str_replace_all("_", " ") |> 
  unique() |>
  sort()  #|> str_subset("CONW", negate = T)


comp <- 
list.files("D:/SPATIAL/RoF_Models/INLA3", ".tif", recursive=T) |> 
  str_subset("^A2/", negate = !run_a2) |>
  str_subset("^A2/\\w{4}/", negate = T) |>
  str_subset("\\s", negate = T) |>
  str_subset("CONI_prev", negate = T) |>
  str_extract(glue::glue("{ifelse(run_a2, '(?<=A2/)','^')}\\w+(_)*\\w+(?=/)")) |>
  str_replace_all("_", " ") |> 
  unique() |>
  sort()

missed <- spp_comp_inla[!spp_comp_inla %in% comp]
# run_predictions("Rusty Blackbird", load_rds = F, gen_map_outputs = TRUE)

rp <- safely(run_predictions)
x <- map(missed, rp,load_rds = F, gen_map_outputs = TRUE)
sink("output_pred_inla.log")
print(x)
sink()



# 
# 
# spp_comp <- 
# list.files(INLA_output_loc_TMP, "_inlabru_model.rds", recursive = T, full.names = F) |> 
#   str_extract("\\w{4}")  |> unique() |> sort() |> str_subset("CONW", negate = T)
# 
# # 
# # pred_comps <- list.files("G:/RoF_Model_Preds_rds/INLA2", "predictions_full.rds", recursive = T, full.names = F) |> 
# #   str_extract("\\w{4}")  |> unique() |> sort() |> str_subset("PAWA", negate = T)
# pred_comps <- list.files("D:/SPATIAL/RoF_Models/INLA2", "mad_", recursive = T, full.names = F) |> 
#   str_extract("\\w{4}")  |> unique() |> sort() 
# # spp_to_run <- c(
# #   "LEYE",
# #   "GRYE",
# #   "PAWA",
# #   "BLPW",
# #   "CONW",
# #   "LISP",
# #   "FOSP",
# #   "LEOW",
# #   "GGOW",
# #   "ATSP",
# #   "OSFL",
# #   "RUBL",
# #   "LCSP",
# #   "NHOW",
# #   "NESP",
# #   "BOCH",
# #   "SBDO",
# #   "CONI"
# #   
# # )
# 
# to_run <- spp_comp[!spp_comp %in% pred_comps]
# 
# # run_predictions("CONI", load_rds = F, gen_map_outputs = TRUE)
# 
# 
# srs <- purrr::safely(run_predictions)
# x <- map(to_run,  srs, load_rds = F, gen_map_outputs = TRUE)
# # x <- purrr::walk(spp_comp,~{
# #                  prep_predictions(.x, save_objects = T)
# #                  srs(.x , load_rds = TRUE) })
# 
# 
# # prep_predictions("PAWA", save_objects = T)
# # run_predictions("PAWA")
# 
# 
# # s_prep <- purrr::safely(prep_predictions)
# # s_run_preds <- purrr::safely(run_predictions)
# # x <- purrr::walk(spp_comp,
# #                  s_prep, save_objects = TRUE )
# # cat("Data prep complete --------------------------\n")
# 
# # x1 <- purrr::walk(spp_comp,
# #                   s_run_preds )
# # s_gm <- safely(gen_maps)
# # # x2 <- purrr::walk(c("GRYE", "GGOW","LEOW", "LEYE", "LISP", "NHOW", "OSFL", "PAWA", "RUBL", "SBDO" ), s_gm)
# # # run_predictions("CONW", load_rds = TRUE)
# # # gen_maps("RUBL")
# # 
# # 
# # pred_comps
# # x <- purrr::walk(pred_comps,s_gm)