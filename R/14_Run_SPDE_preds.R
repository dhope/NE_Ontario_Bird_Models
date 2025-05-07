source("R/09_INLA_preds.R")


spp_comp <- 
list.files(INLA_output_loc_TMP, "_inlabru_model.rds", recursive = T, full.names = F) |> 
  str_extract("\\w{4}")  |> unique() |> sort()#|> str_subset("CONW", negate = T)

# spp_to_run <- c(
#   "LEYE",
#   "GRYE",
#   "PAWA",
#   "BLPW",
#   "CONW",
#   "LISP",
#   "FOSP",
#   "LEOW",
#   "GGOW",
#   "ATSP",
#   "OSFL",
#   "RUBL",
#   "LCSP",
#   "NHOW",
#   "NESP",
#   "BOCH",
#   "SBDO",
#   "CONI"
#   
# )
# srs <- purrr::safely(run_predictions)
# x <- purrr::walk(spp_comp,
#                  srs )
# prep_predictions("PAWA", save_objects = T)
# run_predictions("PAWA")


# s_prep <- purrr::safely(prep_predictions)
# s_run_preds <- purrr::safely(run_predictions)
# x <- purrr::walk(spp_comp,
#                  s_prep, save_objects = TRUE )
# cat("Data prep complete --------------------------\n")

# x1 <- purrr::walk(spp_comp,
#                   s_run_preds )
# s_gm <- safely(gen_maps)
# x2 <- purrr::walk(c("GRYE", "GGOW","LEOW", "LEYE", "LISP", "NHOW", "OSFL", "PAWA", "RUBL", "SBDO" ), s_gm)
gen_maps("RUBL")