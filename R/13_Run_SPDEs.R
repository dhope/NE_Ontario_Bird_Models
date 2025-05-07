source("R/08_INLA_SPDE.R")

spp_to_run <- c(
  "LEYE",
  "GRYE",
  "PAWA",
  "BLPW",
  "CONW",
  "LISP",
  "FOSP",
  "LEOW",
  "GGOW",
  "ATSP",
  "OSFL",
  "RUBL",
  "LCSP",
  "NHOW",
  "NESP",
  "BOCH",
  "SBDO",
  "CONI"

)

srs <- purrr::safely(run_inlabru)
x <- purrr::walk(spp_to_run,
                 srs )
