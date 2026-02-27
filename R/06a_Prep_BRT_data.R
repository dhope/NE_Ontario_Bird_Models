library(tidyverse)
library(sf)
library(tidymodels)
library(offsetreg)
library(patchwork)
library(tictoc)

library(DALEXtra)
run_a2 <- FALSE

# Load in data shared across functions
source("R/__globals.R")
ggplot2::theme_set(ggplot2::theme_minimal(
  base_size = 14,
  base_family = "Roboto Condensed"
))
spp_cov_date <- "2026-02-16"
spp_cov_file <- g("output/rds/{spp_cov_date}_spatial_covariates_data.rds") 
# Load in area of interest from data mesh
aoi <- read_rds(g("output/rds/{date_compiled}_Area_of_focus_mesh.rds"))
# Load data ---------
source("R/__data_prep.R")

all_spp <- unique(counts$species_name_clean)
all_dat <- map(all_spp, prep_brt_data, .progress = T)
df_stds <- all_dat |> transpose() |> pluck("df_std")

names(df_stds) <- all_spp



walk(all_spp, ~{
  df <- df_stds[[.x]]
  if(nrow(df)==0 || sum(df$y)==0) return(NULL)
  spp_df <- na.omit(df$species_name_clean)[[1]]
  spp_ <- spp_df |> str_replace_all(" ", "_")
  if(.x!= spp_df) abort("Bad_name")
  
  
  write_rds(df, g("{brt_spp_dat_loc}/{spp_}.rds")) })
