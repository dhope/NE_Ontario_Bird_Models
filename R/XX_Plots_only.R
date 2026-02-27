library(tidyverse)
library(terra)
library(sf)
source("R/06_BRT_time.R")
ont.proj <- 3161
g <- glue::glue
just_maps <- function(spp){
  cat(g("{spp} -------------------\n"))
out_dir_spatial <- g(
  "{BRT_output_loc_spatial}/{spp}"
)
out_dir <- g(
  "{BRT_output_loc}/{spp}"
)
fit_workflow <- read_rds(
  file = glue::glue("{bundle_locs}/{spp}_time_bundle.rds")
)
df_std <- prep_brt_data(spp)$df_std

t_pred <-
  expand_grid(
    t2se = seq(-60, 240),
    event_gr = c("Dawn", "Dusk"),
    doy = 122:202
  ) |> bind_cols(df_std |> select(-t2se, -doy) |> 
                   summarize(across(where(is.numeric), \(x) median(x, na.rm=T)) )
  ) |> 
  bind_rows(df_std |> select(-t2se, -doy) ) |>
  filter(!is.na(t2se))

pred_t <- predict(fit_workflow, t_pred)
t_pred$.pred <- pred_t$.pred


period_to_use <- t_pred |> summarize(.pred = median(.pred), .by = c(event_gr)) |> 
  slice_max( order_by = .pred, 
             n=1,with_ties = F) |> 
  pull(event_gr)

doy_used <- t_pred |> filter( event_gr == period_to_use) |> 
  summarize(.pred = median(.pred), .by = c(doy)) |> 
  slice_max( order_by = .pred, 
             n=1,with_ties = F) |> 
  pull(doy)

if(!file.exists( glue::glue("{out_dir_spatial}/{spp}_time_brt.tif" ))){
  preds <- read_rds(g(
    "{prediction_layer_loc}/2025-10-22_prediction_rasters_df.rds"
  ))
  preds$d2O <- read_rds(g(
    "{prediction_layer_loc}/distance_2_ocean_prediction.rds"
  ))
  # preds_ez <- preds |>
  #   st_as_sf(coords = c("X", "Y"), crs = ont.proj) |>
  #   select() |>
  #   st_join(select(ontario_ez, on_er = ZONE_NAME))
  # preds$on_er <- preds_ez$on_er
  # rm(preds_ez)
  # xy <- read_rds("{prediction_layer_loc}/2025-02-28_prediction_rasters_xy.rds")
  
  # fit_workflow <- #bundle::unbundle(read_rds(
  #   readRDS(
  #     glue::glue("{bundle_locs}/{spp}_time_bundle.rds")
  #   )
  
  # 
  library(future)
  plan(multisession, workers = 32)
  p <- predict(
    fit_workflow,
    new_data = preds |>
      mutate(
        doy = doy_used,
        t2se = 0, #30,
        Time_period = period_to_use,
        year = factor(2026),
        event_gr = "Dawn",
        source = NA,
        collection = NA,
        type = NA,
        Rec_length = (5),
        SiteN = NA,
        longitude = NA,
        latitude = NA,
        event_id = NA,
        QPAD_offset = NA,
        olcb_99_100 = NA,
        total_100 = NA,
        olcb_99_500 = NA,
        total_500 = NA
      )
  )
  toc()
  
  plan(sequential)
  write_rds(p, glue::glue("{prediction_layer_loc}/{spp}_time_pred_brt.rds"))
  r_pred <- rast(g("{prediction_layer_loc}/2025-10-21_prediction_rasters.nc")) #2025-02-10_prediction_rasters.nc")
  pobs <- 1 - dpois(0, lambda = p$.pred)
  r_pobs <- predicted_raster <- rast(r_pred[[1]])
  r2 <- terra::cellFromXY(predicted_raster, preds[, c("X", "Y")])
  predicted_raster[r2] <- p$.pred
  r_pobs[r2] <- pobs
  predicted_raster <- mask(predicted_raster, aoi)
  r_pobs <- mask(r_pobs, aoi)
  names(r_pobs) <- names(predicted_raster) <- spp
  
  terra::writeRaster(
    predicted_raster,
    glue::glue("{out_dir_spatial}/{spp}_time_brt.tif"),
    overwrite = T
  )
  terra::writeRaster(
    r_pobs,
    glue::glue("{out_dir_spatial}/{spp}_pobs_time_brt.tif"),
    overwrite = T
  )
}


  predicted_raster <- terra::rast(
  glue::glue("{out_dir_spatial}/{spp}_time_brt.tif" ) )
  r_pobs <- rast(
  glue::glue("{out_dir_spatial}/{spp}_pobs_time_brt.tif")
)
  p <- read_rds(glue::glue("{prediction_layer_loc}/{spp}_time_pred_brt.rds"))
  pobs <- 1 - dpois(0, lambda = p$.pred)
min_q <- min(p, na.rm = T)

mf_route <- read_sf(
  g("{mf_route_loc}/MFCAR_Route.shp")
)

mf_r <- crop(
  predicted_raster,
  st_bbox(st_transform(st_buffer(mf_route, 20000), st_crs(predicted_raster)))
)
# quat_loc <- st_crop(quat, st_bbox(st_transform(st_buffer(mf_route, 20000), st_crs(predicted_raster))))

# cairo_pdf(g("{out_dir}/{spp}_time_map_MFCAR.pdf"), fallback_resolution = 300, height = 17, width = 11)
# ggplot( ) + tidyterra::geom_spatraster(data = mf_r) +
#   geom_sf(data = mf_route, alpha = 0.5, linetype =2) +
#   tidyterra::scale_fill_whitebox_c(
#     palette = "muted",#breaks = c(0.05, 0.1, 0.2, 0.3, 0.8),
#     # labels = scales::label_number(suffix = "indiv/ha"),
#     n.breaks = 8,#limits= c(min_q, quantile(p, 0.999, na.rm=T)),
#     guide = guide_legend(reverse = TRUE))  +
#   labs(fill = "individuals/ha", title = spp,
#        colour = "#\nNon-zero\ncounts") +
#   geom_sf(data = quat_loc, fill = NA, colour = 'white')
# dev.off()
blob <- read_sf("E:/SPATIAL/RA_Blob/MiningClaimsRoF_July2025.shp") |>
  st_transform(ont.proj)
# map_plot <- ggplot() +
#   tidyterra::geom_spatraster(data = predicted_raster, maxcell = 2e6) +
#   labs(fill = "individuals/ha", title = spp, colour = "#\nNon-zero\ncounts") +
#   tidyterra::scale_fill_whitebox_c(
#     palette = "muted", #breaks = c(0.05, 0.1, 0.2, 0.3, 0.8),
#     # labels = scales::label_number(suffix = "indiv/ha"),
#     n.breaks = 8,
#     limits = quantile(p, c(0, 0.999), na.rm = T),
#     guide = guide_legend(reverse = TRUE)
#   ) +
#   # geom_sf(data = napken_lake, shape =2 )+
#   # geom_sf(data = locs_in) +
#   geom_sf(data = ra_area, fill = NA, linetype = 2, colour = 'white') +
#   geom_sf(data = mesh_data, fill = NA, linetype = 3, colour = 'black')

can <- st_read("E:/CWS_ONT_local/Base_Data/canada.shp") %>%
  # filter(NAME %in% c("Manitoba", "Québec", )) |>
  st_transform(ont.proj, partial = F, check = T)

inset <- ggplot() +
  geom_sf(data = can, fill = NA) +
  geom_sf(
    data = ra_area_official,
    fill = tidyterra::wiki.colors(1),
    alpha = 0.8
  ) +
  ggthemes::theme_map()

map_plot <-
  ggplot() +
  tidyterra::geom_spatraster(
    data = mask(predicted_raster, ra_area),
    maxcell = 1e6
  ) +
  labs(
    fill = "individuals/ha",
    title = NULL,
    colour = "#\nNon-zero\ncounts"
  ) +
  tidyterra::scale_fill_wiki_c(
    # palette = "muted", #breaks = c(0.05, 0.1, 0.2, 0.3, 0.8),
    # labels = scales::label_number(suffix = "indiv/ha"),
    # n.breaks = 8,
    limits = quantile(p, c(0, 0.999), na.rm = T),
    guide = guide_legend(reverse = TRUE)
  ) +
  geom_sf(data = blob, fill = 'red', alpha = 0.4, colour = NA) +
  ggthemes::theme_map()

# ggsave(
#   "output/CONW_with_inset_rof.jpeg",
#   test_plot +
#     patchwork::inset_element(inset, 0.7, 0.6, 1.2, 1),
#   width = 8,
#   height = 6
# )
# ggsave("output/CONW_no_inset_rof.jpeg", test_plot, width = 8, height = 6)

p_obs_plot <-
  ggplot() +
  tidyterra::geom_spatraster(
    data = mask(r_pobs, ra_area_official),
    maxcell = 2e6
  ) +
  labs(
    fill = "Probability\nof\nobservation",
    title = spp,
    colour = "#\nNon-zero\ncounts"
  ) +
  tidyterra::scale_fill_wiki_c(
    # palette = "muted", #breaks = c(0.05, 0.1, 0.2, 0.3, 0.8),
    # labels = scales::label_number(suffix = "indiv/ha"),
    n.breaks = 8,
    limits = quantile(pobs, c(0, 0.999), na.rm = T),
    guide = guide_legend(reverse = TRUE)
  ) +
  geom_sf(data = BASSr::ontario, fill = NA) +
  # geom_sf(data = napken_lake, shape =2 )+
  # geom_sf(data = locs_in) +
  geom_sf(data = ra_area, fill = NA, linetype = 2, colour = 'white') +
  # geom_sf(data = mesh_data, fill = NA, linetype = 3, colour = 'black') +
  coord_sf(
    xlim = st_bbox(ra_buffer)[c(1, 3)],
    ylim = st_bbox(ra_buffer)[c(2, 4)]
  )

ggsave(
  plot = map_plot,
  g("{out_dir}/{spp}_time_map_{period_to_use}_{doy_used}.jpeg"),
  width = 6.73,
  height = 8.5
)
ggsave(
  plot = p_obs_plot,
  g("{out_dir}/{spp}_time_map_pobs_{period_to_use}_{doy_used}.jpeg"),
  width = 6.73,
  height = 8.5
)
}


spp_completed <- list.files(
  path = bundle_locs,
  "\\w{4}_time",
  recursive = T
) |>
  str_extract("\\w{4}") |>
  unique()

spp_no_map <- list.files(BRT_output_loc, "time_map_D", recursive = T) |> 
  str_subset("Unk", negate=T) |> str_extract("\\w{4}")



walk(spp_completed[!spp_completed %in% spp_no_map],
     just_maps, .progress = T)
