source("R/__globals.R")
full_locs <- locs
locs <- individual_locs

# locations_with_data <- filter(locations, location %in% our_wt_data$location)
# far_north_PSUs_with_data <- locations_with_data |> filter(grepl("^[PQ]\\d+", location))  %>%
#   st_filter(x = far_north_PSUs, y = .,
#             .predicate = st_contains)
prelim_crescent <- read_rds("data/crescent_deployment_2024_prelim.rds") |> 
  st_transform(ont.proj)
locs_split <- full_locs |> group_split(type)


farnorth_boundary <- read_sf(g(
  "{spatial_raster_location}/1_FNLC/FNLC_boundary.shp"
))


exclusion_grid <- st_make_grid(locs, n = c(20, 20)) |>
  st_as_sf() |>
  rename(geometry = x) |>
  st_filter(locs) |>
  st_join(full_locs |> distinct(site_id, type, geometry)) |>
  summarize(n = n(), .by = c(type,geometry)) |>
  arrange(desc(n)) |>
  mutate(ID = row_number())




ggplot(exclusion_grid, aes(fill = n)) + geom_sf() + 
  geom_sf(data = prelim_crescent, fill = NA) + facet_wrap(~type)

avoids <- 
st_filter(exclusion_grid |> filter(type=="ARU recording"), st_buffer(prelim_crescent, 50000))

ex_grd <- filter(exclusion_grid, !ID %in% avoids$ID) |> mutate(p=1/(n))

(exclusion_grid$n == 1) |> sum()
withr::with_seed(56198, {
  test_blocks <- spsurvey::grts(st_centroid(ex_grd), 
                                n_base = c("Point Count" = 10,
                                           "ARU recording" = 10),
                                aux_var = "p", mindis = 40000, stratum_var = "type")

  test_locations <- st_filter(
    locs,
    ex_grd[ex_grd$ID %in% test_blocks$sites_base$ID, ]
  )
}) 
ggplot(ex_grd) + geom_sf() + geom_sf(data = prelim_crescent, fill = NULL) +
  geom_sf(data = test_blocks$sites_base, aes(colour = n)) +
  facet_wrap(~type)
test_locs_ids <- left_join(
  test_locations,
  st_drop_geometry(locs),
  by = join_by(site_id)
)

full_recordings_train <- g("{rds_data_loc}/all_events.rds") |>
  read_rds() |>
  filter(!is.na(site_id) & site_id %in% locs$site_id) |>
  filter(
    collection %in%
      c("WildTrax", "ONATLAS3PC") &
      !site_id %in% test_locs_ids$site_id &
      str_detect(location, "WHA-", negate = T) &
      str_detect(project, "Nocturnal", negate = T) &
      str_detect(project, "Resample", negate = T)
  ) |>
  left_join(
    st_join(locs, select(ontario_ez, on_er = ZONE_NAME)),
    by = join_by(site_id)
  )


withr::with_seed(6894, {
  test_recordings <-
    full_recordings_train |>
    slice_sample(prop = 0.1,
                 by = location)
})


train_recordings <- full_recordings_train |>
  filter(!event_id %in% test_recordings$event_id)


write_rds(
  list(
    exclusion_grid = exclusion_grid,
    test_blocks = test_blocks,
    test_locations = test_locations,
    test_recordings = test_recordings,
    train_recordings = train_recordings
  ),
  g("output/{date_compiled}_test_training_data.rds")
)
