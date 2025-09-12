source("R/__globals.R")
locs <- individual_locs

# locations_with_data <- filter(locations, location %in% our_wt_data$location)
# far_north_PSUs_with_data <- locations_with_data |> filter(grepl("^[PQ]\\d+", location))  %>%
#   st_filter(x = far_north_PSUs, y = .,
#             .predicate = st_contains)

farnorth_boundary <- read_sf(g(
  "{spatial_raster_location}/1_FNLC/FNLC_boundary.shp"
))


exclusion_grid <- st_make_grid(locs, n = c(20, 20)) |>
  st_as_sf() |>
  rename(geometry = x) |>
  st_filter(locs) |>
  st_join(locs) |>
  summarize(n = n(), .by = geometry) |>
  arrange(desc(n)) |>
  mutate(ID = row_number())

ggplot(exclusion_grid, aes(fill = n)) + geom_sf()
(exclusion_grid$n == 1) |> sum()
withr::with_seed(12354, {
  test_blocks <- spsurvey::grts(exclusion_grid, n_base = 30)

  test_locations <- st_filter(
    locs,
    exclusion_grid[exclusion_grid$ID %in% test_blocks$sites_base$ID, ]
  )
})

test_locs_ids <- left_join(
  test_locations,
  st_drop_geometry(aggregated_locs),
  by = join_by(site_id_agg)
)

full_recordings_train <- g("{rds_data_loc}/all_events.rds") |>
  read_rds() |>
  filter(location %in% aggregated_locs$location) |>
  filter(
    collection %in%
      c("WildTrax", "ONATLAS3PC") &
      !location %in% test_locs_ids$location &
      str_detect(location, "WHA-", negate = T) &
      str_detect(project, "Nocturnal", negate = T)
  ) |>
  left_join(
    st_join(aggregated_locs, select(ontario_ez, on_er = ZONE_NAME)),
    by = join_by(project, location, collection)
  )


withr::with_seed(6894, {
  test_recordings <-
    full_recordings_train |>
    slice_sample(prop = 0.2, by = location)
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
