count_na <- function(df) {
  dplyr::summarize(df, dplyr::across(dplyr::everything(), ~ sum(is.na(.x))))
}
count_na(spatial_cov) |> glimpse()


ggplot(spatial_cov, aes(NFIS_age_100)) +
  geom_density() +
  geom_density(
    data = replace_na(
      spatial_cov,
      list = NFIS_age_100 = median(spatial_cov$NFIS_age_100, na.rm = T)
    ),
    fill = 'red',
    alpha = 0.2
  )
