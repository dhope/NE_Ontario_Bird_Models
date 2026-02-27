source("R/06_BRT_time.R")



gen_count_obs <- function(spp){
p <- prep_brt_data(spp)
# ggplot(p$df_std,aes( y)) + geom_histogram(binwidth = 1)

gt_rec <- 
janitor::tabyl(p$df_std, y, type) |> 
  mutate(across(2:3, ~round(.x/sum(.x), 2))) |> 
  gt::gt()|> 
  gt::tab_header(glue::glue("Counts by site - {spp}"))


t2_plot <- 
p$df_std |> ggplot(aes(ymd("2025-01-01") +doy,
                       t2se, colour = y>0 )) +
  geom_point(alpha = 0.2) +
  facet_wrap(~event_gr) +
  rcartocolor::scale_colour_carto_d(palette = 2) +
  labs(x = NULL, y = "Time to sun event",
       colour = "Obs")


df_train_site <- p$df_std |> 
  # filter(y>0) |> 
  slice_max(y, by = c(location), with_ties = F) |> 
  st_as_sf(coords = c("X", "Y"), crs = ont.proj)

gt_site <- 
janitor::tabyl(st_drop_geometry(df_train_site), y, type) |> 
  mutate(across(2:3, ~round(.x/sum(.x), 2))) |> 
  gt::gt() |> 
  gt::tab_header(glue::glue("Max count by site - {spp}"))


gg_obs <- 
  ggplot() + 
  geom_sf(data = st_intersection(ra_area_official, BASSr::ontario)) +
  geom_sf(data = df_train_site[df_train_site$y==0,], colour = 'grey', alpha = 0.5) +
  geom_sf(data =df_train_site[df_train_site$y>0,],
          aes(colour = type)) +
  ggtitle(spp) + ggthemes::theme_map() +
  rcartocolor::scale_colour_carto_d()
  
  
  list(
  p1 = (gg_obs + (patchwork::wrap_table(gt_rec) / patchwork::wrap_table(gt_site))),
  p2 =  t2_plot + labs(title = spp) )
  
  
}


spp_comp <-
  list.files(
    BRT_output_loc,
    "time_model_pred.rds",
    recursive = T,
    full.names = F
  ) |>
  str_extract("\\w{4}") |>
  unique() |>
  sort() 

sgco <- safely(gen_count_obs)
all_spp <- map(spp_comp, sgco, .progress = T)
alspt <- all_spp |> purrr::compact() |> transpose() |> pluck("result") |> compact() |> 
  transpose()
cairo_pdf("species_counts.pdf", width = 8, height = 6)
alspt$p1
dev.off()


cairo_pdf("time_date_obs.pdf", width = 12, height = 6)
alspt$p2
dev.off()


brd_fix<- gen_count_obs("BLJA")
