source("R/__paths.R")
rs <- spatial_raster_location |> list.dirs()
rasts <- map(
  rs,
  ~ {
    list.files(.x, pattern = ".adf", full.names = F, recursive = F) |>
      length()
  }
) |>
  list_c()


rs[rasts > 0]

cairo_pdf("Covariates_zeros.pdf")
walk(
  rs[rasts > 0],
  ~ {
    r <- try({
      terra::rast(.x)
    })
    if (class(r) != "try-error") {
      plot(
        ggplot() +
          tidyterra::geom_spatraster(data = r == 0) +
          labs(
            title = str_remove(
              .x,
              "D:/DELIVERABLES/Spatialworks24_25/DataLayers/"
            ),
            fill = "Is zero?"
          )
      )
    }
  },
  .progress = T
)
dev.off()
