source("R/__globals.R")
locs <- individual_locs
locs_atlas2 <- 
  g("{rds_data_loc}/locations.rds") |> 
  read_rds() |> 
  filter(collection %in% c( "OBBA2PC") &
           str_detect(project, "Nocturnal|Resample", negate=T)) |> 
  st_filter(ra_buffer) 
if(run_a2) locs <- locs_atlas2



ld <- list.dirs(
  spatial_raster_location, recursive = F
)


library(terra)
## List of all layers
ll <- map(ld, list.dirs, recursive=F) |> 
  list_c()

length(ll[str_detect(ll, "info$", negate = T)])








# OLCC --------------------------------------------------------------------


olcc <- terra::rast(ll[ str_detect(ll, "olc", negate = F)&str_detect(ll, "_200|_1km", negate = T)])
olcc_extract <- function(olcc, dist_){
  exactextractr::exact_extract(olcc,
                             st_transform(
                               st_buffer(locs, dist = dist_),
                               st_crs(olcc)), 'sum') |>  
 rowwise() |> 
  mutate(total = sum(c_across(starts_with('sum.')))) |> 
  ungroup() |> 
  mutate(across(starts_with('sum.'), ~{.x/total*100})) |> 
  rename_with(.fn = ~str_remove(.x, "sum."), starts_with('sum.')) |> 
  rename_with(.fn = ~glue::glue("{str_remove(.x, 'cls')}_{dist_}"),
              everything())
}
 
olcc_500 <- olcc_extract(olcc , 500)
olcc_100 <- olcc_extract(olcc , 100)



cat_layers <- terra::rast(
  ll[str_detect(ll, "orig", negate = F) & str_detect(ll, "lc|geo|wat_", negate = F) & str_detect(ll, "asp8|geob|geoq", negate = T)] )



cat_layers_extract <- function(layers, dist_){
  map(names(layers), 
                          ~{exactextractr::exact_extract(layers[[.x]], 
                                          st_transform(
                                            st_buffer(locs, dist = dist_),
                                            st_crs(layers)), 'frac', force_df =T, 
                                          default_value = -69) |> 
                              mutate(across(contains("frac_"), ~{.x*100})) |> 
  rename_with(str_replace, pattern ="frac_(\\d+)",
              replacement = glue::glue("{.x}_\\1_{dist_}" ))  |> as_tibble() 
                          } ) |> 
  list_cbind() |> 
  rename_with(str_remove,pattern = "_orig30m")
}

cat_layers_100 <- cat_layers_extract(cat_layers, 100)
cat_layers_500 <- cat_layers_extract(cat_layers, 500)




mod_layers <- str_subset(ll, "asp8|geob|geoq|/sp_") |> 
  str_subset("1km|200", negate = T) |> rast()

mod_layers_extract <- function(layers, dist_, f = 'mode'){
  map(names(layers), 
      ~{
        n <- glue::glue("{.x}_{dist_}")
        tibble({{n}}  :=
        factor(exactextractr::exact_extract(layers[[.x]], 
                                     st_transform(
                                       st_buffer(locs, dist = dist_),
                                       st_crs(layers)), f) )) }) |> 
    list_cbind() |> 
    rename_with(str_remove,pattern = "_orig30m") |> 
    rename_with(~str_replace(.x, "sp_", "NFIS2_sp_"))
}
  

mod_ex_500 <- mod_layers_extract(mod_layers, 500) 
mod_ex_100 <- mod_layers_extract(mod_layers, 100)

fire_500 <- 
  mod_layers_extract(str_subset(ll, "fire|insct|harv") |> 
                     str_subset("1km|200", negate = T) |> rast(), 500, 'max')
fire_100 <- mod_layers_extract(str_subset(ll, "fire|insct|harv") |> 
                                 str_subset("1km|200", negate = T) |> rast(), 100,
                               'max')


fire_100$fire_rec30m_100 |> janitor::tabyl()
fire_100$insct_rec30m_100 |> janitor::tabyl()
fire_100$harv_rec30m_100 |> janitor::tabyl()


non_lc <- ll[(str_detect(ll, "orig", negate = F) & 
                str_detect(ll, "lc|geo|/sp_|wat_|fire|insct|harv", negate = T))|
               (str_detect(ll, "lc|info|orig|1km|_200|asp8|wat_|fire|insct|harv", negate = T) )
             ] 

get_non_lc <- function(r, dist_){
  map(r, 
                ~{
                  r <- rast(.x)
                  print(r)
                  if(str_detect(.x, "dem")){
                    n <- str_extract(.x,"(dem_\\w+)" )
                  } else{
                    n <- str_extract(.x, "((?<=DataLayers/\\d{1,2}_)[\\w,\\W]+/\\w+(?=_(orig|rec)30m))") |> 
                      str_replace("/", "_") }
                  n <- paste0(c(n, dist_), collapse = "_")
                  exactextractr::exact_extract(r, 
                                               st_transform(
                                                 st_buffer(locs, dist = dist_),
                                                 st_crs(r)), 'mean') %>% 
                    tibble::tibble(
                                   {{n}}:=.)  |> 
                    dplyr::select(2)
                    }) |> bind_cols()
}
   
non_lc_r <- get_non_lc(non_lc, 100)
non_lc_r_500 <- get_non_lc(non_lc, 500)



# To fix ---------
# dem_pslp # perspective (direction)




cov <- bind_cols(cat_layers_100, olcc_100, non_lc_r,
                     cat_layers_500, olcc_500, non_lc_r_500,
                 mod_ex_100, mod_ex_500, fire_100, fire_500) |> 
  dplyr::select(-starts_with('frac_'))
  

if(run_a2){
  write_rds(cov, g("output/rds/{Sys.Date()}_spatial_covariates_Atlas2.rds"))
}else{
write_rds(cov, g("output/rds/{Sys.Date()}_spatial_covariates_data.rds"))
}





