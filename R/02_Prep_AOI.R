source("R/__globals.R")
gen_grids <- TRUE




library(fmesher)

inner_mesh <- fmesher::fm_as_sfc(mesh_inla, format='int') |> st_as_sf() |> 
  sf::st_cast("POLYGON") #|> slice(1) 
outer_mesh <- st_zm( fmesher::fm_as_sfc(mesh_inla, format='bnd') |> st_as_sf() |> 
                  sf::st_cast("POLYGON") )
aoi <- inner_mesh

aoi_bb <- st_bbox(aoi) |> st_as_sfc()

ggplot() +
  geom_sf(data = ra_buffer) + geom_sf(data = inner_mesh, colour = 'red', fill = NA)+
  geom_sf(data = outer_mesh, colour = 'blue', fill = NA) +
  geom_sf(data = locs, alpha = 0.5) +
  geom_sf(data = ra_area, fill = NA,linetype =2) +
  geom_sf(data = aoi_bb, fill = NA)


n_ra <- nrow(st_filter(locs, ra_area)) 
n_tota <- nrow(locs)
(n_tota-n_ra)/n_tota






if(gen_grids){

grid_fine <-st_as_stars(.x=aoi_bb, dx = 200)# |> st_crop(aoi)
grid_coarse <-st_as_stars(aoi_bb, dx = 1000) #|> st_crop(aoi)
write_rds(list(fine = grid_fine, coarse = grid_coarse), "output/rds/grids.rds")



write_rds(aoi, g("output/rds/{date_compiled}_Area_of_focus_mesh.rds"))
}
