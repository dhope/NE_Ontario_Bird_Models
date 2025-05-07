source("R/__globals.R")
library(QPAD)


QPAD_globals <- 
expand_grid(species =QPAD::getBAMspecieslist(), r = c(0.5, 1, Inf), t = dur) |> # dur from NAPOPS_QPAD
  # nest_by(species) |> 
  rowwise() |> 
  mutate(QPAD = list(QPAD::globalBAMcorrections(species, r, t))) |> 
  unnest(QPAD) |> 
  mutate(o=log(A) + log(p) + log(q))

write_rds(QPAD_globals, "output/QPAD_global_offsets.rds")
na <- read_rds("output/na_pops_offsets.rds")

comp <- 
full_join(na, QPAD_globals, by = 
            join_by(
              spp == species, 
              max_dist == r, 
              time_minutes == t
            ))


ggplot(comp, aes(o.x, o.y, colour = factor(max_dist))) + geom_point() +
# ggplot(comp, aes(o.x, o.y, colour = time_minutes)) + geom_point() +
  geom_abline(intercept = 0, slope =1)





