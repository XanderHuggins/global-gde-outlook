### ---------------------\\ 
# Script objective:
# Create a global GDE-terr and GDE-aqua density raster using:
# Huggins et al. 2023 GDE type classification
# and Link et al. 2023 GDE land use suitability
### ---------------------\\
library(here); source(here(("scripts/on_button.R")))

Hu23_classes = terra::rast("D:/Geodatabase/GDEs/Huggins_2023/gde-map.tif")

# Code for GDE classification follows this patter:
# 1xx: terrestrial GDE
# x1x: lentic aquatic GDE
# xx1: lotic aquatic GDE

sub_dict = tibble(
  from = c(1, 10, 11, 100, 101, 110, 111),
  terr = c(0,  0,  0,   1,   1,   1,   1),
  lent = c(0,  1,  1,   0,   0,   1,   1),
  lotc = c(1,  0,  1,   0,   1,   0,   1))

Hu23_Terr = terra::subst(Hu23_classes, from = sub_dict$from, to = sub_dict$terr, others = 0)
Hu23_Lent = terra::subst(Hu23_classes, from = sub_dict$from, to = sub_dict$lent, others = 0)
Hu23_Lotc = terra::subst(Hu23_classes, from = sub_dict$from, to = sub_dict$lotc, others = 0)

Hu23_Aqua = max(Hu23_Lent, Hu23_Lotc)

# area_1km = WGS84_areaRaster(0.5/60) |> rast() # a few minutesto run... 
area_1km = terra::rast("D:/! Project Archive/projects/archetypes/data/wgs-area-ras-30-arcsec.tif")

Hu23_Terr_area = Hu23_Terr * area_1km
Hu23_Lent_area = Hu23_Lent * area_1km
Hu23_Lotc_area = Hu23_Lotc * area_1km
Hu23_Aqua_area = Hu23_Aqua * area_1km

# aggregate to 30 arcmin
grid_30m = WGS84_areaRaster(30/60) |> rast()

# multiply by Link et al. 2023 GDE land cover 
li23 = terra::rasterize(x = terra::vect("D:/Geodatabase/GDEs/Link_2023/GDE_potentials.shp"),
                        y = WGS84_areaRaster(0.5) |> rast(), # this is the base resolution of Link 2023
                        field = "I5", # share of relevant land covers
                        touches = TRUE)

li23_area_30m = li23 * grid_30m

Hu23_Terr_area_30m = (terra::aggregate(x = Hu23_Terr_area, fact = 60, fun = "sum") ) * li23 /1e6
Hu23_Lent_area_30m = (terra::aggregate(x = Hu23_Lent_area, fact = 60, fun = "sum") ) * li23 /1e6 
Hu23_Lotc_area_30m = (terra::aggregate(x = Hu23_Lotc_area, fact = 60, fun = "sum") ) * li23 /1e6
Hu23_Aqua_area_30m = (terra::aggregate(x = Hu23_Aqua_area, fact = 60, fun = "sum") ) * li23 /1e6

names(Hu23_Terr_area_30m) = "terr_GDE_area"
names(Hu23_Aqua_area_30m) = "aqua_GDE_area"

R.gde_types_hu = c(Hu23_Terr_area_30m, Hu23_Aqua_area_30m)
writeRaster(R.gde_types_hu, filename = here("data/GDE_types_hu_for_comparison.tif"), overwrite = T)