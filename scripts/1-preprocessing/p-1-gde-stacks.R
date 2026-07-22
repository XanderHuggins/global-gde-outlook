### ---------------------\\ 
# Script objective:
# Create a raster stack of available global GDE maps at 0.5 degree resolution (coarsest resolution among available datasets)
### ---------------------\\
library(here); source(here(("scripts/on_button.R")))
###

#ro24 = rohde 2024
#hu23 = huggins 2023
#li23 = link 2023

# li23 needs conversion from vector to raster
li23 = terra::rasterize(x = terra::vect("D:/Geodatabase/GDEs/Link_2023/GDEs_at_risk.shp"),
                        y = WGS84_areaRaster(0.5) |> rast(), # this is the base resolution of Link 2023
                        field = "I7", # this is the GDE probability layer
                        touches = TRUE)

hu23 = terra::rast("D:/Geodatabase/GDEs/Huggins_2023/gde-map.tif")

### ----------------------------------------- \\
# create raster stack at 30 arcminute resolution
### ----------------------------------------- \\
ro24 = terra::rast("D:/Geodatabase/GDEs/Rohde_2024/D_Data_GDE_AggregatedLayers/GDE_30arcmin.tif")$GDE_frac_GA

# convert hu23 to area density
  hu23[hu23 > 1] = 1 # create binary representation of hu23
  hu23[is.na(hu23)] = 0
  
  hu23_area = hu23 * (WGS84_areaRaster(0.5/60) |> rast())
  
  # for 5m
  hu23_area_5m = terra::aggregate(x = hu23_area, 
                                  fact = 10,
                                  fun = "sum",
                                  na.rm = T)
  
  hu23_areadens_5m = hu23_area_5m / (WGS84_areaRaster(5/60) |> rast())
  
  writeRaster(x = hu23_areadens_5m, 
              filename = "D:/Geodatabase/GDEs/Huggins_2023/hu23_areadens_5m.tif")
  
  # for 30m
  hu23_area_30m = terra::aggregate(x = hu23_area_5m, 
                                   fact = 6,
                                   fun = "sum",
                                   na.rm = T)
  
  hu23_areadens_30m = hu23_area_30m / (WGS84_areaRaster(30/60) |> rast())
  
  writeRaster(x = hu23_areadens_30m, 
              filename = "D:/Geodatabase/GDEs/Huggins_2023/hu23_areadens_30m.tif")
# end of hu23 conversion

gde_stack_5m = c(li23 |> terra::crop(y= ro24, snap = "near"), 
                 hu23_areadens_30m |> terra::crop(y= ro24, snap = "near"), 
                 ro24) 