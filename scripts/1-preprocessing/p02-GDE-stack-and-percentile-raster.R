### ---------------------\\ 
# Script objective:
# Create GDE comparison raster stack and an average GDE area percentile raster 
### ---------------------\\
library(here); source(here(("scripts/on_button.R")))
###

###
## Import the rasters at 0.5 degree:
###

# Import Rohde et al. GDE/AA fraction at 0.5 degree
ro24_30m = terra::rast("D:/Geodatabase/GDEs/Rohde_2024/GDE_data_deposit_v6/GDE_30arcmin.tif")
ro24_30m = ro24_30m$GDE_frac_GA/1e8
ro24_30m_AAfrac = terra::rast("D:/Geodatabase/GDEs/Rohde_2024/GDE_data_deposit_v6/GDE_30arcmin.tif")$AA_frac_GA/1e8

# Import the Huggins et al. area density at 5 arcmin and resample to 30 arcmin
hu23_30m = terra::resample(x = terra::rast("D:/Geodatabase/GDEs/Huggins_2023/hu23_areadens_5m.tif"),
                           y = WGS84_areaRaster(0.5) |> rast(),
                           method = "bilinear") 

# Link et al. GDE probability at 0.5 degree
li23_30m = terra::rasterize(x = terra::vect("D:/Geodatabase/GDEs/Link_2023/GDEs_at_risk.shp"),
                            y = WGS84_areaRaster(0.5) |> rast(), 
                            field = "I7", 
                            touches = TRUE) 

# Extend RO24 to full global extent
ro24_30m = terra::resample(x = ro24_30m, y = hu23_30m, method = "near")

# raster stack
R.gde_comp = c(ro24_30m, hu23_30m |> terra::crop(y = ro24_30m), li23_30m |> terra::crop(y = ro24_30m))
# R.gde_comp = c(ro24_30m, hu23_30m, li23_30m)
names(R.gde_comp) = c("RO", "HU", "LI")

writeRaster(R.gde_comp, filename = here("data/GDE_comparison_0d5.tif"), overwrite = T)

R.gde_comp_otherlands = c(hu23_30m, li23_30m)
names(R.gde_comp_otherlands) = c("HU", "LI")

# Create a raster to identify places included in Rohde et al. 2024
mask_RO = rast(R.gde_comp$RO)
mask_RO[] = 0
mask_RO[R.gde_comp$RO >= 0] = 1

# Create mask raster for terrestrial land area (i.e., for HU and LI)
mask_HU = rast(R.gde_comp_otherlands$HU)
mask_HU[] = 0
mask_HU[R.gde_comp_otherlands$HU > 0] = 1

mask_LI = rast(R.gde_comp_otherlands$LI)
mask_LI[] = 0
mask_LI[R.gde_comp_otherlands$LI >= 0] = 1

mask_otherlands = terra::mask(x = mask_LI, mask = terra::extend(x=mask_RO, y= mask_LI, fill = 0), 
                              maskvalues = 1, updatevalue = 0)

## remove drylands from otherlands for comparison
mask_HULI = mask_HU*mask_LI
mask_otherlands = terra::mask(x = mask_HULI, mask = terra::extend(x=mask_RO, y= mask_LI, fill = 0), 
                              maskvalues = 1, updatevalue = 0)

R.gde_comp_drylands = R.gde_comp |> terra::mask(mask = mask_RO, maskvalues = 0) 
R.gde_comp_otherlands = R.gde_comp_otherlands |> terra::mask(mask = mask_LI, maskvalues = 0) 

###########################################################################
#### Look at average percentile per pixel
###########################################################################
ptile_stack = c(reclassify_percentile_bins(raster_layer = R.gde_comp$RO, 
                                           area_grid = WGS84_areaRaster(0.5) |> rast() |> 
                                             crop(R.gde_comp_drylands$RO), 
                                           num_bins = 100),
                reclassify_percentile_bins(raster_layer = R.gde_comp$LI, 
                                           area_grid = WGS84_areaRaster(0.5) |> rast() |> 
                                             crop(R.gde_comp_otherlands$LI), 
                                           num_bins = 100),
                reclassify_percentile_bins(raster_layer = terra::mask(R.gde_comp$HU, mask = mask_HU, maskvalues = 0), 
                                           area_grid = WGS84_areaRaster(0.5) |> rast() |> 
                                             crop(R.gde_comp_otherlands$HU), 
                                           num_bins = 100))
names(ptile_stack) = c("RO", "LI", "HU")

ptile_mean = mean(ptile_stack, na.rm = T)
writeRaster(ptile_mean, here("data/GDE_Ptile_mean.tif"), overwrite = T)