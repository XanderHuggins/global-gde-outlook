### ---------------------\\ 
# Script objective:
# Plot areas of percentile agreement or disagreement between rasters 
### ---------------------\\
library(here); source(here(("scripts/on_button.R")))
###

ptile_mean = terra::rast(here("data/GDE_Ptile_mean.tif"))

########### Import necessary additions for map

# 1 Caspain sea
casp = terra::vect("D:/Geodatabase/Admin-ocean-boundaries/worldglwd1.shp")
casp_r = terra::rasterize(x = casp, y = ptile_mean, 1, touches = F)
ptile_mean[casp_r == 1] = NA

# 2 regional call-outs
areas_red = terra::vect("D:/Geodatabase/Admin-ocean-boundaries/boundaries/GDE_panels/study_areas/study_areas.shp")
areas_red = terra::centroids(areas_red)

outline = terra::vect("D:/D_documents/1.projects-scripts/sustainability-puzzles/data/land_mask_polygon.sqlite") |> 
  st_as_sf()

outline_r = terra::vect(outline)
outline_r = terra::rasterize(outline_r, ptile_mean, 1, touches = T)
ptile_mean[is.na(outline_r)] = NA

# Map of P-MEAN of GDE area density globally
map =  
  tm_shape(outline, crs = "+proj=robin") +
  tm_fill(fill ="grey", lwd = 0, border.col = NA) +
  tm_shape(terra::project(x = ptile_mean, y = "+proj=robin", method = "near")) +
  # tm_raster(palette = met.brewer("Isfahan1", n = 100),
  tm_raster(col.scale = tm_scale_continuous(
    values = met.brewer("Paquin", n = 100)[51:100],
    ticks = seq(0, 100)
  )) +
  tm_shape(areas_red) +
  tm_borders(col = "red", lwd = 5) +
  # tm_shape(ramsar.pts) +
  # tm_dots(size = 0.6, shape = 16, fill = "red") + 
  tm_layout(legend.show = F, legend.frame = F, frame = F, bg = FALSE, bg.color = "transparent")
map

tmap_save(map, here("plots/gMAP_gde_AreaDensity_Ptile_MEAN_monotonicscale.png"), dpi = 400, units = "in")

########################
## SI FIGURES BELOW
########################

###########################################################################
#### CREATE BASE MAPS OF GDE DENSITY FOR EACH STUDY
###########################################################################

map =  
  tm_shape(outline, crs = "+proj=robin") +
  tm_fill(col = "grey", lwd = 0, border.col = NA) +
  tm_shape(terra::project(x = R.gde_comp$RO, y = "+proj=robin", method = "near")) +
  tm_raster(palette = met.brewer("Paquin", n = 100)[51:100],
            breaks = seq(0,0.33, length.out = 50)) +
  tm_shape(outline) + 
  tm_borders(col = "black", lwd = 1) +
  tm_layout(legend.show = F, legend.frame = F, frame = F, bg = FALSE, bg.color = "transparent")
map
tmap_save(map, here("plots/gMAP_Rohde_basemap_0_to_0d33.pdf"), dpi = 400, units = "in")

map =  
  tm_shape(outline, crs = "+proj=robin") +
  tm_fill(col = "grey", lwd = 0, border.col = NA) +
  tm_shape(terra::project(x = R.gde_comp$LI, y = "+proj=robin", method = "near")) +
  tm_raster(palette = met.brewer("Paquin", n = 100)[51:100],
            breaks = seq(0,0.33, length.out = 50)) +
  tm_shape(outline) + 
  tm_borders(col = "black", lwd = 1) +
  tm_layout(legend.show = F, legend.frame = F, frame = F, bg = FALSE, bg.color = "transparent")
map
tmap_save(map, here("plots/gMAP_Link_basemap_0_to_0d66.pdf"), dpi = 400, units = "in")

map =  
  tm_shape(outline, crs = "+proj=robin") +
  tm_fill(col = "grey", lwd = 0, border.col = NA) +
  tm_shape(terra::project(x = terra::mask(R.gde_comp$HU, mask = mask_HU, maskvalues = 0), 
                          y = "+proj=robin", method = "near")) +
  tm_raster(palette = met.brewer("Paquin", n = 100)[51:100],
            breaks = seq(0,1, length.out = 50)) +
  tm_shape(outline) + 
  tm_borders(col = "black", lwd = 1) +
  tm_layout(legend.show = F, legend.frame = F, frame = F, bg = FALSE, bg.color = "transparent")
map
tmap_save(map, here("plots/gMAP_Huggins_basemap_0_to_1.pdf"), dpi = 400, units = "in")

###########################################################################
#### CREATE PERCENTILE MAPS OF GDE DENSITY FOR EACH STUDY
###########################################################################
map =  
  tm_shape(outline, crs = "+proj=robin") +
  tm_fill(col ="grey", lwd = 0, border.col = NA) +
  tm_shape(terra::project(x = reclassify_percentile_bins(raster_layer = R.gde_comp$RO, 
                                                         area_grid = WGS84_areaRaster(0.5) |> rast() |> 
                                                           crop(R.gde_comp_drylands$RO), 
                                                         num_bins = 100),
                          y = "+proj=robin", method = "near")) +
  tm_raster(palette = met.brewer("Paquin", n = 100)[51:100],
            breaks = seq(0,100, length.out = 50)) +
  tm_shape(outline) + 
  tm_borders(col = "black", lwd = 1) +
  tm_layout(legend.show = F, legend.frame = F, frame = F, bg = FALSE, bg.color = "transparent")
map
tmap_save(map, here("plots/gMAP_Rohde_percentiles.pdf"), dpi = 400, units = "in")

map =  
  tm_shape(outline, crs = "+proj=robin") +
  tm_fill(col ="grey", lwd = 0, border.col = NA) +
  tm_shape(terra::project(x = reclassify_percentile_bins(raster_layer = R.gde_comp$LI, 
                                                         area_grid = WGS84_areaRaster(0.5) |> rast() |> 
                                                           crop(R.gde_comp_otherlands$LI), 
                                                         num_bins = 100),
                          y = "+proj=robin", method = "near")) +
  tm_raster(palette = met.brewer("Paquin", n = 100)[51:100],
            breaks = seq(0,100, length.out = 50)) +
  tm_shape(outline) + 
  tm_borders(col = "black", lwd = 1) +
  tm_layout(legend.show = F, legend.frame = F, frame = F, bg = FALSE, bg.color = "transparent")
map
tmap_save(map, here("plots/gMAP_Link_percentiles.pdf"), dpi = 400, units = "in")

map =  
  tm_shape(outline, crs = "+proj=robin") +
  tm_fill(col ="grey", lwd = 0, border.col = NA) +
  tm_shape(terra::project(x = reclassify_percentile_bins(raster_layer = terra::mask(R.gde_comp$HU, mask = mask_HU, maskvalues = 0), 
                                                         area_grid = WGS84_areaRaster(0.5) |> rast() |> 
                                                           crop(R.gde_comp_otherlands$HU), 
                                                         num_bins = 100),
                          y = "+proj=robin", method = "near")) +
  tm_raster(palette = met.brewer("Paquin", n = 100)[51:100],
            breaks = seq(0,100, length.out = 50)) +
  tm_shape(outline) + 
  tm_borders(col = "black", lwd = 1) +
  tm_layout(legend.show = F, legend.frame = F, frame = F, bg = FALSE, bg.color = "transparent")
map
tmap_save(map, here("plots/gMAP_Huggins_percentiles.pdf"), dpi = 400, units = "in")
