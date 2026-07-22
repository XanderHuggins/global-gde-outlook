### ---------------------\\ 
# Script objective:
# Plot areas of percentile agreement or disagreement between rasters 
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

ptile_mean = terra::rast(here("data/GDE_Ptile_mean.tif"))

########### Import necessary additions for map

# 1 Caspain sea
casp = terra::vect("D:/Geodatabase/Admin-ocean-boundaries/worldglwd1.shp")
casp_r = terra::rasterize(x = casp, y = ptile_mean, 1, touches = F)
ptile_mean[casp_r == 1] = NA

# 2 regional call-outs
areas_red = terra::vect("D:/Geodatabase/Admin-ocean-boundaries/boundaries/GDE_panels/study_areas/study_areas.shp")
areas_red = terra::centroids(areas_red)

outline = terra::vect("C:/Users/xande/Documents/1.projects-scripts/sustainability-puzzles/data/land_mask_polygon.sqlite") |> 
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
outline = terra::vect("C:/Users/xande/Documents/1.projects-scripts/sustainability-puzzles/data/land_mask_polygon.sqlite") |> 
  st_as_sf()

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

###########################################################################
#### Plot relationship between area density and model variance
###########################################################################
n_nonNA  = sum(!is.na(ptile_stack)) 

id_r = rast(ptile_stack)[[1]]
names(id_r) = "id"
id_r[] = seq(1, ncell(id_r))

ptile_square = c(id_r, mean(ptile_stack, na.rm = T), range(ptile_stack, na.rm = T)[[2]] - range(ptile_stack, na.rm = T)[[1]], n_nonNA)
names(ptile_square) = c("id", "mean", "range", "count")
ptile_square[ptile_square$range == 0] = NA
ptile_square[ptile_square$count < 2] = NA
ptile_square$cov = ptile_square$range / ptile_square$mean
ptile_square$norm_range = ptile_square$range / min(ptile_square$mean, 100 - ptile_square$mean)

df_cv = tibble(id = ptile_square$id[], 
               range = ptile_square$range[],
               mean = ptile_square$mean[],
               count = ptile_square$count[]) |> 
  drop_na() |> 
  set_colnames(c("id", "range", "mean", "count")) |>
  dplyr::filter(count >= 2) 

df_cv$norm_range = rep(NA)
df_cv$norm_range = df_cv$range / pmin(df_cv$mean, 100 - df_cv$mean)
df_cv$norm_range_cv = df_cv$range / df_cv$mean

library(mgcv)

m = gam(range ~ s(mean), data=df_cv, family=gaussian())
# df_cv$resid = residuals(m) |> as.numeric() |> unlist()
# df_cv$resid = residuals(m) |> as.numeric() |> unlist()
df_cv$exp_range = predict(m) |> as.numeric() |> unlist()
df_cv$resid = df_cv$range - df_cv$exp_range
df_cv$rel_range = df_cv$range / df_cv$exp_range

df_cv |> ggplot(aes(x = mean, y= exp_range)) +
  geom_point(col = "black", alpha = 0.05)

df_cv$resid_scale = df_cv$resid / sd(df_cv$resid)
df_cv$resid_scale_p = df_cv$resid_scale
df_cv$resid_scale_p[df_cv$resid_scale_p <= 0] = NA
df_cv$resid_scale_p[df_cv$resid_scale_p <= 1.5] = 1
df_cv$resid_scale_p[df_cv$resid_scale_p >  1.5] = 2


r_plot = terra::subst(id_r, 
                      from = df_cv$id |> as.numeric(), 
                      to = df_cv$resid_scale_p |> as.numeric(), others = NA)
plot(r_plot)

r_plot[is.na(r_plot)] = 0
r_plot = r_plot + 1 


R.density_bins = reclassify_percentile_bins(raster_layer = ptile_mean, 
                                            area_grid = WGS84_areaRaster(0.5) |> rast(), 
                                            num_bins = 3)


R.bivar_bins = R.density_bins*10 + r_plot

writeRaster(R.bivar_bins, here("data/rast_uncert_bivar.tif"), overwrite = T)

# pals::brewer.seqseq1(n=4)
bivir_pal = c(
  "#BEE3F6", "#DEC7AA", "#FFAC5E",  #11 (low dens, cert) to 13 (low dens, uncert)
  "#43ABE2", "#9E93AB", "#FA7B75",  #21 to 23
  "#0888E9", "#5F487B", "#F70009" ) #31  (high dens, cert) to 33 (high dens, uncert)


# pals::brewer.seqseq1(n=4)
bivir_pal = c(
  "#594015", "#B68E52", "#D6B37A",  #11 (low dens, cert) to 13 (low dens, uncert)
  "#54B6C7", "#9E93AB", "#FA7B75",  #21 to 23
  "#115755", "#5F487B", "#F70009" ) #31  (high dens, cert) to 33 (high dens, uncert)



R.bivar_bins_DRY = R.bivar_bins # for future use

map =  
  tm_shape(outline, crs = "+proj=robin") +
  tm_fill(fill ="grey", lwd = 0, border.col = NA) +
  tm_shape(terra::project(x = R.bivar_bins, y = "+proj=robin", method = "near")) +
  tm_raster(palette = bivir_pal,
            breaks = c(
              seq(10.5, 13.5, by = 1),
              seq(21.5, 23.5, by = 1),
              seq(31.5, 33.5, by = 1))) +
  tm_shape(outline) + 
  tm_borders(col ="black", lwd = 1) +
  tm_layout(legend.show = F, legend.frame = F, frame = F)
map
tmap_save(map, here("plots/gMAP_GDEdensity_X_resid_above_exp_range.pdf"), dpi = 400, units = "in")


# Isolate areas with elevated uncertainty 
r_plot[r_plot != 3] = NA

highdesn_uncert = R.bivar_bins_DRY
highdesn_uncert[highdesn_uncert != 33] = NA

u_vect = as.polygons(highdesn_uncert, extent=FALSE)
# u_vect = as.polygons(r_plot, extent=FALSE)

map =  
  tm_shape(outline, crs = "+proj=robin") +
  tm_fill(fill ="grey", lwd = 0, border.col = NA) +
  tm_shape(terra::project(x = ptile_mean, y = "+proj=robin", method = "near")) +
  tm_raster(palette = palette = met.brewer("Paquin", n = 100)[51:100],
            breaks = seq(0,100)) +
  tm_shape(outline) + 
  tm_borders(col = "black", lwd = 1) +
  tm_shape(u_vect) + 
  tm_borders(col = "blue", lwd = 1) +
  tm_layout(legend.show = F, legend.frame = F, frame = F, bg = FALSE, bg.color = "transparent")
map
tmap_save(map, here("plots/gMAP_GDEdensity_outline_uncert.pdf"), dpi = 400, units = "in")