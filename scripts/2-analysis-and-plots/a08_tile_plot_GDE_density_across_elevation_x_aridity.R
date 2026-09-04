### ---------------------// 
# Script objective:
# Evaluate global relationships of:
# aridity, precipitation seasonality, and elevation
# (1) globally and
# (2) within individual dominant landforms and aquifer lithologies
# and only in regions with low human modification gradient (to minimise effects of land use change)
### ---------------------//
library(here); source(here(("scripts/on_button.R")))

# import GDE area fraction rasters
R.gdes = terra::rast(here("data/GDE_types_hu_for_comparison.tif"))
R.gdes$terrFrac = R.gdes$terr_GDE_area / WGS84_areaRaster(0.5) |> rast()
R.gdes$aquaFrac = R.gdes$aqua_GDE_area / WGS84_areaRaster(0.5) |> rast()
R.gdes$allFrac = terra::rast("D:/Geodatabase/GDEs/Huggins_2023/hu23_areadens_30m.tif")
R.gdes$allFrac[is.na(R.gdes)] = NA

# import aridity raster and resample to 30arcmin
aridity_30m = terra::rast("D:/Geodatabase/Climate/Aridityaridity_0d5.tif")

# import topo
topo = terra::rast("D:/Geodatabase/Terrain/topo_5arcmin_gmt_masked_continents.tif")
topo_30m = terra::resample(x = topo, y = R.gdes, method = "bilinear")

# import precipitation seasonality
p_cv_30m = terra::rast("D:/Geodatabase/Climate/WorldClim/precip_cv_0d5.tif")

# human modification index
modgrad = terra::rast("D:/Geodatabase/Land-use/ModGrad_0d5.tif")
# plot(modgrad_low)

# latitude of grid cells
lat_r = rast(R.gdes, nlyr = 1)                 # or rast(R.gdes, nlyr = 1)
lat_r[] = crds(lat_r, df = TRUE)$y

# import the landforms 
R.landform = terra::rast("D:/Geodatabase/Landforms/Global_Landforms_0d5_mode.tif")
R.landform = terra::resample(x = R.landform, y = lat_r, method = "near")

##########
##
## stack of data to plot/analyse
##
##########
stack_ras = c(R.gdes$terr_GDE_area, R.gdes$aqua_GDE_area, 
             R.gdes$terrFrac, R.gdes$aquaFrac, R.gdes$allFrac, 
             aridity_30m, topo_30m, p_cv_30m, lat_r, modgrad, R.landform, 
             WGS84_areaRaster(0.5) |> rast())
names(stack_ras) = c("terrGDEarea", "aquaGDEarea", 
                     "terrFrac", "aquaFrac", "allFrac", 
                     "arid", "topo", "seasonality", "lat", "modgrad",  "landform", 
                     "area")

stack_df = stack_ras |> 
  as_tibble() |> 
  set_colnames(c("terrGDEarea", "aquaGDEarea", 
                 "terrFrac", "aquaFrac", "allFrac", 
                 "arid", "topo", "seasonality", "lat", "modgrad",  "landform", 
                 "area")) |> 
  drop_na()


## plot params
modgrad_threshold = 0.2
stack_df = stack_df  |> filter(modgrad < modgrad_threshold)

##########
##
## make three plots: 
## (1) elevation x aridity - this is the only one used in the manuscript
## (2) elevation x seasonality
## (3) seasonality x aridity
##
##########

################################### 
################################### PLOT 1 - elev X aridity
################################### 
summary_df = stack_df |> 
  mutate(topo_bin = cut(topo,
                        breaks = seq(0, 5000, length.out = 51),
                        include.lowest = TRUE,
                        right = FALSE),
         arid_bin = cut(arid,
                        breaks = (10^seq(
                          log10(0.01),
                          log10(10),
                          length.out = 50 + 1   # +1 because breaks define bin edges
                        )),
                        include.lowest = TRUE,
                        right = FALSE)) |> 
  group_by(topo_bin, arid_bin) |> 
  summarise(
    n = n(),
    area_sum = sum(area, na.rm = TRUE),
    terrFrac_wmed = matrixStats::weightedMedian(terrFrac, w = area, na.rm = TRUE),
    aquaFrac_wmed = matrixStats::weightedMedian(aquaFrac, w = area, na.rm = TRUE),
    allFrac_wmed  = matrixStats::weightedMedian(allFrac,  w = area, na.rm = TRUE),
    .groups = "drop"
  )
plot_df = summary_df |> 
  mutate(
    arid_lo  = as.numeric(str_match(as.character(arid_bin),  "\\[([^,]+),")[,2]),
    arid_hi  = as.numeric(str_match(as.character(arid_bin),  ",([^\\)\\]]+)")[,2]),
    topo_lo = as.numeric(str_match(as.character(topo_bin), "\\[([^,]+),")[,2]),
    topo_hi = as.numeric(str_match(as.character(topo_bin), ",([^\\)\\]]+)")[,2])
  )

ggplot(plot_df) +
  geom_rect(aes(
    xmin = arid_lo, xmax = arid_hi,
    ymin = topo_lo, ymax = topo_hi,
    fill = allFrac_wmed
  )) +
  scale_x_log10(
    limits = c(0.01, 10)   # IMPORTANT: no zeros allowed
  ) + 
  scale_fill_gradientn(
    colours = met.brewer("Paquin", n = 200)[101:200],
    limits  = c(0, 1),
    breaks  = seq(0, 1, by = 0.01),
    oob     = scales::squish
  ) + 
  coord_cartesian(expand = FALSE, xlim = c(0.01, 10), ylim = c(0, 4000)) +
  labs(x = "elevation", 
       y = "aridity") +
  # geom_hline(yintercept = seq(0, 4000, 1000)) +
  # theme_minimal() 
  theme_void() + theme(legend.position = "none")

ggsave(file = paste0(here("plots"), "/tile_allFrac_wmed_Yelev_Xarid.png"),
       plot = last_plot(), device = "png",
       width = 250/4, height = 250/6, units = "mm", dpi = 400)



################################### 
################################### PLOT 2 - elev X seasonality
################################### 
summary_df = stack_df |> 
  mutate(topo_bin = cut(topo,
                        breaks = seq(0, 5000, length.out = 51),
                        include.lowest = TRUE,
                        right = FALSE),
         seasonality_bin = cut(seasonality,
                               breaks = seq(0, 2.5, length.out = 51),
                               include.lowest = TRUE,
                               right = FALSE)) |> 
  group_by(topo_bin, seasonality_bin) |> 
  summarise(
    n = n(),
    area_sum = sum(area, na.rm = TRUE),
    terrFrac_wmed = matrixStats::weightedMedian(terrFrac, w = area, na.rm = TRUE),
    aquaFrac_wmed = matrixStats::weightedMedian(aquaFrac, w = area, na.rm = TRUE),
    allFrac_wmed  = matrixStats::weightedMedian(allFrac,  w = area, na.rm = TRUE),
    .groups = "drop"
  )
plot_df = summary_df |> 
  mutate(
    seasonality_lo  = as.numeric(str_match(as.character(seasonality_bin),  "\\[([^,]+),")[,2]),
    seasonality_high  = as.numeric(str_match(as.character(seasonality_bin),  ",([^\\)\\]]+)")[,2]),
    topo_lo = as.numeric(str_match(as.character(topo_bin), "\\[([^,]+),")[,2]),
    topo_hi = as.numeric(str_match(as.character(topo_bin), ",([^\\)\\]]+)")[,2])
  )

ggplot(plot_df) +
  geom_rect(aes(
    xmin = seasonality_lo, xmax = seasonality_high,
    ymin = topo_lo, ymax = topo_hi,
    fill = allFrac_wmed
  )) + 
  scale_fill_gradientn(
    colours = met.brewer("Paquin", n = 200)[101:200],
    limits  = c(0, 1),
    breaks  = seq(0, 1, by = 0.01),
    oob     = scales::squish
  ) + 
  coord_cartesian(expand = FALSE, xlim = c(0.05, 2), ylim = c(0, 4000)) +
  labs(x = "seasonality", 
       y = "aridity") +
  # geom_vline(xintercept = c(0.5, 1, 1.5, 2)) +
  # theme_minimal()
  theme_void() + theme(legend.position = "none")

ggsave(file = paste0(here("plots"), "/tile_allFrac_wmed_Yelev_Xseasonality.png"), 
       plot = last_plot(), device = "png", 
       width = 250/4, height = 250/6, units = "mm", dpi = 400)


################################### 
################################### PLOT 3 - seasonality x aridity
################################### 
summary_df = stack_df |> 
  mutate(seasonality_bin = cut(seasonality,
                               breaks = seq(0, 2.5, length.out = 51),
                               include.lowest = TRUE,
                               right = FALSE),
         arid_bin = cut(arid,
                        breaks = (10^seq(
                          log10(0.01),
                          log10(10),
                          length.out = 50 + 1   # +1 because breaks define bin edges
                        )),
                        include.lowest = TRUE,
                        right = FALSE)) |> 
  group_by(seasonality_bin, arid_bin) |> 
  summarise(
    n = n(),
    area_sum = sum(area, na.rm = TRUE),
    terrFrac_wmed = matrixStats::weightedMedian(terrFrac, w = area, na.rm = TRUE),
    aquaFrac_wmed = matrixStats::weightedMedian(aquaFrac, w = area, na.rm = TRUE),
    allFrac_wmed  = matrixStats::weightedMedian(allFrac,  w = area, na.rm = TRUE),
    .groups = "drop"
  )
plot_df = summary_df |> 
  mutate(
    seasonality_lo  = as.numeric(str_match(as.character(seasonality_bin),  "\\[([^,]+),")[,2]),
    seasonality_high  = as.numeric(str_match(as.character(seasonality_bin),  ",([^\\)\\]]+)")[,2]),
    arid_lo = as.numeric(str_match(as.character(arid_bin), "\\[([^,]+),")[,2]),
    arid_hi = as.numeric(str_match(as.character(arid_bin), ",([^\\)\\]]+)")[,2])
  )

ggplot(plot_df) +
  geom_rect(aes(
    xmin = arid_lo, xmax = arid_hi,
    ymin = seasonality_lo, ymax = seasonality_high,
    fill = allFrac_wmed
  )) + 
  scale_fill_gradientn(
    colours = met.brewer("Paquin", n = 200)[101:200],
    limits  = c(0, 1),
    breaks  = seq(0, 1, by = 0.01),
    oob     = scales::squish
  ) +
  scale_x_log10(
    limits = c(0.01, 10)   # IMPORTANT: no zeros allowed
  ) + 
  coord_cartesian(expand = FALSE, xlim = c(0.01, 10), ylim = c(0, 2)) +
  labs(x = "aridity", 
       y = "seasonality") +
  # geom_hline(yintercept = c(0.5, 1, 1.5, 2)) +
  # theme_minimal()
  theme_void() + theme(legend.position = "none")

ggsave(file = paste0(here("plots"), "/tile_allFrac_wmed_Yseasonality_Xaridity.png"), 
       plot = last_plot(), device = "png", 
       width = 250/4, height = 250/6, units = "mm", dpi = 400)
