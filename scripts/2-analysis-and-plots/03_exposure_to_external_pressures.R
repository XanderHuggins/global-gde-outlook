### ---------------------\\ 
# Script objective:
# Map co-occurrence of GDE pressures
### ---------------------\\
library(here); source(here(("scripts/on_button.R")))
###

template_r = terra::rast(extent = c(-180, 180, -90, 90),
                         resolution = 5/60)

## variables to cross map with bivariate plots

## 1. Groundwater warming
gw_warm = terra::rast("D:/Geodatabase/Groundwater/Temperature_warming/TempChange_2100to2000_ssp585_GWTable.tif")
gw_warm = gw_warm$ssp585_p50
gw_warm_0d5 = terra::aggregate(x = gw_warm, fact = 6, fun = "modal")
gw_warm_0d5 = terra::resample(x = gw_warm_0d5, y = template_r, "near")

## 4. GW depletion
GW_dep = terra::rast("C:/Users/xande/Documents/1.projects-scripts/GWS_NDVI_resilience/data/grace_gws_ts/gwsa_theilsen_slope.tif")
GW_dep[GW_dep > 0] = 0
GW_dep = abs(GW_dep)
GW_dep_0d5 = terra::resample(x = GW_dep, y = template_r, "near")


## 5. land conversion pressure
CPI = terra::rast(here("data/CPI_5am.tif"))
CPI_0d5 = terra::aggregate(CPI, fact = 6, fun = "mean", na.rm = T)
CPI_0d5 = terra::resample(x = CPI_0d5, y = template_r, "near")

# create binary representations 
gw_warm_0d5_BIN = gw_warm_0d5
gw_warm_0d5_BIN[gw_warm_0d5 < 2] = 0
gw_warm_0d5_BIN[gw_warm_0d5 >= 2] = 1

# create binary representations 
GW_dep_0d5_BIN = GW_dep_0d5
GW_dep_0d5_BIN[GW_dep_0d5 < 5] = 0
GW_dep_0d5_BIN[GW_dep_0d5 >= 5] = 1

# create binary representations 
CPI_0d5_BIN = CPI_0d5
CPI_0d5_BIN[CPI_0d5 < 0.5] = 0
CPI_0d5_BIN[CPI_0d5 >= 0.5] = 1


## create three layered ID
pressure_stack = rast(template_r)
pressure_stack = (100 * GW_dep_0d5_BIN) +
  (10 * gw_warm_0d5_BIN) + 
  CPI_0d5_BIN

pressure_stack[] = terra::as.factor(pressure_stack)

# colour scheme
category_colours = c(
  "0"   = "#d3d3d3",  # None - grey
  "1"   = "#FAB476",  # land use only 
  "10"  = "#FEE08B",  # Warming only 
  "100" = "#91BFDB",  # depletion only 
  "11"  = "#D73027",  # warming and land use 
  "101" = "#762A83",  # depletion and land use 
  "110" = "#3F00B4",  # depletion and warming
  "111" = "#1A1A1A"   # all 
)

writeRaster(pressure_stack, here("data/gde_pressure_stack.tif"))

pressure_stack = terra::rast(here("data/gde_pressure_stack.tif"))

#### map
outline = terra::vect("C:/Users/xande/Documents/1.projects-scripts/sustainability-puzzles/data/land_mask_polygon.sqlite") |> 
  st_as_sf()

# 1 Caspain sea
casp = terra::vect("D:/Geodatabase/Admin-ocean-boundaries/worldglwd1.shp")
casp_r = terra::rasterize(x = casp, y = pressure_stack, 1, touches = F)
pressure_stack[casp_r == 1] = NA

map =  
  tm_shape(outline, crs = "+proj=robin") +
  tm_fill(col ="grey", lwd = 0, border.col = NA) +
  tm_shape(terra::project(x = pressure_stack, y = "+proj=robin", method = "near")) +
  tm_raster(
    col.scale = tm_scale_categorical(
      values = category_colours
    )) + 
  # tm_shape(outline) + 
  # tm_borders(col = "black", lwd = 1) +
  tm_layout(legend.show = F, legend.frame = F, frame = F, bg = FALSE, bg.color = "transparent")
map
tmap_save(map, here("plots/gMAP_pressure_stacking.pdf"), dpi = 400, units = "in")


#######################################
### Create a Voronoi tree to represent the area distribution of these pressures
#######################################
library(WeightedTreemaps)


area_df = c(pressure_stack, WGS84_areaRaster(5/60) |> rast()) |> 
  as_tibble() |> 
  set_colnames(c('id', 'area'))

summary_df = area_df |> 
  group_by(id) |> 
  summarise(
    id_area = sum(area, na.rm = T)
  ) |> 
  filter(id >= 0) |> 
  mutate(
    id_frac = id_area / sum(id_area)
  ) |> 
  mutate(
    id = as.character(id)
  ) |> 
  arrange(-id_frac) |> 
  mutate(order_id = sprintf("%02d", row_number()))


# Compute the additively-weighted Voronoi treemap, clipped to a circle
tm = voronoiTreemap(
  data      = summary_df,
  levels    = "id",
  cell_size = "id_frac",   # "id_area" gives the same result — both get normalised
  shape     = "circle",
  positioning = "clustered",
  # sort        = FALSE,
  seed      = 7           # layout is stochastic; fix the seed for reproducibility
)

cells = get_polygons(tm)   # returns an sf object, one polygon per cell


cells_sf <- st_sf(
  id       = sub("^LEVEL\\d+_", "", names(cells)),  # "LEVEL1_100" -> "100"
  geometry = st_sfc(cells)                           # list of sfg -> one sfc
)

ggplot(cells_sf) +
  geom_sf(aes(fill = id), colour = "white", linewidth = 0.4) +
  scale_fill_manual(values = category_colours) +
  theme_void()

ggsave(file = here("plots/pressures_voronoi_tree.png"),
       plot = last_plot(), device = "png",
       width = 51.8, height = 34.6, units = "mm", dpi = 400)