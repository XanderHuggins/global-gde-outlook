### ---------------------\\ 
# Script objective:
# Prepare raster stack of co-occuring GDE presures
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

## 2. GW depletion
GW_dep = terra::rast("D:/D_documents/1.projects-scripts/GWS_NDVI_resilience/data/grace_gws_ts/gwsa_theilsen_slope.tif")
GW_dep[GW_dep > 0] = 0
GW_dep = abs(GW_dep)
GW_dep_0d5 = terra::resample(x = GW_dep, y = template_r, "near")

## 3. land conversion pressure
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

writeRaster(pressure_stack, here("data/gde_pressure_stack.tif"))