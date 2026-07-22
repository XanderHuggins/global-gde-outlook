### ---------------------\\ 
# Script objective:
# Extract and resample precipitation seasonality to 0d5
### ---------------------\\

# BIO15 = Precipitation Seasonality (Coefficient of Variation)
p_cv = terra::rast("D:/Geodatabase/Climate/WorldClim/wc2.1_10m_bio_15.tif")/100

# reference raster
R.comp_gde = terra::rast(here("data/GDE_comparison_0d5.tif"), lyrs =1 )

p_cv = terra::resample(x = p_cv, y = R.comp_gde, method = "bilinear")
writeRaster(p_cv, "D:/Geodatabase/Climate/WorldClim/precip_cv_0d5.tif")

aridity = terra::rast("D:/Geodatabase/Climate/Aridity/ai_et0/ai_et0.tif") / 10000
aridity_30m = terra::resample(x = aridity, y = R.comp_gde, method = "bilinear")
writeRaster(aridity_30m, "D:/Geodatabase/Climate/Aridityaridity_0d5.tif")