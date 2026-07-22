### ---------------------\\ 
# Script objective:
# reproject CPI to WGS84 and resample to 0d5
### ---------------------\\

gdalwarp(
  srcfile  = "D:/Geodatabase/Land-use/CPI_Ver1_Oakleaf_etal_2024/CPI_Continous_0to1.tif",
  dstfile  = "D:/Geodatabase/Land-use/CPI_Ver1_Oakleaf_etal_2024/CPI_Continous_0to1_WGS84.tif",
  s_srs    = "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs",
  t_srs    = "EPSG:4326",
  te       = c(-180, -90, 180, 90),
  tr       = c(0.5/60, 0.5/60),
  r        = "average",
  dstnodata = -9999,
  multi    = TRUE
)

CPI = terra::rast("D:/Geodatabase/Land-use/CPI_Ver1_Oakleaf_etal_2024/CPI_Continous_0to1_WGS84.tif")

area_rast = WGS84_areaRaster(0.5/60) |> rast()
writeRaster(area_rast, filename = here("data/area_raster_30arcsec.tif"))

CPI_x_area = CPI * area_rast

# half a degree
CPI_x_area_0d5 = terra::aggregate(x = CPI_x_area,
                                  fact = 60,
                                  fun = sum,
                                  cores = 7)

CPI_05deg = CPI_x_area_0d5 / WGS84_areaRaster(0.5) |> rast()
writeRaster(CPI_05deg, filename = here("data/CPI_05deg.tif"))


# 5 arcminute
CPI_x_area_5am = terra::aggregate(x = CPI_x_area,
                                  fact = 10,
                                  fun = sum,
                                  cores = 7)

CPI_5am = CPI_x_area_5am / WGS84_areaRaster(5/60) |> rast()
writeRaster(CPI_5am, filename = here("data/CPI_5am.tif"))