### ---------------------\\
# Script objective:
# Import custom functions
### ---------------------\\

### ---
# NAME: WGS84_areaRaster
# Description: Generate raster with cell values representing grid areas based on WGS84 reference ellipsoid.
# This method is more sophisticated approach than that provided by the e.g. raster::area() function and follows methods outlined in Santini et al. 2010 https://doi.org/10.1111/j.1467-9671.2010.01200.x
WGS84_areaRaster = function(ResT) {
  # Function arguments:
  # ResT: Desired resolution of area raster in decimal degrees
  
  pi.g = 3.14159265358979
  Fl = 0.00335281066474 # Flattening
  SMA = 6378137.0 # Semi major axis
  e = sqrt((2*Fl) - (Fl^2)) # Eccentricity
  RES = ResT
  
  # Check that resolution is a factor of 90
  if ((90/RES) %% 1 > 1e-6) { stop("'ResT' must a factor of 90") }
  
  # Initialize dataframe with geodetic latitudes
  df_a = data.frame(LowLAT = seq(-90, 90-RES, by = RES),
                    UppLAT = seq(-90+RES, 90, by = RES))
  
  # Convert geodetic latitudes from degrees to radians
  df_a$LowLATrad = df_a$LowLAT * pi.g / 180
  df_a$UppLATrad = df_a$UppLAT * pi.g / 180
  
  # Calculate q1 and q2 following Santini et al. (2010)
  df_a$q1 = (1-e*e)* ((sin(df_a$LowLATrad)/(1-e*e*sin(df_a$LowLATrad)^2)) - ((1/(2*e))*log((1-e*sin(df_a$LowLATrad))/(1+e*sin(df_a$LowLATrad)))))
  
  df_a$q2 = (1-e*e)* ((sin(df_a$UppLATrad)/(1-e*e*sin(df_a$UppLATrad)^2)) - ((1/(2*e))*log((1-e*sin(df_a$UppLATrad))/(1+e*sin(df_a$UppLATrad)))))
  
  # Calculate q constant following Santini et al (2010)
  q_const = (1-e*e)* ((sin(pi.g/2)/(1-e*e*sin(pi.g/2)^2)) - ((1/(2*e))*log((1-e*sin(pi.g/2))/(1 + e*sin(pi.g/2)))))
  
  # Calculate authaltic latitudes
  df_a$phi1 = asin(df_a$q1 / q_const)
  df_a$phi2 = asin(df_a$q2 / q_const)
  
  # Calculate authaltic radius
  R.adius = sqrt(SMA*SMA*q_const/2)
  
  # Calculate cell size in radians
  CS = (RES) * pi.g/180
  
  # Calculate cell area in m2
  df_a$area_m2 = R.adius*R.adius*CS*(sin(df_a$phi2)-sin(df_a$phi1))
  
  # Convert to raster, and replicate column across global longitude domain
  WGS84area_km2 = matrix(df_a$area_m2/1e6, nrow = 180/RES, ncol = 360/RES,
                         byrow = FALSE, dimnames = NULL) |> raster()
  extent(WGS84area_km2) = c(-180, 180, -90, 90) # Set extent of raster
  crs(WGS84area_km2) = crs("+proj=longlat") # Set CRS of raster
  
  WGS84area_km2 = raster::flip(WGS84area_km2, direction = "y")
  
  message(paste0("Calculated global surface area at: ", RES,
                 "deg. is ", sum(WGS84area_km2[]), " km2.", sep = ""))
  
  return(WGS84area_km2)
}


### ---
# NAME: reclassify_percentile_bins
# Description: Function to reclassify raster into area-weighted percentile bins
reclassify_percentile_bins = function(raster_layer, area_grid, num_bins = 100) {
  
  D.stack = c(raster_layer, area_grid) |> as_tibble() |> drop_na() |> set_colnames(c("rl", "ga"))
  
  # Calculate percentiles based on the raster values
  percentiles = Hmisc::wtd.quantile(x = D.stack$rl, weights = D.stack$ga,
                                    probs = seq(0, 1, length.out = num_bins + 1), na.rm = TRUE)
  
  # reclassification matrix
  M.reclass = data.frame(low = c(-Inf, percentiles[1:num_bins-1]),
                         high = c(percentiles[2:num_bins], Inf),
                         new = seq(1:num_bins)) |> as.matrix()
  
  # Apply the reclassification to the raster
  R.reclass <- terra::classify(x = raster_layer, rcl = M.reclass)
  
  return(R.reclass)
}

### ---
# NAME: wtd_ptile
# Description: small script to assign a weighted percentile
wtd_pctile = function(x, w) {
  brks = Hmisc::wtd.quantile(x, weights = w, probs = seq(0, 1, 0.01), na.rm = TRUE)
  cut(x, breaks = unique(brks), include.lowest = TRUE, labels = FALSE)
}

### ---
# NAME: wtd_mode
# Description: small script to assign a weighted mode
wtd_mode = function(x, w) {
  ux = unique(x)
  ux[which.max(tapply(w, match(x, ux), sum))]
}
