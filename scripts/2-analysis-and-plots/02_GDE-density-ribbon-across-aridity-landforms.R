### ---------------------\\ 
# Script objective:
# create landform stratified curves of GDE percentile across aridity and seasonality 
### ---------------------\\
library(here); source(here(("scripts/on_button.R")))
###

####################
## Import the Huggins et al. GDE types 
####################
R.gde_types_hu = terra::rast(here("data/GDE_types_hu_for_comparison.tif"))

####################
## now import aridity, seasonality, and landform data layers
####################
aridity = terra::rast("D:/Geodatabase/Climate/Aridity/ai_et0/ai_et0.tif") / 10000
aridity_30m = terra::resample(x = aridity, y = R.gde_types_hu, method = "bilinear")
aridity_30m_r = aridity_30m
aridity_30m_r[aridity_30m > 2] = 2

# Aridity legend
# 0.03 Hyper Arid
# 0.03 – 0.2 Arid
# 0.2 – 0.5 Semi-Arid
# 0.5 – 0.65 Dry sub-humid
# > 0.65 Humid

R.landform = terra::rast("D:/Geodatabase/Landforms/Global_Landforms_0d5_mode.tif")
R.landform = terra::resample(x = R.landform, y = R.gde_types_hu, method = "near")

R.seasonailty = terra::rast("D:/Geodatabase/Climate/WorldClim/precip_cv_0d5.tif")

hlz = terra::vect("D:/Geodatabase/Climate/Holdridge_Life_Zones/holdridge.shp")

hlz_r = terra::rasterize(x = hlz, y = R.gde_types_hu, field = "zone", touches = TRUE)
hlz_r[hlz_r == 0] = NA

coll_df = c(R.gde_types_hu$terr_GDE_area, R.gde_types_hu$aqua_GDE_area, aridity_30m, R.landform, R.seasonailty, hlz_r, WGS84_areaRaster(0.5) |> rast()) |> 
  as_tibble() |> 
  drop_na() |> 
  set_colnames(c("terr_GDE_area", "aqua_GDE_area", "aridity", "landform", "seasonality", "lifezones", "area")) |> 
  mutate(
    aridity_bin = cut(aridity, breaks = c(-Inf, 0.03, 0.2, 0.5, 0.65, Inf), labels = seq(1,5)) |> as.numeric()
  ) |> 
  filter(terr_GDE_area > 0) |> 
  mutate(
    terr_GDE_dens = 100 * (terr_GDE_area / area),
    aqua_GDE_dens = 100 * (aqua_GDE_area / area)
  ) 

# small script to assign a weighted percentile
wtd_pctile = function(x, w) {
  brks = Hmisc::wtd.quantile(x, weights = w, probs = seq(0, 1, 0.01), na.rm = TRUE)
  cut(x, breaks = unique(brks), include.lowest = TRUE, labels = FALSE)
}

coll_df$terr_ptile = wtd_pctile(coll_df$terr_GDE_dens, coll_df$area)
coll_df$aqua_ptile = wtd_pctile(coll_df$aqua_GDE_dens, coll_df$area)


######
## [plot set 1]
## Now summarise GDE-terr and GDE-aqua percentiles across aridity gradient
######

aridity_df = coll_df |> 
  mutate(arid_bin = cut(aridity,
                        breaks = (10^seq(
                          log10(0.01),
                          log10(10),
                          length.out = 50 + 1   # +1 because breaks define bin edges
                        )),
                        include.lowest = TRUE,
                        right = FALSE)) |> 
  group_by(arid_bin, landform) |> 
  summarise(
    n = n(),
    area_sum = sum(area, na.rm = TRUE),
    terr_ptile_wmed = matrixStats::weightedMedian(terr_ptile, w = area, na.rm = TRUE),
    aqua_ptile_wmed = matrixStats::weightedMedian(aqua_ptile, w = area, na.rm = TRUE),
    .groups = "drop"
  ) |> 
  mutate(
    arid_lo  = as.numeric(str_match(as.character(arid_bin),  "\\[([^,]+),")[,2]),
    arid_hi  = as.numeric(str_match(as.character(arid_bin),  ",([^\\)\\]]+)")[,2])
  ) |> 
  mutate(
    arid_med = (arid_lo+arid_hi)/2
  )

library(ggh4x)

# Landform legend
# 1 - mountains 
# 2 - hills
# 3 - plateaus
# 4 - plains

aridity_df |>
  filter(landform == 1) |>
  ggplot(aes(x = arid_med, ymin = terr_ptile_wmed, ymax = aqua_ptile_wmed)) +
  stat_difference(alpha = 0.3) +
  geom_line(aes(y = terr_ptile_wmed), color = "#2C5224", linewidth = 1) +
  geom_line(aes(y = aqua_ptile_wmed), color = "#75BCCD", linewidth = 1) +
  scale_fill_manual(values = c("+" = "#75BCCD", "-" = "#2C5224")) +
  scale_x_log10(limits = c(0.01, 3)) +
  scale_y_continuous(limits = c(0, 100)) +
  coord_cartesian(expand = FALSE) +
  theme_void() + theme(legend.position = "none") 
  # geom_vline(xintercept = c(0.03, 0.2, 0.5, 0.65, 1))

ggsave(file = paste0(here("plots"), "/RIBBON_aridity_MOUNTAINS.png"),
       plot = last_plot(), device = "png",
       width = 250/4, height = 250/4, units = "mm", dpi = 400)


######
## [plot set 2]
## Now summarise GDE-terr and GDE-aqua percentiles across seasonality gradient
######

seasonality_df = coll_df |> 
  mutate(seasonality_bin = cut(seasonality,
                               breaks = seq(0, 2.5, length.out = 51),
                               include.lowest = TRUE,
                               right = FALSE)) |>  
  group_by(seasonality_bin, landform) |> 
  summarise(
    n = n(),
    area_sum = sum(area, na.rm = TRUE),
    terr_ptile_wmed = matrixStats::weightedMedian(terr_ptile, w = area, na.rm = TRUE),
    aqua_ptile_wmed = matrixStats::weightedMedian(aqua_ptile, w = area, na.rm = TRUE),
    .groups = "drop"
  ) |> 
  mutate(
    seasonality_lo  = as.numeric(str_match(as.character(seasonality_bin),  "\\[([^,]+),")[,2]),
    seasonality_hi  = as.numeric(str_match(as.character(seasonality_bin),  ",([^\\)\\]]+)")[,2])
  ) |> 
  mutate(
    seasonality_med = (seasonality_lo+seasonality_hi)/2
  )

# Landform legend
# 1 - mountains 
# 2 - hills
# 3 - plateaus
# 4 - plains

seasonality_df |>
  filter(landform == 4) |>
  ggplot(aes(x = seasonality_med, ymin = terr_ptile_wmed, ymax = aqua_ptile_wmed)) +
  stat_difference(alpha = 0.3) +
  geom_line(aes(y = terr_ptile_wmed), color = "#2C5224", linewidth = 1) +
  geom_line(aes(y = aqua_ptile_wmed), color = "#75BCCD", linewidth = 1) +
  scale_fill_manual(values = c("+" = "#75BCCD", "-" = "#2C5224")) +
  scale_x_continuous(limits = c(0.1, 1.2)) +
  scale_y_continuous(limits = c(0, 100)) +
  coord_cartesian(expand = FALSE) +
  theme_void() + theme(legend.position = "none") 
  # geom_vline(xintercept = c(0.3, 0.6, 0.9, 1.2))

ggsave(file = paste0(here("plots"), "/RIBBON_seasonality_PLAINS.png"),
       plot = last_plot(), device = "png",
       width = 250/4, height = 250/4, units = "mm", dpi = 400)



######
## [plot set 3]
## General relationships between aqua and terrestrial GDE density and density percentiles
######

cor.test(coll_df$aqua_GDE_dens, coll_df$terr_GDE_dens)

cor.test(coll_df$aqua_GDE_dens, coll_df$terr_GDE_dens, method = "spearman")

coll_df |>
  mutate(aqua_bin = cut(aqua_GDE_dens, breaks = seq(0, 100, by = 2), include.lowest = TRUE)) |>
  group_by(aqua_bin) |>
  summarise(
    aqua_mid = mean(aqua_GDE_dens, na.rm = TRUE),
    q25 = quantile(terr_GDE_dens, 0.25, na.rm = TRUE),
    q50 = quantile(terr_GDE_dens, 0.50, na.rm = TRUE),
    q75 = quantile(terr_GDE_dens, 0.75, na.rm = TRUE),
    q05 = quantile(terr_GDE_dens, 0.05, na.rm = TRUE),
    q95 = quantile(terr_GDE_dens, 0.95, na.rm = TRUE),
    .groups = "drop"
  ) |>
  ggplot(aes(x = aqua_mid)) +
  geom_ribbon(aes(ymin = q05, ymax = q95), fill = "grey80", alpha = 0.5) +
  geom_ribbon(aes(ymin = q25, ymax = q75), fill = "grey40", alpha = 0.5) +
  geom_line(aes(y = q50), color = "black", linewidth = 1) +
  theme_void() +
  coord_cartesian(expand = FALSE, ylim = c(0,100), xlim =c(0,100))

ggsave(file = paste0(here("plots"), "/RIBBON_correlation_GDE_t_a.png"),
       plot = last_plot(), device = "png",
       width = 250/4, height = 250/4, units = "mm", dpi = 400)


######
## [plot set 4]
## GDE percentiles across Holdridge life forms
######

hlz_lookup = data.frame(zone = hlz$zone, desc = hlz$desc_) |> unique()

zone_summary = coll_df |>
  filter(lifezones > 1) |>
  bind_rows(
    coll_df |> filter(lifezones == 38) |> mutate(lifezones = 39)
  ) |>
  mutate(
    lifezones = case_when(
      lifezones == 25 ~ 18,
      lifezones == 26 ~ 19,
      lifezones == 27 ~ 20,
      lifezones == 28 ~ 21,
      lifezones == 29 ~ 22,
      lifezones == 30 ~ 23,
      lifezones == 31 ~ 24,
      TRUE ~ lifezones
    )
  ) |>
  group_by(lifezones) |>
  summarise(
    terr_med = matrixStats::weightedMedian(terr_ptile, w = area, na.rm = TRUE),
    aqua_med = matrixStats::weightedMedian(aqua_ptile, w = area, na.rm = TRUE),
    terr_area = sum(terr_GDE_area, na.rm = T),
    aqua_area = sum(aqua_GDE_area, na.rm = T),
    .groups = "drop"
  ) |>
  mutate(
    diff_med = terr_med - aqua_med
  ) |> 
  left_join(hlz_lookup, by = c("lifezones" = "zone")) |>
  mutate(
    desc = case_when(
      lifezones == 18 ~ "Warm temperate / Subtropical desert",
      lifezones == 19 ~ "Warm temperate / Subtropical desert bush",
      lifezones == 20 ~ "Warm temperate / Subtropical thorn steppe",
      lifezones == 21 ~ "Warm temperate / Subtropical dry forest",
      lifezones == 22 ~ "Warm temperate / Subtropical moist forest",
      lifezones == 23 ~ "Warm temperate / Subtropical wet forest",
      lifezones == 24 ~ "Warm temperate / Subtropical rain forest",
      lifezones == 39 ~ "Tropical rain forest (placeholder)",
      TRUE ~ desc
    )
  )

hlz_positions = tribble(
  ~zone, ~row, ~col,
  2,  1, 5,
  3,  2, 4,  4,  2, 5,  5,  2, 6,  6,  2, 7,
  7,  3, 3,  8,  3, 4,  9,  3, 5, 10,  3, 6, 11,  3, 7,
  12,  4, 3, 13,  4, 4, 14,  4, 5, 15,  4, 6, 16,  4, 7, 17,  4, 8,
  18,  5, 2, 19,  5, 3, 20,  5, 4, 21,  5, 5, 22,  5, 6, 23,  5, 7, 24,  5, 8,
  32,  6, 2, 33,  6, 3, 34,  6, 4, 35,  6, 5, 36,  6, 6, 37,  6, 7, 38,  6, 8, 39,  6, 9
)

zone_summary <- zone_summary |>
  left_join(hlz_positions, by = c("lifezones" = "zone")) |>
  mutate(
    x = col + ifelse(row %% 2 == 1, 0.5, 0),
    y = -(row - 1) * sqrt(3) / 2
  ) |>
  filter(lifezones != 2)

zone_summary$GDE_total_med = ((zone_summary$terr_med * zone_summary$terr_area) + 
                                (zone_summary$aqua_med * zone_summary$aqua_area))/ (zone_summary$terr_area + zone_summary$aqua_area) 

full_grid <- tribble(
  ~row, ~col,
  1, 4, 1, 5, 1, 6,
  2, 4, 2, 5, 2, 6, 2, 7,
  3, 3, 3, 4, 3, 5, 3, 6, 3, 7,
  4, 3, 4, 4, 4, 5, 4, 6, 4, 7, 4, 8,
  5, 2, 5, 3, 5, 4, 5, 5, 5, 6, 5, 7, 5, 8,
  6, 2, 6, 3, 6, 4, 6, 5, 6, 6, 6, 7, 6, 8, 6, 9
) |>
  mutate(
    x = col + ifelse(row %% 2 == 1, 0.5, 0),
    y = -(row - 1) * sqrt(3) / 2
  )


# === Plot (change terr_med to aqua_med as needed) ===
ggplot() +
  ggforce::geom_regon(data = zone_summary,
                      aes(x0 = x, y0 = y, sides = 6, r = 1/sqrt(3), angle = pi/2,
                          fill = GDE_total_med),
                      color = NA) +
  ggforce::geom_regon(data = full_grid,
                      aes(x0 = x, y0 = y, sides = 6, r = 1/sqrt(3), angle = pi/2),
                      fill = NA, color = "black", linewidth = 0.4) +
  scale_fill_gradientn(colors = met.brewer("Paquin", n = 100)[51:100], 
                       limits = c(0, 100), oob = scales::squish) +
  geom_text(data = zone_summary,
            aes(x = x, y = y, label = str_wrap(desc, 10)),
            size = 2, lineheight = 0.8) +
  coord_equal() +
  theme_void() +
  theme(legend.position = "none")

ggsave(file = paste0(here("plots"), "/HEX_holdridge_total_med.png"),
       plot = last_plot(), device = "png",
       width = 250/4, height = 250/4, units = "mm", dpi = 400)


ggplot() +
  ggforce::geom_regon(data = zone_summary,
                      aes(x0 = x, y0 = y, sides = 6, r = 1/sqrt(3), angle = pi/2,
                          fill = diff_med),
                      color = NA) +
  ggforce::geom_regon(data = full_grid,
                      aes(x0 = x, y0 = y, sides = 6, r = 1/sqrt(3), angle = pi/2),
                      fill = NA, color = "black", linewidth = 0.4) +
  scale_fill_gradientn(colors = met.brewer("Isfahan1", n = 100, direction = -1), 
                       limits = c(-15, 15), oob = scales::squish) +
  # geom_text(data = zone_summary,
  #           aes(x = x, y = y, label = str_wrap(desc, 10)),
  #           size = 2, lineheight = 0.8) +
  coord_equal() +
  theme_void() +
  theme(legend.position = "none")

ggsave(file = paste0(here("plots"), "/HEX_holdridge_diff.png"),
       plot = last_plot(), device = "png",
       width = 250/4, height = 250/4, units = "mm", dpi = 400)




#### now plot scatter of hexplots lifezones for GDE area density and area total per life zone

hlz_lookup = data.frame(zone = hlz$zone, desc = hlz$desc_) |> unique()

lifezone_df = coll_df |> 
  group_by(lifezones) |> 
  summarise(
    terr_GDE_area = sum(terr_GDE_area, na.rm = T),
    aqua_GDE_area = sum(aqua_GDE_area, na.rm = T),
    total_area = sum(area, na.rm = T)
  ) |> 
  mutate(
    terr_dens = terr_GDE_area / total_area,
    aqua_dens = aqua_GDE_area / total_area
  ) |> 
  left_join(hlz_lookup, by = c("lifezones" = "zone")) |>
  mutate(
    desc = case_when(
      lifezones == 18 ~ "Warm temperate / Subtropical desert",
      lifezones == 19 ~ "Warm temperate / Subtropical desert bush",
      lifezones == 20 ~ "Warm temperate / Subtropical thorn steppe",
      lifezones == 21 ~ "Warm temperate / Subtropical dry forest",
      lifezones == 22 ~ "Warm temperate / Subtropical moist forest",
      lifezones == 23 ~ "Warm temperate / Subtropical wet forest",
      lifezones == 24 ~ "Warm temperate / Subtropical rain forest",
      lifezones == 39 ~ "Tropical rain forest (placeholder)",
      TRUE ~ desc
    )
  )


# Lookup: desc = color
desc_colors <- tribble(
  ~desc,                                            ~color,
  "Polar dry tundra",                               "#C0BFBD",
  "Polar moist tundra",                             "#C0BFBD",
  "Polar wet tundra",                               "#C0BFBD",
  "Polar rain tundra",                              "#C0BFBD",
  "Boreal desert",                                  "#D2DCE5",
  "Boreal dry bush",                                "#D2DCE5",
  "Boreal moist forest",                            "#586393",
  "Boreal wet forest",                              "#586393",
  "Boreal rain forest",                             "#586393",
  "Cool temperate desert",                          "#CCCC76",
  "Cool temperate desert bush",                     "#CCCC76",
  "Cool temperate steppe",                          "#CCCC76",
  "Cool temperate moist forest",                    "#ADC9CA",
  "Cool temperate wet forest",                      "#ADC9CA",
  "Cool temperate rain forest",                     "#ADC9CA",
  "Warm temperate / Subtropical desert",            "#DEDDA5",
  "Warm temperate / Subtropical desert bush",       "#DEDDA5",
  "Warm temperate / Subtropical thorn steppe",      "#DEDDA5",
  "Warm temperate / Subtropical dry forest",        "#968649",
  "Warm temperate / Subtropical moist forest",      "#968649",
  "Warm temperate / Subtropical wet forest",        "#968649",
  "Warm temperate / Subtropical rain forest",       "#968649",
  "Tropical desert",                                "#D3A372",
  "Tropical desert bush",                           "#D3A372",
  "Tropical thorn steppe",                          "#D3A372",
  "Tropical very dry forest",                       "#AE6E6A",
  "Tropical dry forest",                            "#AE6E6A",
  "Tropical moist forest",                          "#804946",
  "Tropical wet forest",                            "#804946",
  "Tropical rain forest (placeholder)",             "#804946"
)

lifezone_df = lifezone_df |> left_join(desc_colors, by = "desc") |> 
  drop_na()

ggplot(lifezone_df) +
  ggforce::geom_regon(aes(x0 = 100*aqua_dens/2, y0 = 100*terr_dens,
                          sides = 6, r = 2.5, angle = pi/2, fill = color),
                      color = "black", linewidth = 0.3) +
  scale_fill_identity() +
  coord_fixed(ratio = 1, xlim = c(0, 50), ylim = c(0, 50), expand = FALSE) +
  theme_void() +
  theme(legend.position = "none") 
  # geom_hline(yintercept = c(0,50)) +
  # geom_vline(xintercept = c(0,100))

ggsave(file = paste0(here("plots"), "/HEXscatter_holdridge_densities.png"),
       plot = last_plot(), device = "png",
       width = 250/4, height = 250/4, units = "mm", dpi = 400)


#########
######### spearmann correlation per life zone
#########

coll_prep = coll_df |>
  filter(lifezones > 1) |>
  bind_rows(
    coll_df |> filter(lifezones == 38) |> mutate(lifezones = 39)
  ) |>
  mutate(
    lifezones = case_when(
      lifezones == 25 ~ 18, lifezones == 26 ~ 19, lifezones == 27 ~ 20,
      lifezones == 28 ~ 21, lifezones == 29 ~ 22, lifezones == 30 ~ 23,
      lifezones == 31 ~ 24, TRUE ~ lifezones
    )
  )

# --- Summarise per lifezone ---
zone_summary = coll_prep |>
  group_by(lifezones) |>
  summarise(
    cor_res   = list(tryCatch(
      cor.test(aqua_GDE_dens, terr_GDE_dens, method = "spearman", exact = FALSE),
      error = function(e) NULL
    )),
    n = n(),
    .groups = "drop"
  ) |>
  mutate(
    rho     = sapply(cor_res, \(x) if (!is.null(x)) unname(x$estimate) else NA_real_),
    p_value = sapply(cor_res, \(x) if (!is.null(x)) x$p.value          else NA_real_),
    sig     = p_value < 0.05
  ) |>
  dplyr::select(-cor_res) |>
  left_join(hlz_lookup, by = c("lifezones" = "zone")) |>
  mutate(
    desc = case_when(
      lifezones == 18 ~ "Warm temperate / Subtropical desert",
      lifezones == 19 ~ "Warm temperate / Subtropical desert bush",
      lifezones == 20 ~ "Warm temperate / Subtropical thorn steppe",
      lifezones == 21 ~ "Warm temperate / Subtropical dry forest",
      lifezones == 22 ~ "Warm temperate / Subtropical moist forest",
      lifezones == 23 ~ "Warm temperate / Subtropical wet forest",
      lifezones == 24 ~ "Warm temperate / Subtropical rain forest",
      lifezones == 39 ~ "Tropical rain forest (placeholder)",
      TRUE ~ desc
    )
  )

# --- Grid positions (unchanged) ---
hlz_positions = tribble(
  ~zone, ~row, ~col,
  2,  1, 5,
  3,  2, 4,  4,  2, 5,  5,  2, 6,  6,  2, 7,
  7,  3, 3,  8,  3, 4,  9,  3, 5, 10,  3, 6, 11,  3, 7,
  12,  4, 3, 13,  4, 4, 14,  4, 5, 15,  4, 6, 16,  4, 7, 17,  4, 8,
  18,  5, 2, 19,  5, 3, 20,  5, 4, 21,  5, 5, 22,  5, 6, 23,  5, 7, 24,  5, 8,
  32,  6, 2, 33,  6, 3, 34,  6, 4, 35,  6, 5, 36,  6, 6, 37,  6, 7, 38,  6, 8, 39,  6, 9
)

full_grid = tribble(
  ~row, ~col,
  1, 4, 1, 5, 1, 6,
  2, 4, 2, 5, 2, 6, 2, 7,
  3, 3, 3, 4, 3, 5, 3, 6, 3, 7,
  4, 3, 4, 4, 4, 5, 4, 6, 4, 7, 4, 8,
  5, 2, 5, 3, 5, 4, 5, 5, 5, 6, 5, 7, 5, 8,
  6, 2, 6, 3, 6, 4, 6, 5, 6, 6, 6, 7, 6, 8, 6, 9
) |>
  mutate(x = col + ifelse(row %% 2 == 1, 0.5, 0), y = -(row - 1) * sqrt(3) / 2)

# --- Attach coordinates ---
zone_summary = zone_summary |>
  left_join(hlz_positions, by = c("lifezones" = "zone")) |>
  mutate(x = col + ifelse(row %% 2 == 1, 0.5, 0), y = -(row - 1) * sqrt(3) / 2) |>
  filter(lifezones != 2)

# --- Plot ---
ggplot() +
  ggforce::geom_regon(data = zone_summary,
                      aes(x0 = x, y0 = y, sides = 6, r = 1/sqrt(3), angle = pi/2,
                          fill = rho),
                      color = NA) +
  ggforce::geom_regon(data = full_grid,
                      aes(x0 = x, y0 = y, sides = 6, r = 1/sqrt(3), angle = pi/2),
                      fill = NA, color = "black", linewidth = 0.4) +
  scale_fill_gradientn(colors = met.brewer("Greek", n = 100, direction = -1),
                       limits = c(0, 1), oob = scales::squish, na.value = "grey80") +
  # geom_text(data = zone_summary,
  #           aes(x = x, y = y, label = str_wrap(desc, 10)),
  #           size = 2, lineheight = 0.8) +
  coord_equal() +
  theme_void() +
  theme(legend.position = "none")

ggsave(file = paste0(here("plots"), "/HEXscatter_holdridge_spearmann_cor.png"),
       plot = last_plot(), device = "png",
       width = 250/4, height = 250/4, units = "mm", dpi = 400)


zone_means <- coll_df |>
  group_by(lifezones) |>
  summarise(
    mean_aqua = mean(aqua_GDE_dens, na.rm = TRUE),
    mean_terr = mean(terr_GDE_dens, na.rm = TRUE)
  )

cor.test(zone_means$mean_aqua, zone_means$mean_terr, method = "spearman")