### ---------------------\\ 
# Script objective:
# and evaluate patterns in GDE distributions across aridity and altitude bands per GDE class
### ---------------------\\
library(here); source(here(("scripts/on_button.R")))
###

####################
## Import the Huggins et al. GDE types 
####################
R.gde_types_hu = terra::rast(here("data/GDE_types_hu_for_comparison.tif"))

####################
## now import aridity and landform data layers
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

# Landform legend
# 1 - mountains 
# 2 - hills
# 3 - plateaus
# 4 - plains

R.seasonailty = terra::rast("D:/Geodatabase/Climate/WorldClim/precip_cv_0d5.tif")

coll_df = c(R.gde_types_hu$terr_GDE_area, R.gde_types_hu$aqua_GDE_area, aridity_30m, R.landform, R.seasonailty, WGS84_areaRaster(0.5) |> rast()) |> 
  as_tibble() |> 
  drop_na() |> 
  set_colnames(c("terr_GDE_area", "aqua_GDE_area", "aridity", "landform", "seasonality", "area")) |> 
  mutate(
    aridity_bin = cut(aridity, breaks = c(-Inf, 0.03, 0.2, 0.5, 0.65, Inf), labels = seq(1,5)) |> as.numeric()
  ) |> 
  filter(terr_GDE_area > 0) |> 
  mutate(
    terr_GDE_dens = 100 * (terr_GDE_area / area),
    aqua_GDE_dens = 100 * (aqua_GDE_area / area)
  )



coll_df = coll_df |> 
  mutate(
    terr_GDE_Ptile = ntile(terr_GDE_dens, 100),
    aqua_GDE_Ptile = ntile(aqua_GDE_dens, 100)) |> 
  group_by(terr_GDE_Ptile, aqua_GDE_Ptile) |> 
  summarise(
    wtd_landform =  wtd_mode(landform, area),
    wtd_aridity  = wtd.quantile(aridity,  weights = area, probs = 0.5),
    wtd_seasonality  = wtd.quantile(seasonality,  weights = area, probs = 0.5),
    total_area = sum(area, na.rm = T),
    .groups = "drop"
  )
coll_df$wtd_landform = as.character(coll_df$wtd_landform)

# now plow these ptile X ptile TILE plots for ARIDITY (and landform)
ggplot(coll_df, aes(x = aqua_GDE_Ptile, 
                    y = terr_GDE_Ptile, 
                    fill = wtd_aridity)) +
  geom_tile(color = "transparent", linewidth = 0) +
  scale_fill_gradientn(
    colours = met.brewer("Veronese", n = 100, direction = 1),
    limits  = c(0, 2),
    breaks  = c(0, 0.03, 0.2, 0.5, 1, 2),
    labels  = c("0", "0.03", "0.2", "0.5", "1", "2"),
    values  = c(0, 0.015, 0.1, 0.5, 0.75, 1),
    oob     = scales::squish
  ) +
  geom_abline(slope = 1, intercept = 0, colour = "red", linewidth = 1, linetype = "dashed") +
  scale_x_continuous(breaks = seq(0,100)) +
  scale_y_continuous(breaks = seq(0,100)) +
  # scale_y_continuous(breaks = c(1,2,4)) +
  labs(x = "aquatic GDE ptile",
       y = "terrestrial GDE ptile") +
  coord_equal() +
  coord_cartesian(expand = 0) +
  theme_void() +
  theme(legend.position = "none")
ggsave(file = here("plots/GDE_terr_VS_aquatic_aridity.png"), 
       plot = last_plot(), device = "png", 
       width = 51.8, height = 34.6, units = "mm", dpi = 400)


# now plow these ptile X ptile TILE plots for (aridity and) LANDFORM
ggplot(coll_df, aes(x = aqua_GDE_Ptile, 
                    y = terr_GDE_Ptile, 
                    fill = wtd_landform)) +
  geom_tile(color = "transparent") +
  scale_fill_manual(
    values = c(
      "1" = "#234657",  # mountains
      "2" = "#6C8F8E",  # hills
      "3" = "#FA6D00",  # plateaus
      "4" = "#ADAE37"   # plains
    ),
    labels = c("1" = "mountains", "2" = "hills", "3" = "plateaus", "4" = "plains")
  ) +
  geom_abline(slope = 1, intercept = 0, colour = "red", linewidth = 1, linetype = "dashed") +
  scale_x_continuous(breaks = seq(0,100)) +
  scale_y_continuous(breaks = seq(0,100)) +
  # scale_y_continuous(breaks = c(1,2,4)) +
  labs(x = "aquatic GDE ptile",
       y = "terrestrial GDE ptile") +
  coord_equal() +
  coord_cartesian(expand = 0) +
  theme_void() +
  theme(legend.position = "none")
ggsave(file = here("plots/GDE_terr_VS_aquatic_landform.png"), 
       plot = last_plot(), device = "png", 
       width = 51.8, height = 34.6, units = "mm", dpi = 400)

# now plow these ptile X ptile TILE plots for SEASONALITY
ggplot(coll_df, aes(x = aqua_GDE_Ptile, 
                    y = terr_GDE_Ptile, 
                    fill = wtd_seasonality)) +
  geom_tile(color = "transparent", linewidth = 0) +
  scale_fill_gradientn(
    colours = met.brewer("Pissaro", n = 100)[50:100],
    limits  = c(0.4, 1.2),
    oob     = scales::squish
  ) + 
  geom_abline(slope = 1, intercept = 0, colour = "red", linewidth = 1, linetype = "dashed") +
  scale_x_continuous(breaks = seq(0,100)) +
  scale_y_continuous(breaks = seq(0,100)) +
  # scale_y_continuous(breaks = c(1,2,4)) +
  labs(x = "aquatic GDE ptile",
       y = "terrestrial GDE ptile") +
  coord_equal() +
  coord_cartesian(expand = 0) #+
  theme_void() +
  theme(legend.position = "none")
ggsave(file = here("plots/GDE_terr_VS_aquatic_seasonality.png"),
       plot = last_plot(), device = "png",
       width = 51.8, height = 34.6, units = "mm", dpi = 400)
  