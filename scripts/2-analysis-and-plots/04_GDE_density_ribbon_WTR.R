### ---------------------\\ 
# Script objective:
# Plot distribution of GDE area density against water table ratio 
### ---------------------\\
library(here); source(here(("scripts/on_button.R")))
###

area = WGS84_areaRaster(5/60) |> rast()
names(area) = "area"

gde_dens = terra::rast(here("data/GDE_types_hu_for_comparison.tif")) / 
  WGS84_areaRaster(0.5) |> rast()
gde_dens = terra::resample(x = gde_dens$terr_GDE_area, y = area, "near")

wtr = terra::rast("D:/Geodatabase/Groundwater/GRT_WTR/LOG_WTR_L_01_5arcmin.tif")
wtr = terra::resample(x = wtr, y = area, "near")

# need to mask-out arid regions 
recharge = terra::rast("D:/Geodatabase/Groundwater/Doell_Recharge2008/r1_doll2008.tif")
recharge = terra::resample(recharge, wtr)
wtr[recharge < 5] = NA

stack_df = c(gde_dens, wtr, area) |> 
  as_tibble() |> 
  set_colnames(c("GDE_frac", "wtr", "area")) |> 
  drop_na()

wtr_gde_df = stack_df |>
  mutate(
    wtr_bin = cut(wtr, 
                  breaks = seq(-5, 5, length.out = 101),
                  include.lowest = TRUE,
                  labels = seq(-5, 5, length.out = 101)[1:100] + 0.05),
    gde_bin = cut(GDE_frac,
                  breaks = seq(0, 1, length.out = 101),
                  include.lowest = TRUE,
                  labels = seq(0, 1, length.out = 101)[1:100] + 0.005)
  ) |>
  mutate(wtr_bin = as.numeric(as.character(wtr_bin)),
         gde_bin = as.numeric(as.character(gde_bin))) |>
  group_by(wtr_bin, gde_bin) |>
  summarise(total_area = sum(area, na.rm = TRUE),
            .groups = "drop")

# Plot
ggplot(wtr_gde_df, aes(x = wtr_bin, y = gde_bin, fill = total_area)) +
  geom_tile() +
  scale_fill_viridis_c(
    trans = "log10",
    na.value = "grey90",
    name = "Total area"
  ) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(
    x = "Water table ratio (WTR)",
    y = "GDE fraction"
  ) +
  theme_minimal()


# Percentile ribbon across WTR bins
ribbon_df = stack_df |>
  mutate(wtr_bin = cut(wtr,
                       breaks = c(-Inf, seq(-3, 3, length.out = 51), Inf),
                       include.lowest = TRUE,
                       labels = c(-3.1, seq(-3, 3, length.out = 50), 3.1))) |>
  mutate(wtr_bin = as.numeric(as.character(wtr_bin))) |>
  group_by(wtr_bin) |>
  summarise(
    p05 = Hmisc::wtd.quantile(x = GDE_frac, weights = area, probs = 0.05, na.rm = T),
    p25 = Hmisc::wtd.quantile(x = GDE_frac, weights = area, probs = 0.25, na.rm = T),
    p50 = Hmisc::wtd.quantile(x = GDE_frac, weights = area, probs = 0.50, na.rm = T),
    p75 = Hmisc::wtd.quantile(x = GDE_frac, weights = area, probs = 0.75, na.rm = T),
    p95 = Hmisc::wtd.quantile(x = GDE_frac, weights = area, probs = 0.95, na.rm = T),
    .groups = "drop"
  )

ggplot(ribbon_df, aes(x = wtr_bin)) +
  geom_hline(yintercept = seq(0, 1, 0.25), color = "grey80", linewidth = 0.5) + 
  geom_vline(xintercept = seq(-2, 2, 1), color = "grey80", linewidth = 0.5) + 
  geom_ribbon(aes(ymin = p05, ymax = p95), 
              fill = "#91bfdb", alpha = 0.3) +
  geom_ribbon(aes(ymin = p25, ymax = p75), 
              fill = "#91bfdb", alpha = 0.5) +
  geom_line(aes(y = p50), 
            color = "#2166ac", linewidth = 0.8) +
  coord_cartesian(xlim = c(-2, 2), ylim = c(0, 1), expand = 0) +
  theme_void()
ggsave(file = here("plots/WTR_vs_TERR_GDE_dens.png"),
       plot = last_plot(), device = "png",
       width = 51.8, height = 34.6, units = "mm", dpi = 400)


######################################
## Map areas that have >50% GDE area frac
## but have negative WTR
######################################

hu23 = terra::rast("D:/Geodatabase/GDEs/Huggins_2023/hu23_areadens_5m.tif")

flag_r = rast(wtr)
flag_r[] = NA
flag_r[wtr < -0.5 & hu23 > 0.50] = 1
flag_r = as.factor(flag_r)

# make the map
outline = terra::vect("C:/Users/xande/Documents/1.projects-scripts/sustainability-puzzles/data/land_mask_polygon.sqlite") |> 
  st_as_sf()
map =  
  tm_shape(outline, crs = "+proj=robin") +
  tm_fill(col = "grey", lwd = 0, border.col = NA) +
  tm_shape(terra::project(x = flag_r, y = "+proj=robin", method = "near")) +
  tm_raster(
    col.scale = tm_scale_categorical(
      values = "red")
    ) +
  tm_shape(outline) + 
  tm_borders(col = "black", lwd = 1) +
  tm_layout(legend.show = F, legend.frame = F, frame = F, bg = FALSE, bg.color = "transparent")
map
tmap_save(map, here("plots/gMAP_GDE_regions_WTR_negative.pdf"), dpi = 400, units = "in")