### ---------------------// 
# Script objective:
# Plot soil organic carbon vs GDE density
### ---------------------//
library(here); source(here(("scripts/on_button.R")))
###

area_r = WGS84_areaRaster(5/60) |> rast()

soc_0_5 = terra::rast("D:/Geodatabase/Carbon/soc_0-5cm_mean_5000.tif") * 5
soc_5_15 = terra::rast("D:/Geodatabase/Carbon/soc_5-15cm_mean_5000.tif") * 10
soc_15_30 = terra::rast("D:/Geodatabase/Carbon/soc_15-30cm_mean_5000.tif") * 15
soc_30_60 = terra::rast("D:/Geodatabase/Carbon/soc_30-60cm_mean_5000.tif") * 30
soc_60_100 = terra::rast("D:/Geodatabase/Carbon/soc_60-100cm_mean_5000.tif") * 40
soc_100_200 = terra::rast("D:/Geodatabase/Carbon/soc_100-200cm_mean_5000.tif") * 100


soc = c(soc_0_5, soc_5_15, soc_15_30, soc_30_60, soc_60_100, soc_100_200) |> sum()
soc = soc/200

soc_wgs = terra::project(x = soc, y = area_r, method = "bilinear")
summary(soc_wgs[])

soc_classified = classify(
  soc_wgs,
  rcl = matrix(c(
    0,   50,  1,
    50,  100,  2,
    100,  200, 3,
    200, 400, 4,
    400, Inf, 5
  ), ncol = 3, byrow = TRUE),
  include.lowest = TRUE
)

# import gde area density
hu23 = terra::rast("D:/Geodatabase/GDEs/Huggins_2023/hu23_areadens_5m.tif")

# multiply by Link area suitability
li23 = terra::rasterize(x = terra::vect("D:/Geodatabase/GDEs/Link_2023/GDE_potentials.shp"),
                        y = area_r, # this is the base resolution of Link 2023
                        field = "I5", # share of relevant land covers
                        touches = TRUE)

hu23 = hu23*li23
hu23_area = hu23 * WGS84_areaRaster(5/60) |> rast()
names(hu23_area) = "gde_area"

hu23_classified = classify(
  hu23,
  rcl = matrix(c(
    0,   0.1,  1,
    0.1,  0.3,  2,
    0.3,  0.5, 3,
    0.5, 0.7, 4,
    0.7, Inf, 5
  ), ncol = 3, byrow = TRUE),
  include.lowest = TRUE
)


gde_x_soc = hu23_classified + (10*soc_classified)


## make map
outline = terra::vect("C:/Users/xande/Documents/1.projects-scripts/sustainability-puzzles/data/land_mask_polygon.sqlite") |> 
  st_as_sf()

map =
  tm_shape(outline, crs = "+proj=robin") +
  tm_fill(col = "grey", lwd = 0, border.col = NA) +
  tm_shape(terra::project(x = gde_x_soc, y = "+proj=robin", method = "near")) +
  tm_raster(
    style   = "cat",
    palette = c(
      "11" = "#C5C2AA",
      "12" = "#7EB0CD",
      "13" = "#54A7D7",
      "14" = "#2A9EE1",
      "15" = "#0096EB",
      "21" = "#D3CB83",
      "22" = "#9BA894",
      "23" = "#808A8E",
      "24" = "#507AAB",
      "25" = "#2B64B7",
      "31" = "#E2D45C",
      "32" = "#AEA36F",
      "33" = "#8C827B",
      "34" = "#696288",
      "35" = "#484294",
      "41" = "#E0C73C",
      "42" = "#C19E49",
      "43" = "#A27457",
      "44" = "#834A64",
      "45" = "#652172",
      "51" = "#F1CC16",
      "52" = "#D59924",
      "53" = "#AB4C3A",
      "54" = "#9D3341",
      "55" = "#820050"
    )
  ) +
  tm_shape(outline) + 
  tm_borders(col = "black", lwd = 1) +
  tm_layout(legend.show = F, legend.frame = F, frame = F, 
            bg = FALSE, bg.color = "transparent")
map


tmap_save(map, here("plots/gMAP_SOC_X_GWDensity.pdf"), dpi = 400, units = "in")



### create empirical CDF of global SOC across GDE density bins
library(zoo)

ecdf_df = c(soc_wgs, 100*hu23, hu23_area) |> 
  as_tibble() |> 
  set_colnames(c('soc', 'gde_dens', 'gde_area')) |> 
  mutate(GDE_dens_bin = cut(gde_dens,
                            breaks = seq(0, 100, by = 1),
                            include.lowest = TRUE,
                            right = FALSE)) |> 
  group_by(GDE_dens_bin) |> 
  summarise(
    local_soc = sum(soc, na.rm = T),
    local_gde_area = sum(gde_area, na.rm = T)
  ) |> 
  drop_na() |> 
  arrange(GDE_dens_bin) |> 
  mutate(
    GDE_bin_mid = seq(0.5, 99.5, by = 1 ),
    soc_rel_gde = local_soc / local_gde_area
  ) |> 
  mutate(
    soc_rel_gde_roll = rollmean(soc_rel_gde, k = 7, fill = NA, align = "center")
  )

ecdf_df$cumsum = cumsum(ecdf_df$local_soc) / sum(ecdf_df$local_soc)


ggplot(ecdf_df) +
  
  # background bands
  # annotate("rect", xmin = 0.01, xmax = 10, ymin =  0, ymax = 10,
  #          fill = "#5C1A0B", alpha = 0.5) +
  # annotate("rect", xmin = 0.01, xmax = 10, ymin = 10, ymax = 30,
  #          fill = "#D4956A", alpha = 0.5) +
  # annotate("rect", xmin = 0.01, xmax = 10, ymin = 30, ymax = 100,
  #          fill = "#2D7A6B", alpha = 0.5) +
  
  geom_line(aes(x = GDE_bin_mid, y = 100*cumsum),  
            colour = "black", linewidth = 2) +  
  
  # geom_line(aes(x = GDE_bin_mid, y = soc_rel_gde_roll),  
  #           colour = "red", linewidth = 1) +   # smoothed
  
  scale_x_continuous(limits = c(0, 100)) +
  scale_y_continuous(limits = c(0,100)) +
  coord_cartesian(expand = 0) +
  # coord_cartesian(ylim = c(0, 50), expand = c(0,0)) +
  labs(x = "GDE area density", y = "SOC cumulative sum") +
  geom_vline(xintercept = seq(0, 100, 20), 
             linewidth = 0.5, lty = "dashed", col=  "grey") +
  geom_hline(yintercept = seq(0, 100, 25), 
             linewidth = 0.5, lty = "dashed", col=  "grey") +
  theme_void()


ggsave(file = here("plots/ECDF_soc_x_gde_dens.png"), 
       plot = last_plot(), device = "png", 
       width = 25.2*2.5, height = 19.2*2.5, units = "mm", dpi = 400)
