### ---------------------// 
# Script objective:
# Evaluate how GDEs are protected 
### ---------------------//
library(here); source(here(("scripts/on_button.R")))
###

# import GDE area estimate
hu23 = terra::rast("D:/Geodatabase/GDEs/Huggins_2023/gde-map.tif")
hu23[hu23>1] = 1

# import area raster, necessary for protected percentages
area_1km = terra::rast("C:/Users/xande/Documents/1.projects-scripts/gde-global-comparison/data/area_raster_30arcsec.tif")

# multiply by Link et al. 2023 GDE land cover 
li23 = terra::rasterize(x = terra::vect("D:/Geodatabase/GDEs/Link_2023/GDE_potentials.shp"),
                        y = area_1km, # this is the base resolution of Link 2023
                        field = "I5", # share of relevant land covers
                        touches = TRUE)

# to calculate GDE area, multiply hu23, li23 area cover, and area raster
gde_area = hu23 * area_1km * li23

gde_area = terra::writeRaster(gde_area, here("data/GDE_area_hugginsXlink_1km.tif"))
gde_area = terra::rast(here("data/GDE_area_hugginsXlink_1km.tif"))

## aggregate to 5 arcmin for easier national computation
gde_area_30am = terra::aggregate(x = gde_area, 
                                fact = 60,
                                fun = "sum", na.rm = T)
names(gde_area_30am) = "gde_area"

# import protected areas estimate
pas = terra::rast("C:/Users/xande/Documents/1.projects-scripts/global-gde-map/dryland-gde-map/data/wdpa_binary_1km.tif")
pas[is.na(pas)] = 0

# calculate protected area (not just binary)
pa_area = pas * area_1km

pa_area_30am = terra::aggregate(x = pa_area, 
                               fact = 60,
                               fun = "sum", na.rm = T)
names(pa_area_30am) = "prot_area"

# resample the GDE data to 1km for comparison with PAs
gde_prot_1km = gde_area * pas


gde_area_in_pa_30am = terra::aggregate(x = gde_prot_1km, 
                                      fact = 60,
                                      fun = "sum", na.rm = T)
names(gde_area_in_pa_30am) = "prot_gde_area"


# import an aridity gradient 
aridity = terra::rast("D:/Geodatabase/Climate/Aridity/ai_et0/ai_et0.tif") / 10000
aridity_30m = terra::resample(x = aridity, y = gde_area_in_pa_30am, method = "bilinear")

gde_area_30am[is.na(gde_area_30am)] = 0
gde_area_in_pa_30am[is.na(gde_area_in_pa_30am)] = 0

gde_density = (100 * gde_area_30am) / WGS84_areaRaster(0.5) |> rast()
names(gde_density) = "gde_density"


# create summary tibble
summary_df = c(gde_area_30am, gde_area_in_pa_30am, gde_density, aridity_30m, WGS84_areaRaster(0.5) |> rast()) |> 
  as_tibble() |> 
  set_colnames(c('gde_area', 'prot_gde_area', 'gde_density', 'aridity', 'area')) |> 
  mutate(GDE_dens_bin = cut(gde_density,
                        breaks = seq(0, 100, by = 1),
                        include.lowest = TRUE,
                        right = FALSE),
         arid_bin = cut(aridity,
                        breaks = (10^seq(
                          log10(0.01),
                          log10(10),
                          length.out = 100 + 1   # +1 because breaks define bin edges
                        )),
                        include.lowest = TRUE,
                        right = FALSE)) |> 
  group_by(GDE_dens_bin, arid_bin) |> 
  summarise(
    n = n(),
    area_sum = sum(area, na.rm = TRUE),
    gde_area_sum = sum(gde_area, na.rm = T),
    gde_prot_sum = sum(prot_gde_area, na.rm = T),
    .groups = "drop"
  )  |> 
  mutate(
    gde_protection = (100 * gde_prot_sum) / gde_area_sum
  )

plot_df = summary_df |> 
  mutate(
    arid_lo  = as.numeric(str_match(as.character(arid_bin),  "\\[([^,]+),")[,2]),
    arid_hi  = as.numeric(str_match(as.character(arid_bin),  ",([^\\)\\]]+)")[,2]),
    gde_d_lo = as.numeric(str_match(as.character(GDE_dens_bin), "\\[([^,]+),")[,2]),
    gde_d_hi = as.numeric(str_match(as.character(GDE_dens_bin), ",([^\\)\\]]+)")[,2])
  )

ggplot(plot_df) +
  geom_hline(yintercept = seq(0, 100, 10), 
             linewidth = 0.5, lty = "dashed", col=  "grey") +
  geom_vline(xintercept = c(0.03, 0.2, 0.5, 1, 5), 
             linewidth = 0.5, lty = "dashed", col=  "grey") +
  
  geom_rect(aes(
    xmin = arid_lo, xmax = arid_hi,
    ymin = gde_d_lo, ymax = gde_d_hi,
    fill = gde_protection
  )) +
  scale_x_log10(limits = c(0.02, 10)) +
  scale_fill_stepsn(
    colours = c("#5C1A0B", "#D4956A", "#2D7A6B"),
    breaks  = c(10, 30),
    limits  = c(0, 100),
    values  = scales::rescale(c(0, 10, 30, 100)),
    oob     = scales::squish,
    name    = "GDE protection"
  ) + 
  coord_cartesian(expand = FALSE, xlim = c(0.01, 10), ylim = c(0, 100)) +
  labs(x = "aridity index", 
       y = "GDE area density") +
 
  
  theme_void() + theme(legend.position = "none")


ggsave(file = here("plots/protection_pct_vs_gde_dens_x_aridity.png"), 
       plot = last_plot(), device = "png", 
       width = 51.8*2, height = 34.6*2, units = "mm", dpi = 400)


## Plot rolling average over aridity and GDE area density
library(zoo)

plot_df$arid_mid = (plot_df$arid_hi + plot_df$arid_lo)/2 

rolling_arid = plot_df |> 
  group_by(arid_mid) |> 
  summarise(
    gde_prot_sum = sum(gde_prot_sum, na.rm = T),
    gde_sum = sum(gde_area_sum, na.rm = T)
  ) |> 
  mutate(
    gde_protection = (100  *gde_prot_sum)/gde_sum
  ) |> 
  arrange(arid_mid) |>
  mutate(
    gde_prot_roll = rollmean(gde_protection, k = 7, fill = NA, align = "center")
  )

ggplot(rolling_arid, aes(x = arid_mid)) +
 
  # background bands
  annotate("rect", xmin = 0.01, xmax = 10, ymin =  0, ymax = 10,
           fill = "#5C1A0B", alpha = 0.5) +
  annotate("rect", xmin = 0.01, xmax = 10, ymin = 10, ymax = 30,
           fill = "#D4956A", alpha = 0.5) +
  annotate("rect", xmin = 0.01, xmax = 10, ymin = 30, ymax = 100,
           fill = "#2D7A6B", alpha = 0.5) +
  
  geom_line(aes(y = gde_protection), colour = "grey80", linewidth = 0.4) +  # raw
  geom_line(aes(y = gde_prot_roll),  colour = "black", linewidth = 1) +   # smoothed
  scale_x_log10(limits = c(0.01, 10)) +
  coord_cartesian(ylim = c(0, 50), expand = c(0,0)) +
  labs(x = "Aridity index", y = "GDE protection (%)") +
  geom_vline(xintercept = c(0.03, 0.2, 0.5, 1, 5), 
             linewidth = 0.5, lty = "dashed", col=  "grey") +
  theme_void()

ggsave(file = here("plots/protection_pct_ROLLING_aridity.png"), 
       plot = last_plot(), device = "png", 
       width = 51.5*2, height = 11.5*2, units = "mm", dpi = 400)



########################
## now repeat for GDE area density


plot_df$gde_mid = (plot_df$gde_d_lo + plot_df$gde_d_hi)/2 

rolling_arid = plot_df |> 
  group_by(gde_mid) |> 
  summarise(
    gde_prot_sum = sum(gde_prot_sum, na.rm = T),
    gde_sum = sum(gde_area_sum, na.rm = T)
  ) |> 
  mutate(
    gde_protection = (100  *gde_prot_sum)/gde_sum
  ) |> 
  arrange(gde_mid) |>
  mutate(
    gde_prot_roll = rollmean(gde_protection, k = 7, fill = NA, align = "center")
  )

ggplot(rolling_arid, aes(x = gde_mid)) +
  
  # background bands
  annotate("rect", xmin = 0, xmax = 100, ymin =  0, ymax = 10,
           fill = "#5C1A0B", alpha = 0.5) +
  annotate("rect", xmin = 0, xmax = 100, ymin = 10, ymax = 30,
           fill = "#D4956A", alpha = 0.5) +
  annotate("rect", xmin = 0, xmax = 100, ymin = 30, ymax = 100,
           fill = "#2D7A6B", alpha = 0.5) +
  
  geom_line(aes(y = gde_protection), colour = "grey80", linewidth = 0.4) +  # raw
  geom_line(aes(y = gde_prot_roll),  colour = "black", linewidth = 1) +   # smoothed
  scale_x_continuous(limits = c(0, 100)) +
  coord_cartesian(ylim = c(0, 50), expand = c(0,0)) +
  labs(x = "Aridity index", y = "GDE protection (%)") +
  geom_vline(xintercept = seq(0, 100, 10), 
             linewidth = 0.5, lty = "dashed", col=  "grey") +
  theme_void()

ggsave(file = here("plots/protection_pct_ROLLING_gde_density.png"), 
       plot = last_plot(), device = "png", 
       width = 34.7*2, height = 11.5*2, units = "mm", dpi = 400)