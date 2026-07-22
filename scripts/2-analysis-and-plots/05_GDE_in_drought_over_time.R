### ---------------------\\ 
# Script objective:
# Evaluate fraction of global GDEs that are under moderate drought from 1970 to 2025
### ---------------------\\
library(here); source(here(("scripts/on_button.R"))); library(terra); library(readr)
###

area = WGS84_areaRaster(5/60) |> rast()
names(area) = "area"

template_0d25 = terra::rast(extent = c(-180, 180, -90, 90),
                            resolution = 0.25)

template_0d5 = terra::rast(extent = c(-180, 180, -90, 90),
                            resolution = 0.5)

R.gde_types_hu = terra::rast(here("data/GDE_types_hu_for_comparison.tif"))

R.gde_terr = R.gde_types_hu$terr_GDE_area
R.gde_aqua = R.gde_types_hu$aqua_GDE_area

sum.gde_terr = sum(R.gde_terr[], na.rm = T)
sum.gde_aqua = sum(R.gde_aqua[], na.rm = T)


frac_tibb = expand.grid(
  year = seq(1970, 2025),
  month = seq(1, 12)) |> 
  as_tibble() |> 
  arrange(year, month) |> 
  mutate(month = sprintf("%02d", month))
frac_tibb$terr_frac = rep(NA)
frac_tibb$aqua_frac = rep(NA)


for (ii in 1:nrow(frac_tibb)) {
  # ii = 1
  
  # import the SPEI12 data
  
  iter_SPEI = terra::rast(paste0("D:/Geodatabase/Climate/SPEI_ERA5/",
                                 "SPEI12_genlogistic_global_era5_moda_ref1991to2020_",
                                 as.character(frac_tibb$year[ii]), 
                                 frac_tibb$month[ii], ".nc"))
  
  message(paste0("D:/Geodatabase/Climate/SPEI_ERA5/",
                 "SPEI12_genlogistic_global_era5_moda_ref1991to2020_",
                 as.character(frac_tibb$year[ii]), 
                 frac_tibb$month[ii], ".nc"))

  iter_SPEI = terra::resample(x = iter_SPEI, y = template_0d25,
                              method = "near")
  iter_SPEI_0d5 = terra::aggregate(x = iter_SPEI, fact = 2, method = "mean", na.rm = T)
  
  BIN_drought = rast(iter_SPEI_0d5)
  BIN_drought[]= NA
  BIN_drought[iter_SPEI_0d5 < -1] = 1
  
  temp_GDE_terr_drought = R.gde_terr * BIN_drought
  temp_GDE_aqua_drought = R.gde_aqua * BIN_drought
  
  frac_tibb$terr_frac[ii] = sum(temp_GDE_terr_drought[], na.rm = T) / sum.gde_terr
  frac_tibb$aqua_frac[ii] = sum(temp_GDE_aqua_drought[], na.rm = T) / sum.gde_aqua
  
  message(ii, " done, which is ", round(100*ii/nrow(frac_tibb), 2), "%")
  
}

frac_tibb_backup = frac_tibb

frac_tibb$monthID = seq(1, nrow(frac_tibb))

# find p20 and p80 range of aquatic and terrestrial 
# areas in drought over 1970-2000
p20 = frac_tibb |> 
  filter(year<2000) |> 
  pull(aqua_frac) |> 
  quantile(0.05)

p80 = frac_tibb |> 
  filter(year<2000) |> 
  pull(aqua_frac) |> 
  quantile(0.95)

ggplot(frac_tibb, aes(x = monthID)) +
  geom_hline(yintercept = seq(0, 0.5, 0.1), color = "grey80", linewidth = 0.5) + 
  geom_vline(xintercept = seq(1, 672, 120), color = "grey80", linewidth = 0.5) + 
  geom_ribbon(aes(ymin = p20, ymax = p80), 
              fill = "#91bfdb", alpha = 0.3) +
  geom_line(aes(y = aqua_frac),
            color = "#2166ac", linewidth = 0.5) +
  geom_line(aes(y = terr_frac ), 
            color = "#779B46", linewidth = 0.5) +
  geom_hline(yintercept = c(p20, p80), color = "black", linewidth = 0.6) + 
  coord_cartesian(xlim = c(1, 672), ylim = c(0, 0.5), expand = 0) +
  theme_void()

ggsave(file = here("plots/GDE_in_drought_ts.png"),
       plot = last_plot(), device = "png",
       width = 54.9, height = 36.6, units = "mm", dpi = 400)