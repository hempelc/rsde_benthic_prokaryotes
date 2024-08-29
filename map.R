library("ggplot2")
library("sf")
library("rnaturalearth")
library("rnaturalearthdata")
library(readxl)
library(dplyr)
library(leaflet)


metadata_file <- "/Users/simplexdna/Desktop/RSDE BAC soil project/RSDE metabarcoding metadata.xlsx"
metadata <- as.data.frame(read_excel(metadata_file))
plot_outdir <- "/Users/simplexdna/Desktop/RSDE BAC soil project/analysis/plots/fancy_plots/"

metadata$Latitude <- as.numeric(metadata$Latitude)
metadata$Longitude <- as.numeric(metadata$Longitude)
metadata$'Red Sea zone' <- factor(metadata$'Red Sea zone', levels=c("Gulf of Aqaba", "North Red Sea", "Central-north Red Sea", "Central-south Red Sea", "South Red Sea"))
metadata <- metadata %>%
  mutate(
    "DSC" = case_when(
      Phase == "DSC" ~ "Discovery site",
      Phase != "DSC" ~ "Not discovery site",
      TRUE ~ NA_character_  # default condition (optional)
    )
  )

theme_set(theme_bw())


world <- ne_countries(scale = "medium", returnclass = "sf")
coords <- metadata[, c("Latitude", "Longitude")]

# By Zone
map_zone <- ggplot(data = world) +
  geom_sf() +
  geom_point(data = metadata, aes(x = Longitude, y = Latitude, shape = DSC, color = `Red Sea zone`), size = 2) +
  scale_shape_manual(values = c(4, 20)) +  # 1 for circle, 3 for cross
  coord_sf(xlim = c(44, 32), ylim = c(30, 12), expand = FALSE) +
  theme_bw()
ggsave(units="px", width=3000, height=3000, dpi=600, file.path(plot_outdir, "map_zone.png"), plot = map_zone)
ggsave(units="px", width=3000, height=3000, dpi=600, file.path(plot_outdir, "map_zone.svg"), plot = map_zone)

# By depth
map_depth <- ggplot(data = world) +
  geom_sf() +
  geom_point(data = metadata, aes(x = Longitude, y = Latitude, shape = DSC, color = `Depth`), size = 2) +
  scale_shape_manual(values = c(4, 20)) +  # 1 for circle, 3 for cross
  coord_sf(xlim = c(44, 32), ylim = c(30, 12), expand = FALSE) +
  theme_bw() +
  scale_colour_viridis_c(direction=-1)
ggsave(units="px", width=3000, height=3000, dpi=600, file.path(plot_outdir, "map_depth.png"), plot = map_depth)
ggsave(units="px", width=3000, height=3000, dpi=600, file.path(plot_outdir, "map_depth.svg"), plot = map_depth)

# Get background
leaflet() %>%
  # Set the initial view (center and zoom level)
  addTiles() %>% 
  addProviderTiles("Esri.WorldImagery") %>%
  fitBounds(44, 32, 30, 12)
