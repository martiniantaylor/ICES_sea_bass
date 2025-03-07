
################################################################################
# Making Maps script
# 
# This script loads lat lons and ices rectangles for feeding spawning and southern samples.
#It produces a feeding map and spawning map with inlay of southern samples.

#It also makes the supp info Sample tables
#
# Author: Martin Taylor (martin.taylor@uea.ac.uk)
################################################################################


# Load necessary libraries
library(rnaturalearth)  # For world map data
library(sf)             # For handling spatial data
library(ggplot2)        # For plotting
library(dplyr)         # For data manipulation
library(tidyr)         # For reshaping data
library(patchwork)     # For combining plots

# Import sample lat/lon datasets
feeding_10 <- read.csv("data/metadata/bass_samples_feeding_10.csv", stringsAsFactors = FALSE)
spawning_10 <- read.csv("data/metadata/bass_samples_spawning_10.csv", stringsAsFactors = FALSE)
southern <- read.csv("data/metadata/bass_samples_southern.csv", stringsAsFactors = FALSE)

# Import ICES area shape file
ices_areas <- read_sf("data/map_shape_files/ICES_Areas_20160601_cut_dense_3857.shp")    
ices_areas2 <- st_simplify(ices_areas, dTolerance = 75)  # Simplify geometry for better performance

# Import Europe shape file
Europe <- ne_countries(continent = "Europe", returnclass = "sf", scale = "large")
Europe2 <- st_simplify(Europe, dTolerance = 75)  # Simplify geometry for better performance

# Create ICES area labels (excluding some subdivisions)
ices_areas_label <- ices_areas %>% 
  unite("sub_division", c("SubArea", "Division"), sep = ".") %>% 
  mutate(sub_division = if_else(sub_division %in% c("5.b", "3.c", "7.j", "2.b", "2.a", "1.a", "1.b", "3.d", "14.a", "14.b"), "", sub_division))

# --- PLOTTING ---

# Feeding samples map
feeding_map <- ggplot(Europe2) +
  geom_sf(data = ices_areas2, color = "#999999", fill = "white") +
  geom_sf(color = "#999999", fill = "#999999") +
  geom_sf_text(data = ices_areas_label, aes(label = sub_division)) +
  coord_sf(xlim = c(-11, 11), ylim = c(49, 63)) +
  geom_point(data = feeding_10, aes(x = lon, y = lat), size = 7, colour = "#FFD700") +
  geom_text(data = feeding_10, aes(label = ices_rect, x = lon, y = lat), size = 2) +
  theme_void() +
  theme(panel.background = element_rect(fill = "#999999", colour = "#999999"),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  labs(title = "Feeding Samples")

# Spawning samples map
spawning_map <- ggplot(Europe2) +
  geom_sf(data = ices_areas2, color = "#999999", fill = "white") +
  geom_sf(color = "#999999", fill = "#999999") +
  geom_sf_text(data = ices_areas_label, aes(label = sub_division)) +
  coord_sf(xlim = c(-11, 11), ylim = c(49, 63)) +
  geom_point(data = spawning_10, aes(x = lon, y = lat), size = 7, colour = "#FFD700") +
  geom_text(data = spawning_10, aes(label = ices_rect, x = lon, y = lat), size = 2) +
  theme_void() +
  theme(panel.background = element_rect(fill = "#999999", colour = "#999999"),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  labs(title = "Spawning Samples")

# Southern sites map (for inset)
southern_map <- ggplot(Europe2) +
  geom_sf(color = "#999999", fill = "#999999") +
  coord_sf(xlim = c(-11, 10), ylim = c(35, 45), expand = FALSE) +
  geom_point(data = southern, aes(x = lon, y = lat), size = 7, colour = "#FFD700") +
  geom_text(data = southern, aes(label = ices_rect, x = lon, y = lat), size = 2) +
  theme_void() +
  theme(panel.background = element_rect(fill = "white", colour = "#505050"),
        panel.border = element_rect(colour = "black", fill = NA, size = 1))

# Use gginlay to manually adjust the inset placement (interactive)
# gginlay(spawning_map, southern_map)

# Combine main map with inset
combined_spawning_map <- spawning_map +
  patchwork::inset_element(southern_map, 0.6, 0.78, 1, 1, align_to = 'panel', on_top = TRUE, ignore_tag = TRUE)

#save plots to file

#spawning map
pdf(file="output/maps/spawning_map.pdf")
spawning_map
dev.off()

#feeding map
pdf(file="output/maps/feeding_map.pdf")
feeding_map
dev.off()

#combined spawning and southern map
pdf(file="output/maps/spawning_map_southern.pdf")
combined_spawning_map
dev.off()

