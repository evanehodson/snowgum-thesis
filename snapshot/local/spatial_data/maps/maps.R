library(tidyverse)
library(sf)
library(rnaturalearth)
library(ggspatial)

# ---- Setup paths ----
base_dir <- "C:/Users/evane/OneDrive/Desktop/Masters Dissertation/genomic_data/land_gen"
ala_dir  <- "C:/Users/evane/OneDrive/Desktop/Masters Dissertation/spatial_data/ALA_occurences"

out_dir  <- file.path(base_dir, "figures")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Set global font family
font_family <- "Times New Roman"

# ---- Step 1: Explicit CSV File Paths (All 6 Subspecies) ----
subspecies_files <- c(
  file.path(ala_dir, "acer_range", "acer_range.csv"),
  file.path(ala_dir, "deb_range", "deb_range.csv"),
  file.path(ala_dir, "hed_range", "hed_range.csv"),
  file.path(ala_dir, "pavi_range", "pavi_range.csv"),
  file.path(ala_dir, "pauci_range", "pauci_range.csv"),
  file.path(ala_dir, "niph_range", "niph_range.csv")
)

existing_files <- subspecies_files[file.exists(subspecies_files)]

read_and_clean_csv <- function(fp) {
  read_csv(fp, col_types = cols(.default = col_character()), show_col_types = FALSE) |>
    mutate(
      decimalLatitude  = as.numeric(decimalLatitude),
      decimalLongitude = as.numeric(decimalLongitude)
    ) |>
    filter(
      is.finite(decimalLatitude), is.finite(decimalLongitude),
      decimalLongitude >= 129.0, decimalLongitude <= 155.0,  # Filter out WA (<129°E) and NZ (>155°E)
      decimalLatitude  >= -45.0, decimalLatitude  <= -10.0
    ) |>
    select(decimalLatitude, decimalLongitude)
}

all_occ_clean <- purrr::map_dfr(existing_files, read_and_clean_csv) |>
  distinct(decimalLatitude, decimalLongitude, .keep_all = TRUE) |>
  st_as_sf(coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)

message(sprintf("Loaded %d unique occurrence points across species range.", nrow(all_occ_clean)))

# ---- Step 2: Administrative Boundaries & sf Land-Clipped Buffer ----
message("Fetching Australian state boundaries...")
aus_states <- ne_download(scale = 50, type = "states", category = "cultural", returnclass = "sf") |>
  filter(admin == "Australia") |>
  st_make_valid()

message("Computing land-clipped range polygon via sf...")
aus_land_union <- st_union(aus_states)

# Buffer points by ~0.35 degrees (~40 km) and intersect with Australian landmass
species_range_sf <- all_occ_clean |>
  st_union() |>
  st_buffer(dist = 0.35) |>
  st_intersection(aus_land_union) |>
  st_sf()

# ---- Step 3: Set Map Limits ----
bbox_occ   <- st_bbox(all_occ_clean)
buffer_deg <- 1.5

map_limits <- c(
  xmin = unname(bbox_occ["xmin"]) - buffer_deg,
  xmax = unname(bbox_occ["xmax"]) + buffer_deg,
  ymin = unname(bbox_occ["ymin"]) - buffer_deg,
  ymax = unname(bbox_occ["ymax"]) + buffer_deg
)

# ---- Step 4: Render Minimalist Black & White Map ----
message("Rendering map...")

p_minimal <- ggplot() +
  # Australian Landmass / State Borders (White fill with solid black outlines)
  geom_sf(data = aus_states, fill = "white", color = "black", linewidth = 0.4) +
  # Land-clipped Buffered Species Range (Light grey fill, dashed black border)
  geom_sf(
    data = species_range_sf,
    fill = "grey90",
    color = "black",
    linetype = "dashed",
    linewidth = 0.5,
    alpha = 0.7
  ) +
  # Occurrence Points (Solid black points)
  geom_sf(
    data = all_occ_clean,
    color = "black",
    size = 0.5,
    alpha = 0.75
  ) +
  # Zoom Bounding Box
  coord_sf(
    xlim = c(map_limits["xmin"], map_limits["xmax"]),
    ylim = c(map_limits["ymin"], map_limits["ymax"]),
    expand = FALSE
  ) +
  # Scale Bar (Times New Roman text)
  annotation_scale(
    location = "bl",
    width_hint = 0.25,
    style = "ticks",
    text_size = 9,
    text_family = font_family
  ) +
  # North Arrow (Times New Roman text)
  annotation_north_arrow(
    location = "tl",
    which_north = "true",
    pad_x = unit(0.25, "in"),
    pad_y = unit(0.25, "in"),
    height = unit(1.2, "cm"),  # Scales up the arrow height (default is 1cm)
    width = unit(1.2, "cm"),   # Scales up the arrow width (default is 1cm)
    style = north_arrow_minimal(
      fill = "black",
      line_col = "black",
      text_family = font_family,
      text_size = 20            # Increases 'N' font size (default is 10)
    )
  ) +
  theme_bw(base_size = 12, base_family = font_family) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA),
    panel.grid       = element_blank(),
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.title       = element_blank(),
    axis.text        = element_text(family = font_family, color = "black", size = 10),
    axis.ticks       = element_line(color = "black"),
    plot.margin      = margin(6, 6, 6, 6)
  )

# ---- Step 5: Save Plot ----
out_filepath <- file.path(out_dir, "fig_species_range_sf_bw.png")
ggsave(out_filepath, p_minimal, width = 7, height = 8, dpi = 400)

message("Finished generating minimalist range map: ", out_filepath)
