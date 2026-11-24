library(tidyverse)
library(sf)
library(rnaturalearth)
library(ggspatial)
library(FNN)
library(ggtext) # Required for markdown legend parsing

# 1. Clear open graphics devices
while (!is.null(dev.list())) { dev.off() }

# 2. Register Times New Roman Font
windowsFonts(Times = windowsFont("Times New Roman"))
font_family <- "Times"

# ---- Paths Setup ----
base_dir <- "C:/Users/evane/OneDrive/Desktop/Masters Dissertation/genomic_data/land_gen"
grid_fp  <- file.path(base_dir, "current_covariate_grid.csv")
meta_fp  <- "C:/Users/evane/OneDrive/Desktop/Masters Dissertation/spatial_data/sample_metadata.csv"

out_dir  <- file.path(base_dir, "figures")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

ala_dir      <- "C:/Users/evane/OneDrive/Desktop/Masters Dissertation/spatial_data/ALA_occurences"
pauci_ala_fp <- file.path(ala_dir, "pauci_range", "pauci_range.csv")
niph_ala_fp  <- file.path(ala_dir, "niph_range", "niph_range.csv")

env_vars <- c("bio15_1981-2010", "bio06_1981-2010", "bio04_1981-2010", "bio17_1981-2010")

grid_lon_col <- "lon"
grid_lat_col <- "lat"

# ---- Colors & Formatting (Color 1: Orange, Color 2: Sky Blue) ----
taxa_labels <- c(
  "subsp. pauciflora" = "subsp. *pauciflora*",
  "subsp. niphophila" = "subsp. *niphophila*"
)

sample_colors <- c(
  "subsp. pauciflora" = "#E69F00",  # Okabe-Ito Color 1 (Orange)
  "subsp. niphophila" = "#56B4E9"   # Okabe-Ito Color 2 (Sky Blue)
)

ala_colors <- c(
  "subsp. pauciflora" = "#E69F00",
  "subsp. niphophila" = "#56B4E9"
)

grid <- read_csv(grid_fp, show_col_types = FALSE)

# ---- Load unified sample metadata ----
samples_raw <- read_csv(meta_fp, show_col_types = FALSE)

stopifnot(
  "sample_metadata.csv is missing expected columns" =
    all(c("sample", "species", "lat", "long", "seq_tech") %in% names(samples_raw))
)

samples_meta <- samples_raw |>
  filter(
    str_to_lower(seq_tech) == "illumina",
    species %in% c("E. pauciflora subsp. pauciflora", "E. pauciflora subsp. niphophila")
  ) |>
  mutate(
    subspecies = case_when(
      str_detect(species, regex("niphophila", ignore_case = TRUE)) ~ "subsp. niphophila",
      str_detect(species, regex("pauciflora", ignore_case = TRUE)) ~ "subsp. pauciflora",
      TRUE ~ "other/unresolved"
    ),
    lat_num  = suppressWarnings(as.numeric(lat)),
    long_num = suppressWarnings(as.numeric(long))
  )

samples_meta <- samples_meta |>
  filter(is.finite(lat_num), is.finite(long_num),
         long_num >= -180, long_num <= 180, lat_num >= -90, lat_num <= 90)

filtered_samples_fp <- file.path(base_dir, "filtered_sample_ids.txt")
writeLines(samples_meta$sample, filtered_samples_fp)

# ---- Extract climate for each sample via nearest-neighbor join ----
grid_coords   <- grid |> select(all_of(c(grid_lon_col, grid_lat_col))) |> as.matrix()
sample_coords <- samples_meta |> select(long_num, lat_num) |> as.matrix()

nn_samples <- FNN::get.knnx(data = grid_coords, query = sample_coords, k = 1)
samples_meta$grid_nn_dist <- nn_samples$nn.dist[, 1]

sample_env <- grid[nn_samples$nn.index[, 1], env_vars] |>
  bind_cols(samples_meta |> select(sample, subspecies, lat_num, long_num))

# ---- PCA fit on full climate grid ----
grid_env <- grid |> select(all_of(env_vars)) |> drop_na()
pca_fit  <- prcomp(grid_env, center = TRUE, scale. = TRUE)

pct_var <- round(100 * summary(pca_fit)$importance[2, 1:2], 1)

loadings <- as_tibble(pca_fit$rotation[, 1:2], rownames = "variable") |>
  mutate(variable = str_remove(variable, "_1981-2010"))

sample_scores <- predict(pca_fit, newdata = sample_env |> select(all_of(env_vars))) |>
  as_tibble() |>
  select(PC1, PC2) |>
  mutate(source = "Sampled", subspecies = sample_env$subspecies)

arrow_scale <- max(abs(sample_scores$PC1), abs(sample_scores$PC2)) * 0.8

# ---- ALA occurrences ----
pauci_ala <- read_csv(pauci_ala_fp, show_col_types = FALSE) |>
  mutate(subspecies_ala = "subsp. pauciflora")

niph_ala <- read_csv(niph_ala_fp, show_col_types = FALSE) |>
  mutate(subspecies_ala = "subsp. niphophila")

ala_occ <- bind_rows(pauci_ala, niph_ala) |>
  filter(!is.na(decimalLatitude), !is.na(decimalLongitude)) |>
  distinct(decimalLatitude, decimalLongitude, subspecies_ala, .keep_all = TRUE)

ala_coords <- ala_occ |> select(decimalLongitude, decimalLatitude) |> as.matrix()
nn_ala <- FNN::get.knnx(data = grid_coords, query = ala_coords, k = 1)

ala_occ$nn_dist <- nn_ala$nn.dist[, 1]

ala_env <- grid[nn_ala$nn.index[, 1], env_vars] |>
  bind_cols(subspecies_ala = ala_occ$subspecies_ala, nn_dist = ala_occ$nn_dist,
            lat = ala_occ$decimalLatitude, lon = ala_occ$decimalLongitude)

ala_scores <- predict(pca_fit, newdata = ala_env |> select(all_of(env_vars))) |>
  as_tibble() |>
  select(PC1, PC2) |>
  mutate(subspecies = ala_env$subspecies_ala, source = "ALA occurrence",
         lat = ala_env$lat, lon = ala_env$lon)

# =====================================================================
# FIGURE 1: Environmental PCA (No Title)
# =====================================================================
p_pca_ala <- ggplot() +
  geom_point(
    data = ala_scores, aes(PC1, PC2, fill = subspecies),
    shape = 21, color = "white", stroke = 0.1, size = 1.8, alpha = 0.5
  ) +
  geom_point(
    data = sample_scores, aes(PC1, PC2, fill = subspecies),
    shape = 21, color = "black", stroke = 0.35, size = 2.2, alpha = 0.95
  ) +
  geom_segment(
    data = loadings,
    aes(x = 0, y = 0, xend = PC1 * arrow_scale, yend = PC2 * arrow_scale),
    arrow = arrow(length = unit(0.2, "cm")), color = "grey20", linewidth = 0.5
  ) +
  geom_text(
    data = loadings,
    aes(x = PC1 * arrow_scale * 1.15, y = PC2 * arrow_scale * 1.15, label = variable),
    size = 3.2, color = "grey20", fontface = "italic", family = font_family
  ) +
  scale_fill_manual(
    values = ala_colors,
    labels = taxa_labels,
    name = "Subspecies", 
    guide = guide_legend(override.aes = list(shape = 21, size = 3, color = "black", stroke = 0.35, alpha = 0.95))
  ) +
  labs(
    x = paste0("PC1 (", pct_var[1], "%)"),
    y = paste0("PC2 (", pct_var[2], "%)")
  ) +
  theme_minimal(base_size = 12, base_family = font_family) +
  theme(
    legend.position = "bottom",
    legend.title    = element_text(size = 10, face = "bold", family = font_family),
    legend.text     = element_markdown(size = 9, family = font_family),
    panel.grid.minor = element_blank(),
    text = element_text(family = font_family)
  )

png(filename = file.path(out_dir, "fig_environmental_pca_with_ALA.png"), width = 7 * 400, height = 7 * 400, res = 400)
print(p_pca_ala)
dev.off()

# =====================================================================
# Shared Spatial Prep (Australia B&W Outlines)
# =====================================================================
pts_sf <- samples_meta |>
  st_as_sf(coords = c("long_num", "lat_num"), crs = 4326)

ala_pts_sf <- ala_occ |>
  st_as_sf(coords = c("decimalLongitude", "decimalLatitude"), crs = 4326) |>
  rename(subspecies = subspecies_ala)

aus_states <- ne_download(scale = 50, type = "states", category = "cultural", returnclass = "sf") |>
  filter(admin == "Australia") |>
  st_make_valid()

aus_land_union <- st_union(aus_states)

bbox_main <- st_bbox(c(
  xmin = min(st_bbox(pts_sf)["xmin"], st_bbox(ala_pts_sf)["xmin"]),
  xmax = max(st_bbox(pts_sf)["xmax"], st_bbox(ala_pts_sf)["xmax"]),
  ymin = min(st_bbox(pts_sf)["ymin"], st_bbox(ala_pts_sf)["ymin"]),
  ymax = max(st_bbox(pts_sf)["ymax"], st_bbox(ala_pts_sf)["ymax"])
), crs = 4326)

buffer_deg_main <- 1.0
ext_main <- c(
  xmin = unname(bbox_main["xmin"]) - buffer_deg_main,
  xmax = unname(bbox_main["xmax"]) + buffer_deg_main,
  ymin = unname(bbox_main["ymin"]) - buffer_deg_main,
  ymax = unname(bbox_main["ymax"]) + buffer_deg_main
)

# =====================================================================
# FIGURE 2: Sampled Sites Overlaid on ALA Occurrences
# =====================================================================
p_map <- ggplot() +
  geom_sf(data = aus_states, fill = "white", color = "black", linewidth = 0.4) +
  geom_sf(
    data = ala_pts_sf, aes(color = subspecies),
    size = 0.8, alpha = 0.4, inherit.aes = FALSE, show.legend = FALSE
  ) +
  scale_color_manual(values = ala_colors, guide = "none") +
  ggnewscale::new_scale_fill() + 
  geom_sf(
    data = pts_sf, aes(fill = subspecies),
    shape = 21, color = "black", stroke = 0.35, size = 2.2, alpha = 0.95,
    inherit.aes = FALSE
  ) +
  scale_fill_manual(
    values = sample_colors,
    labels = taxa_labels,
    name = "Subspecies",
    guide = guide_legend(override.aes = list(size = 3.5, shape = 21, color = "black", stroke = 0.35))
  ) +
  coord_sf(
    xlim = c(ext_main["xmin"], ext_main["xmax"]),
    ylim = c(ext_main["ymin"], ext_main["ymax"]),
    expand = FALSE
  ) +
  annotation_scale(
    location = "bl",
    width_hint = 0.2,
    style = "ticks",
    text_family = font_family
  ) +
  annotation_north_arrow(
    location = "tl",
    which_north = "true",
    pad_x = unit(0.2, "in"),
    pad_y = unit(0.2, "in"),
    height = unit(1.0, "cm"),
    width = unit(1.0, "cm"),
    style = north_arrow_minimal(
      fill = "black",
      line_col = "black",
      text_family = font_family,
      text_size = 14
    )
  ) +
  theme_bw(base_size = 11, base_family = font_family) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA),
    panel.grid       = element_blank(),
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.title       = element_blank(),
    axis.text        = element_text(family = font_family, color = "black", size = 9),
    axis.ticks       = element_line(color = "black"),
    legend.position      = c(0.02, 0.15),
    legend.justification = c(0, 0),
    legend.background    = element_blank(),
    legend.key           = element_blank(),
    legend.text          = element_markdown(size = 9, family = font_family),
    legend.title         = element_text(size = 10, face = "bold", family = font_family),
    plot.margin          = margin(6, 6, 6, 6)
  )

png(filename = file.path(out_dir, "fig_sample_map.png"), width = 7.5 * 400, height = 7 * 400, res = 400)
print(p_map)
dev.off()

# =====================================================================
# FIGURE 3: Standalone pauciflora Map (Color 1: Orange)
# =====================================================================
pauci_config <- list(
  code = "pauci",
  csv_subpath = file.path("pauci_range", "pauci_range.csv"),
  label = "E. pauciflora subsp. pauciflora",
  color = "#E69F00"
)

pauci_occ_df <- pauci_ala |>
  mutate(
    decimalLatitude  = as.numeric(decimalLatitude),
    decimalLongitude = as.numeric(decimalLongitude)
  ) |>
  filter(
    is.finite(decimalLatitude), is.finite(decimalLongitude),
    decimalLongitude >= 129.0, decimalLongitude <= 155.0,
    decimalLatitude  >= -45.0, decimalLatitude  <= -10.0
  ) |>
  distinct(decimalLatitude, decimalLongitude, .keep_all = TRUE)

pauci_occ_sf <- st_as_sf(pauci_occ_df, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)

pauci_range_sf <- suppressWarnings(
  pauci_occ_sf |>
    summarise(geometry = st_union(geometry)) |>
    st_buffer(dist = 0.35) |>
    st_intersection(aus_land_union)
)

bbox_pauci <- st_bbox(pauci_occ_sf)
limits_pauci <- c(
  xmin = unname(bbox_pauci["xmin"]) - 1.0,
  xmax = unname(bbox_pauci["xmax"]) + 1.0,
  ymin = unname(bbox_pauci["ymin"]) - 1.0,
  ymax = unname(bbox_pauci["ymax"]) + 1.0
)

p_pauci <- ggplot() +
  geom_sf(data = aus_states, fill = "white", color = "black", linewidth = 0.4) +
  geom_sf(
    data = pauci_range_sf,
    fill = pauci_config$color,
    color = "black",
    linetype = "dashed",
    linewidth = 0.4,
    alpha = 0.3
  ) +
  geom_sf(
    data = pauci_occ_sf,
    fill = pauci_config$color,
    color = "black",
    shape = 21,
    size = 2.0,
    stroke = 0.3,
    alpha = 0.9
  ) +
  coord_sf(
    xlim = c(limits_pauci["xmin"], limits_pauci["xmax"]),
    ylim = c(limits_pauci["ymin"], limits_pauci["ymax"]),
    expand = FALSE
  ) +
  annotation_scale(
    location = "bl",
    width_hint = 0.2,
    style = "ticks",
    text_family = font_family
  ) +
  annotation_north_arrow(
    location = "tl",
    which_north = "true",
    pad_x = unit(0.2, "in"),
    pad_y = unit(0.2, "in"),
    height = unit(1.0, "cm"),
    width = unit(1.0, "cm"),
    style = north_arrow_minimal(
      fill = "black",
      line_col = "black",
      text_family = font_family,
      text_size = 14
    )
  ) +
  theme_bw(base_size = 11, base_family = font_family) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA),
    panel.grid       = element_blank(),
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.title       = element_blank(),
    axis.text        = element_text(family = font_family, color = "black", size = 9),
    axis.ticks       = element_line(color = "black"),
    legend.position  = "none",
    plot.margin      = margin(6, 6, 6, 6)
  )

png(filename = file.path(out_dir, "fig_pauciflora_clean.png"), width = 8.5 * 400, height = 8.5 * 400, res = 400)
print(p_pauci)
dev.off()

# =====================================================================
# FIGURE 4: Standalone niphophila Map (Color 2: Sky Blue)
# =====================================================================
niph_config <- list(
  code = "niph",
  csv_subpath = file.path("niph_range", "niph_range.csv"),
  label = "E. pauciflora subsp. niphophila",
  color = "#56B4E9"
)

niph_occ_df <- niph_ala |>
  mutate(
    decimalLatitude  = as.numeric(decimalLatitude),
    decimalLongitude = as.numeric(decimalLongitude)
  ) |>
  filter(
    is.finite(decimalLatitude), is.finite(decimalLongitude),
    decimalLongitude >= 129.0, decimalLongitude <= 155.0,
    decimalLatitude  >= -45.0, decimalLatitude  <= -10.0
  ) |>
  distinct(decimalLatitude, decimalLongitude, .keep_all = TRUE)

niph_occ_sf <- st_as_sf(niph_occ_df, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)

niph_range_sf <- suppressWarnings(
  niph_occ_sf |>
    summarise(geometry = st_union(geometry)) |>
    st_buffer(dist = 0.35) |>
    st_intersection(aus_land_union)
)

bbox_niph <- st_bbox(niph_occ_sf)
limits_niph <- c(
  xmin = unname(bbox_niph["xmin"]) - 1.0,
  xmax = unname(bbox_niph["xmax"]) + 1.0,
  ymin = unname(bbox_niph["ymin"]) - 1.0,
  ymax = unname(bbox_niph["ymax"]) + 1.0
)

p_niph <- ggplot() +
  geom_sf(data = aus_states, fill = "white", color = "black", linewidth = 0.4) +
  geom_sf(
    data = niph_range_sf,
    fill = niph_config$color,
    color = "black",
    linetype = "dashed",
    linewidth = 0.4,
    alpha = 0.3
  ) +
  geom_sf(
    data = niph_occ_sf,
    fill = niph_config$color,
    color = "black",
    shape = 21,
    size = 2.0,
    stroke = 0.3,
    alpha = 0.9
  ) +
  coord_sf(
    xlim = c(limits_niph["xmin"], limits_niph["xmax"]),
    ylim = c(limits_niph["ymin"], limits_niph["ymax"]),
    expand = FALSE
  ) +
  annotation_scale(
    location = "bl",
    width_hint = 0.2,
    style = "ticks",
    text_family = font_family
  ) +
  annotation_north_arrow(
    location = "tl",
    which_north = "true",
    pad_x = unit(0.2, "in"),
    pad_y = unit(0.2, "in"),
    height = unit(1.0, "cm"),
    width = unit(1.0, "cm"),
    style = north_arrow_minimal(
      fill = "black",
      line_col = "black",
      text_family = font_family,
      text_size = 14
    )
  ) +
  theme_bw(base_size = 11, base_family = font_family) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA),
    panel.grid       = element_blank(),
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.title       = element_blank(),
    axis.text        = element_text(family = font_family, color = "black", size = 9),
    axis.ticks       = element_line(color = "black"),
    legend.position  = "none",
    plot.margin      = margin(6, 6, 6, 6)
  )

png(filename = file.path(out_dir, "fig_niphophila_clean.png"), width = 8.5 * 400, height = 8.5 * 400, res = 400)
print(p_niph)
dev.off()

message("All figures (PCA, Sample Map, Pauci Map, Niph Map) successfully generated and saved to: ", out_dir)