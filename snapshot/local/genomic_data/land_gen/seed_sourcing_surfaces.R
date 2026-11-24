library(dplyr)
library(gdm)
library(terra)
library(ggplot2)
library(sf)
library(ggspatial)
library(patchwork)

# 1. Load data FIRST ----
load("gdm_outputs.RData")
load("lfmm_outputs.RData")

# 2. Rebuild grid coordinates & completeness mask ----
grid_current <- read.csv("current_covariate_grid.csv")
names(grid_current)[grep("^bio", names(grid_current))] <- gsub("_.*", "", grep("^bio", names(grid_current), value = TRUE))

target_vars <- colnames(grid_curr_clean)
complete_idx <- complete.cases(grid_current[, target_vars])
stopifnot(sum(complete_idx) == nrow(grid_curr_clean))

base_coords <- grid_current[complete_idx, c("lon", "lat")]

# 3. Build a fresh terra template AFTER loading data ----
ux <- sort(unique(grid_current$lon))
uy <- sort(unique(grid_current$lat))
res_x <- min(diff(ux)[diff(ux) > 1e-6])
res_y <- min(diff(uy)[diff(uy) > 1e-6])

r_template <- rast(
  ext(min(ux) - res_x/2, max(ux) + res_x/2, min(uy) - res_y/2, max(uy) + res_y/2),
  res = c(res_x, res_y),
  crs = "EPSG:4326"
)

rast_single <- function(z_vals, coords) {
  pts <- vect(data.frame(x = coords$lon, y = coords$lat, z = z_vals), geom = c("x", "y"))
  rasterize(pts, r_template, field = "z")
}

# Setup shared styling variables
if (.Platform$OS.type == "windows") {
  windowsFonts(Times = windowsFont("Times New Roman"))
}
font_family <- "Times"
pal_colors <- c("darkgreen", "lightgreen", "yellow", "orange", "red", "darkred")
site_order <- c("North", "Central", "South")

clean_predictors <- colnames(grid_curr_clean)
pred_means <- colMeans(grid_curr_clean[, clean_predictors])
pred_sds   <- apply(grid_curr_clean[, clean_predictors], 2, sd)
grid_curr_std <- scale(grid_curr_clean[, clean_predictors], center = pred_means, scale = pred_sds)

lfmm_pred_means <- colMeans(env_train_mat)
lfmm_pred_sds   <- apply(env_train_mat, 2, sd)

gap_from_B_local <- function(new_env, pred_env, means, sds) {
  scale_local <- function(m) scale(as.matrix(m), center = means, scale = sds)
  delta <- scale_local(new_env) - scale_local(pred_env)
  as.vector(rowSums((delta %*% BtB) * delta)) / nrow(B_eff)
}

create_site_map <- function(r, row_range, current_site_row, col_idx, row_idx) {
  df <- as.data.frame(r, xy = TRUE, na.rm = TRUE)
  colnames(df)[3] <- "z_val"
  
  p <- ggplot() +
    geom_raster(data = df, aes(x = x, y = y, fill = z_val)) +
    scale_fill_gradientn(colors = pal_colors, limits = row_range, oob = scales::squish, na.value = "transparent") +
    geom_sf(data = aus_states, fill = NA, color = "black", linewidth = 0.3, inherit.aes = FALSE) +
    geom_point(data = current_site_row, aes(x = lon, y = lat), shape = 21, fill = "white", color = "black", size = 2, stroke = 1.2) +
    coord_sf(xlim = c(ext_main["xmin"], ext_main["xmax"]), ylim = c(ext_main["ymin"], ext_main["ymax"]), expand = FALSE) +
    theme_bw(base_size = 10, base_family = font_family) +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background  = element_rect(fill = "white", color = NA),
      panel.grid       = element_blank(),
      panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.8),
      axis.title       = element_blank(),
      legend.position  = "none",
      plot.margin      = margin(2, 2, 2, 2)
    )
  
  if (col_idx == 1) {
    p <- p + theme(axis.text.y = element_text(family = font_family, color = "black", size = 7), axis.ticks.y = element_line(color = "black"))
  } else {
    p <- p + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  }
  
  if (row_idx == 3) {
    p <- p + theme(axis.text.x = element_text(family = font_family, color = "black", size = 7), axis.ticks.x = element_line(color = "black"))
  } else {
    p <- p + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  }
  
  return(p)
}

text_plot <- function(label, size = 4, face = "bold", angle = 0) {
  ggplot() + 
    annotate("text", x = 0.5, y = 0.5, label = label, family = font_family, size = size, fontface = face, angle = angle) +
    theme_void() +
    theme(plot.margin = margin(0, 0, 0, 0))
}

# Define scenarios including current, and set up aliases for uniform looping
scenarios <- c("current", "126", "370", "585")
trans_env_current  <- grid_trans_env
grid_current_clean <- grid_curr_clean

master_surfaces <- list()

# 4. Generate rasters for ALL scenarios (including current) first ----
for (scen in scenarios) {
  cat("Computing rasters for scenario:", scen, "...\n")
  trans_env_scen  <- get(paste0("trans_env_", scen))
  grid_scen_clean <- get(paste0("grid_", scen, "_clean"))
  
  # GDM site-match surfaces
  gdm_match_rasters <- list()
  for (i in seq_len(nrow(rest_sites))) {
    s_idx  <- rest_sites$idx[i]
    s_name <- rest_sites$site[i]
    target_env <- trans_env_scen[s_idx, , drop = FALSE]
    diff_mat   <- sweep(grid_trans_env, 2, as.numeric(target_env), "-")
    dist_vec   <- sqrt(rowSums(diff_mat^2))
    gdm_match_rasters[[s_name]] <- rast_single(dist_vec, base_coords)
  }
  
  # Climate-matching surfaces
  climate_match_rasters <- list()
  for (i in seq_len(nrow(rest_sites))) {
    s_idx  <- rest_sites$idx[i]
    s_name <- rest_sites$site[i]
    target_raw <- grid_scen_clean[s_idx, clean_predictors, drop = FALSE]
    target_std <- scale(target_raw, center = pred_means, scale = pred_sds)
    dist_vec   <- sqrt(rowSums(sweep(grid_curr_std, 2, as.numeric(target_std), "-")^2))
    climate_match_rasters[[s_name]] <- rast_single(dist_vec, base_coords)
  }
  
  # LFMM site-match surfaces
  lfmm_match_rasters <- list()
  for (i in seq_len(nrow(rest_sites))) {
    s_idx  <- rest_sites$idx[i]
    s_name <- rest_sites$site[i]
    target_env  <- grid_scen_clean[s_idx, colnames(env_train_mat), drop = FALSE]
    pred_rep    <- matrix(as.numeric(target_env), nrow = nrow(grid_curr_clean), ncol = ncol(target_env), byrow = TRUE)
    dist_vec    <- gap_from_B_local(grid_curr_clean[, colnames(env_train_mat), drop = FALSE], pred_rep, lfmm_pred_means, lfmm_pred_sds)
    lfmm_match_rasters[[s_name]] <- rast_single(dist_vec, base_coords)
  }
  
  master_surfaces[[scen]] <- list(
    "Climate Distance" = climate_match_rasters,
    "GDM Offset"       = gdm_match_rasters,
    "LFMM Offset"      = lfmm_match_rasters
  )
}

# 5. Compute global value ranges across ALL scenarios (current + future) for each surface type ----
surface_names <- c("Climate Distance", "GDM Offset", "LFMM Offset")
global_row_ranges <- list()

for (s_name in surface_names) {
  all_vals <- c()
  for (scen in scenarios) {
    site_list <- master_surfaces[[scen]][[s_name]]
    for (site_name in site_order) {
      all_vals <- c(all_vals, values(site_list[[site_name]]))
    }
  }
  global_row_ranges[[s_name]] <- range(all_vals, na.rm = TRUE)
}

# 6. Build and save plots using the global shared ranges ----
for (scen in scenarios) {
  cat("Generating and saving plot for scenario:", scen, "...\n")
  all_surfaces <- master_surfaces[[scen]]
  
  row_list <- list()
  row_idx <- 1
  for (surface_name in surface_names) {
    row_range <- global_row_ranges[[surface_name]]
    col_plots <- list()
    col_idx <- 1
    
    for (site_name in site_order) {
      r <- all_surfaces[[surface_name]][[site_name]]
      s_row <- rest_sites[rest_sites$site == site_name, ]
      col_plots[[col_idx]] <- create_site_map(r, row_range, s_row, col_idx, row_idx)
      col_idx <- col_idx + 1
    }
    
    right_label <- text_plot(surface_name, size = 4, face = "bold", angle = -90)
    row_list[[row_idx]] <- wrap_plots(col_plots[[1]], col_plots[[2]], col_plots[[3]], right_label, ncol = 4) + 
      plot_layout(widths = c(1, 1, 1, 0.15))
    row_idx <- row_idx + 1
  }
  
  header_row <- wrap_plots(text_plot("North", 4, "bold"), text_plot("Central", 4, "bold"), text_plot("South", 4, "bold"), plot_spacer(), ncol = 4) + 
    plot_layout(widths = c(1, 1, 1, 0.15))
  
  final_panel <- header_row / row_list[[1]] / row_list[[2]] / row_list[[3]] +
    plot_layout(heights = c(0.12, 1, 1, 1))
  
  out_filename <- if (scen == "current") {
    "fig_seed_sourcing_surfaces_current_ggplot.png"
  } else {
    paste0("fig_seed_sourcing_surfaces_ssp", scen, "_ggplot.png")
  }
  
  ggsave(filename = out_filename, plot = final_panel, width = 10, height = 11, dpi = 400, device = "png")
  cat("Successfully saved:", out_filename, "\n")
}