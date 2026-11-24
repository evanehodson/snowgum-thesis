library(gdm)
library(readr)
library(dplyr)
library(terra)
library(geodist)
library(stringr)
library(ggplot2)
library(grid)
library(gridExtra)

options(scipen = 999)
RUN_VARIMP    <- TRUE

# 0. Multi-model future grid ensemble processor -----------------------------
model_names <- c("GFDL-ESM4", "IPSL-CM6A-LR", "MPI-ESM1-2-HR", "MRI-ESM2-0", "UKESM1-0-LL")
model_alt   <- paste(model_names, collapse = "|")
id_cols     <- c("grid_id", "lon", "lat", "row", "col")

process_grid_ssp <- function(raw_grid, current_grid, shared_predictors, ssp_target) {
  ssp_cols <- names(raw_grid) %>% str_subset(paste0("_", ssp_target, "$"))
  
  var_names <- ssp_cols %>%
    str_remove(paste0("_(", model_alt, ")_", ssp_target, "$")) %>%
    str_remove("_\\d{4}-\\d{4}$") %>%
    unique()
  
  ens_mean   <- raw_grid %>% select(any_of(id_cols))
  ens_sd     <- raw_grid %>% select(any_of(id_cols))
  per_model  <- setNames(lapply(model_names, function(m) raw_grid %>% select(any_of(id_cols))), model_names)
  stats_list <- vector("list", length(var_names))
  matched_predictors <- character(0)
  
  for (i in seq_along(var_names)) {
    v <- var_names[i]
    v_escaped <- str_replace_all(v, "([.\\-])", "\\\\\\1")
    cols_v <- ssp_cols[str_detect(ssp_cols, paste0("^", v_escaped, "_"))]
    mat_v  <- as.matrix(raw_grid[, cols_v])
    
    target_name <- shared_predictors[str_detect(shared_predictors, paste0("^", v_escaped, "(_|$)"))]
    if (length(target_name) == 0) {
      warning(paste0("process_grid_ssp: no shared_predictors match for '", v, "' (", ssp_target, "); dropping this variable"))
      next
    }
    if (length(target_name) > 1) {
      warning(paste0("process_grid_ssp: multiple shared_predictors matches for '", v, "': ",
                     paste(target_name, collapse = ", "), " - using the first"))
      target_name <- target_name[1]
    }
    
    ens_mean[[target_name]] <- rowMeans(mat_v, na.rm = TRUE)
    ens_sd[[target_name]]   <- apply(mat_v, 1, sd, na.rm = TRUE)
    matched_predictors <- c(matched_predictors, target_name)
    
    for (m in model_names) {
      col_m <- cols_v[str_detect(cols_v, fixed(m))]
      if (length(col_m) == 1) {
        per_model[[m]][[target_name]] <- raw_grid[[col_m]]
      } else if (length(col_m) > 1) {
        warning(paste0("process_grid_ssp: multiple '", m, "' columns matched for '", v, "'; using the first"))
        per_model[[m]][[target_name]] <- raw_grid[[col_m[1]]]
      } else {
        per_model[[m]][[target_name]] <- NA_real_
      }
    }
    
    stats_list[[i]] <- data.frame(
      variable  = v,
      predictor = target_name,
      ssp       = ssp_target,
      n_models  = ncol(mat_v),
      grid_mean = mean(mat_v, na.rm = TRUE),
      grid_sd   = mean(ens_sd[[target_name]], na.rm = TRUE),
      cv_pct    = 100 * mean(ens_sd[[target_name]], na.rm = TRUE) / mean(mat_v, na.rm = TRUE)
    )
  }
  
  static_present <- setdiff(shared_predictors, matched_predictors)
  static_present <- intersect(static_present, colnames(current_grid))
  if (length(static_present) > 0) {
    ens_mean <- ens_mean %>%
      left_join(current_grid %>% select(grid_id, all_of(static_present)), by = "grid_id")
    per_model <- lapply(per_model, function(df) {
      df %>% left_join(current_grid %>% select(grid_id, all_of(static_present)), by = "grid_id")
    })
  }
  
  missing_predictors <- setdiff(shared_predictors, colnames(ens_mean))
  if (length(missing_predictors) > 0) {
    warning(paste0("process_grid_ssp (", ssp_target, "): still missing from output grid: ",
                   paste(missing_predictors, collapse = ", ")))
  }
  
  list(mean = ens_mean, sd = ens_sd, per_model = per_model, stats = bind_rows(stats_list))
}

# 1. Load genetic and spatial data ----------------------------------------
dist_raw <- read.table("genetic_de.mdist", header = FALSE)
ids      <- read.table("genetic_de.mdist.id", header = FALSE)
dist_mat <- as.matrix(dist_raw)
rownames(dist_mat) <- ids$V2
colnames(dist_mat) <- ids$V2

env_present <- read_csv("current_covariates.csv", show_col_types = FALSE)
grid_data   <- read_csv("current_covariate_grid.csv", show_col_types = FALSE)

pauci_samples <- readLines("pauci_samples.txt")
valid_samples <- intersect(pauci_samples, rownames(dist_mat))

site_data <- env_present %>%
  filter(sample_id %in% valid_samples) %>%
  mutate(site_num = as.integer(factor(sample_id))) %>%
  select(site_num, sample_id, long, lat, elev, everything()) %>%
  as.data.frame()

dist_mat_sub <- dist_mat[site_data$sample_id, site_data$sample_id]
dist_df      <- cbind(site_num = site_data$site_num, as.data.frame(dist_mat_sub))
colnames(dist_df) <- c("site_num", site_data$sample_id)

shared_predictors <- intersect(
  setdiff(colnames(env_present), c("sample_id", "long", "lat", "elev")),
  setdiff(colnames(grid_data), c("grid_id", "lon", "lat", "row", "col"))
)

pred_input <- site_data %>% select(site_num, long, lat, all_of(shared_predictors))

# 2. Format site-pairs & filter spatial proximity (> 1km) -----------------
gdm_input_raw <- formatsitepair(
  bioData    = dist_df,
  bioFormat  = 3,
  predData   = pred_input,
  siteColumn = "site_num",
  XColumn    = "long",
  YColumn    = "lat",
  verbose    = TRUE
)

pair_dists <- geodist(
  x = data.frame(lon = gdm_input_raw$s1.xCoord, lat = gdm_input_raw$s1.yCoord),
  y = data.frame(lon = gdm_input_raw$s2.xCoord, lat = gdm_input_raw$s2.yCoord),
  paired = TRUE,
  measure = "haversine"
)

gdm_input_filtered <- gdm_input_raw[pair_dists > 1000, ]

# 3. Fit GDM models & Identify Active Predictors ---------------------------
gdm_mod       <- gdm(gdm_input_filtered, geo = TRUE)

spline_check <- isplineExtract(gdm_mod)
max_y_nogeo  <- apply(spline_check$y, 2, max)
eps <- 1e-6
active_env_vars <- setdiff(names(max_y_nogeo)[max_y_nogeo > eps], "Geographic")

varImp_result <- if (RUN_VARIMP) {
  gdm.varImp(gdm_input_filtered, geo = TRUE)
} else {
  cat("Skipping gdm.varImp (RUN_VARIMP = FALSE)\n")
  NULL
}

if(RUN_VARIMP) print(varImp_result[[2]])

# 4. Transform environmental grids -----------------------------------------
grid_preds_bare  <- grid_data %>% select(long = lon, lat, all_of(shared_predictors)) %>% as.data.frame()
grid_transformed <- gdm.transform(gdm_mod, grid_preds_bare)
grid_trans_mat   <- as.matrix(grid_transformed[, active_env_vars, drop = FALSE])

future_grid_raw <- read_csv("future_covariate_grid.csv", show_col_types = FALSE)

ssp126_out <- process_grid_ssp(future_grid_raw, grid_data, shared_predictors, "ssp126")
ssp370_out <- process_grid_ssp(future_grid_raw, grid_data, shared_predictors, "ssp370")
ssp585_out <- process_grid_ssp(future_grid_raw, grid_data, shared_predictors, "ssp585")

grid_126 <- ssp126_out$mean
grid_370 <- ssp370_out$mean
grid_585 <- ssp585_out$mean

model_agreement <- bind_rows(ssp126_out$stats, ssp370_out$stats, ssp585_out$stats)

trans_mat_126 <- as.matrix(gdm.transform(gdm_mod, grid_126 %>% select(long = lon, lat, all_of(shared_predictors)) %>% as.data.frame()))[, active_env_vars, drop = FALSE]
trans_mat_370 <- as.matrix(gdm.transform(gdm_mod, grid_370 %>% select(long = lon, lat, all_of(shared_predictors)) %>% as.data.frame()))[, active_env_vars, drop = FALSE]
trans_mat_585 <- as.matrix(gdm.transform(gdm_mod, grid_585 %>% select(long = lon, lat, all_of(shared_predictors)) %>% as.data.frame()))[, active_env_vars, drop = FALSE]

complete_rows <- complete.cases(grid_trans_mat) & 
  complete.cases(trans_mat_126) & 
  complete.cases(trans_mat_370) & 
  complete.cases(trans_mat_585)

env_cols  <- active_env_vars

grid_trans_env <- grid_trans_mat[complete_rows, env_cols, drop = FALSE]
trans_env_126  <- trans_mat_126[complete_rows, env_cols, drop = FALSE]
trans_env_370  <- trans_mat_370[complete_rows, env_cols, drop = FALSE]
trans_env_585  <- trans_mat_585[complete_rows, env_cols, drop = FALSE]

# 5. Build robust terra template raster ------------------------------------
ux <- sort(unique(grid_data$lon))
uy <- sort(unique(grid_data$lat))
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

base_coords <- grid_data[complete_rows, c("lon", "lat")]

# 6. Single PCA basis (Grid-derived) & Linear-hue PC1 Mapping -----------
pca_fit <- prcomp(grid_trans_env, center = TRUE, scale. = FALSE)
pca_scores_curr <- pca_fit$x[, 1:3, drop = FALSE]

pc1_vals <- pca_scores_curr[, 1]
pc1_norm <- (pc1_vals - min(pc1_vals)) / (max(pc1_vals) - min(pc1_vals))
hue_vals <- pc1_norm * 270

pixel_colors_pc1 <- hsv(h = hue_vals / 360, s = 0.75, v = 0.85)

xyz_pc1 <- data.frame(x = base_coords$lon, y = base_coords$lat, col_id = seq_along(pixel_colors_pc1))
pts_pc1 <- vect(xyz_pc1, geom = c("x", "y"))
r_pc1 <- rasterize(pts_pc1, r_template, field = "col_id")

# 7. Native GDM Offsets (Isolating Climate Turnover, Suppressing Spatial IBD) ---
env_model_preds <- setdiff(gdm_mod$predictors, "Geographic")

curr_input <- grid_data[complete_rows, c("lon", "lat", env_model_preds)] %>% as.data.frame()
fut_126_in <- grid_126[complete_rows, c("lon", "lat", env_model_preds)] %>% as.data.frame()
fut_370_in <- grid_370[complete_rows, c("lon", "lat", env_model_preds)] %>% as.data.frame()
fut_585_in <- grid_585[complete_rows, c("lon", "lat", env_model_preds)] %>% as.data.frame()

trans_curr_full <- gdm.transform(gdm_mod, curr_input)
trans_126_full  <- gdm.transform(gdm_mod, fut_126_in)
trans_370_full  <- gdm.transform(gdm_mod, fut_370_in)
trans_585_full  <- gdm.transform(gdm_mod, fut_585_in)

trans_curr <- as.matrix(trans_curr_full[, env_model_preds, drop = FALSE])
trans_126  <- as.matrix(trans_126_full[, env_model_preds, drop = FALSE])
trans_370  <- as.matrix(trans_370_full[, env_model_preds, drop = FALSE])
trans_585  <- as.matrix(trans_585_full[, env_model_preds, drop = FALSE])

l1_126 <- rowSums(abs(trans_curr - trans_126))
l1_370 <- rowSums(abs(trans_curr - trans_370))
l1_585 <- rowSums(abs(trans_curr - trans_585))

offset_126 <- 1 - exp(-l1_126)
offset_370 <- 1 - exp(-l1_370)
offset_585 <- 1 - exp(-l1_585)

r_off_126 <- rast_single(offset_126, base_coords)
r_off_370 <- rast_single(offset_370, base_coords)
r_off_585 <- rast_single(offset_585, base_coords)

global_min_off <- min(c(offset_126, offset_370, offset_585), na.rm = TRUE)
global_max_off <- max(c(offset_126, offset_370, offset_585), na.rm = TRUE)

# 8. Extract splines & diagnostic matrices ---------------------------------
splineData <- isplineExtract(gdm_mod)
cor_mat <- cor(site_data[, active_env_vars, drop = FALSE], use = "complete.obs")
cat("Max abs correlation among active predictors:", max(abs(cor_mat[upper.tri(cor_mat)])), "\n")

pal_heat   <- colorRampPalette(c("lightgreen", "yellow", "orange", "red", "darkred"))(100)

# 9. Compile Figure Outputs ------------------------------------------------

# --- Figure 1: Observed vs predicted scatter ---
png("fig1_obs_vs_pred.png", width = 600, height = 600, res = 150)
par(mar = c(4, 4, 2, 1))

set.seed(42)
n_pts <- length(gdm_mod$observed)
sub_idx <- sample(seq_len(n_pts), size = min(n_pts, 300))

plot(gdm_mod$predicted[sub_idx], gdm_mod$observed[sub_idx],
     xlab = "Predicted genomic distance", ylab = "Observed genomic distance",
     col = rgb(0, 0, 0, 0.25), pch = 16, cex = 0.8)
abline(0, 1, col = "black", lwd = 1.2, lty = 2)
dev.off()


# --- Figure 2: 2x2 spline panel (Geographic + 3 active predictors) ---
png("fig2_splines.png", width = 1200, height = 1200, res = 150)
par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
max_y_val <- max(splineData$y)

geo_dist_degrees <- sqrt(
  (gdm_input_filtered$s1.xCoord - gdm_input_filtered$s2.xCoord)^2 +
    (gdm_input_filtered$s1.yCoord - gdm_input_filtered$s2.yCoord)^2
)

if ("Geographic" %in% colnames(splineData$y)) {
  geo_x <- splineData$x[, "Geographic"]
  geo_y <- splineData$y[, "Geographic"]
  geo_ord <- order(geo_x)
  
  plot(geo_x[geo_ord], geo_y[geo_ord], type = "l", col = "black", lwd = 1.8,
       xlab = "Geographic distance (decimal degrees)", ylab = "genomic distance",
       main = "Geographic distance", ylim = c(0, max_y_val))
  
  lines(geo_x[geo_ord], geo_y[geo_ord], col = "black", lwd = 1.8)
  points(geo_dist_degrees[sub_idx], gdm_mod$observed[sub_idx], col = rgb(0, 0, 0, 0.25), pch = 16, cex = 0.4)
}

max_vals <- apply(splineData$y[, setdiff(colnames(splineData$y), "Geographic"), drop = FALSE], 2, max)
max_vals <- max_vals[max_vals > eps] 
sorted_vars <- names(sort(max_vals, decreasing = TRUE))[1:3]

# --- UPDATED CONVERSION FUNCTION ---
chelsa_convert <- function(varname, x_vals) {
  var_root <- gsub("_.*", "", varname)
  abs_temp_vars <- c("bio01", "bio02", "bio05", "bio06", "bio07", "bio08", "bio09", "bio10", "bio11")
  precip_vars   <- c("bio12", "bio13", "bio14", "bio15", "bio16", "bio17", "bio18", "bio19")
  
  if (var_root == "bio04") {
    return(x_vals / 1000)            # CHELSA v2 bio04 scale (dK * 100 to °C SD)
  } else if (var_root %in% abs_temp_vars) {
    return((x_vals / 10) - 273.15)  # CHELSA v2 temperature scale (K/10 to °C)
  } else if (var_root %in% precip_vars) {
    return(x_vals / 10)            # CHELSA v2 precipitation scale (1/10 mm to mm)
  } else {
    return(x_vals)
  }
}

axis_labels <- c(
  "bio06_1981-2010" = "Min temp of coldest month (\u00b0C)",
  "bio04_1981-2010" = "Temperature seasonality (\u00b0C SD)",
  "bio17_1981-2010" = "Precipitation of driest quarter (mm)"
)

title_labels <- c(
  "bio06_1981-2010" = "Bioclim 6",
  "bio04_1981-2010" = "Bioclim 4",
  "bio17_1981-2010" = "Bioclim 17"
)

for (v in sorted_vars) {
  if (!is.na(v)) {
    x_vals_raw <- splineData$x[, v]
    y_vals     <- splineData$y[, v]
    ord        <- order(x_vals_raw)
    
    x_vals_converted <- chelsa_convert(v, x_vals_raw)
    
    plot(x_vals_converted[ord], y_vals[ord], type = "n",
         xlab = if (v %in% names(axis_labels)) axis_labels[[v]] else v,
         ylab = "genomic distance", 
         main = if (v %in% names(title_labels)) title_labels[[v]] else v, 
         ylim = c(0, max_y_val))
    lines(x_vals_converted[ord], y_vals[ord], col = "black", lwd = 1.8)
  }
}
dev.off()

# --- Figure 3: Baseline composition PC1 ---
library(ggplot2)
library(sf)
library(ggspatial)
library(tidyterra)

# 2. Extract non-NA cell coordinates from r_pc1
pc1_df <- as.data.frame(r_pc1, xy = TRUE, na.rm = TRUE)

# 3. Map pixel_colors_pc1 directly to the non-NA grid cells
pc1_df$color_val <- pixel_colors_pc1

# 4. Convert site data to sf object
sites_sf <- st_as_sf(site_data, coords = c("long", "lat"), crs = 4326)

# ---- Figure 3: Baseline Composition PC1 Map ----
p_pc1_map <- ggplot() +
  # 1. PC1 Raster Surface using non-NA coordinates
  geom_raster(data = pc1_df, aes(x = x, y = y, fill = color_val)) +
  scale_fill_identity() +
  
  # 2. State Outlines
  geom_sf(data = aus_states, fill = NA, color = "black", linewidth = 0.4, inherit.aes = FALSE) +
  
  # 3. Sampled Sites Overlay
  geom_sf(
    data = sites_sf,
    shape = 21, fill = "white", color = "black", stroke = 0.35, size = 2.0, alpha = 0.9,
    inherit.aes = FALSE
  ) +
  
  # 4. Extent Alignment (Matches ext_main / limits_pauci)
  coord_sf(
    xlim = c(ext_main["xmin"], ext_main["xmax"]),
    ylim = c(ext_main["ymin"], ext_main["ymax"]),
    expand = FALSE
  ) +
  
  # 5. Scale Bar
  annotation_scale(
    location = "bl",
    width_hint = 0.2,
    style = "ticks",
    text_family = font_family
  ) +
  
  # 6. North Arrow
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
  
  # 7. Theme Formatting
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

# Save high-resolution PNG
png(filename = file.path(out_dir, "fig3_baseline_composition_PC1_final.png"), 
    width = 8.5 * 400, height = 8.5 * 400, res = 400)
print(p_pc1_map)
dev.off()

# 10. Save Outputs ----------------------------------------------------------
write_csv(model_agreement, "climate_variable_model_agreement.csv")

rest_sites <- data.frame(
  site = c("North", "South", "Central"),
  lon = c(151.6856, 148.6164, 149.5631),
  lat = c(-30.7212, -36.3881, -34.2418)
)

save(
  gdm_mod,
  active_env_vars, 
  grid_trans_env, 
  trans_env_126, 
  trans_env_370, 
  trans_env_585,
  grid_126,
  grid_370,
  grid_585,
  rest_sites,
  site_data,
  shared_predictors,
  dist_mat_sub,
  model_agreement,
  file = "gdm_outputs.RData"
)
