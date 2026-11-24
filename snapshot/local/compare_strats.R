library(grid)
library(png)
suppressMessages(library(httr2))

set.seed(42)

# ---- Geometry & Spatial Setup ------------------------------------------------
RES <- 250

x_vec <- seq(0, 1, length.out = RES)
y_vec <- seq(1, 0, length.out = RES)

X <- matrix(rep(x_vec, times = RES), nrow = RES, ncol = RES, byrow = TRUE)
Y <- matrix(rep(y_vec, times = RES), nrow = RES, ncol = RES, byrow = FALSE)

site_xy <- c(0.50, 0.48)

# ---- Fetch real satellite imagery (Esri World Imagery) ----------------------
kosci_bbox <- c(xmin = 148.28, ymin = -36.48, xmax = 148.42, ymax = -36.38)

fetch_satellite_png <- function(bbox, out_path, size_px = 800) {
  bbox_str <- paste(bbox["xmin"], bbox["ymin"], bbox["xmax"], bbox["ymax"], sep = ",")
  
  req <- request("https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/export") |>
    req_url_query(
      bbox     = bbox_str,
      bboxSR   = 4326,
      imageSR  = 4326,
      size     = paste0(size_px, ",", size_px),
      format   = "png32",
      f        = "image"
    ) |>
    req_timeout(60)
  
  resp <- req_perform(req)
  writeBin(resp_body_raw(resp), out_path)
  invisible(out_path)
}

sat_dir <- tempdir()
sat_path <- file.path(sat_dir, "sat_tile1.png")
fetch_satellite_png(kosci_bbox, sat_path)
tile1 <- readPNG(sat_path)

sat_rasters <- list(Local = tile1, Composite = tile1, Admixture = tile1)

# ---- Fractal noise field generator -------------------------------------------
upsample_bilinear <- function(coarse_mat, RES) {
  n <- nrow(coarse_mat)
  coarse_coords <- seq(0, 1, length.out = n)
  fine_coords   <- seq(0, 1, length.out = RES)
  
  tmp <- matrix(0, nrow = n, ncol = RES)
  for (i in 1:n) {
    tmp[i, ] <- approx(coarse_coords, coarse_mat[i, ], xout = fine_coords, rule = 2)$y
  }
  out <- matrix(0, nrow = RES, ncol = RES)
  for (j in 1:RES) {
    out[, j] <- approx(coarse_coords, tmp[, j], xout = fine_coords, rule = 2)$y
  }
  out
}

fractal_noise <- function(RES, seed, octaves = 2, base_freq = 3,
                          persistence = 0.55, lacunarity = 2.1) {
  set.seed(seed)
  Z <- matrix(0, nrow = RES, ncol = RES)
  amplitude <- 1
  freq <- base_freq
  total_amp <- 0
  
  for (o in 1:octaves) {
    n <- max(2, round(freq))
    coarse <- matrix(rnorm(n * n), nrow = n, ncol = n)
    layer <- upsample_bilinear(coarse, RES)
    Z <- Z + amplitude * layer
    
    total_amp <- total_amp + amplitude
    amplitude <- amplitude * persistence
    freq <- freq * lacunarity
  }
  
  Z <- Z / total_amp
  (Z - min(Z)) / (max(Z) - min(Z))
}

env_temp   <- fractal_noise(RES, seed = 101, octaves = 2, base_freq = 3, persistence = 0.6, lacunarity = 1.0)
env_precip <- fractal_noise(RES, seed = 202, octaves = 2, base_freq = 4, persistence = 0.58, lacunarity = 1.0)

fine_noise_1 <- fractal_noise(RES, seed = 303, octaves = 2, base_freq = 5, persistence = 0.55, lacunarity = 1.0)
fine_noise_2 <- fractal_noise(RES, seed = 404, octaves = 2, base_freq = 6, persistence = 0.5, lacunarity = 1.0)
fine_noise_3 <- fractal_noise(RES, seed = 505, octaves = 2, base_freq = 7, persistence = 0.5, lacunarity = 1.0)

# ---- Complex Organic Population Boundaries -----------------------------------
make_complex_blob <- function(cx, cy, r_base, seed) {
  set.seed(seed)
  n_pts <- 16
  angles <- seq(0, 2 * pi, length.out = n_pts + 1)[-(n_pts + 1)]
  
  radii <- r_base + runif(n_pts, -r_base * 0.3, r_base * 0.5)
  radii <- (radii + c(radii[-1], radii[1]) + c(radii[n_pts], radii[1:(n_pts - 1)])) / 3
  
  s <- spline(c(angles, angles[1] + 2 * pi), c(radii, radii[1]), n = 150)
  
  list(x = cx + s$y[1:150] * cos(s$x[1:150]),
       y = cy + s$y[1:150] * sin(s$x[1:150]) * 0.85)
}

pop_A <- make_complex_blob(0.15, 0.62, 0.11, seed = 11)
pop_B <- make_complex_blob(0.50, 0.48, 0.12, seed = 22)
pop_C <- make_complex_blob(0.85, 0.38, 0.10, seed = 33)
populations <- list(A = pop_A, B = pop_B, C = pop_C)

# ---- Seed Sources (Coordinates) ---------------------------------------------
make_sources <- function(strategy) {
  switch(strategy,
         # Local: 6 donors in central range (minimum distance from site enforced)
         "Local" = data.frame(
           x = site_xy[1] + c(-0.045, 0.045, -0.035, 0.035, -0.015, 0.025),
           y = site_xy[2] + c(0.030,  0.035, -0.040, -0.035,  0.050, -0.050)
         ),
         # Composite: 2 in left, 4 in central, 2 in right
         "Composite" = data.frame(
           x = c(0.13, 0.16, site_xy[1] + c(-0.04, 0.04, -0.03, 0.03), 0.84, 0.87),
           y = c(0.63, 0.59, site_xy[2] + c(0.03, 0.03, -0.03, -0.03), 0.39, 0.36)
         ),
         # Admixture: 3 in each range
         "Admixture" = data.frame(
           x = c(0.13, 0.16, 0.18, site_xy[1] + c(-0.04, 0.04, 0.00), 0.83, 0.85, 0.88),
           y = c(0.65, 0.60, 0.64, site_xy[2] + c(0.03, -0.03, 0.045), 0.41, 0.35, 0.39)
         ),
         "Climate matching" = data.frame(
           x = c(0.13, 0.15, 0.17, 0.14, 0.16),
           y = c(0.60, 0.65, 0.58, 0.66, 0.62)
         ),
         "GDM-optimized" = data.frame(
           x = c(0.15, 0.17, 0.13, site_xy[1] - 0.035, site_xy[1] + 0.035, site_xy[1] + 0.04),
           y = c(0.60, 0.65, 0.58, site_xy[2] + 0.035, site_xy[2] - 0.035, site_xy[2] + 0.03)
         ),
         # LFMM-optimized: 5 in left, 1 in right
         "LFMM-optimized" = data.frame(
           x = c(0.13, 0.15, 0.17, 0.14, 0.16, 0.85),
           y = c(0.62, 0.66, 0.64, 0.58, 0.60, 0.38)
         )
  )
}

# ---- Space Transforms (Bottom row) -------------------------------------------
normalize <- function(mat) { (mat - min(mat)) / (max(mat) - min(mat)) }

t_target <- env_temp[round(0.38 * RES), round(0.15 * RES)]
p_target <- env_precip[round(0.38 * RES), round(0.15 * RES)]
clim_dist <- sqrt((env_temp - t_target)^2 + (env_precip - p_target)^2)
clim_dist <- normalize(clim_dist)
clim_dist <- normalize(clim_dist * 0.65 + fine_noise_1 * 0.2 + fine_noise_2 * 0.15)

gdm_t <- 1 / (1 + exp(-8 * (env_temp - 0.5)))
gdm_p <- sin(env_precip * 1.5 * pi)
gt_target <- gdm_t[round(0.52 * RES), round(0.50 * RES)]
gp_target <- gdm_p[round(0.52 * RES), round(0.50 * RES)]
gdm_dist <- sqrt((gdm_t - gt_target)^2 + (gdm_p - gp_target)^2)
gdm_dist <- normalize(gdm_dist)
gdm_dist <- normalize(gdm_dist * 0.6 + fine_noise_2 * 0.2 + fine_noise_3 * 0.2)

hotspot <- exp(-((X - 0.15)^2 + (Y - 0.62)^2) / 0.02)
lfmm_dist <- normalize(1 - hotspot)
lfmm_dist <- normalize(lfmm_dist * 0.55 + fine_noise_1 * 0.25 + fine_noise_3 * 0.2)

# ---- Color Palette -----------------------------------------------------------
dist_ramp <- colorRamp(c("#FFFFFF", "#FEE08B", "#F46D43", "#B2182B", "#6B1119"))

matrix_to_hex <- function(mat, ramp_func) {
  rgb_mat <- ramp_func(as.vector(mat))
  hex_vec <- rgb(rgb_mat[, 1], rgb_mat[, 2], rgb_mat[, 3], maxColorValue = 255)
  matrix(hex_vec, nrow = RES, ncol = RES)
}

# ---- Panel Rendering ---------------------------------------------------------
draw_panel <- function(strategy) {
  grid.rect(gp = gpar(fill = "#FFFFFF", col = "#CCCCCC", lwd = 1))
  
  grid.text(strategy, x = 0.04, y = 0.92, just = "left",
            gp = gpar(fontsize = 12, fontface = "bold", col = "#111111", fontfamily = "sans"))
  
  pushViewport(viewport(x = 0.5, y = 0.43, width = 0.92, height = 0.78))
  
  if (strategy %in% names(sat_rasters)) {
    grid.raster(sat_rasters[[strategy]], width = 1, height = 1)
    grid.rect(gp = gpar(fill = "#FFFFFF", alpha = 0.12, col = NA))
  } else if (strategy == "Climate matching") {
    grid.raster(matrix_to_hex(clim_dist, dist_ramp), width = 1, height = 1)
  } else if (strategy == "GDM-optimized") {
    grid.raster(matrix_to_hex(gdm_dist, dist_ramp), width = 1, height = 1)
  } else if (strategy == "LFMM-optimized") {
    grid.raster(matrix_to_hex(lfmm_dist, dist_ramp), width = 1, height = 1)
  }
  
  grid.rect(gp = gpar(fill = NA, col = "#888888", lwd = 0.8))
  
  # Thicker single dashed boundary line without white outline
  for (p in populations) {
    grid.polygon(p$x, p$y, gp = gpar(fill = NA, col = "#111111", lwd = 1.8, lty = "dashed"))
  }
  
  pts <- make_sources(strategy)
  
  # Seed donor points (no extra outline)
  grid.points(pts$x, pts$y, pch = 21,
              gp = gpar(fill = "#FFFFFF", col = "#000000", lwd = 1.2),
              size = unit(0.13, "inches"))
  
  # Restoration site point
  grid.points(site_xy[1], site_xy[2], pch = 23,
              gp = gpar(fill = "#FF3311", col = "#000000", lwd = 1.2),
              size = unit(0.19, "inches"))
  
  popViewport()
}

# ---- Main Layout Assembly ----------------------------------------------------
render_figure <- function() {
  grid.newpage()
  grid.rect(gp = gpar(fill = "#FAFAFA", col = NA))
  
  pushViewport(viewport(y = 0.12, height = 0.85, just = "bottom",
                        layout = grid.layout(2, 3,
                                             widths = unit(rep(1, 3), "null"),
                                             heights = unit(rep(1, 2), "null"))))
  
  strats <- c("Local", "Composite", "Admixture", "Climate matching", "GDM-optimized", "LFMM-optimized")
  
  for (idx in seq_along(strats)) {
    row <- ifelse(idx <= 3, 1, 2)
    col <- ifelse(idx %% 3 == 0, 3, idx %% 3)
    
    pushViewport(viewport(layout.pos.row = row, layout.pos.col = col))
    pushViewport(viewport(width = 0.95, height = 0.95))
    draw_panel(strats[idx])
    popViewport(2)
  }
  popViewport()
  
  # ---- CLEAN LEGEND ----
  pushViewport(viewport(y = 0, height = 0.12, just = "bottom"))
  
  grid.points(0.22, 0.5, pch = 23, gp = gpar(fill = "#FF3311", col = "#000000", lwd = 1.2), size = unit(0.17, "inches"))
  grid.text("Restoration site", x = 0.24, y = 0.5, just = "left", gp = gpar(fontsize = 11, col = "#111111"))
  
  grid.points(0.44, 0.5, pch = 21, gp = gpar(fill = "#FFFFFF", col = "#000000", lwd = 1.2), size = unit(0.13, "inches"))
  grid.text("Seed donor", x = 0.46, y = 0.5, just = "left", gp = gpar(fontsize = 11, col = "#111111"))
  
  grid.lines(x = c(0.62, 0.65), y = c(0.5, 0.5), gp = gpar(col = "#111111", lwd = 1.8, lty = "dashed"))
  grid.text("Population bounds", x = 0.66, y = 0.5, just = "left", gp = gpar(fontsize = 11, col = "#111111"))
  
  grid.text("Satellite imagery: Esri, Maxar, Earthstar Geographics (World Imagery)",
            x = 0.98, y = 0.08, just = "right", gp = gpar(fontsize = 7, col = "#999999"))
  
  popViewport()
}

# ---- Render High-Res PNG -----------------------------------------------------
png("seed_sourcing_strategies_v14.png", width = 11, height = 7.5, units = "in", res = 300)
render_figure()
dev.off()
