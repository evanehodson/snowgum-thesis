library(dplyr)
library(matrixStats)
library(LEA)
library(terra)
library(stringr)

options(scipen = 999)
K_LATENT <- 9

# 1. Load Genotypes & Sample IDs ------------------------------------------
dt <- read.delim("candidates_geno.raw", sep = "\t", header = TRUE)
geno_mat <- as.matrix(dt[, -(1:6)])
rownames(geno_mat) <- dt$IID

pauci_samples <- readLines("pauci_samples.txt")

# 2. Load Environment Data ------------------------------------------------
env_present     <- read.csv("current_covariates.csv")
grid_current    <- read.csv("current_covariate_grid.csv")
grid_future_raw <- read.csv("future_covariate_grid.csv")

names(env_present)[grep("^bio", names(env_present))] <- gsub("_.*", "", grep("^bio", names(env_present), value = TRUE))
names(grid_current)[grep("^bio", names(grid_current))] <- gsub("_.*", "", grep("^bio", names(grid_current), value = TRUE))

target_vars <- grep("^bio", names(grid_current), value = TRUE)

average_ssp_grid <- function(future_df, vars_clean, ssp_target) {
  res <- data.frame(grid_id = future_df$grid_id, lon = future_df$lon, lat = future_df$lat)
  for (v in vars_clean) {
    match_cols <- grep(paste0("^", v, "_.*_", ssp_target, "$"), names(future_df), value = TRUE)
    res[[v]] <- rowMeans(as.matrix(future_df[, match_cols]), na.rm = TRUE)
  }
  return(res)
}

grid_126 <- average_ssp_grid(grid_future_raw, target_vars, "ssp126")
grid_370 <- average_ssp_grid(grid_future_raw, target_vars, "ssp370")
grid_585 <- average_ssp_grid(grid_future_raw, target_vars, "ssp585")


# genetic.gap() stops on any NA, so cells must be complete in current AND all futures
complete_idx <- complete.cases(grid_current[, target_vars]) &
  complete.cases(grid_126[, target_vars]) &
  complete.cases(grid_370[, target_vars]) &
  complete.cases(grid_585[, target_vars])

grid_curr_clean <- as.matrix(grid_current[complete_idx, target_vars])
grid_126_clean  <- as.matrix(grid_126[complete_idx, target_vars])
grid_370_clean  <- as.matrix(grid_370[complete_idx, target_vars])
grid_585_clean  <- as.matrix(grid_585[complete_idx, target_vars])

future_grids <- list("126" = grid_126_clean, "370" = grid_370_clean, "585" = grid_585_clean)


# Align genotypes to available environmental samples
common_samples <- intersect(intersect(pauci_samples, rownames(geno_mat)), env_present$sample_id)
geno_mat       <- geno_mat[common_samples, ]
env_train      <- env_present[match(common_samples, env_present$sample_id), ]
stopifnot(identical(rownames(geno_mat), env_train$sample_id))


# Imputation
col_means <- round(colMeans2(geno_mat, na.rm = TRUE))
for (j in seq_len(ncol(geno_mat))) {
  geno_mat[is.na(geno_mat[, j]), j] <- col_means[j]
}


# Scale training environment
env_train_mat <- as.matrix(env_train[, target_vars])
rownames(env_train_mat) <- env_train$sample_id
train_means   <- colMeans(env_train_mat)
train_sds     <- apply(env_train_mat, 2, sd)
stopifnot(all(train_sds > 0))

scale_env <- function(m) scale(as.matrix(m), center = train_means, scale = train_sds)


# 3. Fit LFMM2 Model ------------------------------------------------------
set.seed(1)
mod_lfmm <- lfmm2(input = geno_mat, env = scale_env(env_train_mat), K = K_LATENT, effect.sizes = TRUE)
B_eff    <- as.matrix(mod_lfmm@B)

rownames(B_eff) <- colnames(geno_mat)


# Genetic gap in the training-scaled environmental space (Gain et al. 2023).
# Confirmed against LEA::genetic.gap():
#   gg = rowSums(((X.new - X.pred) %*% t(B))^2) / nrow(B)
# which is algebraically identical to rowSums((delta %*% crossprod(B)) * delta) / nrow(B).
# This reformulation avoids ever forming the full pixels x loci matrix
# (tens of GB here) since BtB is only d x d (d = number of env variables).
BtB <- crossprod(B_eff)

gap_from_B <- function(new_env, pred_env) {
  delta <- scale_env(new_env) - scale_env(pred_env)
  as.vector(rowSums((delta %*% BtB) * delta)) / nrow(B_eff)
}


# 4. Temporal Offsets (current -> future, per pixel) ----------------------
# Check against genetic.gap() on a pixel subsample.
sub_idx <- seq_len(min(2000, nrow(grid_curr_clean)))

set.seed(1)
gap_canonical <- genetic.gap(
  input    = geno_mat,
  env      = env_train_mat,
  new.env  = grid_curr_clean[sub_idx, , drop = FALSE],
  pred.env = future_grids[["126"]][sub_idx, , drop = FALSE],
  K        = K_LATENT,
  scale    = TRUE
)$offset

gap_check <- gap_from_B(grid_curr_clean[sub_idx, , drop = FALSE],
                        future_grids[["126"]][sub_idx, , drop = FALSE])

max_abs_diff <- max(abs(gap_canonical - gap_check))
max_rel_diff <- max(abs(gap_canonical - gap_check) / pmax(abs(gap_canonical), 1e-12))
gap_cor      <- cor(gap_canonical, gap_check)

cat(sprintf("[genetic.gap validation] max abs diff = %.3g | max rel diff = %.3g | correlation = %.8f\n",
            max_abs_diff, max_rel_diff, gap_cor))

if (gap_cor < 0.999) {
  stop("gap_from_B and genetic.gap() disagree beyond floating-point noise (correlation < 0.999). ",
       "Do not proceed to the full 400k-pixel run until this is resolved.")
} else if (max_abs_diff > 1e-3) {
  warning("gap_from_B matches genetic.gap() closely (correlation > 0.999) but with larger-than-expected ",
          "absolute differences. Likely SVD/factor-estimation noise in lfmm2() at this sample size -- ",
          "safe to proceed, but worth a second look if downstream offset magnitudes look off.")
} else {
  cat("[genetic.gap validation] OK -- gap_from_B matches canonical genetic.gap() within floating-point noise.\n")
}


offset_126 <- gap_from_B(grid_curr_clean, future_grids[["126"]])
offset_370 <- gap_from_B(grid_curr_clean, future_grids[["370"]])
offset_585 <- gap_from_B(grid_curr_clean, future_grids[["585"]])


# 4b. Heterozygosity (inbreeding-depression proxy)
het_file_path <- "pauci_stellulata_het.het"

if (file.exists(het_file_path)) {
  het_raw <- readr::read_tsv(het_file_path, show_col_types = FALSE) |>
    janitor::clean_names() |>
    dplyr::mutate(iid = str_remove(iid, "^(.*)_\\1$"))
  het_filtered <- het_raw |>
    dplyr::filter(iid %in% common_samples)
  f_obs_temp <- setNames(het_filtered$f, het_filtered$iid)
  het_obs <- f_obs_temp[common_samples]
} else {
  warning("Heterozygosity file not found at '", het_file_path,
          "' -- het_obs is a placeholder of NA for all common_samples. ",
          "Monte Carlo heterozygosity comparisons will be meaningless until this is real.")
  het_obs <- setNames(rep(NA_real_, length(common_samples)), common_samples)
}


# 5. Build Terra Template Raster & Spatial Mapping ------------------------
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

rast_rgb <- function(score_mat, coords) {
  xyz <- data.frame(x = coords$lon, y = coords$lat, score_mat)
  names(xyz) <- c("x", "y", paste0("PC", seq_len(ncol(score_mat))))
  pts <- vect(xyz, geom = c("x", "y"))
  rasterize(pts, r_template, field = paste0("PC", seq_len(ncol(score_mat))))
}


base_coords <- grid_current[complete_idx, c("lon", "lat")]


# 6. Baseline PCA of Predicted Current Composition ------------------------
curr_scaled <- scale_env(grid_curr_clean)
eig_BtB <- eigen(BtB)
eig_BtB$values[eig_BtB$values < 0] <- 0
pca_input <- curr_scaled %*% eig_BtB$vectors %*% diag(sqrt(eig_BtB$values))
pca_fit         <- prcomp(pca_input, center = TRUE, scale. = FALSE)
pca_scores_curr <- pca_fit$x[, 1:3, drop = FALSE]

pc1_vals <- pca_scores_curr[, 1]
pc1_norm <- (pc1_vals - min(pc1_vals)) / (max(pc1_vals) - min(pc1_vals))  # true linear 0-1
hue_vals <- pc1_norm * 270

pixel_colors_pc1 <- hsv(h = hue_vals / 360, s = 0.75, v = 0.85)

xyz_pc1 <- data.frame(x = base_coords$lon, y = base_coords$lat, col_id = seq_along(pixel_colors_pc1))
pts_pc1 <- vect(xyz_pc1, geom = c("x", "y"))
r_pc1   <- rasterize(pts_pc1, r_template, field = "col_id")


# 7. Restoration Seed Sourcing (spatiotemporal / forward-looking offsets) ----
rest_sites <- data.frame(
  site = c("North", "South", "Central"),
  lon = c(151.6856, 148.6164, 149.5631),
  lat = c(-30.7212, -36.3881, -34.2418)
)

valid_idx <- which(complete_idx)
rest_sites$idx <- sapply(1:nrow(rest_sites), function(i) {
  d <- (grid_current$lon[valid_idx] - rest_sites$lon[i])^2 + (grid_current$lat[valid_idx] - rest_sites$lat[i])^2
  valid_idx[which.min(d)]
})

match_rasters     <- list()
site_match_scores <- list()

for (i in 1:nrow(rest_sites)) {
  s_idx  <- rest_sites$idx[i]
  s_name <- rest_sites$site[i]
  
  for (ssp in names(future_grids)) {
    # 1. Get the restoration site's FUTURE climate under this specific SSP
    fut_grid_ssp <- future_grids[[ssp]]
    target_env   <- fut_grid_ssp[s_idx, ]
    
    # 2. Create a reference matrix where every row is the target site's future climate
    pred_rep_pix <- matrix(target_env, nrow = nrow(grid_curr_clean), ncol = length(target_env), byrow = TRUE)
    
    # 3. Calculate the genetic gap from each pixel's CURRENT climate to the site's FUTURE climate
    dist_vec <- gap_from_B(new_env = grid_curr_clean, pred_env = pred_rep_pix)
    
    match_rasters[[paste(s_name, ssp)]] <- rast_single(dist_vec, base_coords)
    
    # Per-training-sample gaps (training environments today vs. site future climate)
    pred_rep_sites <- matrix(target_env, nrow = nrow(env_train_mat), ncol = length(target_env), byrow = TRUE)
    site_match_scores[[paste(s_name, ssp)]] <- setNames(
      gap_from_B(new_env = env_train_mat, pred_env = pred_rep_sites),
      rownames(env_train_mat)
    )
  }
}


# 8. Compile PDF Visual Outputs -------------------------------------------
pdf("lea_genomic_offsets_and_future.pdf", width = 14, height = 10)


pal_heat  <- colorRampPalette(c("lightgreen", "yellow", "orange", "red", "darkred"))(100)
pal_match <- colorRampPalette(c("darkgreen", "lightgreen", "yellow", "orange", "red", "darkred"))(100)

global_min_off <- min(c(offset_126, offset_370, offset_585), na.rm = TRUE)
global_max_off <- max(c(offset_126, offset_370, offset_585), na.rm = TRUE)

all_match_vals  <- unlist(lapply(match_rasters, values))
global_site_min <- min(all_match_vals, na.rm = TRUE)
global_site_max <- max(all_match_vals, na.rm = TRUE)


par(mfrow = c(2, 2), oma = c(0, 0, 2, 0))


par(mar = c(2, 2, 3, 1))
plot(r_pc1, col = pixel_colors_pc1, legend = FALSE, axes = FALSE, box = FALSE)
points(env_train$long, env_train$lat, col = "black", pch = 20, cex = 0.5)
title(main = "A) Baseline Current Adaptive Composition", line = -1, adj = 0.5, cex.main = 1.2)


plot_off_panel <- function(r, panel_title) {
  par(mar = c(2, 2, 3, 4))
  plot(r, col = pal_heat, range = c(global_min_off, global_max_off), main = panel_title, axes = FALSE)
}

r_off_126 <- rast_single(offset_126, base_coords)
r_off_370 <- rast_single(offset_370, base_coords)
r_off_585 <- rast_single(offset_585, base_coords)

plot_off_panel(r_off_126, "B) LEA Genomic Offset: SSP1-2.6")
plot_off_panel(r_off_370, "C) LEA Genomic Offset: SSP3-7.0")
plot_off_panel(r_off_585, "D) LEA Genomic Offset: SSP5-8.5")

mtext("LEA Genomic Offset Vulnerability Maps", outer = TRUE, cex = 1.5, font = 2)


for (i in 1:nrow(rest_sites)) {
  s_name <- rest_sites$site[i]
  s_lon  <- rest_sites$lon[i]
  s_lat  <- rest_sites$lat[i]
  par(mfrow = c(2, 2), mar = c(2, 2, 3, 4), oma = c(0, 0, 2, 0))
  for (ssp in names(future_grids)) {
    r_match <- match_rasters[[paste(s_name, ssp)]]
    plot(r_match, col = pal_match, range = c(global_site_min, global_site_max),
         main = paste("SSP", ssp), axes = FALSE)
    points(s_lon, s_lat, col = "black", bg = "white", pch = 21, cex = 1.5, lwd = 1.5)
  }
  mtext(paste("Seed Sourcing Match:", s_name), outer = TRUE, cex = 1.5, font = 2)
}

dev.off()


save(
  mod_lfmm, 
  B_eff, 
  BtB, 
  offset_126, 
  offset_370, 
  offset_585,
  rest_sites, 
  match_rasters, 
  site_match_scores,
  env_train_mat,
  grid_curr_clean,
  grid_126_clean,
  grid_370_clean,
  grid_585_clean,
  het_obs,
  file = "lfmm_outputs.RData"
)