library(tidyr)
library(dplyr)
library(readr)
library(ggplot2)
library(gdm)
library(patchwork)

set.seed(2026)

N_DRAW     <- 10
K_POOL     <- 20
N_ITER     <- 100
N_CLUSTERS <- 3

load("gdm_outputs.RData")
load("lfmm_outputs.RData")

# Replace pooled heterozygosity/inbreeding with within-taxon pauciflora data
het_new <- readr::read_tsv("het_pauciflora.het", show_col_types = FALSE) |>  janitor::clean_names()
het_obs <- setNames(het_new$f, het_new$iid)
length(intersect(names(het_obs), rownames(env_train_mat)))   # expect about 258

haversine_km <- function(lon1, lat1, lon2, lat2) {
  R <- 6371; rad <- pi / 180
  dlat <- (lat2 - lat1) * rad
  dlon <- (lon2 - lon1) * rad
  a <- sin(dlat/2)^2 + cos(lat1*rad)*cos(lat2*rad)*sin(dlon/2)^2
  R * 2 * atan2(sqrt(a), sqrt(1-a))
}

site_trans_env <- gdm.transform(
  gdm_mod, 
  site_data %>% dplyr::select(long, lat, dplyr::all_of(shared_predictors))
)
site_trans_env_mat <- as.matrix(site_trans_env[, active_env_vars, drop = FALSE])
rownames(site_trans_env_mat) <- site_data$sample_id

common_ids <- Reduce(intersect, list(
  site_data$sample_id,
  rownames(env_train_mat),
  names(het_obs)
))
n_common <- length(common_ids)

cat(sprintf(
  "[sample alignment] site_data: %d | lfmm training: %d | heterozygosity: %d | common: %d\n",
  nrow(site_data), nrow(env_train_mat), length(het_obs), n_common
))

sample_ids          <- common_ids
site_data_mc        <- site_data[match(common_ids, site_data$sample_id), ]
site_trans_env_mat_mc <- site_trans_env_mat[match(common_ids, rownames(site_trans_env_mat)), , drop = FALSE]
het_obs_mc          <- het_obs[common_ids]

site_match_scores_mc <- lapply(site_match_scores, function(v) v[common_ids])
stopifnot(!anyNA(unlist(site_match_scores_mc)), !anyNA(het_obs_mc))

# Climate-adjusted
names(site_data_mc)[grep("^bio", names(site_data_mc))] <- gsub("_.*", "", grep("^bio", names(site_data_mc), value = TRUE))

clean_predictors <- colnames(grid_curr_clean)
stopifnot(all(clean_predictors %in% colnames(site_data_mc)))

pred_means         <- colMeans(grid_curr_clean[, clean_predictors])
pred_sds           <- apply(grid_curr_clean[, clean_predictors], 2, sd)
site_covariate_mat <- as.matrix(site_data_mc[, clean_predictors])
site_covariate_std <- scale(site_covariate_mat, center = pred_means, scale = pred_sds)

# Admixture
site_geo_space <- scale(site_data_mc[, c("long", "lat")])
km_res      <- kmeans(site_geo_space, centers = N_CLUSTERS, nstart = 25)
cluster_ids <- km_res$cluster

get_admixture_pool <- function() {
  draws_per_cluster <- ceiling(K_POOL / N_CLUSTERS)
  selected_indices <- c()
  for (cl in unique(cluster_ids)) {
    cluster_members <- which(cluster_ids == cl)
    sampled_members <- sample(cluster_members, size = min(length(cluster_members), draws_per_cluster))
    selected_indices <- c(selected_indices, sampled_members)
  }
  unique(selected_indices)
}

# Composite
get_composite_pool <- function(site_lon, site_lat) {
  dist_km <- haversine_km(site_lon, site_lat, site_data_mc$long, site_data_mc$lat)
  weights <- 1 / (dist_km + 1)
  sample(seq_along(weights), size = K_POOL, prob = weights, replace = FALSE)
}

run_monte_for_site <- function(site_name, site_idx_in_grid, site_lon, site_lat) {
  
  dist_km    <- haversine_km(site_lon, site_lat, site_data_mc$long, site_data_mc$lat)
  local_pool <- order(dist_km)[1:K_POOL]
  resident_idx <- which(dist_km <= 100)
  admixture_pool <- get_admixture_pool()
  composite_pool <- get_composite_pool(site_lon, site_lat)
  
  target_trans_by_ssp <- list(
    "126" = trans_env_126[site_idx_in_grid, ],
    "370" = trans_env_370[site_idx_in_grid, ],
    "585" = trans_env_585[site_idx_in_grid, ]
  )
  target_future_raw_by_ssp <- list(
    "126" = grid_126_clean[site_idx_in_grid, clean_predictors, drop = FALSE],
    "370" = grid_370_clean[site_idx_in_grid, clean_predictors, drop = FALSE],
    "585" = grid_585_clean[site_idx_in_grid, clean_predictors, drop = FALSE]
  )
  
  ssp_scenarios     <- c("126", "370", "585")
  site_results_list <- list()
  
  for (ssp_draw in ssp_scenarios) {
    target_vec <- target_trans_by_ssp[[ssp_draw]]
    target_raw <- target_future_raw_by_ssp[[ssp_draw]]
    
    target_std       <- scale(target_raw, center = pred_means, scale = pred_sds)
    raw_climate_dist <- sqrt(rowSums(sweep(site_covariate_std, 2, as.numeric(target_std), "-")^2))
    climate_pool <- head(order(raw_climate_dist), K_POOL)
    
    gdm_offsets_all <- sqrt(rowSums(sweep(site_trans_env_mat_mc, 2, target_vec, "-")^2))
    go_pool <- head(order(gdm_offsets_all), K_POOL)
    
    lfmm_target_offsets <- site_match_scores_mc[[paste(site_name, ssp_draw)]]
    lfmm_pool <- head(order(lfmm_target_offsets), K_POOL)
    
    pools <- list(local = local_pool, climate = climate_pool, go = go_pool, lfmm = lfmm_pool,
                  admixture = admixture_pool, composite = composite_pool)
    iter_rows <- vector("list", N_ITER)
    
    for (i in seq_len(N_ITER)) {
      draws <- lapply(pools, function(p) sample(p, size = N_DRAW, replace = length(p) < N_DRAW))
      
      strat_rows <- lapply(names(draws), function(strat) {
        idx <- draws[[strat]]
        data.frame(
          iter        = i, site = site_name, ssp = ssp_draw, strategy = strat,
          gdm_offset        = mean(gdm_offsets_all[idx]),
          lfmm_offset       = mean(lfmm_target_offsets[idx]),
          inbreeding        = mean(het_obs_mc[idx], na.rm = TRUE),
          outbreeding_dist  = mean(dist_mat_sub[sample_ids[idx], sample_ids[resident_idx]]),
          stringsAsFactors = FALSE
        )
      })
      iter_rows[[i]] <- bind_rows(strat_rows)
    }
    site_results_list[[ssp_draw]] <- bind_rows(iter_rows)
  }
  
  bind_rows(site_results_list)
}

all_results <- bind_rows(lapply(seq_len(nrow(rest_sites)), function(i) {
  run_monte_for_site(rest_sites$site[i], rest_sites$idx[i], rest_sites$lon[i], rest_sites$lat[i])
}))

metric_labels <- c(
  gdm_offset        = "GDM Offset",
  lfmm_offset       = "LFMM Genetic Gap",
  inbreeding        = "Inbreeding Coefficient (F)",
  outbreeding_dist  = "Outbreeding Distance (to resident population)"
)

plot_data_base <- all_results %>%
  pivot_longer(cols = names(metric_labels), names_to = "metric", values_to = "value") %>%
  mutate(
    strategy = factor(strategy, levels = c("local", "climate", "admixture", "composite", "go", "lfmm"),
                      labels = c("Local", "Climate-matched", "Admixture", "Composite", "GDM-Optimized", "LFMM-Optimized")),
    site     = factor(site, levels = c("North", "Central", "South")),
    metric   = factor(metric, levels = names(metric_labels), labels = unname(metric_labels))
  )

ssp_scenarios <- c("126", "370", "585")
metrics_vec <- unname(metric_labels)

# Loop over each SSP scenario and generate separate files
for (ssp_val in ssp_scenarios) {
  plot_data_ssp <- filter(plot_data_base, ssp == ssp_val)
  
  plot_list <- lapply(seq_along(metrics_vec), function(idx) {
    m <- metrics_vec[idx]
    is_rightmost <- (idx == length(metrics_vec))
    
    p <- ggplot(filter(plot_data_ssp, metric == m), aes(x = strategy, y = value, fill = strategy)) +
      geom_violin(alpha = 0.5, trim = FALSE, scale = "width", width = 0.75) +
      geom_boxplot(width = 0.12, outlier.shape = NA, alpha = 0.8, color = "black") +
      facet_grid(site ~ ., scales = "free_y") +
      labs(
        title = m,
        x = NULL, y = NULL
      ) +
      theme_bw(base_size = 9) +
      theme(
        legend.position = "none",
        plot.title = element_text(face = "bold", size = 10, hjust = 0.5),
        axis.text.x = element_text(angle = 35, hjust = 1, size = 7.5),
        panel.grid.minor = element_blank()
      )
    
    if (!is_rightmost) {
      p <- p + theme(
        strip.text.y = element_blank(),
        strip.background = element_blank()
      )
    } else {
      p <- p + theme(
        strip.text.y = element_text(face = "bold", size = 8),
        strip.background = element_rect(fill = "grey90", color = "grey50")
      )
    }
    p
  })
  
  p_combined <- wrap_plots(plot_list, nrow = 1) + 
    plot_annotation(
      title = paste0("Seed Sourcing Strategy Comparison (SSP", ssp_val, ")"),
      subtitle = "Evaluation across independent criteria (columns) and restoration sites (rows)",
      theme = theme(
        plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
        plot.subtitle = element_text(size = 9.5, hjust = 0.5)
      )
    )
  
  legend_plot <- ggplot(filter(plot_data_ssp, metric == metrics_vec[1]), aes(x = strategy, y = value, fill = strategy)) +
    geom_boxplot() +
    theme_bw() +
    theme(legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = 8.5))
  shared_legend <- cowplot::get_legend(legend_plot)
  
  p_final <- p_combined / shared_legend + plot_layout(heights = c(12, 1))
  
  filename <- paste0("monte_carlo_comparison_ssp", ssp_val, ".png")
  ggsave(filename, plot = p_final, width = 14, height = 6.5, dpi = 300)
}

write_csv(all_results, "monte_carlo_results.csv")