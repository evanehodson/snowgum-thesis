library(ggplot2)
library(sf)
library(ggspatial)
library(terra)
library(dplyr)
library(patchwork)

# 1. Register font & load workspace environments
windowsFonts(Times = windowsFont("Times New Roman"))
font_family <- "Times"

# load("gdm_outputs.RData")
# load("lfmm_outputs.RData")

strategy_meta <- data.frame(
  strategy = c("local", "climate", "admixture", "composite", "go", "lfmm"),
  label    = c("Local", "Climate-matched", "Admixture", "Composite", "GDM-Optimized", "LFMM-Optimized"),
  color    = c("#F8766D", "#B79F00", "#00BA38", "#00BFC4", "#619CFF", "#F564E3"),
  stringsAsFactors = FALSE
)

sites_sf <- st_as_sf(env_train, coords = c("long", "lat"), crs = 4326)
if (!"sample_id" %in% names(sites_sf)) {
  sites_sf$sample_id <- env_train$sample_id
}

# Define fixed scenario and sites to iterate over
ssp_val <- "370"
focal_sites <- c("North", "Central", "South")

for (focal_site_name in focal_sites) {
  
  focal_idx <- which(rest_sites$site == focal_site_name)
  if (length(focal_idx) == 0) {
    warning(paste("Site", focal_site_name, "not found in rest_sites. Skipping."))
    next
  }
  focal_row <- rest_sites[focal_idx, ]
  
  dist_km_focal <- haversine_km(focal_row$lon, focal_row$lat, site_data$long, site_data$lat)
  local_pool    <- order(dist_km_focal)[1:K_POOL]
  
  target_raw <- switch(ssp_val,
                       "126" = grid_126_clean[focal_row$idx, clean_predictors, drop = FALSE],
                       "370" = grid_370_clean[focal_row$idx, clean_predictors, drop = FALSE],
                       "585" = grid_585_clean[focal_row$idx, clean_predictors, drop = FALSE]
  )
  target_std <- scale(target_raw, center = pred_means, scale = pred_sds)
  raw_climate_dist <- sqrt(rowSums(sweep(site_covariate_std, 2, as.numeric(target_std), "-")^2))
  climate_pool <- head(order(raw_climate_dist), K_POOL)
  
  admixture_pool <- get_admixture_pool()
  composite_pool <- get_composite_pool(focal_row$lon, focal_row$lat)
  
  target_vec <- switch(ssp_val,
                       "126" = trans_env_126[focal_row$idx, ],
                       "370" = trans_env_370[focal_row$idx, ],
                       "585" = trans_env_585[focal_row$idx, ]
  )
  gdm_offsets_all <- sqrt(rowSums(sweep(site_trans_env_mat_mc, 2, target_vec, "-")^2))
  go_pool <- head(order(gdm_offsets_all), K_POOL)
  
  lfmm_target_offsets <- site_match_scores_mc[[paste(focal_site_name, ssp_val)]]
  lfmm_pool <- head(order(lfmm_target_offsets), K_POOL)
  
  pools_list <- list(
    local     = sample_ids[local_pool],
    climate   = sample_ids[climate_pool],
    admixture = sample_ids[admixture_pool],
    composite = sample_ids[composite_pool],
    go        = sample_ids[go_pool],
    lfmm      = sample_ids[lfmm_pool]
  )
  
  panel_plots <- lapply(seq_len(nrow(strategy_meta)), function(i) {
    strat_key   <- strategy_meta$strategy[i]
    strat_label <- strategy_meta$label[i]
    strat_col   <- strategy_meta$color[i]
    
    selected_ids <- pools_list[[strat_key]]
    sites_sf$selected <- sites_sf$sample_id %in% selected_ids
    
    # Updated for 3 columns x 2 rows indexing
    col_idx <- ((i - 1) %% 3) + 1
    row_idx <- floor((i - 1) / 3) + 1
    
    p <- ggplot() +
      geom_sf(data = aus_states, fill = "grey97", color = "grey70", linewidth = 0.3, inherit.aes = FALSE) +
      geom_sf(
        data = filter(sites_sf, !selected),
        shape = 21, fill = "grey80", color = "grey50", stroke = 0.2, size = 1.2, alpha = 0.5,
        inherit.aes = FALSE
      ) +
      geom_sf(
        data = filter(sites_sf, selected),
        shape = 21, fill = strat_col, color = "black", stroke = 0.35, size = 2.0, alpha = 0.95,
        inherit.aes = FALSE
      ) +
      coord_sf(
        xlim = c(ext_main["xmin"], ext_main["xmax"]),
        ylim = c(ext_main["ymin"], ext_main["ymax"]),
        expand = FALSE
      ) +
      labs(title = strat_label) +
      theme_bw(base_size = 9, base_family = font_family) +
      theme(
        panel.background = element_rect(fill = "white", color = NA),
        plot.background  = element_rect(fill = "white", color = NA),
        panel.grid       = element_blank(),
        panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.7),
        axis.title       = element_blank(),
        plot.title       = element_text(family = font_family, face = "bold", size = 9.5, hjust = 0.5),
        plot.margin      = margin(2, 2, 2, 2)
      )
    
    # X-axis ticks only on the bottom row (Row 2)
    if (row_idx == 2) {
      p <- p + theme(
        axis.text.x = element_text(family = font_family, color = "black", size = 7.5),
        axis.ticks.x = element_line(color = "black")
      )
    } else {
      p <- p + theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
      )
    }
    
    # Y-axis ticks only on the leftmost column (Column 1)
    if (col_idx == 1) {
      p <- p + theme(
        axis.text.y = element_text(family = font_family, color = "black", size = 7.5),
        axis.ticks.y = element_line(color = "black")
      )
    } else {
      p <- p + theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()
      )
    }
    
    # Scale bar in bottom-right panel (Panel 6: col 3, row 2)
    if (i == 6) {
      p <- p + annotation_scale(location = "bl", width_hint = 0.25, style = "ticks", text_family = font_family, text_size = 5.5)
    }
    
    # North arrow in top-left panel (Panel 1: col 1, row 1)
    if (i == 1) {
      p <- p + annotation_north_arrow(
        location = "tl", which_north = "true",
        pad_x = unit(0.12, "in"), pad_y = unit(0.12, "in"),
        height = unit(0.7, "cm"), width = unit(0.7, "cm"),
        style = north_arrow_minimal(fill = "black", line_col = "black", text_family = font_family, text_size = 10)
      )
    }
    p
  })
  
  # Combine into a 3 columns x 2 rows grid
  p_combined_map <- wrap_plots(panel_plots, ncol = 3, nrow = 2)
  
  # Dynamic output filename based on site and fixed ssp370 scenario
  out_filename <- paste0("fig_seed_sourcing_map_", tolower(focal_site_name), "_ssp", ssp_val, ".png")
  ggsave(
    filename = out_filename,
    plot = p_combined_map,
    width = 10.5,   # Adjusted wider for 3 columns
    height = 6.8,   # Adjusted shorter for 2 rows
    dpi = 400
  )
}
