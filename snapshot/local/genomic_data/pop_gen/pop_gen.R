library(SNPRelate)
library(tidyverse)
library(patchwork)
library(janitor)
library(geosphere)
library(vegan)
library(dendextend)
library(ggdendro)

species_colors <- c(
  "Eucalyptus pauciflora subsp. niphophila"   = "#56B4E9", 
  "Eucalyptus pauciflora subsp. pauciflora"   = "#E69F00", 
  "Eucalyptus stellulata"                    = "#CC79A7"  
)

abbrev_species <- function(x) str_replace(x, "^Eucalyptus\\s+", "E. ")

euc_meta <- read_csv("https://docs.google.com/spreadsheets/d/1QI4KE4M7i8ST6mLNoG8rkkfkUio_3AGQRQY6PQ3Sg24/export?gid=4473427&format=csv")

euc_meta <- euc_meta %>%
  mutate(fastq_sample_name = paste(fastq_sample_name, fastq_sample_name, sep = "_"))

pruned_gds_path <- "snowgums_200k_pruned.gds"
if (!file.exists(pruned_gds_path)) {
  snpgdsVCF2GDS("snowgums_200k_pruned.vcf.gz", pruned_gds_path)
}
pruned_gds <- snpgdsOpen(pruned_gds_path)

pruned_sampinfo <- snpgdsSummary(pruned_gds)
initial_good_samples <- pruned_sampinfo$sample.id[
  pruned_sampinfo$sample.id %in% euc_meta$fastq_sample_name
]

meta <- euc_meta[match(initial_good_samples, euc_meta$fastq_sample_name), ] %>%
  mutate(species_clean = species %>% str_replace("X.*", "") %>% str_replace(" x .*", ""))

# Kinship (KING-robust)
king <- snpgdsIBDKING(pruned_gds, sample.id = initial_good_samples, autosome.only = FALSE, verbose = FALSE)
kinship_mat <- king$kinship
rownames(kinship_mat) <- colnames(kinship_mat) <- king$sample.id

kinship_long <- kinship_mat %>%
  as.data.frame() %>%
  rownames_to_column("sample_1") %>%
  pivot_longer(-sample_1, names_to = "sample_2", values_to = "kinship") %>%
  filter(sample_1 != sample_2) %>%
  left_join(meta %>% select(sample_id = fastq_sample_name, species_clean), by = join_by(sample_1 == sample_id)) %>%
  left_join(meta %>% select(sample_id = fastq_sample_name, species_clean), by = join_by(sample_2 == sample_id), suffix = c("_1", "_2"))

# flag pairs at or above the first-degree-relative threshold (kinship > 0.177)
close_relatives <- kinship_long %>%
  filter(sample_1 < sample_2, kinship > 0.177) %>%
  arrange(desc(kinship))

print(close_relatives)
write_csv(close_relatives, "king_close_relative_pairs.csv")

# per-sample missingness used to decide which member of a close-relative pair to drop
pruned_missrate <- snpgdsSampMissRate(pruned_gds, sample.id = initial_good_samples)
names(pruned_missrate) <- initial_good_samples

related_exclusions <- close_relatives %>%
  mutate(
    miss_1 = pruned_missrate[sample_1],
    miss_2 = pruned_missrate[sample_2],
    exclude = if_else(miss_1 >= miss_2, sample_1, sample_2)
  ) %>%
  pull(exclude) %>%
  unique()

initial_good_samples <- setdiff(initial_good_samples, related_exclusions)
meta <- meta %>% filter(!fastq_sample_name %in% related_exclusions)

gds <- pruned_gds

# PCA
pca_to_df <- function(pca) {
  pcadf <- pca$eigenvect %>% as.data.frame()
  colnames(pcadf) <- sprintf("PC%02d", 1:ncol(pcadf))
  pcadf <- pcadf %>% mutate(sample_id = as.character(pca$sample.id))
  return(pcadf)
}

run_iterative_pca_cleanup <- function(gds, initial_samples, filter_strings) {
  current_samples <- initial_samples
  for (filter_str in filter_strings) {
    pca <- snpgdsPCA(gds, sample.id = current_samples, autosome.only = FALSE, verbose = FALSE)
    df  <- pca_to_df(pca)
    condition_expr <- parse(text = filter_str)
    outliers       <- pca$sample.id[eval(condition_expr, envir = df)]
    current_samples <- setdiff(current_samples, outliers)
  }
  final_pca <- snpgdsPCA(gds, sample.id = current_samples, autosome.only = FALSE, verbose = FALSE)
  final_df  <- pca_to_df(final_pca)
  return(list(pca = final_pca, df = final_df))
}

global_filters <- c(
  "PC01 < -0.2"
)

pca_results <- run_iterative_pca_cleanup(gds, initial_good_samples, global_filters)

pca_data <- pca_results$df %>%
  left_join(meta, by = join_by(sample_id == fastq_sample_name))

plot_pca <- ggplot(pca_data, aes(x = PC01, y = PC02, color = species_clean)) +
  geom_point(size = 2.5, alpha = 0.85) +
  scale_colour_manual(values = species_colors) +
  scale_x_continuous(breaks = seq(-1, 1, 0.1)) +
  scale_y_continuous(breaks = seq(-1, 1, 0.1)) +
  labs(
    x = paste0("PC01 (", round(pca_results$pca$varprop[1] * 100, 1), "%)"),
    y = paste0("PC02 (", round(pca_results$pca$varprop[2] * 100, 1), "%)"),
    color = "Taxa"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 8),
    plot.title = element_text(face = "bold")
  )

print(plot_pca)
ggsave("clean_global_pca.png", plot = plot_pca, width = 8, height = 6.5, dpi = 600)

writeLines(pca_results$pca$sample.id, "final_pca_samples.txt")

clean_sample_ids <- pca_data$sample_id

# IBS Outlier
technical_outlier_exclusions <- c("MON0667_MON0667")

clean_sample_ids <- setdiff(clean_sample_ids, technical_outlier_exclusions)
pca_data <- pca_data %>% filter(!sample_id %in% technical_outlier_exclusions)

# Dendrogram
ibs <- snpgdsIBS(gds, sample.id = clean_sample_ids, autosome.only = FALSE, verbose = FALSE)
ibs_mat <- 1 - ibs$ibs
rownames(ibs_mat) <- colnames(ibs_mat) <- ibs$sample.id
dist_mat <- as.dist(ibs_mat)
hc <- hclust(dist_mat)
species_lookup <- setNames(pca_data$species_clean, pca_data$sample_id)
dend <- as.dendrogram(hc)

collapse_pure_subtrees <- function(dend, species_lookup, min_size = 2) {
  collapse_node <- function(node) {
    labs <- labels(node)
    node_species <- species_lookup[labs]
    
    if (length(labs) >= min_size && length(unique(node_species)) == 1) {
      taxon <- unique(node_species)[1]
      
      new_leaf <- labs[1]
      class(new_leaf) <- "dendrogram"
      attr(new_leaf, "leaf") <- TRUE
      attr(new_leaf, "members") <- 1L
      attr(new_leaf, "height") <- 0
      attr(new_leaf, "label") <- labs[1]
      attr(new_leaf, "species") <- taxon
      attr(new_leaf, "n_members") <- length(labs)
      attr(new_leaf, "is_collapsed") <- TRUE
      return(new_leaf)
    }
    
    if (!is.leaf(node)) {
      for (i in seq_along(node)) {
        node[[i]] <- collapse_node(node[[i]])
      }
    }
    return(node)
  }
  collapse_node(dend)
}

dend_collapsed <- collapse_pure_subtrees(dend, species_lookup, min_size = 2)

tag_singletons <- function(node, species_lookup) {
  if (is.leaf(node)) {
    if (is.null(attr(node, "is_collapsed"))) {
      sample_id <- labels(node)[1]
      
      attr(node, "label") <- sample_id
      attr(node, "species") <- unname(species_lookup[sample_id])
      attr(node, "n_members") <- 1L
      attr(node, "is_collapsed") <- FALSE
      attr(node, "height") <- 0
    }
    return(node)
  }
  for (i in seq_along(node)) {
    node[[i]] <- tag_singletons(node[[i]], species_lookup)
  }
  return(node)
}

dend_collapsed <- tag_singletons(dend_collapsed, species_lookup)

rebuild_tree_geometry <- function(node) {
  if (is.leaf(node)) {
    attr(node, "members") <- 1L
    attr(node, "midpoint") <- 0
    return(node)
  }
  
  m <- 0L
  for (i in seq_along(node)) {
    node[[i]] <- rebuild_tree_geometry(node[[i]])
    m <- m + attr(node[[i]], "members")
  }
  attr(node, "members") <- m
  
  if (length(node) == 2) {
    m1 <- attr(node[[1]], "midpoint")
    m2 <- attr(node[[2]], "midpoint")
    if (is.null(m1)) m1 <- 0
    if (is.null(m2)) m2 <- 0
    mem1 <- attr(node[[1]], "members")
    attr(node, "midpoint") <- (m1 + (mem1 + m2)) / 2
  } else {
    attr(node, "midpoint") <- 0
  }
  return(node)
}

dend_collapsed <- rebuild_tree_geometry(dend_collapsed)

extract_leaf_meta <- function(node) {
  meta <- list()
  collect <- function(n) {
    if (is.leaf(n)) {
      lab <- attr(n, "label")
      if (is.null(lab)) lab <- labels(n)[1]
      
      spec <- attr(n, "species")
      if (is.null(spec)) spec <- unname(species_lookup[lab])
      
      meta[[length(meta) + 1]] <<- data.frame(
        label = lab,
        species = ifelse(is.na(spec) || is.null(spec), NA_character_, as.character(spec)),
        n_members = ifelse(is.null(attr(n, "n_members")), 1L, attr(n, "n_members")),
        is_collapsed = ifelse(is.null(attr(n, "is_collapsed")), FALSE, attr(n, "is_collapsed")),
        stringsAsFactors = FALSE
      )
    } else {
      for (i in seq_along(n)) collect(n[[i]])
    }
  }
  collect(node)
  do.call(rbind, meta)
}

leaf_meta <- extract_leaf_meta(dend_collapsed)

hcdata <- dendro_data(dend_collapsed, type = "rectangle")

leaves <- hcdata$labels
leaves$species <- leaf_meta$species
leaves$n_members <- leaf_meta$n_members
leaves$is_collapsed <- leaf_meta$is_collapsed
leaves$display_label <- sprintf("(n = %d)", leaves$n_members)

dend_plot <- ggplot() +
  geom_segment(data = segment(hcdata),
               aes(x = x, y = y, xend = xend, yend = yend),
               linewidth = 0.4, color = "grey30") +
  geom_point(data = leaves,
             aes(x = x, y = 0, color = species, size = n_members),
             alpha = 0.85) +
  geom_text(data = leaves,
            aes(x = x, y = -0.012, label = display_label),
            hjust = 0, size = 3.0, color = "grey20") +
  scale_color_manual(values = species_colors, guide = "none") +
  scale_size_continuous(range = c(2.0, 5.5), guide = "none") +
  scale_y_reverse(expand = expansion(mult = c(0.15, 0.25))) +
  labs(title = "Hierarchical Clustering (1 - IBS)", y = "Distance (1 - IBS)", x = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5)
  ) +
  coord_flip()

legend_data <- data.frame(species = names(species_colors))
legend_plot <- ggplot(legend_data, aes(x = 1, y = species, color = species)) +
  geom_point(size = 3.5) +
  geom_text(aes(label = abbrev_species(species)), hjust = 0, nudge_x = 0.08, size = 3.5) +
  scale_color_manual(values = species_colors) +
  xlim(0.9, 2.5) +
  theme_void() +
  theme(legend.position = "none")

final_plot <- dend_plot / legend_plot + patchwork::plot_layout(heights = c(10, 1))
ggsave("ibs_hclust_dendrogram.png", plot = final_plot, width = 8, height = 12, dpi = 600, bg = "white")

# ADMIXTURE
extract_cv_error <- function(log_path) {
  lines <- readLines(log_path, warn = FALSE)
  cv_line <- lines[grepl("CV error", lines)]
  if (length(cv_line) == 0) {
    warning(paste("No CV error line found in", log_path))
    return(NA_real_)
  }
  as.numeric(str_extract(cv_line, "[0-9.]+$"))
}

cv_best <- tibble(K = 1:7) %>%
  mutate(
    log_file = sprintf("snowgums_clean_pruned.K%d.log", K),
    q_file   = sprintf("snowgums_clean_pruned.%d.Q", K)
  )

stopifnot(all(file.exists(cv_best$log_file)))

cv_best <- cv_best %>%
  mutate(cv_error = purrr::map_dbl(log_file, extract_cv_error))

print(cv_best)

elbow_plot <- ggplot(cv_best, aes(x = K, y = cv_error)) +
  geom_line(linewidth = 0.8, color = "grey40") +
  geom_point(size = 3, color = "black") +
  scale_x_continuous(breaks = cv_best$K) +
  labs(
    x = "K",
    y = "Cross-validation error",
  ) +
  theme_bw()

print(elbow_plot)
ggsave("admixture_cv_elbow.png", plot = elbow_plot, width = 6, height = 4.5, dpi = 600)


# Ancestry bar plots
fam <- read_table(
  "snowgums_clean_pruned.fam",
  col_names = c("fid", "sample_id", "father", "mother", "sex", "phenotype"),
  col_types = cols(.default = "c")
)

fam <- fam %>%
  mutate(
    sample_id_doubled = paste(sample_id, sample_id, sep = "_"),
    species_clean = species_lookup[sample_id_doubled]
  )

taxon_order <- c("Eucalyptus stellulata",
                 "Eucalyptus pauciflora subsp. pauciflora",
                 "Eucalyptus pauciflora subsp. niphophila")

sample_order <- fam %>%
  filter(!is.na(species_clean)) %>%
  mutate(species_clean = factor(species_clean, levels = taxon_order)) %>%
  arrange(species_clean, sample_id) %>%
  pull(sample_id)

read_q_long_fixed <- function(K, q_file, fam_df, sample_order) {
  q_mat <- read_table(q_file, col_names = sprintf("cluster_%d", 1:K), col_types = cols(.default = "d"))
  q_df <- bind_cols(fam_df %>% select(sample_id, species_clean), q_mat)
  
  q_df %>%
    filter(!is.na(species_clean)) %>%
    mutate(sample_id = factor(sample_id, levels = sample_order)) %>%
    pivot_longer(starts_with("cluster_"), names_to = "cluster", values_to = "ancestry") %>%
    mutate(
      K = K,
      cluster = paste0("K", K, "_", cluster)
    )
}

cv_subset <- cv_best %>% filter(K %in% c(2, 3))

admixture_long <- purrr::pmap_dfr(
  list(cv_subset$K, cv_subset$q_file),
  function(K, q_file) read_q_long_fixed(K, q_file, fam, sample_order)
)

map_cluster_colors_fixed <- function(df) {
  mean_df <- df %>%
    group_by(K, species_clean, cluster) %>%
    summarise(mean_anc = mean(ancestry), .groups = "drop")
  
  res_list <- list()
  
  sub_k2 <- filter(mean_df, K == 2)
  stell_c2 <- sub_k2 %>% filter(species_clean == "Eucalyptus stellulata") %>% slice_max(mean_anc, n = 1) %>% pull(cluster)
  pauci_c2 <- setdiff(unique(sub_k2$cluster), stell_c2)[1]
  
  res_list[["2"]] <- data.frame(
    cluster = c(stell_c2, pauci_c2),
    color = c(species_colors["Eucalyptus stellulata"], species_colors["Eucalyptus pauciflora subsp. pauciflora"])
  )
  
  sub_k3 <- filter(mean_df, K == 3)
  niph_c3  <- sub_k3 %>% filter(species_clean == "Eucalyptus pauciflora subsp. niphophila") %>% slice_max(mean_anc, n = 1) %>% pull(cluster)
  pauci_c3 <- sub_k3 %>% filter(species_clean == "Eucalyptus pauciflora subsp. pauciflora") %>% slice_max(mean_anc, n = 1) %>% pull(cluster)
  stell_c3 <- sub_k3 %>% filter(species_clean == "Eucalyptus stellulata") %>% slice_max(mean_anc, n = 1) %>% pull(cluster)
  
  res_list[["3"]] <- data.frame(
    cluster = c(niph_c3, pauci_c3, stell_c3),
    color = c(
      species_colors["Eucalyptus pauciflora subsp. niphophila"],
      species_colors["Eucalyptus pauciflora subsp. pauciflora"],
      species_colors["Eucalyptus stellulata"]
    )
  )
  
  bind_rows(res_list)
}

cluster_color_map <- map_cluster_colors_fixed(admixture_long)
cluster_color_vector <- setNames(cluster_color_map$color, cluster_color_map$cluster)

ancestry_plot <- ggplot(admixture_long, aes(x = sample_id, y = ancestry, fill = cluster)) +
  geom_col(width = 1, color = NA) +
  scale_fill_manual(values = cluster_color_vector, guide = "none") +
  facet_grid(K ~ species_clean, scales = "free_x", space = "free_x",
             labeller = labeller(K = function(k) paste0("K=", k),
                                 species_clean = abbrev_species)) +
  labs(
    x = NULL, y = "Ancestry proportion",
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.spacing.x = unit(0.15, "lines"),
    panel.grid = element_blank(),
    strip.text.y = element_text(face = "bold"),
    strip.text.x = element_text(face = "italic", size = 8)
  )

print(ancestry_plot)
ggsave("admixture_ancestry_k2_k3.png", plot = ancestry_plot, width = 10, height = 6, dpi = 600)

# Stats
taxon_files <- c(
  "Eucalyptus pauciflora subsp. niphophila" = "het_niphophila.het",
  "Eucalyptus pauciflora subsp. pauciflora" = "het_pauciflora.het",
  "Eucalyptus stellulata"                  = "het_stellulata.het"
)

het_raw <- purrr::map_dfr(names(taxon_files), function(tax_name) {
  file_path <- taxon_files[[tax_name]]
  
  read_tsv(file_path, show_col_types = FALSE) %>%
    clean_names() %>%
    mutate(
      iid           = paste(iid, iid, sep = "_"),
      obs_het       = (obs_ct - o_hom) / obs_ct,
      exp_het       = (obs_ct - e_hom) / obs_ct,
      species_clean = tax_name
    )
})

pop_stats <- het_raw %>%
  filter(!is.na(species_clean)) %>%
  group_by(species_clean) %>%
  summarize(
    N_Samples   = n(),
    Mean_Ho     = mean(obs_het, na.rm = TRUE),
    SD_Ho       = sd(obs_het, na.rm = TRUE),
    Mean_He     = mean(exp_het, na.rm = TRUE),
    Ho_He_ratio = mean(obs_het, na.rm = TRUE) / mean(exp_het, na.rm = TRUE),
    Mean_F      = 1 - Ho_He_ratio,
    SD_F        = sd(f, na.rm = TRUE),
    .groups     = "drop"
  )

print(pop_stats)
write_csv(pop_stats, "ch3_descriptive_stats.csv")

# Global Fst
short_to_full_species <- c(
  "niphophila"  = "Eucalyptus pauciflora subsp. niphophila",
  "pauciflora"  = "Eucalyptus pauciflora subsp. pauciflora",
  "stellulata"  = "Eucalyptus stellulata"
)

fst_pairwise <- read_tsv("pairwise_fst_results.fst.summary", show_col_types = FALSE,
                         col_names = c("pop1", "pop2", "wc_fst"), skip = 1) %>%
  mutate(
    pop1 = unname(short_to_full_species[pop1]),
    pop2 = unname(short_to_full_species[pop2])
  )

pop_names_filtered <- unname(short_to_full_species)

matrix_fst_data <- fst_pairwise %>%
  transmute(Pop1 = pop1, Pop2 = pop2, FST = pmax(wc_fst, 0))

full_matrix_data <- bind_rows(
  matrix_fst_data,
  matrix_fst_data %>% rename(Pop1 = Pop2, Pop2 = Pop1),
  tibble(Pop1 = pop_names_filtered, Pop2 = pop_names_filtered, FST = 0)
) %>%
  mutate(across(c(Pop1, Pop2), \(x) factor(x, levels = pop_names_filtered)))

fst_heatmap_plot <- ggplot(full_matrix_data, aes(x = Pop2, y = Pop1, fill = FST)) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(aes(label = sprintf("%.4f", FST), color = FST > 0.05), fontface = "bold", size = 4.2) +
  scale_fill_viridis_c(option = "rocket", direction = -1, name = expression(italic(F)[ST])) +
  scale_color_manual(values = c("TRUE" = "white", "FALSE" = "black"), guide = "none") +
  scale_x_discrete(labels = abbrev_species) +
  scale_y_discrete(labels = abbrev_species) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 12, base_family = "Times New Roman") +
  theme(
    text = element_text(family = "Times New Roman"),
    axis.text.x = element_text(angle = 35, hjust = 1, face = "bold", color = "black"),
    axis.text.y = element_text(face = "bold", color = "black"),
    panel.grid = element_blank()
  )

print(fst_heatmap_plot)
ggsave("clean_fst_matrix_heatmap.png", plot = fst_heatmap_plot, width = 7.5, height = 6.5, dpi = 600)

# IBD (whole-genome PLINK 1-IBS distance matrix; already computed)
mdist_ids <- read_tsv("pauci_stellulata_1minusIBS.mdist.id", show_col_types = FALSE,
                      col_names = c("fid", "iid")) %>%
  mutate(sample_id = paste(iid, iid, sep = "_"))

gen_mat <- as.matrix(read_tsv("pauci_stellulata_1minusIBS.mdist", show_col_types = FALSE,
                              col_names = FALSE))
rownames(gen_mat) <- colnames(gen_mat) <- mdist_ids$sample_id

valid_ibd_samples <- intersect(pca_data$sample_id, mdist_ids$sample_id)

geo_coords <- pca_data %>% 
  filter(sample_id %in% valid_ibd_samples, !is.na(latitude), !is.na(longitude))

geo_mat <- distm(geo_coords[, c("longitude", "latitude")], fun = distHaversine) / 1000
rownames(geo_mat) <- colnames(geo_mat) <- geo_coords$sample_id

valid_samples <- intersect(rownames(gen_mat), rownames(geo_mat))
species_map   <- setNames(pca_data$species_clean, pca_data$sample_id)
all_pairs     <- combn(valid_samples, 2, simplify = FALSE)

ibd_pair_data <- tibble(
  s1 = purrr::map_chr(all_pairs, 1),
  s2 = purrr::map_chr(all_pairs, 2)
) %>%
  mutate(
    sp1          = species_map[s1],
    sp2          = species_map[s2],
    geo_dist_km  = geo_mat[cbind(s1, s2)],
    genetic_dist = gen_mat[cbind(s1, s2)]
  )

shared_xlim <- c(0, 150)
shared_ylim <- range(ibd_pair_data$genetic_dist, na.rm = TRUE)

make_ibd_panel <- function(rowSp, colSp, show_x = FALSE, show_y = FALSE) {
  panel_data <- ibd_pair_data %>%
    filter(sp1 %in% c(rowSp, colSp) & sp2 %in% c(rowSp, colSp)) %>%
    mutate(
      pair_type = case_when(
        sp1 == rowSp & sp2 == rowSp ~ "within_row",
        sp1 == colSp & sp2 == colSp ~ "within_col",
        TRUE ~ "between"
      )
    )
  
  ggplot(panel_data, aes(x = geo_dist_km, y = genetic_dist, colour = pair_type)) +
    geom_point(size = 0.6, alpha = 0.3) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.9) +
    scale_colour_manual(
      values = c("within_row" = "red", "within_col" = "blue", "between" = "purple"), 
      guide = "none"
    ) +
    coord_cartesian(xlim = shared_xlim, ylim = shared_ylim) +
    labs(
      x = if (show_x) "Geographic distance (km)" else NULL,
      y = if (show_y) "Genetic distance (1 - IBS)" else NULL
    ) +
    theme_bw(base_size = 8) +
    theme(
      axis.text.x = if (show_x) element_text(size = 6) else element_blank(),
      axis.text.y = if (show_y) element_text(size = 6) else element_blank(),
      axis.ticks  = element_blank(),
      plot.margin = margin(2, 2, 2, 2)
    )
}

make_label <- function(label, angle = 0, bold = TRUE) {
  ggplot() + 
    theme_void() +
    annotate("text", x = 0, y = 0, label = abbrev_species(label), 
             angle = angle, fontface = ifelse(bold, "bold", "plain"), size = 3) +
    theme(plot.margin = margin(1, 1, 1, 1))
}

pops <- pop_names_filtered
n <- length(pops)
grid_cells <- list()

for (i in 2:n) {
  for (j in 1:(n-1)) {
    if (j == 1) {
      grid_cells <- append(grid_cells, list(make_label(pops[i], angle = 90)))
    }
    if (i > j) {
      grid_cells <- append(grid_cells, list(
        make_ibd_panel(pops[i], pops[j], show_x = (i == n), show_y = (j == 1))
      ))
    } else {
      grid_cells <- append(grid_cells, list(plot_spacer()))
    }
  }
}

grid_cells <- append(grid_cells, list(plot_spacer()))
for (j in 1:(n-1)) {
  grid_cells <- append(grid_cells, list(make_label(pops[j])))
}

n_data <- n - 1

ibd_matrix_plot <- wrap_plots(grid_cells, ncol = n_data + 1, nrow = n_data + 1) +
  plot_layout(
    widths  = c(0.12, rep(1, n_data)),
    heights = c(rep(1, n_data), 0.1)
  ) +
  plot_annotation(
    theme = theme(plot.title = element_text(face = "bold"))
  )

print(ibd_matrix_plot)
ggsave("ibd_scatterplot_matrix.png", plot = ibd_matrix_plot,
       width = 3.2 * n_data + 1, height = 3.2 * n_data + 0.5, dpi = 300)

taxa <- unique(c(ibd_pair_data$sp1, ibd_pair_data$sp2))

mantel_results <- list()

for (i in 1:length(taxa)) {
  for (j in i:length(taxa)) {
    t1 <- taxa[i]
    t2 <- taxa[j]
    samples_t1 <- names(species_map[species_map == t1])
    samples_t2 <- names(species_map[species_map == t2])
    
    if (t1 == t2) {
      sub_samples <- intersect(valid_samples, samples_t1)
      if (length(sub_samples) < 3) next
      
      gen_sub <- as.dist(gen_mat[sub_samples, sub_samples])
      geo_sub <- as.dist(geo_mat[sub_samples, sub_samples])
      
      m_test <- mantel(gen_sub, geo_sub, method = "pearson", permutations = 999)
      
      mantel_results[[length(mantel_results) + 1]] <- tibble(
        Comparison_Type = "Within Taxon",
        Taxon_1         = t1,
        Taxon_2         = t2,
        Mantel_r        = m_test$statistic,
        p_value         = m_test$signif
      )
    } else {
      sub_samples <- intersect(valid_samples, c(samples_t1, samples_t2))
      if (length(sub_samples) < 3) next
      
      gen_sub <- as.dist(gen_mat[sub_samples, sub_samples])
      geo_sub <- as.dist(geo_mat[sub_samples, sub_samples])
      
      m_test <- mantel(gen_sub, geo_sub, method = "pearson", permutations = 999)
      
      mantel_results[[length(mantel_results) + 1]] <- tibble(
        Comparison_Type = "Between Taxa",
        Taxon_1         = t1,
        Taxon_2         = t2,
        Mantel_r        = m_test$statistic,
        p_value         = m_test$signif
      )
    }
  }
}

mantel_df <- bind_rows(mantel_results)

elevation_summary <- ibd_pair_data %>%
  mutate(
    pair_label = case_when(
      sp1 == sp2 ~ paste0("Within ", sp1),
      TRUE ~ paste0("Between ", pmin(sp1, sp2), " & ", pmax(sp1, sp2))
    )
  ) %>%
  group_by(pair_label) %>%
  group_modify(~ {
    fit <- lm(genetic_dist ~ geo_dist_km, data = .x)
    tibble(
      Mean_Genetic_Dist = mean(.x$genetic_dist, na.rm = TRUE),
      Intercept_at_0km  = coef(fit)[1],
      Slope_per_km      = coef(fit)[2],
      R_squared         = summary(fit)$r.squared
    )
  })

print(mantel_df)
print(elevation_summary)
write_csv(mantel_df, "ibd_mantel_test_results.csv")
write_csv(elevation_summary, "ibd_elevation_regression_summary.csv")

# Windowed Fst
process_fst <- function(file_path, comparison_label, min_variants = 10) {
  read_tsv(file_path, show_col_types = FALSE) %>%
    filter(N_VARIANTS >= min_variants) %>%
    mutate(
      WEIGHTED_FST = pmax(WEIGHTED_FST, 0),
      bin_mid_mb   = (BIN_START + BIN_END) / 2 / 1e6,
      CHROM        = factor(CHROM),
      Comparison   = comparison_label
    )
}

plot_fst_facets <- function(data, title_str, point_color = "#2b5c8f", smooth_color = "#e41a1c") {
  num_comparisons <- length(unique(data$Comparison))
  facet_formula   <- if (num_comparisons > 1) Comparison ~ CHROM else ~ CHROM
  
  ggplot(data, aes(x = bin_mid_mb, y = WEIGHTED_FST)) +
    geom_point(alpha = 0.5, size = 0.8, color = point_color) +
    geom_smooth(method = "loess", span = 0.05, color = smooth_color, se = FALSE, linewidth = 0.7) +
    facet_grid(facet_formula, scales = "free_x") +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
    labs(
      x     = "Genomic Position (Mb)",
      y     = expression(Weighted~italic(F)[ST])
    ) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey92"),
      strip.text       = element_text(face = "bold", size = 8),
      strip.text.y     = element_text(angle = 0, face = "bold.italic") # keeps stacked labels horizontal
    )
}

fst_niph_pauc <- process_fst("fst_niph_vs_pauc_50kb.windowed.weir.fst", "niphophila vs pauciflora")
fst_niph_stel <- process_fst("fst_niph_vs_stel_50kb.windowed.weir.fst", "niphophila vs stellulata")
fst_pauc_stel <- process_fst("fst_pauc_vs_stel_50kb.windowed.weir.fst", "pauciflora vs stellulata")

plot_subsp <- plot_fst_facets(
  data      = fst_niph_pauc,
)
print(plot_subsp)
ggsave("fst_scan_niph_vs_pauc_50kb.png", plot = plot_subsp, width = 16, height = 3.5, dpi = 300)

fst_stell_combined <- bind_rows(fst_niph_stel, fst_pauc_stel)

plot_stell <- plot_fst_facets(
  data      = fst_stell_combined,
)
print(plot_stell)
ggsave("fst_scan_stellulata_stacked_50kb.png", plot = plot_stell, width = 16, height = 6.5, dpi = 300)

# Chr6 elevated region summary stats
target_chr_id <- "CM074670.1"
region_start  <- 38e6
region_end    <- 48.5e6

background_mean <- mean(fst_niph_pauc$WEIGHTED_FST, na.rm = TRUE)
n_windows_total <- nrow(fst_niph_pauc)

chr6_region <- fst_niph_pauc %>%
  filter(CHROM == target_chr_id, BIN_START >= region_start, BIN_START <= region_end)

region_mean <- mean(chr6_region$WEIGHTED_FST, na.rm = TRUE)
region_max  <- max(chr6_region$WEIGHTED_FST, na.rm = TRUE)

mean_excl_region <- fst_niph_pauc %>%
  filter(!(CHROM == target_chr_id & BIN_START >= region_start & BIN_START <= region_end)) %>%
  pull(WEIGHTED_FST) %>%
  mean(na.rm = TRUE)

pct_drop_region <- (background_mean - mean_excl_region) / background_mean * 100

# Compare against removing the single highest-peak window (+/-50kb buffer) on every other chromosome
other_chrom_peak_effects <- fst_niph_pauc %>%
  filter(CHROM != target_chr_id) %>%
  group_by(CHROM) %>%
  group_modify(~ {
    peak_row   <- .x %>% slice_max(WEIGHTED_FST, n = 1, with_ties = FALSE)
    peak_start <- peak_row$BIN_START - 50000
    peak_end   <- peak_row$BIN_END + 50000
    mean_excl_peak <- fst_niph_pauc %>%
      filter(!(CHROM == .y$CHROM & BIN_START >= peak_start & BIN_START <= peak_end)) %>%
      pull(WEIGHTED_FST) %>%
      mean(na.rm = TRUE)
    tibble(
      peak_max        = peak_row$WEIGHTED_FST,
      peak_max_ratio  = peak_row$WEIGHTED_FST / background_mean,
      pct_drop_removal = (background_mean - mean_excl_peak) / background_mean * 100
    )
  })

cat("Genome-wide background mean Fst:", round(background_mean, 4), "\n")
cat("Chr6 region mean Fst:", round(region_mean, 4), "  (", round(region_mean / background_mean, 2), "x background )\n")
cat("Chr6 region max Fst:", round(region_max, 4), "  (", round(region_max / background_mean, 2), "x background )\n")
cat("Genome-wide mean Fst excluding chr6 region:", round(mean_excl_region, 4), "\n")
cat("Percent drop from removing chr6 region:", round(pct_drop_region, 2), "%\n")
print(other_chrom_peak_effects)

write_csv(other_chrom_peak_effects, "fst_other_chrom_peak_effects.csv")

# Local PCA (Chr6 elevated region; PLINK2 --pca and --het, unpruned, already computed on Gadi)
chr6_pca <- read_tsv("chr6_region_pca.eigenvec", show_col_types = FALSE,
                     col_names = c("fid", "iid", paste0("PC", 1:10)), skip = 1) %>%
  mutate(sample_id = paste(iid, iid, sep = "_"))

chr6_het <- read_tsv("chr6_region_het.het", show_col_types = FALSE,
                     col_names = c("fid", "iid", "o_hom", "e_hom", "obs_ct", "f"), skip = 1) %>%
  mutate(
    local_het = (obs_ct - o_hom) / obs_ct,
    sample_id = paste(iid, iid, sep = "_")
  )

chr6_local <- chr6_pca %>%
  left_join(chr6_het %>% select(sample_id, local_het), by = "sample_id") %>%
  left_join(meta %>% select(sample_id = fastq_sample_name, species_clean), by = "sample_id")

chr6_local_plot <- ggplot(chr6_local, aes(x = PC1, y = local_het, colour = species_clean)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_colour_manual(values = species_colors) +
  labs(
    x = "PC1 (chr6 38-48.5Mb region)",
    y = "Local Heterozygosity",
    colour = "Taxa"
  ) +
  theme_bw()

print(chr6_local_plot)
ggsave("chr6_local_pca_het.png", plot = chr6_local_plot, width = 8, height = 6, dpi = 600)

# Latent factors
env_data <- read_csv("current_covariates.csv") %>%
  mutate(fastq_sample_name = paste(sample_id, sample_id, sep = "_"))

env_pca_input <- pca_data %>%
  filter(species_clean == "Eucalyptus pauciflora subsp. pauciflora") %>%
  select(sample_id) %>%
  left_join(env_data, by = c("sample_id" = "fastq_sample_name")) %>%
  column_to_rownames("sample_id") %>%
  select(lat, long, elev, starts_with("bio")) %>%
  drop_na()

env_pca <- prcomp(env_pca_input, center = TRUE, scale = TRUE)
env_pc1 <- env_pca$x[, 1]

vp <- data.frame(
  varprop = pca_results$pca$varprop, 
  pc = sprintf("PC%02d", 1:length(pca_results$pca$varprop))
) %>% filter(is.finite(varprop))

pear_pca_cor <- pca_data %>%
  filter(species_clean == "Eucalyptus pauciflora subsp. pauciflora") %>%
  pivot_longer(starts_with("PC"), names_to = "genetic_pc", values_to = "pc_score") %>%
  left_join(as.data.frame(env_pc1) %>% rownames_to_column("sample_id") %>% rename(env_pc1 = 2), by = "sample_id") %>%
  group_by(genetic_pc) %>% 
  summarize(
    correlation_to_environment = cor(pc_score, env_pc1, method = "pearson", use = "pairwise.complete.obs"),
    .groups = "drop"
  ) %>%
  left_join(vp, by = join_by(genetic_pc == pc)) %>%
  mutate(scaled_cor = correlation_to_environment * varprop)

pear_pca_cor %>%
  pivot_longer(-genetic_pc) %>%
  mutate(name = factor(name,
                       levels = c("varprop", "correlation_to_environment", "scaled_cor"),
                       labels = c("Variance Explained", "Correlation to Env PC1", "Scaled Correlation")
  )) %>%
  ggplot(aes(x = genetic_pc, y = value, group = name, color = name)) + 
  geom_line() + 
  geom_point() +
  facet_grid(name ~ ., scales = "free_y") +
  labs(
    title = "Environmental PC1 Correlation & Scaling (E. pauciflora subsp. pauciflora)",
    x = "Genetic PC",
    y = NULL
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5),
    legend.position = "none"
  )

ggsave("env_pc1_correlation_scaling.png", width = 10, height = 8, dpi = 300)
