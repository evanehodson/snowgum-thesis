library(tidyverse)
library(duckplyr)
library(DBI)
library(duckdb)

# 1. Connect & Unpivot Parquet Data ---------------------------------------
con <- DBI::dbConnect(duckdb::duckdb())

DBI::dbExecute(con, "
    CREATE OR REPLACE TEMP VIEW pvals AS
    FROM read_parquet('lfmm_pauci_k9_nonmf.parquet')
    UNPIVOT (p FOR envvar IN (COLUMNS('^p_.*')))
    SELECT chr, pos, envvar, p
")

# Extract distinct environmental variable names
envvars <- DBI::dbGetQuery(con, "SELECT DISTINCT envvar FROM pvals ORDER BY envvar")[[1]]

# 2. Generate Plots per Variable ------------------------------------------
for (ev in envvars) {
  message("Processing: ", ev)
  
  # Clean variable name for file naming
  clean_ev <- gsub("^p_", "", ev)
  
  # --- Manhattan Plot ---
  # Filter points with p < 1e-2 to keep rendering fast and file sizes manageable
  man_df <- tbl(con, "pvals") |>
    filter(envvar == !!ev, p < 1e-2) |>
    collect()
  
  p_man <- ggplot(man_df, aes(x = pos, y = -log10(p))) +
    geom_point(size = 0.6, alpha = 0.7, color = "darkblue") +
    facet_grid(~ chr, scales = "free_x", space = "free_x") +
    labs(
      title = paste("Manhattan Plot -", clean_ev),
      x = "Chromosome Position",
      y = expression(-log[10](p))
    ) +
    theme_bw() +
    theme(
      panel.spacing = unit(0.1, "lines"),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )
  
  ggsave(
    filename = paste0("manhattan_", clean_ev, ".png"),
    plot = p_man,
    width = 12,
    height = 5,
    dpi = 300
  )
  
  # --- QQ Plot ---
  # Binned p-values aggregated in DuckDB for quick QQ calculation
  qq_hist <- tbl(con, sql(sprintf("
      FROM pvals
      SELECT ROUND(-LN(p)/LN(10.0), 3) AS logp_bin,
             COUNT(*) AS n
      WHERE envvar = '%s' AND p > 0
      GROUP BY ROUND(-LN(p)/LN(10.0), 3)
      ORDER BY logp_bin DESC
  ", ev))) |> 
    collect()
  
  qq_data <- qq_hist |>
    arrange(desc(logp_bin)) |>
    mutate(
      cum_n = cumsum(n),
      total_n = sum(n),
      expected_logp = -log10(cum_n / total_n)
    )
  
  p_qq <- ggplot(qq_data, aes(x = expected_logp, y = logp_bin)) +
    geom_line(color = "blue", linewidth = 0.8) +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
    labs(
      title = paste("QQ Plot -", clean_ev),
      x = expression(Expected ~ -log[10](p)),
      y = expression(Observed ~ -log[10](p))
    ) +
    theme_bw()
  
  ggsave(
    filename = paste0("qq_", clean_ev, ".png"),
    plot = p_qq,
    width = 6,
    height = 6,
    dpi = 300
  )
}

DBI::dbDisconnect(con, shutdown = TRUE)