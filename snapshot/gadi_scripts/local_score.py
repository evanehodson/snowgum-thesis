import duckdb

con = duckdb.connect()

# 1. Load data into DuckDB view
con.execute("""
    CREATE OR REPLACE TEMP VIEW pvals AS
    FROM read_parquet('/scratch/xe2/eh4663/LFMM/lfmm_pauci_k9_nonmf.parquet')
    UNPIVOT (p FOR envvar IN (COLUMNS('^p_.*')))
    SELECT chr, pos, envvar, p
""")

A_COEF = {
    1: (-5.50, 6.76, -5.66, -2.51),
    2: (2.47, -4.16, -1.82, -4.58),
    3: (2.04, -5.76, 1.04, -6.95),
    4: (0.22, -4.08, 1.16, -9.16),
}
B_COEF = {
    1: (-1.22, 3.17, -1.99),
    2: (0.37, 2.14, -2.35),
    3: (2.55, -0.02, -2.31),
    4: (3.45, -0.98, -2.33),
}


def local_score_duckdb(source, con, xi=2, alpha=0.05):
  a, b = A_COEF[xi], B_COEF[xi]

  con.execute(f"""
        CREATE OR REPLACE TEMP TABLE _ls_thresh AS
        WITH chr_stats AS (
            SELECT envvar, chr,
                   COUNT(*) AS nmarker,
                   COALESCE(ABS(CORR(p, lag_p)), 0.0) AS rho
            FROM (
                SELECT envvar, chr, p,
                       LAG(p) OVER (PARTITION BY envvar, chr ORDER BY pos) AS lag_p
                FROM {source}
            )
            GROUP BY envvar, chr
        )
        SELECT envvar, chr, nmarker, rho,
               (LN(-LN(1.0 - {alpha})) - (
                   LN(nmarker) + {a[0]}*POW(rho,3) + {a[1]}*POW(rho,2) + {a[2]}*rho + {a[3]}
               )) / (
                   {b[0]}*POW(rho,2) + {b[1]}*rho + {b[2]}
               ) AS lindley_thresh
        FROM chr_stats
    """)

  con.execute("""
        CREATE OR REPLACE TEMP TABLE _ls_result (
            envvar VARCHAR, chr VARCHAR, pos BIGINT, p DOUBLE,
            score DOUBLE, lindley_score DOUBLE, lindley_thresh DOUBLE,
            lindley_window INTEGER
        )
    """)

  envvars = [
      r[0]
      for r in con.execute(
          "SELECT DISTINCT envvar FROM _ls_thresh ORDER BY envvar"
      ).fetchall()
  ]

  for ev in envvars:
    con.execute(f"""
            INSERT INTO _ls_result
            WITH scores AS (
                SELECT p.envvar, p.chr, p.pos, p.p,
                       LN(1.0/p.p)/LN(10.0) - {xi} AS score,
                       t.lindley_thresh
                FROM {source} p
                JOIN _ls_thresh t ON p.envvar = t.envvar AND p.chr = t.chr
                WHERE p.envvar = '{ev}'
            ),
            cumul AS (
                SELECT *, SUM(score) OVER (PARTITION BY chr ORDER BY pos) AS cumsum_score
                FROM scores
            ),
            lindley_calc AS (
                SELECT *,
                       cumsum_score - LEAST(0.0,
                           MIN(cumsum_score) OVER (PARTITION BY chr ORDER BY pos ROWS UNBOUNDED PRECEDING)
                       ) AS lindley_score
                FROM cumul
            ),
            transitions AS (
                SELECT *,
                       CASE WHEN lindley_score > lindley_thresh
                            AND COALESCE(LAG(lindley_score) OVER (PARTITION BY chr ORDER BY pos) <= lindley_thresh, TRUE)
                            THEN 1 ELSE 0 END AS is_new_window
                FROM lindley_calc
            )
            SELECT envvar, chr, pos, p, score, lindley_score, lindley_thresh,
                   CASE WHEN lindley_score > lindley_thresh
                        THEN SUM(is_new_window) OVER (PARTITION BY chr ORDER BY pos)
                        ELSE NULL END AS lindley_window
            FROM transitions
            WHERE lindley_score > 1
            ORDER BY chr, pos
        """)

  return con.table("_ls_result")


# 1. Run local score calculation using standard alpha=0.05
local_score_duckdb("pvals", con, xi=2, alpha=0.05)

# 2. Export filtered and merged candidate regions with peak score threshold >= 100
con.execute("""
    COPY (
        WITH raw_intervals AS (
            SELECT chr, MIN(pos) AS start, MAX(pos) AS stop,
                   MAX(lindley_score) AS peak_score
            FROM _ls_result
            WHERE lindley_window IS NOT NULL
            GROUP BY envvar, chr, lindley_window
        ),
        filtered_intervals AS (
            -- Filter to keep only high-confidence peaks (peak score >= 100)
            SELECT chr, start, stop
            FROM raw_intervals
            WHERE peak_score >= 100
        ),
        ordered_intervals AS (
            SELECT chr, start, stop,
                   ROW_NUMBER() OVER (PARTITION BY chr ORDER BY start, stop) as rn
            FROM filtered_intervals
        ),
        merged_groups AS (
            SELECT chr, start, stop,
                   SUM(CASE WHEN start > max_prev_stop THEN 1 ELSE 0 END) OVER (PARTITION BY chr ORDER BY start, stop) as grp
            FROM (
                SELECT chr, start, stop,
                       COALESCE(MAX(stop) OVER (PARTITION BY chr ORDER BY start, stop ROWS BETWEEN UNBOUNDED PRECEDING AND 1 PRECEDING), start) as max_prev_stop
                FROM ordered_intervals
            )
        )
        SELECT chr, MIN(start) AS start, MAX(stop) AS stop
        FROM merged_groups
        GROUP BY chr, grp
        ORDER BY chr, start
    ) TO 'candidates_merged.bed' (HEADER false, DELIMITER '\t')
""")
print("Successfully generated 'candidates_merged.bed'.")