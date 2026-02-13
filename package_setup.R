# Prevent scientific notation globally
options(scipen = 999)

# Install vctrs first (do not load any other packages yet)
if (!requireNamespace("vctrs", quietly = TRUE)) {
  install.packages("vctrs")
}

# Install required packages for setup + city inspection
required_packages <- c("data.table", "dbscan", "ggplot2", "dplyr", "mco")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  install.packages(missing_packages)
}

# Load packages AFTER installation
suppressPackageStartupMessages({
  library(vctrs)
  library(data.table)
  library(dbscan)
  library(ggplot2)
  library(dplyr)
  library(mco)
})

############################################################
# Step 1 — Load + inspect all cities (no combining)
#
# WHY:
# - Quick schema check across cities in one run
# - Confirms consistency before you start clustering
# - Does NOT merge/stack datasets (each city stays separate)
############################################################

city_paths <- list(
  berlin  = "/content/drive/MyDrive/berlin_flows_utm.rds",
  munich  = "/content/drive/MyDrive/munich_flows_utm.rds",
  cologne = "/content/drive/MyDrive/cologne_flows_utm.rds"
)

# Store each city separately (no combining)
city_data <- lapply(names(city_paths), function(city_id) {
  dt <- readRDS(city_paths[[city_id]])
  dt <- as.data.table(dt)
  dt[, city := city_id]  # traceability tag
  dt
})
names(city_data) <- names(city_paths)

# Print inspection for each city
for (city_id in names(city_data)) {
  dt <- city_data[[city_id]]

  cat("\n=============================\n")
  cat("CITY:", toupper(city_id), "\n")
  cat("=============================\n")

  cat("class():\n")
  print(class(dt))

  cat("\ndim():\n")
  print(dim(dt))

  cat("\nnames():\n")
  print(names(dt))
}

############################################################
# Step 2 — Structural correction + scaling (per city)
#
# WHAT this step produces:
# - X_scaled: scaled feature matrix for clustering
# - city_clean: filtered dataset aligned with X_scaled rows
#
# Key decisions (explicit WHY):
# - We drop rows with NA in clustering features (complete.cases)
#   because distance-based clustering cannot handle missing values reliably.
# - We scale features because DBSCAN/HDBSCAN use distance computations;
#   scaling prevents any one dimension from dominating due to units/magnitude.
############################################################

make_X_and_clean <- function(city_dt,
                             feature_cols = c("x_o", "y_o", "dx", "dy")) {

  # WHY: fail fast if a city has a different schema (examiner-friendly)
  missing_cols <- setdiff(feature_cols, names(city_dt))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  # Build feature matrix (same as your berlin snippet, just generalized)
  X <- city_dt %>%
    dplyr::select(dplyr::all_of(feature_cols)) %>%
    as.matrix()

  # Drop rows with any NA in features
  ok <- complete.cases(X)

  # Scale only the valid rows
  X_scaled <- scale(X[ok, , drop = FALSE])

  # Keep dataset aligned with X_scaled rows
  city_clean <- city_dt[ok, ]

  list(X = X[ok, , drop = FALSE], ok = ok, X_scaled = X_scaled, city_clean = city_clean)
}

# ---- Run for all cities (no combining)
# Assumes you already created: city_data$berlin, city_data$munich, city_data$cologne
prepared <- lapply(names(city_data), function(city_id) {
  res <- make_X_and_clean(city_data[[city_id]])
  res$city_id <- city_id
  res
})
names(prepared) <- names(city_data)

############################################################
# Step 2b — Store total trips per city (used by eval_solution)
############################################################

N_total_by_city <- sapply(names(prepared), function(city_id) {
  nrow(prepared[[city_id]]$X)   # filtered X, aligned with X_scaled
})

print(N_total_by_city)

# ---- Optional: keep Berlin-style objects for your existing code compatibility
berlin_clean  <- prepared$berlin$city_clean
munich_clean  <- prepared$munich$city_clean
cologne_clean <- prepared$cologne$city_clean

X_scaled_berlin  <- prepared$berlin$X_scaled
X_scaled_munich  <- prepared$munich$X_scaled
X_scaled_cologne <- prepared$cologne$X_scaled

# Quick audit print (helps you catch NA issues early)
for (city_id in names(prepared)) {
  cat("\n---", toupper(city_id), "---\n")
  cat("Rows original:", nrow(city_data[[city_id]]), "\n")
  cat("Rows kept:",     nrow(prepared[[city_id]]$city_clean), "\n")
  cat("Rows dropped:",  sum(!prepared[[city_id]]$ok), "\n")
  cat("X_scaled dim:",  paste(dim(prepared[[city_id]]$X_scaled), collapse = " x "), "\n")
}

############################################################
# Step 3 — kNN-distance plot configuration
#
# WHY this step:
# - kNN-distance plots are used to visually assess a reasonable
#   neighborhood scale (eps) for density-based clustering.
# - We inspect multiple k values to understand sensitivity
#   across different minimum neighborhood sizes.
#
# IMPORTANT:
# - These k values are diagnostic only (not optimization yet).
# - They are aligned with minPts values later used in DBSCAN/HDBSCAN
#   for consistency and interpretability.
############################################################

# Chosen k values (same across all cities)
# WHY these values:
# - 10   : captures local dense patterns
# - 20-30: moderate neighborhood size (often stable)
# - 50   : stress-test for over-smoothing
k_values <- c(10, 20, 30, 50)

############################################################
# Step 3a — Plot layout setup (kNN-distance diagnostics)
#
# WHY this step:
# - We visualize multiple kNN-distance plots side-by-side
#   to compare neighborhood scales for different k values.
# - A fixed 2x2 layout is sufficient for the chosen k grid (4 values).
#
# SAFETY:
# - We store the original graphical parameters and restore them later
#   to avoid side effects on subsequent plots.
############################################################

# Save current graphical parameters
old_par <- par(no.readonly = TRUE)

# Set plotting layout (one panel per k value)
par(mfrow = c(2, 2))

############################################################
# Step 3b — kNN-distance plots
# WHY: visual support for choosing a reasonable eps scale
############################################################

for (city_id in names(prepared)) {

  X_scaled <- prepared[[city_id]]$X_scaled

  cat("\nPlotting kNN-distance curves for:", toupper(city_id), "\n")

  old_par_city <- par(no.readonly = TRUE)
  par(mfrow = c(2, 2))

  for (k in k_values) {
    kNN <- kNNdist(X_scaled, k = k)

    plot(sort(kNN), type = "l",
         xlab = "Points sorted by distance",
         ylab = paste0(k, "-NN distance"),
         main = paste(city_id, "| kNN distance (k =", k, ")"))

    grid()
  }

  par(old_par_city)
}

# Restore original graphics settings
par(old_par)



############################################################
# Step 4 — DBSCAN parameter grids
# WHY: explore sensitivity of density scale (eps) and minPts
############################################################

# Neighborhood radius candidates
eps_grid <- c(0.2, 0.3, 0.4, 0.5, 0.6)

# minPts candidates (aligned with kNN diagnostics)
k_grid <- c(20, 30, 50)

############################################################
# Step 5 — Initialize results container
# WHY: collect clustering outputs across cities and parameters
############################################################

results <- list()

############################################################
# Step 5a — Result index counter
# WHY: sequentially store results in a list
############################################################

idx <- 1

############################################################
# Step 5b — Sanity checks / debugging prints
############################################################

# Index counter
if (exists("idx")) print(idx)

# Parameter grids
print(k_values)
print(eps_grid)
print(k_grid)

# Results container
if (exists("results")) print(results)

# Check feature matrices per city
for (city_id in names(prepared)) {
  cat("\n---", toupper(city_id), "---\n")
  cat("X head:\n")
  print(head(prepared[[city_id]]$X))
  cat("X_scaled head:\n")
  print(head(prepared[[city_id]]$X_scaled))
}

############################################################
# Step 6 — DBSCAN grid search
# WHY: evaluate clustering structure across (eps, minPts) values
############################################################

results <- list()
idx <- 1

for (city_id in names(prepared)) {

  X_scaled <- prepared[[city_id]]$X_scaled
  n_total  <- nrow(X_scaled)

  for (k in k_grid) {
    for (eps_val in eps_grid) {

      db <- dbscan(X_scaled, eps = eps_val, minPts = k)
      cl <- db$cluster

      n_noise   <- sum(cl == 0)
      noise_pct <- 100 * n_noise / length(cl)

      n_clusters <- length(unique(cl[cl != 0]))

      max_cluster_size <- if (n_clusters > 0) {
        max(table(cl[cl != 0]))
      } else {
        0
      }

      max_cluster_pct <- if (n_clusters > 0) {
        100 * max_cluster_size / length(cl)
      } else {
        0
      }

      results[[idx]] <- data.frame(
        city = city_id,
        algorithm = "DBSCAN",
        minPts = k,
        eps = eps_val,
        n_clusters = n_clusters,
        noise_pct = noise_pct,
        max_cluster_pct = max_cluster_pct,
        n_total = n_total
      )

      idx <- idx + 1
    }
  }
}

grid_df <- dplyr::bind_rows(results)
print(grid_df)

############################################################
# Step 7 — Final DBSCAN parameters (per city)
############################################################

dbscan_params <- list(
  berlin  = list(eps = 0.35, minPts = 20),  # DBCV = 0.062
  munich  = list(eps = 0.35, minPts = 20),  # DBCV = 0.087
  cologne = list(eps = 0.45, minPts = 20)   # updated due to DBCV = -0.329
)

############################################################
# Step 7a — Final DBSCAN run (all cities)
############################################################

for (city_id in names(prepared)) {

  params <- dbscan_params[[city_id]]

  db_final <- dbscan(
    prepared[[city_id]]$X_scaled,
    eps    = params$eps,
    minPts = params$minPts
  )

  stopifnot(
    length(db_final$cluster) == nrow(prepared[[city_id]]$city_clean)
  )

  prepared[[city_id]]$city_clean$cluster <- db_final$cluster

  cat(
    "DBSCAN done |",
    city_id,
    "| eps =", params$eps,
    "| minPts =", params$minPts,
    "\n"
  )
}

############################################################
# Step 8 — Verify final clustering (cluster size check)
############################################################

for (city_id in names(prepared)) {

  cat("\n---", toupper(city_id), "---\n")

  cl_tab <- table(prepared[[city_id]]$city_clean$cluster)

  # Print first clusters (including noise = 0)
  print(head(cl_tab, 65))
}


############################################################
# Step 8b — DBCV (internal validation, expert configuration)
# (works with dbcv() that returns a LIST with $score)
############################################################

set.seed(1)

dbcv_scores <- list()

N_MAX_TOTAL <- 20000
N_MAX_PERCLUSTER <- 2000
MIN_CLUSTER_FOR_DBCV <- 3

for (city_id in names(prepared)) {

  X_scaled <- prepared[[city_id]]$X_scaled
  cl_full <- prepared[[city_id]]$city_clean$cluster

  # Exclude noise (0)
  keep <- cl_full != 0
  Xc <- X_scaled[keep, , drop = FALSE]
  cl <- cl_full[keep]

  if (length(unique(cl)) < 2) {
    dbcv_scores[[city_id]] <- NA_real_
    cat("CITY:", toupper(city_id), "| DBCV = NA (fewer than 2 clusters)\n")
    next
  }

  # Stratified sampling
  dt <- data.table(idx = seq_along(cl), cluster = cl)
  samp_idx <- dt[, .(idx = sample(idx, size = min(.N, N_MAX_PERCLUSTER))), by = cluster]$idx
  if (length(samp_idx) > N_MAX_TOTAL) samp_idx <- sample(samp_idx, N_MAX_TOTAL)

  X_s <- Xc[samp_idx, , drop = FALSE]
  cl_s <- cl[samp_idx]

  # Drop tiny clusters (diagnostic stability)
  tab <- table(cl_s)
  keep_clusters <- names(tab)[tab >= MIN_CLUSTER_FOR_DBCV]
  keep2 <- cl_s %in% keep_clusters
  X_s <- X_s[keep2, , drop = FALSE]
  cl_s <- cl_s[keep2]

  # Relabel clusters to 1..K
  cl_s <- as.integer(factor(cl_s))

  if (length(unique(cl_s)) < 2) {
    dbcv_scores[[city_id]] <- NA_real_
    cat("CITY:", toupper(city_id), "| DBCV = NA (after filtering tiny clusters)\n")
    next
  }

  # Compute DBCV — IMPORTANT: pass d = ncol(X_s)
  dbcv_raw <- tryCatch(
    dbcv(X_s, cl_s, d = ncol(X_s), metric = "euclidean"),
    error = function(e) e
  )

  if (inherits(dbcv_raw, "error")) {
    dbcv_scores[[city_id]] <- NA_real_
    cat("CITY:", toupper(city_id), "| DBCV = NA (dbcv() error)\n")
    next
  }

  # Extract numeric score (your version returns a list)
  dbcv_val <- dbcv_raw$score
  if (!is.numeric(dbcv_val) || length(dbcv_val) != 1 || !is.finite(dbcv_val)) {
    dbcv_val <- NA_real_
  }

  dbcv_scores[[city_id]] <- dbcv_val

  cat(
    "CITY:", toupper(city_id),
    "| DBCV(sample n=", nrow(X_s), ", K=", length(unique(cl_s)), ") =",
    if (is.na(dbcv_val)) "NA" else sprintf("%.3f", dbcv_val),
    "\n"
  )
}

############################################################
# Step 9 — Cluster summary (size + means)
############################################################

cluster_summaries <- list()

for (city_id in names(prepared)) {

  city_dt <- prepared[[city_id]]$city_clean

  # Optional fields may differ by dataset source; keep summary robust.
  has_len <- "len" %in% names(city_dt)
  has_duration <- "duration" %in% names(city_dt)

  cluster_summary <- city_dt %>%
    dplyr::filter(cluster != 0) %>%            # exclude noise (DBSCAN label 0)
    dplyr::group_by(cluster) %>%
    dplyr::summarise(
      size = dplyr::n(),
      mean_x_o = mean(x_o, na.rm = TRUE),
      mean_y_o = mean(y_o, na.rm = TRUE),
      mean_dx = mean(dx, na.rm = TRUE),
      mean_dy = mean(dy, na.rm = TRUE),
      mean_len = if (has_len) mean(len, na.rm = TRUE) else NA_real_,
      mean_duration = if (has_duration) mean(duration, na.rm = TRUE) else NA_real_,
      .groups = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(size))          # sort by cluster size only

  cluster_summaries[[city_id]] <- cluster_summary

  cat("\n---", toupper(city_id), "---\n")
  if (!has_len || !has_duration) {
    cat("Note: len/duration missing; summary columns filled with NA where needed.\n")
  }
  print(head(cluster_summary, 15))
}

############################################################
# Step 10 — Analyze dominant cluster (but do NOT exclude)
############################################################

moo_clusters <- list()

for (city_id in names(cluster_summaries)) {

  cs <- cluster_summaries[[city_id]]

  dominant_cluster <- if (nrow(cs) > 0) cs$cluster[which.max(cs$size)] else NA_integer_

  cs_moo <- cs  # no filtering

  moo_clusters[[city_id]] <- cs_moo

  # Attach N_total once (so eval_solution can read it later)
  attr(moo_clusters[[city_id]], "N_total") <- nrow(prepared[[city_id]]$X)

  cat(
    city_id,
    "| dominant cluster:", if (is.na(dominant_cluster)) "NA" else dominant_cluster,
    "| clusters kept:", nrow(cs_moo),
    "| N_total =", attr(moo_clusters[[city_id]], "N_total"),
    "\n"
  )
}

############################################################
# Step 10a — Reference cluster count (context for MOO comparison)
############################################################

for (city_id in names(moo_clusters)) {

  cluster_ids <- moo_clusters[[city_id]]$cluster
  K <- length(cluster_ids)

  cat("CITY:", toupper(city_id), "| K =", K, "\n")
}

############################################################
# Step 11 — Directional similarity matrix (per city)
############################################################

direction_vectors <- list()
direction_similarity <- list()

for (city_id in names(moo_clusters)) {

  cs <- moo_clusters[[city_id]]

  # Unit direction vectors
  v <- as.matrix(cs[, c("mean_dx", "mean_dy")])
  norms <- sqrt(rowSums(v^2))
  norms[norms == 0] <- 1
  u <- v / norms

  # Cosine similarity matrix (K x K)
  S <- u %*% t(u)
  diag(S) <- NA

  direction_vectors[[city_id]] <- u
  direction_similarity[[city_id]] <- S

  cat(
    "CITY:", toupper(city_id),
    "| K =", nrow(S),
    "| similarity matrix computed\n"
  )
}

############################################################
# Step 11a — Unweighted direction concentration
# Measures pattern-level directional similarity:
# do selected clusters point in similar directions?
############################################################

dir_concentration_unweighted <- function(sel, u) {
  if (length(sel) < 2) return(0)

  m <- colMeans(u[sel, , drop = FALSE])  # equal weight per cluster
  sqrt(sum(m^2))                         # value in [0, 1]
}

############################################################
# Step 12 — Two objectives with data-normalized K penalty
# N_total is optional: read from df attribute if not provided
############################################################

eval_solution <- function(x, df, u, N_total = NULL) {

  if (is.null(N_total)) {
    N_total <- attr(df, "N_total")
  }
  if (is.null(N_total) || !is.numeric(N_total) || length(N_total) != 1 || N_total <= 0) {
    stop("eval_solution: N_total missing/invalid. Provide N_total or set attr(df,'N_total').")
  }

  sel <- which(x == 1)
  K <- length(sel)

  if (K == 0L) {
    return(c(
      f1_neg_size = 1,
      f2_neg_dir_concentration = 1,
      K = 0,
      coverage_raw = 0,
      penalty = 0
    ))
  }

  total_size <- sum(df$size[sel], na.rm = TRUE)
  coverage_raw <- total_size / N_total

  K_norm <- log(K + 1) / log(N_total + 1)
  penalty <- (1 - K_norm)

  coverage_eff <- coverage_raw * penalty

  f1 <- -coverage_eff
  dir_conc <- dir_concentration_unweighted(sel, u)
  f2 <- if (K < 2) 1 else -dir_conc

  c(
    f1_neg_size = f1,
    f2_neg_dir_concentration = f2,
    K = K,
    coverage_raw = coverage_raw,
    penalty = penalty
  )
}

############################################################
# Step 12a — Sanity checks (none/all) per city
############################################################

for (city_id in names(moo_clusters)) {

  df <- moo_clusters[[city_id]]
  u <- direction_vectors[[city_id]]
  K <- nrow(df)

  cat("\nCITY:", toupper(city_id), "\n")
  cat("K candidates in df:", K, "\n")

  # select none
  print(eval_solution(rep(0, K), df, u))

  # select all
  print(eval_solution(rep(1, K), df, u))
}

############################################################
# Step 12b — Random sanity checks (per city)
############################################################

set.seed(1)

for (city_id in names(moo_clusters)) {

  df <- moo_clusters[[city_id]]
  u <- direction_vectors[[city_id]]
  K <- nrow(df)

  cat("\nCITY:", toupper(city_id), "\n")

  if (K == 0) {
    cat("No clusters available; skipping random sanity checks.\n")
    next
  }

  for (i in 1:5) {
    x <- rbinom(K, 1, 0.2)

    # avoid trivial all-zero case
    if (sum(x) == 0) x[sample.int(K, 1)] <- 1

    res <- eval_solution(x, df, u)

    cat(
      "draw", i,
      "| selected =", sum(x),
      "| f1 =", round(res["f1_neg_size"], 6),
      "| f2 =", round(res["f2_neg_dir_concentration"], 6),
      "\n"
    )
  }
}

############################################################
# Step 13 — Scatter plot of random solutions (per city)
# Rebuild random samples using UPDATED eval_solution()
############################################################

set.seed(1)

vals_samples_new <- list()

for (city_id in names(moo_clusters)) {

  df <- moo_clusters[[city_id]]
  u <- direction_vectors[[city_id]]
  K <- nrow(df)

  # generate random solutions (adjust n_draws if you want)
  n_draws <- 300
  out <- data.frame(
    coverage_eff = numeric(n_draws),
    dir_concentration = numeric(n_draws),
    selected = integer(n_draws)
  )

  if (K == 0) {
    vals_samples_new[[city_id]] <- out
    cat("CITY:", toupper(city_id), "| no clusters; scatter skipped\n")
    next
  }

  for (i in 1:n_draws) {
    x <- rbinom(K, 1, 0.2)
    if (sum(x) == 0) x[sample.int(K, 1)] <- 1

    res <- eval_solution(x, df, u)

    out$coverage_eff[i] <- -as.numeric(res["f1_neg_size"])
    out$dir_concentration[i] <- -as.numeric(res["f2_neg_dir_concentration"])
    out$selected[i] <- sum(x)
  }

  vals_samples_new[[city_id]] <- out

  plot(
    out$coverage_eff,
    out$dir_concentration,
    xlab = "Effective coverage (coverage_raw × penalty)",
    ylab = "Direction concentration",
    main = paste("Random solutions |", toupper(city_id))
  )
}

############################################################
# Step 13a — SANITY CHECK: objective update consistency
############################################################

cat("\n=== Sanity check: N_total lookup ===\n")

# 1) N_total_by_city must exist and match prepared
if (!exists("N_total_by_city")) stop("N_total_by_city does not exist. Define it after Step 2.")
stopifnot(all(names(N_total_by_city) %in% names(prepared)))

for (city_id in names(prepared)) {
  n_now <- nrow(prepared[[city_id]]$X)
  n_ref <- as.numeric(N_total_by_city[[city_id]])
  cat("CITY:", toupper(city_id), "| prepared N_total =", n_now, "| lookup =", n_ref, "\n")
  stopifnot(n_now == n_ref)
}

cat("\n=== Sanity check: moo_clusters attributes + eval_solution ===\n")

for (city_id in names(moo_clusters)) {

  df0 <- moo_clusters[[city_id]]
  u <- direction_vectors[[city_id]]
  K <- nrow(df0)

  cat("\n---", toupper(city_id), "---\n")
  cat("K clusters:", K, "\n")

  # 2) Check attribute existence
  n_attr <- attr(df0, "N_total")
  cat("attr(df,'N_total'):", if (is.null(n_attr)) "NULL" else n_attr, "\n")

  # 3) Try eval_solution without passing N_total
  x_all <- rep(1, K)
  ok1 <- TRUE
  res1 <- tryCatch(eval_solution(x_all, df0, u), error = function(e) {
    ok1 <<- FALSE
    e
  })
  cat("eval_solution(x_all, df0, u) ->", if (ok1) "OK" else "ERROR", "\n")
  if (ok1) {
    cat(
      "  f1 =", res1["f1_neg_size"],
      "| f2 =", res1["f2_neg_dir_concentration"],
      "| coverage_raw =", res1["coverage_raw"],
      "| penalty =", res1["penalty"],
      "\n"
    )
  } else {
    cat("  message:", conditionMessage(res1), "\n")
  }

  # 4) Try eval_solution with explicit N_total from lookup (must always work)
  N_total <- as.numeric(N_total_by_city[[city_id]])
  ok2 <- TRUE
  res2 <- tryCatch(eval_solution(x_all, df0, u, N_total = N_total), error = function(e) {
    ok2 <<- FALSE
    e
  })
  cat("eval_solution(..., N_total=lookup) ->", if (ok2) "OK" else "ERROR", "\n")
  if (!ok2) stop("Even explicit N_total failed for ", city_id, ": ", conditionMessage(res2))

  # 5) Check if attribute gets stripped by common conversions
  df_plain <- as.data.frame(df0)  # may drop attributes
  n_attr2 <- attr(df_plain, "N_total")
  cat("After as.data.frame(): attr(df,'N_total'):", if (is.null(n_attr2)) "NULL" else n_attr2, "\n")

  ok3 <- TRUE
  res3 <- tryCatch(eval_solution(x_all, df_plain, u), error = function(e) {
    ok3 <<- FALSE
    e
  })
  cat("eval_solution(x_all, as.data.frame(df), u) ->", if (ok3) "OK" else "ERROR", "\n")
  if (!ok3) {
    cat("  (Expected if attributes were stripped) message:", conditionMessage(res3), "\n")
    cat("  Fix: pass N_total explicitly OR reattach attr(df,'N_total') before calling.\n")
  }

  # 6) Objective scale check (effective coverage should be between 0 and 1)
  f1_pos <- -as.numeric(res2["f1_neg_size"])
  if (f1_pos < -1e-9 || f1_pos > 1 + 1e-9) {
    cat("WARNING: effective coverage out of [0,1]:", f1_pos, "\n")
  } else {
    cat("effective coverage (=-f1) in [0,1]:", f1_pos, "\n")
  }
}

cat("\n=== All sanity checks completed ===\n")

############################################################
# Step 14' — Objective function for DBSCAN parameter optimization (UPDATED)
# Uses: effective coverage (clustered share × soft K-penalty) + directional cohesion
############################################################

eval_dbscan_params <- function(eps, minPts, X_scaled, city_dt) {
  # city_dt: data.table aligned with X_scaled and containing dx, dy

  db <- dbscan(X_scaled, eps = eps, minPts = minPts)
  cl <- db$cluster

  N_total <- length(cl)

  # clustered share (exclude noise)
  n_clustered <- sum(cl != 0)
  coverage_raw <- if (N_total > 0) n_clustered / N_total else 0

  # number of clusters produced (exclude noise)
  K <- length(unique(cl[cl != 0]))

  # same soft, data-normalized penalty used in previous objective block
  penalty <- if (K == 0) 0 else (1 - log(K + 1) / log(N_total + 1))
  coverage_eff <- coverage_raw * penalty

  # Need at least 2 non-noise clusters to define cohesion
  if (K < 2 || n_clustered == 0) {
    return(c(f1_neg_coverage = -coverage_eff, f2_neg_dir = 1))
  }

  # Attach cluster labels (in-place; avoids copying)
  city_dt[, cluster := cl]

  # Mean displacement per cluster (exclude noise)
  agg <- city_dt[cluster != 0 & is.finite(dx) & is.finite(dy),
                 .(mean_dx = mean(dx), mean_dy = mean(dy)),
                 by = cluster]

  # Safety: if aggregation collapses too much
  if (nrow(agg) < 2) {
    return(c(f1_neg_coverage = -coverage_eff, f2_neg_dir = 1))
  }

  v <- as.matrix(agg[, .(mean_dx, mean_dy)])
  norms <- sqrt(rowSums(v^2))
  norms[norms == 0] <- 1
  u <- v / norms

  # Use the same directional concentration helper style as previous blocks
  dir_cohesion <- dir_concentration_unweighted(seq_len(nrow(u)), u)

  c(
    f1_neg_coverage = -coverage_eff,
    f2_neg_dir = -dir_cohesion
  )
}
make_objfun_city_cached <- function(city_id, prepared,
                                    eps_min, eps_max, minPts_min, minPts_max,
                                    eps_step = 0.05,
                                    minPts_allowed = c(20L, 30L, 50L),
                                    eps_digits = 3) {

  X_scaled <- prepared[[city_id]]$X_scaled

  # Convert once (important for speed)
  city_dt <- as.data.table(prepared[[city_id]]$city_clean)

  # quick schema check
  if (!all(c("dx", "dy") %in% names(city_dt))) {
    stop("city_clean for ", city_id, " must contain dx and dy")
  }

  cache <- new.env(parent = emptyenv())

  function(x) {
    eps_raw <- max(eps_min, min(x[1], eps_max))
    eps <- round(eps_min + round((eps_raw - eps_min) / eps_step) * eps_step, eps_digits)
    eps <- max(eps_min, min(eps, eps_max))

    minPts_raw <- as.integer(round(x[2]))
    minPts_raw <- max(minPts_min, min(minPts_raw, minPts_max))
    minPts <- minPts_allowed[which.min(abs(minPts_allowed - minPts_raw))]

    key <- paste0(sprintf(paste0("%.", eps_digits, "f"), eps), "_", minPts)

    if (exists(key, envir = cache, inherits = FALSE)) {
      return(get(key, envir = cache, inherits = FALSE))
    }

    val <- eval_dbscan_params(eps, minPts, X_scaled, city_dt)
    assign(key, val, envir = cache)
    val
  }
}

cat("OK: Updated DBSCAN objective (effective coverage + cohesion) is defined.\n")

if (!requireNamespace("mco", quietly = TRUE)) {
  stop("Package 'mco' is required for NSGA-II. Install it with install.packages('mco').")
}

############################################################
# Step 14a — NSGA-II search on subsampled data (per city)
############################################################

set.seed(1)
N_SUB <- 10000  # 8k–15k is a good Colab range

prepared_sub <- lapply(prepared, function(res) {
  n <- nrow(res$X_scaled)
  idx <- if (n > N_SUB) sample.int(n, N_SUB) else seq_len(n)

  res$X_scaled <- res$X_scaled[idx, , drop = FALSE]
  res$city_clean <- res$city_clean[idx, ]
  res
})
data_for_opt <- prepared_sub

set.seed(1)

eps_min <- 0.30
eps_max <- 0.60
minPts_min <- 20
minPts_max <- 50
# align with earlier grid for fair comparison
eps_step <- 0.05
minPts_allowed <- c(20L, 30L, 50L)

popsize <- 40
generations <- 25

objfun_params <- list()
nsga_param_results <- list()

for (city_id in names(data_for_opt)) {

  cat("\nStarting:", toupper(city_id), "\n")
  t0 <- Sys.time()

  objfun_params[[city_id]] <- make_objfun_city_cached(
    city_id, data_for_opt,
    eps_min, eps_max, minPts_min, minPts_max,
    eps_step = eps_step,
    minPts_allowed = minPts_allowed
  )

  res <- mco::nsga2(
    fn = objfun_params[[city_id]],
    idim = 2,
    odim = 2,
    lower.bounds = c(eps_min, minPts_min),
    upper.bounds = c(eps_max, minPts_max),
    popsize = popsize,
    generations = generations
  )

  nsga_param_results[[city_id]] <- res

  cat(
    "Finished:", toupper(city_id),
    "| minutes =", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2),
    "\n"
  )
}

############################################################
# Step 14b — Pareto table + knee-point selection (per city)
############################################################

pareto_tables <- list()

for (city_id in names(nsga_param_results)) {

  res <- nsga_param_results[[city_id]]

  eps_raw <- res$par[, 1]
  minPts_raw <- as.integer(round(res$par[, 2]))

  # snap back to the same discrete search grid for fair comparison
  eps_snap <- round(eps_min + round((eps_raw - eps_min) / eps_step) * eps_step, 3)
  eps_snap <- pmax(eps_min, pmin(eps_snap, eps_max))
  minPts_raw <- pmax(minPts_min, pmin(minPts_raw, minPts_max))
  minPts_snap <- vapply(minPts_raw, function(v) minPts_allowed[which.min(abs(minPts_allowed - v))], integer(1))

  pareto_tables[[city_id]] <- data.frame(
    city = city_id,
    eps = eps_snap,
    minPts = as.integer(minPts_snap),
    coverage = -res$value[, 1],
    dir_cohesion = -res$value[, 2]
  )

  # filter fallback objective values (f2 = 1 => dir_cohesion = -1) before knee selection
  pareto_tables[[city_id]] <- pareto_tables[[city_id]][pareto_tables[[city_id]]$dir_cohesion >= 0, ]

  # remove duplicate parameter pairs
  pareto_tables[[city_id]] <- pareto_tables[[city_id]][
    !duplicated(pareto_tables[[city_id]][, c("eps", "minPts")]),
  ]

  cat("\n---", toupper(city_id), "Pareto solutions (head) ---\n")
  print(head(pareto_tables[[city_id]], 10))
}

knee_solutions <- list()

for (city_id in names(pareto_tables)) {

  df <- pareto_tables[[city_id]]

  if (nrow(df) == 0) {
    warning("No valid Pareto rows with non-negative dir_cohesion for ", city_id)
    knee_solutions[[city_id]] <- NA
    next
  }

  # Normalize to [0,1]
  cov_n <- (df$coverage - min(df$coverage)) / (max(df$coverage) - min(df$coverage) + 1e-9)
  coh_n <- (df$dir_cohesion - min(df$dir_cohesion)) / (max(df$dir_cohesion) - min(df$dir_cohesion) + 1e-9)

  # Distance to ideal (1,1)
  d <- sqrt((1 - cov_n)^2 + (1 - coh_n)^2)

  knee_solutions[[city_id]] <- df[which.min(d), ]

  cat("\nKNEE:", toupper(city_id), "\n")
  print(knee_solutions[[city_id]])
}
