# Prevent scientific notation globally
options(scipen = 999)

# Install vctrs first (do not load any other packages yet)
if (!requireNamespace("vctrs", quietly = TRUE)) {
  install.packages("vctrs")
}

# Install required packages for setup + city inspection
required_packages <- c("data.table", "dbscan", "ggplot2", "dplyr")
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
