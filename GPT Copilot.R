############################################################
# Step 0 — Integrated OD-prep pipeline for Berlin + Cologne + Munich
############################################################

pkgs <- c("data.table", "sf")
to_install <- pkgs[!vapply(pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
if (length(to_install)) install.packages(to_install)

library(data.table)
library(sf)

city_files <- list(
  berlin  = "dt_bolt_berlin_06_05.rds",
  cologne = "dt_voi_cologne_06_05.rds",
  munich  = "dt_voi_munich_06_05.rds"
)

needed <- c("start_loc_lon", "start_loc_lat", "dest_loc_lon", "dest_loc_lat")

prep_city_od <- function(rds_path,
                         needed = c("start_loc_lon","start_loc_lat","dest_loc_lon","dest_loc_lat"),
                         crs_in = 4326,
                         crs_out = 25832,
                         drop_zero_len = TRUE) {
  
  dt <- readRDS(rds_path)
  setDT(dt)
  
  if (!all(needed %in% names(dt))) {
    miss <- setdiff(needed, names(dt))
    stop("Missing columns in ", rds_path, ": ", paste(miss, collapse = ", "))
  }
  
  n_before <- nrow(dt)
  
  dt <- dt[
    !is.na(start_loc_lat) & !is.na(start_loc_lon) &
      !is.na(dest_loc_lat) & !is.na(dest_loc_lon)
  ]
  n_after_na <- nrow(dt)
  
  dt <- dt[
    start_loc_lat >= -90 & start_loc_lat <= 90 &
      dest_loc_lat  >= -90 & dest_loc_lat  <= 90 &
      start_loc_lon >= -180 & start_loc_lon <= 180 &
      dest_loc_lon  >= -180 & dest_loc_lon  <= 180
  ]
  n_after_range <- nrow(dt)
  
  orig_sf <- st_as_sf(
    dt,
    coords = c("start_loc_lon", "start_loc_lat"),
    crs = crs_in,
    remove = FALSE
  )
  
  dest_sf <- st_as_sf(
    dt,
    coords = c("dest_loc_lon", "dest_loc_lat"),
    crs = crs_in,
    remove = FALSE
  )
  
  orig_xy <- st_coordinates(st_transform(orig_sf, crs_out))
  dest_xy <- st_coordinates(st_transform(dest_sf, crs_out))
  
  dt[, `:=`(
    x_o = orig_xy[, 1],
    y_o = orig_xy[, 2],
    x_d = dest_xy[, 1],
    y_d = dest_xy[, 2]
  )]
  
  dt[, `:=`(dx = x_d - x_o, dy = y_d - y_o)]
  dt[, `:=`(
    angle = atan2(dy, dx),
    len   = sqrt(dx^2 + dy^2)
  )]
  
  n_zero_len <- dt[len == 0, .N]
  if (drop_zero_len) dt <- dt[len > 0]
  n_after_zero <- nrow(dt)
  
  info <- list(
    file = rds_path,
    n_before = n_before,
    n_after_na = n_after_na,
    n_after_range = n_after_range,
    n_zero_len = n_zero_len,
    n_after_zero = n_after_zero
  )
  
  list(data = dt, info = info)
}

prepared <- lapply(city_files, prep_city_od, needed = needed, drop_zero_len = TRUE)

cleaning_summary <- rbindlist(lapply(names(prepared), function(city) {
  i <- prepared[[city]]$info
  data.table(
    city = city,
    file = i$file,
    n_before = i$n_before,
    n_after_na = i$n_after_na,
    n_after_range = i$n_after_range,
    zero_len_trips = i$n_zero_len,
    n_after_zero_len_removed = i$n_after_zero
  )
}), fill = TRUE)

print(cleaning_summary)

berlin  <- prepared$berlin$data
cologne <- prepared$cologne$data
munich  <- prepared$munich$data

cat("\nZero-length remaining (should be 0):\n")
cat("Berlin :", sum(berlin$len == 0), "\n")
cat("Cologne:", sum(cologne$len == 0), "\n")
cat("Munich :", sum(munich$len == 0), "\n")

############################################################
# Step 1 — Load + inspect all cities (no combining)
# BRIDGE: build city_data from Step 0 output (prepared)
############################################################

# Keep a stable pointer because Step 2 will reassign `prepared`
prepared_od <- prepared

# Store each city separately (no combining) — SAME OBJECT NAME: city_data
city_data <- lapply(names(prepared_od), function(city_id) {
  dt <- prepared_od[[city_id]]$data
  setDT(dt)
  if (!("city" %in% names(dt))) dt[, city := city_id]  # traceability tag
  dt
})
names(city_data) <- names(prepared_od)

# Print inspection for each city (your original loop)
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
############################################################

make_X_and_clean <- function(city_dt,
                             feature_cols = c("x_o", "y_o", "dx", "dy")) {
  
  missing_cols <- setdiff(feature_cols, names(city_dt))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  X <- as.matrix(city_dt[, ..feature_cols])
  
  ok <- complete.cases(X)
  X_scaled <- scale(X[ok, , drop = FALSE])
  city_clean <- city_dt[ok, ]
  
  list(X = X[ok, , drop = FALSE], ok = ok, X_scaled = X_scaled, city_clean = city_clean)
}

# Run for all cities (no combining) — SAME OBJECT NAME: prepared (Step 2 output)
prepared <- lapply(names(city_data), function(city_id) {
  res <- make_X_and_clean(city_data[[city_id]])
  res$city_id <- city_id
  res
})
names(prepared) <- names(city_data)

N_total_by_city <- sapply(names(prepared), function(city_id) {
  nrow(prepared[[city_id]]$X)
})
print(N_total_by_city)

berlin_clean  <- prepared$berlin$city_clean
munich_clean  <- prepared$munich$city_clean
cologne_clean <- prepared$cologne$city_clean

X_scaled_berlin  <- prepared$berlin$X_scaled
X_scaled_munich  <- prepared$munich$X_scaled
X_scaled_cologne <- prepared$cologne$X_scaled

for (city_id in names(prepared)) {
  cat("\n---", toupper(city_id), "---\n")
  cat("Rows original:", nrow(city_data[[city_id]]), "\n")
  cat("Rows kept:",     nrow(prepared[[city_id]]$city_clean), "\n")
  cat("Rows dropped:",  sum(!prepared[[city_id]]$ok), "\n")
  cat("X_scaled dim:",  paste(dim(prepared[[city_id]]$X_scaled), collapse = " x "), "\n")
}


############################################################
# Search-space justification diagnostics (run after `prepared`)
############################################################

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(dbscan)
  library(ggplot2)
})

# -----------------------------
# Preconditions
# -----------------------------
if (!exists("prepared") || !is.list(prepared) || length(prepared) == 0) {
  stop("Object `prepared` not found or not a non-empty list. Run data preparation first.")
}
for (city_id in names(prepared)) {
  if (is.null(prepared[[city_id]]$X_scaled)) stop("Missing X_scaled for city: ", city_id)
}

# -----------------------------
# Helpers
# -----------------------------

knee_elbow <- function(d_sorted) {
  d <- as.numeric(d_sorted)
  n <- length(d)
  if (n < 10) return(NA_real_)
  x <- seq_len(n)
  
  x1 <- 1; y1 <- d[1]
  x2 <- n; y2 <- d[n]
  
  num <- abs((y2 - y1) * x - (x2 - x1) * d + x2*y1 - y2*x1)
  den <- sqrt((y2 - y1)^2 + (x2 - x1)^2)
  dist <- num / den
  
  d[which.max(dist)]
}

knn_summary_one <- function(X_scaled, k,
                            q = c(0.01,0.05,0.10,0.25,0.50,0.75,0.90,0.95,0.99)) {
  d <- dbscan::kNNdist(X_scaled, k = k)
  ds <- sort(d)
  
  out <- data.frame(
    k = k,
    n = length(d),
    elbow_eps = knee_elbow(ds),
    t(quantile(d, probs = q, na.rm = TRUE, names = FALSE))
  )
  names(out) <- c("k","n","elbow_eps", paste0("q", sprintf("%02d", round(q*100))))
  out
}

snn_sharednn_diagnostics <- function(X_scaled, k = 30, n_probe = 4000, seed = 1) {
  set.seed(seed)
  n <- nrow(X_scaled)
  idx <- if (n > n_probe) sample.int(n, n_probe) else seq_len(n)
  
  knn <- dbscan::kNN(X_scaled, k = k)
  nn_list <- lapply(idx, function(i) knn$id[i, ])
  
  m <- length(nn_list)
  checks <- min(8000, m * 2)
  shared_counts <- integer(checks)
  
  for (t in seq_len(checks)) {
    a <- sample.int(m, 1)
    b <- sample.int(m, 1)
    shared_counts[t] <- length(intersect(nn_list[[a]], nn_list[[b]]))
  }
  
  data.frame(
    k = k,
    probe_n = m,
    shared_mean = mean(shared_counts),
    shared_median = median(shared_counts),
    shared_q90 = unname(quantile(shared_counts, 0.90)),
    shared_q95 = unname(quantile(shared_counts, 0.95)),
    shared_q99 = unname(quantile(shared_counts, 0.99))
  )
}

hdbscan_mcs_candidates <- function(N_eff, n = 14, frac_min = 0.005, frac_max = 0.10) {
  lo <- max(10L, as.integer(round(frac_min * N_eff)))
  hi <- max(lo + 1L, as.integer(round(frac_max * N_eff)))
  vals <- unique(as.integer(round(exp(seq(log(lo), log(hi), length.out = n)))))
  sort(vals)
}

# -----------------------------
# Config
# -----------------------------
k_values_dbscan <- c(10, 20, 30, 50)

minPts_allowed <- c(6L, 8L, 10L, 12L, 15L, 20L, 25L, 30L, 40L, 50L, 60L, 70L)

k_allowed_snn <- c(15L, 20L, 25L, 30L, 35L)
n_probe_snn <- 4000

min_samples_allowed_hdb <- c(10, 15, 20, 25, 30, 35, 40, 50, 60)

# -----------------------------
# 1) DBSCAN eps range justification (kNN distances)
# -----------------------------
knn_tables <- list()
eps_reco <- list()

for (city_id in names(prepared)) {
  X_scaled <- prepared[[city_id]]$X_scaled
  
  tab_city <- bind_rows(lapply(k_values_dbscan, function(k) knn_summary_one(X_scaled, k))) %>%
    mutate(city = city_id) %>%
    relocate(city)
  
  knn_tables[[city_id]] <- tab_city
  
  eps_min <- min(tab_city$q05, na.rm = TRUE)
  eps_max <- max(tab_city$q99, na.rm = TRUE)
  elbow_med <- median(tab_city$elbow_eps, na.rm = TRUE)
  
  eps_reco[[city_id]] <- data.frame(
    city = city_id,
    eps_min_q05 = eps_min,
    eps_max_q99 = eps_max,
    elbow_median = elbow_med
  )
}

knn_df <- bind_rows(knn_tables)
eps_reco_df <- bind_rows(eps_reco)

cat("\n=== kNN distance summaries (per city, per k) ===\n")
print(knn_df)

cat("\n=== Recommended eps bounds from kNN distances ===\n")
print(eps_reco_df)

plot_knn_curves <- function(city_id, k_values = k_values_dbscan) {
  X_scaled <- prepared[[city_id]]$X_scaled
  dt <- rbindlist(lapply(k_values, function(k) {
    d <- sort(dbscan::kNNdist(X_scaled, k = k))
    data.table(city = city_id, k = k, rank = seq_along(d), dist = d)
  }))
  ggplot(dt, aes(x = rank, y = dist, color = factor(k))) +
    geom_line(alpha = 0.8) +
    labs(
      title = paste0("kNN-distance curves — ", toupper(city_id)),
      x = "Points sorted by distance",
      y = "kNN distance",
      color = "k"
    ) +
    theme_minimal()
}

# Example:
# print(plot_knn_curves("berlin"))

# -----------------------------
# 2) SNN bounds sanity check
# -----------------------------
snn_diag <- list()
for (city_id in names(prepared)) {
  X_scaled <- prepared[[city_id]]$X_scaled
  
  diag_city <- bind_rows(lapply(k_allowed_snn, function(k) {
    snn_sharednn_diagnostics(X_scaled, k = k, n_probe = n_probe_snn, seed = 1) %>%
      mutate(city = city_id)
  })) %>%
    relocate(city)
  
  snn_diag[[city_id]] <- diag_city
}

snn_diag_df <- bind_rows(snn_diag)

cat("\n=== SNN shared-neighbor diagnostics (probe sample) ===\n")
print(snn_diag_df)

cat("\nInterpretation guide:\n")
cat("- eps_snn should not be near 0 (over-connect) or near k (over-fragment).\n")
cat("- shared_q90/shared_q95 suggest where meaningful overlap begins per city.\n")

# -----------------------------
# 3) HDBSCAN bounds justification (from N_SUB used in optimization)
# -----------------------------
N_SUB <- 10000L
mcs_vals <- hdbscan_mcs_candidates(N_eff = N_SUB, n = 14, frac_min = 0.005, frac_max = 0.10)

hdb_bounds_df <- data.frame(
  N_eff = N_SUB,
  frac_min = 0.005,
  frac_max = 0.10,
  mcs_min = min(mcs_vals),
  mcs_max = max(mcs_vals),
  mcs_values = paste(mcs_vals, collapse = ", "),
  min_samples_allowed = paste(min_samples_allowed_hdb, collapse = ", ")
)

cat("\n=== HDBSCAN search-space implied by fraction rule (on subsample) ===\n")
print(hdb_bounds_df)

# -----------------------------
# 4) Optional export
# -----------------------------
# write.csv(knn_df, "knn_distance_summaries.csv", row.names = FALSE)
# write.csv(eps_reco_df, "dbscan_eps_recommendations.csv", row.names = FALSE)
# write.csv(snn_diag_df, "snn_sharednn_diagnostics.csv", row.names = FALSE)
# write.csv(hdb_bounds_df, "hdbscan_bounds_summary.csv", row.names = FALSE)
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
# - 20–30: moderate neighborhood size (often stable)
# - 50   : stress-test for over-smoothing
k_values <- c(10, 20, 30, 50)

# Save current graphical parameters
old_par <- par(no.readonly = TRUE)

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

# Restore original graphical parameters
par(old_par)

suppressPackageStartupMessages({
  library(dbscan)   # kNN()
})

k_values <- c(10, 20, 30, 50)

# Allocation-free intersection size for sorted integer vectors
# (two-pointer scan; MUCH faster than intersect())
intersect_count_sorted <- function(a, b) {
  i <- 1L; j <- 1L
  na <- length(a); nb <- length(b)
  cnt <- 0L
  while (i <= na && j <= nb) {
    av <- a[i]; bv <- b[j]
    if (av == bv) { cnt <- cnt + 1L; i <- i + 1L; j <- j + 1L }
    else if (av < bv) i <- i + 1L
    else j <- j + 1L
  }
  cnt
}

# SNN diagnostic curve:
# For each point i: shared(i, k-th NN of i)
# Sort across points to get a curve for choosing snn_threshold.
snn_kth_shared_curve <- function(X_scaled, k) {
  knn <- dbscan::kNN(X_scaled, k = k)
  nn_id <- knn$id  # n x k

  # Ensure rows are sorted to make intersection counting valid/fast.
  # (kNN() usually returns neighbors ordered by distance already,
  #  but sorting is cheap and guarantees correctness.)
  nn_id <- t(apply(nn_id, 1, sort))

  n <- nrow(nn_id)
  shared_kth <- integer(n)

  for (i in seq_len(n)) {
    Ni <- nn_id[i, ]
    j  <- Ni[k]               # k-th nearest neighbor of i
    Nj <- nn_id[j, ]
    shared_kth[i] <- intersect_count_sorted(Ni, Nj)
  }

  sort(shared_kth)
}

# Plotting
old_par <- par(no.readonly = TRUE)

for (city_id in names(prepared)) {

  X_scaled <- prepared[[city_id]]$X_scaled  # enforce same input as DBSCAN

  cat("\nPlotting SNN shared-neighbor curves for:", toupper(city_id), "\n")

  old_par_city <- par(no.readonly = TRUE)
  par(mfrow = c(2, 2))

  for (k in k_values) {
    s_curve <- snn_kth_shared_curve(X_scaled, k = k)

    plot(s_curve, type = "l",
         xlab = "Points sorted by shared neighbors",
         ylab = paste0("shared neighbors with ", k, "-th NN (k=", k, ")"),
         main = paste(city_id, "| SNN diagnostic (k =", k, ")"))
    grid()
  }

  par(old_par_city)
}

par(old_par)
                             
############################################################
# Step 4 — DBSCAN parameter grids (FINAL, MOO-fair)
# DBSCAN params:
# - eps    : neighborhood radius (distance scale)
# - minPts : core density requirement (aligned with Step 3 k diagnostics)
############################################################

# Shared locality scales (use same values diagnosed in Step 3)
k_grid_core <- c(10, 20, 30)

# eps candidates (based on your k-distance knee range)
eps_grid <- seq(0.08, 0.40, by = 0.02)   # 17 values

# minPts candidates aligned with Step 3 k diagnostics
minPts_grid_dbscan <- k_grid_core       # 4 values

# Total DBSCAN configs = length(eps_grid) * length(minPts_grid_dbscan) = 68

############################################################
# Step 4 — SNN parameter grids (UPDATED, MOO-fair + covers final configs)
# SNN params:
# - kNN size             : k
# - shared-neighbor min  : snn_threshold (integer, 0..k-1)
# - minPts               : core requirement in SNN clustering
#
# Rationale:
# - Stage-1 runs showed two failure modes:
#   (i) over-merge (few clusters, ~0% noise) at large k / low threshold
#   (ii) fragmentation / high noise at too strict settings
# - This grid includes a stable mid-regime observed in full-data runs.
############################################################

# k values: include 15 (needed for dense Cologne), keep moderate for Colab
k_grid_snn <- c(15, 20, 30)

# Threshold grid depends on k (bounded by k-1).
# Use a wider band that includes the empirically stable region (~40–55% of k)
# while still exploring stricter values.
snn_threshold_grid <- function(k) {
  thr <- round(seq(0.35 * k, 0.70 * k, length.out = 7))
  thr <- unique(pmax(2L, pmin(k - 1L, as.integer(thr))))
  thr
}

# minPts: include 10 and 12 (needed for stable Cologne/Berlin),
# but keep upper values for robustness checks.
minPts_grid_snn <- c(10, 12, 15, 18, 20)

############################################################
# Step 5 — Initialize results container (DBSCAN)
############################################################

results_dbscan <- list()
idx_dbscan <- 1

############################################################
# Helper — direction head (prefer columns in X_scaled / X)
############################################################
get_dir_head <- function(prep_city) {
  if (!is.null(prep_city$X_scaled)) {
    cn <- colnames(prep_city$X_scaled)
    if (!is.null(cn) && all(c("dx","dy") %in% cn)) {
      return(head(as.data.frame(prep_city$X_scaled[, c("dx","dy"), drop = FALSE])))
    }
    if (!is.null(cn) && all(c("ux","uy") %in% cn)) {
      return(head(as.data.frame(prep_city$X_scaled[, c("ux","uy"), drop = FALSE])))
    }
  }
  if (!is.null(prep_city$X)) {
    cn <- colnames(prep_city$X)
    if (!is.null(cn) && all(c("dx","dy") %in% cn)) {
      return(head(as.data.frame(prep_city$X[, c("dx","dy"), drop = FALSE])))
    }
    if (!is.null(cn) && all(c("ux","uy") %in% cn)) {
      return(head(as.data.frame(prep_city$X[, c("ux","uy"), drop = FALSE])))
    }
  }
  if (!is.null(prep_city$ux) && !is.null(prep_city$uy)) {
    return(head(data.frame(ux = prep_city$ux, uy = prep_city$uy)))
  }
  if (!is.null(prep_city$dx) && !is.null(prep_city$dy)) {
    return(head(data.frame(dx = prep_city$dx, dy = prep_city$dy)))
  }
  NULL
}

############################################################
# Step 5a — Sanity checks / debugging prints (DBSCAN)
############################################################

cat("\n[DBSCAN] idx =", idx_dbscan, "\n")

cat("\n[DBSCAN] Parameter grids:\n")
print(eps_grid)
print(minPts_grid_dbscan)

cat("\n[DBSCAN] results_dbscan length =", length(results_dbscan), "\n")

# Check feature matrices per city (ENFORCE X_scaled)
for (city_id in names(prepared)) {
  cat("\n---", toupper(city_id), "---\n")

  if (is.null(prepared[[city_id]]$X_scaled)) {
    stop("[DBSCAN] Missing X_scaled for city: ", city_id)
  }

  X_scaled <- prepared[[city_id]]$X_scaled
  cat("X_scaled dims:", paste(dim(X_scaled), collapse = " x "), "\n")
  cat("X_scaled head:\n")
  print(head(X_scaled))

  dir_head <- get_dir_head(prepared[[city_id]])
  if (!is.null(dir_head)) {
    cat("Direction head:\n")
    print(dir_head)
  } else {
    cat("WARNING: direction vectors not found (dx/dy or ux/uy missing).\n")
  }
}

############################################################
# Step 5 — Initialize results container (SNN)
############################################################

results_snn <- list()
idx_snn <- 1

############################################################
# Helper — direction head (prefer columns in X_scaled / X)
############################################################
get_dir_head <- function(prep_city) {
  if (!is.null(prep_city$X_scaled)) {
    cn <- colnames(prep_city$X_scaled)
    if (!is.null(cn) && all(c("dx","dy") %in% cn)) {
      return(head(as.data.frame(prep_city$X_scaled[, c("dx","dy"), drop = FALSE])))
    }
    if (!is.null(cn) && all(c("ux","uy") %in% cn)) {
      return(head(as.data.frame(prep_city$X_scaled[, c("ux","uy"), drop = FALSE])))
    }
  }
  if (!is.null(prep_city$X)) {
    cn <- colnames(prep_city$X)
    if (!is.null(cn) && all(c("dx","dy") %in% cn)) {
      return(head(as.data.frame(prep_city$X[, c("dx","dy"), drop = FALSE])))
    }
    if (!is.null(cn) && all(c("ux","uy") %in% cn)) {
      return(head(as.data.frame(prep_city$X[, c("ux","uy"), drop = FALSE])))
    }
  }
  if (!is.null(prep_city$ux) && !is.null(prep_city$uy)) {
    return(head(data.frame(ux = prep_city$ux, uy = prep_city$uy)))
  }
  if (!is.null(prep_city$dx) && !is.null(prep_city$dy)) {
    return(head(data.frame(dx = prep_city$dx, dy = prep_city$dy)))
  }
  NULL
}

############################################################
# Step 5a — Sanity checks / debugging prints (SNN)
############################################################

cat("\n[SNN] idx =", idx_snn, "\n")

cat("\n[SNN] Parameter grids:\n")
print(k_grid_snn)
print(snn_threshold_grid)
print(minPts_grid_snn)

cat("\n[SNN] results_snn length =", length(results_snn), "\n")

# Check feature matrices per city (ENFORCE X_scaled)
for (city_id in names(prepared)) {
  cat("\n---", toupper(city_id), "---\n")

  if (is.null(prepared[[city_id]]$X_scaled)) {
    stop("[SNN] Missing X_scaled for city: ", city_id)
  }

  X_scaled <- prepared[[city_id]]$X_scaled
  cat("X_scaled dims:", paste(dim(X_scaled), collapse = " x "), "\n")
  cat("X_scaled head:\n")
  print(head(X_scaled))

  dir_head <- get_dir_head(prepared[[city_id]])
  if (!is.null(dir_head)) {
    cat("Direction head:\n")
    print(dir_head)
  } else {
    cat("WARNING: direction vectors not found (dx/dy or ux/uy missing).\n")
  }
}

############################################################
# Step 6 — DBSCAN grid search (FULL DATA, X_scaled)
############################################################

suppressPackageStartupMessages({
  library(dbscan)
  library(dplyr)
})

results_dbscan <- list()
idx_dbscan <- 1

for (city_id in names(prepared)) {

  if (is.null(prepared[[city_id]]$X_scaled)) {
    stop("[DBSCAN] Missing X_scaled for city: ", city_id)
  }

  X_scaled <- prepared[[city_id]]$X_scaled
  n_total  <- nrow(X_scaled)

  for (minPts in minPts_grid_dbscan) {
    for (eps_val in eps_grid) {

      fit <- dbscan::dbscan(X_scaled, eps = eps_val, minPts = minPts)
      cl  <- fit$cluster  # 0 = noise

      n_noise   <- sum(cl == 0)
      noise_pct <- 100 * n_noise / length(cl)

      # number of non-noise clusters
      n_clusters <- length(unique(cl[cl != 0]))

      max_cluster_size <- if (n_clusters > 0) max(table(cl[cl != 0])) else 0
      max_cluster_pct  <- if (n_clusters > 0) 100 * max_cluster_size / length(cl) else 0

      results_dbscan[[idx_dbscan]] <- data.frame(
        city = city_id,
        algorithm = "DBSCAN",
        minPts = minPts,
        eps = eps_val,
        n_clusters = n_clusters,
        noise_pct = noise_pct,
        max_cluster_pct = max_cluster_pct,
        n_total = n_total,
        n_used = n_total,
        stringsAsFactors = FALSE
      )

      idx_dbscan <- idx_dbscan + 1
    }
  }
}

grid_df_dbscan <- dplyr::bind_rows(results_dbscan)
print(grid_df_dbscan)

############################################################
# Step 6 — SNN grid search (FULL DATA; cached sNN per (city,k); jp=FALSE)
############################################################

suppressPackageStartupMessages({
  library(dbscan)
  library(dplyr)
})

# Colab-friendly k values (drop 50)
k_grid_snn <- c(10, 20, 30)

# Threshold grid scaled with k (shared neighbors in [0, k-1])
make_snn_threshold_grid <- function(k) {
  thr <- round(seq(0.35 * k, 0.65 * k, length.out = 6))
  thr <- unique(pmax(2L, pmin(k - 1L, as.integer(thr))))
  thr
}

# Your useful range
minPts_grid_snn <- c(15, 18, 20, 25)

results_snn <- list()
idx_snn <- 1

for (city_id in names(prepared)) {

  if (is.null(prepared[[city_id]]$X_scaled)) {
    stop("[SNN] Missing X_scaled for city: ", city_id)
  }

  X_scaled <- prepared[[city_id]]$X_scaled
  n_total  <- nrow(X_scaled)

  cat("\n=== ", toupper(city_id), " | n = ", n_total, " (FULL DATA) ===\n", sep = "")

  for (k in k_grid_snn) {

    snn_threshold_grid <- make_snn_threshold_grid(k)
    cat("  k =", k, " | snn_threshold_grid =", paste(snn_threshold_grid, collapse = ","), "\n")

    # Build once per (city,k) and reuse (FAST)
    snn_obj <- dbscan::sNN(X_scaled, k = k, sort = FALSE, jp = FALSE)

    for (snn_thr in snn_threshold_grid) {
      for (minPts in minPts_grid_snn) {

        fit <- dbscan::sNNclust(
          x = snn_obj,
          k = k,
          eps = snn_thr,        # shared-neighbor threshold (integer)
          minPts = minPts,
          borderPoints = TRUE
        )

        cl <- fit$cluster  # 0 = noise

        n_noise   <- sum(cl == 0)
        noise_pct <- 100 * n_noise / length(cl)

        n_clusters <- length(unique(cl[cl != 0]))

        max_cluster_size <- if (n_clusters > 0) max(table(cl[cl != 0])) else 0
        max_cluster_pct  <- if (n_clusters > 0) 100 * max_cluster_size / length(cl) else 0

        results_snn[[idx_snn]] <- data.frame(
          city = city_id,
          algorithm = "SNN",
          jp = FALSE,
          k = k,
          snn_threshold = snn_thr,
          minPts = minPts,
          eps = NA_real_,              # alignment column for DBSCAN
          n_clusters = n_clusters,
          noise_pct = noise_pct,
          max_cluster_pct = max_cluster_pct,
          n_total = n_total,
          n_used = n_total,
          stringsAsFactors = FALSE
        )

        idx_snn <- idx_snn + 1
      }
    }
  }
}

grid_df_snn <- dplyr::bind_rows(results_snn)
print(grid_df_snn)

############################################################
# Step 7 — Final DBSCAN parameters (per city)
############################################################

dbscan_params <- list(
  berlin  = list(minPts = 20, eps = 0.200),
  munich  = list(minPts = 20, eps = 0.200),
  cologne = list(minPts = 10, eps = 0.250)
)

############################################################
# Step 7a — Final DBSCAN run (all cities) — X_scaled only
############################################################

suppressPackageStartupMessages(library(dbscan))

for (city_id in names(prepared)) {

  if (is.null(prepared[[city_id]]$X_scaled)) {
    stop("[DBSCAN] Missing X_scaled for city: ", city_id)
  }
  if (is.null(prepared[[city_id]]$city_clean)) {
    stop("[DBSCAN] Missing city_clean for city: ", city_id)
  }

  X_scaled <- prepared[[city_id]]$X_scaled
  city_tbl <- prepared[[city_id]]$city_clean
  stopifnot(nrow(X_scaled) == nrow(city_tbl))

  params <- dbscan_params[[city_id]]

  fit <- dbscan::dbscan(X_scaled, eps = params$eps, minPts = params$minPts)

  stopifnot(length(fit$cluster) == nrow(city_tbl))
  prepared[[city_id]]$city_clean$cluster <- fit$cluster

  cat(
    "DBSCAN done |", city_id,
    "| eps =", params$eps,
    "| minPts =", params$minPts,
    "| clusters =", length(unique(fit$cluster[fit$cluster != 0])),
    "| noise% =", round(100 * mean(fit$cluster == 0), 2),
    "\n"
  )
}

############################################################
# Step 7 — Final SNN parameters (per city) — jp=FALSE
############################################################
snn_params <- list(
  berlin  = list(k = 20, snn_threshold = 8,  minPts = 12),
  munich  = list(k = 30, snn_threshold = 11, minPts = 18),
  cologne = list(k = 15, snn_threshold = 6,  minPts = 10)
)
############################################################
# Step 7a — Final SNN run (all cities) — X_scaled only, jp=FALSE
############################################################

suppressPackageStartupMessages(library(dbscan))

for (city_id in names(prepared)) {

  if (is.null(prepared[[city_id]]$X_scaled)) {
    stop("[SNN] Missing X_scaled for city: ", city_id)
  }
  if (is.null(prepared[[city_id]]$city_clean)) {
    stop("[SNN] Missing city_clean for city: ", city_id)
  }

  X_scaled <- prepared[[city_id]]$X_scaled
  city_tbl <- prepared[[city_id]]$city_clean
  stopifnot(nrow(X_scaled) == nrow(city_tbl))

  params <- snn_params[[city_id]]

  snn_obj <- dbscan::sNN(X_scaled, k = params$k, sort = FALSE, jp = FALSE)

  fit <- dbscan::sNNclust(
    x = snn_obj,
    k = params$k,
    eps = params$snn_threshold,
    minPts = params$minPts,
    borderPoints = TRUE
  )

  stopifnot(length(fit$cluster) == nrow(city_tbl))
  prepared[[city_id]]$city_clean$cluster <- fit$cluster

  cat(
    "SNN done |", city_id,
    "| jp = FALSE",
    "| k =", params$k,
    "| snn_threshold =", params$snn_threshold,
    "| minPts =", params$minPts,
    "| clusters =", length(unique(fit$cluster[fit$cluster != 0])),
    "| noise% =", round(100 * mean(fit$cluster == 0), 2),
    "\n"
  )
}

############################################################
# Step 8 — Verify final clustering (cluster size check)
# (Unified block for DBSCAN / SNN / HDBSCAN)
############################################################

for (city_id in names(prepared)) {

  cat("\n---", toupper(city_id), "---\n")

  if (is.null(prepared[[city_id]]$city_clean)) {
    stop("[Step 8] Missing city_clean for city: ", city_id)
  }
  if (is.null(prepared[[city_id]]$city_clean$cluster)) {
    stop("[Step 8] Missing cluster labels in city_clean for city: ", city_id)
  }

  cl_tab <- table(prepared[[city_id]]$city_clean$cluster)

  # Print first clusters (including noise = 0)
  print(head(cl_tab, 65))
}

############################################################
# Step 9 — Cluster summary (size + means) — UNIFIED
# - Works for DBSCAN/SNN/HDBSCAN labels (0 = noise)
# - Assumes final cluster labels are stored in prepared[[city]]$city_clean$cluster
# - Robust to optional columns (len, duration) and minor naming differences
############################################################

suppressPackageStartupMessages({
  library(dplyr)
})

# helper: pick first existing column name from a set
pick_col <- function(df, candidates) {
  hit <- candidates[candidates %in% names(df)]
  if (length(hit) == 0) return(NA_character_)
  hit[[1]]
}

cluster_summaries <- list()

for (city_id in names(prepared)) {

  if (is.null(prepared[[city_id]]$city_clean)) {
    stop("[Step 9] Missing city_clean for city: ", city_id)
  }

  city_dt <- prepared[[city_id]]$city_clean

  # need cluster labels
  if (!("cluster" %in% names(city_dt))) {
    stop("[Step 9] Missing cluster column in city_clean for city: ", city_id)
  }

  # column mapping (supports cleaned or raw naming)
  col_xo <- pick_col(city_dt, c("x_o", "origin_lon"))
  col_yo <- pick_col(city_dt, c("y_o", "origin_lat"))
  col_dx <- pick_col(city_dt, c("dx", "ux", "delta_x", "d_lon"))
  col_dy <- pick_col(city_dt, c("dy", "uy", "delta_y", "d_lat"))

  # optional fields
  has_len <- "len" %in% names(city_dt)
  has_duration <- "duration" %in% names(city_dt)

  cluster_summary <- city_dt %>%
    filter(cluster != 0) %>%                 # exclude noise
    group_by(cluster) %>%
    summarise(
      size = n(),
      mean_x_o = if (!is.na(col_xo)) mean(.data[[col_xo]], na.rm = TRUE) else NA_real_,
      mean_y_o = if (!is.na(col_yo)) mean(.data[[col_yo]], na.rm = TRUE) else NA_real_,
      mean_dx  = if (!is.na(col_dx)) mean(.data[[col_dx]], na.rm = TRUE) else NA_real_,
      mean_dy  = if (!is.na(col_dy)) mean(.data[[col_dy]], na.rm = TRUE) else NA_real_,
      mean_len = if (has_len) mean(len, na.rm = TRUE) else NA_real_,
      mean_duration = if (has_duration) mean(duration, na.rm = TRUE) else NA_real_,
      .groups = "drop"
    ) %>%
    arrange(desc(size))

  cluster_summaries[[city_id]] <- cluster_summary

  cat("\n---", toupper(city_id), "---\n")
  if (is.na(col_xo) || is.na(col_yo) || is.na(col_dx) || is.na(col_dy)) {
    cat("Note: some expected coordinate/direction columns missing; filled with NA where needed.\n")
  }
  if (!has_len || !has_duration) {
    cat("Note: len/duration missing; summary columns filled with NA where needed.\n")
  }
  print(head(cluster_summary, 15))
}

############################################################
# Step 10 — Analyze dominant cluster (but do NOT exclude) — UNIFIED
# - uses cluster_summaries created from city_clean
# - attaches N_total consistently using X_scaled row count
############################################################

moo_clusters <- list()

for (city_id in names(cluster_summaries)) {

  cs <- cluster_summaries[[city_id]]

  dominant_cluster <- if (nrow(cs) > 0) cs$cluster[which.max(cs$size)] else NA_integer_

  cs_moo <- cs  # keep all clusters (no filtering)
  moo_clusters[[city_id]] <- cs_moo

  # N_total must match the feature matrix used for clustering/evaluation
  if (is.null(prepared[[city_id]]$X_scaled)) {
    stop("[Step 10] Missing X_scaled for city: ", city_id)
  }
  attr(moo_clusters[[city_id]], "N_total") <- nrow(prepared[[city_id]]$X_scaled)

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
  K <- nrow(moo_clusters[[city_id]])  # one row per cluster (noise already excluded in Step 9)
  cat("CITY:", toupper(city_id), "| K =", K, "\n")
}

############################################################
# Step 11 — Trip-level unit direction vectors (FULL, aligned) — UNIFIED
############################################################

trip_dirs_full <- lapply(names(prepared), function(city_id) {

  # Prefer X_scaled direction columns (same matrix used in clustering)
  if (!is.null(prepared[[city_id]]$X_scaled)) {
    Xs <- prepared[[city_id]]$X_scaled
    cn <- colnames(Xs)
    if (!is.null(cn) && all(c("dx","dy") %in% cn)) {
      dx <- as.numeric(Xs[, "dx"])
      dy <- as.numeric(Xs[, "dy"])
    } else {
      # fall back to table
      dt <- prepared[[city_id]]$city_clean
      dx <- as.numeric(dt$dx)
      dy <- as.numeric(dt$dy)
    }
  } else {
    dt <- prepared[[city_id]]$city_clean
    dx <- as.numeric(dt$dx)
    dy <- as.numeric(dt$dy)
  }

  norm <- sqrt(dx^2 + dy^2)
  norm[norm == 0] <- NA_real_

  cbind(ux = dx / norm, uy = dy / norm)
})
names(trip_dirs_full) <- names(prepared)

# sanity check: must match labels length
for (city_id in names(prepared)) {
  stopifnot(nrow(prepared[[city_id]]$X_scaled) == nrow(trip_dirs_full[[city_id]]))
}

############################################################
# Step 11a — Size-weighted within-cluster directional cohesion
############################################################
dir_cohesion_within_weighted <- function(labels, u,
                                         noise_label = 0L,
                                         min_cluster_size = 2L,
                                         map_to_01 = TRUE) {

  labels <- as.integer(labels)
  ux <- as.numeric(u[, "ux"])
  uy <- as.numeric(u[, "uy"])

  ok <- !is.na(labels) & labels != noise_label & !is.na(ux) & !is.na(uy)
  if (!any(ok)) return(0)

  labels <- labels[ok]
  ux <- ux[ok]
  uy <- uy[ok]

  cl_ids <- sort(unique(labels))
  if (length(cl_ids) == 0) return(0)

  num <- 0
  den <- 0

  for (cid in cl_ids) {
    idx <- which(labels == cid)
    n <- length(idx)
    if (n < min_cluster_size) next

    mx <- mean(ux[idx]); my <- mean(uy[idx])
    nm <- sqrt(mx^2 + my^2)
    if (!is.finite(nm) || nm < 1e-12) next
    mx <- mx / nm; my <- my / nm

    cos_i <- ux[idx] * mx + uy[idx] * my
    cos_i <- pmax(-1, pmin(1, cos_i))  # clamp

    coh_c <- mean(cos_i)

    num <- num + n * coh_c
    den <- den + n
  }

  if (den == 0) return(0)

  coh <- num / den  # [-1, 1]
  if (map_to_01) (coh + 1) / 2 else coh
}

############################################################
# Step 12 — Two objectives (parameter-configuration evaluation)
# Objective 1: maximize coverage (clustered trips share)
# Objective 2: maximize within-cluster directional cohesion (size-weighted)
#
# NOTE (NSGA-II minimization form):
# We return f1 = -coverage and f2 = -dir_cohesion.
# In this objective space, "worst" is (0,0) and "best" is (-1,-1).
############################################################

eval_solution <- function(labels, u_trips,
                          noise_label = 0L,
                          min_cluster_size = 2L) {

  labels <- as.integer(labels)

  if (!is.matrix(u_trips) && !is.data.frame(u_trips)) {
    stop("eval_solution: u_trips must be a matrix/data.frame with columns ux, uy")
  }
  if (!all(c("ux", "uy") %in% colnames(u_trips))) {
    stop("eval_solution: u_trips must have columns named 'ux' and 'uy'")
  }

  N_total <- length(labels)
  if (N_total <= 0) stop("eval_solution: empty labels")

  if (nrow(u_trips) != N_total) {
    stop("eval_solution: length(labels) must match nrow(u_trips)")
  }

  clustered <- !is.na(labels) & labels != noise_label
  n_clustered <- sum(clustered)
  coverage <- n_clustered / N_total

  K <- if (n_clustered == 0) 0L else length(unique(labels[clustered]))

  coh <- dir_cohesion_within_weighted(
    labels = labels,
    u = u_trips,
    noise_label = noise_label,
    min_cluster_size = min_cluster_size,
    map_to_01 = TRUE          # ensures coh in [0,1]
  )

  # NSGA-II minimizes => negate objectives
  f1 <- if (coverage == 0) 0 else -coverage   # in [-1, 0], worst = 0
  f2 <- if (coverage == 0) 0 else -coh        # in [-1, 0], worst = 0

  c(
    f1_neg_coverage = f1,
    f2_neg_dir_cohesion = f2,
    coverage = coverage,
    dir_cohesion = coh,
    K = K,
    clustered_n = n_clustered,
    N_total = N_total
  )
}

############################################################
# Step 12a — Sanity checks (per city) using FULL directions
############################################################

set.seed(1)

for (city_id in names(trip_dirs_full)) {

  u_trips <- trip_dirs_full[[city_id]]
  N <- nrow(u_trips)

  cat("\nCITY:", toupper(city_id), "\n")
  cat("N trips =", N, "\n")

  labels_all_noise <- rep(0L, N)
  print(eval_solution(labels_all_noise, u_trips, noise_label = 0L))

  labels_one_cluster <- rep(1L, N)
  print(eval_solution(labels_one_cluster, u_trips, noise_label = 0L))

  labels_toy <- rep(0L, N)
  clustered_idx <- sample.int(N, size = floor(0.7 * N))
  labels_toy[clustered_idx] <- sample.int(20, length(clustered_idx), replace = TRUE)
  print(eval_solution(labels_toy, u_trips, noise_label = 0L))
}

############################################################
# Step 12b — Random sanity checks (per city) using FULL directions
############################################################

set.seed(1)

for (city_id in names(trip_dirs_full)) {

  u_trips <- trip_dirs_full[[city_id]]
  N <- nrow(u_trips)

  cat("\nCITY:", toupper(city_id), "\n")

  for (i in 1:5) {

    clustered_share <- runif(1, 0.3, 0.9)
    K_toy <- sample(c(5, 10, 20, 40), 1)

    labels <- rep(0L, N)
    clustered_idx <- sample.int(N, size = floor(clustered_share * N))
    labels[clustered_idx] <- sample.int(K_toy, length(clustered_idx), replace = TRUE)

    res <- eval_solution(labels, u_trips, noise_label = 0L)

    cat(
      "draw", i,
      "| clustered_share ~", round(clustered_share, 2),
      "| K_toy =", K_toy,
      "| f1 =", round(res["f1_neg_coverage"], 6),
      "| f2 =", round(res["f2_neg_dir_cohesion"], 6),
      "| coverage =", round(res["coverage"], 4),
      "| cohesion =", round(res["dir_cohesion"], 4),
      "\n"
    )
  }
}

############################################################
# Step 14 — NSGA-II for DBSCAN (ALIGNED + HARDENED)
# - Minimizes: (-coverage, -dir_cohesion)
# - Uses shared snap helpers in BOTH objective wrapper & Pareto table
# - Reproducible: seed for subsample + per-city nsga2 seed
############################################################

suppressPackageStartupMessages({
  library(dbscan)
  library(dplyr)
  library(ggplot2)
  library(mco)
})

# -----------------------------
# Helpers (single source of truth)
# -----------------------------
snap_to_allowed <- function(v, allowed) allowed[which.min(abs(allowed - v))]

snap_eps <- function(eps_raw, eps_min, eps_max, eps_step, eps_digits = 3) {
  eps_raw <- max(eps_min, min(eps_raw, eps_max))
  eps <- eps_min + round((eps_raw - eps_min) / eps_step) * eps_step
  eps <- round(eps, eps_digits)
  eps <- max(eps_min, min(eps, eps_max))
  eps
}

# ==========================================================
# Objective function for DBSCAN
# ==========================================================
eval_dbscan_params <- function(eps, minPts, X_scaled, u_trips,
                               noise_label = 0L, min_cluster_size = 2L) {

  fit <- dbscan::dbscan(X_scaled, eps = eps, minPts = minPts)
  cl <- as.integer(fit$cluster)  # 0 = noise

  N_total <- length(cl)
  if (N_total == 0) return(c(f1 = 1, f2 = 1))

  clustered <- (cl != noise_label) & !is.na(cl)
  n_clustered <- sum(clustered)
  coverage <- n_clustered / N_total

  if (n_clustered == 0) return(c(f1 = 1, f2 = 1))

  coh <- dir_cohesion_within_weighted(
    labels = cl,
    u = u_trips,
    noise_label = noise_label,
    min_cluster_size = min_cluster_size,
    map_to_01 = TRUE
  )

  c(f1 = -coverage, f2 = -coh)
}

# ==========================================================
# Cached wrapper (snaps eps + minPts; caches evaluations)
# ==========================================================
make_objfun_city_cached_dbscan <- function(city_id, data_for_opt, trip_dirs_for_opt,
                                           eps_min, eps_max, eps_step,
                                           minPts_allowed,
                                           eps_digits = 3,
                                           noise_label = 0L,
                                           min_cluster_size = 2L) {

  X_scaled <- data_for_opt[[city_id]]$X_scaled
  u_trips  <- trip_dirs_for_opt[[city_id]]
  stopifnot(nrow(X_scaled) == nrow(u_trips))

  cache <- new.env(parent = emptyenv())

  function(x) {
    eps <- snap_eps(x[1], eps_min, eps_max, eps_step, eps_digits)
    mp  <- as.integer(round(x[2]))
    minPts <- snap_to_allowed(mp, minPts_allowed)

    key <- paste0(sprintf(paste0("%.", eps_digits, "f"), eps), "_", minPts)
    if (exists(key, envir = cache, inherits = FALSE)) {
      return(get(key, envir = cache, inherits = FALSE))
    }

    val <- eval_dbscan_params(
      eps = eps,
      minPts = minPts,
      X_scaled = X_scaled,
      u_trips = u_trips,
      noise_label = noise_label,
      min_cluster_size = min_cluster_size
    )

    assign(key, val, envir = cache)
    val
  }
}

# ==========================================================
# Step 14a — Subsample for optimization (shared indices)
# ==========================================================
set.seed(1)
N_SUB <- 10000

sub_idx <- lapply(names(prepared), function(city_id) {
  n <- nrow(prepared[[city_id]]$X_scaled)
  if (n > N_SUB) sample.int(n, N_SUB) else seq_len(n)
})
names(sub_idx) <- names(prepared)

data_for_opt <- lapply(names(prepared), function(city_id) {
  idx <- sub_idx[[city_id]]
  list(X_scaled = prepared[[city_id]]$X_scaled[idx, , drop = FALSE])
})
names(data_for_opt) <- names(prepared)

trip_dirs_for_opt <- lapply(names(prepared), function(city_id) {
  idx <- sub_idx[[city_id]]
  trip_dirs_full[[city_id]][idx, , drop = FALSE]
})
names(trip_dirs_for_opt) <- names(prepared)

for (city_id in names(prepared)) {
  stopifnot(nrow(data_for_opt[[city_id]]$X_scaled) == nrow(trip_dirs_for_opt[[city_id]]))
}

cat("\n[DBSCAN] Subsample sizes:\n")
for (city_id in names(prepared)) {
  cat(city_id, "FULL =", nrow(prepared[[city_id]]$X_scaled),
      "| SUB =", nrow(data_for_opt[[city_id]]$X_scaled), "\n")
}

# ==========================================================
# Step 14b — NSGA-II run (per city) [ALIGNED BUDGET]
# ==========================================================
eps_min  <- 0.08
eps_max  <- 0.45
eps_step <- 0.02
eps_digits <- 3

minPts_allowed <- c(6L, 8L, 10L, 12L, 15L, 20L, 25L, 30L, 40L, 50L, 60L, 70L)

popsize <- 100
generations <- 50

objfun_dbscan <- list()
nsga_dbscan <- list()

for (i in seq_along(names(prepared))) {
  city_id <- names(prepared)[i]
  cat("\n[DBSCAN] Starting NSGA-II for:", toupper(city_id), "\n")
  t0 <- Sys.time()

  objfun_dbscan[[city_id]] <- make_objfun_city_cached_dbscan(
    city_id = city_id,
    data_for_opt = data_for_opt,
    trip_dirs_for_opt = trip_dirs_for_opt,
    eps_min = eps_min, eps_max = eps_max, eps_step = eps_step,
    minPts_allowed = minPts_allowed,
    eps_digits = eps_digits,
    noise_label = 0L,
    min_cluster_size = 2L
  )

  set.seed(1000 + i)  # reproducible per city
  res <- mco::nsga2(
    fn = objfun_dbscan[[city_id]],
    idim = 2,
    odim = 2,
    lower.bounds = c(eps_min, min(minPts_allowed)),
    upper.bounds = c(eps_max, max(minPts_allowed)),
    popsize = popsize,
    generations = generations
  )

  nsga_dbscan[[city_id]] <- res

  cat("[DBSCAN] Finished:", toupper(city_id),
      "| minutes =", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2), "\n")
}

# ==========================================================
# Step 14c — Pareto tables + knee selection (ALIGNED)
# ==========================================================
pareto_dbscan <- list()

for (city_id in names(nsga_dbscan)) {

  res <- nsga_dbscan[[city_id]]

  eps_raw   <- res$par[, 1]
  minPts_raw <- as.integer(round(res$par[, 2]))

  eps_snap <- vapply(eps_raw, snap_eps, numeric(1),
                     eps_min = eps_min, eps_max = eps_max,
                     eps_step = eps_step, eps_digits = eps_digits)

  minPts_snap <- vapply(minPts_raw, snap_to_allowed, integer(1), allowed = minPts_allowed)

  df <- data.frame(
    city = city_id,
    eps = as.numeric(eps_snap),
    minPts = as.integer(minPts_snap),
    coverage = -res$value[, 1],
    dir_cohesion = -res$value[, 2]
  )

  df <- df %>%
    filter(is.finite(coverage), is.finite(dir_cohesion),
           coverage >= 0, dir_cohesion >= 0) %>%
    distinct(eps, minPts, .keep_all = TRUE)

  pareto_dbscan[[city_id]] <- df

  cat("[DBSCAN] CITY:", toupper(city_id), "| Pareto candidates:", nrow(df), "\n")
  print(head(df, 8))
}

knee_dbscan <- list()

for (city_id in names(pareto_dbscan)) {

  df <- pareto_dbscan[[city_id]]
  if (is.null(df) || nrow(df) == 0) {
    knee_dbscan[[city_id]] <- NULL
    next
  }

  cov_n <- (df$coverage - min(df$coverage)) / (max(df$coverage) - min(df$coverage) + 1e-9)
  coh_n <- (df$dir_cohesion - min(df$dir_cohesion)) / (max(df$dir_cohesion) - min(df$dir_cohesion) + 1e-9)
  d <- sqrt((1 - cov_n)^2 + (1 - coh_n)^2)

  knee_dbscan[[city_id]] <- df[which.min(d), , drop = FALSE]

  cat("\n[DBSCAN] KNEE:", toupper(city_id), "\n")
  print(knee_dbscan[[city_id]])
}

plot_pareto <- function(df, knee = NULL, expert = NULL, city_name, title_prefix = "Pareto front (subsample)") {
  p <- ggplot(df, aes(x = coverage, y = dir_cohesion)) +
    geom_point(size = 2, alpha = 0.8) +
    labs(
      title = paste(title_prefix, "—", toupper(city_name)),
      x = "Coverage (clustered share)",
      y = "Directional cohesion (within-cluster, size-weighted)"
    ) +
    theme_minimal()

  if (!is.null(knee) && is.data.frame(knee) && nrow(knee) > 0) {
    p <- p + geom_point(data = knee, aes(x = coverage, y = dir_cohesion),
                        color = "red", size = 4)
  }
  if (!is.null(expert) && is.data.frame(expert) && nrow(expert) > 0) {
    p <- p + geom_point(data = expert, aes(x = coverage, y = dir_cohesion),
                        shape = 4, color = "black", size = 4, stroke = 1.2)
  }
  p
}

for (city_id in names(pareto_dbscan)) {
  df <- pareto_dbscan[[city_id]]
  if (!is.null(df) && nrow(df) > 0) {
    print(plot_pareto(df, knee = knee_dbscan[[city_id]], city_name = city_id, title_prefix = "[DBSCAN] Pareto front (subsample)"))
  }
}

############################################################
# Step 14 — NSGA-II for SNN (ALIGNED to DBSCAN + to your Step 7)
# - Minimizes: (-coverage, -dir_cohesion)
# - SAME subsample indices/budget/filtering/knee logic as DBSCAN
# - Caches sNN object per k for speed (important)
############################################################

suppressPackageStartupMessages({
  library(dbscan)
  library(dplyr)
  library(ggplot2)
  library(mco)
})

# -----------------------------
# Helpers (same as DBSCAN)
# -----------------------------
snap_to_allowed <- function(v, allowed) allowed[which.min(abs(allowed - v))]

# ==========================================================
# Objective function for SNN (matches your Step 7)
# ==========================================================
eval_snn_params <- function(k, snn_threshold, minPts,
                            X_scaled, u_trips,
                            noise_label = 0L, min_cluster_size = 2L,
                            jp = FALSE, sort = FALSE, borderPoints = TRUE) {

  snn_obj <- dbscan::sNN(X_scaled, k = k, sort = sort, jp = jp)

  fit <- dbscan::sNNclust(
    x = snn_obj,
    k = k,
    eps = snn_threshold,
    minPts = minPts,
    borderPoints = borderPoints
  )

  cl <- as.integer(fit$cluster)  # 0 = noise

  N_total <- length(cl)
  if (N_total == 0) return(c(f1 = 1, f2 = 1))

  clustered <- (cl != noise_label) & !is.na(cl)
  n_clustered <- sum(clustered)
  coverage <- n_clustered / N_total

  if (n_clustered == 0) return(c(f1 = 1, f2 = 1))

  coh <- dir_cohesion_within_weighted(
    labels = cl,
    u = u_trips,
    noise_label = noise_label,
    min_cluster_size = min_cluster_size,
    map_to_01 = TRUE
  )

  c(f1 = -coverage, f2 = -coh)
}

# ==========================================================
# Cached wrapper for SNN
# - snaps (k, threshold, minPts) to allowed sets
# - enforces constraints: threshold <= k-1, minPts <= k
# - caches sNN object per k (big speed gain)
# ==========================================================
make_objfun_city_cached_snn <- function(city_id, data_for_opt, trip_dirs_for_opt,
                                        k_allowed, thr_allowed, minPts_allowed,
                                        noise_label = 0L, min_cluster_size = 2L,
                                        jp = FALSE, sort = FALSE, borderPoints = TRUE) {

  X_scaled <- data_for_opt[[city_id]]$X_scaled
  u_trips  <- trip_dirs_for_opt[[city_id]]
  stopifnot(nrow(X_scaled) == nrow(u_trips))

  cache_val <- new.env(parent = emptyenv())  # caches objective values
  cache_snn <- new.env(parent = emptyenv())  # caches snn_obj by k

  function(x) {
    k_raw   <- as.integer(round(x[1]))
    th_raw  <- as.integer(round(x[2]))
    mp_raw  <- as.integer(round(x[3]))

    k  <- snap_to_allowed(k_raw,  k_allowed)
    th <- snap_to_allowed(th_raw, thr_allowed)
    mp <- snap_to_allowed(mp_raw, minPts_allowed)

    # enforce SNN constraints (consistent everywhere)
    th <- max(1L, min(th, k - 1L))
    mp <- max(2L, min(mp, k))

    key <- paste(k, th, mp, sep = "_")
    if (exists(key, envir = cache_val, inherits = FALSE)) {
      return(get(key, envir = cache_val, inherits = FALSE))
    }

    # cache sNN object per k for speed
    k_key <- as.character(k)
    if (exists(k_key, envir = cache_snn, inherits = FALSE)) {
      snn_obj <- get(k_key, envir = cache_snn, inherits = FALSE)
    } else {
      snn_obj <- dbscan::sNN(X_scaled, k = k, sort = sort, jp = jp)
      assign(k_key, snn_obj, envir = cache_snn)
    }

    fit <- dbscan::sNNclust(
      x = snn_obj,
      k = k,
      eps = th,
      minPts = mp,
      borderPoints = borderPoints
    )

    cl <- as.integer(fit$cluster)

    N_total <- length(cl)
    if (N_total == 0) {
      val <- c(f1 = 1, f2 = 1)
    } else {
      clustered <- (cl != noise_label) & !is.na(cl)
      n_clustered <- sum(clustered)
      coverage <- n_clustered / N_total

      if (n_clustered == 0) {
        val <- c(f1 = 1, f2 = 1)
      } else {
        coh <- dir_cohesion_within_weighted(
          labels = cl,
          u = u_trips,
          noise_label = noise_label,
          min_cluster_size = min_cluster_size,
          map_to_01 = TRUE
        )
        val <- c(f1 = -coverage, f2 = -coh)
      }
    }

    assign(key, val, envir = cache_val)
    val
  }
}

# ==========================================================
# Step 14a — Subsample for optimization (REUSE DBSCAN subsample!)
# IMPORTANT: Do NOT resample here. Use the SAME:
#   sub_idx, data_for_opt, trip_dirs_for_opt
# created in the DBSCAN block above.
# ==========================================================

# ==========================================================
# Step 14b — NSGA-II run (per city) [ALIGNED BUDGET]  ✅ FIXED BOUNDS TYPE
# ==========================================================
k_allowed      <- c(15L, 20L, 25L, 30L, 35L)
thr_allowed    <- 4L:14L
minPts_allowed <- c(6L, 8L, 10L, 12L, 15L, 18L, 20L)

popsize <- 100
generations <- 50

objfun_snn <- list()
nsga_snn <- list()

lb <- as.numeric(c(min(k_allowed), min(thr_allowed), min(minPts_allowed)))
ub <- as.numeric(c(max(k_allowed), max(thr_allowed), max(minPts_allowed)))

for (i in seq_along(names(prepared))) {
  city_id <- names(prepared)[i]
  cat("\n[SNN] Starting NSGA-II for:", toupper(city_id), "\n")
  t0 <- Sys.time()

  objfun_snn[[city_id]] <- make_objfun_city_cached_snn(
    city_id = city_id,
    data_for_opt = data_for_opt,
    trip_dirs_for_opt = trip_dirs_for_opt,
    k_allowed = k_allowed,
    thr_allowed = thr_allowed,
    minPts_allowed = minPts_allowed,
    noise_label = 0L,
    min_cluster_size = 2L,
    jp = FALSE, sort = FALSE, borderPoints = TRUE
  )

  set.seed(2000 + i)
  res <- mco::nsga2(
    fn = objfun_snn[[city_id]],
    idim = 3,
    odim = 2,
    lower.bounds = lb,
    upper.bounds = ub,
    popsize = popsize,
    generations = generations
  )

  nsga_snn[[city_id]] <- res

  cat("[SNN] Finished:", toupper(city_id),
      "| minutes =", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2), "\n")
}

# ==========================================================
# Step 14c — Pareto tables + knee selection (ALIGNED)
# ==========================================================
pareto_snn <- list()

for (city_id in names(nsga_snn)) {

  res <- nsga_snn[[city_id]]

  k_raw  <- as.integer(round(res$par[, 1]))
  th_raw <- as.integer(round(res$par[, 2]))
  mp_raw <- as.integer(round(res$par[, 3]))

  k  <- vapply(k_raw,  snap_to_allowed, integer(1), allowed = k_allowed)
  th <- vapply(th_raw, snap_to_allowed, integer(1), allowed = thr_allowed)
  mp <- vapply(mp_raw, snap_to_allowed, integer(1), allowed = minPts_allowed)

  # enforce constraints same as objective wrapper
  th <- pmax(1L, pmin(th, k - 1L))
  mp <- pmax(2L, pmin(mp, k))

  df <- data.frame(
    city = city_id,
    k = k,
    snn_threshold = th,
    minPts = mp,
    coverage = -res$value[, 1],
    dir_cohesion = -res$value[, 2]
  )

  df <- df %>%
    filter(is.finite(coverage), is.finite(dir_cohesion),
           coverage >= 0, dir_cohesion >= 0) %>%
    distinct(k, snn_threshold, minPts, .keep_all = TRUE)

  pareto_snn[[city_id]] <- df

  cat("[SNN] CITY:", toupper(city_id), "| Pareto candidates:", nrow(df), "\n")
  print(head(df, 8))
}

knee_snn <- list()

for (city_id in names(pareto_snn)) {

  df <- pareto_snn[[city_id]]
  if (is.null(df) || nrow(df) == 0) {
    knee_snn[[city_id]] <- NULL
    next
  }

  cov_n <- (df$coverage - min(df$coverage)) / (max(df$coverage) - min(df$coverage) + 1e-9)
  coh_n <- (df$dir_cohesion - min(df$dir_cohesion)) / (max(df$dir_cohesion) - min(df$dir_cohesion) + 1e-9)
  d <- sqrt((1 - cov_n)^2 + (1 - coh_n)^2)

  knee_snn[[city_id]] <- df[which.min(d), , drop = FALSE]

  cat("\n[SNN] KNEE:", toupper(city_id), "\n")
  print(knee_snn[[city_id]])
}

for (city_id in names(pareto_snn)) {
  df <- pareto_snn[[city_id]]
  if (!is.null(df) && nrow(df) > 0) {
    print(plot_pareto(df, knee = knee_snn[[city_id]], city_name = city_id, title_prefix = "[SNN] Pareto front (subsample)"))
  }
}

# ==========================================================
# Print ALL DBSCAN Pareto points (per city)
# ==========================================================
options(max.print = 1e6)  # prevent truncation

for (city_id in names(pareto_dbscan)) {

  cat("\n========================================\n")
  cat("[DBSCAN] CITY:", toupper(city_id), "\n")
  cat("N Pareto points:", nrow(pareto_dbscan[[city_id]]), "\n")
  cat("========================================\n\n")

  df <- pareto_dbscan[[city_id]]

  if (is.null(df) || nrow(df) == 0) {
    cat("No valid Pareto rows.\n")
    next
  }

  # Optional: sort by coverage for readability
  df <- df[order(df$coverage), ]

  print(df, row.names = FALSE)
}

# ==========================================================
# Print ALL SNN Pareto points (per city)
# ==========================================================
options(max.print = 1e6)

for (city_id in names(pareto_snn)) {

  cat("\n========================================\n")
  cat("[SNN] CITY:", toupper(city_id), "\n")
  cat("N Pareto points:", nrow(pareto_snn[[city_id]]), "\n")
  cat("========================================\n\n")

  df <- pareto_snn[[city_id]]

  if (is.null(df) || nrow(df) == 0) {
    cat("No valid Pareto rows.\n")
    next
  }

  # Optional: sort by coverage for readability
  df <- df[order(df$coverage), ]

  print(df, row.names = FALSE)
}

############################################################
# Step 15 — Validate selected DBSCAN Pareto configs on FULL data
# Selection strategy (per city):
#   - knee (from subsample Pareto)
#   - 2 minimally higher-coverage points
#   - 2 minimally higher-cohesion points
#
# Uses FULL data objects:
#   prepared_full <- prepared
############################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(dbscan)
})

# ---- Freeze FULL data object
prepared_full <- prepared

# ---- Quick sanity: print full sizes
for (city_id in names(prepared_full)) {
  cat("[DBSCAN] FULL:", toupper(city_id), "N =", nrow(prepared_full[[city_id]]$X_scaled), "\n")
}
if (exists("data_for_opt")) {
  for (city_id in names(data_for_opt)) {
    cat("[DBSCAN] SUB :", toupper(city_id), "N =", nrow(data_for_opt[[city_id]]$X_scaled), "\n")
  }
}

# ----------------------------------------------------------
# Clean + deduplicate Pareto table (from subsample NSGA-II)
# ----------------------------------------------------------
pareto_unique_dbscan <- lapply(names(pareto_dbscan), function(city_id) {
  df <- pareto_dbscan[[city_id]]
  if (is.null(df) || nrow(df) == 0) return(df)

  df %>%
    select(city, eps, minPts, coverage, dir_cohesion) %>%
    mutate(
      eps = as.numeric(eps),
      minPts = as.integer(minPts),
      coverage = as.numeric(coverage),
      dir_cohesion = as.numeric(dir_cohesion)
    ) %>%
    filter(is.finite(coverage), is.finite(dir_cohesion),
           coverage >= 0, dir_cohesion >= 0) %>%
    distinct(eps, minPts, .keep_all = TRUE)
})
names(pareto_unique_dbscan) <- names(pareto_dbscan)

# ----------------------------------------------------------
# Select ~5 configs around knee
# ----------------------------------------------------------
select_configs_dbscan <- function(df, knee_row, n_each = 2) {

  if (is.null(df) || nrow(df) == 0) return(df)

  # fallback if knee is missing
  if (is.null(knee_row) || !is.data.frame(knee_row) || nrow(knee_row) == 0) {
    return(df[order(-df$coverage, -df$dir_cohesion), ][1:min(5, nrow(df)), , drop = FALSE])
  }

  # match knee by parameters (tolerant float match on eps)
  knee_match <- df %>%
    filter(abs(eps - knee_row$eps[1]) < 1e-6, minPts == knee_row$minPts[1])

  if (nrow(knee_match) == 0) {
    # fallback: match by objective values
    knee_match <- df %>%
      filter(abs(coverage - knee_row$coverage[1]) < 1e-9,
             abs(dir_cohesion - knee_row$dir_cohesion[1]) < 1e-9)
  }

  if (nrow(knee_match) == 0) {
    knee_match <- df[order(-df$coverage, -df$dir_cohesion), ][1, , drop = FALSE]
  } else {
    knee_match <- knee_match[1, , drop = FALSE]
  }

  higher_cov <- df %>%
    filter(coverage > knee_match$coverage) %>%
    arrange(coverage, desc(dir_cohesion)) %>%
    head(n_each)

  higher_coh <- df %>%
    filter(dir_cohesion > knee_match$dir_cohesion) %>%
    arrange(dir_cohesion, desc(coverage)) %>%
    head(n_each)

  bind_rows(knee_match, higher_cov, higher_coh) %>%
    distinct(eps, minPts, .keep_all = TRUE) %>%
    head(1 + 2*n_each)
}

# ----------------------------------------------------------
# FULL-data evaluation of one DBSCAN config
# ----------------------------------------------------------
full_eval_one_dbscan <- function(city_id, eps, minPts, prepared_full,
                                 noise_label = 0L, min_cluster_size = 2L) {

  X_full <- prepared_full[[city_id]]$X_scaled
  city_clean <- prepared_full[[city_id]]$city_clean
  stopifnot(nrow(X_full) == nrow(city_clean))

  # unit directions (aligned to X_full rows)
  dx <- as.numeric(city_clean$dx)
  dy <- as.numeric(city_clean$dy)
  norm <- sqrt(dx^2 + dy^2)
  norm[norm == 0] <- NA_real_
  u_full <- cbind(ux = dx / norm, uy = dy / norm)

  fit <- dbscan::dbscan(X_full, eps = eps, minPts = minPts)
  cl <- as.integer(fit$cluster)

  N_total <- length(cl)
  clustered <- (cl != noise_label) & !is.na(cl)
  n_clustered <- sum(clustered)

  coverage <- if (N_total > 0) n_clustered / N_total else 0
  K <- if (n_clustered == 0) 0L else length(unique(cl[clustered]))

  coh <- if (n_clustered == 0) 0 else dir_cohesion_within_weighted(
    labels = cl,
    u = u_full,
    noise_label = noise_label,
    min_cluster_size = min_cluster_size,
    map_to_01 = TRUE
  )

  data.frame(
    city = city_id,
    eps = as.numeric(eps),
    minPts = as.integer(minPts),
    coverage = as.numeric(coverage),
    clustered_pct = 100 * coverage,
    dir_cohesion = as.numeric(coh),
    n_clusters = as.integer(K),
    noise_pct = if (N_total > 0) 100 * sum(cl == noise_label, na.rm = TRUE) / N_total else NA_real_,
    n_total = as.integer(N_total)
  )
}

# ----------------------------------------------------------
# Run FULL validation for selected Pareto configs + expert baseline
# ----------------------------------------------------------
full_validated_dbscan <- list()

for (city_id in names(pareto_unique_dbscan)) {

  df <- pareto_unique_dbscan[[city_id]]
  knee <- knee_dbscan[[city_id]]  # <-- updated name

  if (is.null(df) || nrow(df) == 0) {
    cat("[DBSCAN] Skipping", toupper(city_id), "(no Pareto rows).\n")
    next
  }

  chosen <- select_configs_dbscan(df, knee_row = knee, n_each = 2)

  cat("\n[DBSCAN] Validating", nrow(chosen), "Pareto configs on FULL data for", toupper(city_id), "\n")

  pareto_out <- do.call(rbind, lapply(seq_len(nrow(chosen)), function(i) {
    full_eval_one_dbscan(city_id, chosen$eps[i], chosen$minPts[i], prepared_full)
  }))
  pareto_out$config <- "pareto_local"

  # ---- Expert baseline (must exist as dbscan_params list)
  if (exists("dbscan_params") && !is.null(dbscan_params[[city_id]])) {
    ex <- dbscan_params[[city_id]]
    expert_out <- full_eval_one_dbscan(city_id, ex$eps, ex$minPts, prepared_full)
    expert_out$config <- "expert"
  } else {
    expert_out <- NULL
  }

  full_validated_dbscan[[city_id]] <- bind_rows(pareto_out, expert_out)
}

full_validated_dbscan_df <- if (length(full_validated_dbscan) > 0) do.call(rbind, full_validated_dbscan) else data.frame()
print(full_validated_dbscan_df)

# ----------------------------------------------------------
# Re-select knee on FULL data (among validated Pareto configs only)
# ----------------------------------------------------------
knee_full_dbscan <- NULL
if (nrow(full_validated_dbscan_df) > 0) {

  knee_full_dbscan <- full_validated_dbscan_df %>%
    filter(config == "pareto_local") %>%
    group_by(city) %>%
    mutate(
      cov_n = (coverage - min(coverage)) / (max(coverage) - min(coverage) + 1e-9),
      coh_n = (dir_cohesion - min(dir_cohesion)) / (max(dir_cohesion) - min(dir_cohesion) + 1e-9),
      dist_ideal = sqrt((1 - cov_n)^2 + (1 - coh_n)^2)
    ) %>%
    slice_min(dist_ideal, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(config = "knee_full") %>%
    select(city, config, eps, minPts, coverage, clustered_pct, dir_cohesion, n_clusters, noise_pct)

  cat("\n[DBSCAN] KNEE (FULL data, among validated Pareto configs):\n")
  print(knee_full_dbscan)

  compare_df <- full_validated_dbscan_df %>%
    filter(config %in% c("expert")) %>%
    select(city, config, eps, minPts, coverage, clustered_pct, dir_cohesion, n_clusters, noise_pct) %>%
    bind_rows(knee_full_dbscan) %>%
    arrange(city, factor(config, levels = c("expert", "knee_full")))

  cat("\n[DBSCAN] EXPERT vs KNEE_FULL (FULL data):\n")
  print(compare_df)
}

############################################################
# Step 15 — Validate selected SNN Pareto configs on FULL data
# Same selection strategy as DBSCAN
############################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(dbscan)
})

prepared_full <- prepared

# ---- Quick sanity: print full sizes
for (city_id in names(prepared_full)) {
  cat("[SNN] FULL:", toupper(city_id), "N =", nrow(prepared_full[[city_id]]$X_scaled), "\n")
}

# ----------------------------------------------------------
# Clean + deduplicate Pareto table (from subsample NSGA-II)
# ----------------------------------------------------------
pareto_unique_snn <- lapply(names(pareto_snn), function(city_id) {
  df <- pareto_snn[[city_id]]
  if (is.null(df) || nrow(df) == 0) return(df)

  df %>%
    select(city, k, snn_threshold, minPts, coverage, dir_cohesion) %>%
    mutate(
      k = as.integer(k),
      snn_threshold = as.integer(snn_threshold),
      minPts = as.integer(minPts),
      coverage = as.numeric(coverage),
      dir_cohesion = as.numeric(dir_cohesion)
    ) %>%
    filter(is.finite(coverage), is.finite(dir_cohesion),
           coverage >= 0, dir_cohesion >= 0) %>%
    distinct(k, snn_threshold, minPts, .keep_all = TRUE)
})
names(pareto_unique_snn) <- names(pareto_snn)

# ----------------------------------------------------------
# Select ~5 configs around knee
# ----------------------------------------------------------
select_configs_snn <- function(df, knee_row, n_each = 2) {

  if (is.null(df) || nrow(df) == 0) return(df)

  if (is.null(knee_row) || !is.data.frame(knee_row) || nrow(knee_row) == 0) {
    return(df[order(-df$coverage, -df$dir_cohesion), ][1:min(5, nrow(df)), , drop = FALSE])
  }

  # match knee by parameters
  knee_match <- df %>%
    filter(k == knee_row$k[1],
           snn_threshold == knee_row$snn_threshold[1],
           minPts == knee_row$minPts[1])

  if (nrow(knee_match) == 0) {
    knee_match <- df %>%
      filter(abs(coverage - knee_row$coverage[1]) < 1e-9,
             abs(dir_cohesion - knee_row$dir_cohesion[1]) < 1e-9)
  }

  if (nrow(knee_match) == 0) {
    knee_match <- df[order(-df$coverage, -df$dir_cohesion), ][1, , drop = FALSE]
  } else {
    knee_match <- knee_match[1, , drop = FALSE]
  }

  higher_cov <- df %>%
    filter(coverage > knee_match$coverage) %>%
    arrange(coverage, desc(dir_cohesion)) %>%
    head(n_each)

  higher_coh <- df %>%
    filter(dir_cohesion > knee_match$dir_cohesion) %>%
    arrange(dir_cohesion, desc(coverage)) %>%
    head(n_each)

  bind_rows(knee_match, higher_cov, higher_coh) %>%
    distinct(k, snn_threshold, minPts, .keep_all = TRUE) %>%
    head(1 + 2*n_each)
}

# ----------------------------------------------------------
# FULL-data evaluation of one SNN config
# (recomputes u_full from FULL city_clean)
# ----------------------------------------------------------
full_eval_one_snn <- function(city_id, k, snn_threshold, minPts, prepared_full,
                              noise_label = 0L, min_cluster_size = 2L,
                              jp = FALSE, sort = FALSE, borderPoints = TRUE) {

  X_full <- prepared_full[[city_id]]$X_scaled
  city_clean <- prepared_full[[city_id]]$city_clean
  stopifnot(nrow(X_full) == nrow(city_clean))

  # unit directions
  dx <- as.numeric(city_clean$dx)
  dy <- as.numeric(city_clean$dy)
  norm <- sqrt(dx^2 + dy^2)
  norm[norm == 0] <- NA_real_
  u_full <- cbind(ux = dx / norm, uy = dy / norm)

  # build SNN object then cluster (matches your Step 7)
  snn_obj <- dbscan::sNN(X_full, k = k, sort = sort, jp = jp)
  fit <- dbscan::sNNclust(
    x = snn_obj,
    k = k,
    eps = snn_threshold,
    minPts = minPts,
    borderPoints = borderPoints
  )

  cl <- as.integer(fit$cluster)

  N_total <- length(cl)
  clustered <- (cl != noise_label) & !is.na(cl)
  n_clustered <- sum(clustered)

  coverage <- if (N_total > 0) n_clustered / N_total else 0
  K <- if (n_clustered == 0) 0L else length(unique(cl[clustered]))

  coh <- if (n_clustered == 0) 0 else dir_cohesion_within_weighted(
    labels = cl,
    u = u_full,
    noise_label = noise_label,
    min_cluster_size = min_cluster_size,
    map_to_01 = TRUE
  )

  data.frame(
    city = city_id,
    k = as.integer(k),
    snn_threshold = as.integer(snn_threshold),
    minPts = as.integer(minPts),
    coverage = as.numeric(coverage),
    clustered_pct = 100 * coverage,
    dir_cohesion = as.numeric(coh),
    n_clusters = as.integer(K),
    noise_pct = if (N_total > 0) 100 * sum(cl == noise_label, na.rm = TRUE) / N_total else NA_real_,
    n_total = as.integer(N_total)
  )
}

# ----------------------------------------------------------
# Run FULL validation for selected Pareto configs + expert baseline
# ----------------------------------------------------------
full_validated_snn <- list()

for (city_id in names(pareto_unique_snn)) {

  df <- pareto_unique_snn[[city_id]]
  knee <- knee_snn[[city_id]]   # <-- updated name

  if (is.null(df) || nrow(df) == 0) {
    cat("[SNN] Skipping", toupper(city_id), "(no Pareto rows).\n")
    next
  }

  chosen <- select_configs_snn(df, knee_row = knee, n_each = 2)

  cat("\n[SNN] Validating", nrow(chosen), "Pareto configs on FULL data for", toupper(city_id), "\n")

  pareto_out <- do.call(rbind, lapply(seq_len(nrow(chosen)), function(i) {
    full_eval_one_snn(city_id,
                      k = chosen$k[i],
                      snn_threshold = chosen$snn_threshold[i],
                      minPts = chosen$minPts[i],
                      prepared_full = prepared_full)
  }))
  pareto_out$config <- "pareto_local"

  # ---- Expert baseline (must exist as snn_params list from your Step 7)
  if (exists("snn_params") && !is.null(snn_params[[city_id]])) {
    ex <- snn_params[[city_id]]
    expert_out <- full_eval_one_snn(city_id,
                                    k = ex$k,
                                    snn_threshold = ex$snn_threshold,
                                    minPts = ex$minPts,
                                    prepared_full = prepared_full)
    expert_out$config <- "expert"
  } else {
    expert_out <- NULL
  }

  full_validated_snn[[city_id]] <- bind_rows(pareto_out, expert_out)
}

full_validated_snn_df <- if (length(full_validated_snn) > 0) do.call(rbind, full_validated_snn) else data.frame()
print(full_validated_snn_df)

# ----------------------------------------------------------
# Re-select knee on FULL data (among validated Pareto configs only)
# ----------------------------------------------------------
knee_full_snn <- NULL
if (nrow(full_validated_snn_df) > 0) {

  knee_full_snn <- full_validated_snn_df %>%
    filter(config == "pareto_local") %>%
    group_by(city) %>%
    mutate(
      cov_n = (coverage - min(coverage)) / (max(coverage) - min(coverage) + 1e-9),
      coh_n = (dir_cohesion - min(dir_cohesion)) / (max(dir_cohesion) - min(dir_cohesion) + 1e-9),
      dist_ideal = sqrt((1 - cov_n)^2 + (1 - coh_n)^2)
    ) %>%
    slice_min(dist_ideal, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(config = "knee_full") %>%
    select(city, config, k, snn_threshold, minPts, coverage, clustered_pct, dir_cohesion, n_clusters, noise_pct)

  cat("\n[SNN] KNEE (FULL data, among validated Pareto configs):\n")
  print(knee_full_snn)

  compare_df <- full_validated_snn_df %>%
    filter(config %in% c("expert")) %>%
    select(city, config, k, snn_threshold, minPts, coverage, clustered_pct, dir_cohesion, n_clusters, noise_pct) %>%
    bind_rows(knee_full_snn) %>%
    arrange(city, factor(config, levels = c("expert", "knee_full")))

  cat("\n[SNN] EXPERT vs KNEE_FULL (FULL data):\n")
  print(compare_df)
}

############################################################
# Step 16 — Hypervolume analysis for Pareto fronts (coverage + cohesion)
# - Works for BOTH algorithms: DBSCAN + SNN
# - Uses columns: coverage, dir_cohesion
# - Global normalization across ALL cities AND BOTH algorithms
# - Converts to MIN space: f1 = 1 - coverage_norm, f2 = 1 - cohesion_norm
############################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

# -----------------------------
# Non-dominated filter for MINIMIZATION
# -----------------------------
nondominated_min <- function(F) {
  n <- nrow(F)
  keep <- rep(TRUE, n)
  for (i in seq_len(n)) {
    if (!keep[i]) next
    for (j in seq_len(n)) {
      if (i == j || !keep[i]) next
      if (all(F[j, ] <= F[i, ]) && any(F[j, ] < F[i, ])) {
        keep[i] <- FALSE
      }
    }
  }
  F[keep, , drop = FALSE]
}

# -----------------------------
# Exact 2D dominated hypervolume for MINIMIZATION
# -----------------------------
hv2d_min <- function(F, ref) {
  if (nrow(F) == 0) return(0)

  F <- nondominated_min(F)
  o <- order(F[, 1], F[, 2])
  F <- F[o, , drop = FALSE]

  hv <- 0
  h <- ref[2]
  for (i in seq_len(nrow(F))) {
    f1 <- F[i, 1]
    f2 <- F[i, 2]
    if (f2 < h) {
      hv <- hv + (ref[1] - f1) * (h - f2)
      h <- f2
    }
  }
  hv
}

# Leave-one-out HV contributions
hv2d_contrib_min <- function(F, ref) {
  hv_all <- hv2d_min(F, ref)
  n <- nrow(F)
  sapply(seq_len(n), function(i) hv_all - hv2d_min(F[-i, , drop = FALSE], ref))
}

# Safe normalization
safe_norm <- function(z, zmin, zmax) {
  if (!is.finite(zmin) || !is.finite(zmax) || isTRUE(all.equal(zmax, zmin))) {
    return(rep(0, length(z)))
  }
  (z - zmin) / (zmax - zmin)
}

# Gini coefficient for contribution concentration
gini <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  x <- sort(pmax(0, x))
  s <- sum(x)
  if (s <= 0) return(0)
  n <- length(x)
  (2 * sum(seq_len(n) * x) / (n * s)) - (n + 1) / n
}

# -----------------------------
# Collect Pareto tables from BOTH algorithms
# -----------------------------
pareto_by_algo <- list(
  DBSCAN = pareto_dbscan,
  SNN    = pareto_snn
)

# -----------------------------
# Build global pool for normalization (ALL cities + BOTH algos)
# -----------------------------
all_pts <- do.call(rbind, lapply(names(pareto_by_algo), function(algo) {
  pt_list <- pareto_by_algo[[algo]]
  if (is.null(pt_list) || length(pt_list) == 0) return(NULL)

  do.call(rbind, lapply(names(pt_list), function(city_id) {
    df <- pt_list[[city_id]]
    if (is.null(df) || nrow(df) == 0) return(NULL)

    df2 <- df[, c("coverage", "dir_cohesion")]
    df2 <- df2[is.finite(df2$coverage) & is.finite(df2$dir_cohesion), , drop = FALSE]
    if (nrow(df2) == 0) return(NULL)

    df2$algorithm <- algo
    df2$city <- city_id
    df2
  }))
}))

if (is.null(all_pts) || nrow(all_pts) == 0) {
  hv_results <- data.frame()
  cat("No valid Pareto points available for HV analysis.\n")
} else {

  cov_min <- min(all_pts$coverage, na.rm = TRUE)
  cov_max <- max(all_pts$coverage, na.rm = TRUE)
  coh_min <- min(all_pts$dir_cohesion, na.rm = TRUE)
  coh_max <- max(all_pts$dir_cohesion, na.rm = TRUE)

  # Reference point in MIN space (slightly worse than worst = 1)
  ref_min <- c(1.05, 1.05)
  hv_max <- prod(ref_min - c(0, 0))  # max possible HV if front reached (0,0)

  # -----------------------------
  # Compute HV per (algorithm, city)
  # -----------------------------
  hv_results <- do.call(rbind, lapply(names(pareto_by_algo), function(algo) {

    pt_list <- pareto_by_algo[[algo]]
    if (is.null(pt_list) || length(pt_list) == 0) return(NULL)

    do.call(rbind, lapply(names(pt_list), function(city_id) {

      df <- pt_list[[city_id]]
      if (is.null(df) || nrow(df) == 0) {
        return(data.frame(
          algorithm = algo,
          city = city_id,
          n_points = 0,
          hv_raw = NA_real_,
          hv_unit = NA_real_,
          contrib_gini = NA_real_,
          contrib_top10_share = NA_real_,
          contrib_cv = NA_real_
        ))
      }

      df <- df[, c("coverage", "dir_cohesion")]
      df <- df[is.finite(df$coverage) & is.finite(df$dir_cohesion), , drop = FALSE]
      if (nrow(df) == 0) {
        return(data.frame(
          algorithm = algo,
          city = city_id,
          n_points = 0,
          hv_raw = NA_real_,
          hv_unit = NA_real_,
          contrib_gini = NA_real_,
          contrib_top10_share = NA_real_,
          contrib_cv = NA_real_
        ))
      }

      cov_n <- safe_norm(df$coverage, cov_min, cov_max)
      coh_n <- safe_norm(df$dir_cohesion, coh_min, coh_max)

      # MIN space for HV
      Fmin <- cbind(f1 = 1 - cov_n, f2 = 1 - coh_n)

      hv_raw <- hv2d_min(Fmin, ref_min)
      hv_unit <- hv_raw / hv_max

      contrib <- hv2d_contrib_min(Fmin, ref_min)
      contrib <- contrib[is.finite(contrib)]

      top10_share <- if (length(contrib) == 0 || sum(contrib) == 0) NA_real_ else {
        n_top <- max(1, floor(0.1 * length(contrib)))
        sum(sort(contrib, decreasing = TRUE)[seq_len(n_top)]) / sum(contrib)
      }

      data.frame(
        algorithm = algo,
        city = city_id,
        n_points = nrow(df),
        hv_raw = hv_raw,
        hv_unit = hv_unit,
        contrib_gini = gini(contrib),
        contrib_top10_share = top10_share,
        contrib_cv = if (length(contrib) > 1 && mean(contrib) > 0) sd(contrib) / mean(contrib) else NA_real_
      )
    }))
  }))

  hv_results <- hv_results %>%
    arrange(algorithm, desc(hv_unit), city)

  cat("\n===== Hypervolume results (global normalization) =====\n")
  print(hv_results)

  # Convenience summary for plots
  hv_summary_df <- hv_results %>%
    transmute(
      algorithm,
      city,
      hv_unit,
      gini = contrib_gini
    )

  # -----------------------------
  # Optional plots
  # -----------------------------
  if (nrow(hv_summary_df) > 0) {

    print(
      ggplot(hv_summary_df, aes(x = city, y = hv_unit, fill = algorithm)) +
        geom_col(position = "dodge", width = 0.7) +
        labs(
          title = "Unit-scaled hypervolume by city and algorithm",
          x = "City",
          y = "Hypervolume (0–1)"
        ) +
        theme_minimal()
    )

    print(
      ggplot(hv_summary_df, aes(x = city, y = gini, fill = algorithm)) +
        geom_col(position = "dodge", width = 0.7) +
        labs(
          title = "Hypervolume shape (Gini of contributions) by city and algorithm",
          x = "City",
          y = "Gini (0 = even, 1 = concentrated)"
        ) +
        theme_minimal()
    )

    print(
      ggplot(hv_results, aes(x = hv_unit, y = contrib_gini, label = paste0(toupper(city), " ", algorithm))) +
        geom_point(size = 3) +
        geom_text(nudge_y = 0.02, size = 3) +
        labs(
          title = "Global vs shape of Pareto fronts (all cities, both algorithms)",
          x = "Global hypervolume (unit-scaled)",
          y = "Shape (Gini of HV contributions)"
        ) +
        theme_minimal()
    )
  }
}

############################################################
# STEP 14.N — NULL-TEST (DBSCAN): REAL vs RANDOM
# Place AFTER Step 14 (DBSCAN NSGA-II on real) and BEFORE Step 15.
#
# Null model:
# - keep origins fixed
# - permute destination PAIRS (x_d, y_d) with ONE permutation
# - recompute dx, dy
# - rebuild features using the SAME make_X_and_clean() pipeline
# Then:
# - run same NSGA-II on subsample for REAL and RAND
# - compare with unit-scaled HV (global normalization across both tags)
############################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(dbscan)
  library(mco)
  library(ggplot2)
  library(tidyr)
})

# ----------------------------------------------------------
# 0) Preconditions
# ----------------------------------------------------------
prepared_real <- if (exists("prepared_full")) prepared_full else prepared
stopifnot(is.list(prepared_real), length(prepared_real) > 0)

if (!exists("make_X_and_clean")) {
  stop("Function make_X_and_clean() not found. Use your Step 2 feature builder for REAL and RAND.")
}
if (!exists("dir_cohesion_within_weighted")) {
  stop("Function dir_cohesion_within_weighted() not found. Load your cohesion helper before running null-test.")
}

# Use EXACT same feature columns as REAL pipeline
feature_cols_used <- c("x_o", "y_o", "dx", "dy")  # <-- keep identical to your REAL Step 2

# ----------------------------------------------------------
# 1) Null-model generator: permute destination PAIRS
# ----------------------------------------------------------
make_randomized_od <- function(city_clean,
                               origin_cols = c("x_o", "y_o"),
                               dest_cols   = c("x_d", "y_d"),
                               recompute_dxdy = TRUE,
                               seed = 1L) {

  need <- c(origin_cols, dest_cols)
  miss <- setdiff(need, names(city_clean))
  if (length(miss) > 0) {
    stop("Null model needs missing columns in city_clean: ", paste(miss, collapse = ", "))
  }

  set.seed(seed)
  n <- nrow(city_clean)
  perm <- sample.int(n)

  city_rand <- city_clean
  city_rand[[dest_cols[1]]] <- city_clean[[dest_cols[1]]][perm]
  city_rand[[dest_cols[2]]] <- city_clean[[dest_cols[2]]][perm]

  if (recompute_dxdy) {
    city_rand$dx <- city_rand[[dest_cols[1]]] - city_rand[[origin_cols[1]]]
    city_rand$dy <- city_rand[[dest_cols[2]]] - city_rand[[origin_cols[2]]]
  }

  city_rand
}

# ----------------------------------------------------------
# 2) Build prepared_rand using SAME feature pipeline
# ----------------------------------------------------------
prepared_rand <- lapply(names(prepared_real), function(city_id) {

  cc_real <- prepared_real[[city_id]]$city_clean
  stopifnot(all(c("x_o","y_o","x_d","y_d") %in% names(cc_real)))

  cc_rand <- make_randomized_od(cc_real, seed = 100 + match(city_id, names(prepared_real)))

  # rebuild features exactly the same way as REAL
  res <- make_X_and_clean(cc_rand, feature_cols = feature_cols_used)

  # keep same structure as prepared_real (must contain city_clean + X_scaled)
  # adapt here if your make_X_and_clean returns different field names
  if (is.null(res$city_clean)) res$city_clean <- cc_rand
  if (is.null(res$X_scaled) && !is.null(res$X)) res$X_scaled <- res$X  # fallback if your object uses $X

  res
})
names(prepared_rand) <- names(prepared_real)

# Quick sanity: row counts match
for (city_id in names(prepared_real)) {
  if (nrow(prepared_real[[city_id]]$city_clean) != nrow(prepared_rand[[city_id]]$city_clean)) {
    stop("Row mismatch REAL vs RAND for city: ", city_id)
  }
}

# ----------------------------------------------------------
# 3) Sanity checks (optional but recommended)
# ----------------------------------------------------------
sanity_check_randomization <- function(city_real_clean, city_rand_clean,
                                       origin_cols = c("x_o","y_o"),
                                       dest_cols   = c("x_d","y_d"),
                                       n_show = 8) {

  stopifnot(nrow(city_real_clean) == nrow(city_rand_clean))
  n <- nrow(city_real_clean)

  origin_same <- all(city_real_clean[[origin_cols[1]]] == city_rand_clean[[origin_cols[1]]]) &&
                 all(city_real_clean[[origin_cols[2]]] == city_rand_clean[[origin_cols[2]]])

  same_dest_pair <- (city_real_clean[[dest_cols[1]]] == city_rand_clean[[dest_cols[1]]]) &
                    (city_real_clean[[dest_cols[2]]] == city_rand_clean[[dest_cols[2]]])
  match_rate <- mean(same_dest_pair)

  perm_ok_x <- identical(sort(city_real_clean[[dest_cols[1]]]), sort(city_rand_clean[[dest_cols[1]]]))
  perm_ok_y <- identical(sort(city_real_clean[[dest_cols[2]]]), sort(city_rand_clean[[dest_cols[2]]]))

  cat("\n============================\n")
  cat("SANITY CHECK — Randomized OD\n")
  cat("============================\n")
  cat("Rows:", n, "\n")
  cat("Origins identical row-by-row:", origin_same, "\n")
  cat("Destination pair match rate (should be ~0):", round(match_rate * 100, 3), "%\n")
  cat("Dest X is a permutation of original:", perm_ok_x, "\n")
  cat("Dest Y is a permutation of original:", perm_ok_y, "\n")

  show_idx <- seq_len(min(n_show, n))
  preview <- data.frame(
    i = show_idx,
    x_o = city_real_clean[[origin_cols[1]]][show_idx],
    y_o = city_real_clean[[origin_cols[2]]][show_idx],
    x_d_real = city_real_clean[[dest_cols[1]]][show_idx],
    y_d_real = city_real_clean[[dest_cols[2]]][show_idx],
    x_d_rand = city_rand_clean[[dest_cols[1]]][show_idx],
    y_d_rand = city_rand_clean[[dest_cols[2]]][show_idx]
  )
  cat("\nPreview (first rows):\n")
  print(preview, row.names = FALSE)

  invisible(list(
    origin_same = origin_same,
    dest_pair_match_rate = match_rate,
    perm_ok_x = perm_ok_x,
    perm_ok_y = perm_ok_y
  ))
}

for (city_id in names(prepared_real)) {
  cat("\n\n### CITY:", toupper(city_id), "###\n")
  sanity_check_randomization(prepared_real[[city_id]]$city_clean,
                             prepared_rand[[city_id]]$city_clean,
                             n_show = 8)
}

# ----------------------------------------------------------
# 4) Directions (must be built separately for REAL and RAND)
# ----------------------------------------------------------
make_trip_dirs_from_clean <- function(city_clean) {
  stopifnot(all(c("dx","dy") %in% names(city_clean)))
  v <- as.matrix(city_clean[, c("dx","dy")])
  nrm <- sqrt(rowSums(v^2)) + 1e-12
  u <- v / nrm
  colnames(u) <- c("ux","uy")
  u
}

trip_dirs_real <- lapply(names(prepared_real), function(city_id) make_trip_dirs_from_clean(prepared_real[[city_id]]$city_clean))
names(trip_dirs_real) <- names(prepared_real)

trip_dirs_rand <- lapply(names(prepared_rand), function(city_id) make_trip_dirs_from_clean(prepared_rand[[city_id]]$city_clean))
names(trip_dirs_rand) <- names(prepared_rand)

# ----------------------------------------------------------
# 5) Shared subsample indices (same for REAL and RAND)
# ----------------------------------------------------------
make_sub_idx <- function(prepared_in, N_SUB = 10000, seed = 1) {
  set.seed(seed)
  idx <- lapply(names(prepared_in), function(city_id) {
    n <- nrow(prepared_in[[city_id]]$X_scaled)
    if (n > N_SUB) sample.int(n, N_SUB) else seq_len(n)
  })
  names(idx) <- names(prepared_in)
  idx
}

apply_subsample <- function(prepared_in, trip_dirs_in, sub_idx) {
  prep_sub <- lapply(names(prepared_in), function(city_id) {
    ii <- sub_idx[[city_id]]
    list(
      X_scaled = prepared_in[[city_id]]$X_scaled[ii, , drop = FALSE],
      city_clean = prepared_in[[city_id]]$city_clean[ii, , drop = FALSE]
    )
  })
  names(prep_sub) <- names(prepared_in)

  dirs_sub <- lapply(names(trip_dirs_in), function(city_id) {
    ii <- sub_idx[[city_id]]
    trip_dirs_in[[city_id]][ii, , drop = FALSE]
  })
  names(dirs_sub) <- names(trip_dirs_in)

  for (city_id in names(prep_sub)) {
    stopifnot(nrow(prep_sub[[city_id]]$X_scaled) == nrow(dirs_sub[[city_id]]))
  }

  list(prepared_sub = prep_sub, trip_dirs_sub = dirs_sub)
}

N_SUB_NULL <- 10000
sub_idx_null <- make_sub_idx(prepared_real, N_SUB = N_SUB_NULL, seed = 1)

sub_real <- apply_subsample(prepared_real, trip_dirs_real, sub_idx_null)
sub_rand <- apply_subsample(prepared_rand, trip_dirs_rand, sub_idx_null)

# ----------------------------------------------------------
# 6) Objective + NSGA-II runner (DBSCAN)
# ----------------------------------------------------------
eval_dbscan_params <- function(eps, minPts, X_scaled, u_trips,
                               noise_label = 0L, min_cluster_size = 2L) {

  fit <- dbscan::dbscan(X_scaled, eps = eps, minPts = minPts)
  cl <- as.integer(fit$cluster)

  N_total <- length(cl)
  if (N_total == 0) return(c(f1 = 1, f2 = 1))

  clustered <- (cl != noise_label) & !is.na(cl)
  n_clustered <- sum(clustered)
  coverage <- n_clustered / N_total
  if (n_clustered == 0) return(c(f1 = 1, f2 = 1))

  coh <- dir_cohesion_within_weighted(
    labels = cl,
    u = u_trips,
    noise_label = noise_label,
    min_cluster_size = min_cluster_size,
    map_to_01 = TRUE
  )

  c(f1 = -coverage, f2 = -coh)
}

snap_to_allowed <- function(v, allowed) allowed[which.min(abs(allowed - v))]
snap_eps <- function(eps_raw, eps_min, eps_max, eps_step, eps_digits = 3) {
  eps_raw <- max(eps_min, min(eps_raw, eps_max))
  eps <- eps_min + round((eps_raw - eps_min) / eps_step) * eps_step
  eps <- round(eps, eps_digits)
  eps <- max(eps_min, min(eps, eps_max))
  eps
}

make_objfun_city_cached_dbscan <- function(city_id, prepared_sub, trip_dirs_sub,
                                           eps_min, eps_max, eps_step,
                                           minPts_allowed,
                                           eps_digits = 3) {

  X_scaled <- prepared_sub[[city_id]]$X_scaled
  u_trips  <- trip_dirs_sub[[city_id]]
  stopifnot(nrow(X_scaled) == nrow(u_trips))

  cache <- new.env(parent = emptyenv())

  function(x) {
    eps <- snap_eps(x[1], eps_min, eps_max, eps_step, eps_digits)
    mp  <- as.integer(round(x[2]))
    minPts <- snap_to_allowed(mp, minPts_allowed)

    key <- paste0(sprintf(paste0("%.", eps_digits, "f"), eps), "_", minPts)
    if (exists(key, envir = cache, inherits = FALSE)) return(get(key, envir = cache, inherits = FALSE))

    val <- eval_dbscan_params(eps, minPts, X_scaled, u_trips)
    assign(key, val, envir = cache)
    val
  }
}

run_nsga2_dbscan_variant <- function(prepared_sub, trip_dirs_sub, tag,
                                     popsize = 100, generations = 50, seed0 = 1) {

  eps_min  <- 0.08
  eps_max  <- 0.45
  eps_step <- 0.02
  eps_digits <- 3
  minPts_allowed <- c(6L, 8L, 10L, 12L, 15L, 20L, 25L, 30L, 40L, 50L, 60L, 70L)

  pareto <- list()

  for (i in seq_along(names(prepared_sub))) {
    city_id <- names(prepared_sub)[i]
    cat("\n[NULL-TEST][", tag, "] NSGA-II:", toupper(city_id), "\n")
    t0 <- Sys.time()

    objfun <- make_objfun_city_cached_dbscan(
      city_id, prepared_sub, trip_dirs_sub,
      eps_min, eps_max, eps_step,
      minPts_allowed, eps_digits
    )

    set.seed(seed0 + i)
    res <- mco::nsga2(
      fn = objfun,
      idim = 2, odim = 2,
      lower.bounds = as.numeric(c(eps_min, min(minPts_allowed))),
      upper.bounds = as.numeric(c(eps_max, max(minPts_allowed))),
      popsize = popsize,
      generations = generations
    )

    eps_raw <- res$par[, 1]
    mp_raw  <- as.integer(round(res$par[, 2]))

    eps_snap <- vapply(eps_raw, snap_eps, numeric(1),
                       eps_min = eps_min, eps_max = eps_max,
                       eps_step = eps_step, eps_digits = eps_digits)
    mp_snap  <- vapply(mp_raw, snap_to_allowed, integer(1), allowed = minPts_allowed)

    df <- data.frame(
      city = city_id, tag = tag,
      eps = as.numeric(eps_snap),
      minPts = as.integer(mp_snap),
      coverage = -res$value[, 1],
      dir_cohesion = -res$value[, 2]
    ) %>%
      filter(is.finite(coverage), is.finite(dir_cohesion),
             coverage >= 0, dir_cohesion >= 0) %>%
      distinct(eps, minPts, .keep_all = TRUE)

    pareto[[city_id]] <- df

    cat("[NULL-TEST][", tag, "] Done:", toupper(city_id),
        "| pareto_n =", nrow(df),
        "| minutes =", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2), "\n")
  }

  list(pareto = pareto)
}

out_real <- run_nsga2_dbscan_variant(sub_real$prepared_sub, sub_real$trip_dirs_sub, tag = "real",
                                    popsize = 100, generations = 50, seed0 = 10)

out_rand <- run_nsga2_dbscan_variant(sub_rand$prepared_sub, sub_rand$trip_dirs_sub, tag = "rand",
                                    popsize = 100, generations = 50, seed0 = 20)

# ----------------------------------------------------------
# 7) Hypervolume comparison (REAL vs RAND)
# ----------------------------------------------------------
nondominated_min <- function(F) {
  n <- nrow(F); keep <- rep(TRUE, n)
  for (i in seq_len(n)) {
    if (!keep[i]) next
    for (j in seq_len(n)) {
      if (i == j || !keep[i]) next
      if (all(F[j, ] <= F[i, ]) && any(F[j, ] < F[i, ])) keep[i] <- FALSE
    }
  }
  F[keep, , drop = FALSE]
}

hv2d_min <- function(F, ref) {
  if (nrow(F) == 0) return(0)
  F <- nondominated_min(F)
  o <- order(F[, 1], F[, 2])
  F <- F[o, , drop = FALSE]
  hv <- 0; h <- ref[2]
  for (i in seq_len(nrow(F))) {
    f1 <- F[i, 1]; f2 <- F[i, 2]
    if (f2 < h) { hv <- hv + (ref[1] - f1) * (h - f2); h <- f2 }
  }
  hv
}

safe_norm <- function(z, zmin, zmax) {
  if (!is.finite(zmin) || !is.finite(zmax) || isTRUE(all.equal(zmax, zmin))) return(rep(0, length(z)))
  (z - zmin) / (zmax - zmin)
}

pareto_by_tag <- list(real = out_real$pareto, rand = out_rand$pareto)

all_pts <- do.call(rbind, lapply(names(pareto_by_tag), function(tag) {
  pts <- pareto_by_tag[[tag]]
  do.call(rbind, lapply(names(pts), function(city_id) {
    df <- pts[[city_id]]
    if (is.null(df) || nrow(df) == 0) return(NULL)
    df <- df[, c("coverage", "dir_cohesion")]
    df <- df[is.finite(df$coverage) & is.finite(df$dir_cohesion), , drop = FALSE]
    if (nrow(df) == 0) return(NULL)
    df$city <- city_id; df$tag <- tag
    df
  }))
}))

if (!is.null(all_pts) && nrow(all_pts) > 0) {

  cov_min <- min(all_pts$coverage, na.rm = TRUE)
  cov_max <- max(all_pts$coverage, na.rm = TRUE)
  coh_min <- min(all_pts$dir_cohesion, na.rm = TRUE)
  coh_max <- max(all_pts$dir_cohesion, na.rm = TRUE)

  ref_min <- c(1.05, 1.05)
  hv_max  <- prod(ref_min - c(0, 0))

  hv_results_tag <- do.call(rbind, lapply(names(pareto_by_tag), function(tag) {
    pts <- pareto_by_tag[[tag]]
    do.call(rbind, lapply(names(pts), function(city_id) {
      df <- pts[[city_id]]
      df <- df[, c("coverage", "dir_cohesion")]
      df <- df[is.finite(df$coverage) & is.finite(df$dir_cohesion), , drop = FALSE]
      if (nrow(df) == 0) return(data.frame(city=city_id, tag=tag, hv_unit=NA_real_))

      cov_n <- safe_norm(df$coverage, cov_min, cov_max)
      coh_n <- safe_norm(df$dir_cohesion, coh_min, coh_max)
      Fmin  <- cbind(f1 = 1 - cov_n, f2 = 1 - coh_n)

      hv_unit <- hv2d_min(Fmin, ref_min) / hv_max
      data.frame(city = city_id, tag = tag, hv_unit = hv_unit)
    }))
  }))

  hv_results_tag <- hv_results_tag %>%
    mutate(tag = factor(tag, levels = c("real","rand"))) %>%
    arrange(city, tag)

  print(hv_results_tag)

  p_null <- ggplot(hv_results_tag, aes(x = city, y = hv_unit, fill = tag)) +
    geom_col(position = "dodge", width = 0.65) +
    labs(
      title = "Null-test: Unit-scaled hypervolume (REAL vs RANDOM)",
      subtitle = paste0("DBSCAN NSGA-II on subsample (N_SUB=", N_SUB_NULL, "), ref=", paste(ref_min, collapse = ",")),
      x = "City", y = "Hypervolume (0–1)", fill = "Dataset"
    ) +
    theme_minimal()

  print(p_null)
} else {
  cat("No valid Pareto points available for HV analysis.\n")
}

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(dbscan)
})

# Use full-data prepared if available, else subsample prepared
prepared_plot <- if (exists("prepared_full")) prepared_full else prepared
stopifnot(is.list(prepared_plot), length(prepared_plot) > 0)

# ----------------------------------------------------------
# Plot: top-N clusters (optional grey background)
# ----------------------------------------------------------
plot_flow_topN <- function(city_id, prepared, labels,
                           algo_name = "ALGO",
                           noise_label = 0L,
                           out_dir = "plots_flows_top3",
                           top_n_clusters = 3,
                           sample_n_segments_top = 25000,
                           sample_n_segments_bg  = 30000,
                           show_background = FALSE,
                           seed = 1,
                           show_plots = TRUE) {

  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  dt <- as.data.table(prepared[[city_id]]$city_clean)
  req <- c("x_o","y_o","x_d","y_d")
  if (!all(req %in% names(dt))) stop("Missing columns in city_clean: ", paste(req, collapse=", "))

  if (length(labels) != nrow(dt)) {
    stop("labels length != nrow(city_clean) for ", city_id,
         " | labels=", length(labels), " rows=", nrow(dt))
  }

  dt[, cluster := as.integer(labels)]

  # Top clusters by size (exclude noise/NA)
  cl_nonnoise <- dt$cluster[!is.na(dt$cluster) & dt$cluster != noise_label]
  tab <- sort(table(cl_nonnoise), decreasing = TRUE)
  top_ids <- as.integer(names(head(tab, top_n_clusters)))

  if (length(top_ids) == 0) {
    warning("No non-noise clusters for ", city_id, " (", algo_name, ").")
    return(NULL)
  }

  p <- ggplot()

  # Optional background: all clustered flows in grey
  if (isTRUE(show_background)) {
    dt_bg <- dt[!is.na(cluster) & cluster != noise_label]
    if (nrow(dt_bg) > 0) {
      set.seed(seed)
      idx_bg <- if (nrow(dt_bg) > sample_n_segments_bg) sample.int(nrow(dt_bg), sample_n_segments_bg) else seq_len(nrow(dt_bg))
      dt_bg <- dt_bg[idx_bg, ]
      p <- p + geom_segment(
        data = dt_bg,
        aes(x = x_o, y = y_o, xend = x_d, yend = y_d),
        alpha = 0.05, linewidth = 0.20
      )
    }
  }

  # Top clusters in color
  dt_top <- dt[cluster %in% top_ids]
  set.seed(seed)
  idx_top <- if (nrow(dt_top) > sample_n_segments_top) sample.int(nrow(dt_top), sample_n_segments_top) else seq_len(nrow(dt_top))
  dt_top <- dt_top[idx_top, ]
  dt_top[, cluster_f := factor(cluster, levels = top_ids)]

  p <- p +
    geom_segment(
      data = dt_top,
      aes(x = x_o, y = y_o, xend = x_d, yend = y_d, color = cluster_f),
      alpha = 0.28, linewidth = 0.30
    ) +
    coord_equal() +
    labs(
      title = paste0(algo_name, " — Top ", top_n_clusters, " clusters — ", toupper(city_id)),
      subtitle = if (show_background)
        "Colored = top clusters | Grey = all clustered flows (sampled)"
      else
        "Colored = top clusters (sampled)",
      x = "UTM Easting (m)", y = "UTM Northing (m)", color = "Cluster"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "right"
    )

  f_out <- file.path(out_dir, paste0("flows_", algo_name, "_", city_id, "_top", top_n_clusters,
                                    ifelse(show_background, "_bg", ""), ".png"))
  ggsave(f_out, plot = p, width = 9, height = 7.5, dpi = 300)

  if (isTRUE(show_plots)) print(p)

  list(file = f_out, top_cluster_ids = top_ids, plot = p)
}

# ----------------------------------------------------------
# Select final configuration per city (knee_full preferred)
# ----------------------------------------------------------
get_dbscan_selected_params <- function(preferred = "knee_full") {
  out <- list()
  for (city_id in names(prepared_plot)) {

    if (preferred == "knee_full" && exists("knee_full_dbscan") &&
        is.data.frame(knee_full_dbscan) && nrow(knee_full_dbscan) > 0) {
      row <- knee_full_dbscan %>% filter(city == city_id) %>% slice(1)
      if (nrow(row) == 1) {
        out[[city_id]] <- list(eps = as.numeric(row$eps), minPts = as.integer(row$minPts))
        next
      }
    }

    if (exists("dbscan_params") && !is.null(dbscan_params[[city_id]])) {
      out[[city_id]] <- list(eps = as.numeric(dbscan_params[[city_id]]$eps),
                             minPts = as.integer(dbscan_params[[city_id]]$minPts))
      next
    }

    stop("[DBSCAN] No selected parameters found for city: ", city_id,
         " (need knee_full_dbscan or dbscan_params).")
  }
  out
}

get_snn_selected_params <- function(preferred = "knee_full") {
  out <- list()
  for (city_id in names(prepared_plot)) {

    if (preferred == "knee_full" && exists("knee_full_snn") &&
        is.data.frame(knee_full_snn) && nrow(knee_full_snn) > 0) {
      row <- knee_full_snn %>% filter(city == city_id) %>% slice(1)
      if (nrow(row) == 1) {
        out[[city_id]] <- list(
          k = as.integer(row$k),
          snn_threshold = as.integer(row$snn_threshold),
          minPts = as.integer(row$minPts)
        )
        next
      }
    }

    if (exists("snn_params") && !is.null(snn_params[[city_id]])) {
      out[[city_id]] <- list(
        k = as.integer(snn_params[[city_id]]$k),
        snn_threshold = as.integer(snn_params[[city_id]]$snn_threshold),
        minPts = as.integer(snn_params[[city_id]]$minPts)
      )
      next
    }

    stop("[SNN] No selected parameters found for city: ", city_id,
         " (need knee_full_snn or snn_params).")
  }
  out
}

selected_params_dbscan <- get_dbscan_selected_params(preferred = "knee_full")
selected_params_snn    <- get_snn_selected_params(preferred = "knee_full")

# ----------------------------------------------------------
# Run + plot DBSCAN
# ----------------------------------------------------------
for (city_id in names(prepared_plot)) {
  p <- selected_params_dbscan[[city_id]]
  X <- prepared_plot[[city_id]]$X_scaled

  fit <- dbscan::dbscan(X, eps = p$eps, minPts = p$minPts)
  labels <- as.integer(fit$cluster)

  plot_flow_topN(
    city_id = city_id,
    prepared = prepared_plot,
    labels = labels,
    algo_name = "DBSCAN",
    out_dir = "plots_flows_top3",
    top_n_clusters = 3,
    show_background = FALSE,
    show_plots = TRUE
  )
}

# ----------------------------------------------------------
# Run + plot SNN
# ----------------------------------------------------------
for (city_id in names(prepared_plot)) {
  p <- selected_params_snn[[city_id]]
  X <- prepared_plot[[city_id]]$X_scaled

  snn_obj <- dbscan::sNN(X, k = p$k, sort = FALSE, jp = FALSE)
  fit <- dbscan::sNNclust(
    x = snn_obj,
    k = p$k,
    eps = p$snn_threshold,
    minPts = p$minPts,
    borderPoints = TRUE
  )
  labels <- as.integer(fit$cluster)

  plot_flow_topN(
    city_id = city_id,
    prepared = prepared_plot,
    labels = labels,
    algo_name = "SNN",
    out_dir = "plots_flows_top3",
    top_n_clusters = 3,
    show_background = FALSE,
    show_plots = TRUE
  )
}                                                        
