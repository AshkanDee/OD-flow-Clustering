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
