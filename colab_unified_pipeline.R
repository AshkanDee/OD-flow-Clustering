# Unified OD-flow pipeline scaffold (DBSCAN + HDBSCAN + SNN)
# This block keeps city datasets separate and validates schema before clustering.

# ---- Package setup ----
options(repos = c(CRAN = "https://cloud.r-project.org"))
options(scipen = 999)  # Prevent scientific notation globally

required_pkgs <- c(
  "vctrs",
  "data.table",
  "dbscan",
  "ggplot2",
  "dplyr",
  "cluster",
  "tibble",
  "tidyr",
  "purrr"
)

missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) install.packages(missing_pkgs)
invisible(lapply(required_pkgs, library, character.only = TRUE))

# ---- Step 1: Load + inspect all cities (no combining) ----
city_paths <- list(
  berlin  = "/content/drive/MyDrive/berlin_flows_utm.rds",
  munich  = "/content/drive/MyDrive/munich_flows_utm.rds",
  cologne = "/content/drive/MyDrive/cologne_flows_utm.rds"
)

required_columns <- c("x_o", "y_o", "x_d", "y_d", "dx", "dy", "len", "angle")

city_data <- lapply(names(city_paths), function(city_id) {
  path <- city_paths[[city_id]]
  if (!file.exists(path)) {
    stop(sprintf("Missing file for %s: %s", city_id, path))
  }

  dt <- as.data.table(readRDS(path))
  dt[, city := city_id]

  missing_cols <- setdiff(required_columns, names(dt))
  if (length(missing_cols) > 0) {
    warning(sprintf("%s missing columns: %s", city_id, paste(missing_cols, collapse = ", ")))
  }

  dt
})
names(city_data) <- names(city_paths)

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

  cat("\nmissing required columns:\n")
  print(setdiff(required_columns, names(dt)))

  cat("\nNA counts in key fields:\n")
  key_cols <- intersect(required_columns, names(dt))
  if (length(key_cols) > 0) {
    print(dt[, lapply(.SD, function(x) sum(is.na(x))), .SDcols = key_cols])
  }
}
