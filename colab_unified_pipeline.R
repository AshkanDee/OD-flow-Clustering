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
# WHY:
# - Quick schema check across cities in one run
# - Confirms consistency before clustering
# - Keeps each city dataset separate (no merge/stack)

city_paths <- list(
  berlin  = "/content/drive/MyDrive/berlin_flows_utm.rds",
  munich  = "/content/drive/MyDrive/munich_flows_utm.rds",
  cologne = "/content/drive/MyDrive/cologne_flows_utm.rds"
)

city_data <- lapply(names(city_paths), function(city_id) {
  path <- city_paths[[city_id]]
  if (!file.exists(path)) {
    stop(sprintf("Missing file for %s: %s", city_id, path))
  }

  dt <- as.data.table(readRDS(path))
  dt[, city := city_id]  # traceability tag
  dt
})
names(city_data) <- names(city_paths)

print_city_inspection <- function(city_id, dt) {
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

invisible(Map(print_city_inspection, names(city_data), city_data))
