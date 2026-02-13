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
