# Prevent scientific notation globally
options(scipen = 999)

# Install vctrs first (do not load any other packages yet)
if (!requireNamespace("vctrs", quietly = TRUE)) {
  install.packages("vctrs")
}

# Install other dependent packages if missing
deps <- c("dbscan", "ggplot2", "dplyr")
missing <- deps[!vapply(deps, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  install.packages(missing)
}

# Load packages AFTER installation
library(vctrs)
library(dbscan)
library(ggplot2)
library(dplyr)
