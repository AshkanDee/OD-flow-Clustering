# Please if dataset files are not in your current working directory, change them to your path, codes in line 163 and 2109
# Tidy OD-flow clustering pipeline (DBSCAN + SNN)
# ------------------------------------------------
# This script consolidates shared logic that was previously duplicated
# across Colab blocks and keeps city datasets separate throughout.

# -------------------------
# 0) Package setup
# -------------------------
setup_packages <- function(pkgs = c("data.table", "sf", "dbscan", "dplyr", "ggplot2"),
                           repos = "https://cloud.r-project.org") {
  options(repos = c(CRAN = repos))

  missing <- pkgs[!vapply(pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
  if (length(missing) > 0) {
    message("Installing missing packages: ", paste(missing, collapse = ", "))
    install.packages(missing, dependencies = TRUE)
  }

  still_missing <- pkgs[!vapply(pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
  if (length(still_missing) > 0) {
    stop(
      "Failed to install/load required packages: ", paste(still_missing, collapse = ", "),
      "\nInstall them manually, then rerun."
    )
  }

  for (pkg in pkgs) {
    suppressPackageStartupMessages(library(pkg, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE))
  }

  invisible(TRUE)
}

# -------------------------
# 1) OD preprocessing
# -------------------------
prep_city_od <- function(rds_path,
                         needed = c("start_loc_lon", "start_loc_lat", "dest_loc_lon", "dest_loc_lat"),
                         crs_in = 4326,
                         crs_out = 25832,
                         drop_zero_len = TRUE) {
  dt <- data.table::as.data.table(readRDS(rds_path))

  missing <- setdiff(needed, names(dt))
  if (length(missing) > 0) {
    stop("Missing columns in ", rds_path, ": ", paste(missing, collapse = ", "))
  }

  n_before <- nrow(dt)

  dt <- dt[
    !is.na(start_loc_lat) & !is.na(start_loc_lon) &
      !is.na(dest_loc_lat) & !is.na(dest_loc_lon)
  ]
  n_after_na <- nrow(dt)

  dt <- dt[
    start_loc_lat >= -90 & start_loc_lat <= 90 &
      dest_loc_lat >= -90 & dest_loc_lat <= 90 &
      start_loc_lon >= -180 & start_loc_lon <= 180 &
      dest_loc_lon >= -180 & dest_loc_lon <= 180
  ]
  n_after_range <- nrow(dt)

  orig_sf <- sf::st_as_sf(dt, coords = c("start_loc_lon", "start_loc_lat"), crs = crs_in, remove = FALSE)
  dest_sf <- sf::st_as_sf(dt, coords = c("dest_loc_lon", "dest_loc_lat"), crs = crs_in, remove = FALSE)

  orig_xy <- sf::st_coordinates(sf::st_transform(orig_sf, crs_out))
  dest_xy <- sf::st_coordinates(sf::st_transform(dest_sf, crs_out))

  dt[, `:=`(
    x_o = orig_xy[, 1],
    y_o = orig_xy[, 2],
    x_d = dest_xy[, 1],
    y_d = dest_xy[, 2]
  )]

  dt[, `:=`(dx = x_d - x_o, dy = y_d - y_o)]
  dt[, `:=`(angle = atan2(dy, dx), len = sqrt(dx^2 + dy^2))]

  n_zero_len <- dt[len == 0, .N]
  if (drop_zero_len) dt <- dt[len > 0]

  list(
    data = dt,
    info = list(
      file = rds_path,
      n_before = n_before,
      n_after_na = n_after_na,
      n_after_range = n_after_range,
      n_zero_len = n_zero_len,
      n_after_zero = nrow(dt)
    )
  )
}

prepare_all_cities <- function(city_files,
                               needed = c("start_loc_lon", "start_loc_lat", "dest_loc_lon", "dest_loc_lat")) {
  prepared <- lapply(city_files, prep_city_od, needed = needed, drop_zero_len = TRUE)

  cleaning_summary <- data.table::rbindlist(lapply(names(prepared), function(city) {
    i <- prepared[[city]]$info
    data.table::data.table(
      city = city,
      file = i$file,
      n_before = i$n_before,
      n_after_na = i$n_after_na,
      n_after_range = i$n_after_range,
      zero_len_trips = i$n_zero_len,
      n_after_zero_len_removed = i$n_after_zero
    )
  }), fill = TRUE)

  list(prepared = prepared, cleaning_summary = cleaning_summary)
}

# -------------------------
# 2) Feature matrix prep
# -------------------------
build_city_data <- function(prepared_od) {
  city_data <- lapply(names(prepared_od), function(city_id) {
    dt <- prepared_od[[city_id]]$data
    data.table::setDT(dt)
    if (!("city" %in% names(dt))) dt[, city := city_id]
    dt
  })
  names(city_data) <- names(prepared_od)
  city_data
}

make_X_and_clean <- function(city_dt, feature_cols = c("x_o", "y_o", "dx", "dy")) {
  missing_cols <- setdiff(feature_cols, names(city_dt))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  X <- as.matrix(city_dt[, ..feature_cols])
  ok <- complete.cases(X)

  list(
    X = X[ok, , drop = FALSE],
    ok = ok,
    X_scaled = scale(X[ok, , drop = FALSE]),
    city_clean = city_dt[ok, ]
  )
}

prepare_feature_matrices <- function(city_data, feature_cols = c("x_o", "y_o", "dx", "dy")) {
  out <- lapply(names(city_data), function(city_id) {
    res <- make_X_and_clean(city_data[[city_id]], feature_cols = feature_cols)
    res$city_id <- city_id
    res
  })
  names(out) <- names(city_data)
  out
}





# Default dataset mapping used when city_files is omitted
default_city_files <- function(base_dir = ".") {
  list(
    berlin = file.path(base_dir, "dt_bolt_berlin_06_05.rds"),
    cologne = file.path(base_dir, "dt_voi_cologne_06_05.rds"),
    munich = file.path(base_dir, "dt_voi_munich_06_05.rds")
  )
}

resolve_city_files <- function(city_files = NULL, base_dir = ".", extra_dirs = NULL) {
  if (!is.null(city_files)) return(city_files)

  candidate_dirs <- unique(c(base_dir, getwd(), extra_dirs))
  candidate_dirs <- candidate_dirs[!is.na(candidate_dirs) & nzchar(candidate_dirs)]

  for (d in candidate_dirs) {
    guessed <- default_city_files(base_dir = d)
    if (all(file.exists(unlist(guessed)))) {
      message("Using default city_files from: ", normalizePath(d, winslash = "/", mustWork = FALSE))
      return(guessed)
    }
  }

  stop(
    "city_files not provided and default files were not all found.
",
    "Checked directories: ", paste(normalizePath(candidate_dirs, winslash = "/", mustWork = FALSE), collapse = " | "), "
",
    "Provide explicit paths, e.g.:
",
    "city_files <- list(berlin='.../dt_bolt_berlin_06_05.rds', cologne='.../dt_voi_cologne_06_05.rds', munich='.../dt_voi_munich_06_05.rds')"
  )
}

# -------------------------
# Quick start: load city datasets first
# -------------------------
# 1) Define your city files (edit paths):
# city_files <- list(
#   berlin  = "dt_bolt_berlin_06_05.rds",
#   cologne = "dt_voi_cologne_06_05.rds",
#   munich  = "dt_voi_munich_06_05.rds"
# )
#
# 2) Initialize core objects used by all later steps:
# init <- initialize_pipeline_state(city_files)  # or initialize_pipeline_state() if default files exist in wd
# prepared_od <- init$prepared_od
# city_data <- init$city_data
# prepared <- init$prepared
# cleaning_summary <- init$cleaning_summary
# print(cleaning_summary)

initialize_pipeline_state <- function(city_files = NULL,
                                      base_dir = ".",
                                      needed = c("start_loc_lon", "start_loc_lat", "dest_loc_lon", "dest_loc_lat"),
                                      feature_cols = c("x_o", "y_o", "dx", "dy")) {
  city_files <- resolve_city_files(city_files, base_dir = base_dir)

  if (is.null(names(city_files)) || any(names(city_files) == "")) {
    stop("initialize_pipeline_state: `city_files` must be named (e.g., berlin, cologne, munich).")
  }

  missing_paths <- city_files[!file.exists(unlist(city_files))]
  if (length(missing_paths) > 0) {
    stop("initialize_pipeline_state: missing files -> ", paste(unname(missing_paths), collapse = ", "))
  }

  stage0 <- prepare_all_cities(city_files = city_files, needed = needed)
  prepared_od <- stage0$prepared
  city_data <- build_city_data(prepared_od)
  prepared <- prepare_feature_matrices(city_data, feature_cols = feature_cols)

  list(
    prepared_od = prepared_od,
    city_data = city_data,
    prepared = prepared,
    cleaning_summary = stage0$cleaning_summary
  )
}

# -------------------------
# 3) Diagnostics helpers
# -------------------------
knee_elbow <- function(d_sorted) {
  d <- as.numeric(d_sorted)
  n <- length(d)
  if (n < 10) return(NA_real_)

  x <- seq_len(n)
  x1 <- 1
  y1 <- d[1]
  x2 <- n
  y2 <- d[n]

  num <- abs((y2 - y1) * x - (x2 - x1) * d + x2 * y1 - y2 * x1)
  den <- sqrt((y2 - y1)^2 + (x2 - x1)^2)
  d[which.max(num / den)]
}

knn_summary_one <- function(X_scaled, k,
                            q = c(0.01, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99)) {
  d <- dbscan::kNNdist(X_scaled, k = k)
  out <- data.frame(
    k = k,
    n = length(d),
    elbow_eps = knee_elbow(sort(d)),
    t(stats::quantile(d, probs = q, na.rm = TRUE, names = FALSE))
  )
  names(out) <- c("k", "n", "elbow_eps", paste0("q", sprintf("%02d", round(q * 100))))
  out
}

snn_sharednn_diagnostics <- function(X_scaled, k = 30, n_probe = 4000, seed = 1) {
  set.seed(seed)
  n <- nrow(X_scaled)
  idx <- if (n > n_probe) sample.int(n, n_probe) else seq_len(n)

  knn <- dbscan::kNN(X_scaled, k = k)
  nn_list <- lapply(idx, function(i) knn$id[i, ])

  checks <- min(8000, length(nn_list) * 2)
  shared_counts <- integer(checks)

  for (t in seq_len(checks)) {
    a <- sample.int(length(nn_list), 1)
    b <- sample.int(length(nn_list), 1)
    shared_counts[t] <- length(intersect(nn_list[[a]], nn_list[[b]]))
  }

  data.frame(
    k = k,
    probe_n = length(nn_list),
    shared_mean = mean(shared_counts),
    shared_median = stats::median(shared_counts),
    shared_q90 = unname(stats::quantile(shared_counts, 0.90)),
    shared_q95 = unname(stats::quantile(shared_counts, 0.95)),
    shared_q99 = unname(stats::quantile(shared_counts, 0.99))
  )
}

intersect_count_sorted <- function(a, b) {
  i <- 1L
  j <- 1L
  na <- length(a)
  nb <- length(b)
  cnt <- 0L

  while (i <= na && j <= nb) {
    av <- a[i]
    bv <- b[j]
    if (av == bv) {
      cnt <- cnt + 1L
      i <- i + 1L
      j <- j + 1L
    } else if (av < bv) {
      i <- i + 1L
    } else {
      j <- j + 1L
    }
  }
  cnt
}

snn_kth_shared_curve <- function(X_scaled, k) {
  knn <- dbscan::kNN(X_scaled, k = k)
  nn_id <- t(apply(knn$id, 1, sort))

  n <- nrow(nn_id)
  shared_kth <- integer(n)
  for (i in seq_len(n)) {
    Ni <- nn_id[i, ]
    j <- Ni[k]
    Nj <- nn_id[j, ]
    shared_kth[i] <- intersect_count_sorted(Ni, Nj)
  }

  sort(shared_kth)
}

# -------------------------
# 4) Shared evaluation helpers
# -------------------------
build_metrics_row <- function(city_id, algorithm, cl, n_total, extra = list()) {
  n_noise <- sum(cl == 0)
  n_clusters <- length(unique(cl[cl != 0]))
  max_cluster_size <- if (n_clusters > 0) max(table(cl[cl != 0])) else 0

  base <- list(
    city = city_id,
    algorithm = algorithm,
    n_clusters = n_clusters,
    noise_pct = 100 * n_noise / length(cl),
    max_cluster_pct = if (n_clusters > 0) 100 * max_cluster_size / length(cl) else 0,
    n_total = n_total,
    n_used = n_total
  )

  as.data.frame(c(base, extra), stringsAsFactors = FALSE)
}

run_dbscan_grid_search <- function(prepared, eps_grid, minPts_grid_dbscan) {
  rows <- list()
  idx <- 1L

  for (city_id in names(prepared)) {
    X_scaled <- prepared[[city_id]]$X_scaled
    if (is.null(X_scaled)) stop("[DBSCAN] Missing X_scaled for city: ", city_id)

    n_total <- nrow(X_scaled)
    for (minPts in minPts_grid_dbscan) {
      for (eps_val in eps_grid) {
        fit <- dbscan::dbscan(X_scaled, eps = eps_val, minPts = minPts)
        rows[[idx]] <- build_metrics_row(
          city_id = city_id,
          algorithm = "DBSCAN",
          cl = fit$cluster,
          n_total = n_total,
          extra = list(minPts = minPts, eps = eps_val)
        )
        idx <- idx + 1L
      }
    }
  }

  dplyr::bind_rows(rows)
}

run_snn_grid_search <- function(prepared, k_grid_snn, make_snn_threshold_grid, minPts_grid_snn) {
  rows <- list()
  idx <- 1L

  for (city_id in names(prepared)) {
    X_scaled <- prepared[[city_id]]$X_scaled
    if (is.null(X_scaled)) stop("[SNN] Missing X_scaled for city: ", city_id)

    n_total <- nrow(X_scaled)

    for (k in k_grid_snn) {
      snn_obj <- dbscan::sNN(X_scaled, k = k, sort = FALSE, jp = FALSE)
      thresholds <- make_snn_threshold_grid(k)

      for (snn_thr in thresholds) {
        for (minPts in minPts_grid_snn) {
          fit <- dbscan::sNNclust(
            x = snn_obj,
            k = k,
            eps = snn_thr,
            minPts = minPts,
            borderPoints = TRUE
          )

          rows[[idx]] <- build_metrics_row(
            city_id = city_id,
            algorithm = "SNN",
            cl = fit$cluster,
            n_total = n_total,
            extra = list(jp = FALSE, k = k, snn_threshold = snn_thr, minPts = minPts, eps = NA_real_)
          )
          idx <- idx + 1L
        }
      }
    }
  }

  dplyr::bind_rows(rows)
}

# -------------------------
# 5) End-to-end runner
# -------------------------
run_pipeline <- function(city_files = NULL, base_dir = ".") {
  setup_packages()

  init <- initialize_pipeline_state(city_files, base_dir = base_dir)
  prepared <- init$prepared

  # Suggested search space from your original pipeline
  eps_grid <- seq(0.08, 0.40, by = 0.02)
  minPts_grid_dbscan <- c(10, 20, 30)

  k_grid_snn <- c(10, 20, 30)
  make_snn_threshold_grid <- function(k) {
    thr <- round(seq(0.35 * k, 0.65 * k, length.out = 6))
    unique(pmax(2L, pmin(k - 1L, as.integer(thr))))
  }
  minPts_grid_snn <- c(15, 18, 20, 25)

  grid_df_dbscan <- run_dbscan_grid_search(prepared, eps_grid, minPts_grid_dbscan)
  grid_df_snn <- run_snn_grid_search(prepared, k_grid_snn, make_snn_threshold_grid, minPts_grid_snn)

  list(
    cleaning_summary = init$cleaning_summary,
    prepared_od = init$prepared_od,
    city_data = init$city_data,
    prepared = prepared,
    grid_df_dbscan = grid_df_dbscan,
    grid_df_snn = grid_df_snn
  )
}

# Example usage:
# city_files <- list(
#   berlin  = "dt_bolt_berlin_06_05.rds",
#   cologne = "dt_voi_cologne_06_05.rds",
#   munich  = "dt_voi_munich_06_05.rds"
# )
# out <- run_pipeline(city_files)
# print(out$cleaning_summary)
# print(out$grid_df_dbscan)
# print(out$grid_df_snn)


# -------------------------
# 6) Final per-city runs + post-clustering analysis (Steps 7-12)
# -------------------------

apply_final_dbscan <- function(prepared,
                               dbscan_params = list(
                                 berlin = list(minPts = 20, eps = 0.200),
                                 munich = list(minPts = 20, eps = 0.200),
                                 cologne = list(minPts = 10, eps = 0.250)
                               ),
                               label_col = "cluster") {
  for (city_id in names(prepared)) {
    if (is.null(prepared[[city_id]]$X_scaled)) stop("[DBSCAN] Missing X_scaled for city: ", city_id)
    if (is.null(prepared[[city_id]]$city_clean)) stop("[DBSCAN] Missing city_clean for city: ", city_id)
    if (is.null(dbscan_params[[city_id]])) stop("[DBSCAN] Missing final params for city: ", city_id)

    X_scaled <- prepared[[city_id]]$X_scaled
    city_tbl <- prepared[[city_id]]$city_clean
    stopifnot(nrow(X_scaled) == nrow(city_tbl))

    params <- dbscan_params[[city_id]]
    fit <- dbscan::dbscan(X_scaled, eps = params$eps, minPts = params$minPts)

    stopifnot(length(fit$cluster) == nrow(city_tbl))
    prepared[[city_id]]$city_clean[[label_col]] <- fit$cluster

    cat(
      "DBSCAN done |", city_id,
      "| eps =", params$eps,
      "| minPts =", params$minPts,
      "| clusters =", length(unique(fit$cluster[fit$cluster != 0])),
      "| noise% =", round(100 * mean(fit$cluster == 0), 2),
      "\n"
    )
  }

  prepared
}

apply_final_snn <- function(prepared,
                            snn_params = list(
                              berlin = list(k = 20, snn_threshold = 8, minPts = 12),
                              munich = list(k = 30, snn_threshold = 11, minPts = 18),
                              cologne = list(k = 15, snn_threshold = 6, minPts = 10)
                            ),
                            label_col = "cluster") {
  for (city_id in names(prepared)) {
    if (is.null(prepared[[city_id]]$X_scaled)) stop("[SNN] Missing X_scaled for city: ", city_id)
    if (is.null(prepared[[city_id]]$city_clean)) stop("[SNN] Missing city_clean for city: ", city_id)
    if (is.null(snn_params[[city_id]])) stop("[SNN] Missing final params for city: ", city_id)

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
    prepared[[city_id]]$city_clean[[label_col]] <- fit$cluster

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

  prepared
}

verify_cluster_sizes <- function(prepared, label_col = "cluster", n_head = 65) {
  for (city_id in names(prepared)) {
    cat("\n---", toupper(city_id), "---\n")

    if (is.null(prepared[[city_id]]$city_clean)) stop("[verify_cluster_sizes] Missing city_clean for city: ", city_id)
    city_dt <- prepared[[city_id]]$city_clean
    if (!(label_col %in% names(city_dt))) stop("[verify_cluster_sizes] Missing label column for city: ", city_id)

    cl_tab <- table(city_dt[[label_col]])
    print(utils::head(cl_tab, n_head))
  }

  invisible(TRUE)
}

pick_col <- function(df, candidates) {
  hit <- candidates[candidates %in% names(df)]
  if (length(hit) == 0) return(NA_character_)
  hit[[1]]
}

build_cluster_summaries <- function(prepared, label_col = "cluster") {
  cluster_summaries <- list()

  for (city_id in names(prepared)) {
    if (is.null(prepared[[city_id]]$city_clean)) stop("[build_cluster_summaries] Missing city_clean for city: ", city_id)
    city_dt <- prepared[[city_id]]$city_clean
    if (!(label_col %in% names(city_dt))) stop("[build_cluster_summaries] Missing label column for city: ", city_id)

    col_xo <- pick_col(city_dt, c("x_o", "origin_lon"))
    col_yo <- pick_col(city_dt, c("y_o", "origin_lat"))
    col_dx <- pick_col(city_dt, c("dx", "ux", "delta_x", "d_lon"))
    col_dy <- pick_col(city_dt, c("dy", "uy", "delta_y", "d_lat"))

    has_len <- "len" %in% names(city_dt)
    has_duration <- "duration" %in% names(city_dt)

    cluster_summary <- city_dt %>%
      dplyr::filter(.data[[label_col]] != 0) %>%
      dplyr::group_by(.data[[label_col]]) %>%
      dplyr::summarise(
        size = dplyr::n(),
        mean_x_o = if (!is.na(col_xo)) mean(.data[[col_xo]], na.rm = TRUE) else NA_real_,
        mean_y_o = if (!is.na(col_yo)) mean(.data[[col_yo]], na.rm = TRUE) else NA_real_,
        mean_dx = if (!is.na(col_dx)) mean(.data[[col_dx]], na.rm = TRUE) else NA_real_,
        mean_dy = if (!is.na(col_dy)) mean(.data[[col_dy]], na.rm = TRUE) else NA_real_,
        mean_len = if (has_len) mean(.data[["len"]], na.rm = TRUE) else NA_real_,
        mean_duration = if (has_duration) mean(.data[["duration"]], na.rm = TRUE) else NA_real_,
        .groups = "drop"
      ) %>%
      dplyr::arrange(dplyr::desc(.data$size))

    names(cluster_summary)[1] <- "cluster"

    cluster_summaries[[city_id]] <- cluster_summary

    cat("\n---", toupper(city_id), "---\n")
    if (is.na(col_xo) || is.na(col_yo) || is.na(col_dx) || is.na(col_dy)) {
      cat("Note: some expected coordinate/direction columns missing; filled with NA where needed.\n")
    }
    if (!has_len || !has_duration) {
      cat("Note: len/duration missing; summary columns filled with NA where needed.\n")
    }
    print(utils::head(cluster_summary, 15))
  }

  cluster_summaries
}

build_moo_clusters <- function(cluster_summaries, prepared) {
  moo_clusters <- list()

  for (city_id in names(cluster_summaries)) {
    cs <- cluster_summaries[[city_id]]
    dominant_cluster <- if (nrow(cs) > 0) cs$cluster[which.max(cs$size)] else NA_integer_

    moo_clusters[[city_id]] <- cs
    if (is.null(prepared[[city_id]]$X_scaled)) stop("[build_moo_clusters] Missing X_scaled for city: ", city_id)
    attr(moo_clusters[[city_id]], "N_total") <- nrow(prepared[[city_id]]$X_scaled)

    cat(
      city_id,
      "| dominant cluster:", if (is.na(dominant_cluster)) "NA" else dominant_cluster,
      "| clusters kept:", nrow(cs),
      "| N_total =", attr(moo_clusters[[city_id]], "N_total"),
      "\n"
    )
  }

  moo_clusters
}

print_cluster_count_context <- function(moo_clusters) {
  for (city_id in names(moo_clusters)) {
    K <- nrow(moo_clusters[[city_id]])
    cat("CITY:", toupper(city_id), "| K =", K, "\n")
  }
  invisible(TRUE)
}

build_trip_dirs_full <- function(prepared) {
  trip_dirs_full <- lapply(names(prepared), function(city_id) {
    if (!is.null(prepared[[city_id]]$X_scaled)) {
      Xs <- prepared[[city_id]]$X_scaled
      cn <- colnames(Xs)
      if (!is.null(cn) && all(c("dx", "dy") %in% cn)) {
        dx <- as.numeric(Xs[, "dx"])
        dy <- as.numeric(Xs[, "dy"])
      } else {
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

  for (city_id in names(prepared)) {
    stopifnot(nrow(prepared[[city_id]]$X_scaled) == nrow(trip_dirs_full[[city_id]]))
  }

  trip_dirs_full
}

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

    mx <- mean(ux[idx])
    my <- mean(uy[idx])
    nm <- sqrt(mx^2 + my^2)
    if (!is.finite(nm) || nm < 1e-12) next
    mx <- mx / nm
    my <- my / nm

    cos_i <- ux[idx] * mx + uy[idx] * my
    cos_i <- pmax(-1, pmin(1, cos_i))

    coh_c <- mean(cos_i)
    num <- num + n * coh_c
    den <- den + n
  }

  if (den == 0) return(0)

  coh <- num / den
  if (map_to_01) (coh + 1) / 2 else coh
}

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
  if (nrow(u_trips) != N_total) stop("eval_solution: length(labels) must match nrow(u_trips)")

  clustered <- !is.na(labels) & labels != noise_label
  n_clustered <- sum(clustered)
  coverage <- n_clustered / N_total

  K <- if (n_clustered == 0) 0L else length(unique(labels[clustered]))

  coh <- dir_cohesion_within_weighted(
    labels = labels,
    u = u_trips,
    noise_label = noise_label,
    min_cluster_size = min_cluster_size,
    map_to_01 = TRUE
  )

  f1 <- if (coverage == 0) 0 else -coverage
  f2 <- if (coverage == 0) 0 else -coh

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

run_eval_sanity_checks <- function(trip_dirs_full, seed = 1) {
  set.seed(seed)

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

  invisible(TRUE)
}

run_eval_random_checks <- function(trip_dirs_full, seed = 1, draws = 5) {
  set.seed(seed)

  for (city_id in names(trip_dirs_full)) {
    u_trips <- trip_dirs_full[[city_id]]
    N <- nrow(u_trips)

    cat("\nCITY:", toupper(city_id), "\n")

    for (i in seq_len(draws)) {
      clustered_share <- stats::runif(1, 0.3, 0.9)
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

  invisible(TRUE)
}

# -------------------------
# 7) NSGA-II search + FULL-data validation (Steps 14-15)
# -------------------------

snap_to_allowed <- function(v, allowed) allowed[which.min(abs(allowed - v))]

snap_eps <- function(eps_raw, eps_min, eps_max, eps_step, eps_digits = 3) {
  eps_raw <- max(eps_min, min(eps_raw, eps_max))
  eps <- eps_min + round((eps_raw - eps_min) / eps_step) * eps_step
  eps <- round(eps, eps_digits)
  max(eps_min, min(eps, eps_max))
}

build_opt_subsample <- function(prepared, trip_dirs_full, n_sub = 10000L, seed = 1) {
  set.seed(seed)

  sub_idx <- lapply(names(prepared), function(city_id) {
    n <- nrow(prepared[[city_id]]$X_scaled)
    if (n > n_sub) sample.int(n, n_sub) else seq_len(n)
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

  list(sub_idx = sub_idx, data_for_opt = data_for_opt, trip_dirs_for_opt = trip_dirs_for_opt)
}

eval_dbscan_params <- function(eps, minPts, X_scaled, u_trips,
                               noise_label = 0L, min_cluster_size = 2L) {
  fit <- dbscan::dbscan(X_scaled, eps = eps, minPts = minPts)
  cl <- as.integer(fit$cluster)
  n_total <- length(cl)
  if (n_total == 0) return(c(f1 = 1, f2 = 1))

  clustered <- (cl != noise_label) & !is.na(cl)
  n_clustered <- sum(clustered)
  coverage <- n_clustered / n_total
  if (n_clustered == 0) return(c(f1 = 1, f2 = 1))

  coh <- dir_cohesion_within_weighted(cl, u_trips, noise_label, min_cluster_size, map_to_01 = TRUE)
  c(f1 = -coverage, f2 = -coh)
}

make_objfun_city_cached_dbscan <- function(city_id, data_for_opt, trip_dirs_for_opt,
                                           eps_min, eps_max, eps_step,
                                           minPts_allowed, eps_digits = 3,
                                           noise_label = 0L, min_cluster_size = 2L) {
  X_scaled <- data_for_opt[[city_id]]$X_scaled
  u_trips <- trip_dirs_for_opt[[city_id]]
  stopifnot(nrow(X_scaled) == nrow(u_trips))

  cache <- new.env(parent = emptyenv())

  function(x) {
    eps <- snap_eps(x[1], eps_min, eps_max, eps_step, eps_digits)
    mp <- snap_to_allowed(as.integer(round(x[2])), minPts_allowed)

    key <- paste0(sprintf(paste0("%.", eps_digits, "f"), eps), "_", mp)
    if (exists(key, envir = cache, inherits = FALSE)) return(get(key, envir = cache, inherits = FALSE))

    val <- eval_dbscan_params(eps, mp, X_scaled, u_trips, noise_label, min_cluster_size)
    assign(key, val, envir = cache)
    val
  }
}

run_nsga_dbscan <- function(prepared, data_for_opt, trip_dirs_for_opt,
                            eps_min = 0.08, eps_max = 0.45, eps_step = 0.02, eps_digits = 3,
                            minPts_allowed = c(6L, 8L, 10L, 12L, 15L, 20L, 25L, 30L, 40L, 50L, 60L, 70L),
                            popsize = 100, generations = 50, seed_base = 1000) {
  objfun_dbscan <- list()
  nsga_dbscan <- list()

  for (i in seq_along(names(prepared))) {
    city_id <- names(prepared)[i]
    objfun_dbscan[[city_id]] <- make_objfun_city_cached_dbscan(
      city_id, data_for_opt, trip_dirs_for_opt,
      eps_min, eps_max, eps_step, minPts_allowed, eps_digits
    )

    set.seed(seed_base + i)
    nsga_dbscan[[city_id]] <- mco::nsga2(
      fn = objfun_dbscan[[city_id]],
      idim = 2,
      odim = 2,
      lower.bounds = c(eps_min, min(minPts_allowed)),
      upper.bounds = c(eps_max, max(minPts_allowed)),
      popsize = popsize,
      generations = generations
    )
  }

  list(objfun_dbscan = objfun_dbscan, nsga_dbscan = nsga_dbscan,
       eps_min = eps_min, eps_max = eps_max, eps_step = eps_step,
       eps_digits = eps_digits, minPts_allowed = minPts_allowed)
}

build_pareto_dbscan <- function(nsga_dbscan, eps_min, eps_max, eps_step, eps_digits, minPts_allowed) {
  pareto_dbscan <- list()

  for (city_id in names(nsga_dbscan)) {
    res <- nsga_dbscan[[city_id]]
    eps_snap <- vapply(res$par[, 1], snap_eps, numeric(1), eps_min = eps_min, eps_max = eps_max,
                       eps_step = eps_step, eps_digits = eps_digits)
    mp_snap <- vapply(as.integer(round(res$par[, 2])), snap_to_allowed, integer(1), allowed = minPts_allowed)

    df <- data.frame(
      city = city_id,
      eps = as.numeric(eps_snap),
      minPts = as.integer(mp_snap),
      coverage = -res$value[, 1],
      dir_cohesion = -res$value[, 2]
    )

    pareto_dbscan[[city_id]] <- df %>%
      dplyr::filter(is.finite(.data$coverage), is.finite(.data$dir_cohesion), .data$coverage >= 0, .data$dir_cohesion >= 0) %>%
      dplyr::distinct(.data$eps, .data$minPts, .keep_all = TRUE)
  }

  pareto_dbscan
}

select_knee_points <- function(pareto_list) {
  knee <- list()
  for (city_id in names(pareto_list)) {
    df <- pareto_list[[city_id]]
    if (is.null(df) || nrow(df) == 0) {
      knee[[city_id]] <- NULL
      next
    }

    cov_n <- (df$coverage - min(df$coverage)) / (max(df$coverage) - min(df$coverage) + 1e-9)
    coh_n <- (df$dir_cohesion - min(df$dir_cohesion)) / (max(df$dir_cohesion) - min(df$dir_cohesion) + 1e-9)
    d <- sqrt((1 - cov_n)^2 + (1 - coh_n)^2)
    knee[[city_id]] <- df[which.min(d), , drop = FALSE]
  }
  knee
}

plot_pareto <- function(df, knee = NULL, expert = NULL, city_name, title_prefix = "Pareto front (subsample)") {
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$coverage, y = .data$dir_cohesion)) +
    ggplot2::geom_point(size = 2, alpha = 0.8) +
    ggplot2::labs(
      title = paste(title_prefix, "—", toupper(city_name)),
      x = "Coverage (clustered share)",
      y = "Directional cohesion (within-cluster, size-weighted)"
    ) +
    ggplot2::theme_minimal()

  if (!is.null(knee) && is.data.frame(knee) && nrow(knee) > 0) {
    p <- p + ggplot2::geom_point(data = knee, ggplot2::aes(x = .data$coverage, y = .data$dir_cohesion), color = "red", size = 4)
  }
  if (!is.null(expert) && is.data.frame(expert) && nrow(expert) > 0) {
    p <- p + ggplot2::geom_point(data = expert, ggplot2::aes(x = .data$coverage, y = .data$dir_cohesion), shape = 4, color = "black", size = 4, stroke = 1.2)
  }

  p
}

make_objfun_city_cached_snn <- function(city_id, data_for_opt, trip_dirs_for_opt,
                                        k_allowed, thr_allowed, minPts_allowed,
                                        noise_label = 0L, min_cluster_size = 2L,
                                        jp = FALSE, sort = FALSE, borderPoints = TRUE) {
  X_scaled <- data_for_opt[[city_id]]$X_scaled
  u_trips <- trip_dirs_for_opt[[city_id]]
  stopifnot(nrow(X_scaled) == nrow(u_trips))

  cache_val <- new.env(parent = emptyenv())
  cache_snn <- new.env(parent = emptyenv())

  function(x) {
    k <- snap_to_allowed(as.integer(round(x[1])), k_allowed)
    th <- snap_to_allowed(as.integer(round(x[2])), thr_allowed)
    mp <- snap_to_allowed(as.integer(round(x[3])), minPts_allowed)

    th <- max(1L, min(th, k - 1L))
    mp <- max(2L, min(mp, k))

    key <- paste(k, th, mp, sep = "_")
    if (exists(key, envir = cache_val, inherits = FALSE)) return(get(key, envir = cache_val, inherits = FALSE))

    k_key <- as.character(k)
    if (exists(k_key, envir = cache_snn, inherits = FALSE)) {
      snn_obj <- get(k_key, envir = cache_snn, inherits = FALSE)
    } else {
      snn_obj <- dbscan::sNN(X_scaled, k = k, sort = sort, jp = jp)
      assign(k_key, snn_obj, envir = cache_snn)
    }

    fit <- dbscan::sNNclust(x = snn_obj, k = k, eps = th, minPts = mp, borderPoints = borderPoints)
    cl <- as.integer(fit$cluster)

    n_total <- length(cl)
    if (n_total == 0) {
      val <- c(f1 = 1, f2 = 1)
    } else {
      clustered <- (cl != noise_label) & !is.na(cl)
      n_clustered <- sum(clustered)
      if (n_clustered == 0) {
        val <- c(f1 = 1, f2 = 1)
      } else {
        coverage <- n_clustered / n_total
        coh <- dir_cohesion_within_weighted(cl, u_trips, noise_label, min_cluster_size, map_to_01 = TRUE)
        val <- c(f1 = -coverage, f2 = -coh)
      }
    }

    assign(key, val, envir = cache_val)
    val
  }
}

run_nsga_snn <- function(prepared, data_for_opt, trip_dirs_for_opt,
                         k_allowed = c(15L, 20L, 25L, 30L, 35L),
                         thr_allowed = 4L:14L,
                         minPts_allowed = c(6L, 8L, 10L, 12L, 15L, 18L, 20L),
                         popsize = 100, generations = 50, seed_base = 2000) {
  objfun_snn <- list()
  nsga_snn <- list()

  lb <- as.numeric(c(min(k_allowed), min(thr_allowed), min(minPts_allowed)))
  ub <- as.numeric(c(max(k_allowed), max(thr_allowed), max(minPts_allowed)))

  for (i in seq_along(names(prepared))) {
    city_id <- names(prepared)[i]
    objfun_snn[[city_id]] <- make_objfun_city_cached_snn(
      city_id, data_for_opt, trip_dirs_for_opt, k_allowed, thr_allowed, minPts_allowed
    )

    set.seed(seed_base + i)
    nsga_snn[[city_id]] <- mco::nsga2(
      fn = objfun_snn[[city_id]], idim = 3, odim = 2,
      lower.bounds = lb, upper.bounds = ub,
      popsize = popsize, generations = generations
    )
  }

  list(objfun_snn = objfun_snn, nsga_snn = nsga_snn,
       k_allowed = k_allowed, thr_allowed = thr_allowed, minPts_allowed = minPts_allowed)
}

build_pareto_snn <- function(nsga_snn, k_allowed, thr_allowed, minPts_allowed) {
  pareto_snn <- list()

  for (city_id in names(nsga_snn)) {
    res <- nsga_snn[[city_id]]

    k <- vapply(as.integer(round(res$par[, 1])), snap_to_allowed, integer(1), allowed = k_allowed)
    th <- vapply(as.integer(round(res$par[, 2])), snap_to_allowed, integer(1), allowed = thr_allowed)
    mp <- vapply(as.integer(round(res$par[, 3])), snap_to_allowed, integer(1), allowed = minPts_allowed)

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

    pareto_snn[[city_id]] <- df %>%
      dplyr::filter(is.finite(.data$coverage), is.finite(.data$dir_cohesion), .data$coverage >= 0, .data$dir_cohesion >= 0) %>%
      dplyr::distinct(.data$k, .data$snn_threshold, .data$minPts, .keep_all = TRUE)
  }

  pareto_snn
}

print_all_pareto <- function(pareto_list, algo = "DBSCAN") {
  options(max.print = 1e6)
  for (city_id in names(pareto_list)) {
    cat("\n========================================\n")
    cat("[", algo, "] CITY:", toupper(city_id), "\n", sep = "")
    cat("N Pareto points:", nrow(pareto_list[[city_id]]), "\n")
    cat("========================================\n\n")

    df <- pareto_list[[city_id]]
    if (is.null(df) || nrow(df) == 0) {
      cat("No valid Pareto rows.\n")
      next
    }

    print(df[order(df$coverage), ], row.names = FALSE)
  }
}

select_configs_dbscan <- function(df, knee_row, n_each = 2) {
  if (is.null(df) || nrow(df) == 0) return(df)

  if (is.null(knee_row) || !is.data.frame(knee_row) || nrow(knee_row) == 0) {
    return(df[order(-df$coverage, -df$dir_cohesion), ][1:min(5, nrow(df)), , drop = FALSE])
  }

  knee_match <- df %>% dplyr::filter(abs(.data$eps - knee_row$eps[1]) < 1e-6, .data$minPts == knee_row$minPts[1])
  if (nrow(knee_match) == 0) {
    knee_match <- df %>% dplyr::filter(abs(.data$coverage - knee_row$coverage[1]) < 1e-9,
                                       abs(.data$dir_cohesion - knee_row$dir_cohesion[1]) < 1e-9)
  }
  if (nrow(knee_match) == 0) knee_match <- df[order(-df$coverage, -df$dir_cohesion), ][1, , drop = FALSE] else knee_match <- knee_match[1, , drop = FALSE]

  higher_cov <- df %>% dplyr::filter(.data$coverage > knee_match$coverage) %>% dplyr::arrange(.data$coverage, dplyr::desc(.data$dir_cohesion)) %>% utils::head(n_each)
  higher_coh <- df %>% dplyr::filter(.data$dir_cohesion > knee_match$dir_cohesion) %>% dplyr::arrange(.data$dir_cohesion, dplyr::desc(.data$coverage)) %>% utils::head(n_each)

  dplyr::bind_rows(knee_match, higher_cov, higher_coh) %>% dplyr::distinct(.data$eps, .data$minPts, .keep_all = TRUE) %>% utils::head(1 + 2 * n_each)
}

full_eval_one_dbscan <- function(city_id, eps, minPts, prepared_full,
                                 noise_label = 0L, min_cluster_size = 2L) {
  X_full <- prepared_full[[city_id]]$X_scaled
  city_clean <- prepared_full[[city_id]]$city_clean
  stopifnot(nrow(X_full) == nrow(city_clean))

  dx <- as.numeric(city_clean$dx)
  dy <- as.numeric(city_clean$dy)
  norm <- sqrt(dx^2 + dy^2)
  norm[norm == 0] <- NA_real_
  u_full <- cbind(ux = dx / norm, uy = dy / norm)

  fit <- dbscan::dbscan(X_full, eps = eps, minPts = minPts)
  cl <- as.integer(fit$cluster)

  n_total <- length(cl)
  clustered <- (cl != noise_label) & !is.na(cl)
  n_clustered <- sum(clustered)
  coverage <- if (n_total > 0) n_clustered / n_total else 0
  K <- if (n_clustered == 0) 0L else length(unique(cl[clustered]))

  coh <- if (n_clustered == 0) 0 else dir_cohesion_within_weighted(cl, u_full, noise_label, min_cluster_size, map_to_01 = TRUE)

  data.frame(
    city = city_id,
    eps = as.numeric(eps),
    minPts = as.integer(minPts),
    coverage = as.numeric(coverage),
    clustered_pct = 100 * coverage,
    dir_cohesion = as.numeric(coh),
    n_clusters = as.integer(K),
    noise_pct = if (n_total > 0) 100 * sum(cl == noise_label, na.rm = TRUE) / n_total else NA_real_,
    n_total = as.integer(n_total)
  )
}

validate_pareto_dbscan_full <- function(prepared, pareto_dbscan, knee_dbscan, dbscan_params = NULL, n_each = 2) {
  pareto_unique_dbscan <- lapply(names(pareto_dbscan), function(city_id) {
    df <- pareto_dbscan[[city_id]]
    if (is.null(df) || nrow(df) == 0) return(df)
    df %>%
      dplyr::select(.data$city, .data$eps, .data$minPts, .data$coverage, .data$dir_cohesion) %>%
      dplyr::mutate(eps = as.numeric(.data$eps), minPts = as.integer(.data$minPts),
                    coverage = as.numeric(.data$coverage), dir_cohesion = as.numeric(.data$dir_cohesion)) %>%
      dplyr::filter(is.finite(.data$coverage), is.finite(.data$dir_cohesion), .data$coverage >= 0, .data$dir_cohesion >= 0) %>%
      dplyr::distinct(.data$eps, .data$minPts, .keep_all = TRUE)
  })
  names(pareto_unique_dbscan) <- names(pareto_dbscan)

  full_validated_dbscan <- list()
  for (city_id in names(pareto_unique_dbscan)) {
    df <- pareto_unique_dbscan[[city_id]]
    knee <- knee_dbscan[[city_id]]
    if (is.null(df) || nrow(df) == 0) next

    chosen <- select_configs_dbscan(df, knee, n_each = n_each)
    pareto_out <- do.call(rbind, lapply(seq_len(nrow(chosen)), function(i) {
      full_eval_one_dbscan(city_id, chosen$eps[i], chosen$minPts[i], prepared)
    }))
    pareto_out$config <- "pareto_local"

    expert_out <- NULL
    if (!is.null(dbscan_params) && !is.null(dbscan_params[[city_id]])) {
      ex <- dbscan_params[[city_id]]
      expert_out <- full_eval_one_dbscan(city_id, ex$eps, ex$minPts, prepared)
      expert_out$config <- "expert"
    }

    full_validated_dbscan[[city_id]] <- dplyr::bind_rows(pareto_out, expert_out)
  }

  full_df <- if (length(full_validated_dbscan) > 0) do.call(rbind, full_validated_dbscan) else data.frame()

  knee_full <- NULL
  compare_df <- NULL
  if (nrow(full_df) > 0) {
    knee_full <- full_df %>%
      dplyr::filter(.data$config == "pareto_local") %>%
      dplyr::group_by(.data$city) %>%
      dplyr::mutate(
        cov_n = (.data$coverage - min(.data$coverage)) / (max(.data$coverage) - min(.data$coverage) + 1e-9),
        coh_n = (.data$dir_cohesion - min(.data$dir_cohesion)) / (max(.data$dir_cohesion) - min(.data$dir_cohesion) + 1e-9),
        dist_ideal = sqrt((1 - .data$cov_n)^2 + (1 - .data$coh_n)^2)
      ) %>%
      dplyr::slice_min(.data$dist_ideal, n = 1, with_ties = FALSE) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(config = "knee_full") %>%
      dplyr::select(.data$city, .data$config, .data$eps, .data$minPts, .data$coverage,
                    .data$clustered_pct, .data$dir_cohesion, .data$n_clusters, .data$noise_pct)

    compare_df <- full_df %>%
      dplyr::filter(.data$config %in% c("expert")) %>%
      dplyr::select(.data$city, .data$config, .data$eps, .data$minPts, .data$coverage,
                    .data$clustered_pct, .data$dir_cohesion, .data$n_clusters, .data$noise_pct) %>%
      dplyr::bind_rows(knee_full) %>%
      dplyr::arrange(.data$city, factor(.data$config, levels = c("expert", "knee_full")))
  }

  list(full_validated_dbscan_df = full_df, knee_full_dbscan = knee_full, compare_dbscan = compare_df)
}

select_configs_snn <- function(df, knee_row, n_each = 2) {
  if (is.null(df) || nrow(df) == 0) return(df)

  if (is.null(knee_row) || !is.data.frame(knee_row) || nrow(knee_row) == 0) {
    return(df[order(-df$coverage, -df$dir_cohesion), ][1:min(5, nrow(df)), , drop = FALSE])
  }

  knee_match <- df %>% dplyr::filter(.data$k == knee_row$k[1], .data$snn_threshold == knee_row$snn_threshold[1], .data$minPts == knee_row$minPts[1])
  if (nrow(knee_match) == 0) {
    knee_match <- df %>% dplyr::filter(abs(.data$coverage - knee_row$coverage[1]) < 1e-9,
                                       abs(.data$dir_cohesion - knee_row$dir_cohesion[1]) < 1e-9)
  }
  if (nrow(knee_match) == 0) knee_match <- df[order(-df$coverage, -df$dir_cohesion), ][1, , drop = FALSE] else knee_match <- knee_match[1, , drop = FALSE]

  higher_cov <- df %>% dplyr::filter(.data$coverage > knee_match$coverage) %>% dplyr::arrange(.data$coverage, dplyr::desc(.data$dir_cohesion)) %>% utils::head(n_each)
  higher_coh <- df %>% dplyr::filter(.data$dir_cohesion > knee_match$dir_cohesion) %>% dplyr::arrange(.data$dir_cohesion, dplyr::desc(.data$coverage)) %>% utils::head(n_each)

  dplyr::bind_rows(knee_match, higher_cov, higher_coh) %>% dplyr::distinct(.data$k, .data$snn_threshold, .data$minPts, .keep_all = TRUE) %>% utils::head(1 + 2 * n_each)
}

full_eval_one_snn <- function(city_id, k, snn_threshold, minPts, prepared_full,
                              noise_label = 0L, min_cluster_size = 2L,
                              jp = FALSE, sort = FALSE, borderPoints = TRUE) {
  X_full <- prepared_full[[city_id]]$X_scaled
  city_clean <- prepared_full[[city_id]]$city_clean
  stopifnot(nrow(X_full) == nrow(city_clean))

  dx <- as.numeric(city_clean$dx)
  dy <- as.numeric(city_clean$dy)
  norm <- sqrt(dx^2 + dy^2)
  norm[norm == 0] <- NA_real_
  u_full <- cbind(ux = dx / norm, uy = dy / norm)

  snn_obj <- dbscan::sNN(X_full, k = k, sort = sort, jp = jp)
  fit <- dbscan::sNNclust(x = snn_obj, k = k, eps = snn_threshold, minPts = minPts, borderPoints = borderPoints)
  cl <- as.integer(fit$cluster)

  n_total <- length(cl)
  clustered <- (cl != noise_label) & !is.na(cl)
  n_clustered <- sum(clustered)

  coverage <- if (n_total > 0) n_clustered / n_total else 0
  K <- if (n_clustered == 0) 0L else length(unique(cl[clustered]))

  coh <- if (n_clustered == 0) 0 else dir_cohesion_within_weighted(cl, u_full, noise_label, min_cluster_size, map_to_01 = TRUE)

  data.frame(
    city = city_id,
    k = as.integer(k),
    snn_threshold = as.integer(snn_threshold),
    minPts = as.integer(minPts),
    coverage = as.numeric(coverage),
    clustered_pct = 100 * coverage,
    dir_cohesion = as.numeric(coh),
    n_clusters = as.integer(K),
    noise_pct = if (n_total > 0) 100 * sum(cl == noise_label, na.rm = TRUE) / n_total else NA_real_,
    n_total = as.integer(n_total)
  )
}

validate_pareto_snn_full <- function(prepared, pareto_snn, knee_snn, snn_params = NULL, n_each = 2) {
  pareto_unique_snn <- lapply(names(pareto_snn), function(city_id) {
    df <- pareto_snn[[city_id]]
    if (is.null(df) || nrow(df) == 0) return(df)
    df %>%
      dplyr::select(.data$city, .data$k, .data$snn_threshold, .data$minPts, .data$coverage, .data$dir_cohesion) %>%
      dplyr::mutate(k = as.integer(.data$k), snn_threshold = as.integer(.data$snn_threshold),
                    minPts = as.integer(.data$minPts), coverage = as.numeric(.data$coverage),
                    dir_cohesion = as.numeric(.data$dir_cohesion)) %>%
      dplyr::filter(is.finite(.data$coverage), is.finite(.data$dir_cohesion), .data$coverage >= 0, .data$dir_cohesion >= 0) %>%
      dplyr::distinct(.data$k, .data$snn_threshold, .data$minPts, .keep_all = TRUE)
  })
  names(pareto_unique_snn) <- names(pareto_snn)

  full_validated_snn <- list()
  for (city_id in names(pareto_unique_snn)) {
    df <- pareto_unique_snn[[city_id]]
    knee <- knee_snn[[city_id]]
    if (is.null(df) || nrow(df) == 0) next

    chosen <- select_configs_snn(df, knee, n_each = n_each)
    pareto_out <- do.call(rbind, lapply(seq_len(nrow(chosen)), function(i) {
      full_eval_one_snn(city_id, chosen$k[i], chosen$snn_threshold[i], chosen$minPts[i], prepared)
    }))
    pareto_out$config <- "pareto_local"

    expert_out <- NULL
    if (!is.null(snn_params) && !is.null(snn_params[[city_id]])) {
      ex <- snn_params[[city_id]]
      expert_out <- full_eval_one_snn(city_id, ex$k, ex$snn_threshold, ex$minPts, prepared)
      expert_out$config <- "expert"
    }

    full_validated_snn[[city_id]] <- dplyr::bind_rows(pareto_out, expert_out)
  }

  full_df <- if (length(full_validated_snn) > 0) do.call(rbind, full_validated_snn) else data.frame()

  knee_full <- NULL
  compare_df <- NULL
  if (nrow(full_df) > 0) {
    knee_full <- full_df %>%
      dplyr::filter(.data$config == "pareto_local") %>%
      dplyr::group_by(.data$city) %>%
      dplyr::mutate(
        cov_n = (.data$coverage - min(.data$coverage)) / (max(.data$coverage) - min(.data$coverage) + 1e-9),
        coh_n = (.data$dir_cohesion - min(.data$dir_cohesion)) / (max(.data$dir_cohesion) - min(.data$dir_cohesion) + 1e-9),
        dist_ideal = sqrt((1 - .data$cov_n)^2 + (1 - .data$coh_n)^2)
      ) %>%
      dplyr::slice_min(.data$dist_ideal, n = 1, with_ties = FALSE) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(config = "knee_full") %>%
      dplyr::select(.data$city, .data$config, .data$k, .data$snn_threshold, .data$minPts,
                    .data$coverage, .data$clustered_pct, .data$dir_cohesion, .data$n_clusters, .data$noise_pct)

    compare_df <- full_df %>%
      dplyr::filter(.data$config %in% c("expert")) %>%
      dplyr::select(.data$city, .data$config, .data$k, .data$snn_threshold, .data$minPts,
                    .data$coverage, .data$clustered_pct, .data$dir_cohesion, .data$n_clusters, .data$noise_pct) %>%
      dplyr::bind_rows(knee_full) %>%
      dplyr::arrange(.data$city, factor(.data$config, levels = c("expert", "knee_full")))
  }

  list(full_validated_snn_df = full_df, knee_full_snn = knee_full, compare_snn = compare_df)
}

# -------------------------
# 8) Hypervolume analysis + null-test + flow plotting utilities (Step 16)
# -------------------------

nondominated_min_points <- function(F) {
  n <- nrow(F)
  keep <- rep(TRUE, n)
  for (i in seq_len(n)) {
    if (!keep[i]) next
    for (j in seq_len(n)) {
      if (i == j || !keep[i]) next
      if (all(F[j, ] <= F[i, ]) && any(F[j, ] < F[i, ])) keep[i] <- FALSE
    }
  }
  F[keep, , drop = FALSE]
}

hv2d_min_exact <- function(F, ref) {
  if (nrow(F) == 0) return(0)
  F <- nondominated_min_points(F)
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

hv2d_contrib_min_exact <- function(F, ref) {
  hv_all <- hv2d_min_exact(F, ref)
  n <- nrow(F)
  vapply(seq_len(n), function(i) hv_all - hv2d_min_exact(F[-i, , drop = FALSE], ref), numeric(1))
}

safe_norm_vec <- function(z, zmin, zmax) {
  if (!is.finite(zmin) || !is.finite(zmax) || isTRUE(all.equal(zmax, zmin))) return(rep(0, length(z)))
  (z - zmin) / (zmax - zmin)
}

gini_coeff <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  x <- sort(pmax(0, x))
  s <- sum(x)
  if (s <= 0) return(0)
  n <- length(x)
  (2 * sum(seq_len(n) * x) / (n * s)) - (n + 1) / n
}

compute_hv_results <- function(pareto_by_group, ref_min = c(1.05, 1.05),
                               group_col = "algorithm", city_col = "city") {
  all_pts <- do.call(rbind, lapply(names(pareto_by_group), function(grp) {
    pt_list <- pareto_by_group[[grp]]
    if (is.null(pt_list) || length(pt_list) == 0) return(NULL)

    do.call(rbind, lapply(names(pt_list), function(city_id) {
      df <- pt_list[[city_id]]
      if (is.null(df) || nrow(df) == 0) return(NULL)
      df2 <- df[, c("coverage", "dir_cohesion")]
      df2 <- df2[is.finite(df2$coverage) & is.finite(df2$dir_cohesion), , drop = FALSE]
      if (nrow(df2) == 0) return(NULL)
      df2[[group_col]] <- grp
      df2[[city_col]] <- city_id
      df2
    }))
  }))

  if (is.null(all_pts) || nrow(all_pts) == 0) {
    return(list(hv_results = data.frame(), hv_summary_df = data.frame(), global_pool = all_pts))
  }

  cov_min <- min(all_pts$coverage, na.rm = TRUE)
  cov_max <- max(all_pts$coverage, na.rm = TRUE)
  coh_min <- min(all_pts$dir_cohesion, na.rm = TRUE)
  coh_max <- max(all_pts$dir_cohesion, na.rm = TRUE)

  hv_max <- prod(ref_min - c(0, 0))

  hv_results <- do.call(rbind, lapply(names(pareto_by_group), function(grp) {
    pt_list <- pareto_by_group[[grp]]
    if (is.null(pt_list) || length(pt_list) == 0) return(NULL)

    do.call(rbind, lapply(names(pt_list), function(city_id) {
      df <- pt_list[[city_id]]
      if (is.null(df) || nrow(df) == 0) {
        return(data.frame(
          group = grp,
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
          group = grp,
          city = city_id,
          n_points = 0,
          hv_raw = NA_real_,
          hv_unit = NA_real_,
          contrib_gini = NA_real_,
          contrib_top10_share = NA_real_,
          contrib_cv = NA_real_
        ))
      }

      cov_n <- safe_norm_vec(df$coverage, cov_min, cov_max)
      coh_n <- safe_norm_vec(df$dir_cohesion, coh_min, coh_max)
      Fmin <- cbind(f1 = 1 - cov_n, f2 = 1 - coh_n)

      hv_raw <- hv2d_min_exact(Fmin, ref_min)
      hv_unit <- hv_raw / hv_max

      contrib <- hv2d_contrib_min_exact(Fmin, ref_min)
      contrib <- contrib[is.finite(contrib)]

      top10_share <- if (length(contrib) == 0 || sum(contrib) == 0) NA_real_ else {
        n_top <- max(1, floor(0.1 * length(contrib)))
        sum(sort(contrib, decreasing = TRUE)[seq_len(n_top)]) / sum(contrib)
      }

      data.frame(
        group = grp,
        city = city_id,
        n_points = nrow(df),
        hv_raw = hv_raw,
        hv_unit = hv_unit,
        contrib_gini = gini_coeff(contrib),
        contrib_top10_share = top10_share,
        contrib_cv = if (length(contrib) > 1 && mean(contrib) > 0) stats::sd(contrib) / mean(contrib) else NA_real_
      )
    }))
  }))

  hv_results <- hv_results %>% dplyr::arrange(.data$group, dplyr::desc(.data$hv_unit), .data$city)
  hv_summary_df <- hv_results %>% dplyr::transmute(group = .data$group, city = .data$city, hv_unit = .data$hv_unit, gini = .data$contrib_gini)

  list(hv_results = hv_results, hv_summary_df = hv_summary_df, global_pool = all_pts)
}

plot_hv_summaries <- function(hv_results, hv_summary_df) {
  if (nrow(hv_summary_df) == 0) return(invisible(NULL))

  p1 <- ggplot2::ggplot(hv_summary_df, ggplot2::aes(x = .data$city, y = .data$hv_unit, fill = .data$group)) +
    ggplot2::geom_col(position = "dodge", width = 0.7) +
    ggplot2::labs(title = "Unit-scaled hypervolume by city and group", x = "City", y = "Hypervolume (0–1)") +
    ggplot2::theme_minimal()

  p2 <- ggplot2::ggplot(hv_summary_df, ggplot2::aes(x = .data$city, y = .data$gini, fill = .data$group)) +
    ggplot2::geom_col(position = "dodge", width = 0.7) +
    ggplot2::labs(title = "Hypervolume shape (Gini of contributions)", x = "City", y = "Gini (0=even, 1=concentrated)") +
    ggplot2::theme_minimal()

  p3 <- ggplot2::ggplot(hv_results, ggplot2::aes(x = .data$hv_unit, y = .data$contrib_gini,
                                                 label = paste0(toupper(.data$city), " ", .data$group))) +
    ggplot2::geom_point(size = 3) +
    ggplot2::geom_text(nudge_y = 0.02, size = 3) +
    ggplot2::labs(title = "Global vs shape of Pareto fronts", x = "Global hypervolume (unit-scaled)", y = "Shape (Gini of HV contributions)") +
    ggplot2::theme_minimal()

  print(p1); print(p2); print(p3)
  invisible(list(hv_plot = p1, gini_plot = p2, scatter_plot = p3))
}

make_randomized_od <- function(city_clean,
                               origin_cols = c("x_o", "y_o"),
                               dest_cols = c("x_d", "y_d"),
                               recompute_dxdy = TRUE,
                               seed = 1L) {
  need <- c(origin_cols, dest_cols)
  miss <- setdiff(need, names(city_clean))
  if (length(miss) > 0) stop("Null model missing columns: ", paste(miss, collapse = ", "))

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

build_prepared_randomized <- function(prepared_real, feature_cols_used = c("x_o", "y_o", "dx", "dy"), seed_base = 100) {
  prepared_rand <- lapply(names(prepared_real), function(city_id) {
    cc_real <- prepared_real[[city_id]]$city_clean
    stopifnot(all(c("x_o", "y_o", "x_d", "y_d") %in% names(cc_real)))

    cc_rand <- make_randomized_od(cc_real, seed = seed_base + match(city_id, names(prepared_real)))
    res <- make_X_and_clean(cc_rand, feature_cols = feature_cols_used)

    if (is.null(res$city_clean)) res$city_clean <- cc_rand
    if (is.null(res$X_scaled) && !is.null(res$X)) res$X_scaled <- res$X
    res
  })
  names(prepared_rand) <- names(prepared_real)

  for (city_id in names(prepared_real)) {
    if (nrow(prepared_real[[city_id]]$city_clean) != nrow(prepared_rand[[city_id]]$city_clean)) {
      stop("Row mismatch REAL vs RAND for city: ", city_id)
    }
  }

  prepared_rand
}

sanity_check_randomization <- function(city_real_clean, city_rand_clean,
                                       origin_cols = c("x_o", "y_o"),
                                       dest_cols = c("x_d", "y_d"),
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

  list(
    origin_same = origin_same,
    dest_pair_match_rate = match_rate,
    perm_ok_x = perm_ok_x,
    perm_ok_y = perm_ok_y,
    preview = preview
  )
}

make_trip_dirs_from_clean <- function(city_clean) {
  stopifnot(all(c("dx", "dy") %in% names(city_clean)))
  v <- as.matrix(city_clean[, c("dx", "dy")])
  nrm <- sqrt(rowSums(v^2)) + 1e-12
  u <- v / nrm
  colnames(u) <- c("ux", "uy")
  u
}

make_sub_idx <- function(prepared_in, n_sub = 10000L, seed = 1) {
  set.seed(seed)
  idx <- lapply(names(prepared_in), function(city_id) {
    n <- nrow(prepared_in[[city_id]]$X_scaled)
    if (n > n_sub) sample.int(n, n_sub) else seq_len(n)
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

run_nsga2_dbscan_variant <- function(prepared_sub, trip_dirs_sub, tag,
                                     popsize = 100, generations = 50, seed0 = 1,
                                     eps_min = 0.08, eps_max = 0.45, eps_step = 0.02,
                                     eps_digits = 3,
                                     minPts_allowed = c(6L, 8L, 10L, 12L, 15L, 20L, 25L, 30L, 40L, 50L, 60L, 70L)) {
  pareto <- list()

  for (i in seq_along(names(prepared_sub))) {
    city_id <- names(prepared_sub)[i]

    objfun <- make_objfun_city_cached_dbscan(
      city_id = city_id,
      data_for_opt = prepared_sub,
      trip_dirs_for_opt = trip_dirs_sub,
      eps_min = eps_min,
      eps_max = eps_max,
      eps_step = eps_step,
      minPts_allowed = minPts_allowed,
      eps_digits = eps_digits
    )

    set.seed(seed0 + i)
    res <- mco::nsga2(
      fn = objfun,
      idim = 2,
      odim = 2,
      lower.bounds = as.numeric(c(eps_min, min(minPts_allowed))),
      upper.bounds = as.numeric(c(eps_max, max(minPts_allowed))),
      popsize = popsize,
      generations = generations
    )

    eps_snap <- vapply(res$par[, 1], snap_eps, numeric(1), eps_min = eps_min, eps_max = eps_max,
                       eps_step = eps_step, eps_digits = eps_digits)
    mp_snap <- vapply(as.integer(round(res$par[, 2])), snap_to_allowed, integer(1), allowed = minPts_allowed)

    pareto[[city_id]] <- data.frame(
      city = city_id,
      tag = tag,
      eps = as.numeric(eps_snap),
      minPts = as.integer(mp_snap),
      coverage = -res$value[, 1],
      dir_cohesion = -res$value[, 2]
    ) %>%
      dplyr::filter(is.finite(.data$coverage), is.finite(.data$dir_cohesion), .data$coverage >= 0, .data$dir_cohesion >= 0) %>%
      dplyr::distinct(.data$eps, .data$minPts, .keep_all = TRUE)
  }

  list(pareto = pareto)
}

run_dbscan_null_test <- function(prepared_real,
                                 feature_cols_used = c("x_o", "y_o", "dx", "dy"),
                                 n_sub = 10000L,
                                 popsize = 100,
                                 generations = 50,
                                 sub_seed = 1,
                                 real_seed0 = 10,
                                 rand_seed0 = 20) {
  prepared_rand <- build_prepared_randomized(prepared_real, feature_cols_used = feature_cols_used)

  trip_dirs_real <- lapply(names(prepared_real), function(city_id) make_trip_dirs_from_clean(prepared_real[[city_id]]$city_clean))
  names(trip_dirs_real) <- names(prepared_real)
  trip_dirs_rand <- lapply(names(prepared_rand), function(city_id) make_trip_dirs_from_clean(prepared_rand[[city_id]]$city_clean))
  names(trip_dirs_rand) <- names(prepared_rand)

  sub_idx <- make_sub_idx(prepared_real, n_sub = n_sub, seed = sub_seed)
  sub_real <- apply_subsample(prepared_real, trip_dirs_real, sub_idx)
  sub_rand <- apply_subsample(prepared_rand, trip_dirs_rand, sub_idx)

  out_real <- run_nsga2_dbscan_variant(sub_real$prepared_sub, sub_real$trip_dirs_sub, tag = "real",
                                       popsize = popsize, generations = generations, seed0 = real_seed0)
  out_rand <- run_nsga2_dbscan_variant(sub_rand$prepared_sub, sub_rand$trip_dirs_sub, tag = "rand",
                                       popsize = popsize, generations = generations, seed0 = rand_seed0)

  hv <- compute_hv_results(pareto_by_group = list(real = out_real$pareto, rand = out_rand$pareto),
                           ref_min = c(1.05, 1.05), group_col = "tag", city_col = "city")

  list(prepared_rand = prepared_rand, sub_idx = sub_idx,
       out_real = out_real, out_rand = out_rand,
       hv_results_tag = hv$hv_results, hv_summary_tag = hv$hv_summary_df)
}

plot_flow_topN <- function(city_id, prepared, labels,
                           algo_name = "ALGO",
                           noise_label = 0L,
                           out_dir = "plots_flows_top3",
                           top_n_clusters = 3,
                           sample_n_segments_top = 25000,
                           sample_n_segments_bg = 30000,
                           show_background = FALSE,
                           seed = 1,
                           show_plots = TRUE) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  dt <- data.table::as.data.table(prepared[[city_id]]$city_clean)
  req <- c("x_o", "y_o", "x_d", "y_d")
  if (!all(req %in% names(dt))) stop("Missing columns in city_clean: ", paste(req, collapse = ", "))
  if (length(labels) != nrow(dt)) stop("labels length != nrow(city_clean) for ", city_id)

  dt[, cluster := as.integer(labels)]

  cl_nonnoise <- dt$cluster[!is.na(dt$cluster) & dt$cluster != noise_label]
  tab <- sort(table(cl_nonnoise), decreasing = TRUE)
  top_ids <- as.integer(names(utils::head(tab, top_n_clusters)))
  if (length(top_ids) == 0) return(NULL)

  p <- ggplot2::ggplot()

  if (isTRUE(show_background)) {
    dt_bg <- dt[!is.na(cluster) & cluster != noise_label]
    if (nrow(dt_bg) > 0) {
      set.seed(seed)
      idx_bg <- if (nrow(dt_bg) > sample_n_segments_bg) sample.int(nrow(dt_bg), sample_n_segments_bg) else seq_len(nrow(dt_bg))
      dt_bg <- dt_bg[idx_bg, ]
      p <- p + ggplot2::geom_segment(data = dt_bg,
                                     ggplot2::aes(x = .data$x_o, y = .data$y_o, xend = .data$x_d, yend = .data$y_d),
                                     alpha = 0.05, linewidth = 0.20)
    }
  }

  dt_top <- dt[cluster %in% top_ids]
  set.seed(seed)
  idx_top <- if (nrow(dt_top) > sample_n_segments_top) sample.int(nrow(dt_top), sample_n_segments_top) else seq_len(nrow(dt_top))
  dt_top <- dt_top[idx_top, ]
  dt_top[, cluster_f := factor(cluster, levels = top_ids)]

  p <- p +
    ggplot2::geom_segment(data = dt_top,
                          ggplot2::aes(x = .data$x_o, y = .data$y_o, xend = .data$x_d, yend = .data$y_d, color = .data$cluster_f),
                          alpha = 0.28, linewidth = 0.30) +
    ggplot2::coord_equal() +
    ggplot2::labs(
      title = paste0(algo_name, " — Top ", top_n_clusters, " clusters — ", toupper(city_id)),
      subtitle = if (show_background) "Colored = top clusters | Grey = all clustered flows (sampled)" else "Colored = top clusters (sampled)",
      x = "UTM Easting (m)", y = "UTM Northing (m)", color = "Cluster"
    ) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"), legend.position = "right")

  f_out <- file.path(out_dir, paste0("flows_", algo_name, "_", city_id, "_top", top_n_clusters,
                                     ifelse(show_background, "_bg", ""), ".png"))
  ggplot2::ggsave(f_out, plot = p, width = 9, height = 7.5, dpi = 300)
  if (isTRUE(show_plots)) print(p)

  list(file = f_out, top_cluster_ids = top_ids, plot = p)
}

get_dbscan_selected_params <- function(prepared_plot, preferred = "knee_full", knee_full_dbscan = NULL, dbscan_params = NULL) {
  out <- list()
  for (city_id in names(prepared_plot)) {
    if (preferred == "knee_full" && is.data.frame(knee_full_dbscan) && nrow(knee_full_dbscan) > 0) {
      row <- knee_full_dbscan %>% dplyr::filter(.data$city == city_id) %>% dplyr::slice(1)
      if (nrow(row) == 1) {
        out[[city_id]] <- list(eps = as.numeric(row$eps), minPts = as.integer(row$minPts))
        next
      }
    }

    if (!is.null(dbscan_params) && !is.null(dbscan_params[[city_id]])) {
      out[[city_id]] <- list(eps = as.numeric(dbscan_params[[city_id]]$eps), minPts = as.integer(dbscan_params[[city_id]]$minPts))
      next
    }

    stop("[DBSCAN] No selected parameters for city: ", city_id)
  }
  out
}

get_snn_selected_params <- function(prepared_plot, preferred = "knee_full", knee_full_snn = NULL, snn_params = NULL) {
  out <- list()
  for (city_id in names(prepared_plot)) {
    if (preferred == "knee_full" && is.data.frame(knee_full_snn) && nrow(knee_full_snn) > 0) {
      row <- knee_full_snn %>% dplyr::filter(.data$city == city_id) %>% dplyr::slice(1)
      if (nrow(row) == 1) {
        out[[city_id]] <- list(k = as.integer(row$k), snn_threshold = as.integer(row$snn_threshold), minPts = as.integer(row$minPts))
        next
      }
    }

    if (!is.null(snn_params) && !is.null(snn_params[[city_id]])) {
      out[[city_id]] <- list(k = as.integer(snn_params[[city_id]]$k), snn_threshold = as.integer(snn_params[[city_id]]$snn_threshold), minPts = as.integer(snn_params[[city_id]]$minPts))
      next
    }

    stop("[SNN] No selected parameters for city: ", city_id)
  }
  out
}

run_and_plot_selected_flows <- function(prepared_plot,
                                        selected_params_dbscan,
                                        selected_params_snn,
                                        out_dir = "plots_flows_top3",
                                        top_n_clusters = 3,
                                        show_background = FALSE,
                                        show_plots = TRUE) {
  out_dbscan <- list()
  out_snn <- list()

  for (city_id in names(prepared_plot)) {
    p <- selected_params_dbscan[[city_id]]
    X <- prepared_plot[[city_id]]$X_scaled
    fit <- dbscan::dbscan(X, eps = p$eps, minPts = p$minPts)
    out_dbscan[[city_id]] <- plot_flow_topN(
      city_id = city_id,
      prepared = prepared_plot,
      labels = as.integer(fit$cluster),
      algo_name = "DBSCAN",
      out_dir = out_dir,
      top_n_clusters = top_n_clusters,
      show_background = show_background,
      show_plots = show_plots
    )
  }

  for (city_id in names(prepared_plot)) {
    p <- selected_params_snn[[city_id]]
    X <- prepared_plot[[city_id]]$X_scaled
    snn_obj <- dbscan::sNN(X, k = p$k, sort = FALSE, jp = FALSE)
    fit <- dbscan::sNNclust(x = snn_obj, k = p$k, eps = p$snn_threshold, minPts = p$minPts, borderPoints = TRUE)
    out_snn[[city_id]] <- plot_flow_topN(
      city_id = city_id,
      prepared = prepared_plot,
      labels = as.integer(fit$cluster),
      algo_name = "SNN",
      out_dir = out_dir,
      top_n_clusters = top_n_clusters,
      show_background = show_background,
      show_plots = show_plots
    )
  }

  list(dbscan = out_dbscan, snn = out_snn)
}

# -------------------------
# 9) End-to-end verbose runner (prints tables/plots like notebook workflow)
# -------------------------
run_pipeline_verbose <- function(city_files = NULL,
                                 base_dir = ".",
                                 dbscan_params = list(
                                   berlin = list(minPts = 20, eps = 0.200),
                                   munich = list(minPts = 20, eps = 0.200),
                                   cologne = list(minPts = 10, eps = 0.250)
                                 ),
                                 snn_params = list(
                                   berlin = list(k = 20, snn_threshold = 8, minPts = 12),
                                   munich = list(k = 30, snn_threshold = 11, minPts = 18),
                                   cologne = list(k = 15, snn_threshold = 6, minPts = 10)
                                 ),
                                 popsize = 100,
                                 generations = 50,
                                 n_sub = 10000L,
                                 show_plots = TRUE,
                                 print_all_tables = TRUE) {
  setup_packages(c("data.table", "sf", "dbscan", "dplyr", "ggplot2", "mco"))

  # Start state
  init <- initialize_pipeline_state(city_files, base_dir = base_dir)
  prepared <- init$prepared

  cat("\n===== CLEANING SUMMARY =====\n")
  print(init$cleaning_summary)

  # Grid searches
  eps_grid <- seq(0.08, 0.40, by = 0.02)
  minPts_grid_dbscan <- c(10, 20, 30)
  k_grid_snn <- c(10, 20, 30)
  make_snn_threshold_grid <- function(k) {
    thr <- round(seq(0.35 * k, 0.65 * k, length.out = 6))
    unique(pmax(2L, pmin(k - 1L, as.integer(thr))))
  }
  minPts_grid_snn <- c(15, 18, 20, 25)

  grid_df_dbscan <- run_dbscan_grid_search(prepared, eps_grid, minPts_grid_dbscan)
  grid_df_snn <- run_snn_grid_search(prepared, k_grid_snn, make_snn_threshold_grid, minPts_grid_snn)

  cat("\n===== DBSCAN GRID RESULTS =====\n")
  print(grid_df_dbscan)
  cat("\n===== SNN GRID RESULTS =====\n")
  print(grid_df_snn)

  # Final labels + summaries
  prepared_db <- apply_final_dbscan(prepared, dbscan_params = dbscan_params)
  prepared_snn <- apply_final_snn(prepared, snn_params = snn_params)

  cat("\n===== DBSCAN CLUSTER SIZE CHECK =====\n")
  verify_cluster_sizes(prepared_db)
  cat("\n===== SNN CLUSTER SIZE CHECK =====\n")
  verify_cluster_sizes(prepared_snn)

  cat("\n===== DBSCAN CLUSTER SUMMARIES =====\n")
  cluster_summaries_db <- build_cluster_summaries(prepared_db)
  cat("\n===== SNN CLUSTER SUMMARIES =====\n")
  cluster_summaries_snn <- build_cluster_summaries(prepared_snn)

  # NSGA-II on shared subsample
  trip_dirs_full <- build_trip_dirs_full(prepared)
  opt <- build_opt_subsample(prepared, trip_dirs_full, n_sub = n_sub, seed = 1)

  db_nsga <- run_nsga_dbscan(prepared, opt$data_for_opt, opt$trip_dirs_for_opt,
                             popsize = popsize, generations = generations)
  pareto_dbscan <- build_pareto_dbscan(db_nsga$nsga_dbscan,
                                       eps_min = db_nsga$eps_min,
                                       eps_max = db_nsga$eps_max,
                                       eps_step = db_nsga$eps_step,
                                       eps_digits = db_nsga$eps_digits,
                                       minPts_allowed = db_nsga$minPts_allowed)
  knee_dbscan <- select_knee_points(pareto_dbscan)

  snn_nsga <- run_nsga_snn(prepared, opt$data_for_opt, opt$trip_dirs_for_opt,
                           popsize = popsize, generations = generations)
  pareto_snn <- build_pareto_snn(snn_nsga$nsga_snn,
                                 k_allowed = snn_nsga$k_allowed,
                                 thr_allowed = snn_nsga$thr_allowed,
                                 minPts_allowed = snn_nsga$minPts_allowed)
  knee_snn <- select_knee_points(pareto_snn)

  if (isTRUE(print_all_tables)) {
    cat("\n===== ALL DBSCAN PARETO TABLES =====\n")
    print_all_pareto(pareto_dbscan, algo = "DBSCAN")
    cat("\n===== ALL SNN PARETO TABLES =====\n")
    print_all_pareto(pareto_snn, algo = "SNN")
  }

  # Full-data validations
  valid_db <- validate_pareto_dbscan_full(prepared, pareto_dbscan, knee_dbscan, dbscan_params = dbscan_params)
  valid_snn <- validate_pareto_snn_full(prepared, pareto_snn, knee_snn, snn_params = snn_params)

  cat("\n===== DBSCAN FULL VALIDATION =====\n")
  print(valid_db$full_validated_dbscan_df)
  cat("\n===== DBSCAN EXPERT vs KNEE_FULL =====\n")
  print(valid_db$compare_dbscan)

  cat("\n===== SNN FULL VALIDATION =====\n")
  print(valid_snn$full_validated_snn_df)
  cat("\n===== SNN EXPERT vs KNEE_FULL =====\n")
  print(valid_snn$compare_snn)

  # Hypervolume DBSCAN vs SNN
  hv <- compute_hv_results(pareto_by_group = list(DBSCAN = pareto_dbscan, SNN = pareto_snn),
                           ref_min = c(1.05, 1.05),
                           group_col = "algorithm",
                           city_col = "city")
  cat("\n===== HYPERVOLUME RESULTS (DBSCAN vs SNN) =====\n")
  print(hv$hv_results)

  if (isTRUE(show_plots) && nrow(hv$hv_summary_df) > 0) {
    plot_hv_summaries(hv$hv_results, hv$hv_summary_df)
    for (city_id in names(pareto_dbscan)) {
      if (!is.null(pareto_dbscan[[city_id]]) && nrow(pareto_dbscan[[city_id]]) > 0) {
        print(plot_pareto(pareto_dbscan[[city_id]], knee = knee_dbscan[[city_id]], city_name = city_id,
                          title_prefix = "[DBSCAN] Pareto front (subsample)"))
      }
    }
    for (city_id in names(pareto_snn)) {
      if (!is.null(pareto_snn[[city_id]]) && nrow(pareto_snn[[city_id]]) > 0) {
        print(plot_pareto(pareto_snn[[city_id]], knee = knee_snn[[city_id]], city_name = city_id,
                          title_prefix = "[SNN] Pareto front (subsample)"))
      }
    }
  }

  # Optional null-test + flow plots
  null_out <- run_dbscan_null_test(prepared_real = prepared, n_sub = n_sub, popsize = popsize, generations = generations)
  cat("\n===== NULL TEST HV (REAL vs RAND) =====\n")
  print(null_out$hv_results_tag)

  selected_params_dbscan <- get_dbscan_selected_params(prepared, knee_full_dbscan = valid_db$knee_full_dbscan, dbscan_params = dbscan_params)
  selected_params_snn <- get_snn_selected_params(prepared, knee_full_snn = valid_snn$knee_full_snn, snn_params = snn_params)
  flow_plots <- run_and_plot_selected_flows(
    prepared_plot = prepared,
    selected_params_dbscan = selected_params_dbscan,
    selected_params_snn = selected_params_snn,
    out_dir = "plots_flows_top3",
    top_n_clusters = 3,
    show_background = FALSE,
    show_plots = show_plots
  )

  list(
    init = init,
    prepared = prepared,
    grid_df_dbscan = grid_df_dbscan,
    grid_df_snn = grid_df_snn,
    cluster_summaries_db = cluster_summaries_db,
    cluster_summaries_snn = cluster_summaries_snn,
    pareto_dbscan = pareto_dbscan,
    pareto_snn = pareto_snn,
    knee_dbscan = knee_dbscan,
    knee_snn = knee_snn,
    valid_db = valid_db,
    valid_snn = valid_snn,
    hv = hv,
    null_out = null_out,
    flow_plots = flow_plots
  )
}

# Example full-report run (explicit city_files):
# city_files <- list(
#   berlin  = "dt_bolt_berlin_06_05.rds",
#   cologne = "dt_voi_cologne_06_05.rds",
#   munich  = "dt_voi_munich_06_05.rds"
# )
# report <- run_pipeline_verbose(city_files, show_plots = TRUE, print_all_tables = TRUE)


# -------------------------
# -------------------------
# 10) RUN COMMAND (non-commented): shareable execution block
# -------------------------
# Your exact files (as requested):
CITY_FILES <- list(
  berlin = "dt_bolt_berlin_06_05.rds",
  cologne = "dt_voi_cologne_06_05.rds",
  munich = "dt_voi_munich_06_05.rds"
)

RUN_PIPELINE_NOW <- TRUE  # set FALSE if you only want function definitions

if (isTRUE(RUN_PIPELINE_NOW)) {
  init <- initialize_pipeline_state(city_files = CITY_FILES)
  report <- run_pipeline_verbose(
    city_files = CITY_FILES,
    show_plots = interactive(),
    print_all_tables = TRUE
  )
}

