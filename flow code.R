 install.packages(c("data.table"))  # run once if not installed
library(data.table)

berlin <- readRDS("dt_bolt_berlin_06_05.rds")

class(berlin)
dim(berlin)
names(berlin)

needed <- c("start_loc_lon", "start_loc_lat", "dest_loc_lon", "dest_loc_lat")
needed %in% names(berlin)

library(data.table)

setDT(berlin)

berlin <- berlin[
  !is.na(start_loc_lat) & !is.na(start_loc_lon) &
    !is.na(dest_loc_lat) & !is.na(dest_loc_lon)
]

berlin <- berlin[
  start_loc_lat >= -90 & start_loc_lat <= 90 &
    dest_loc_lat   >= -90 & dest_loc_lat   <= 90 &
    start_loc_lon >= -180 & start_loc_lon <= 180 &
    dest_loc_lon   >= -180 & dest_loc_lon   <= 180
]

class(berlin)
dim(berlin)
needed %in% names(berlin)

install.packages("sf")   # only if not installed
library(sf)

orig_sf <- st_as_sf(
  berlin,
  coords = c("start_loc_lon", "start_loc_lat"),
  crs = 4326,
  remove = FALSE
)

dest_sf <- st_as_sf(
  berlin,
  coords = c("dest_loc_lon", "dest_loc_lat"),
  crs = 4326,
  remove = FALSE
)

orig_utm <- st_transform(orig_sf, 25832)
dest_utm <- st_transform(dest_sf, 25832)

orig_xy <- st_coordinates(orig_utm)
dest_xy <- st_coordinates(dest_utm)

berlin[, `:=`(
  x_o = orig_xy[, 1],
  y_o = orig_xy[, 2],
  x_d = dest_xy[, 1],
  y_d = dest_xy[, 2]
)]

summary(berlin$x_o)
summary(berlin$y_o)
summary(berlin$x_d)
summary(berlin$y_d)

berlin[, `:=`(
  dx = x_d - x_o,
  dy = y_d - y_o
)]

berlin[, `:=`(
  angle = atan2(dy, dx),     # radians [-pi, pi]
  len = sqrt(dx^2 + dy^2)    # meters
)]

summary(berlin$len)
range(berlin$angle)
sum(berlin$len == 0)

saveRDS(berlin, "C:/Users/User/Downloads/berlin_flows_utm.rds")
file.exists("C:/Users/User/Downloads/berlin_flows_utm.rds")

names(berlin)
berlin <- readRDS("berlin_flows_utm.rds")

install.packages(c("dbscan", "ggplot2", "dplyr"))
library(dbscan)
library(ggplot2)
library(dplyr)

# berlin is your data.frame

# new code for structural correction
X <- berlin %>%
  dplyr::select(x_o, y_o, dx, dy) %>%
  as.matrix()

ok <- complete.cases(X)
X_scaled <- scale(X[ok, , drop = FALSE])
berlin_clean <- berlin[ok, ]

# standardize (mean 0, sd 1)
X_scaled <- scale(X)

X <- X[ok, , drop = FALSE]
berlin_clean <- berlin[ok, , drop = FALSE]
X_scaled <- scale(X)
# kNN-distance plots (safe version)
k_values <- c(10, 20, 30, 50)

old_par <- par(no.readonly = TRUE)
par(mfrow = c(2, 2))

for (k in k_values) {
  kNN <- kNNdist(X_scaled, k = k)
  plot(sort(kNN), type = "l",
       xlab = "Points sorted by distance",
       ylab = paste0(k, "-NN distance"),
       main = paste("kNN distance plot (k =", k, ")"))
  grid()
}

par(old_par)
eps_grid <- c(0.2, 0.3, 0.4, 0.5, 0.6)
k_grid   <- c(20, 30, 50)

# new parameter search with adjusted eps
results <- list()
idx <- 1

for (k in k_grid) {
  for (eps_val in eps_grid) {

    db <- dbscan(X_scaled, eps = eps_val, minPts = k)
    cl <- db$cluster

    n_noise <- sum(cl == 0)
    noise_pct <- 100 * n_noise / length(cl)
    n_clusters <- length(unique(cl[cl != 0]))

    max_cluster_size <- if (n_clusters > 0) max(table(cl[cl != 0])) else 0
    max_cluster_pct  <- if (n_clusters > 0) 100 * max_cluster_size / length(cl) else 0

    results[[idx]] <- data.frame(
      minPts = k,
      eps = eps_val,
      n_clusters = n_clusters,
      noise_pct = noise_pct,
      max_cluster_pct = max_cluster_pct
    )

    idx <- idx + 1
  }
}

grid_df <- dplyr::bind_rows(results)
print(grid_df)

# Recommended Berlin choice 
minPts = 30
eps    = 0.2

# Cologne choice
minPts = 30
eps    = 0.3


# Munich choice
minPts = 20
eps    = 0.2

eps_star <- 0.2
minPts_star <- 30
# Final DBSCAN run
db_final <- dbscan(X_scaled, eps = eps_star, minPts = minPts_star)

berlin_clean$cluster <- db_final$cluster  # 0 = noise

# Verify the final clustering
table(berlin_clean$cluster)[1:10]

# cluster summary
cluster_summary <- berlin_clean %>%
  filter(cluster != 0) %>%          # exclude noise
  group_by(cluster) %>%
  summarise(
    size = n(),
    mean_x_o = mean(x_o),
    mean_y_o = mean(y_o),
    mean_dx  = mean(dx),
    mean_dy  = mean(dy),
    mean_len = mean(len, na.rm = TRUE),
    mean_duration = mean(duration, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(size))

print(head(cluster_summary, 15))

# excluding cluster 1
moo_clusters <- cluster_summary %>%
  dplyr::filter(cluster >= 2)

nrow(moo_clusters)

# vector of cluster IDs
cluster_ids <- moo_clusters$cluster
K <- length(cluster_ids)
K

# Unit direction vectors
v <- as.matrix(moo_clusters[, c("mean_dx", "mean_dy")])
norms <- sqrt(rowSums(v^2))
# avoid divide-by-zero (shouldn't happen if len>0, but safe)
norms[norms == 0] <- 1
u <- v / norms

# Cosine similarity matrix (K x K)
S <- u %*% t(u)
diag(S) <- NA  # ignore self-similarity

# weighted direction concentration
dir_concentration_weighted <- function(sel, u, w) {
  if (length(sel) < 2) return(0)
  ww <- w[sel] / sum(w[sel])
  m <- colSums(u[sel, , drop = FALSE] * ww)
  sqrt(sum(m^2))
}
# evaluation function
eval_solution <- function(x, df, u) {
  sel <- which(x == 1)

  total_size <- sum(df$size[sel])
  f1 <- -total_size

  dir_conc <- dir_concentration_weighted(sel, u, df$size)

  # Penalize selecting fewer than 2 clusters
  if (length(sel) < 2) {
    f2 <- 1   # worst (since we minimize)
  } else {
    f2 <- -dir_conc
  }

  c(f1_neg_size = f1,
    f2_neg_dir_concentration = f2)
}

eval_solution(x, moo_clusters, u)

eval_solution(rep(0, K), moo_clusters, u)
eval_solution(rep(1, K), moo_clusters, u)

K <- nrow(moo_clusters)

# select none
eval_solution(rep(0, K), moo_clusters, u)

# select all
eval_solution(rep(1, K), moo_clusters, u)
# sanity check
set.seed(1)
for (i in 1:5) {
  x <- rbinom(K, 1, 0.2)
  print(eval_solution(x, moo_clusters, u))
}

# random sampling check
set.seed(1)
N <- 500

vals <- matrix(NA, nrow = N, ncol = 2)

for (i in 1:N) {
  x <- rbinom(K, 1, 0.2)
  vals[i, ] <- eval_solution(x, moo_clusters, u)
}

vals_df <- data.frame(
  covered_size = -vals[,1],
  dir_concentration = -vals[,2]
)

summary(vals_df)
# ploting it
plot(vals_df$covered_size, vals_df$dir_concentration,
     xlab = "Covered demand (sum size)",
     ylab = "Direction concentration",
     main = "Random solutions: coverage vs direction concentration")
# rebuild pareto_df
F <- res$value
pareto_df <- data.frame(
  covered_demand = -F[,1],
  dir_concentration = -F[,2]
)

summary(pareto_df)

# knee detection
x <- pareto_df$covered_demand
y <- pareto_df$dir_concentration

x_n <- (x - min(x)) / (max(x) - min(x))
y_n <- (y - min(y)) / (max(y) - min(y))

dist <- sqrt((1 - x_n)^2 + (1 - y_n)^2)
knee_idx <- which.min(dist)

knee_idx
pareto_df[knee_idx, ]

# selected clusters
x_knee <- as.integer(res$par[knee_idx, ] > 0.5)
selected_clusters <- moo_clusters$cluster[x_knee == 1]

length(selected_clusters)
selected_clusters


# same as above
x_knee <- as.integer(res$par[knee_idx, ] > 0.5)
clusters_knee <- moo_clusters$cluster[x_knee == 1]

clusters_knee
length(clusters_knee)

# gives extrems and knee solution
best_dir_idx <- which.max(pareto_df$dir_concentration)
best_cov_idx <- which.max(pareto_df$covered_demand)

pareto_df[best_dir_idx, ]
pareto_df[best_cov_idx, ]
pareto_df[knee_idx, ]

# pareto solutions for abjectives
decode_solution <- function(sol_idx, res, moo_clusters) {
  x_bin <- as.integer(res$par[sol_idx, ] > 0.5)
  moo_clusters$cluster[x_bin == 1]
}

clusters_best_dir <- decode_solution(best_dir_idx, res, moo_clusters)
clusters_best_cov <- decode_solution(best_cov_idx, res, moo_clusters)
clusters_knee     <- decode_solution(knee_idx, res, moo_clusters)

length(clusters_best_dir); clusters_best_dir
length(clusters_best_cov); clusters_best_cov
length(clusters_knee);     clusters_knee

# Create a comparison table for the three solutions
make_solution_table <- function(selected, cluster_summary) {
  cluster_summary %>%
    dplyr::filter(cluster %in% selected) %>%
    dplyr::select(cluster, size, mean_dx, mean_dy, mean_len, mean_duration) %>%
    dplyr::arrange(dplyr::desc(size))
}

tab_best_dir <- make_solution_table(clusters_best_dir, cluster_summary)
tab_knee     <- make_solution_table(clusters_knee, cluster_summary)
tab_best_cov <- make_solution_table(clusters_best_cov, cluster_summary)

head(tab_best_dir, 10)
head(tab_knee, 10)
head(tab_best_cov, 10)

# map the knee solution
library(ggplot2)

plot_df <- cluster_summary %>%
  dplyr::filter(cluster %in% clusters_knee) %>%
  mutate(
    x1 = mean_x_o,
    y1 = mean_y_o,
    x2 = mean_x_o + mean_dx,
    y2 = mean_y_o + mean_dy
  )

ggplot(plot_df) +
  geom_segment(
    aes(x = x1, y = y1, xend = x2, yend = y2, linewidth = size),
    arrow = arrow(length = unit(0.15, "cm"))
  ) +
  scale_linewidth_continuous(range = c(0.3, 1.5)) +
  labs(
    title = "Knee solution: selected flow corridors",
    x = "x", y = "y"
  )

# better map
# Normalize arrow length (optional)
scale_factor <- 0.5

plot_df <- plot_df %>%
  mutate(
    x2 = mean_x_o + scale_factor * mean_dx,
    y2 = mean_y_o + scale_factor * mean_dy
  )

#Fix aspect ratio (important for maps)
ggplot(plot_df) +
  geom_segment(
    aes(x = x1, y = y1, xend = x2, yend = y2, linewidth = size),
    arrow = arrow(length = unit(0.15, "cm"))
  ) +
  coord_equal() +
  scale_linewidth_continuous(range = c(0.3, 1.5)) +
  labs(
    title = "Knee solution: selected flow corridors",
    x = "Easting (m)", y = "Northing (m)"
  )

