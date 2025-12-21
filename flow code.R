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

X <- berlin %>%
  dplyr::select(x_o, y_o, x_d, y_d) %>%
  as.matrix()

# standardize (mean 0, sd 1)
X_scaled <- scale(X)

ok <- stats::complete.cases(X)
if (!all(ok)) {
  message("Removed ", sum(!ok), " rows with NA in (x_o,y_o,x_d,y_d).")
}
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
eps_grid <- c(0.5, 0.6, 0.7)
k_grid   <- c(15, 20, 30)
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
