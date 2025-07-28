# This script is to plot the map
# Author: Haitao Gao
# Date: 21st July 2025

## Load the package
library(sf)
library(ggplot2)
library(ggmap)
library(terra)
library(data.table)
library(dplyr)
library(evd)
library(parallel)

## read the file
load("dataset/application_florida.RData")



# ## must include other DOPGrid file except shp file like .shx and .dbf
grid_sf <- read_sf("dataset/DOPGrid.shp")
# 
all_IDs <- names(read.csv("dataset/PixelRain15min_1995.csv", header = TRUE, nrows = 1))[-1]

all_IDs_num <-  as.numeric(stringi::stri_extract_first(all_IDs, regex = "[0-9]+"))
# 
# fill_grid = grid_sf$PIXEL %in% all_IDs_num
# 
## create a Google Map static key
register_google(key = "")
map <- get_googlemap(center=c(-82.273528,28.209394),zoom=9,maptype = "terrain",style = "feature:all|element:all|saturation:-100|lightness:50")
# 
# ## Transform the grid_sf to latitude and longitude in decimal degrees
# grid_transform <- st_transform(grid_sf, crs = 4326)
# 
# 
# 
# ggmap(map) + theme_void() +
# ggtitle("Tampa Bay") + theme(plot.title = element_text(hjust = 0.5))  +
# geom_sf(data=grid_transform, aes(colour=as.factor(fill_grid)),alpha=0.5,inherit.aes = FALSE) + scale_color_brewer(palette = "Set1",name="Grid")
# 
# 
# coord_geo <- as.data.frame(st_coordinates(grid_transform))
# idx_pixel <- unlist(lapply(all_IDs_num,function(x){which(grid_transform$PIXEL==x)}))
# coord_geo <- matrix(unlist(lapply(idx_pixel, function(i){apply(coord_geo[coord_geo$L2==i,1:2][1:4,],2,mean)})),ncol=2,byrow=TRUE)


grid_transform <- st_transform(grid_sf, crs = 4326)
coord <- as.data.frame(st_coordinates(grid_transform))

idx_pixel <- unlist(lapply(all_IDs_num,function(x){which(grid_sf$PIXEL==x)}))
# 
# coord_grid <- matrix(unlist(lapply(idx_pixel, function(i){apply(coord[coord$L2==i,1:2][1:4,],2,mean)})),ncol=2,byrow=TRUE)
# 
# coord_df <- as.data.frame(coord_grid)
# names(coord_df) <- c("x", "y")
# 
# coord_df$y_rounded <- round(coord_df$y, 4)
# 
# length(unique(coord_df$y_rounded))
# 
# coord_df$row_id <- as.numeric(factor(coord_df$y_rounded))
# 
# 
# coord_trimmed <- coord_df %>%
#   group_by(row_id) %>%
#   arrange(x) %>%
#   slice_tail(n = min_count) %>%
#   ungroup()
# 
# coord_final <- as.matrix(coord_trimmed[, c("x", "y")])


centers_df <- do.call(rbind, lapply(idx_pixel, function(i) {
  center <- colMeans(coord[coord$L2 == i, 1:2][1:4, ])
  data.frame(x = center[1], y = center[2], pixel_idx = i)
}))

centers_df$y_rounded <- round(centers_df$y, 4)
centers_df$row_id <- as.numeric(factor(centers_df$y_rounded))

row_lengths <- centers_df %>%
  group_by(row_id) %>%
  summarise(count = n())

min_count <- min(row_lengths$count)

coord_trimmed <- centers_df %>%
  group_by(row_id) %>%
  arrange(x) %>%
  slice_tail(n = min_count) %>%
  ungroup()


selected_sf <- grid_transform[coord_trimmed$pixel_idx, ]


H <- 75
W <- 55
kernel <- 3  # kernel size and stride

coord_mat_x <- matrix(coord_trimmed$x, nrow = H, byrow = TRUE)
coord_mat_y <- matrix(coord_trimmed$y, nrow = H, byrow = TRUE)
coord_mat_id <- matrix(coord_trimmed$pixel_idx, nrow = H, byrow = TRUE)

rain_mat <- tmp[, -1]  # pixel columns only

H_out <- floor((H - kernel) / kernel + 1)
W_out <- floor((W - kernel) / kernel + 1)

down_coord_proj <- matrix(NA, nrow = H_out * W_out, ncol = 2)  # x1, y1
down_coord_geo  <- matrix(NA, nrow = H_out * W_out, ncol = 2)  # x2, y2
down_pixel_idx  <- rep(NA, H_out * W_out)

idx <- 1
for (i in 1:H_out) {
  for (j in 1:W_out) {
    row_range <- (kernel * (i - 1) + 1):(kernel * (i - 1) + kernel)
    col_range <- (kernel * (j - 1) + 1):(kernel * (j - 1) + kernel)
    
    block_ids <- coord_mat_id[row_range, col_range]
    
    # Get PIXEL name strings
    pixel_names <- paste0("P", grid_sf$PIXEL[block_ids])
    
    # Compute total rainfall for each pixel in this block
    block_sums <- colSums(rain_mat[, pixel_names, drop = FALSE], na.rm = TRUE)
    
    # Get the pixel with the maximum total rainfall
    max_idx <- which.max(block_sums)
    selected_id <- block_ids[max_idx]
    
    # Save projected coords (x1, y1)
    # Get geometry before transformation
    selected_geom_proj <- st_geometry(grid_sf[selected_id, ])
    
    # Compute projected centroid: x1, y1 (in meters)
    center_proj <- st_coordinates(st_centroid(selected_geom_proj))
    down_coord_proj[idx, ] <- center_proj
    
    # Transform to lat-lon for x2, y2 (in degrees)
    selected_geom_latlon <- st_transform(selected_geom_proj, crs = 4326)
    center_latlon <- st_coordinates(st_centroid(selected_geom_latlon))
    down_coord_geo[idx, ] <- center_latlon
    
    
    # Transform to lat-lon and compute centroid (x2, y2)
    selected_geom <- st_geometry(grid_sf[selected_id, ])
    selected_geom_latlon <- st_transform(selected_geom, crs = 4326)
    center_latlon <- st_coordinates(st_centroid(selected_geom_latlon))
    
    down_coord_geo[idx, ] <- center_latlon
    down_pixel_idx[idx] <- selected_id
    
    idx <- idx + 1
  }
}


down_df <- data.frame(
  x1 = down_coord_proj[, 1],
  y1 = down_coord_proj[, 2],
  x2 = down_coord_geo[, 1],
  y2 = down_coord_geo[, 2],
  pixel_idx = down_pixel_idx
)




rain_mat <- tmp[, -1]

# Build mapping from pixel index (grid_sf row index) → column name
pixel_idx_to_colname <- function(pixel_idx) {
  paste0("P", grid_sf$PIXEL[pixel_idx])
}

# Convert all selected pixel indices to their corresponding column names
col_ids <- sapply(down_df$pixel_idx, pixel_idx_to_colname)

# Subset rain_mat to retain only selected columns (pixels)
rain_down <- rain_mat[, col_ids, drop = FALSE]



# Step 1: Rank transform (column-wise) to uniform scale
rank_colwise <- function(j) {
  x <- rain_down[, j]
  ind <- x > 0
  x[ind] <- rank(x[ind]) / (sum(ind) + 1)
  return(x)
}

# Apply rank transform in parallel to each column
system.time({
  rain_down_uniform <- do.call(cbind, mclapply(1:ncol(rain_down), rank_colwise, mc.cores = 4))
})

# Step 2: Identify and keep only non-zero rows
nonzeros_row <- apply(rain_down_uniform, 1, function(x) any(x > 0))
rain_down_uniform_filtered <- rain_down_uniform[nonzeros_row, ]

# Optional: how many non-zero rows?
num.nonzeros_row <- apply(rain_down_uniform, 1, function(x) sum(x > 0))

# Step 3: Reformat each row into list: (non-zero indices, values)
func_extract_nonzeros <- function(i) {
  x <- rain_down_uniform_filtered[i, ]
  ind <- which(x > 0)
  list(ind, x[ind])
}

system.time({
  data_down_raw <- mclapply(1:nrow(rain_down_uniform_filtered), func_extract_nonzeros, mc.cores = 4)
})

# Step 4: Apply GPD transform
data_down_list <- mclapply(data_down_raw, function(x) {
  list(x[[1]], qgpd(x[[2]], scale = 1, shape = 1))
}, mc.cores = 4)

# Step 5: Summary statistics for filtering
data_down_mean <- sapply(data_down_list, function(x) if (length(x[[2]]) > 0) mean(x[[2]]) else NA)
data_down_max  <- sapply(data_down_list, function(x) if (length(x[[2]]) > 0) max(x[[2]])  else NA)



data_down_fit_sum <- data_down_list[data_down_mean > quantile(data_down_mean, 0.9995, na.rm = TRUE)]
data_down_fit_max <- data_down_list[data_down_max  > quantile(data_down_max,  0.9995, na.rm = TRUE)]

down_coord_grid <- as.matrix(down_df[, c("x1", "y1")]) 


ggmap(map) + 
  theme_void() +
  ggtitle("Tampa Bay: Full Grid (blue), Trimmed Grid (red), Downsampled Grid (yellow)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  
  # Full grid in light gray
  geom_sf(data = grid_transform, fill = NA, color = "blue", size = 0.1, inherit.aes = FALSE) +
  
  # Trimmed grid in red (polygons)
  geom_sf(data = grid_transform[coord_trimmed$pixel_idx, ], fill = "red", color = "black", size = 0.2, inherit.aes = FALSE) +
  
  # Downsampled centers in yellow (use transformed coords)
  geom_point(data = down_df, aes(x = x2, y = y2), color = "gold", size = 1.5, inherit.aes = FALSE)



# 
# load("dataset/application_florida_list.RData")
# 
save(data_down_raw, down_coord_grid, file = "dataset/Florida_new.RData")



                                                        