# helper functions for manually calculating Voronoi areas


make_round_mask = function(radius){
  # makes a raster round mask. It assumes things are pixelized, hence +1
  width = radius * 2 + 1
  height = radius * 2 + 1
  center = radius + 1
  mask = matrix(FALSE, nrow = width, ncol = height)
  x_locs = matrix(rep(1:width, height), byrow = TRUE, nrow = width)
  y_locs = matrix(rep(1:height, width), byrow = FALSE, nrow = width)
  mask[((x_locs - center) ^ 2 + (y_locs - center) ^ 2) <= radius ^ 2] = TRUE
  mask
}

make_rect_mask = function(width, height){
  mask = matrix(TRUE, nrow = width + 1, ncol = height + 1)
}



get_squared_distance_array = function(x, y, mask){
  width = dim(mask)[1]
  height = dim(mask)[2]
  y_locs = matrix(rep(1:width, height), byrow = TRUE, nrow = width)
  x_locs = matrix(rep(1:height, width), byrow = FALSE, nrow = width)
  dist_array = array(0, dim = c(length(x), width, height))
  for (i in 1:length(x)){
    dist_array[i,,] = (x_locs - x[i])^2 + (y_locs - y[i]) ^ 2
    dist_array[i,,] = dist_array[i,,] * mask
  }
  dist_array
}

apply_weight = function(dist_array, weights){
  # multiples each distance matrix in dist_array by the 1 / weights in weights
  # multiples by the inverse so higher weights -> smaller distances
  for (i in 1:length(weights)){
    dist_array[i,,] = dist_array[i,,] / weights[i]
  }
  dist_array
}

get_closest_map = function(dist_array, mask){
  # takes an array of length i,x,y where each sub matrix for each i is a 
  # distance matrix. Then, finds which distance is the minimum, and returns this
  # "map"
  new_map = apply(dist_array, MARGIN = c(2,3), FUN = function(x){
    winners = which(x == min(x))
    if (length(winners) > 1){
      winners = sample(winners, 1)
    }
    winners
  })
  new_map * mask
}

get_areas = function(area_map, n){
  areas = sapply(1:n, FUN = function(x) sum(area_map == x))
  areas
}



