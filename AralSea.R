###################################################
##### Preprocessing Aral Sea Chlorophyll Data #####
###################################################

rm(list=ls())
library(ggplot2)

# set working directory if necessary
setwd('D:/Documents/_TAMU/Projects/Github/BAST')
source('ComplexDomainFun.R')

### load data ----

data('aral', package = 'gamair')
data('aral.bnd', package = 'gamair')


### preprocess boundary and observed locations ----

# make boundary
lon = aral.bnd[[1]]
lat = aral.bnd[[2]]
lon[8:36] = aral.bnd[[1]][8:36] + 0.08
lon[43:46] = aral.bnd[[1]][43:46] - 0.02
lon[47:52] = aral.bnd[[1]][47:52] - 0.08
lon[58:60] = aral.bnd[[1]][58:60] + 0.08
lon[61:63] = aral.bnd[[1]][61:63] + 0.02
lon[64:65] = aral.bnd[[1]][64:65] + 0.05
lon[72:93] = aral.bnd[[1]][72:93] - 0.08
lon[100:107] = aral.bnd[[1]][100:107] - 0.08
lat[1:7] = aral.bnd[[2]][1:7] + 0.08
lat[37:42] = aral.bnd[[2]][37:42] - 0.08
lat[53:57] = aral.bnd[[2]][53:57] - 0.08
lat[94:99] = aral.bnd[[2]][94:99] + 0.08
lat[67:71] = aral.bnd[[2]][67:71] - 0.08

rmp = -c(99)
bnd = cbind(lon[rmp], lat[rmp])
bnd = rbind(bnd, bnd[1, ])

# refine boundary
bnd = rbind(refineLine(bnd[1, ], bnd[2, ], 5), bnd[3:nrow(bnd), ])
bnd = rbind(bnd[1:4, ], refineLine(bnd[5, ], bnd[6, ], 5), bnd[7:nrow(bnd), ])
bnd = rbind(bnd[1:12, ], refineLine(bnd[13, ], bnd[14, ], 3), bnd[15:nrow(bnd), ])
bnd = rbind(bnd[1:14, ], refineLine(bnd[15, ], bnd[16, ], 3), bnd[17:nrow(bnd), ])
bnd = rbind(bnd[1:18, ], refineLine(bnd[19, ], bnd[20, ], 3), bnd[21:nrow(bnd), ])
bnd = rbind(bnd[1:31, ], refineLine(bnd[32, ], bnd[33, ], 3), bnd[34:nrow(bnd), ])
bnd = rbind(bnd[1:45, ], refineLine(bnd[46, ], bnd[47, ], 7), bnd[48:nrow(bnd), ])
bnd = rbind(bnd[1:51, ], refineLine(bnd[52, ], bnd[53, ], 3), bnd[54:nrow(bnd), ])
bnd = rbind(bnd[1:56, ], refineLine(bnd[57, ], bnd[58, ], 8), bnd[59:nrow(bnd), ])
bnd = rbind(bnd[1:63, ], refineLine(bnd[64, ], bnd[65, ], 3), bnd[66:nrow(bnd), ])
bnd = rbind(bnd[1:67, ], refineLine(bnd[68, ], bnd[69, ], 3), bnd[70:nrow(bnd), ])
bnd = rbind(bnd[1:70, ], refineLine(bnd[71, ], bnd[72, ], 3), bnd[73:nrow(bnd), ])
bnd = rbind(bnd[1:76, ], refineLine(bnd[77, ], bnd[78, ], 3), bnd[79:nrow(bnd), ])
bnd = rbind(bnd[1:81, ], refineLine(bnd[82, ], bnd[83, ], 4), bnd[84:nrow(bnd), ])
bnd = rbind(bnd[1:86, ], refineLine(bnd[87, ], bnd[88, ], 3), bnd[89:nrow(bnd), ])
bnd = rbind(bnd[1:88, ], refineLine(bnd[89, ], bnd[90, ], 3), bnd[91:nrow(bnd), ])
bnd = rbind(bnd[1:96, ], refineLine(bnd[97, ], bnd[98, ], 4), bnd[99:nrow(bnd), ])
bnd = rbind(bnd[1:102, ], refineLine(bnd[103, ], bnd[104, ], 4), bnd[105:nrow(bnd), ])
bnd = rbind(bnd[1:107, ], refineLine(bnd[108, ], bnd[109, ], 3), bnd[110:nrow(bnd), ])
bnd = rbind(bnd[1:126, ], refineLine(bnd[127, ], bnd[128, ], 3), bnd[129:nrow(bnd), ])
bnd = rbind(bnd[1:128, ], refineLine(bnd[129, ], bnd[130, ], 3), bnd[131:nrow(bnd), ])
bnd = rbind(bnd[1:130, ], refineLine(bnd[131, ], bnd[132, ], 3), bnd[133:nrow(bnd), ])
bnd = rbind(bnd[1:143, ], refineLine(bnd[144, ], bnd[145, ], 3))

# rescale boundary and observed locations
shrink = 10

bnd[, 1] = (bnd[, 1] - mean(aral[, 1]) ) * shrink
bnd[, 2] = (bnd[, 2] - mean(aral[, 2]) ) * shrink
bnd = list('x' = bnd[, 1], 'y' = bnd[, 2])

bnd_original = aral.bnd
bnd_original$lon = (bnd_original$lon - mean(aral[, 1]) ) * shrink
bnd_original$lat = (bnd_original$lat - mean(aral[, 2]) ) * shrink

coords = as.matrix(aral[, 1:2])
coords[, 1] = (coords[, 1] - mean(aral[, 1]) ) * shrink
coords[, 2] = (coords[, 2] - mean(aral[, 2]) ) * shrink
colnames(coords) = c('lon', 'lat')

use = which( !is.na(aral[, 3]) )
coords = coords[use, ]


### center responses ----

Y = aral[use, 3]
mean_Y = mean(Y)
Y = as.numeric(scale(Y, center = TRUE, scale = FALSE))

# plot observed data
ggplot() + 
  geom_boundary(bnd) +
  geom_point(aes(x = lon, y = lat, col = Y), data = as.data.frame(coords)) +
  scale_color_gradientn(colours = rainbow(5), name = 'Y') +
  labs(x = 'Scaled Lon.', y = 'Scaled Lat.') + 
  ggtitle('Observed Data')


### make grid points ----

lon_unique = sort(unique(coords[, 1]))
lat_unique = sort(unique(coords[, 2]))
lon_gap = lon_unique[2] - lon_unique[1]
lat_gap = lat_unique[2] - lat_unique[1]
lon_min = lon_unique[1]; lon_max = lon_unique[length(lon_unique)]
lat_min = lat_unique[1]; lat_max = lat_unique[length(lat_unique)]

coords_grid = expand.grid(seq(lon_min - lon_gap, lon_max + lon_gap, lon_gap / 2),
                          seq(lat_min - lat_gap, lat_max + lat_gap, lat_gap / 2))
coords_grid = as.matrix(coords_grid)
inside_bnd = apply(coords_grid, 1, function(coord) {
  lon = coord[1]; lat = coord[2]
  mgcv::inSide(bnd_original, lon, lat)
})
coords_grid = coords_grid[inside_bnd, ]
colnames(coords_grid) = c('lon', 'lat')

# plot grid points
ggplot() + 
  geom_boundary(bnd) +
  geom_point(aes(x = lon, y = lat), data = as.data.frame(coords_grid)) +
  labs(x = 'Scaled Lon.', y = 'Scaled Lat.') + 
  ggtitle('Grid Points')


### save data ----
save(bnd, coords, Y, coords_grid, mean_Y, file = 'aral_data.RData')
