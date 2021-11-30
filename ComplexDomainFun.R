### Functions related to complex domains ###

library(mgcv)
library(fdaPDE)
library(igraph)

### Functions for 2-d constrained domains -----

# function to create a 2d mesh
gen2dMesh <- function(coords, bnd, ...) {
  # note the first and last boundary nodes are the same
  n = nrow(coords); n_bnd = length(bnd$x) - 1
  coords_all = rbind(coords, cbind(bnd$x, bnd$y)[1:n_bnd, ])
  
  # get boundary segments
  segments = cbind( (n+1):(n+n_bnd), c((n+2):(n+n_bnd), n+1) )
  
  mesh = create.mesh.2D(coords_all, segments = segments, ...)
  mesh$n_int = n  # number of interior nodes
  # edges that connect boundary nodes
  mesh$bnd_edges = apply(mesh$edges, 1, FUN = function(x) any(x > n))
  return(mesh)
}

# function to obtain a constrained Delaunay triangulation graph from a mesh
constrainedDentri <- function(n, mesh, threshold = 1000) {
  coords = mesh$nodes[1:n, ]
  
  # drop edges that connect boundary nodes
  rid_drop = mesh$bnd_edges
  edge_list = mesh$edges[!rid_drop, ]
  
  # compute edge length
  distance = sqrt( rowSums((coords[edge_list[, 1], ] - coords[edge_list[, 2], ]) ^ 2) )
  
  rid_drop = distance > threshold
  edge_list = edge_list[!rid_drop, ]
  distance = distance[!rid_drop]
  
  graph0 = graph_from_edgelist(edge_list, directed = F)
  E(graph0)$weight = distance
  return(graph0)
}

# function to generate equally spaced points along a given line
refineLine <- function(start, end, n) {
  grids_x = seq(start[1], end[1], length.out = n)
  grids_y = seq(start[2], end[2], length.out = n)
  grids = cbind(grids_x, grids_y)
  return(grids)
}

### Constrained KNN -----

# function to get a KNN graph given a distance matrix
KNNGraph <- function(dist_mat, k_nn = 5, cross_dist = F, return_graph = T) {
  
  n1 = nrow(dist_mat); n2 = ncol(dist_mat)
  if(cross_dist) {
    # dist_mat is cross distance matrix
    adj_list = apply(dist_mat, 1, function(x) order(x)[1:k_nn])
  } else {
    adj_list = apply(dist_mat, 1, function(x) order(x)[2:(k_nn+1)])
  }
  adj_list = t(adj_list)
  
  if(return_graph & n1 == n2) {
    require(igraph)
    require(Matrix)
    
    i_all = c(); j_all = c(); x_all = c()
    for(i in 1:n1) {
      for(cidx in 1:k_nn) {
        i_all = c(i_all, i)
        j_all = c(j_all, adj_list[i, cidx])
        x_all = c(x_all, dist_mat[i, adj_list[i, cidx] ])
      }
    }
    adj_mat = sparseMatrix(i = i_all, j = j_all, x = x_all, dims = c(n1, n2))
    # get knn graph
    knngraph = graph_from_adjacency_matrix(adj_mat, mode = 'directed', weighted = T)
    knngraph = as.undirected(knngraph, mode = 'collapse', edge.attr.comb = 'first')
    return(knngraph)
  } else {
    return(lapply(1:n1, function(i) adj_list[i, ]))
  }
}


### Plotting functions -----

# function to plot complex domain
geom_boundary <- function(bnd, ...) {
  n = length(bnd$x)
  segments = cbind(bnd$x[-n], bnd$y[-n], bnd$x[-1], bnd$y[-1])
  segments = data.frame(segments)
  names(segments) = c('x1', 'y1', 'x2', 'y2')
  return(geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2), data = segments, ...))
}

# plotting functions from SCC
plotGraph <- function(coords, graph, title = NULL){
  require(ggplot2)
  edgelist = get.edgelist(graph) 
  edgedata = data.frame(coords[edgelist[,1 ], ], coords[edgelist[, 2], ])
  colnames(edgedata) = c("X1", "Y1", "X2", "Y2")
  
  ggplot() + geom_segment(aes(x=X1, y=Y1, xend = X2, yend = Y2), data = edgedata, size = 0.5, colour = "grey") +
    labs(title = title, x = "lon", y = "lat")+
    theme(plot.title = element_text(hjust = 0.5))
}

# function to plot spatial field
plotField <- function(coords, Data, col_lim = NULL, legend_name = NULL, title = NULL, colors = rainbow(10)){
  if(missing(col_lim)) {col_lim = range(Data)}
  ggplot() + 
    geom_tile(data = data.frame(coords), aes(lon, lat, fill = Data)) +
    scale_fill_gradientn(colours = colors, limits = col_lim, name = legend_name, na.value = 'white') +
    ggtitle(title) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          # legend.title=element_blank(), 
          plot.title = element_text(hjust = 0.5, size = 10),
          legend.text = element_text(size = 9, hjust = 0, margin = margin(l = 3)))
}
