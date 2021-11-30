### Bayesian Additive Spanning Trees ###
library(igraph)
library(fields)
library(Matrix)
library(FNN)
source('FEMFun.R')

### Main functions -----

#' Train a BAST model
#'
#' @param Y Standardized vector of responses
#' @param graph0 Spatial graph. Should be an \code{igraph} object containing \code{length(Y)} vertices
#' @param init_val Named list of initial values. Should include the following items.
#'                 \code{'tree'}: A list of M initial spanning trees, where M is the number
#'                                of weak learners. Each should be an \code{igraph} object.
#'                 \code{'cluster'}: Initial cluster membership. Should be an integer matrix of size
#'                                   \code{length(Y) * M}. Each column should contain consecutive 
#'                                   integers from 1 to the number of clusters.
#'                 \code{'mu'}: A list of initial values of \eqn{\mu}. Each item should be a vector of
#'                              length equal to the number of clusters.
#'                 \code{'sigmasq_y'}: Initial value of noise variance.
#' @param hyperpar Named list of hyperparameters. Should include the following items.
#'                 \code{'M'}: Number of weak learners.
#'                 \code{'sigmasq_mu'}: Variance of the Gaussian prior for \eqn{\mu}.
#'                 \code{'lambda_s'}: Scale parameter of the prior for noise variance \eqn{\sigma^2_y}.
#'                 \code{'nu'}: Degree of freedom of the prior for noise variance \eqn{\sigma^2_y}.
#'                 \code{'lambda_k'}: Mean parameter of the prior for the number of clusters. 
#'                 \code{'k_max'}: Maximum number of clusters in each weak learner.
#' @param MCMC Number of MCMC iterations
#' @param BURNIN Number of burnin iterations
#' @param THIN Length thinning intervals. Will retain samples every \code{THIN} iterations
#' @param seed Random seed
#'
#' @return A list of MCMC samples containing the following items. 
#'         \code{'cluster_out'}: A list of \code{length(Y) * M} cluster membership matrices.
#'         \code{'mu_out'}: A list of posterior sample of \eqn{\mu}. \code{mu_out[[i]][[j]][k]} 
#'                          is the constant value for the k-th cluster in the j-th weak learner
#'                          in the i-th posterior sample.
#'         \code{'sigmasq_y_out'}: A vector of posterior draws of noise variance.
#'         \code{'log_post_out'}: A vector of log posterior densities.
#'          
fitBAST <- function(Y, graph0, init_val, hyperpar, MCMC, BURNIN, THIN, seed = 1234) {
  set.seed(seed)
  
  n = vcount(graph0)
  # hyper-parameter
  M = hyperpar['M']  # number of trees
  sigmasq_mu = hyperpar['sigmasq_mu']
  lambda_s = hyperpar['lambda_s']
  nu = hyperpar['nu']
  lambda_k = hyperpar['lambda_k']
  k_max = hyperpar['k_max']
  hyper = c(sigmasq_mu, lambda_s, nu, lambda_k)
  
  if('name' %in% names(vertex_attr(graph0))) {
    graph0 = delete_vertex_attr(graph0, 'name')
  }
  inc_mat = get.edgelist(graph0, names = F)
  adj_list = lapply(as_adj_list(graph0), FUN = function(x) {x$vid})
  adj_edge_list = lapply(as_adj_edge_list(graph0), FUN = function(x) {x$eid})
  
  
  # initial values
  mstgraph_lst = init_val[['trees']]
  mu = init_val[['mu']]
  cluster = init_val[['cluster']]  # n*M matrix
  sigmasq_y = init_val[['sigmasq_y']]
  k = as.numeric(apply(cluster, 2, max))  # number of clusters
  
  csize = list() # cluster size
  subgraphs = list()
  eid_btw_mst = list()
  g = matrix(0, nrow = n, ncol = M)  # n*M matrix of fitted mu's
  for(m in 1:M) {
    cluster_m = cluster[, m]
    g[, m] = mu[[m]][cluster_m]
    
    csize[[m]] = Rfast::Table(cluster_m)
    
    mstgraph_m = mstgraph_lst[[m]]
    clust_vid_m = split(1:n, cluster_m)
    subgraphs[[m]] = lapply(clust_vid_m, function(vids, mstgraph) {
      induced_subgraph(mstgraph, vids)
    }, mstgraph_m)
    
    inc_mat_mst = get.edgelist(mstgraph_m, names = F)
    # idx for bewteen edges
    c1_m = cluster_m[inc_mat_mst[, 1]]; c2_m = cluster_m[inc_mat_mst[, 2]]
    idx_btw = which(c1_m != c2_m)
    eid_btw_mst[[m]] = (E(mstgraph_m)$eid)[idx_btw]
  }
  
  # whether an edge in graph0 is within a cluster or bewteen two clusters
  # n*M matrix
  edge_status = apply(cluster, 2, FUN = getEdgeStatus, inc_mat)
  
  ################# MCMC ####################
  
  ## MCMC results
  mu_out = list()
  sigmasq_y_out = numeric((MCMC-BURNIN)/THIN)
  cluster_out = array(0, dim = c((MCMC-BURNIN)/THIN, n, M))
  tree_out = list()
  log_post_out = numeric((MCMC-BURNIN)/THIN)
  
  ## MCMC iteration
  for(iter in 1:MCMC) {
    for(m in 1:M) {
      k_m = k[m]
      mstgraph_m = mstgraph_lst[[m]]
      edge_status_m = edge_status[, m]
      cluster_m = cluster[, m]
      csize_m = csize[[m]]
      subgraphs_m = subgraphs[[m]]
      eid_btw_mst_m = eid_btw_mst[[m]]
      
      e_m = Y - Rfast::rowsums(g[, -m])
      
      if(k_m == 1) {rb = 0.9; rd = 0; rc = 0; rhy = 0.1
      } else if(k_m == min(k_max, n)) {rb = 0; rd = 0.6; rc = 0.3; rhy = 0.1
      } else {rb = 0.3; rd = 0.3; rc = 0.3; rhy = 0.1}
      move = sample(4, 1, prob = c(rb, rd, rc, rhy))
      
      if(move == 1) { ## Birth move
        # split an existing cluster
        split_res = splitCluster(mstgraph_m, k_m, subgraphs_m, csize_m)
        vid_new = split_res$vid_new; vid_old = split_res$vid_old
        
        # compute log-prior ratio
        log_A = log(lambda_k) - log(k_m + 1)
        # compute log-proposal ratio
        if(k_m == min(k_max, n)-1) {
          rd_new = 0.6
        } else {rd_new = 0.3}
        log_P = log(rd_new) - log(rb)
        # compute log-likelihood ratio
        log_L = evalLogLikeRatio('split', e_m, vid_old, vid_new, sigmasq_y, sigmasq_mu)
        
        #acceptance probability
        acc_prob = min(0, log_A + log_P + log_L)
        acc_prob = exp(acc_prob)
        if(runif(1) < acc_prob){
          # accept
          update_res = updateSplit(split_res, subgraphs_m, k_m, csize_m, eid_btw_mst_m, 
                                   cluster_m, edge_status_m, adj_list, adj_edge_list)
          subgraphs[[m]] = update_res$subgraphs
          csize[[m]] = update_res$csize
          eid_btw_mst[[m]] = update_res$eid_btw_mst
          cluster[, m] = update_res$cluster
          k[m] = k[m] + 1
          edge_status[, m] = update_res$estatus
        }
      }
      
      if(move == 2) { ## Death move
        # merge two existing clusters (c1, c2) -> c2
        merge_res = mergeCluster(mstgraph_m, eid_btw_mst_m, subgraphs_m, csize_m, 
                                 cluster_m, inc_mat)
        vid_old = merge_res$vid_old; vid_new = merge_res$vid_new
        
        # compute log-prior ratio
        log_A = -log(lambda_k) + log(k_m)
        # # compute log-proposal ratio
        if(k_m == 2) {rb_new = 0.9
        }else {rb_new = 0.3}
        log_P = -(log(rd) - log(rb_new))
        # compute log-likelihood ratio
        log_L = evalLogLikeRatio('merge', e_m, vid_old, vid_new, sigmasq_y, sigmasq_mu)
        
        # acceptance probability
        acc_prob = min(0, log_A + log_P + log_L)
        acc_prob = exp(acc_prob)
        if(runif(1) < acc_prob){
          # accept
          update_res = updateMerge(merge_res, subgraphs_m, csize_m, eid_btw_mst_m, 
                                   cluster_m, edge_status_m, adj_list, adj_edge_list, mstgraph_m)
          subgraphs[[m]] = update_res$subgraphs
          csize[[m]] = update_res$csize
          eid_btw_mst[[m]] = update_res$eid_btw_mst
          cluster[, m] = update_res$cluster
          k[m] = k[m] - 1
          edge_status[, m] = update_res$estatus
        }
      }
      
      if(move == 3) { ## change move
        # first perform death move: (c1, c2) -> c2
        merge_res = mergeCluster(mstgraph_m, eid_btw_mst_m, subgraphs_m, csize_m, 
                                 cluster_m, inc_mat, change = T)
        # then perform birth move
        split_res = splitCluster(mstgraph_m, k_m-1, merge_res$subgraphs, merge_res$csize)
        
        # compute log-likelihood ratio
        log_L = evalLogLikeRatio('merge', e_m, merge_res$vid_old, merge_res$vid_new, sigmasq_y, sigmasq_mu) + 
          evalLogLikeRatio('split', e_m, split_res$vid_old, split_res$vid_new, sigmasq_y, sigmasq_mu)
        
        # acceptance probability
        acc_prob = min(0, log_L)
        acc_prob = exp(acc_prob)
        if(runif(1) < acc_prob){
          # accept
          update_res = updateMerge(merge_res, subgraphs_m, csize_m, eid_btw_mst_m, 
                                   cluster_m, edge_status_m, adj_list, adj_edge_list, mstgraph_m)
          update_res = updateSplit(split_res, update_res$subgraphs, k_m-1, update_res$csize, 
                                   update_res$eid_btw_mst, update_res$cluster, 
                                   update_res$estatus, adj_list, adj_edge_list)
          
          subgraphs[[m]] = update_res$subgraphs
          csize[[m]] = update_res$csize
          eid_btw_mst[[m]] = update_res$eid_btw_mst
          cluster[, m] = update_res$cluster
          edge_status[, m] = update_res$estatus
        }
      }
      
      if(move == 4) {
        # update MST
        mstgraph_m = proposeMST(graph0, edge_status[, m], subgraphs_m)
        mstgraph_lst[[m]] = mstgraph_m$mstgraph
        
        # update eid_btw_mst
        eid_btw_mst[[m]] = mstgraph_m$eid_btw_mst
        # update subgraphs
        subgraphs[[m]] = mstgraph_m$subgraphs
      }
      
      # update mu_m
      k_m = k[m]; cluster_m = cluster[, m]; csize_m = csize[[m]]
      Qinv_diag = 1 / (csize_m/sigmasq_y + 1/sigmasq_mu)
      b = Qinv_diag * Rfast::group.sum(e_m, cluster_m) / sigmasq_y
      mu[[m]] = rnorm(k_m, b, sqrt(Qinv_diag))
      
      g[, m] = mu[[m]][cluster_m]
    }
    
    # update sigmasq_y
    Y_hat = g[, M] + Y - e_m
    rate = 0.5*(nu*lambda_s + sum((Y - Y_hat)^2))
    sigmasq_y = 1/rgamma(1, shape = (n+nu)/2, rate = rate)
    
    
    ## save result
    if(iter > BURNIN & (iter - BURNIN) %% THIN == 0) {
      mu_out[[(iter-BURNIN)/THIN]] = mu
      sigmasq_y_out[(iter-BURNIN)/THIN] = sigmasq_y
      #tree_out[[(iter-BURNIN)/THIN]] = mstgraph_lst
      cluster_out[(iter-BURNIN)/THIN, , ] = cluster
      
      log_post_out[(iter-BURNIN)/THIN] = evalLogPost(mu, g, sigmasq_y, k, Y, hyper)
    }
    
    if(iter %% 100 == 0)
      cat('Iteration', iter, 'done\n')
  }
  
  mode(cluster_out) = 'integer'  # to save memory
  return(list('mu_out' = mu_out,
              'sigmasq_y_out' = sigmasq_y_out,
              'cluster_out' = cluster_out, 'log_post_out' = log_post_out))
}


#' Prediction from BAST
#'
#' @param mcmc_res MCMC results returned by \code{fitBAST}.
#' @param coords A n*2 matrix of coordinates, where n is the number of training locations.
#' @param coords_new A n_new*2 matrix of coordinates, where n_new is the number of testing locations.
#' @param method Prediction method. Should be one of the followings.
#'               \code{'soft-knn'}: Soft prediction using K nearest neighbors
#'               \code{'soft-mesh'}: Soft prediction using constrained Delaunay triangulation mesh,
#'                                   for domains in R^2 only.
#' @param mesh Constrained Delaunay triangulation mesh. Should be a \code{mesh.2D} object in package \code{fdaPDE}.
#'             Required when \code{method == 'soft-mesh'}.
#' @param return_type Type of returned values. Should be one of the followings.
#'               \code{'mean'}: Return posterior mean prediction only, stored in a vector of length n_new
#'               \code{'all'}: Return all posterior prediction samples, stored in a n_new * n_post matrix,
#'                             where n_post is the number of posterior samples.
#' @param weighting Method to obtain soft prediction weights. Should be one of the followings.
#'               \code{'uniform'}: Uniform weights for all neighbors.
#'               \code{'distance'}: Inverse distance weighting.
#' @param cdist_mat Cross geodesic distance matrix of size n_new * n. Required when \code{method == 'soft-knn'}.
#' @param k_nn Number of neighbors for K nearest neighbors. Required when \code{method == 'soft-knn'}.
#' @param seed Random seed
#'
#' @return A vector of length n_new when \code{return_type == 'mean'}. A matrix of n_new * n_post
#'         when \code{return_type == 'all'}.
#'         
predictBAST <- function(mcmc_res, coords, coords_new, method = 'soft-knn', 
                        mesh = NULL, return_type = 'mean', weighting = 'uniform', 
                        cdist_mat = NULL, k_nn = 5, seed = 12345) {
  set.seed(seed)
  
  # get posterior mean in-sample prediction
  npost = length(mcmc_res$log_post_out)
  cluster_out = mcmc_res$cluster_out
  mu_out = mcmc_res$mu_out
  
  M = dim(cluster_out)[3]
  n = dim(cluster_out)[2]
  
  mu_all = array(0, dim = c(npost, n, M))
  for(i in 1:npost) {
    cluster_i = cluster_out[i, , , drop = F]
    mu_i = mu_out[[i]]
    mu_all_i = matrix(0, ncol = M, nrow = n) # n*M matrix
    for(m in 1:M) {
      mu_all_i[, m] = mu_i[[m]][cluster_i[1, , m]]
    }
    mu_all[i, , ] = mu_all_i
  }
  Y_hat_all = t(apply(mu_all, c(1, 2), sum))
  
  if(method == 'soft-mesh' | method == 'soft-knn') {
    ## soft prediction
    ## method == 'soft-mesh': soft prediction using FEM mesh (for 2-d domains only)
    ## method == 'soft-knn': soft prediction based on KNN graph
    
    Y_new_hat = softPredict2(mu_all, coords, coords_new, method, mesh,
                             return_type, weighting, cdist_mat, k_nn)
  } else {
    stop('Unsupported method')
  }
  
  return(Y_new_hat)
}


### Functions for fitting BAST -----

# function to get whether an edge is within a cluster or bewteen two clusters
getEdgeStatus <- function(membership, inc_mat) {
  membership_head = membership[inc_mat[, 1]]
  membership_tail = membership[inc_mat[, 2]]
  edge_status = rep('w', nrow(inc_mat))
  edge_status[membership_head != membership_tail] = 'b'
  return(edge_status)
}

# function to split an existing cluster given MST
# subgraphs: cluster_id -> subgraph
# cluster: vid -> cluster_id
splitCluster <- function(mstgraph, k, subgraphs, csize) { 
  clust_split = sample.int(k, 1, prob = csize - 1)
  mst_subgraph = subgraphs[[clust_split]]
  
  edge_cutted = E(mst_subgraph)[sample.int(csize[clust_split]-1, 1)]
  eid_cutted = edge_cutted$eid
  mst_subgraph = delete.edges(mst_subgraph, edge_cutted)
  connect_comp = components(mst_subgraph)
  idx_new = (connect_comp$membership == 2)
  vid_new = V(mst_subgraph)$vid[idx_new]
  vid_old = V(mst_subgraph)$vid[!idx_new]
  
  return(list(vid_old = vid_old, vid_new = vid_new, eid_cutted = eid_cutted,
              clust_old = clust_split, idx_new = idx_new))
}

# function to update if a split move is accepted
updateSplit <- function(split_res, subgraphs, k, csize, eid_btw_mst, cluster, edge_status,
                        adj_list, adj_edge_list) {
  clust_split = split_res$clust_old
  vid_old = split_res$vid_old; vid_new = split_res$vid_new
  
  subgraph_split = subgraphs[[clust_split]]
  idx_new = split_res$idx_new
  subgraphs[[clust_split]] = induced_subgraph(subgraph_split, !idx_new)  # subgraph of old cluster
  subgraphs[[k+1]] = induced_subgraph(subgraph_split, idx_new) # subgraph of new cluster
  
  csize[clust_split] = length(vid_old)
  csize[k+1] = length(vid_new)
  
  cluster[vid_new] = k + 1
  eid_btw_mst = c(eid_btw_mst, split_res$eid_cutted)
  
  # update edge status
  adj_vid_old = unlist(adj_list[vid_old])
  adj_eid_old = unlist(adj_edge_list[vid_old])
  clust_adj_old = cluster[adj_vid_old]
  idx_btw = which(clust_adj_old != clust_split)
  eid_btw = adj_eid_old[idx_btw]
  edge_status[eid_btw] = 'b'
  
  return(list(subgraphs = subgraphs, csize = csize, cluster = cluster,
              eid_btw_mst = eid_btw_mst, estatus = edge_status))
}

# function to merge two existing clusters
mergeCluster <- function(mstgraph, eid_btw_mst, subgraphs, csize, cluster, edge_list, 
                         change = F) {
  # edge for merging
  edge_merge = sample.int(length(eid_btw_mst), 1)
  # update cluster information
  # clusters of endpoints of edge_merge
  eid_merge = eid_btw_mst[edge_merge]
  clusters_merge = cluster[edge_list[eid_merge, ]]
  clusters_merge = sort(clusters_merge)
  c1 = clusters_merge[1]; c2 = clusters_merge[2] # note c1 < c2
  # merge c2 to c1
  
  # vid of vertices in c2
  vid_old = V(subgraphs[[c2]])$vid
  # vid in merged cluster
  vid_new = c(V(subgraphs[[c1]])$vid, vid_old)
  
  csize_new = NULL; subgraphs_new = NULL
  if(change) {
    subgraphs_new = subgraphs
    subgraphs_new[[c1]] = induced_subgraph(mstgraph, vid_new)
    subgraphs_new[[c2]] = NULL
    
    csize_new = csize
    csize_new[c1] = length(vid_new)
    csize_new = csize_new[-c2]
  }
  
  # now drop c2
  return(list(vid_old = vid_old, vid_new = vid_new, clust_old = c2, clust_new = c1,
              edge_merge = edge_merge, subgraphs = subgraphs_new, csize = csize_new))
}

# function to update if a merge move is accepted
updateMerge <- function(res_merge, subgraphs, csize, eid_btw_mst, cluster, edge_status,
                        adj_list, adj_edge_list, mstgraph) {
  clust_old = res_merge$clust_old; clust_new = res_merge$clust_new
  vid_old = V(subgraphs[[clust_old]])$vid
  vid_new = c(V(subgraphs[[clust_new]])$vid, vid_old)
  subgraphs[[clust_new]] = induced_subgraph(mstgraph, vid_new)
  subgraphs[[clust_old]] = NULL
  
  csize[clust_new] = length(vid_new)
  csize = csize[-clust_old]
  
  cluster[vid_old] = clust_new
  idx = which(cluster > clust_old)
  cluster[idx] = cluster[idx] - 1
  
  eid_btw_mst = eid_btw_mst[-res_merge$edge_merge]
  
  # update edge status
  adj_vid_old = unlist(adj_list[vid_old])
  adj_eid_old = unlist(adj_edge_list[vid_old])
  clust_adj_old = cluster[adj_vid_old]
  idx_within = which(clust_adj_old == clust_new)
  eid_within = adj_eid_old[idx_within]
  edge_status[eid_within] = 'w'
  
  return(list(subgraphs = subgraphs, csize = csize, cluster = cluster,
              eid_btw_mst = eid_btw_mst, estatus = edge_status))
}

# function to propose a new MST
proposeMST <- function(graph0, edge_status, subgraphs) {
  nedge = length(edge_status)
  nb = sum(edge_status == 'b')
  nw = nedge - nb
  weight = numeric(nedge)
  weight[edge_status == 'w'] = runif(nw)
  weight[edge_status == 'b'] = runif(nb, 10, 20)
  E(graph0)$weight = weight
  mstgraph = mst(graph0, weights = E(graph0)$weight)
  
  # update subgraphs
  subgraphs_new = lapply(subgraphs, function(g, mstgraph) {induced_subgraph(mstgraph, V(g)$vid)},
                         mstgraph)
  # update eid_btw_mst
  eid_btw_mst = E(mstgraph)$eid[E(mstgraph)$weight >= 10]
  
  return(list(mstgraph = mstgraph, subgraphs = subgraphs_new, eid_btw_mst = eid_btw_mst))
}

# function to get log likelihood ratio
evalLogLikeRatio <- function(move, e_m, vid_old, vid_new, sigmasq_y, sigmasq_mu) {
  sigma_ratio = sigmasq_y / sigmasq_mu
  csize_old = length(vid_old); csize_new = length(vid_new)
  sum_e_old = sum(e_m[vid_old]); sum_e_new = sum(e_m[vid_new])
  
  if(move == 'split') {
    logdetdiff = -0.5 * (log(csize_old+sigma_ratio) + log(csize_new+sigma_ratio) - 
                           log(csize_old+csize_new+sigma_ratio) - log(sigma_ratio))
    quaddiff = 0.5 * (sum_e_old^2/(csize_old+sigma_ratio) + sum_e_new^2/(csize_new+sigma_ratio) - 
                        (sum_e_old+sum_e_new)^2/(csize_old+csize_new+sigma_ratio)) / sigmasq_y
  }
  
  if(move == 'merge') {
    logdetdiff = -0.5 * (log(csize_new+sigma_ratio) - log(csize_old+sigma_ratio) - 
                           log(csize_new-csize_old+sigma_ratio) + log(sigma_ratio))
    quaddiff = 0.5 * (sum_e_new^2/(csize_new+sigma_ratio) - sum_e_old^2/(csize_old+sigma_ratio) - 
                        (sum_e_new-sum_e_old)^2/(csize_new-csize_old+sigma_ratio)) / sigmasq_y
  }
  return(logdetdiff + quaddiff)
}

# function to get log posterior density (up to a constant)
evalLogPost <- function(mu, mu_all, sigmasq_y, k, Y, hyper) {
  n = length(Y); M = length(k)
  sigmasq_mu = hyper[1]; lambda_s = hyper[2]; nu = hyper[3]; lambda_k = hyper[4]
  log_prior =  -(nu/2+1)*log(sigmasq_y) - nu*lambda_s/(2*sigmasq_y)
  log_prior = log_prior + sum(-lchoose(n-1, k-1) + k*log(lambda_k) - lfactorial(k))
  log_prior = log_prior - sum(unlist(mu)^2) / (2*sigmasq_mu)
  
  Y_hat = rowSums(mu_all)
  log_like = -n/2*log(sigmasq_y) - sum((Y - Y_hat)^2) / (2*sigmasq_y)
  return(log_like + log_prior)
}

# function to standardize Y
standardize <- function(x) {
  xmean = mean(x)
  x = x - xmean
  xscale = 2 * max(abs(x))
  x = x / xscale
  param = c('mean' = xmean, 'scale' = xscale)
  return(list(x = x, std_par = param))
}

# function to unstandardize Y
unstandardize <- function(x, std_par, nomean = F, s2 = F) {
  if(s2) {
    x = x * std_par['scale'] ^ 2
  } else {
    x = x * std_par['scale']
  }
  if(!nomean) x = x + std_par['mean']
  return(x)
}


### Prediction functions -----

# function to get adjacency list for boundary mesh nodes
getBndAdjList <- function(mesh) {
  n_int = mesh$n_int  # number of interior nodes
  n_bnd = nrow(mesh$nodes) - n_int  # number of boundary nodes
  edge_list = mesh$edges[mesh$bnd_edges, ]
  
  adj_list = rep(list(c()), n_bnd)  # empty list
  names(adj_list) = as.character((n_int + 1):(n_int + n_bnd))
  for(i in 1:nrow(edge_list)) {
    endpoints = edge_list[i, ]
    v_bnd = max(endpoints); v_int = min(endpoints)
    if(v_int > n_int) next
    v_bnd_char = as.character(v_bnd)
    if(is.null(adj_list[[v_bnd_char]])) {
      adj_list[[v_bnd_char]] = c(v_int)
    } else {
      adj_list[[v_bnd_char]] = c(adj_list[[v_bnd_char]], v_int)
    }
  }
  return(adj_list)
}

# function to get neighbor probability matrix for new locations given a mesh
# nn_prob[i, j]: probability that node j is used for prediction at new location i
# ncol(nn_prob) == number of nodes in mesh (including boundary nodes)
getNNProb <- function(mesh, coords_new, type = 'uniform', b = 1) {
  if(type == 'uniform' | type == 'distance') {
    tri_idx = apply(coords_new, 1, function(loc) R_insideIndex2(mesh, loc, all_tri = T))
    nn_prob = matrix(0, nrow = nrow(coords_new), ncol = nrow(mesh$nodes))
    
    n_nodes = nrow(mesh$nodes); n_int = mesh$n_int
    if(type == 'distance') {
      weights = 1 / as.matrix(rdist(coords_new, mesh$nodes)) ^ b
    } else {
      weights = matrix(1, nrow = nrow(coords_new), ncol = nrow(mesh$nodes))
    }
    
    if(class(tri_idx) == 'list') {
      for(i in 1:nrow(coords_new)) {
        nn_idx = mesh$triangles[tri_idx[[i]], ]
        nn_idx = unique(c(nn_idx))
        # use boundary nodes only if all neighbors are on boundary
        if(!all(nn_idx > n_int))
          nn_idx = nn_idx[nn_idx <= n_int]
        nn_prob[i, nn_idx] = weights[i, nn_idx]
      }
      
    } else { # tri_idx is a vector
      for(i in 1:nrow(coords_new)) {
        nn_idx = mesh$triangles[tri_idx[i], ]
        # use boundary nodes only if all neighbors are on boundary
        if(!all(nn_idx > n_int))
          nn_idx = nn_idx[nn_idx <= n_int]
        nn_prob[i, nn_idx] = weights[i, nn_idx]
      }
    }
    
    # normalize nn_prob
    nn_prob = nn_prob / rowSums(nn_prob)
    
  }  else {
    stop('Unsupported weighting types')
  }
  
  return(nn_prob)
}

# extend prediction to boundary nodes of a mesh
getBndPred <- function(Y_hat, mesh, adj_list_bnd = NULL) {
  if(is.null(adj_list_bnd))
    adj_list_bnd = getBndAdjList(mesh)
  n = length(Y_hat)  # number of interior nodes
  n_bnd = length(adj_list_bnd)  # number of boundary nodes
  Y_hat_bnd = rep(0, n_bnd)
  for(i in (n + 1):(n + n_bnd))
    Y_hat_bnd[i - n] = mean(Y_hat[ adj_list_bnd[[as.character(i)]] ])
  return(Y_hat_bnd)
}


# helper function for soft prediction
softPredict2 <- function(mu_all, coords, coords_new, method = 'soft-knn', 
                         mesh = NULL, return_type = 'mean', 
                         weighting = 'uniform', cdist_mat = NULL, k_nn = 5) {
  extend_bnd = F  # do we have to extend prediction to boundary nodes in mesh
  if(method == 'soft-mesh') {
    # soft prediction using mesh
    # require(fdaPDE)
    if(is.null(mesh)) 
      stop("Missing 'mesh'.")
    
    # get neighbor choices for new locations
    nn_prob = getNNProb(mesh, coords_new, weighting)
    nn_prob = round(nn_prob, 8)
    
    # check if boundary nodes are useful for prediction
    n_int = mesh$n_int
    if(any(nn_prob[, (n_int+1):ncol(nn_prob)] > 0))
      extend_bnd = T
    
    # pre-compute adjacency list for boundary mesh nodes
    if(extend_bnd)
      adj_list_bnd = getBndAdjList(mesh)
    
  } else if(method == 'soft-knn') {
    # soft prediction using constrained KNN
    require(fields)
    
    if(is.null(cdist_mat))
      cdist_mat = as.matrix(rdist(coords_new, coords))
    nn_list = KNNGraph(cdist_mat, k_nn, cross_dist = T, return_graph = F)

    # get probability for each neighbor
    n_new = nrow(coords_new)
    nn_prob = matrix(0, nrow = n_new, ncol = nrow(coords))
    if(weighting == 'uniform') {
      for(i in 1:n_new)
        nn_prob[i, nn_list[[i]] ] = 1
    } else if(weighting == 'distance') {
      # deal with zero distance
      cdist_mat[cdist_mat == 0] = 1e-9
      for(i in 1:n_new)
        nn_prob[i, nn_list[[i]] ] = 1 / cdist_mat[i, nn_list[[i]] ]
    } else {
      stop('Unsupported weighting types')
    }
    # normalize nn_prob
    nn_prob = nn_prob / rowSums(nn_prob)
  }
  
  n_nodes = ncol(nn_prob)
  npost = dim(mu_all)[1]
  M = dim(mu_all)[3]
  Y_new_all = matrix(0, nrow = nrow(coords_new), ncol = npost)
  
  
  for(i in 1:npost) {
    # sample neighbors
    nn = t(apply(nn_prob, 1, function(p) sample.int(n_nodes, M, replace = T, prob = p)))
    
    mu_i = mu_all[i, , ]
    
    if(method == 'soft-mesh' & extend_bnd) {
      # extend mu to boundary nodes in FEM mesh
      mu_i_bnd = apply(mu_i, 2, getBndPred, mesh = mesh, adj_list_bnd = adj_list_bnd)
      mu_i = rbind(mu_i, mu_i_bnd)
    }
    
    mu_new_i = matrix(0, nrow = nrow(coords_new), ncol = M)
    for(m in 1:M)
      mu_new_i[, m] = mu_i[nn[, m], m]
    Y_new_all[, i] = rowSums(mu_new_i)
  }
  if(return_type == 'mean') {
    Y_new_hat = rowMeans(Y_new_all)
  } else {
    Y_new_hat = Y_new_all
  }
  return(Y_new_hat)
}