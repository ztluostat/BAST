####################################################################
##### Demo Code of Bayesian Additive Regression Spanning Trees #####
####################################################################

rm(list=ls())
library(igraph)
library(fields)

# set working directory if necessary
setwd('D:/Documents/_TAMU/Projects/Github/BAST')
source('ComplexDomainFun.R')
source('BASTFun.R')

### load data ----

load('aral_data.RData')
n = length(Y)

# plot observed data
ggplot() + 
  geom_boundary(bnd) +
  geom_point(aes(x = lon, y = lat, col = Y), data = as.data.frame(coords)) +
  scale_color_gradientn(colours = rainbow(5), name = 'Y') +
  labs(x = 'Scaled Lon.', y = 'Scaled Lat.') + 
  ggtitle('Observed Data')


### initializing parameters ----

# get mesh and triangulation
mesh = gen2dMesh(coords, bnd)
graph0 = constrainedDentri(n, mesh)
E(graph0)$eid = c(1:ecount(graph0))  # edge id
V(graph0)$vid = c(1:vcount(graph0))  # vertex id
mstgraph = mst(graph0)  # initial spanning tree
graph0 = delete_edge_attr(graph0, 'weight')
mstgraph0 = delete_edge_attr(mstgraph, 'weight')

# plot spatial graph
plotGraph(coords, graph0) + 
  geom_boundary(bnd) +
  labs(x = 'Scaled Lon.', y = 'Scaled Lat.') + 
  ggtitle('Spatial Graph')

M = 30      # number of weak learners
k_max = 5   # maximum number of clusters per weak learner
mu = list() # initial values of mu
mstgraph_lst = list()  # initial spanning trees
cluster = matrix(1, nrow = n, ncol = M)  # initial cluster memberships
for(m in 1:M) {
  mu[[m]] = c(0)
  mstgraph_lst[[m]] = mstgraph0
}

init_val = list()
init_val[['trees']] = mstgraph_lst
init_val[['mu']] = mu
init_val[['cluster']] = cluster
init_val[['sigmasq_y']] = 1

# standardize Y
std_res = standardize(Y)
Y_std = std_res$x; std_par = std_res$std_par

# find lambda_s
nu = 3; q = 0.9
quant = qchisq(1-q, nu)
lambda_s = quant * var(Y_std) / nu

hyperpar = c()
hyperpar['sigmasq_mu'] = (0.5/(2*sqrt(M)))^2
hyperpar['lambda_s'] = lambda_s
hyperpar['nu'] = nu
hyperpar['lambda_k'] = 4
hyperpar['M'] = M
hyperpar['k_max'] = k_max

# MCMC parameters
# number of posterior samples = (MCMC - BURNIN) / THIN
MCMC = 30000    # MCMC iterations
BURNIN = 15000  # burnin period length
THIN = 5        # thinning intervals


### Run BAST ----
# warning: this step is time consuming
# for quick results, try smaller MCMC and BURNIN
# API documentations can be found in BASTFun.R

# train BAST
mcmc_res = fitBAST(Y_std, graph0, init_val, hyperpar, MCMC, BURNIN, THIN, seed = 1234)

# soft prediction
Y_grid_all = predictBAST(mcmc_res, coords, coords_grid, method = 'soft-mesh', mesh = mesh, 
                         weighting = 'uniform', return_type = 'all', seed = 12345)
Y_grid_all = apply(Y_grid_all, 2, unstandardize, std_par = std_par)

# use posterior mean as predictor
Y_grid = rowMeans(Y_grid_all)

# back to original scale
Y_grid = Y_grid + mean_Y

# plot predictive field
plotField(coords_grid, Y_grid, title = 'Predictive Field') +
  xlab('Scaled Long.') + ylab('Scaled Lat.')


