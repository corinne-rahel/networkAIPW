############################################################################
#### AIPW with ML on Networks ##############################################
############################################################################

library("igraph")
library("Matrix")
library("parallel")
library("expm")
library("ranger")
library("Rcpp")
library("tidyverse")
library("ggplot2")
library("wesanderson")
library("RColorBrewer")
library("rootSolve")
library("inferference")

source("graphs_features_models.R")
source("plotting.R")
simul_file <- "simul.R" 
source(simul_file)

# prename for file names
prename <- "DML_"
do_product <- FALSE
model <- "1"
conf_model <- gen_Y_conf
data_model <- data_model
g1_foo <- g1 
g0_foo <- g0 

# set parameters
mc.cores_nrep <- 77
mc.cores_nsim <- mc.cores_nrep
nrep_graph <- 1000 
nsim <- 5000
nsim_inner <- 1000
typeofgraph <-"WS" #c("rand_npfix", "WS")
error.type <- "runif" 

start_time_1 <- Sys.time()

# for rand_npfix
prob <- 0.01 # 0.01 probability of an edge in random graphs
const <- 3 # constant for random graphs with np = const

# for WS
prob.rewiring <- 0.05  # rewiring probability of an edge in WS
num.nei <- 2 

sigma_Y <- 0.1 # sigma for error term for Y
sigma_C <- sqrt(1 / 12) # sigma for error term for C
K <- 5 # number of sample splittings
S <- 1 # number of times the sample splitting is repeated on one generated dataset
R <- 300 # number of bootstrap repetitions 
do_bootscheme <- TRUE
nvals <- 5000 * 2 ^ c(-3:6) 
eta <- 1 
cluster_Ik <- TRUE
alpha = 0.05

# rf settings
num.trees <- 500 
min.node.size <- 5 

# function to compute sd
power_foo <- function(nval) {
  #nval ^ (1 / 9)
  #1
  1
}

# save the setting 
save_setting(simul_file = simul_file, nrep_graph = nrep_graph, 
             nsim = nsim, # for computing theta0N
             nsim_inner = nsim_inner, # for computing theta0N
             typeofgraph = typeofgraph, 
             error.type = error.type, 
             model = model, 
             prob = prob, 
             const = const, 
             prob.rewiring = prob.rewiring, 
             sigma_Y = sigma_Y, 
             sigma_C = sigma_C, 
             K = K, 
             S = S, 
             #B = B, 
             R = R, 
             nvals = nvals, 
             eta = eta, 
             prename = prename, 
             num.trees = num.trees, 
             min.node.size = min.node.size, 
             fun_model = data_model, 
             feat_foo = list(feat_X1_C = feat_X1_C, feat_C = feat_C), 
             conf_model = conf_model, 
             g1 = g1_foo, 
             g0 = g0_foo, 
             get_sigma2 = get_sigma2, 
             power_foo = power_foo, 
             cluster_Ik = cluster_Ik, 
             num.nei = num.nei, 
             alpha = alpha, 
             get_theta0N = get_theta0N, 
             simulation_foo = simulation_foo, 
             do_bootscheme = do_bootscheme, 
             h_foo = h_foo, 
             get_features_X = get_features_X)

set.seed(1)
seeds <- sample(c(0:10000), nrep_graph, replace = FALSE)


for (nval in nvals) {
  print(paste0("nval: ", nval, ", ", Sys.time()))
  
  Results_nval <- 
    simulation_foo(typeofgraph = typeofgraph, nval = nval, mc.cores_nrep = mc.cores_nrep, 
                   mc.cores_nsim = mc.cores_nsim, nsim = nsim, nsim_inner = nsim_inner,
                   nrep_graph = nrep_graph, error.type = error.type, sigma_C = sigma_C, 
                   eta = eta, sigma_Y = sigma_Y, S = S, K = K, B = B, 
                   data_model = data_model, num.trees = num.trees, 
                   prob = prob, const = const, prob.rewiring = prob.rewiring, 
                   seeds = seeds, power_foo = power_foo, cluster_Ik = cluster_Ik, 
                   model = model, R = R, only_netAIPW = TRUE, 
                   alpha = alpha)
  
  filename <- paste0(prename, nval, "_", typeofgraph)
  
  filename <- gsub("[.]","",filename)
  assign(filename, Results_nval)
  save(Results_nval, file = paste(filename, ".RData",sep=""))
}

(time_elapsed <- Sys.time() - start_time_1)
plot_results()
