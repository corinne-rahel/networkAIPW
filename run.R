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

# set parameters
nrep_graph <- 1000 
typeofgraph <-"WS" #c("rand_npfix", "WS")
error.type <- "runif" 

start_time_1 <- Sys.time()

# for rand_npfix
prob <- 0.01 # probability of an edge in random graphs
const <- 3 # constant for random graphs with np = const

# for WS
prob.rewiring <- 0.05  # rewiring probability of an edge in WS
num.nei <- 2 

sigma_Y <- 0.1 # sigma for error term for Y
sigma_C <- sqrt(1 / 12) # sigma for error term for C
K <- 10 # number of sample splittings
S <- 20 # number of times the sample splitting is repeated on one generated dataset
nvals <- c(2500 / 4, 2500 / 2, 2500, 5000, 10000)
eta <- 1 
cluster_Ik <- TRUE

# rf settings
num.trees <- 500 
min.node.size <- 5 

fun_model <- get(paste0("data_model", model))

# function to compute sd
get_sigma2 <- get_sigma2
power_foo <- function(nval) {
  #nval ^ (1 / 9)
  1
}

# save the setting 
save_setting(nrep_graph = nrep_graph, 
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
             nvals = nvals, 
             eta = eta, 
             prename = prename, 
             num.trees = num.trees, 
             min.node.size = min.node.size, 
             fun_model = fun_model, 
             feat_foo = feat_X1_C, 
             conf_model = conf_model, 
             g1 = g1, 
             g0 = g0, 
             get_sigma2 = get_sigma2, 
             power_foo = power_foo, 
             cluster_Ik = cluster_Ik, 
             num.nei = num.nei)

set.seed(1)
seeds <- sample(c(0:10000), nrep_graph, replace = FALSE)

mc.cores_nrep <- 100 


for (nval in nvals) {
  print(paste0("nval: ", nval, ", ", Sys.time()))
  
  Results_nval <- 
    simulation_foo(typeofgraph = typeofgraph, nval = nval, mc.cores_nrep = mc.cores_nrep, 
                   nrep_graph = nrep_graph, error.type = error.type, sigma_C = sigma_C, 
                   eta = eta, sigma_Y = sigma_Y, S = S, K = K, 
                   fun_model = fun_model, num.trees = num.trees, 
                   prob = prob, const = const, prob.rewiring = prob.rewiring, 
                   seeds = seeds, power_foo = power_foo, cluster_Ik = cluster_Ik)
  
  filename <- paste0(prename, nval, "_", typeofgraph)
  
  filename <- gsub("[.]","",filename)
  assign(filename, Results_nval)
  save(Results_nval, file = paste(filename, ".RData",sep=""))
}

(time_elapset <- Sys.time() - start_time_1)
plot_results()
