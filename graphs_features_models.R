# Script storing functions that generate different graphs, features, outcome generating models

#############################################################################################################
################################### Graphs  #################################################################
#############################################################################################################
## Different Graph Models
## Generation of Adacency Matrix A
normalize_A <- function(A){
  weights <- 1 / rowSums(A)
  weights[which(weights == "Inf")] <- 0
  weight_matrix <- matrix(rep(weights, dim(A)[1]), byrow = F, ncol = dim(A)[1])
  A_norm <- A * weight_matrix
  return(A_norm)
}

## Erdoes-Renyi G(n,p) random graph
A_random_graph <- function(n, prob){
  # construct Adjacency matrix A with uniformly with bernoulli distr.
  A_entries_nonsym <- rbinom(n * n,size=1, prob = prob)
  A_matrix_nonsym <- matrix(A_entries_nonsym, ncol = n)
  sym_structure <- lower.tri(A_matrix_nonsym, diag = FALSE)
  sym_structure <- matrix(as.numeric(sym_structure), ncol = n, nrow = n)
  ind <- which(sym_structure == 0)
  A_matrix_nonsym[ind] <- 0
  A <- A_matrix_nonsym + t(A_matrix_nonsym) # this is now a symmetric adj. matrix
}

## Watts-Strogatz small-world graph
## Has size^dim many nodes
## eg. n =1000 nodes is received by dim=3, size=10
A_Watts_Strogatz_new <- function(size, p, nei){
  net <- 
    sample_smallworld(dim = 1, size = size, 
                      nei = nei, p = p, loops = FALSE, multiple = FALSE)
  A <- get.adjacency(net, type = "both",attr = NULL, names = F, sparse = T)
  return(A)
}

#############################################################################################################
########################################### Features  #######################################################
##############################################################################################################

#------------------------------------------------------------------------------------------------------------
# fraction of treated neighbors, weightd by C
#------------------------------------------------------------------------------------------------------------

feat_X1_C <- function(A, W, C){
  W[W == 0] <- -1
  A_tilde <- as.matrix(normalize_A(A))
  feat <- c(A_tilde %*% (W * C))
  feat[which(c(feat) == "NaN")] <- 0
  return(feat)
}


##############################################################################################
##############################################################################################
## New Models, not the one of the notes anymore (angepasst: more combinations of the different possible
## mechanisms, no weights for the beginning)
##############################################################################################
##############################################################################################

sigm <- function(x){
  return(1/(1+exp(-x)))
}

# g functions
g1 <- function(x, c) {
  1.5 * (x[1] >= 0.5) * (c >= -0.2) * (x[1] < 0.7) + 4 * 
    (x[1] >= 0.5) * (c >= -0.2) * (x[1] >= 0.7) + 0.5 * 
    (x[1] >= 0.5) * (c < -0.2) + 3.5 * (x[1] < 0.5) * 
    (c >= -0.2) + 2.5 * (x[1] < 0.5) * (c < -0.2)
}

g0 <- function(x, c) {
  0.5 * (x[1] >= 0.4) * (c >= 0.2) - 0.75 * (x[1] >= 0.4) * 
    (c < 0.2) + 0.25 * (x[1] < 0.4) * (c >= 0.2) - 0.5 * 
    (x[1] < 0.4) * (c < 0.2)
}

gen_Y_conf <- function(feat, W, g1, g0, C,sigma_Y, error.type){
  n <- length(W)
  
  error.Y <- if(error.type=="rnorm"){
    rnorm(n,0,sigma_Y)
  } else if (error.type == "runif"){
    # make sure that the sd is sigma_Y, that is sqrt(1/12(b-a)^2)
    runif(n, -sqrt(12*sigma_Y^2)/2, sqrt(12*sigma_Y^2)/2)
  } else if (error.type=="rt"){
    # make sure that the sd is sigma_Y, that is sqrt(df/(df-2))
    rt(n, df = (2 * sigma_Y ^ 2 / (sigma_Y ^ 2 - 1)))
  }
  
  g1.eval <- unname(sapply(1:n, function(i){
    g1(feat[i], C[i])
  }))
  
  g0.eval <- unname(sapply(1:n, function(i){
    g0(feat[i], C[i])
  }))
  
  Y <- W * g1.eval + (1 - W) * g0.eval + error.Y  
  
  theta <- mean(g1.eval) - mean(g0.eval)
  
  return(list(Y = Y, theta = theta))
  
}

# Model 1: 
data_model1 <- function(A,W, g1, g0,C, sigma_Y, error.type){
  feat <- feat_X1_C(A = A, W = W, C = C)
  dat <- gen_Y_conf(feat = feat, W = W, g1 = g1, g0 = g0, C = C, 
                    sigma_Y = sigma_Y, error.type = error.type)
  return(list(Y = as.numeric(dat$Y),theta = dat$theta,featX1 = feat))
}
