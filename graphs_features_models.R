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

get_A_tilde <- function(A, neighbor_degree) {
  A_prod <- if (neighbor_degree == 1) {
    A
  } else if (neighbor_degree == 2) {
    A %*% A
  } else if (neighbor_degree == 3) {
    A %*% A %*% A
  }
  diag(A_prod) <- 0
  as.matrix(normalize_A(A_prod))
}

feat_X1_C <- function(A_tilde, W, C){
  W[W == 0] <- -1
  feat <- c(A_tilde %*% (W * C))
  feat[which(c(feat) == "NaN")] <- 0
  return(feat)
}

#------------------------------------------------------------------------------------------------------------
# fraction of neighbor's C's
#------------------------------------------------------------------------------------------------------------

feat_C <- function(A_tilde, W, C){
  feat <- c(A_tilde %*% C)
  feat[which(c(feat) == "NaN")] <- 0
  return(feat)
}

##############################################################################################
##############################################################################################
## New Models, not the one of the notes anymore (angepasst: more combinations of the different possible
## mechanisms, no weights for the beginning)
##############################################################################################
##############################################################################################

# data generating functions
get_C <- function(error.type, sigma_C, nval.eff) {
  if(error.type == "rnorm"){
    rnorm(nval.eff, 0, sigma_C)
  } else if (error.type == "runif"){
    # make sure that the sd is sigma_C, that is sqrt(1 / 12 * (b - a) ^ 2)
    runif(nval.eff, 0, sqrt(12 * sigma_C ^ 2))
  } else if (error.type == "rt"){
    # make sure that the sd is sigma_C, that is sqrt( df / (df - 2))
    rt(nval.eff, df = (2 * sigma_C ^ 2 / (sigma_C ^ 2 - 1)))
  }
}

h_foo <- function(eta, C, feat_Z) {
  0.15 * (C < 0.33) + 0.5 * (0.33 <= C) * (C < 0.66) + 0.85 * (0.66 <= C)
}

get_W <- function(eta, C, nval.eff, feat_Z) {
  stopifnot(ncol(feat_Z) == 1)
  feat_Z[is.na(feat_Z[, 1]), 1] <- 0
  sigm_eval <- h_foo(eta, C, feat_Z)
  rbinom(n = nval.eff, size = 1, prob = sigm_eval) 
}

get_error.Y <- function(error.type, sigma_Y, n) {
  if(error.type == "rnorm"){
    rnorm(n, 0, sigma_Y)
  } else if (error.type == "runif"){
    # make sure that the sd is sigma_Y, that is sqrt(1/12(b-a)^2)
    runif(n, -sqrt(12 * sigma_Y ^ 2) / 2, sqrt(12 * sigma_Y ^ 2) / 2)
  } else if (error.type=="rt"){
    # make sure that the sd is sigma_Y, that is sqrt(df/(df-2))
    rt(n, df = (2 * sigma_Y ^ 2 / (sigma_Y ^ 2 - 1)))
  }
}

get_features_X <- function(W, C, model, A_tilde_deg1 = NULL, A_tilde_deg2 = NULL) {
  if (model == "1") {
    cbind(
      X1 = feat_X1_C(A_tilde = A_tilde_deg1, W = W, C = C)
    )
  } else if (model == "2") {
    cbind(X1 = feat_C(A_tilde = A_tilde_deg1, W = NULL, C = W), 
          X2 = feat_C(A_tilde = A_tilde_deg2, W = NULL, C = C)) 
  }
}

get_features_Z <- function(C, model) {
  if (model == "1") {
    cbind(Z1 = NA_real_)
  } else if (model == "2") {
    cbind(Z1 = feat_C(A_tilde = A_tilde_deg1, W = NULL, C = C))
  }
}

sigm <- function(x){
  return(1/(1+exp(-x)))
}

# g functions
g1 <- function(x, c) {
  1.5 * (x[1] >= 0.5) * (c >= -0.2) * (x[1] < 0.7) + 
    4 * (x[1] >= 0.5) * (c >= -0.2) * (x[1] >= 0.7) + 
    0.5 * (x[1] >= 0.5) * (c < -0.2) + 
    3.5 * (x[1] < 0.5) * (c >= -0.2) + 
    2.5 * (x[1] < 0.5) * (c < -0.2)
}

g0 <- function(x, c) {
  0.5 * (x[1] >= 0.4) * (c >= 0.2) - 0.75 * (x[1] >= 0.4) * 
    (c < 0.2) + 0.25 * (x[1] < 0.4) * (c >= 0.2) - 0.5 * 
    (x[1] < 0.4) * (c < 0.2)
}

gen_Y_conf <- function(feat, W, g1_foo, g0_foo, C, sigma_Y, error.type){
  n <- length(W)
  
  error.Y <- get_error.Y(error.type = error.type, sigma_Y = sigma_Y, n = n)
  
  g1.eval <- unname(sapply(1:n, function(i){
    #g1_2(feat[i, 1], feat[i, 2], C[i])
    g1_foo(feat[i, ], C[i])
  }))
  
  g0.eval <- unname(sapply(1:n, function(i){
    #g0_2(feat[i, 1], feat[i, 2], C[i])
    g0_foo(feat[i, ], C[i])
  }))
  
  Y <- W * g1.eval + (1 - W) * g0.eval + error.Y  
  
  theta <- mean(g1.eval) - mean(g0.eval)
  
  return(list(Y = Y, theta = theta, error.Y = error.Y
  ))
}

data_model <- function(W, g1_foo, g0_foo, C, sigma_Y, error.type, feat_X, feat_Z){
  dat <- 
    gen_Y_conf(feat = feat_X, W = W, g1_foo = g1_foo, g0_foo = g0_foo, C = C, 
               sigma_Y = sigma_Y, error.type = error.type)
  return(list(Y = as.numeric(dat$Y), feat_X = cbind(feat_X), feat_Z = cbind(feat_Z), 
              theta = dat$theta, error.Y = dat$error.Y))
}

get_theta0N <- function(error.type, sigma_C, nval.eff, eta, C, 
                        sigma_Y, 
                        nsim, mc.cores_nsim, nsim_inner, model, g1_foo, g0_foo, 
                        A_tilde_deg1) {
  seeds <- 1:nsim
  res <- mclapply(seq_len(nsim), function(nsims) {
    set.seed(seeds[nsims])
    
    C <- get_C(error.type = error.type, sigma_C = sigma_C, nval.eff = nval.eff)
    feat_Z <- get_features_Z(C = C, model = model)
    W <- get_W(eta = eta, C = C, nval.eff = nval.eff, feat_Z = feat_Z)
    feat_X <- get_features_X(W = W, C = C, model = model, A_tilde_deg1 = A_tilde_deg1)
    eps_Y <- get_error.Y(error.type = error.type, sigma_Y = sigma_Y, n = nval.eff)
    
    g1.eval <- 
      unname(sapply(seq_len(nval.eff), function(i){
        g1_foo(feat_X[i, ], C[i])
      }))
    g0.eval <- 
      unname(sapply(seq_len(nval.eff), function(i){
        g0_foo(feat_X[i, ], C[i])
      }))
    h.eval <- sapply(seq_len(nval.eff), function(i) {
      get_W(eta = eta, C = rep(C[i], nsim_inner), nval.eff = nsim_inner, feat_Z = feat_Z) %>%
        mean()
    })
    
    theta0N <- mean(g1.eval - g0.eval) #/ (nval.eff)
    correct_one <- W / h.eval * eps_Y
    correct_one[W == 0 | is.infinite(correct_one)] <- 0
    correct_zero <- (1 - W) / (1 - h.eval) * eps_Y
    correct_zero[W == 1 | is.infinite(correct_zero)] <- 0
    theta0N_full <- mean(g1.eval - g0.eval + correct_one - correct_zero)
    list(theta0N = theta0N, theta0N_full = theta0N_full)
  }, mc.cores = mc.cores_nsim)  %>%
    do.call(what = rbind)
  
  return(list(theta0N = res[, "theta0N"] %>% do.call(what = c) %>% mean(), 
              var0N = res[, "theta0N_full"] %>% do.call(what = c) %>% var(), 
              var0N_new = res[, "theta0N"] %>% do.call(what = c) %>% var()))
}
