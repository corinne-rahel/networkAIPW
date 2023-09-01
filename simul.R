# save the simulation setting
save_setting <- function(nrep_graph, nsim, nsim_inner, typeofgraph, error.type, model, 
                         prob, const, prob.rewiring, sigma_Y, sigma_C, K, S, nvals, 
                         eta, prename, 
                         num.trees, min.node.size,
                         data_model, feat_foo, conf_model, g1, g0, 
                         subdir = "Results", get_sigma2, power_foo, cluster_Ik, 
                         num.nei, R, fun_model, alpha,  
                         simul_file,
                         get_theta0N, 
                         simulation_foo, 
                         do_bootscheme, 
                         h_foo,
                         get_features_X) {
  now <- Sys.time()
  print(now)
  now.date <- format(now, "%d_%B_%Y")
  now.time <- format(now, "%H_%M_%S")
  dirname <- if (grepl("product", prename)) {
    paste(subdir, "/pics__", now.date, "__", now.time, "__", typeofgraph, "__product", sep = "")
  } else {
    paste(subdir, "/pics__", now.date, "__", now.time, "__", typeofgraph, sep = "")
  }
  
  if (!dir.exists(dirname)) {
    dir.create(dirname)
  }
  setwd(dirname)
  
  setup <- list(simul_file = simul_file, 
                nrep_graph = nrep_graph, 
                nsim = nsim, 
                nsim_inner = nsim_inner,
                typeofgraph = typeofgraph, 
                error.type = error.type, 
                model = model, 
                prob = prob, 
                const = const, 
                prob.rewiring = prob.rewiring, 
                num.nei = num.nei,
                sigma_Y = sigma_Y, 
                sigma_C = sigma_C, 
                K = K, S = S, R = R, 
                nvals = nvals, 
                eta = eta, 
                cluster_Ik = cluster_Ik, 
                num.trees = num.trees, 
                min.node.size = min.node.size, 
                fun_model = fun_model, 
                feat_foo = feat_foo,
                conf_model = conf_model, 
                g1 = g1, g0 = g0, 
                get_sigma2 = get_sigma2, 
                power_foo = power_foo, 
                alpha = alpha, 
                get_theta0N = get_theta0N, 
                simulation_foo = simulation_foo, 
                do_bootscheme = do_bootscheme, 
                h_foo = h_foo, 
                get_features_X = get_features_X)
  save(setup, file = "preliminary.RData")
  dput(setup, file = "preliminary.txt")
}

# compute variance
get_sigma2 <- function(pred.g1.all, pred.g0.all, pred.h.all, nval, W, Y, 
                       D, degrees_D, len_D, B, degrees_B, degrees_B_indices) {
  ind1 <- which(W == 1)
  ind0 <- which(W == 0)
  
  # initialize loss function psi and build it up
  psi <- 
    pred.g1.all - pred.g0.all #- theta.hat
  psi[ind1] <- psi[ind1] + 1 / pred.h.all[ind1] * (Y[ind1] - pred.g1.all[ind1])
  psi[ind0] <-  psi[ind0] - 1 / (1 - pred.h.all[ind0]) * (Y[ind0] - pred.g0.all[ind0])
  psi[is.na(psi) | is.infinite(psi)] <- 0
  
  # compute theta.hat for each realized degree
  theta_D <- sapply(degrees_D, function(i) {
    mean(psi[which(D == i)])
  })
  # subtract the above computed degree-dependent theta.hat
  # from the current psi to obtain the complete psi,
  # and compute successively the expectation of psi^2 
  Epsi2 <- 0
  for (i in degrees_D) {
    inds_update <- which(D == i)
    psi[inds_update] <- 
      psi[inds_update] - theta_D[which(i == degrees_D)]
    Epsi2 <- Epsi2 + 
      sum(psi[inds_update] ^ 2) / nval 
    # = mean(() ^ 2) * length(inds_update) / nval
  }
  
  # cov_hat[[l]][i, j] = Cov(psi(S_1), psi(S_2)) * A_i_j_l / nval,
  # where |X_1 intersect X_2| = l and deg(S_1) = i and deg(S_2) = j > i
  cov_hat <- 
    lapply(degrees_B, function(l){
      # get the entries of B that equal l. 
      # The first column of ind_B_deg_l represents row indices, 
      # the second column of ind_B_deg_l represents column indices. 
      ind_B_deg_l <- #which(B == l, arr.ind = TRUE)
        degrees_B_indices[[which(l == degrees_B)]]
      lapply(seq_len(len_D - 1), function(i){
        sapply(c((i + 1):len_D), function(j) {
          di <- degrees_D[i]
          dj <- degrees_D[j]
          # extract rows of ind_B_deg_l where the nodes in both columns
          # have degree di or dj
          ind_B_deg_i_j_l <- 
            ind_B_deg_l[(D[ind_B_deg_l[, 1]] == di | D[ind_B_deg_l[, 1]] == dj) & 
                          (D[ind_B_deg_l[, 2]] == di | D[ind_B_deg_l[, 2]] == dj), , drop = FALSE]
          sum(psi[ind_B_deg_i_j_l[, 1]] * psi[ind_B_deg_i_j_l[, 2]]) / nval
          
        })
      })
    })
  
  cov_hat <- unlist(cov_hat)
  
  # estimated variance
  var_est <- (Epsi2 + 2 * sum(cov_hat)) / nval
  var_est
}

# perform simulation
simulation_foo <- function(typeofgraph, nval, mc.cores_nrep, mc.cores_nsim, nrep_graph, 
                           error.type, sigma_C, eta, sigma_Y, S, K, B, 
                           data_model, num.trees, nsim, nsim_inner, 
                           prob, const, prob.rewiring, seeds, power_foo, cluster_Ik, 
                           model, R, only_netAIPW, alpha = 0.05) {
  
  ########################################################################################
  ############################### Data Generation ####################################
  ########################################################################################
  z <- qnorm(1 - alpha / 2)
  
  # generate adjacency matrix of the network
  A <- if (typeofgraph == "rand_npfix"){
    A_random_graph(nval, const / nval * power_foo(nval))
  } else if (typeofgraph == "WS"){
    A_Watts_Strogatz_new(nval, prob.rewiring, num.nei * round(power_foo(nval))) %>%
      as.matrix(A)
  } else {print("graph type not known")}
  A <- Matrix(A, sparse = T)
  nval.eff <- nrow(A) #effective nval (not always possible to construct certain graph of size nval)
  # plot(graph_from_adjacency_matrix(adjmatrix = A, mode = "undirected"))
  
  # dependency graph adjacency matrix
  identity_matrix <- diag(nrow = nval.eff)
  A.depgr <- if (model == "1") {
    #A %*% A + A
    A %*% (A + identity_matrix)
  } else if (model =="2") {
    #A %*% A %*% A %*% A + A %*% A %*% A + A %*% A + A
    A %*% (A %*% (A %*% (A + identity_matrix) + identity_matrix) + identity_matrix)
  }
  # A.depgr should only encode neighbors, so delete diagonal
  diag(A.depgr) <- 0
  # remove double edges
  A.depgr[A.depgr > 0] <- 1
  
  #plot(graph_from_adjacency_matrix(adjmatrix = A.depgr, mode = "undirected"))
  # inds_depgr_neighb[[i]] = all indices j such j--i in the dependency graph
  inds_depgr_neighb <- lapply(seq_len(nval.eff), function(i){
    inds <- which(A.depgr[i, ] >= 1)
    return(inds)
  }) # end lapply
  
  # inds_others.A[[i]] = all indices j such that Wj -> Yi
  inds_others.A <- lapply(seq_len(nval.eff), function(i){
    inds <- which(A[i, ] >= 1)
    return(inds)
  }) # end lapply
  
  # D[i] = |Xi| = number of W's in Xi 
  # = deg(i) the degree of i in the dependency graph
  D <- sapply(seq_len(nval.eff), function(i) {
    length(inds_others.A[[i]])
  })
  # compute theta.hat for each realized degree
  degrees_D <- unique(D)
  len_D <- length(degrees_D)
  
  # B[i, j] = |Xi intersect Xj| = number of W's Xi and Xj share
  # only fill upper triangular, and leave diagonal blank
  B <- Matrix(0, nrow = nval.eff, ncol = nval.eff, sparse = TRUE)
  for (i in seq_len(nval.eff - 1)) {
    ind_others_i <- c(inds_others.A[[i]], i)
    B[i, c((i + 1) : nval.eff)] <- sapply(c((i + 1) : nval.eff), function(j) {
      length(intersect(ind_others_i, c(inds_others.A[[j]], j)))
    })
  }
  # unique degee entries of B
  degrees_B <- unique(B@x)
  # degrees_B_indices[[i]] = all indices (row, col) of B such that B[row, col] = i
  degrees_B_indices <- lapply(degrees_B, function(l){
    which(B == l, arr.ind = TRUE)
  })
  
  # compute estimate of theta0N (is conditional on N and network given by
  # the adjacency matrix A)
  A_tilde_deg1 <- get_A_tilde(A = A, neighbor_degree = 1)
  A_tilde_deg2 <- get_A_tilde(A = A, neighbor_degree = 2)
  res_theta0N <-
    get_theta0N(error.type = error.type, sigma_C = sigma_C,
                nval.eff = nval.eff, eta = eta, C = C,
                sigma_Y = sigma_Y, nsim = nsim, mc.cores_nsim = mc.cores_nsim,
                nsim_inner = nsim_inner,
                model = model, g1_foo = g1_foo, g0_foo = g0_foo, A_tilde_deg1 = A_tilde_deg1)
  var0N <- res_theta0N$var0N
  var0N_new <- res_theta0N$var0N_new
  cat(paste0("ended computing theta0N at ", Sys.time(), "\n"))
  
  save(A, file = paste0("adjacencymatrix_DML_n", nval, "_", typeofgraph, ".RData"))
  
  ########################################################################################
  ############################### Select tree depth for propensity model #################
  ########################################################################################
  
  max.dep.opt <- 0
  
  get_subnet <- function(inds.Ik.choose_from, num.dat.per.fold, 
                         cluster_Ik, inds_depgr_neighb, 
                         forbidden_Ik = vector(mode = "list", length = 1)) {
    inds.Ik.choose_from_orig <- inds.Ik.choose_from
    if (cluster_Ik) {
      # choose indices Ik as randomly as possible and such that they 
      # are connected by a tie in the dependency graph
      Ik.our_selection <- vector(mode = "numeric")
      Ik.choose_from <- vector(mode = "numeric")
      while (length(Ik.our_selection) < num.dat.per.fold) {
        if (is_empty(Ik.choose_from)) {
          Ik.choose_from <- sample(inds.Ik.choose_from, 1)
        }
        missing_number_of_inds <- num.dat.per.fold - length(Ik.our_selection)
        fillup <- Ik.choose_from[seq_len(min(length(Ik.choose_from), missing_number_of_inds))]
        inds.Ik.choose_from <- setdiff(inds.Ik.choose_from, fillup)
        Ik.our_selection <- c(Ik.our_selection, fillup)
        Ik.choose_from <- c(Ik.choose_from, 
                            intersect(do.call(c, inds_depgr_neighb[fillup]), inds.Ik.choose_from_orig)) %>%
          unique() 
        Ik.choose_from <- 
          setdiff(Ik.choose_from, 
                  c(Ik.our_selection, 
                    do.call(c, forbidden_Ik)))
      }
      list(Ik.our_selection = Ik.our_selection, 
           inds.Ik.choose_from = inds.Ik.choose_from)
    } else {
      list(Ik.our_selection = sample(inds.Ik.choose_from, size = num.dat.per.fold, replace = TRUE))
    }
  }
  
  mclapply(1:nrep_graph, function(iter){
    set.seed(seeds[iter])
    
    ########################################################################################
    ############################### Data Generation ####################################
    ########################################################################################
    
    C <- runif(nval.eff, 0, sqrt(12 * sigma_C ^ 2))
    feat_Z <- cbind(Z1 = NA_real_)
    W <- rbinom(n = nval.eff, size = 1, prob = h_foo(eta, C, feat_Z = 0)) 
    feat_X <- cbind(X1 = feat_X1_C(A_tilde = A_tilde_deg1, W = W, C = C))
    error.Y <- runif(nval.eff, -sqrt(12 * sigma_Y ^ 2) / 2, sqrt(12 * sigma_Y ^ 2) / 2)
    g1.eval <- unname(sapply(1:nval.eff, function(i){
      g1_foo(feat_X[i, ], C[i])
    }))
    g0.eval <- unname(sapply(1:nval.eff, function(i){
      g0_foo(feat_X[i, ], C[i])
    }))
    theta0N <- mean(g1.eval) - mean(g0.eval)
    Y <- W * g1.eval + (1 - W) * g0.eval + error.Y 
    
    table_W <- table(W)
    num0 <- table_W[1]
    num1 <- table_W[2]
    
    
    ######################################################################################
    ############################# Estimation #############################################
    ######################################################################################
    
    #stopifnot(S == 1)
    theta.hat.S <- matrix(0, nrow = S, ncol = K)
    theta.IPW.S <- matrix(0, nrow = S, ncol = K)
    sigma2.hat.S <- rep(0, S)
    len_Ikc <- matrix(0, nrow = S, ncol = K)
    mse.g1 <- matrix(0, nrow = S, ncol = K)
    mse.g0 <- matrix(0, nrow = S, ncol = K)
    mse.h <- matrix(0, nrow = S, ncol = K)
    
    for (s in seq_len(S)) {
      # reshuffle data 
      num.dat.per.fold <- rep(floor(nval.eff / K), K)
      num.dat.per.fold[seq_len(nval.eff - K * floor(nval.eff / K))] <- 
        num.dat.per.fold[seq_len(nval.eff - K * floor(nval.eff / K))] + 1
      
      inds.Ik.choose_from <- seq(nval.eff) %>%
        sample(., size = nval.eff, replace = FALSE)
      Ik <- vector(mode = "list", length = K)
      Ik_res <- list(inds.Ik.choose_from = inds.Ik.choose_from)
      for (k in seq_len(K - 1)) {
        Ik[[k]] <-
          (Ik_res <- get_subnet(inds.Ik.choose_from = Ik_res$inds.Ik.choose_from,
                                num.dat.per.fold = num.dat.per.fold[[k]],
                                cluster_Ik = cluster_Ik,
                                inds_depgr_neighb = inds_depgr_neighb,
                                forbidden_Ik = Ik))$Ik.our_selection
      }
      Ik[[K]] <- setdiff(seq_len(nval.eff), do.call(c, Ik))
      stopifnot(do.call(c, Ik) %>% unique() %>% length() == nval.eff)
      stopifnot(do.call(c, Ik) %>% length() == do.call(c, Ik) %>% unique() %>% length())
      
      theta.hat.k.boot <- vector(mode = "list", length = K)
      theta.hat.k <- rep(0, K)
      theta.IPW.hat.k <- rep(0, K)
      var.hat.k.boot <- vector(mode = "list", length = K)
      inds_Ikc <- vector(mode = "list", length = K)
      var_formula <- lapply(seq_len(K), function(i) rep(0, nval.eff))
      pred.g1.all <- rep(0, nval.eff)
      pred.g0.all <- rep(0, nval.eff)
      pred.h.all <- rep(0, nval.eff)
      fit.g1.all <- vector(mode = "list", length = K)
      fit.g0.all <- vector(mode = "list", length = K)
      fit.h.all <- vector(mode = "list", length = K)
      for (k in 1:K){
        S_Ik <- data.frame(Y = Y[Ik[[k]]], 
                           feat_X = feat_X[Ik[[k]], ], 
                           C = C[Ik[[k]]], 
                           W = W[Ik[[k]]])
        
        inds_Ikc[[k]] <- setdiff(1:nval.eff, unique(c(do.call(c, inds_depgr_neighb[Ik[[k]]]), Ik[[k]])))
        len_Ikc[s, k] <- length(inds_Ikc[[k]])
        S_Ik_complement <-
          data.frame(Y = Y[inds_Ikc[[k]]],
                     feat_X = feat_X[inds_Ikc[[k]], ],
                     C = C[inds_Ikc[[k]]],
                     W = W[inds_Ikc[[k]]])
        
        # estimation of nuisance parameters on complement(S_Ik)
        # remove the W from the fit because we subset already to W==1 or W==0
        fit.g1.all[[k]] <-
          ranger(formula = Y~.,
                 data =  S_Ik_complement[S_Ik_complement$W == 1, ][, colnames(S_Ik_complement) != "W"],
                 num.trees = num.trees,
                 min.node.size = min.node.size#, 
                 #max.depth = 2,
          )
        fit.g0.all[[k]] <-
          ranger(formula = Y~.,
                 data =  S_Ik_complement[S_Ik_complement$W == 0, ][, colnames(S_Ik_complement) != "W"],
                 num.trees = num.trees,
                 min.node.size = min.node.size#, 
                 #max.depth = 2, 
          )
        pred.g1.all[Ik[[k]]] <-
          predict(fit.g1.all[[k]],
                  data = S_Ik[, colnames(S_Ik) != "W"])$predictions
        pred.g0.all[Ik[[k]]] <-
          predict(fit.g0.all[[k]],
                  data = S_Ik[, colnames(S_Ik) != "W"])$predictions
        
        dat.h.est <- 
          if (nrow(feat_Z) == 1 & is.na(feat_Z[1, 1])) { # no Z-features
            data.frame(W = factor(S_Ik_complement$W, levels = c(0, 1)),
                       C = S_Ik_complement$C)
          } else { # Z-features
            data.frame(W = factor(S_Ik_complement$W, levels = c(0, 1)),
                       feat_Z = feat_Z[inds_Ikc[[k]], ],
                       C = S_Ik_complement$C)
          }
        fit.h.all[[k]] <-
          ranger(probability = TRUE,
                 class.weights = table(dat.h.est$W) / length(dat.h.est$W),
                 formula = W ~ ., data = dat.h.est,
                 num.trees = num.trees,
                 min.node.size = min.node.size,
                 max.depth = 2 #max.dep.opt
          )
        dat.h.predict <-
          if (nrow(feat_Z) == 1 & is.na(feat_Z[1, 1])) { # no Z-features
            data.frame(W = factor(S_Ik$W, levels = c(0, 1)),
                       C = S_Ik$C)
          } else { # Z-features
            data.frame(W = factor(S_Ik$W, levels = c(0, 1)),
                       feat_Z = feat_Z[Ik[[k]], ],
                       C = S_Ik$C)
          }
        pred.h.all[Ik[[k]]] <-
          predict(fit.h.all[[k]], data = dat.h.predict)$predictions[, "1"]
        
        #pred.g1.all[Ik[[k]]] <- 
        pred.g1.0 <- unname(sapply(1:nrow(S_Ik), function(i){
          g1_foo(S_Ik$feat_X[[i]], S_Ik$C[i])
        }))
        #pred.g0.all[Ik[[k]]] <- 
        pred.g0.0 <- unname(sapply(1:nrow(S_Ik), function(i){
          g0_foo(S_Ik$feat_X[[i]], S_Ik$C[i])
        }))
        #pred.h.all[Ik[[k]]] <- 
        pred.h.0 <- unname(sapply(1:nrow(S_Ik), function(i){
          h_foo(eta = eta, C = S_Ik$C[i], feat_Z = 0)
        }))
        mse.g1[s, k] <- mean((pred.g0.0 - pred.g1.all[Ik[[k]]]) ^ 2)
        mse.g0[s, k] <- mean((pred.g0.0 - pred.g0.all[Ik[[k]]]) ^ 2)
        mse.h[s, k] <- mean((pred.h.0 - pred.h.all[Ik[[k]]]) ^ 2)
        
        inds_zero <- S_Ik$W == 0 
        inds_one <- !inds_zero
        part_one <- S_Ik$W / pred.h.all[Ik[[k]]] * (S_Ik$Y - pred.g1.all[Ik[[k]]])
        part_one[inds_zero] <- 0
        part_zero <- (1 - S_Ik$W) / (1 - pred.h.all[Ik[[k]]]) * (S_Ik$Y - pred.g0.all[Ik[[k]]])
        part_zero[inds_one] <- 0
        theta.hat.Ik <- pred.g1.all[Ik[[k]]] - pred.g0.all[Ik[[k]]] + part_one - part_zero
        
        # potentially take median here instead of mean
        theta.hat.k[k] <- mean(theta.hat.Ik[!is.infinite(theta.hat.Ik) & !is.na(theta.hat.Ik)])
        
        # IPW estimator
        IPW_one <- S_Ik$W * S_Ik$Y / pred.h.all[Ik[[k]]]
        IPW_one[inds_zero] <- 0
        IPW_zero <- (1 - S_Ik$W) * S_Ik$Y / (1 - pred.h.all[Ik[[k]]])
        IPW_zero[inds_one] <- 0
        IPW_onezero <- IPW_one - IPW_zero
        theta.IPW.hat.k[k] <- 
          mean(IPW_onezero[!is.infinite(IPW_onezero) & !is.na(IPW_onezero)])
        
      } # end for k in 1:K
      
      theta.hat.S[s, ] <- theta.hat.k
      theta.IPW.S[s, ] <- theta.IPW.hat.k
      sigma2.hat.S[s] <-
        get_sigma2(pred.g1.all = pred.g1.all, pred.g0.all = pred.g0.all,
                   pred.h.all = pred.h.all, nval = nval.eff, W = W, Y = Y,
                   D = D, degrees_D = degrees_D, len_D = len_D, B = B,
                   degrees_B = degrees_B, degrees_B_indices = degrees_B_indices)
      
    } # end for (s in 1:S)
    
    theta.hat.interm <- apply(theta.hat.S, 1, function(x) {xx = x[!is.na(x) & !is.infinite(x)]; mean(xx)})
    theta.hat <- median(theta.hat.interm)
    
    sigma2.hat.old <- median(sigma2.hat.S + (theta.hat.interm - theta.hat) ^ 2)
    sigma2.hat.old <- ifelse (sigma2.hat.old < 0, NA, sigma2.hat.old)
    
    # compute Horvitz-Thompson estimator
    pi.hat <- mean(W)
    theta.HT <- 
      sum(W * Y) / (pi.hat * nval.eff) - sum((1 - W) * Y) / ((1 - pi.hat) * nval.eff)
    
    
    ##############################################################
    # Bootstrap
    ##############################################################
    
    if (do_bootscheme) {
      Y.hat <- W * pred.g1.all + (1 - W) * pred.g0.all
      eps.hat <- Y - Y.hat
      eps.hat.c <- eps.hat - mean(eps.hat)
      Ik_full <- Ik
      lookup_table_Ik <- matrix(0, nrow = nval.eff, ncol = 2)
      lookup_table_Ik[, 1] <- 1:nval.eff
      lookup_table_Ik[, 2] <- sapply(1:nval.eff, function(i) {
        k <- 1
        found <- FALSE
        while(!found) {
          if (is.element(i, Ik[[k]])) {
            found <- TRUE
            return(k)
          } else {
            k <- k + 1
          }
        }
      })
      
      m <- nval.eff 
      theta.hat.boot.peter <- sapply(1:(R), function(j) {
        inds.boot <- sample(1:nval.eff, m, replace = TRUE) %>%
          sort()
        lookup_table_Ik.boot <- lookup_table_Ik[inds.boot, ]
        C.boot <- C[inds.boot]
        feat_Z.boot <- cbind(Z1 = NA_real_)
        
        g1.eval.boot <- rep(0, nval.eff)
        g0.eval.boot <- rep(0, nval.eff)
        h.eval.boot <- rep(0, nval.eff)
        for (k in 1:K) {
          inds.to.eval <- which(lookup_table_Ik.boot[, 2] == k)
          h.eval.boot[inds.to.eval] <- 
            predict(fit.h.all[[k]], 
                    data = data.frame(C = C.boot[inds.to.eval]))$predictions[, "1"]
        }
        W.boot <- rbinom(m, 1, h.eval.boot)
        feat_X.boot <- cbind(X1 = #feat_C(A_tilde = A_tilde_deg1, W = NULL, C = W.boot)
                               feat_X1_C(A_tilde = A_tilde_deg1, W = W.boot, C = C.boot))
        colnames(feat_X.boot) <- NULL
        dat.boot <- data.frame(feat_X = feat_X.boot, C = C.boot)
        for (k in 1:K) {
          inds.to.eval <- which(lookup_table_Ik.boot[, 2] == k)
          g1.eval.boot[inds.to.eval] <- 
            predict(fit.g1.all[[k]],
                    data = dat.boot[inds.to.eval, ])$predictions
          g0.eval.boot[inds.to.eval] <- 
            predict(fit.g0.all[[k]],
                    data = dat.boot[inds.to.eval, ])$predictions
          
        }
        Y.boot <- W.boot * g1.eval.boot + (1 - W.boot) * g0.eval.boot + 
          sample(eps.hat.c, m, replace = TRUE)
        
        inds.Ik.choose_from <- seq(nval.eff) %>%
          sample(., size = nval.eff, replace = FALSE)
        Ik <- vector(mode = "list", length = K)
        Ik_res <- list(inds.Ik.choose_from = inds.Ik.choose_from)
        filled <- vector()
        for (k in 1:(K - 1)) {
          Ik[[k]] <-
            (Ik_res <- get_subnet(inds.Ik.choose_from = Ik_res$inds.Ik.choose_from,
                                  num.dat.per.fold = num.dat.per.fold[[k]],
                                  cluster_Ik = cluster_Ik,
                                  inds_depgr_neighb = inds_depgr_neighb,
                                  forbidden_Ik = Ik))$Ik.our_selection
        }
        Ik[[K]] <- setdiff(seq_len(nval.eff), do.call(c, Ik))
        stopifnot(do.call(c, Ik) %>% unique() %>% length() == nval.eff)
        stopifnot(do.call(c, Ik) %>% length() == do.call(c, Ik) %>% unique() %>% length())
        
        theta.hat.k.bootscheme <- rep(0, K)
        for (k in 1:K){
          S_Ik.boot <- 
            data.frame(Y = Y.boot[Ik[[k]]], 
                       feat_X = feat_X.boot[Ik[[k]], ], 
                       C = C.boot[Ik[[k]]], 
                       W = W.boot[Ik[[k]]]
            )
          inds_Ikc <- setdiff(1:nval.eff, unique(c(do.call(c, inds_depgr_neighb[Ik[[k]]]), Ik[[k]])))
          S_Ik_complement.boot <-
            data.frame(Y = Y.boot[inds_Ikc],
                       feat_X = feat_X.boot[inds_Ikc, ],
                       C = C.boot[inds_Ikc],
                       W = W.boot[inds_Ikc])
          
          # estimation of nuisance parameters on complement(S_Ik)
          # remove the W from the fit because we subset already to W==1 or W==0
          g1.est <-
            ranger(formula = Y~.,
                   data =  S_Ik_complement.boot[S_Ik_complement.boot$W == 1, ][, colnames(S_Ik_complement.boot) != "W"],
                   num.trees = num.trees,
                   min.node.size = min.node.size#, 
                   #max.depth = 3
            )
          g0.est <-
            ranger(formula = Y~.,
                   data =  S_Ik_complement.boot[S_Ik_complement.boot$W == 0, ][, colnames(S_Ik_complement.boot) != "W"],
                   num.trees = num.trees,
                   min.node.size = min.node.size#, 
                   #max.depth = 3
            )
          pred.g1.boot <-
            predict(g1.est,
                    data = S_Ik.boot[, colnames(S_Ik.boot) != "W"])$predictions
          pred.g0.boot <-
            predict(g0.est,
                    data = S_Ik.boot[, colnames(S_Ik.boot) != "W"])$predictions
          
          dat.h.est.boot <- 
            if (nrow(feat_Z.boot) == 1 & is.na(feat_Z.boot[1, 1])) { # no Z-features
              data.frame(W = factor(S_Ik_complement.boot$W, levels = c(0, 1)),
                         C = S_Ik_complement.boot$C)
            } else { # Z-features
              data.frame(W = factor(S_Ik_complement.boot$W, levels = c(0, 1)),
                         feat_Z = feat_Z.boot[inds_Ikc, ],
                         C = S_Ik_complement.boot$C)
            }
          h.est <-
            ranger(probability = TRUE,
                   class.weights = table(dat.h.est.boot$W) / length(dat.h.est.boot$W),
                   formula = W ~ ., data = dat.h.est.boot,
                   num.trees = num.trees,
                   min.node.size = min.node.size,
                   max.depth = 2 #max.dep.opt
            )
          dat.h.predict.boot <-
            if (nrow(feat_Z.boot) == 1 & is.na(feat_Z.boot[1, 1])) { # no Z-features
              data.frame(W = factor(S_Ik.boot$W, levels = c(0, 1)),
                         C = S_Ik.boot$C)
            } else { # Z-features
              data.frame(W = factor(S_Ik.boot$W, levels = c(0, 1)),
                         feat_Z = feat_Z.boot[Ik[[k]], ],
                         C = S_Ik.boot$C)
            }
          pred.h.boot <-
            predict(h.est, data = dat.h.predict.boot)$predictions[, "1"]
          
          inds_zero <- S_Ik.boot$W == 0 
          inds_one <- !inds_zero
          part_one <- S_Ik.boot$W / pred.h.boot * (S_Ik.boot$Y - pred.g1.boot)
          part_one[inds_zero] <- 0
          part_zero <- (1 - S_Ik.boot$W) / (1 - pred.h.boot) * (S_Ik.boot$Y - pred.g0.boot)
          part_zero[inds_one] <- 0
          theta.hat.Ik.boot <- pred.g1.boot - pred.g0.boot + part_one - part_zero
          
          # potentially take median here instead of mean
          theta.hat.k.bootscheme[k] <- mean(theta.hat.Ik.boot[!is.infinite(theta.hat.Ik.boot) & !is.na(theta.hat.Ik.boot)])
          
        } # end for k in 1:K
        
        mean(theta.hat.k.bootscheme)
      })
      var.bootscheme <- var(theta.hat.boot.peter) 
    } else {
      var.bootscheme <- 0
    }
    
    list(nval.eff = nval.eff,
         num0 = 1,
         num1 = 1,
         theta0N = theta0N,
         var0N = var0N, 
         var0N_new = var0N_new, 
         theta.hat = theta.hat,
         theta.hat.S = theta.hat.S,
         theta.HT = theta.HT, 
         theta.IPW.hat.S = theta.IPW.S,
         sigma2.hat.S = sigma2.hat.S, 
         sigma2.hat = 1, 
         sigma2.hat_formula = 1,
         sigma2.hat_old = sigma2.hat.old,
         sigma2.hat_bootscheme = var.bootscheme, 
         ci.lwr = 1, 
         ci.upp = 1, 
         ci.lwr.formula = 1, 
         ci.upp.formula = 1, 
         ci.lwr.old = theta.hat - z * sqrt(sigma2.hat.old), 
         ci.upp.old = theta.hat + z * sqrt(sigma2.hat.old), 
         ci.lwr.var0N = theta.hat - z * sqrt(var0N), 
         ci.upp.var0N = theta.hat + z * sqrt(var0N), 
         ci.lwr.bootscheme = theta.hat - z * sqrt(var.bootscheme), 
         ci.upp.bootscheme = theta.hat + z * sqrt(var.bootscheme), 
         mse.g1 = rowMeans(mse.g1) %>% mean(), 
         mse.g0 = rowMeans(mse.g0) %>% mean(), 
         mse.h = rowMeans(mse.h) %>% mean(), 
         max.dep.opt = 1)
  }, mc.cores = mc.cores_nrep)
}
