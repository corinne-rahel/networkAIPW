# save the simulation setting
save_setting <- function(nrep_graph, typeofgraph, error.type, model, 
                         prob, const, prob.rewiring, sigma_Y, sigma_C, K, S, nvals, 
                         eta, prename, 
                         num.trees, min.node.size,
                         fun_model, feat_foo, conf_model, g1, g0, 
                         subdir = "Results", get_sigma2, power_foo, cluster_Ik, 
                         num.nei) {
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
  
  setup <- list(nrep_graph = nrep_graph, 
                typeofgraph = typeofgraph, 
                error.type = error.type, 
                model = model, 
                prob = prob, 
                const = const, 
                prob.rewiring = prob.rewiring, 
                num.nei = num.nei,
                sigma_Y = sigma_Y, 
                sigma_C = sigma_C, 
                K = K, S = S, 
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
                power_foo = power_foo)
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
  # make sure return value is positive ???
  # if (var_est >= 0) {
  #   var_est
  # } else {
  #   mean(psi) ^ 2 / nval
  # }
  var_est
}

# perform simulation
simulation_foo <- function(typeofgraph, nval, mc.cores_nrep, nrep_graph, 
                           error.type, sigma_C, eta, sigma_Y, S, K, 
                           fun_model, num.trees, 
                           prob, const, prob.rewiring, seeds, power_foo, cluster_Ik) {
  
  ########################################################################################
  ############################### Data Generation ####################################
  ########################################################################################
  
  if (typeofgraph == "rand_npfix"){
    A <- A_random_graph(nval, const / nval * power_foo(nval))
    A <- Matrix(A, sparse=T)
  } else if (typeofgraph == "WS"){
    A <- A_Watts_Strogatz_new(nval, prob.rewiring, num.nei * round(power_foo(nval)))
    A <- as.matrix(A)
    A <- Matrix(A, sparse = T)
  } else {print("graph type not known")}
  
  # plot(graph_from_adjacency_matrix(adjmatrix = A, mode = "undirected"))
  
  # define dependency adjacency matrix (corresponding to dependency graph) 
  # for model 1 with the two features this is given by: 
  A.all.connections <- A %*% A+ A
  # Diagonale raus, da man neighbors bestimmen will
  diag(A.all.connections) <- 0
  # remove double edges
  A.all.connections[A.all.connections > 0] <- 1
  #plot(graph_from_adjacency_matrix(adjmatrix = A.all.connections, mode = "undirected"))
  if (FALSE) {
    par(mfrow = c(1,2))
    col <- rep("white", nval)
    col[inds_Ikc] <- "yellow"
    col[Ik] <- "pink"
    plot(graph_from_adjacency_matrix(adjmatrix = A, mode = "undirected"), 
         vertex.color = W)
    # vertex.color = to set color of vertices
    plot(graph_from_adjacency_matrix(adjmatrix = A.all.connections, mode = "undirected"), 
         vertex.color = col)
  }
  
  nval.eff <- nrow(A) #effective nval (not always possible to construct certain graph of size nval)
  
  # inds_others.all[[i]] = all indices j such j--i in the dependency graph
  inds_others.all <- lapply(seq_len(nval.eff), function(i){
    inds <- which(A.all.connections[i, ] >= 1)
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
  
  mclapply(1:nrep_graph, function(iter){
    #print(paste0("graph nr: ", iter))
    set.seed(seeds[iter])
    
    ########################################################################################
    ############################### Data Generation ####################################
    ########################################################################################
    
    C <- if(error.type == "rnorm"){
      rnorm(nval.eff, 0, sigma_C)
    } else if (error.type == "runif"){
      # make sure that the sd is sigma_C, that is sqrt(1 / 12 * (b - a) ^ 2)
      runif(nval.eff, 0, sqrt(12 * sigma_C ^ 2))
    } else if (error.type == "rt"){
      # make sure that the sd is sigma_C, that is sqrt( df / (df - 2))
      rt(nval.eff, df = (2 * sigma_C ^ 2 / (sigma_C ^ 2 - 1)))
    }
    
    sigm_eval <- sigm(eta * (C - 0.25))  
    W <- rbinom(n = nval.eff, size = 1, prob = sigm_eval)
    
    table_W <- table(W)
    num0 <- table_W[1]
    num1 <- table_W[2]
    
    dat <- 
      fun_model(A = A, W = W, g1 = g1, g0 = g0, C = C, 
                sigma_Y = sigma_Y, error.type = error.type)
    
    
    ######################################################################################
    ############################# Estimation #############################################
    ######################################################################################
    
    theta.hat.S <- matrix(0, nrow = S, ncol = K)
    theta.IPW.S <- matrix(0, nrow = S, ncol = K)
    sigma2.hat.S <- rep(0, S)
    len_Ikc <- matrix(0, nrow = S, ncol = K)
    for (s in 1:S) {
      # reshuffle data 
      inds.dat <- sample(nval.eff, replace = F)
      # num.dat.per.fold[k] = number of datapoints i belonging to kth fold
      num.dat.per.fold <- rep(floor(nval.eff / K), K)
      num.dat.per.fold[seq_len(nval.eff - K * floor(nval.eff / K))] <- 
        num.dat.per.fold[seq_len(nval.eff - K * floor(nval.eff / K))] + 1
      
      # initialize data containers
      theta.hat.k <- rep(0, K)
      theta.IPW.hat.k <- rep(0, K)
      pred.g1.all <- rep(NA_real_, nval.eff)
      pred.g0.all <- rep(NA_real_, nval.eff)
      pred.h.all <- rep(NA_real_, nval.eff)
      
      inds.Ik.choose_from <- seq(nval.eff)
      
      for (k in 1:K){
        Ik <- if (cluster_Ik) {
          # choose indices Ik as randomly as possible and such that they 
          # are connected by a tie in the dependency graph
          Ik.our_selection <- vector(mode = "numeric")
          Ik.choose_from <- vector(mode = "numeric")
          while (length(Ik.our_selection) < num.dat.per.fold[[k]]) {
            if (is_empty(Ik.choose_from)) {
              Ik.update <- sample(inds.Ik.choose_from, 1)
              Ik.choose_from <- c(Ik.update, inds_others.all[[Ik.update]])
              Ik.our_selection <- c(Ik.our_selection, Ik.update)
            }
            
            Ik.update <- Ik.choose_from[1]
            tmp <- setdiff(inds_others.all[[Ik.update]], Ik.our_selection)
            Ik.our_selection <-
              c(Ik.our_selection, tmp)
            Ik.our_selection <-
              Ik.our_selection[seq_len(min(length(Ik.our_selection), num.dat.per.fold[k]))]
            Ik.choose_from <- Ik.choose_from[-1]
            Ik.choose_from <- unique(c(Ik.choose_from, tmp))
          }
          inds.Ik.choose_from <- setdiff(inds.Ik.choose_from, Ik.our_selection)
          
          Ik.our_selection
        } else {
          inds.dat[(1 + sum(num.dat.per.fold[seq_len(k - 1)])):(sum(num.dat.per.fold[seq_len(k)]))]
        }
        
        S_Ik <- data.frame(Y = dat$Y[Ik], 
                           featX1 = dat$featX1[Ik], 
                           C = C[Ik], 
                           W = W[Ik])
        
        inds_Ikc <- setdiff(1:nval.eff, unique(c(do.call(c, inds_others.all[Ik]), Ik)))
        len_Ikc[s, k] <- length(inds_Ikc)
        S_Ik_complement <- 
          data.frame(Y = dat$Y[inds_Ikc], 
                     featX1 = dat$featX1[inds_Ikc],
                     C = C[inds_Ikc], 
                     W = W[inds_Ikc])
        
        if (FALSE) {
          par(mfrow = c(1,2))
          col <- rep("white", nval.eff)
          col[inds_Ikc] <- "yellow"
          col[Ik] <- "pink"
          plot(graph_from_adjacency_matrix(adjmatrix = A, mode = "undirected"), 
               vertex.color = col)
          # vertex.color = to set color of vertices
          plot(graph_from_adjacency_matrix(adjmatrix = A.all.connections, mode = "undirected"), 
               vertex.color = col)
        }
        
        
        # estimation of nuisance parameters on complement(S_Ik)
        # remove the W from the fit because we subset already to W==1 or W==0
        g1.est <-
          ranger(formula = Y~.,
                 data =  S_Ik_complement[S_Ik_complement$W == 1, ][, colnames(S_Ik_complement) != "W"],
                 num.trees = num.trees,
                 min.node.size = min.node.size
          )
        
        g0.est <-
          ranger(formula = Y~.,
                 data =  S_Ik_complement[S_Ik_complement$W == 0, ][, colnames(S_Ik_complement) != "W"],
                 num.trees = num.trees,
                 min.node.size = min.node.size
          )
        
        dat.h.est <- #S_Ik_complement[, c("W", "C")]
          data.frame(W = factor(S_Ik_complement$W, levels = c(0, 1)), 
                     C = S_Ik_complement$C)
        
        h.est <- 
          ranger(probability = TRUE,
                 class.weights = table(dat.h.est$W) / length(dat.h.est$W),
                 formula = W ~ ., data = dat.h.est ,
                 num.trees = num.trees, 
                 min.node.size = min.node.size
          )
        
        pred.g1.all[Ik] <- 
          predict(g1.est, 
                  data = S_Ik[, colnames(S_Ik) != "W"])$predictions
        pred.g0.all[Ik] <- 
          predict(g0.est, 
                  data = S_Ik[, colnames(S_Ik) != "W"])$predictions
        
        pred.h.all[Ik] <- 
          predict(h.est, data = data.frame(W = factor(S_Ik$W, levels = c(0, 1)), C = S_Ik$C)
          )$predictions[, "1"]
        
        inds_zero <- S_Ik$W == 0
        inds_one <- !inds_zero
        part_one <- S_Ik$W / pred.h.all[Ik] * (S_Ik$Y - pred.g1.all[Ik])
        part_one[inds_zero] <- 0
        part_zero <- (1 - S_Ik$W) / (1 - pred.h.all[Ik]) * (S_Ik$Y - pred.g0.all[Ik])
        part_zero[inds_one] <- 0
        theta.hat.Ik <- pred.g1.all[Ik] - pred.g0.all[Ik] + part_one - part_zero
        
        # potentially take median here instead of mean
        theta.hat.k[k] <- mean(theta.hat.Ik[!is.infinite(theta.hat.Ik) & !is.na(theta.hat.Ik)])
        
        # IPW estimator
        IPW_one <- S_Ik$W * S_Ik$Y / pred.h.all[Ik]
        IPW_one[inds_zero] <- 0
        IPW_zero <- (1 - S_Ik$W) * S_Ik$Y / (1 - pred.h.all[Ik])
        IPW_zero[inds_one] <- 0
        IPW_onezero <- IPW_one - IPW_zero
        theta.IPW.hat.k[k] <- 
          mean(IPW_onezero[!is.infinite(IPW_onezero) & !is.na(IPW_onezero)])
      } # end for k in 1:K
      
      # theta.hat.S[s] <- mean(theta.hat.k[!is.infinite(theta.hat.k) & !is.na(theta.hat.k)])
      theta.IPW.S[s, ] <- theta.IPW.hat.k
      theta.hat.S[s, ] <- theta.hat.k
      sigma2.hat.S[s] <- 
        get_sigma2(pred.g1.all = pred.g1.all, pred.g0.all = pred.g0.all, 
                   pred.h.all = pred.h.all, nval = nval.eff, W = W, Y = dat$Y, 
                   D = D, degrees_D = degrees_D, len_D = len_D, B = B, 
                   degrees_B = degrees_B, degrees_B_indices = degrees_B_indices)
    } # end for s in 1:S

    theta.hat.interm <- apply(theta.hat.S, 1, function(x) {xx = x[!is.na(x) & !is.infinite(x)]; mean(xx)})
    theta.hat <- median(theta.hat.interm)
    sigma2.hat <- median(sigma2.hat.S + (theta.hat.interm - theta.hat) ^ 2)
  
    # compute Horvitz-Thompson estimator
    pi.hat <- mean(W)
    theta.HT <- 
      sum(W * dat$Y) / (pi.hat * nval.eff) - sum((1 - W) * dat$Y) / ((1 - pi.hat) * nval.eff)
    
    # compute IPW under interference
    # perform clustering of the network
    kc <- cluster_fast_greedy(graph_from_adjacency_matrix(adjmatrix = A, mode = "undirected"))
    # extract membership of nodes
    members <- as.numeric(membership(kc))
    # prepare data
    data.IPC.interf <- data.frame(Y = dat$Y, C = C, W = W, group = members)
    mean_W <- mean(W)
    IPW_fit_interf <- 
      interference(formula = Y | W ~ C | group, 
                   data = data.IPC.interf, 
                   allocations = c((0 + mean_W) / 2, mean_W, (mean_W + 1) / 2), 
                   causal_estimation_options = list(variance_estimation = "naive"), 
                   model_method = "glm")
    IPW_interf <- direct_effect(IPW_fit_interf, mean_W)[, "estimate"]
    IPW_std.error <- direct_effect(IPW_fit_interf, mean_W)[, "std.error"]
    IPW_CI <- 
      c(direct_effect(IPW_fit_interf, mean_W)[, "conf.low"], 
        direct_effect(IPW_fit_interf, mean_W)[, "conf.high"])
    IPW_alpha <- direct_effect(IPW_fit_interf, mean_W)[, "alpha1"]
    
    return(list(
      nval.eff = nval.eff,
      num0 = num0,
      num1 = num1,
      theta = dat$theta, 
      theta.hat = theta.hat,
      sigma2.hat = sigma2.hat, 
      sigma2.hat.S = sigma2.hat.S, 
      theta.hat.S = theta.hat.S, 
      theta.HT = theta.HT, 
      theta.IPW.hat.S = theta.IPW.S,
      IPW_interf = IPW_interf, 
      no_components = length(kc), 
      IPW_std.error = IPW_std.error, 
      IPW_CI = IPW_CI, 
      IPW_alpha = IPW_alpha, 
      len_Ikc = len_Ikc
    )
    ) 
    
  }, mc.cores=mc.cores_nrep)
}