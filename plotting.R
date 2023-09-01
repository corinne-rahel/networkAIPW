library("igraph")
library("Matrix")
library("parallel")
library("expm")
library("ranger")
library("Rcpp")
library("tidyverse")
library("ggplot2")
library("ggpubr")
library("wesanderson")
library("RColorBrewer")
library("rootSolve")
library("inferference")
library(ggh4x) # to force panel size

###################### This includes the p-value aggregation

xfoo=function(x){
  if ((4<x) & (x<=6)) {
    1/4*x-1
  } else if ((6<x) & (x<=8)) {
    2-0.25*x
  } else {
    0
  }
}

# compute confidence interval for estimated probability parameter (p.hat)
# of Bernoulli distribution that has been computed with rep.times
# many repetitions.
get.CIs <- function(p.hat, rep.times, z) {
  if (p.hat == 0) {
    ci.est <- c(0, 3 / rep.times)
  } else if (p.hat == 1) {
    ci.est <- c(1 - 3 / rep.times, 1)
  } else {
    ci.est <- 
      c(p.hat - z * sqrt(p.hat * (1 - p.hat) / rep.times),
        p.hat + z * sqrt(p.hat * (1 - p.hat) / rep.times))
  }
  ci.est[1] <- pmax(ci.est[1], 0)
  ci.est[2] <- pmin(ci.est[2], 1)
  ci.est
}

plot_results <- function(folders = NULL, do_product = FALSE, alpha = 0.05, 
                         show_empirical = FALSE, show_blow = FALSE, 
                         final.graph = NULL, nvals = NULL, do_mse = TRUE) {
  # folders = NULL; do_product = FALSE; alpha = 0.05; show_empirical = FALSE; show_blow = FALSE; final.graph = NULL; nvals = NULL; do_mse = TRUE
  if (is.null(final.graph)) {
    if (is.null(folders)) { 
      folders <- tail(str_split(getwd(), "/")[[1]], 1)
    } else {
      setwd(folders[1])
    }
  } 
  
  # set colors
  myColors <- c("#117DE1", "#FC8D62", "#A6D854", "#E6AB02", "red", "cyan")
  names(myColors) <- c("netAIPW", "Hajek", "IPW", "IPW interference", "netAIPW_blow", "netAIPW_emp")
  
  # set line type
  myLTY <- c(1, 3)
  my_alpha <- 0.15
  my_cex <- 1.5
  my_lwd <- 0.5 
  my_lwd_small <- 0.5
  textsize <- 10
  hadjust <- 0.5
  no_lines <- 1
  
  load("preliminary.RData")
  nrep_graph <- setup$nrep_graph
  if (is.null(nvals)) {
    nvals <- setup$nvals
    update_nvals <- TRUE
  }  else {
    update_nvals <- FALSE
  }
  z <- qnorm(1 - alpha / 2, 0, 1)
  K <- if (do_product) {
    1
  } else {
    setup$K
  }
  S <- setup$S
  
  res <- tibble()
  res_mse <- tibble()
  for (dir in folders) {
    setwd(paste0("../", dir))
    all_files <- list.files()
    load("preliminary.RData")
    typeofgraph <- setup$typeofgraph
    if (update_nvals) {
      nvals <- setup$nvals
    }
    
    if (grepl("neighb", dir)) {
      ordered_levels <- c("correct", "wrong")
      name_growth <- "Assumed spillover"
      power_foo <- if (grepl("mis", dir)) {
        "wrong"
      } else {
        "correct"
      }
    } else {
      ordered_levels <- c("const", "growing")
      name_growth <- "Network growth"
      power_foo <- setup$power_foo
      power_foo <- as.character(body(power_foo))[2]
      power_foo <- if(grepl("nval^(", power_foo, fixed = TRUE)) {
        #"N ^ (1 / 15)"
        "growing"
      } else if (power_foo == "1") {
        #power_foo
        "const"
      }
    }
    names(myLTY) <- ordered_levels #c("const", "growing")
    
    long_graphname <- 
      if (typeofgraph == "family2") {
        "Disconnected components"
      } else if (typeofgraph == "rand_npfix") {
        "Erdos-Renyi"
      } else if (typeofgraph == "WS") {
        "Watts-Strogatz"
      } else {
        stop("unknown graph type")
      }
    
    print(c(long_graphname, power_foo))
    for (n in nvals) {
      file_to_load <- 
        grep(paste0("_", n, "_", setup$typeofgraph, ".RData"), all_files, value = TRUE)
      stopifnot(length(file_to_load) == 1)
      load(file_to_load)
      Results_nval_cols <- do.call(rbind, Results_nval) 
      keep_ind <- lapply(Results_nval_cols[, "sigma2.hat_old"], function(x) !is.na(x)) %>%
        do.call(what = c)
      print(paste0(nrep_graph - sum(keep_ind), " out of ", nrep_graph, 
                   " repetitons with estimated sigma were < 0."))
      
      theta.hat <- do.call(c, Results_nval_cols[, "theta.hat"])
      
      sigma.hat_bootscheme <- do.call(c, Results_nval_cols[, "sigma2.hat_bootscheme"]) %>%
        sqrt()
      theta0N <- do.call(c, Results_nval_cols[, "theta0N"]) %>%
        mean()
      var0N <- do.call(c, Results_nval_cols[, "var0N"])
      
      ci.lwr.bootscheme <- do.call(c, Results_nval_cols[, "ci.lwr.bootscheme"])
      ci.upp.bootscheme <- do.call(c, Results_nval_cols[, "ci.upp.bootscheme"])
      
      max.dep.opt <- do.call(c, Results_nval_cols[, "max.dep.opt"])
      
      theta.IPW.hat.S <- 
        t(t(sapply(1:nrep_graph, function(i) {
          rowSums(do.call(rbind, Results_nval)[, "theta.IPW.hat.S"][[i]]) / K
        })))
      theta.IPW.median <- if (S == 1) {
        apply(theta.IPW.hat.S, 1, median)
      } else {
        apply(theta.IPW.hat.S, 2, median)
      }
      
      # compute Horvitz-Thompson estimator
      theta.HT <- do.call(c, do.call(rbind, Results_nval)[, "theta.HT"])
      
      new_res <- 
        tibble(theta.hat = theta.hat, 
               theta0N = theta0N, 
               var0N = var0N, 
               N = n, 
               graph = long_graphname, 
               power_foo = power_foo, 
               max.dep.opt = max.dep.opt, 
               sigma.hat = sigma.hat_bootscheme, 
               ci.l = ci.lwr.bootscheme, ci.u = ci.upp.bootscheme, 
               method = "netAIPW") %>%
        add_row(theta.hat = theta.HT, 
                theta0N = theta0N, 
                var0N = var0N, 
                N = n, 
                graph = long_graphname, 
                power_foo = power_foo, 
                max.dep.opt = max.dep.opt, 
                sigma.hat = sd(theta.HT), 
                ci.l = theta.HT - z * sd(theta.HT), 
                ci.u = theta.HT + z * sd(theta.HT), 
                method = "Hajek") %>%
        add_row(theta.hat = theta.IPW.median, 
                theta0N = theta0N, 
                var0N = var0N, 
                N = n, 
                graph = long_graphname, 
                power_foo = power_foo, 
                max.dep.opt = max.dep.opt, 
                sigma.hat = sd(theta.IPW.median), 
                ci.l = theta.IPW.median - z * sd(theta.IPW.median), 
                ci.u = theta.IPW.median + z * sd(theta.IPW.median), 
                method = "IPW") %>%
        group_by(method) %>%
        mutate(is.in = ci.l <= theta0N & theta0N <= ci.u, 
               coverage = mean(is.in), 
               ci.length = ci.u - ci.l, 
               CI_len_mean = mean(ci.length), 
               bias = theta0N - theta.hat, 
               mean_bias = mean(bias))
      res <- rbind(res, new_res)
      
      if (do_mse) {
        mse.g1 <- do.call(c, Results_nval_cols[, "mse.g1"])
        mse.g0 <- do.call(c, Results_nval_cols[, "mse.g0"])
        mse.h <- do.call(c, Results_nval_cols[, "mse.h"])
        
        res_mse_new <- 
          tibble(N = n, 
                 mse = c(mse.g1, mse.g0, mse.h), 
                 foo = c(rep("g1", length(mse.g1)), rep("g0", length(mse.g0)), rep("h", length(mse.h))))
        res_mse <- rbind(res_mse, res_mse_new)
      }
      
    }
  }
  
  res_max_dep <- res %>%
    group_by(N) %>%
    select(N, max.dep.opt) %>%
    slice(1)
  res_max_dep %>%
    ggplot(aes(x = N, y = max.dep.opt)) + 
    geom_point() + 
    geom_line() 
  nam <- paste0("maxdep_", long_graphname, ".pdf")
  ggsave(nam, width = 7, height = 3)
  
  
  res <- res %>%
    mutate(
      graph = factor(graph, levels = c("Erdos-Renyi", "Watts-Strogatz"), ordered = TRUE), 
      power_foo = factor(power_foo, levels =  ordered_levels, ordered = TRUE))
  
  res_cov <- res %>%
    filter(method != "netAIPW_old_FWER") %>% 
    group_by(power_foo, N, method) %>%
    slice(1) %>%
    select(-c(theta.hat, sigma.hat, is.in, ci.length, bias)) %>%
    ungroup()
  CI_cov <- sapply(seq_len(nrow(res_cov)), function(i) {
    get.CIs(as.numeric(res_cov[i, "coverage"]), nrep_graph, z)
  })
  res_cov <- res_cov %>%
    mutate(CI_cov_low = as.numeric(CI_cov[1, ]),
           CI_cov_upp = as.numeric(CI_cov[2, ]))
  
  p_cov <-
    res_cov %>%
    group_by(power_foo, N, method) %>%
    ggplot(aes(x = N)) +
    geom_ribbon(aes(ymin = CI_cov_low, ymax = CI_cov_upp, fill = method, pch = power_foo), alpha = my_alpha) +
    geom_hline(aes(yintercept = 1 - alpha), lwd = my_lwd_small) +
    geom_line(aes(y = coverage, col = method, lty = power_foo), lwd = my_lwd) +
    geom_point(aes(y = coverage, col = method, pch = power_foo), cex = my_cex) + 
    scale_colour_manual(name = "Method", values = myColors) +
    scale_fill_manual(name = "Method", values = myColors) +
    scale_linetype_manual(name = name_growth, values = myLTY) + 
    guides(colour = guide_legend(override.aes = list(alpha = 1))) +
    theme(plot.title = element_text(face = "bold", size = textsize, hjust = hadjust),
          legend.background = element_rect(fill = NA, size = 4, colour = NA),
          legend.key = element_blank(),
          plot.background = element_rect(fill = "white", colour = "white"),
          panel.background = element_rect(fill = "white", colour = "white"),
          axis.ticks = element_line(colour = "grey70", size = 0.2),
          panel.grid.major = element_line(colour = "grey70", size = 0.2),
          panel.grid.minor = element_line(colour = "grey70", size = 0.05),
          axis.text.x = element_text(angle = 90),
          panel.spacing = unit(no_lines, "lines"),
          strip.text.y = element_blank()) +
    labs(color = "Method", fill = "Method", lty = name_growth, pch = name_growth,
         x = quote(N), title = "Coverage", y = NULL)
  
  dat_len <- res %>%
    group_by(power_foo, N, method) %>%
    slice(1) %>%
    ungroup() %>%
    select(N, method, power_foo, CI_len_mean) 
  p_len <-
    dat_len %>%
    ggplot(aes(x = log10(N))) +
    scale_y_continuous(trans = 'log10') +
    geom_line(aes(y = CI_len_mean, col = method, lty = power_foo), lwd = my_lwd) +
    geom_point(aes(y = CI_len_mean, col = method, pch = power_foo), cex = my_cex) +
    scale_color_manual(name = "Method", values = myColors) + 
    scale_fill_manual(name = "Method", values = myColors) + 
    scale_linetype_manual(name = name_growth, values = myLTY) + 
    guides(colour = guide_legend(override.aes = list(alpha = 1))) +
    theme(plot.title = element_text(face = "bold", size = textsize, hjust = hadjust),
          legend.background = element_rect(fill = NA, size = 4, colour = NA),
          legend.key=element_blank(),
          plot.background = element_rect(fill = "white", colour = "white"),
          panel.background = element_rect(fill = "white", colour = "white"),
          axis.ticks = element_line(colour = "grey70", size = 0.2),
          panel.grid.major = element_line(colour = "grey70", size = 0.2),
          panel.grid.minor = element_line(colour = "grey70", size = 0.05),
          axis.text.x = element_text(angle = 90), 
          strip.text.y = element_blank(), 
          panel.spacing = unit(no_lines, "lines")) +
    labs(color = "Method", fill = "Method", lty = name_growth, pch = name_growth,
         x = quote(log(N)), title = "Log mean CI length", y = NULL)
  
  which_IPW <- sapply(seq_len(nrow(res)), function (i) {
    if (res$method[i] == "IPW" & res$power_foo[i] %in% c("growing", "correct", "wrong")) {
      0
    } else {
      1
    }
  })
  
  which_IPW <- factor(which_IPW, ordered = TRUE, levels = c(0, 1, 2, 3, 4, 5))
  res_bias <- res %>%
    ungroup() %>%
    mutate(is_IPW = which_IPW) %>%
    group_by(power_foo, N, method) %>%
    slice(1) %>%
    select(power_foo, method, N, theta0N, coverage, mean_bias, is_IPW) %>%
    ungroup() 
  p_bias <-
    res_bias %>%
    ggplot(aes(x = N)) +
    geom_hline(aes(yintercept = 0), lwd = my_lwd_small) +
    geom_line(aes(y = mean_bias, lty = power_foo, col = method), lwd = my_lwd) +
    geom_point(aes(y = mean_bias, col = method, pch = power_foo), cex = my_cex) + 
    facet_grid(is_IPW ~ ., scales = "free_y") + 
    force_panelsizes(rows = c(1, 2), 
                     TRUE) + 
    scale_color_manual(name = "Method", values = myColors) + 
    scale_fill_manual(name = "Method", values = myColors) + 
    scale_linetype_manual(name = name_growth, values = myLTY) + 
    guides(colour = guide_legend(override.aes = list(alpha = 1))) +
    theme(plot.title = element_text(face = "bold", size = textsize, hjust = hadjust),
          legend.background = element_rect(fill = NA, size = 4, colour = NA),
          legend.key=element_blank(), 
          plot.background = element_rect(fill = "white", colour = "white"),
          panel.background = element_rect(fill = "white", colour = "white"),
          axis.ticks = element_line(colour = "grey70", size = 0.2),
          panel.grid.major = element_line(colour = "grey70", size = 0.2),
          panel.grid.minor = element_line(colour = "grey70", size = 0.05),
          axis.text.x = element_text(angle = 90), 
          strip.text.y = element_blank(), 
          panel.spacing = unit(no_lines, "lines")) +
    labs(color = "Method", fill = "Method", lty = name_growth, pch = name_growth,
         x = quote(N), title = "Mean bias", y = NULL)
  
  fol <- sapply(seq_len(length(folders)), function (i) {
    a <- strsplit(folders[i], split = "__")[[1]]
    paste0(a[2], "__", a[3])
  })
  fol <- paste(fol[1], fol[2], fol[3], sep = "---")
  nam <- paste0("coverage_", long_graphname, "_", fol, ".pdf")
  
  p_final <- ggarrange(p_cov, p_len, p_bias, ncol = 3, nrow = 1, common.legend = TRUE, 
                       legend = "right")
  annotate_figure(p_final, left = text_grob(paste(long_graphname, "\n"), face = "bold", size = 17, rot = 90))
  
  
  ggsave(nam, width = 7, height = 3)
  
  if (do_mse) {
    res_mse %>%
      group_by(foo, N) %>%
      mutate(mean_mse = mean(mse)) %>%
      slice(1) %>%
      ungroup() %>%
      group_by(foo) %>%
      ggplot(aes(x = N, y = mean_mse, color = foo)) +
      #facet_wrap(~ foo) + 
      geom_point() + 
      geom_line()
    nam <- paste0("mse_", long_graphname, "_", fol, ".pdf")
    ggsave(nam, width = 7, height = 3)
  }
  
  res_sigma <- res %>%
    group_by(N, method) %>%
    mutate(sigma.hat.mean = mean(sigma.hat)) %>%
    slice(1) %>%
    select(N, graph, power_foo, sigma.hat.mean, method)
  
  res_sigma %>%
    filter(method != "netAIPW_old_FWER") %>% 
    ggplot(aes(x = log10(N), y = sigma.hat.mean, color = method)) + 
    scale_y_continuous(trans = 'log10') + 
    geom_point() + 
    geom_line() + 
    force_panelsizes(rows = c(1, 2), 
                     TRUE) + 
    scale_color_manual(name = "Method", values = myColors) + 
    scale_fill_manual(name = "Method", values = myColors) + 
    labs(color = "Method", fill = "Method", lty = name_growth, pch = name_growth,
         x = quote(N), title = "Mean bias", y = NULL)
  nam <- paste0("sigmahat_", long_graphname, "_", fol, ".pdf")
  ggsave(nam, width = 7, height = 3)
}
