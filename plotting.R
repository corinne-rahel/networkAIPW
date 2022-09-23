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
                         final.graph = NULL) {
  if (is.null(final.graph)) {
    if (is.null(folders)) { 
      folders <- tail(str_split(getwd(), "/")[[1]], 1)
    } else {
      setwd(folders[1])
    }
  } 
  
  # set colors
  myColors <- c("#117DE1", "#FC8D62", "#A6D854", "red", "cyan") # "#E6AB02"
  names(myColors) <- c("netAIPW", "Hajek", "IPW", "netAIPW_blow", "netAIPW_emp") # "IPW interference"
  if (!show_empirical) {
    myColors <- myColors[-6]
  }
  if (!show_blow) {
    myColors <- myColors[-5]
  }
  
  # set line type
  myLTY <- c(1, 3)
  #names(myLTY) <- c("Erdos-Renyi", "Watts-Strogatz")
  names(myLTY) <- c("const", "N ^ (1 / 9)")
  
  my_alpha <- 0.15
  my_cex <- 1.5
  my_lwd <- 0.5 #0.75
  my_lwd_small <- 0.5
  textsize <- 10
  hadjust <- 0.5
  no_lines <- 1
  
  load("preliminary.RData")
  nrep_graph <- setup$nrep_graph
  nvals <- setup$nvals
  z <- qnorm(1 - alpha / 2, 0, 1)
  K <- if (do_product) {
    1
  } else {
    setup$K
  }
  S <- setup$S
  
  res <- tibble()
  for (dir in folders) {
    setwd(paste0("../", dir))
    all_files <- list.files()
    load("preliminary.RData")
    typeofgraph <- setup$typeofgraph
    nvals <- setup$nvals
    power_foo <- setup$power_foo
    power_foo <- as.character(body(power_foo))[2]
    power_foo <- if(power_foo == "nval^(1/9)") {
      "N ^ (1 / 9)"
    } else if (power_foo == "1") {
      #power_foo
      "const"
    }
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
      print(n)
      file_to_load <- 
        grep(paste0("_", n, "_", setup$typeofgraph, ".RData"), all_files, value = TRUE)
      stopifnot(length(file_to_load) == 1)
      load(file_to_load)
      
      sigma2.hat.S.help <- do.call(rbind, Results_nval)[, "sigma2.hat.S"]
      a <- sapply(seq_len(nrep_graph), function(i) length(sigma2.hat.S.help[[i]]))
      delete_ind <- which(a == 1)
      
      sigma2.hat.S <- do.call(rbind, sigma2.hat.S.help[setdiff(1:nrep_graph, delete_ind)])
      theta.hat.S <- 
        t(t(sapply(setdiff(c(1:nrep_graph), delete_ind), function(i) {
          rowSums(do.call(rbind, Results_nval)[, "theta.hat.S"][[i]]) / K
        })))
      theta.IPW.hat.S <- 
        t(t(sapply(setdiff(c(1:nrep_graph), delete_ind), function(i) {
          rowSums(do.call(rbind, Results_nval)[, "theta.IPW.hat.S"][[i]]) / K
        })))
      
      sigma2.aggregate <-
        sapply(seq_len(nrep_graph - length(delete_ind)), function(i) {
          median(sigma2.hat.S[i, ][sigma2.hat.S[i, ] >= 0])
        })
      theta.median <- if (S == 1) {
        apply(theta.hat.S, 1, median)
      } else {
        apply(theta.hat.S, 2, median)
      }
      theta.IPW.median <- if (S == 1) {
        apply(theta.IPW.hat.S, 1, median)
      } else {
        apply(theta.IPW.hat.S, 2, median)
      }
      sd.median <- sd(theta.median)
      sigma.blow <- sqrt(do.call(c, do.call(rbind, Results_nval)[setdiff(1:nrep_graph, delete_ind), "sigma2.hat"]))
      
      # compute p-value from aggregated method
      a <- pnorm(abs(theta.median), sd = sqrt(sigma2.hat.S), lower.tail = FALSE)
      pval.median.aggregate <- 
        2 * sapply(seq_len(nrep_graph - length(delete_ind)), function(i) median(a[i, !is.na(a[i, ])]))
      
      # compute confidence interval from aggregated method
      alpha.tilde <- 0.05
      gamma <- qnorm(1 - alpha.tilde / 4)
      g <- function(x, th, sig) {
        median(abs(th - x) / sig) - gamma
      }
      no_violations <- 0
      takeout_inds <- c()
      endpoints <- sapply(seq_len(nrep_graph - length(delete_ind)), function(i) {
        #print(i)
        valid_inds <- which(sigma2.hat.S[i, ] >= 0)
        theta_valid_inds <- theta.hat.S[, i][valid_inds]
        sig_valid_inds <- sqrt(sigma2.hat.S[i, ][valid_inds])
        gg <- function(x) {
          g(x = x, th = theta_valid_inds, sig = sig_valid_inds)
        }
        search_interval <- 
          if (length(valid_inds) > 0) {
            median(theta_valid_inds) + median(sig_valid_inds) * c(-10, 10)
          } else {
            takeout_inds <<- c(takeout_inds, i)
            median(theta_valid_inds) + c(-10, 10)
          }
        roots <- if (length(valid_inds) > 0) {
          uniroot.all(Vectorize(gg), search_interval)
        } else {
        }
        if (length(roots) > 2) {
          no_violations <<- no_violations + 1
        }
        # only proceed if at least one zero was found
        if (length(roots) >= 1) {
          # If the CI consists of several intervals
          if (length(roots) >= 2) {
            #c(min(roots), max(roots))
            mid_intervals <- sapply(seq_len(length(roots) - 1), function(k) {
              midpoint <- (roots[k] + roots[k + 1]) / 2
              if (gg(midpoint) < 0) {
                if (length(roots) == 2) {
                  list(c(roots[k], roots[k + 1]))
                } else {
                  c(roots[k], roots[k + 1])
                }
              }
            })
            # check if interval (-Inf, min(roots)) is valid
            if (gg(min(roots) - 1) < 0) {
              mid_intervals <- c(mid_intervals, list(c(-Inf, min(roots))))
            }
            # check if interval (max(roots), Inf) is valid
            if (gg(max(roots) + 1) < 0) {
              mid_intervals <- c(mid_intervals, list(c(max(roots), Inf)))
            }
            #if (length(mid_intervals)>1) {print(i)}
            mid_intervals
          } else if (length(roots) == 1) {
            if (gg(roots - 1) < 0) { # roots marks right endpoint of interval
              c(-Inf, roots)
            } else { # roots marks left endpoint of interval
              c(roots, Inf)
            }
          }
        }
        #roots
      })
      print(no_violations)
      keep_inds <- setdiff(seq_len(nrep_graph - length(delete_ind)), takeout_inds)
      CI_len_aggregate <- sapply(seq_len(nrep_graph - length(delete_ind)), function(l) {
        endpt <- endpoints[[l]]
        if (!is.list(endpt)) {
          endpt <- list(endpt)
        }
        len <- 0
        for (j in seq_len(length(endpt))) {
          if (!is.null(endpt[[j]])) {
            len <- len + endpt[[j]][2] - endpt[[j]][1]
          }
        }
        len
      })
      theta0N.est <- mean(unlist(do.call(rbind, Results_nval)[setdiff(1:nrep_graph, delete_ind), "theta"]))
      theta0_is_in_aggregate <- sapply(seq_len(nrep_graph - length(delete_ind)), function(l) {
        endpt <- endpoints[[l]]
        if (!is.list(endpt)) {
          endpt <- list(endpt)
        }
        is_in <- 0
        for (j in seq_len(length(endpt))) {
          if (!is.null(endpt[[j]])) {
            is_in <- is_in + ((endpt[[j]][1] <= theta0N.est) & (theta0N.est <= endpt[[j]][2]))
          }
        }
        is_in > 0
      })
      
      # compute Horvitz-Thompson estimator
      theta.HT <- do.call(c, do.call(rbind, Results_nval)[setdiff(1:nrep_graph, delete_ind), "theta.HT"])
      # extract IPW estimator for interference
      theta.IPW.interf <- -do.call(c, do.call(rbind, Results_nval)[setdiff(1:nrep_graph, delete_ind), "IPW_interf"])
      
      new_res <- 
        tibble(theta = theta.median[keep_inds], 
               sigma = sqrt(pmax(sigma2.aggregate, 0))[keep_inds], 
               method = "netAIPW", 
               sigma_type = "aggregate", 
               theta_type = "median", 
               pvalue = pval.median.aggregate[keep_inds]) %>%
        add_row(theta = theta.median, 
                sigma = sigma.blow, 
                method = "netAIPW_blow", 
                sigma_type = "blow", 
                theta_type = "median", 
                pvalue = 2 * pnorm(abs(theta.median), sd = sigma.blow, lower.tail = FALSE)) %>%
        add_row(theta = theta.median, 
                sigma = sd.median, 
                method = "netAIPW_emp", 
                sigma_type = "empirical", 
                theta_type = "median", 
                pvalue = 2 * pnorm(abs(theta.median), sd = sd.median, lower.tail = FALSE)) %>%
        add_row(theta = theta.HT, 
                sigma = sd(theta.HT), 
                method = "Hajek", 
                sigma_type = "emp.HAJ", 
                theta_type = "Hajek", 
                pvalue = 2 * pnorm(abs(theta.HT), sd = sd(theta.HT), lower.tail = FALSE)) %>%
        add_row(theta = theta.IPW.median, 
                sigma = sd(theta.IPW.median), 
                method = "IPW", 
                sigma_type = "emp.IPW", 
                theta_type = "IPW", 
                pvalue = 2 * pnorm(abs(theta.IPW.median), sd = sd(theta.IPW.median), lower.tail = FALSE)) %>%
        mutate(N = n, 
               graph = long_graphname, 
               power_foo = power_foo, 
               CI.lower = theta - z * sigma, 
               CI.upper = theta + z * sigma, 
               theta0N = mean(do.call(c, do.call(rbind, Results_nval)[setdiff(1:nrep_graph, delete_ind), "theta"])),
               is.in = (CI.lower <= theta0N) & (theta0N <= CI.upper), 
               CI_len = CI.upper - CI.lower, 
               bias = theta - theta0N)
      
      new_res[new_res$method == "netAIPW", "CI_len"] <- 
        CI_len_aggregate[keep_inds]
      new_res[new_res$method == "netAIPW", "is.in"] <- 
        theta0_is_in_aggregate[keep_inds]
      
      new_res <- new_res %>%
        group_by(N, method) %>%
        filter(!is.na(sigma)) %>%
        mutate(coverage = mean(is.in), 
               sigma_mean = mean(sigma), 
               theta_mean = mean(theta), 
               CI_len_med = median(CI_len), 
               bias_med = median(bias)) %>%
        ungroup() %>%
        mutate(center_scale = (theta - theta0N) / sigma_mean) %>%
        add_row(theta = mean(do.call(c, do.call(rbind, Results_nval)[setdiff(1:nrep_graph, delete_ind), "theta"])), 
                sigma_type = "0", 
                pvalue = NA, 
                theta_type = "truth", 
                method = "0",
                graph = long_graphname, 
                theta_mean = mean(do.call(c, do.call(rbind, Results_nval)[setdiff(1:nrep_graph, delete_ind), "theta"])), 
                N = n, 
                power_foo = power_foo)
      
      res <- rbind(res, new_res)
    }
  }
  
  if (!show_empirical) {
    res <- res %>%
      filter(method != "netAIPW_emp")
  }
  if (!show_blow) {
    res <- res %>%
      filter(method != "netAIPW_blow")
  }
  
  res <- res %>%
    mutate(method = factor(method, levels = c("netAIPW", "Hajek", "IPW", "0", "netAIPW_blow", "netAIPW_emp"), ordered = TRUE), 
           graph = factor(graph, levels = c("Erdos-Renyi", "Watts-Strogatz"), ordered = TRUE), 
           power_foo = factor(power_foo, levels =  c("const", "N ^ (1 / 9)"), ordered = TRUE))
  
  res_cov <- res %>%
    group_by(power_foo, N, method) %>%
    slice(1) %>%
    select(power_foo, method, sigma_type, theta_type, N, theta0N, coverage, sigma_mean, theta_mean) %>%
    filter(theta_type != "truth") %>%
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
    geom_jitter(aes(y = coverage, col = method, pch = power_foo), cex = my_cex, width = 0.02, height = 0.01) + 
    scale_colour_manual(name = "Method", values = myColors) +
    scale_fill_manual(name = "Method", values = myColors) +
    scale_linetype_manual(name = "Network growth", values = myLTY) + 
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
    labs(color = "Method", fill = "Method", lty = "Network growth", pch = "Network growth",
         x = quote(N), title = "Coverage", y = NULL)
  
  dat_pval <- res %>%
    group_by(power_foo, N, method) %>%
    slice(1) %>%
    ungroup() %>%
    select(N, method, power_foo, sigma_type, theta_type, pvalue) %>%
    filter(theta_type != "truth")
  p_pval <- dat_pval %>%
    group_by(power_foo, N, method) %>%
    ggplot(aes(x = N)) +
    geom_line(aes(y = pvalue, col = method, lty = power_foo), lwd = my_lwd) +
    #geom_jitter(aes(y = coverage, col = datafun, pch = method), cex = my_cex, width = 0, height = 0.000165) +
    geom_point(aes(y = pvalue, col = method, pch = power_foo), cex = my_cex) + 
    #facet_grid(cov_splitter~., scales="free_y") +
    #ylim(ylim_cov) +
    scale_colour_manual(name = "Method", values = myColors) +
    scale_fill_manual(name = "Method", values = myColors) +
    scale_linetype_manual(name = "Network growth", values = myLTY) + 
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
    labs(color = "Method", fill = "Method", lty = "Network growth", pch = "Network growth",
         x = quote(N), title = "p-value", y = NULL)
  
  dat_len <- res %>%
    group_by(power_foo, N, method) %>%
    slice(1) %>%
    ungroup() %>%
    select(N, method, power_foo, sigma_type, theta_type, CI_len_med) %>%
    filter(theta_type != "truth") %>%
    mutate(sqrt_N = 1 / sqrt(N) * 50)
  p_len <-
    dat_len %>%
    ggplot(aes(x = log10(N))) +
    scale_y_continuous(trans = 'log10') +
    geom_line(aes(y = CI_len_med, col = method, lty = power_foo), lwd = my_lwd) +
    geom_point(aes(y = CI_len_med, col = method, pch = power_foo), cex = my_cex) +
    scale_color_manual(name = "Method", values = myColors) + 
    scale_fill_manual(name = "Method", values = myColors) + 
    scale_linetype_manual(name = "Network growth", values = myLTY) + 
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
          #panel.spacing = unit(0.5, "lines"),
          strip.text.y = element_blank(), 
          panel.spacing = unit(no_lines, "lines")) +
    labs(color = "Method", fill = "Method", lty = "Network growth", pch = "Network growth",
         x = quote(log(N)), title = "Log median CI length", y = NULL)
  
  which_IPW <- sapply(seq_len(nrow(res)), function (i) {
    if (res$theta_type[i] == "IPW") {
      0
    } else {
      1
    }
  })

  which_IPW <- factor(which_IPW, ordered = TRUE, levels = c(0, 1, 2))
  res_bias <- res %>%
    mutate(is_IPW = which_IPW) %>%
    group_by(power_foo, N, method) %>%
    slice(1) %>%
    select(is_IPW, power_foo, method, sigma_type, theta_type, N, theta0N, coverage, sigma_mean, bias_med, theta_mean) %>%
    filter(theta_type != "truth") %>%
    ungroup() 
  p_bias <-
    res_bias %>%
    ggplot(aes(x = N)) +
    geom_hline(aes(yintercept = 0), lwd = my_lwd_small) +
    geom_line(aes(y = bias_med, lty = power_foo, col = method), lwd = my_lwd) +
    geom_point(aes(y = bias_med, col = method, pch = power_foo), cex = my_cex) + #, height = 0.02, width = 0) +
    facet_grid(is_IPW ~ ., scales = "free_y") + 
    force_panelsizes(rows = c(1, 2), 
                     TRUE) + 
    scale_color_manual(name = "Method", values = myColors) + 
    scale_fill_manual(name = "Method", values = myColors) + 
    scale_linetype_manual(name = "Network growth", values = myLTY) + 
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
    labs(color = "Method", fill = "Method", lty = "Network growth", pch = "Network growth",
         x = quote(N), title = "Median bias", y = NULL)
  
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
  
}
