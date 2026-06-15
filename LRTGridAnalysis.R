library(RTMB)


## =========================================================
## 1. Hjelpefunksjoner for miksfordelingen
## =========================================================

pmix_chisq01 <- function(x) {
  ifelse(
    x < 0,
    0,
    ifelse(x == 0, 0.5, 0.5 + 0.5 * pchisq(x, df = 1))
  )
}

qmix_chisq01 <- function(u) {
  ifelse(u <= 0.5, 0, qchisq(2 * u - 1, df = 1))
}

pval_mix_chisq01 <- function(lrt) {
  ifelse(lrt <= 0, 1, 0.5 * pchisq(lrt, df = 1, lower.tail = FALSE))
}

## =========================================================
## 2. Cache for parindekser
## =========================================================

.make_pair_cache <- local({
  cache <- new.env(parent = emptyenv())
  
  function(n) {
    key <- as.character(n)
    if (!exists(key, envir = cache, inherits = FALSE)) {
      pairs <- utils::combn(n, 2)
      assign(
        key,
        list(
          i = pairs[1, ],
          j = pairs[2, ],
          m = ncol(pairs)
        ),
        envir = cache
      )
    }
    get(key, envir = cache, inherits = FALSE)
  }
})

## =========================================================
## 3. Simulering under nullmodellen med konstant s
## =========================================================

sim_data_null_const_s <- function(n, log_r, s, d = 2) {
  stopifnot(d == 2)
  
  pc <- .make_pair_cache(n)
  i <- pc$i
  j <- pc$j
  m <- pc$m
  
  X <- matrix(rnorm(n * d), nrow = n, ncol = d)
  r <- exp(log_r)
  
  eta <- r * (X[i, 1] - X[j, 1])
  p <- plogis(eta)
  z <- rbinom(m, size = s, prob = p)
  
  list(
    z = z,
    i = i,
    j = j,
    s = rep.int(s, m),
    n = n,
    X = X,
    eta = eta,
    p = p,
    true_log_r = log_r
  )
}

## =========================================================
## 4. Model
## =========================================================
invlogit_ad <- function(eta) {
  1 / (1 + exp(-eta))
}

lid_eta <- function(x, i, j, log_r, a) {
  r <- exp(log_r)
  
  r * (x[i, 1] - x[j, 1]) +
    a * (x[i, 2] * x[j, 1] - x[i, 1] * x[j, 2])
}

bt_eta <- function(x, i, j, log_r) {
  r <- exp(log_r)
  
  r * (x[i, 1] - x[j, 1])
}

f_null <- function(parms, data) {
  z <- OBS(data$z)
  
  x <- parms$x
  eta <- bt_eta(x, data$i, data$j, parms$log_r)
  
  nll <- 0
  
  # Do not use as.vector(x) here if you want checkConsistency()
  nll <- nll - sum(dnorm(x, mean = 0, sd = 1, log = TRUE))
  
  nll <- nll - sum(dbinom(
    z,
    size = data$s,
    prob = invlogit_ad(eta),
    log = TRUE
  ))
  
  nll
}

f_full <- function(parms, data) {
  z <- OBS(data$z)
  
  x <- parms$x
  eta <- lid_eta(x, data$i, data$j, parms$log_r, parms$a)
  
  nll <- 0
  
  # Direct density on x, not as.vector(x)
  nll <- nll - sum(dnorm(x, mean = 0, sd = 1, log = TRUE))
  
  nll <- nll - sum(dbinom(
    z,
    size = data$s,
    prob = invlogit_ad(eta),
    log = TRUE
  ))
  
  nll
}

f_only_theta <- function(parms, data) {
  z <- OBS(data$z)
  
  x <- parms$x
  theta_raw <- parms$theta_raw
  theta <- exp(parms$log_sigma_theta) * theta_raw
  
  eta <- bt_eta(x, data$i, data$j, parms$log_r) + theta
  
  nll <- 0
  
  nll <- nll - sum(dnorm(x, mean = 0, sd = 1, log = TRUE))
  nll <- nll - sum(dnorm(theta_raw, mean = 0, sd = 1, log = TRUE))
  
  nll <- nll - sum(dbinom(
    z,
    size = data$s,
    prob = invlogit_ad(eta),
    log = TRUE
  ))
  
  nll
}
cmb <- function(f, d) function(p) f(p, d)

## =========================================================
## 5. Tilpasning
## =========================================================

fit_model <- function(data,
                      func,
                      parms,
                      random = "x",
                      lower = NULL,
                      upper = NULL) {
  obj <- RTMB::MakeADFun(
    cmb(func, data),
    parameters = parms,
    random = random,
    silent = TRUE
  )
  
  if (is.null(lower)) {
    lower <- rep(-Inf, length(obj$par))
  }
  
  if (is.null(upper)) {
    upper <- rep(Inf, length(obj$par))
  }
  
  stopifnot(length(lower) == length(obj$par))
  stopifnot(length(upper) == length(obj$par))
  
  opt <- try(
    nlminb(
      start = obj$par,
      objective = obj$fn,
      gradient = obj$gr,
      lower = lower,
      upper = upper,
      control = list(eval.max = 1000, iter.max = 1000)
    ),
    silent = TRUE
  )
  
  if (inherits(opt, "try-error")) {
    return(list(
      converged = FALSE,
      logLik = NA_real_,
      obj = obj,
      opt = NULL,
      sdr = NULL
    ))
  }
  
  logLik <- try(-obj$fn(opt$par), silent = TRUE)
  
  if (inherits(logLik, "try-error") || !is.finite(logLik)) {
    return(list(
      converged = FALSE,
      logLik = NA_real_,
      obj = obj,
      opt = opt,
      sdr = NULL
    ))
  }
  
  sdr <- try(sdreport(obj), silent = TRUE)
  
  list(
    converged = isTRUE(opt$convergence == 0),
    logLik = as.numeric(logLik),
    obj = obj,
    opt = opt,
    sdr = if (inherits(sdr, "try-error")) NULL else sdr
  )
}

fit_models_once <- function(data,
                            init_log_r_null = 0,
                            init_log_r_full = 0,
                            init_a_full = 0.1,
                            init_log_r_theta = 0,
                            init_sigma_theta = 0.1,
                            tol = 1e-7) {
  n <- data$n
  m <- length(data$z)
  
  parms_null <- list(
    x = matrix(0, n, 2),
    log_r = init_log_r_null
  )
  
  parms_full <- list(
    x = matrix(0, n, 2),
    log_r = init_log_r_full,
    a = init_a_full
  )
  
  parms_theta <- list(
    x = matrix(0, n, 2),
    theta_raw = rep(0, m),
    log_r = init_log_r_theta,
    log_sigma_theta = log(init_sigma_theta)
  )
  
  ## Null model: BT
  fit0 <- fit_model(
    data = data,
    func = f_null,
    parms = parms_null,
    random = "x"
  )
  
  ## Full model: latent intransitive model
  fit1 <- fit_model(
    data = data,
    func = f_full,
    parms = parms_full,
    random = "x",
    lower = c(log_r = -Inf, a = 0),
    upper = c(log_r =  Inf, a = Inf)
  )
  
  ## Dyad-augmented BT model
  fit_theta <- fit_model(
    data = data,
    func = f_only_theta,
    parms = parms_theta,
    random = c("x", "theta_raw"),
  )
  
  l0 <- fit0$logLik
  l1 <- fit1$logLik
  lt <- fit_theta$logLik
  
  lrt_full_vs_null <- if (is.finite(l0) && is.finite(l1)) {
    2 * (l1 - l0)
  } else {
    NA_real_
  }
  
  lrt_theta_vs_null <- if (is.finite(l0) && is.finite(lt)) {
    2 * (lt - l0)
  } else {
    NA_real_
  }
  
  lrt_full_vs_null <- ifelse(lrt_full_vs_null < tol, 0, lrt_full_vs_null)
  lrt_theta_vs_null <- ifelse(lrt_theta_vs_null < tol, 0, lrt_theta_vs_null)
  
  full_a_est <- NA_real_
  full_log_r_est <- NA_real_
  null_log_r_est <- NA_real_
  theta_log_r_est <- NA_real_
  theta_sigma_est <- NA_real_
  
  if (!is.null(fit1$opt)) {
    if ("a" %in% names(fit1$opt$par)) {
      full_a_est <- fit1$opt$par[["a"]]
    }
    if ("log_r" %in% names(fit1$opt$par)) {
      full_log_r_est <- fit1$opt$par[["log_r"]]
    }
  }
  
  if (!is.null(fit0$opt)) {
    if ("log_r" %in% names(fit0$opt$par)) {
      null_log_r_est <- fit0$opt$par[["log_r"]]
    }
  }
  
  if (!is.null(fit_theta$opt)) {
    if ("log_r" %in% names(fit_theta$opt$par)) {
      theta_log_r_est <- fit_theta$opt$par[["log_r"]]
    }
    if ("sigma_theta" %in% names(fit_theta$opt$par)) {
      theta_sigma_est <- fit_theta$opt$par[["sigma_theta"]]
    }
  }
  
  list(
    logLik_null = l0,
    logLik_full = l1,
    logLik_theta = lt,
    
    lrt_full_vs_null = lrt_full_vs_null,
    lrt_theta_vs_null = lrt_theta_vs_null,
    
    null_converged = fit0$converged,
    full_converged = fit1$converged,
    theta_converged = fit_theta$converged,
    
    null_opt = fit0$opt,
    full_opt = fit1$opt,
    theta_opt = fit_theta$opt,
    
    full_a_est = full_a_est,
    full_log_r_est = full_log_r_est,
    null_log_r_est = null_log_r_est,
    theta_log_r_est = theta_log_r_est,
    theta_sigma_est = theta_sigma_est
  )
}

## =========================================================
## 6. Simuler mange repetisjoner for ??n (n, s)-celle
## =========================================================

simulate_lrt_null <- function(nsim, n, log_r, s,
                              seed = NULL,
                              verbose = TRUE) {
  if (!is.null(seed)) set.seed(seed)
  
  out <- vector("list", nsim)
  
  for (b in seq_len(nsim)) {
    dat <- sim_data_null_const_s(n = n, log_r = log_r, s = s)
    fit <- fit_models_once(dat)
    
    out[[b]] <- data.frame(
      sim = b,
      n = n,
      s = s,
      true_log_r = log_r,
      
      null_converged = fit$null_converged,
      full_converged = fit$full_converged,
      theta_converged = fit$theta_converged,
      
      logLik_null = fit$logLik_null,
      logLik_full = fit$logLik_full,
      logLik_theta = fit$logLik_theta,
      
      lrt_full_vs_null = fit$lrt_full_vs_null,
      lrt_theta_vs_null = fit$lrt_theta_vs_null,
      
      pval_mix_full_vs_null = pval_mix_chisq01(fit$lrt_full_vs_null),
      pval_mix_theta_vs_null = pval_mix_chisq01(fit$lrt_theta_vs_null),
      
      full_a_est = fit$full_a_est,
      full_log_r_est = fit$full_log_r_est,
      null_log_r_est = fit$null_log_r_est,
      theta_log_r_est = fit$theta_log_r_est,
      theta_sigma_est = fit$theta_sigma_est
    )
    
    if (verbose && (b %% 50 == 0)) {
      cat(sprintf("n=%d, s=%d: done %d of %d\n", n, s, b, nsim))
    }
  }
  
  do.call(rbind, out)
}

## =========================================================
## 7. Plotfunksjon: histogram + QQ + PP
## =========================================================

plot_lrt_diagnostics <- function(lrt, n, s, label = "LRT") {
  lrt <- lrt[is.finite(lrt)]
  lrt <- sort(lrt)
  N <- length(lrt)
  
  if (N == 0) {
    plot.new()
    title(main = sprintf("n=%d, s=%d: No finite %s values", n, s, label))
    return(invisible(NULL))
  }
  
  probs <- (seq_len(N) - 0.5) / N
  theo_q <- qmix_chisq01(probs)
  theo_cdf <- pmix_chisq01(lrt)
  emp_cdf <- probs
  
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  
  par(mfrow = c(1, 3))
  
  hist(
    lrt,
    breaks = 40,
    probability = TRUE,
    border = "white",
    main = sprintf("n=%d, s=%d: %s under H0", n, s, label),
    xlab = label
  )
  
  xx <- seq(0, max(lrt), length.out = 500)
  lines(xx, 0.5 * dchisq(xx, df = 1), lwd = 2)
  abline(v = 0, lty = 2)
  mtext(sprintf("Empirical mass at 0: %.3f", mean(lrt == 0)), side = 3, line = -2)
  
  plot(
    theo_q, lrt,
    pch = 19, cex = 0.45,
    xlab = expression("Theoretical quantiles: " ~ 0.5 * chi^2[0] + 0.5 * chi^2[1]),
    ylab = paste("Empirical", label, "quantiles"),
    main = sprintf("n=%d, s=%d: QQ-plot", n, s)
  )
  abline(0, 1, col = 2, lwd = 2)
  
  plot(
    theo_cdf, emp_cdf,
    type = "l",
    lwd = 2,
    xlab = "Theoretical CDF",
    ylab = "Empirical CDF",
    main = sprintf("n=%d, s=%d: PP-plot", n, s)
  )
  abline(0, 1, col = 2, lwd = 2)
}

## =========================================================
## 8. Oppsummeringsm??l for ??n celle
## =========================================================

summarise_lrt_col <- function(res_sim,
                              lrt_col,
                              pval_col,
                              est_col = NULL,
                              test_name = lrt_col) {
  lrt <- res_sim[[lrt_col]]
  lrt_fin <- lrt[is.finite(lrt)]
  
  if (length(lrt_fin) == 0) {
    return(data.frame(
      test = test_name,
      n = unique(res_sim$n),
      s = unique(res_sim$s),
      nsim = nrow(res_sim),
      n_finite_lrt = 0,
      prop_null_converged = mean(res_sim$null_converged, na.rm = TRUE),
      prop_alt_converged = NA_real_,
      mass_at_0 = NA_real_,
      mean_lrt = NA_real_,
      median_lrt = NA_real_,
      q90_lrt = NA_real_,
      q95_lrt = NA_real_,
      prop_sig_005 = NA_real_,
      ks_mix = NA_real_,
      cvm_mix = NA_real_,
      mean_abs_est = NA_real_
    ))
  }
  
  lrt_sorted <- sort(lrt_fin)
  N <- length(lrt_sorted)
  emp_cdf <- (seq_len(N) - 0.5) / N
  theo_cdf <- pmix_chisq01(lrt_sorted)
  
  ks_mix <- max(abs(emp_cdf - theo_cdf))
  cvm_mix <- mean((emp_cdf - theo_cdf)^2)
  
  alt_conv <- if (test_name == "full_vs_null") {
    mean(res_sim$full_converged, na.rm = TRUE)
  } else if (test_name == "theta_vs_null") {
    mean(res_sim$theta_converged, na.rm = TRUE)
  } else {
    NA_real_
  }
  
  mean_abs_est <- if (!is.null(est_col)) {
    mean(abs(res_sim[[est_col]]), na.rm = TRUE)
  } else {
    NA_real_
  }
  
  data.frame(
    test = test_name,
    n = unique(res_sim$n),
    s = unique(res_sim$s),
    nsim = nrow(res_sim),
    n_finite_lrt = length(lrt_fin),
    prop_null_converged = mean(res_sim$null_converged, na.rm = TRUE),
    prop_alt_converged = alt_conv,
    mass_at_0 = mean(lrt_fin == 0),
    mean_lrt = mean(lrt_fin),
    median_lrt = median(lrt_fin),
    q90_lrt = unname(quantile(lrt_fin, 0.90)),
    q95_lrt = unname(quantile(lrt_fin, 0.95)),
    prop_sig_005 = mean(res_sim[[pval_col]] < 0.05, na.rm = TRUE),
    ks_mix = ks_mix,
    cvm_mix = cvm_mix,
    mean_abs_est = mean_abs_est
  )
}

summarise_res_sim <- function(res_sim) {
  rbind(
    summarise_lrt_col(
      res_sim = res_sim,
      lrt_col = "lrt_full_vs_null",
      pval_col = "pval_mix_full_vs_null",
      est_col = "full_a_est",
      test_name = "full_vs_null"
    ),
    summarise_lrt_col(
      res_sim = res_sim,
      lrt_col = "lrt_theta_vs_null",
      pval_col = "pval_mix_theta_vs_null",
      est_col = "theta_sigma_est",
      test_name = "theta_vs_null"
    )
  )
}
## =========================================================
## 9. Gridanalyse
## =========================================================

run_lrt_grid_analysis <- function(n_grid,
                                  s_grid,
                                  nsim,
                                  log_r,
                                  outdir = "lrt_grid_results",
                                  seed = 1,
                                  verbose = TRUE) {
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(outdir, "res_sim"), showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(outdir, "plots"), showWarnings = FALSE, recursive = TRUE)
  
  if (!is.null(seed)) set.seed(seed)
  
  grid <- expand.grid(
    n = n_grid,
    s = s_grid,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  
  summary_list <- vector("list", nrow(grid))
  
  for (k in seq_len(nrow(grid))) {
    n_val <- grid$n[k]
    s_val <- grid$s[k]
    
    if (verbose) {
      cat(sprintf("\n========== Cell %d of %d: n=%d, s=%d ==========\n",
                  k, nrow(grid), n_val, s_val))
    }
    
    ## Egen celle-seed for reproduksjon
    cell_seed <- if (is.null(seed)) NULL else seed + 1000 * k
    
    t0 <- proc.time()[3]
    
    res_sim <- simulate_lrt_null(
      nsim = nsim,
      n = n_val,
      log_r = log_r,
      s = s_val,
      seed = cell_seed,
      verbose = verbose
    )
    
    runtime_sec <- proc.time()[3] - t0
    
    ## Lagre res_sim
    res_file <- file.path(outdir, "res_sim",
                          sprintf("res_sim_n%d_s%d.rds", n_val, s_val))
    saveRDS(res_sim, res_file)
    
    ## Lagre plot
    plot_file_full <- file.path(
      outdir, "plots",
      sprintf("lrt_full_vs_null_diag_n%d_s%d.png", n_val, s_val)
    )
    
    png(filename = plot_file_full, width = 1800, height = 700, res = 150)
    plot_lrt_diagnostics(
      res_sim$lrt_full_vs_null,
      n = n_val,
      s = s_val,
      label = "LRT: LID vs BT"
    )
    dev.off()
    
    plot_file_theta <- file.path(
      outdir, "plots",
      sprintf("lrt_theta_vs_null_diag_n%d_s%d.png", n_val, s_val)
    )
    
    png(filename = plot_file_theta, width = 1800, height = 700, res = 150)
    plot_lrt_diagnostics(
      res_sim$lrt_theta_vs_null,
      n = n_val,
      s = s_val,
      label = "LRT: BT + dyad RE vs BT"
    )
    dev.off()
    
    ## Oppsummering
    sum_row <- summarise_res_sim(res_sim)
    sum_row$runtime_sec <- runtime_sec
    sum_row$res_file <- res_file
    sum_row$plot_file_full <- plot_file_full
    sum_row$plot_file_theta <- plot_file_theta
    
    summary_list[[k]] <- sum_row
    
    if (verbose) {
      cat(sprintf("Saved: %s\n", res_file))
      cat(sprintf("Saved: %s\n", plot_file_full))
      cat(sprintf("Saved: %s\n", plot_file_theta))
      cat(sprintf("Runtime: %.2f sec\n", runtime_sec))
    }
  }
  
  summary_df <- do.call(rbind, summary_list)
  
  summary_file_rds <- file.path(outdir, "grid_summary.rds")
  summary_file_csv <- file.path(outdir, "grid_summary.csv")
  
  saveRDS(summary_df, summary_file_rds)
  write.csv(summary_df, summary_file_csv, row.names = FALSE)
  
  if (verbose) {
    cat("\n========== DONE ==========\n")
    cat("Summary saved to:\n")
    cat(summary_file_rds, "\n")
    cat(summary_file_csv, "\n")
  }
  
  invisible(summary_df)
}


summary_grid <- run_lrt_grid_analysis(
  n_grid = c(150), #approx quantiles c(0.5,0.95,0.99) and max
  s_grid = c(10), #approx quantiles c(0.25,0.5,0.75,0.99)
  nsim = 1000,
  log_r = log(2),
  outdir = "lrt_null_dist_analysis/lrt_grid_analysis_thesis",
  seed = 123,
  verbose = TRUE
)

res_n150 <- readRDS("~/lrt_null_dist_analysis/lrt_grid_analysis_thesis/res_sim/res_sim_n150_s10.rds")
res_n50 <- readRDS("~/lrt_null_dist_analysis/lrt_grid_analysis_thesis/res_sim/res_sim_n50_s10.rds")

library(dplyr)
library(ggplot2)
library(rlang)

# ------------------------------------------------------------
# 0.5 chi_0^2 + 0.5 chi_1^2 mixture distribution
# ------------------------------------------------------------

pmix_chisq <- function(x, w0 = 0.5, df = 1) {
  ifelse(
    x < 0,
    0,
    w0 + (1 - w0) * pchisq(x, df = df)
  )
}

qmix_chisq <- function(p, w0 = 0.5, df = 1) {
  ifelse(
    p <= w0,
    0,
    qchisq((p - w0) / (1 - w0), df = df)
  )
}

dmix_chisq_cont <- function(x, w0 = 0.5, df = 1) {
  ifelse(
    x <= 0,
    NA_real_,
    (1 - w0) * dchisq(x, df = df)
  )
}

plot_lrt_hist_mix <- function(data,
                              lrt_col,
                              w0 = 0.5,
                              df = 1,
                              bins = 40,
                              title = NULL) {
  
  lrt_col <- enquo(lrt_col)
  
  df_plot <- data %>%
    transmute(lrt = !!lrt_col) %>%
    filter(is.finite(lrt), lrt >= 0)
  
  x_max <- max(df_plot$lrt, na.rm = TRUE)
  
  curve_df <- tibble(
    x = seq(1e-8, x_max, length.out = 1000),
    density = dmix_chisq_cont(x, w0 = w0, df = df)
  )
  
  p_zero_emp <- mean(df_plot$lrt < 1e-8)
  
  ggplot(df_plot, aes(x = lrt)) +
    geom_histogram(
      aes(y = after_stat(density)),
      bins = bins,
      boundary = 0,
      fill = "grey80",
      color = "white"
    ) +
    geom_line(
      data = curve_df,
      aes(x = x, y = density),
      linewidth = 1
    ) +
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      linewidth = 0.8
    ) +
    annotate(
      "text",
      x = 0,
      y = Inf,
      label = paste0("Empirical mass near 0: ", round(p_zero_emp, 3),
                     "\nAsymptotic mass at 0: ", w0),
      hjust = -0.05,
      vjust = 1.2,
      size = 3.5
    ) +
    labs(
      x = "LRT statistic",
      y = "Density",
      title = title %||% expression("LRT histogram vs. " * 0.5 * chi[0]^2 + 0.5 * chi[1]^2)
    ) +
    theme_bw(base_size = 20)
}
res_n50 <- readRDS("~/lrt_null_dist_analysis/lrt_grid_analysis_thesis/res_sim/res_sim_n50_s10.rds")
res_n50
res_n25 <- readRDS("~/lrt_null_dist_analysis/lrt_grid_analysis_thesis/res_sim/res_sim_n25_s10.rds")
plot_lrt_hist_mix(
  res_n150,
  lrt_full_vs_null,
  title = "Full vs null: LRT histogram"
)

plot_lrt_hist_mix(
  res_n150,
  lrt_theta_vs_null,
  title = "Theta vs null: LRT histogram"
)
hist(res_n150$lrt_full_vs_null)

plot_lrt_qq_mix <- function(data,
                            lrt_col,
                            w0 = 0.5,
                            df = 1,
                            title = NULL) {
  
  lrt_col <- enquo(lrt_col)
  
  lrt <- data %>%
    transmute(lrt = !!lrt_col) %>%
    filter(is.finite(lrt), lrt >= 0) %>%
    pull(lrt) %>%
    sort()
  
  n <- length(lrt)
  p <- ppoints(n)
  
  qq_df <- tibble(
    theoretical = qmix_chisq(p, w0 = w0, df = df),
    empirical = lrt
  )
  
  max_val <- max(c(qq_df$theoretical, qq_df$empirical), na.rm = TRUE)
  
  ggplot(qq_df, aes(x = theoretical, y = empirical)) +
    geom_point(alpha = 0.6, size = 3.6) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    coord_equal(xlim = c(0, max_val), ylim = c(0, max_val)) +
    labs(
      x = "Theoretical quantiles",
      y = "Empirical quantiles",
      title = title %||% expression("QQ plot against " * 0.5 * chi[0]^2 + 0.5 * chi[1]^2)
    ) +
    theme_bw(base_size = 24)
}
results <- readRDS("~/DomArchiveResults/four_models_comparison_filtered_betterdiag_results.rds")
Shimoji_2014c <- results %>%
  filter(fileid == "Shimoji_2014c")
n <- length(Shimoji_2014c$lrt)
qq_df <- tibble(
  theoretical = qmix_chisq(p, w0 = w0, df = df),
  empirical = lrt
)

max_val <- max(c(qq_df$theoretical, qq_df$empirical), na.rm = TRUE)

ggplot(qq_df, aes(x = theoretical, y = empirical)) +
  geom_point(alpha = 0.6, size = 3.6) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  coord_equal(xlim = c(0, max_val), ylim = c(0, max_val)) +
  labs(
    x = "Theoretical quantiles",
    y = "Empirical quantiles",
    title = title %||% expression("QQ plot against " * 0.5 * chi[0]^2 + 0.5 * chi[1]^2)
  ) +
  theme_bw(base_size = 24)


plot_lrt_qq_mix(
  res_n25,
  lrt_full_vs_null,
  title = ""
)

plot_lrt_qq_mix(
  res_n25,
  lrt_theta_vs_null,
  title = ""
)
plot_lrt_pp_mix <- function(data,
                            lrt_col,
                            w0 = 0.5,
                            df = 1,
                            title = NULL) {
  
  lrt_col <- enquo(lrt_col)
  
  lrt <- data %>%
    transmute(lrt = !!lrt_col) %>%
    filter(is.finite(lrt), lrt >= 0) %>%
    pull(lrt) %>%
    sort()
  
  n <- length(lrt)
  
  pp_df <- tibble(
    empirical_p = ppoints(n),
    theoretical_p = pmix_chisq(lrt, w0 = w0, df = df),
    lrt = lrt
  )
  
  ggplot(pp_df, aes(x = theoretical_p, y = empirical_p)) +
    geom_point(alpha = 0.6, size = 3.6) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    coord_equal(xlim = c(0.5, 1), ylim = c(0, 1)) +
    labs(
      x = "Asymptotic mixture CDF",
      y = "Empirical CDF",
      title = title %||% expression("PP plot against " * 0.5 * chi[0]^2 + 0.5 * chi[1]^2)
    ) +
    theme_bw(base_size = 24) +
    theme(aspect.ratio = 1)
}
plot_lrt_pp_mix(
  res_n25,
  lrt_full_vs_null,
  title = ""
)
quantile(results$interactions_per_observed_dyad, c(0.25,0.75))
results %>%
  select(fileid, proportion_unknown) %>%
  arrange(desc(proportion_unknown))
plot_lrt_pp_mix(
  res_n25,
  lrt_theta_vs_null,
  title = ""
)
