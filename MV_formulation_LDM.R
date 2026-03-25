library(RTMB)
#install.packages("BradleyTerry2")
library(BradleyTerry2)

set.seed(42)
d <- 2
n <- 50
beta <- c(1,1)

sim_data <- function(n, beta, pois = TRUE, lambda = 0.5, d = 2, a = 2, s = 1){
  
  stopifnot(length(beta) == d)
  if (d != 2) stop("Som skrevet: A er 2x2, s?? d m?? v??re 2.")
  
  X <- matrix(rnorm(n * d), nrow = n)
  beta <- as.numeric(beta)
  
  A <- matrix(c(0, -a,
                a,  0), 2, 2)
  
  # Full eta
  bX <- drop(X %*% beta)
  eta <- outer(bX, bX, "-") + X %*% A %*% t(X)
  diag(eta) <- 0
  prob <- plogis(eta)
  
  # Alle kombinasjoner i<j
  pairs <- t(combn(n, 2))
  
  i_vec <- pairs[,1]
  j_vec <- pairs[,2]
  
  if(pois){
    size <- rpois(nrow(pairs), lambda)
  } else {
    size <- rep(s, nrow(pairs))
  }
  # Trekker kampresultater
  z <- rbinom(nrow(pairs), size = size,
              prob = prob[cbind(i_vec, j_vec)])
  
  list(
    z = z,
    i = i_vec,
    j = j_vec,
    X = X,
    eta = eta,
    prob = prob[cbind(i_vec, j_vec)],
    s = size,
    n = n
  )
}

parameters_plain <- list(beta = rep(0.1,d), x    = matrix(0,n,d), a    = 0.1)


f_plain <- function(parms, data) {
  getAll(data, parms, warn = FALSE)
  z <- OBS(z)
  # Forutsetter at data inneholder i, j, s
  # x: N x 2 (latent/random)
  # beta: length 2 (fixed)
  # a: scalar (fixed)
  
  # Bygg A (2x2 antisymmetrisk) fra a
  A <- matrix(0, 2, 2)
  A[1,2] <- -a
  A[2,1] <-  a
  
  # NLL
  nll <- 0
  
  # Prior p?? x (isotrop; juster ved behov)
  nll <- nll - sum(dnorm(x, 0, 1, log = TRUE))
  
  # AD-sikker init av eta (lengde = length(z))
  eta <- (x[,1] * 0)   # lengde N, riktig AD-type
  eta <- eta[i] * 0    # lengde length(i)=length(z)
  
  for (k in seq_along(z)) {
    ii <- i[k]
    jj <- j[k]
    
    # line??r del: beta^T (x_i - x_j)
    dx1 <- x[ii, 1] - x[jj, 1]
    dx2 <- x[ii, 2] - x[jj, 2]
    lin <- beta[1] * dx1 + beta[2] * dx2
    
    # biline??r del: x_i^T A x_j
    # (1x2) %*% (2x2) %*% (2x1) -> 1x1
    xi <- matrix(x[ii, ], nrow = 1)
    xj <- matrix(x[jj, ], nrow = 2)     # kolonne
    bil <- (xi %*% A %*% xj)[1,1]
    
    eta[k] <- lin + bil
  }
  
  p <- plogis(eta)
  
  # Binomial likelihood (s kan v??re skalar eller vektor)
  nll <- nll - sum(dbinom(z, size = s, prob = p, log = TRUE))
  
  nll
}
f_null <- function(parms, data) {
  z <- OBS(data$z)
  i <- data$i
  j <- data$j
  s <- data$s
  
  x <- parms$x
  log_r <- parms$log_r
  
  r <- exp(log_r)
  
  nll <- 0
  nll <- nll - sum(dnorm(x, 0, 1, log = TRUE))
  
  eta <- r * (x[i, 1] - x[j, 1])
  p <- plogis(eta)
  
  nll <- nll - sum(dbinom(z, size = s, prob = p, log = TRUE))
  
  nll
}
f_full <- function(parms, data) {
  z <- OBS(data$z)
  i <- data$i
  j <- data$j
  s <- data$s
  
  x <- parms$x
  log_r <- parms$log_r
  a <- parms$a
  
  r <- exp(log_r)
  beta_rot <- c(r, 0)
  A <- matrix(c(0, -a, a, 0), 2, 2, byrow = TRUE)
  
  nll <- 0
  nll <- nll - sum(dnorm(x, 0, 1, log = TRUE))
  
  Xi <- x[i, , drop = FALSE]
  Xj <- x[j, , drop = FALSE]
  
  lin <- drop((Xi - Xj) %*% beta_rot)
  bil <- rowSums((Xi %*% A) * Xj)
  
  eta <- lin + bil
  p <- plogis(eta)
  
  nll <- nll - sum(dbinom(z, size = s, prob = p, log = TRUE))
  nll
}

cmb <- function(f, d) function(p) f(p, d)
fit_model <- function(data, func, parms) {
  obj <- MakeADFun(cmb(func, data),
                   parms,
                   random = "x")
  opt <- nlminb(obj$par, obj$fn, obj$gr)
  logLik <- -obj$fn(opt$par)
  
  if (!is.finite(logLik)) {
    stop("Non-finite logLik")
  }
  sdr <- sdreport(obj)
  
  list(
    obj = obj,
    opt = opt,
    sdr = sdr,
    logLik = logLik
  )
}
get_estimates <- function(fit) {
  #rep <- fit$obj$env$report(fit$opt$par, fit$obj$env$last.par.best)
  par <- fit$obj$env$parList(fit$opt$par)
  par
}

lrt <- function(loglik0, loglik1) {
  d <- -2 * (loglik0 - loglik1)
  d
}
lrt_sim <- function(n, nsim = 40, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  vals <- rep(NA_real_, nsim)
  beta <- c(1, 1)
  
  for (b in seq_len(nsim)) {
    cat("simulation:", b, "\n")
    
    data <- sim_data(n = n, beta = beta, lambda = 1.0, a = 0)
    
    fit_null <- tryCatch(
      fit_model(data, f_null, params0),
      warning = function(w) w,
      error = function(e) e
    )
    
    if (inherits(fit_null, "warning") || inherits(fit_null, "error")) {
      cat("null model issue at simulation", b, "\n")
      print(fit_null)
      next
    }
    
    fit_full <- tryCatch(
      fit_model(data, f_full, params),
      warning = function(w) w,
      error = function(e) e
    )
    
    if (inherits(fit_full, "warning") || inherits(fit_full, "error")) {
      cat("full model issue at simulation", b, "\n")
      print(fit_full)
      next
    }
    
    if (is.finite(fit_null$logLik) && is.finite(fit_full$logLik)) {
      if (fit_full$logLik < fit_null$logLik) {
        cat("negative LRT candidate at simulation", b, "\n")
      } else {
        vals[b] <- lrt(fit_null$logLik, fit_full$logLik)
      }
    } else {
      cat("non-finite logLik at simulation", b, "\n")
    }
  }
  
  vals
}
res0 <- readRDS("lrt_null_grid_seed123/lrt_grid_all_results.rds")

safe_fit_model <- function(data, func, parms) {
  res <- tryCatch({
    fit <- fit_model(data, func, parms)
    
    if (!is.finite(fit$logLik)) {
      return(list(
        ok = FALSE,
        fit = NULL,
        logLik = NA_real_,
        message = "Non-finite logLik"
      ))
    }
    
    list(
      ok = TRUE,
      fit = fit,
      logLik = fit$logLik,
      message = NA_character_
    )
  }, warning = function(w) {
    list(
      ok = FALSE,
      fit = NULL,
      logLik = NA_real_,
      message = paste("Warning:", conditionMessage(w))
    )
  }, error = function(e) {
    list(
      ok = FALSE,
      fit = NULL,
      logLik = NA_real_,
      message = paste("Error:", conditionMessage(e))
    )
  })
  res
}
run_lrt_scenario <- function(
    n,
    lambda,
    nsim,
    beta = c(1, 1),
    a_true = 0,
    save_path = NULL
) {
  out <- vector("list", nsim)
  
  params_full <- list(
    x = matrix(0, n, 2),
    log_r = 0,
    a = 0.1
  )
  
  params_null <- list(
    x = matrix(0, n, 2),
    log_r = 0
  )
  
  for (b in seq_len(nsim)) {
    cat("n =", n, "lambda =", lambda, "sim =", b, "\n")
    
    row <- list(
      n = n,
      lambda = lambda,
      sim = b,
      logLik_null = NA_real_,
      logLik_full = NA_real_,
      lrt = NA_real_,
      null_ok = FALSE,
      full_ok = FALSE,
      null_msg = NA_character_,
      full_msg = NA_character_
    )
    
    data <- tryCatch({
      sim_data(n = n, beta = beta, lambda = lambda, a = a_true)
    }, error = function(e) e)
    
    if (inherits(data, "error")) {
      row$null_msg <- paste("sim_data failed:", conditionMessage(data))
      row$full_msg <- row$null_msg
      out[[b]] <- as.data.frame(row, stringsAsFactors = FALSE)
      next
    }
    
    fit_null <- safe_fit_model(data, f_null, params_null)
    fit_full <- safe_fit_model(data, f_full, params_full)
    
    row$logLik_null <- fit_null$logLik
    row$logLik_full <- fit_full$logLik
    row$null_ok <- fit_null$ok
    row$full_ok <- fit_full$ok
    row$null_msg <- fit_null$message
    row$full_msg <- fit_full$message
    
    if (isTRUE(fit_null$ok) && isTRUE(fit_full$ok)) {
      row$lrt <- lrt(fit_null$logLik, fit_full$logLik)
    }
    
    out[[b]] <- as.data.frame(row, stringsAsFactors = FALSE)
    
    if (!is.null(save_path) && b %% 25 == 0) {
      tmp <- do.call(rbind, out[seq_len(b)])
      saveRDS(tmp, save_path)
    }
  }
  
  res <- do.call(rbind, out)
  
  if (!is.null(save_path)) {
    saveRDS(res, save_path)
  }
  
  res
}
run_lrt_grid <- function(
    n_values,
    lambda_values,
    nsim,
    beta = c(1, 1),
    a_true = 0,
    out_dir = "lrt_grid_results"
) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  all_results <- list()
  idx <- 1
  
  for (n in n_values) {
    for (lambda in lambda_values) {
      file_name <- paste0(
        "lrt_null_n_", n,
        "_lambda_", gsub("\\.", "_", as.character(lambda)),
        ".rds"
      )
      save_path <- file.path(out_dir, file_name)
      
      res <- run_lrt_scenario(
        n = n,
        lambda = lambda,
        nsim = nsim,
        beta = beta,
        a_true = a_true,
        save_path = save_path
      )
      
      all_results[[idx]] <- res
      idx <- idx + 1
    }
  }
  
  combined <- do.call(rbind, all_results)
  saveRDS(combined, file.path(out_dir, "lrt_grid_all_results.rds"))
  combined
}

set.seed(12345)
res <- run_lrt_grid(
  n_values = c(8, 10, 12, 15, 20),
  lambda_values = c(0.5, 1., 1.5, 2., 3.),
  nsim = 500,
  beta = c(1, 1),
  a_true = 0,
  out_dir = "lrt_null_grid"
)

res_ok <- subset(res, null_ok & full_ok & is.finite(lrt) & lrt >= 0)
library(ggplot2)
res
ggplot(res_ok, aes(x = lrt)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 20,
                 color = "white") +
  stat_function(
    fun = dchisq,
    args = list(df = 1),
    linewidth = 0.8
  ) +
  facet_grid(lambda ~ n, scales = "fixed") +
  labs(
    title = "Empirical null distribution of LRT by n and lambda",
    x = "LRT statistic",
    y = "Density"
  ) +
  theme_bw()

plot_nonfinite_full_loglik <- function(res, digits = 2) {
  if (!all(c("n", "lambda", "logLik_full") %in% names(res))) {
    stop("res must contain columns: n, lambda, logLik_full")
  }
  
  # indikator for non-finite logLik i full modell
  tmp <- res
  tmp$nonfinite_full <- !is.finite(tmp$logLik_full)
  
  # oppsummering per (n, lambda)
  summary_df <- aggregate(
    nonfinite_full ~ n + lambda,
    data = tmp,
    FUN = mean
  )
  
  summary_df$percent_nonfinite <- summary_df$nonfinite_full
  
  # lag matrise til base-R plot
  n_vals <- sort(unique(summary_df$n))
  lambda_vals <- sort(unique(summary_df$lambda))
  
  mat <- matrix(
    NA_real_,
    nrow = length(lambda_vals),
    ncol = length(n_vals),
    dimnames = list(
      paste0("lambda=", lambda_vals),
      paste0("n=", n_vals)
    )
  )
  
  for (k in seq_len(nrow(summary_df))) {
    i <- match(summary_df$lambda[k], lambda_vals)
    j <- match(summary_df$n[k], n_vals)
    mat[i, j] <- summary_df$percent_nonfinite[k]
  }
  
  # plot
  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  
  par(mar = c(5, 5, 4, 6))
  
  image(
    x = seq_along(n_vals),
    y = seq_along(lambda_vals),
    z = t(mat),
    axes = FALSE,
    xlab = "n",
    ylab = expression(lambda),
    main = "Percent non-finite logLik (full model)",
    col = hcl.colors(20, "YlOrRd", rev = FALSE)
  )
  
  axis(1, at = seq_along(n_vals), labels = n_vals)
  axis(2, at = seq_along(lambda_vals), labels = lambda_vals)
  box()
  
  # skriv prosent i hver celle
  for (i in seq_along(n_vals)) {
    for (j in seq_along(lambda_vals)) {
      if (!is.na(mat[j, i])) {
        text(i, j, labels = round(mat[j, i], digits))
      }
    }
  }
  
  # enkel fargeskala
  usr <- par("usr")
  x_right <- usr[2] + 0.6
  y_seq <- seq(usr[3], usr[4], length.out = 21)
  z_seq <- seq(min(mat, na.rm = TRUE), max(mat, na.rm = TRUE), length.out = 20)
  cols <- hcl.colors(20, "YlOrRd", rev = TRUE)
  
  for (k in seq_along(cols)) {
    rect(x_right, y_seq[k], x_right + 0.25, y_seq[k + 1], col = cols[k], border = NA, xpd = TRUE)
  }
  text(x_right + 0.35, usr[4], labels = round(max(mat, na.rm = TRUE), digits), adj = 0, xpd = TRUE)
  text(x_right + 0.35, usr[3], labels = round(min(mat, na.rm = TRUE), digits), adj = 0, xpd = TRUE)
  text(x_right, usr[4] + 0.05 * diff(usr[3:4]), labels = "%", adj = 0, xpd = TRUE)
  
  invisible(list(
    summary = summary_df[order(summary_df$lambda, summary_df$n), ],
    matrix = mat
  ))
}
out <- plot_nonfinite_full_loglik(res0)
out$matrix
params <- list(log_r = 0.1, x = matrix(0,n,d), a = 0.1)
params0 <- list(log_r = 0.1, x = matrix(0,n,d))
lrt_vals <- lrt_sim(n = n, nsim = 40, seed = 123)
hist(lrt_vals,
     probability = TRUE,
     breaks = 20,
     main = "Empirical null distribution of LRT",
     xlab = "LRT statistic",
     border = "white")

curve(dchisq(x, df = 1),
      add = TRUE,
      lwd = 2,
      lty = 1)

Tdata <- sim_data(n=10,beta,lambda=1,a=0)
fit_full <- fit_model(data, f_full, params)
fit_null <- fit_model(data, f_null, params0)





predict_prob_bt_intrans <- function(data, est, output = "df") {
  prob <- with(c(data, est), {
    A <- matrix(c(0, -a,
                  a,  0), 2, 2)
    beta <- as.numeric(c(exp(log_r), 0))
    bx <- drop(x %*% beta)
    eta <- outer(bx, bx, "-") + x %*% A %*% t(x)
    diag(eta) <- 0
    prob <- plogis(eta)
    })
    if(output == "matrix"){
      return(prob)
    }
    
    p_obs <- prob[cbind(data$i, data$j)]
    
    if(output == "vector"){
      return(p_obs)
    }
    
    data.frame(
      i = data$i,
      j = data$j,
      p = p_obs
    )
}

make_count_matrix <- function(z, s, i, j, n) {
  M <- matrix(0, n, n)
  for(k in seq_along(z)){
    M[i[k],j[k]] <- M[i[k],j[k]] + z[k]
    M[j[k],i[k]] <- M[j[k],i[k]] + (s[k] - z[k])
  }
  diag(M) <- 0
  M
}
binary_dominance_matrix <- function(M) {
  B <- matrix(0, nrow(M), ncol(M))
  B[M > t(M)] <- 1
  diag(B) <- 0
  B
}
complete_binary_matrix_random <- function(M) {
  n <- nrow(M)
  stopifnot(n == ncol(M))
  
  B <- matrix(0L, n, n)
  upper <- upper.tri(M)
  
  m1 <- M[upper]
  m2 <- t(M)[upper]
  
  vals <- integer(length(m1))
  vals[m1 > m2] <- 1L
  vals[m1 < m2] <- 0L
  
  ties <- (m1 == m2)
  vals[ties] <- rbinom(sum(ties), size = 1, prob = 0.5)
  
  B[upper] <- vals
  B[lower.tri(B)] <- 1L - t(B)[lower.tri(B)]
  
  diag(B) <- 0L
  B
}

proportion_unknown <- function(M) {
  N <- M + t(M)
  n <- ncol(N)
  idx <- N == 0
  unknown = (sum(idx) - n) / 2
  prop_unknwn <- unknown / sum(upper.tri(N))
  prop_unknwn
}
dci <- function(M) {
  n <- nrow(M)
  vals <- c()
  for (a in 1:(n-1)) {
    for (b in (a+1):n) {
      H <- max(M[a,b], M[b,a])
      L <- min(M[a,b], M[b,a])
      if (H + L > 0) {
        vals <- c(vals, (H - L) / (H + L))
      }
    }
  }
  mean(vals)
}
ttri <- function(M) {
  B <- binary_dominance_matrix(M)
  n <- nrow(B)
  transitive <- 0
  total_defined <- 0
  
  for (a in 1:(n-2)) {
    for (b in (a+1):(n-1)) {
      for (c in (b+1):n) {
        pairs_defined <- (
          (B[a,b] + B[b,a] == 1) &&
            (B[a,c] + B[c,a] == 1) &&
            (B[b,c] + B[c,b] == 1)
        )
        if (pairs_defined) {
          total_defined <- total_defined + 1
          outdeg <- c(
            B[a,b] + B[a,c],
            B[b,a] + B[b,c],
            B[c,a] + B[c,b]
          )
          if (all(sort(outdeg) == c(0,1,2))) {
            transitive <- transitive + 1
          }
        }
      }
    }
  }
  
  if (total_defined == 0) return(NA_real_)
  pt <- transitive / total_defined
  4 * (pt - 0.75)
}
davids_scores <- function(M, normalize = FALSE) {
  n <- nrow(M)
  
  # Antall m??ter per par
  N <- M +t(M)
  
  # Proporsjon seire
  P <- matrix(0, n, n)
  idx <- N > 0
  P[idx] <- M[idx] / N[idx]
  diag(P) <- 0
  
  w <- rowSums(P)
  l <- colSums(P)
  
  w2 <- as.vector(P %*% w)
  l2 <- as.vector(t(P) %*% l)
  
  ds <- w + w2 - l - l2
  
  if (!normalize){
    return(ds)
  }
  # Enkel normalisering til [0, n-1] ved rangbevaring
  # Dette er ikke eneste mulige normalisering, men er praktisk og stabil.
  r <- rank(ds, ties.method = "average")
  ds_norm <- (r - 1)
  return(ds_norm)
}
steepness <- function(M) {
  ds_norm <- davids_scores(M, normalize = TRUE)
  ds_sorted <- sort(ds_norm, decreasing = TRUE)
  ranks <- seq_along(ds_sorted)
  fit <- lm(ds_sorted ~ ranks)
  abs(coef(fit)[2])
}
landau_h <- function(B) {
  #B <- binary_dominance_matrix(M)
  n <- nrow(B)
  R <- rowSums(B)
  
  h <- (12 / (n^3 - n)) * sum((R - (n - 1) / 2)^2)
  h
}
modified_landau_h <- function(M, nrep = 1000) {
  vals <- numeric(nrep)
  for(i in seq_len(nrep)){
    B <- complete_binary_matrix_random(M)
    vals[i] <- landau_h(B)
  }
  mean(vals)
}

stat_prop_unknown <- function(z, s, i, j, n) {
  M <- make_count_matrix(z, s, i, j, n)
  proportion_unknown(M)
}
stat_dci <- function(z, s, i, j, n) {
  M <- make_count_matrix(z, s, i, j, n)
  dci(M)
}
stat_ttri <- function(z, s, i, j, n) {
  M <- make_count_matrix(z, s, i, j, n)
  ttri(M)
}
stat_davids_score <- function(z, s, i, j, n, normalize = FALSE) {
  M <- make_count_matrix(z, s, i, j, n)
  davids_scores(M, normalize = normalize)
}
stat_steepness <- function(z, s, i, j, n) {
  M <- make_count_matrix(z, s, i, j, n)
  steepness(M)
}
stat_landau_h <- function(z, s, i, j, n) {
  M <- make_count_matrix(z, s, i, j, n)
  B <- binary_dominance_matrix(M)
  landau_h(B)
}
stat_modified_landau_h <- function(z, s, i, j, n, nrep = 100) {
  M <- make_count_matrix(z, s, i, j, n)
  modified_landau_h(M, nrep = nrep)
}

posterior_predictive_check <- function(fit, data, nsim = 1000, stat_fun) {
  est <- get_estimates(fit)
  p <- predict_prob_bt_intrans(data, est, output = "vector")
  
  z_obs <- data$z
  s <- data$s
  stat_obs <- stat_fun(z_obs, s, data$i, data$j, data$n)
  
  stat_rep <- numeric(nsim)
  
  for(l in seq_len(nsim)){
    z_rep <- rbinom(length(p), size = s, prob = p)
    stat_rep[l] <- stat_fun(z_rep, s, data$i, data$j, data$n)
  }
  
  ppp <- mean(stat_rep >= stat_obs)
  
  list(
    observed = stat_obs,
    replicated = stat_rep,
    ppp = ppp,
    p = p
  )
}
data <- sim_data(n=10,beta,pois=TRUE, lambda=2)
m <- make_count_matrix(data$z,data$s,data$i,data$j,data$n)
modified_landau_h(m)
fit <- fit_bt_intrans(data, f_constrained, params)
params <- list(
  log_r = 0.1,
  x = matrix(0,n,d),
  a = 0.1
)
ppc <- posterior_predictive_check(
  fit = fit,
  data = data,
  nsim = 1000,
  stat_fun = stat_modified_landau_h
)
ppc$observed
ppc$ppp
ppc$replicated
hist(ppc$replicated)
abline(v = ppc$observed, lwd = 2)
##################################################################
# Data set from DomArchive
data
install.packages("devtools")
devtools::install_github("DomArchive/DomArchive", build_vignettes = FALSE)
dom.data$Adcock_2015a$matrix
dom.data$Adcock_2015a$metadata
str(dom.data[[1]]$metadata)
dom.data$
library(purrr)
library(dplyr)

metadata_df <- map_dfr(
  dom.data,
  ~ as.data.frame(.x$metadata),
  .id = "id"
)
metadata_df
###################################################################