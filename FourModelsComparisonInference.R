library(DomArchive)
library(dplyr)
library(purrr)
library(tidyr)
library(ggplot2)
library(plotly)
library(RTMB)
library(tidyverse)
library(patchwork)
library(scales)
library(gt)
library(knitr)

# ============================================================
# 1. Data preparation
# ============================================================
set.seed(123)
make_pair_data <- function(M, drop_zero = TRUE) {
  stopifnot(is.matrix(M), nrow(M) == ncol(M))
  diag(M) <- 0
  
  if (anyNA(M)) stop("Matrix contains NA values outside the diagonal.")
  
  n <- nrow(M)
  pairs <- t(combn(n, 2))
  i <- pairs[, 1]
  j <- pairs[, 2]
  
  z <- as.numeric(M[cbind(i, j)])
  s <- as.numeric(M[cbind(i, j)] + M[cbind(j, i)])
  
  if (drop_zero) {
    keep <- s > 0
    i <- i[keep]
    j <- j[keep]
    z <- z[keep]
    s <- s[keep]
  }
  
  list(
    z = z,
    i = i,
    j = j,
    s = s,
    n = n,
    n_dyads_observed = length(s),
    n_interactions = sum(s)
  )
}

cmb <- function(f, d) {
  function(p) f(p, d)
}

# ============================================================
# 2. Model components
# ============================================================

nll_x_prior <- function(x) {
  -sum(RTMB::dnorm(as.vector(x), mean = 0, sd = 1, log = TRUE))
}

nll_theta_prior <- function(theta_raw) {
  -sum(RTMB::dnorm(theta_raw, mean = 0, sd = 1, log = TRUE))
}
invlogit_ad <- function(eta) {
  1 / (1 + exp(-eta))
}

binom_nll <- function(z, s, eta) {
  p <- invlogit_ad(eta)
  -sum(RTMB::dbinom(z, size = s, prob = p, log = TRUE))
}
# ============================================================
# 2. Linear predictors
# ============================================================

lid_eta <- function(x, i, j, log_r, a) {
  r <- exp(log_r)
  
  r * (x[i, 1] - x[j, 1]) +
    a * (x[i, 2] * x[j, 1] - x[i, 1] * x[j, 2])
}

bt_eta <- function(x, i, j, log_r) {
  r <- exp(log_r)
  
  r * (x[i, 1] - x[j, 1])
}


# ============================================================
# 3. Models
# ============================================================

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

f_full_theta <- function(parms, data) {
  z <- OBS(data$z)
  x <- parms$x
  theta_raw <- parms$theta_raw
  theta <- exp(parms$log_sigma_theta) * theta_raw
  
  eta <- lid_eta(x, data$i, data$j, parms$log_r, parms$a) + theta
  
  nll <- 0
  nll <- nll - sum(dnorm(x, mean = 0, sd = 1, log = TRUE))
  nll <- nll - sum(dnorm(theta_raw, mean = 0, sd = 1, log = TRUE))
  nll <- nll - sum(dbinom(
    z,
    size = data$s,
    prob = plogis(eta),
    log = TRUE
  ))
  nll
}
# ============================================================
# 4. Parameter constructors
# ============================================================

make_pars <- function(data, model,
                      a_start = 0.1,
                      log_r_start = 0,
                      log_sigma_theta_start = -2) {
  
  n <- data$n
  n_dyads <- length(data$z)
  
  switch(
    model,
    
    null = list(
      x = matrix(0, n, 1),
      log_r = log_r_start
    ),
    
    full = list(
      x = matrix(0, n, 2),
      log_r = log_r_start,
      a = a_start
    ),
    
    only_theta = list(
      x = matrix(0, n, 1),
      log_r = log_r_start,
      log_sigma_theta = log_sigma_theta_start,
      theta_raw = rep(0, n_dyads)
    ),
    
    full_theta = list(
      x = matrix(0, n, 2),
      log_r = log_r_start,
      a = a_start,
      log_sigma_theta = log_sigma_theta_start,
      theta_raw = rep(0, n_dyads)
    ),
    
    stop("Unknown model: ", model)
  )
}

# ============================================================
# 5. Model specifications
# ============================================================

model_specs <- tibble(
  model = c("null", "full", "only_theta", "full_theta"),
  func = list(f_null, f_full, f_only_theta, f_full_theta),
  random = list("x", "x", c("x", "theta_raw"), c("x", "theta_raw")),
  constrain_a = c(FALSE, TRUE, FALSE, TRUE),
  k_aic = c(1, 2, 2, 3),
  start_grid = list(
    tibble(a_start = NA_real_, log_sigma_theta_start = NA_real_),
    crossing(a_start = c(0.01, 0.1, 1),
             log_sigma_theta_start = NA_real_),
    crossing(a_start = NA_real_,
             log_sigma_theta_start = c(-2, -1, 0)),
    crossing(a_start = c(0.01, 0.1, 1),
             log_sigma_theta_start = c(-2, -1, 0))
  )
)

# ============================================================
# 6. Fitting helpers
# ============================================================

make_bounds <- function(par,
                        constrain_a = TRUE,
                        lower_a = 0,
                        lower_log_sigma_theta = -Inf,
                        upper_log_sigma_theta = Inf,
                        lower_log_r = -Inf,
                        upper_log_r = Inf) {
  lower <- rep(-Inf, length(par))
  upper <- rep( Inf, length(par))
  names(lower) <- names(upper) <- names(par)
  
  # a is optimized on the natural scale, so positivity must be imposed directly.
  if (constrain_a && "a" %in% names(par)) {
    lower["a"] <- lower_a
  }
  
  # log_r is optimized on the log scale. No positivity bound is needed for r,
  # but optional numerical bounds can be imposed on log_r if desired.
  if ("log_r" %in% names(par)) {
    lower["log_r"] <- lower_log_r
    upper["log_r"] <- upper_log_r
  }
  
  # log_sigma_theta is optimized on the log scale
  if ("log_sigma_theta" %in% names(par)) {
    lower["log_sigma_theta"] <- lower_log_sigma_theta
    upper["log_sigma_theta"] <- upper_log_sigma_theta
  }
  
  list(lower = lower, upper = upper)
}

get_named <- function(x, name, default = NA_real_) {
  if (!is.null(x) && name %in% names(x)) unname(x[name]) else default
}

safe_hessian_diagnostics <- function(par, fn, gr, eig_tol = 1e-8) {
  H <- tryCatch(
    optimHess(par, fn, gr),
    error = function(e) e
  )
  
  if (inherits(H, "error")) {
    return(list(
      H = NULL,
      ok = FALSE,
      message = conditionMessage(H),
      n_par = length(par),
      eig = numeric(0),
      min_eig = NA_real_,
      max_eig = NA_real_,
      min_abs_eig = NA_real_,
      max_abs_eig = NA_real_,
      abs_cond = NA_real_,
      signed_cond_ratio = NA_real_,
      log10_abs_cond = NA_real_,
      n_negative_eig = NA_integer_,
      n_near_zero_eig = NA_integer_,
      positive_definite = NA,
      near_singular = NA,
      indefinite = NA,
      status = "Hessian failed"
    ))
  }
  
  if (is.null(H) || any(!is.finite(H))) {
    return(list(
      H = H,
      ok = FALSE,
      message = "Hessian is NULL or contains non-finite values",
      n_par = length(par),
      eig = numeric(0),
      min_eig = NA_real_,
      max_eig = NA_real_,
      min_abs_eig = NA_real_,
      max_abs_eig = NA_real_,
      abs_cond = NA_real_,
      signed_cond_ratio = NA_real_,
      log10_abs_cond = NA_real_,
      n_negative_eig = NA_integer_,
      n_near_zero_eig = NA_integer_,
      positive_definite = NA,
      near_singular = NA,
      indefinite = NA,
      status = "Hessian non-finite"
    ))
  }
  
  eig <- tryCatch(
    eigen(H, symmetric = TRUE, only.values = TRUE)$values,
    error = function(e) rep(NA_real_, nrow(H))
  )
  
  if (any(!is.finite(eig))) {
    return(list(
      H = H,
      ok = FALSE,
      message = "Eigenvalue computation failed or returned non-finite values",
      n_par = length(par),
      eig = eig,
      min_eig = NA_real_,
      max_eig = NA_real_,
      min_abs_eig = NA_real_,
      max_abs_eig = NA_real_,
      abs_cond = NA_real_,
      signed_cond_ratio = NA_real_,
      log10_abs_cond = NA_real_,
      n_negative_eig = NA_integer_,
      n_near_zero_eig = NA_integer_,
      positive_definite = NA,
      near_singular = NA,
      indefinite = NA,
      status = "Eigen failed"
    ))
  }
  
  min_eig <- min(eig)
  max_eig <- max(eig)
  abs_eig <- abs(eig)
  min_abs_eig <- min(abs_eig)
  max_abs_eig <- max(abs_eig)
  
  abs_cond <- if (min_abs_eig > 0) {
    max_abs_eig / min_abs_eig
  } else {
    Inf
  }
  
  signed_cond_ratio <- if (min_eig != 0) {
    max_eig / min_eig
  } else {
    NA_real_
  }
  
  n_negative_eig <- sum(eig < -eig_tol)
  n_near_zero_eig <- sum(abs(eig) <= eig_tol)
  
  positive_definite <- all(eig > eig_tol)
  near_singular <- any(abs(eig) <= eig_tol)
  indefinite <- any(eig < -eig_tol)
  
  status <- dplyr::case_when(
    indefinite ~ "Indefinite",
    near_singular ~ "Near-singular",
    positive_definite ~ "Positive definite",
    TRUE ~ "Borderline"
  )
  
  list(
    H = H,
    ok = TRUE,
    message = NA_character_,
    n_par = length(par),
    eig = eig,
    min_eig = min_eig,
    max_eig = max_eig,
    min_abs_eig = min_abs_eig,
    max_abs_eig = max_abs_eig,
    abs_cond = abs_cond,
    signed_cond_ratio = signed_cond_ratio,
    log10_abs_cond = if (is.finite(abs_cond) && abs_cond > 0) log10(abs_cond) else NA_real_,
    n_negative_eig = n_negative_eig,
    n_near_zero_eig = n_near_zero_eig,
    positive_definite = positive_definite,
    near_singular = near_singular,
    indefinite = indefinite,
    status = status
  )
}

safe_sdreport <- function(obj, opt_par) {
  sdr <- tryCatch(sdreport(obj, par.fixed = opt_par), error = function(e) e)
  
  list(
    sdr = if (!inherits(sdr, "error")) sdr else NULL,
    ok = !inherits(sdr, "error"),
    message = if (!inherits(sdr, "error")) NA_character_ else conditionMessage(sdr)
  )
}

infer_theta_hat <- function(theta_raw_hat, sigma_theta_hat) {
  if (is.null(theta_raw_hat)) return(NULL)
  if (!is.finite(sigma_theta_hat)) return(rep(NA_real_, length(theta_raw_hat)))
  sigma_theta_hat * theta_raw_hat
}

predict_eta_from_inference <- function(model, x_hat, log_r_hat, a_hat = NA_real_,
                                       theta_hat = NULL, data) {
  if (is.null(x_hat)) {
    return(rep(NA_real_, length(data$z)))
  }
  
  eta <- switch(
    model,
    
    null = {
      bt_eta(
        x = x_hat,
        i = data$i,
        j = data$j,
        log_r = log_r_hat
      )
    },
    
    full = {
      lid_eta(
        x = x_hat,
        i = data$i,
        j = data$j,
        log_r = log_r_hat,
        a = a_hat
      )
    },
    
    only_theta = {
      base_eta <- bt_eta(
        x = x_hat,
        i = data$i,
        j = data$j,
        log_r = log_r_hat
      )
      
      if (!is.null(theta_hat)) {
        base_eta + theta_hat
      } else {
        base_eta
      }
    },
    
    full_theta = {
      base_eta <- lid_eta(
        x = x_hat,
        i = data$i,
        j = data$j,
        log_r = log_r_hat,
        a = a_hat
      )
      
      if (!is.null(theta_hat)) {
        base_eta + theta_hat
      } else {
        base_eta
      }
    },
    
    stop("Unknown model: ", model)
  )
  
  as.numeric(eta)
}

safe_random_effects <- function(obj, opt_par, random_names) {
  tryCatch({
    
    # Force update of conditional modes at opt_par
    invisible(obj$fn(opt_par))
    
    # Full parameter vector: fixed + conditional modes of random effects
    full_par <- obj$env$last.par.best
    
    # Structured parameter list matching make_pars()
    par_list <- obj$env$parList(full_par)
    
    # Extract only declared random effects
    random_list <- par_list[random_names]
    
    list(
      ok = TRUE,
      full_par = full_par,
      par_list = par_list,
      random = random_list,
      message = NA_character_
    )
    
  }, error = function(e) {
    list(
      ok = FALSE,
      full_par = NULL,
      par_list = NULL,
      random = NULL,
      message = conditionMessage(e)
    )
  })
}


build_inference_object <- function(model, data, opt_par, re) {
  log_r_hat <- get_named(opt_par, "log_r")
  a_hat <- get_named(opt_par, "a")
  log_sigma_theta_hat <- get_named(opt_par, "log_sigma_theta")
  
  r_hat <- if (is.finite(log_r_hat)) exp(log_r_hat) else NA_real_
  sigma_theta_hat <- if (is.finite(log_sigma_theta_hat)) {
    exp(log_sigma_theta_hat)
  } else {
    NA_real_
  }
  
  x_hat <- if (!is.null(re$random) && "x" %in% names(re$random)) {
    re$random$x
  } else {
    NULL
  }
  
  theta_raw_hat <- if (!is.null(re$random) && "theta_raw" %in% names(re$random)) {
    as.numeric(re$random$theta_raw)
  } else {
    NULL
  }
  
  theta_hat <- infer_theta_hat(
    theta_raw_hat = theta_raw_hat,
    sigma_theta_hat = sigma_theta_hat
  )
  
  eta_hat <- predict_eta_from_inference(
    model = model,
    x_hat = x_hat,
    log_r_hat = log_r_hat,
    a_hat = a_hat,
    theta_hat = theta_hat,
    data = data
  )
  
  p_hat <- plogis(eta_hat)
  
  dyad_data <- tibble(
    i = data$i,
    j = data$j,
    z = data$z,
    s = data$s,
    p_sat = data$z / data$s,
    eta_hat = eta_hat,
    p_hat = p_hat
  )
  
  if (!is.null(theta_raw_hat)) {
    dyad_data <- dyad_data |>
      mutate(
        theta_raw_hat = theta_raw_hat,
        theta_hat = theta_hat
      )
  }
  
  list(
    model = model,
    
    fixed = list(
      log_r_hat = log_r_hat,
      r_hat = r_hat,
      a_hat = a_hat,
      log_sigma_theta_hat = log_sigma_theta_hat,
      sigma_theta_hat = sigma_theta_hat
    ),
    
    random_ok = re$ok,
    random_message = re$message,
    
    x_hat = x_hat,
    theta_raw_hat = theta_raw_hat,
    theta_hat = theta_hat,
    
    eta_hat = eta_hat,
    p_hat = p_hat,
    dyad_data = dyad_data,
    
    par_full_hat = re$full_par,
    par_list_hat = re$par_list,
    random_hat = re$random
  )
}

################################################
# Fitting
################################################

fit_one_start <- function(data, spec, start_row,
                          control = list(eval.max = 1000, iter.max = 1000)) {
  tryCatch({
    
    parms <- make_pars(
      data = data,
      model = spec$model,
      a_start = coalesce(start_row$a_start, 0.1),
      log_sigma_theta_start = coalesce(start_row$log_sigma_theta_start, -2)
    )
    
    obj <- RTMB::MakeADFun(
      cmb(spec$func[[1]], data),
      parms,
      random = spec$random[[1]]
    )
    
    bounds <- make_bounds(
      par = obj$par,
      constrain_a = spec$constrain_a
    )
    
    opt <- nlminb(
      start = obj$par,
      objective = obj$fn,
      gradient = obj$gr,
      lower = bounds$lower,
      upper = bounds$upper,
      control = control
    )
    
    nll <- obj$fn(opt$par)
    logLik <- -nll
    grad <- obj$gr(opt$par)
    
    hd <- safe_hessian_diagnostics(opt$par, obj$fn, obj$gr)
    sr <- safe_sdreport(obj, opt$par)
    
    re <- safe_random_effects(
      obj = obj,
      opt_par = opt$par,
      random_names = spec$random[[1]]
    )
    
    inference <- build_inference_object(
      model = spec$model,
      data = data,
      opt_par = opt$par,
      re = re
    )
    
    log_r_hat <- inference$fixed$log_r_hat
    r_hat <- inference$fixed$r_hat
    a_hat <- inference$fixed$a_hat
    log_sigma_theta_hat <- inference$fixed$log_sigma_theta_hat
    sigma_theta_hat <- inference$fixed$sigma_theta_hat
    
    list(
      ok = is.finite(logLik),
      model = spec$model,
      
      obj = obj,
      opt = opt,
      sdr = sr$sdr,
      sdr_ok = sr$ok,
      sdr_message = sr$message,
      
      inference = inference,
      
      # Convenience aliases
      x_hat = inference$x_hat,
      theta_raw_hat = inference$theta_raw_hat,
      theta_hat = inference$theta_hat,
      eta_hat = inference$eta_hat,
      p_hat = inference$p_hat,
      dyad_data = inference$dyad_data,
      
      random_ok = inference$random_ok,
      random_message = inference$random_message,
      par_full_hat = inference$par_full_hat,
      par_list_hat = inference$par_list_hat,
      random_hat = inference$random_hat,
      
      logLik = logLik,
      nll = nll,
      AIC = -2 * logLik + 2 * spec$k_aic,
      k_aic = spec$k_aic,
      
      log_r_hat = log_r_hat,
      r_hat = r_hat,
      a_hat = a_hat,
      a_on_boundary = is.finite(a_hat) && abs(a_hat) < 1e-8,
      log_sigma_theta_hat = log_sigma_theta_hat,
      sigma_theta_hat = sigma_theta_hat,
      
      max_grad = max(abs(grad)),
      grad_log_r = get_named(grad, "log_r"),
      grad_a = get_named(grad, "a"),
      grad_log_sigma_theta = get_named(grad, "log_sigma_theta"),
      
      fixed_hessian = hd$H,
      hess_ok = hd$ok,
      hess_message = hd$message,
      hess_n_par = hd$n_par,
      hess_eigenvalues = hd$eig,
      
      hess_min_eig = hd$min_eig,
      hess_max_eig = hd$max_eig,
      hess_min_abs_eig = hd$min_abs_eig,
      hess_max_abs_eig = hd$max_abs_eig,
      
      hess_abs_cond = hd$abs_cond,
      hess_signed_cond_ratio = hd$signed_cond_ratio,
      hess_log10_abs_cond = hd$log10_abs_cond,
      
      hess_n_negative_eig = hd$n_negative_eig,
      hess_n_near_zero_eig = hd$n_near_zero_eig,
      hess_positive_definite = hd$positive_definite,
      hess_near_singular = hd$near_singular,
      hess_indefinite = hd$indefinite,
      hess_status = hd$status,
      
      convergence = opt$convergence,
      message = opt$message,
      
      a_start = start_row$a_start,
      log_sigma_theta_start = start_row$log_sigma_theta_start
    )
    
  }, error = function(e) {
    list(
      ok = FALSE,
      model = spec$model,
      
      obj = NULL,
      opt = NULL,
      sdr = NULL,
      sdr_ok = FALSE,
      sdr_message = NA_character_,
      
      inference = NULL,
      
      x_hat = NULL,
      theta_raw_hat = NULL,
      theta_hat = NULL,
      eta_hat = NULL,
      p_hat = NULL,
      dyad_data = NULL,
      
      random_ok = FALSE,
      random_message = conditionMessage(e),
      par_full_hat = NULL,
      par_list_hat = NULL,
      random_hat = NULL,
      
      logLik = NA_real_,
      nll = NA_real_,
      AIC = NA_real_,
      k_aic = spec$k_aic,
      
      log_r_hat = NA_real_,
      r_hat = NA_real_,
      a_hat = NA_real_,
      a_on_boundary = NA,
      log_sigma_theta_hat = NA_real_,
      sigma_theta_hat = NA_real_,
      
      max_grad = NA_real_,
      grad_log_r = NA_real_,
      grad_a = NA_real_,
      grad_log_sigma_theta = NA_real_,
      
      fixed_hessian = NULL,
      hess_ok = FALSE,
      hess_message = conditionMessage(e),
      hess_n_par = NA_integer_,
      hess_eigenvalues = numeric(0),
      
      hess_min_eig = NA_real_,
      hess_max_eig = NA_real_,
      hess_min_abs_eig = NA_real_,
      hess_max_abs_eig = NA_real_,
      
      hess_abs_cond = NA_real_,
      hess_signed_cond_ratio = NA_real_,
      hess_log10_abs_cond = NA_real_,
      
      hess_n_negative_eig = NA_integer_,
      hess_n_near_zero_eig = NA_integer_,
      hess_positive_definite = NA,
      hess_near_singular = NA,
      hess_indefinite = NA,
      hess_status = "Fit failed",
      
      convergence = NA_integer_,
      message = conditionMessage(e),
      
      a_start = start_row$a_start,
      log_sigma_theta_start = start_row$log_sigma_theta_start
    )
  })
}
fit_model_multistart <- function(data, spec,
                                 control = list(eval.max = 1000, iter.max = 1000)) {
  
  fits <- pmap(
    spec$start_grid[[1]],
    ~ fit_one_start(
      data = data,
      spec = spec,
      start_row = tibble(
        a_start = ..1,
        log_sigma_theta_start = ..2
      ),
      control = control
    )
  )
  
  ok_fits <- keep(fits, ~ isTRUE(.x$ok))
  
  best <- if (length(ok_fits) > 0) {
    ok_fits[[which.max(map_dbl(ok_fits, "logLik"))]]
  } else {
    NULL
  }
  
  start_diag <- map_dfr(fits, \(fit) {
    tibble(
      model = fit$model,
      ok = fit$ok,
      logLik = fit$logLik,
      AIC = fit$AIC,
      convergence = fit$convergence,
      message = fit$message,
      max_grad = fit$max_grad,
      hess_min_eig = fit$hess_min_eig,
      hess_cond = fit$hess_cond,
      random_ok = fit$random_ok,
      random_message = fit$random_message,
      a_start = fit$a_start,
      log_sigma_theta_start = fit$log_sigma_theta_start,
      log_r_hat = fit$log_r_hat,
      r_hat = fit$r_hat,
      a_hat = fit$a_hat,
      sigma_theta_hat = fit$sigma_theta_hat
    )
  }) |>
    arrange(desc(logLik))
  
  list(
    ok = !is.null(best),
    model = spec$model,
    
    best = best,
    tries = fits,
    start_diag = start_diag,
    
    best_inference = if (!is.null(best)) best$inference else NULL,
    best_dyad_data = if (!is.null(best)) best$dyad_data else NULL,
    
    logLik = if (!is.null(best)) best$logLik else NA_real_,
    AIC = if (!is.null(best)) best$AIC else NA_real_,
    
    n_starts = length(fits),
    n_ok = length(ok_fits),
    logLik_range = if (length(ok_fits) > 1) {
      diff(range(map_dbl(ok_fits, "logLik"), na.rm = TRUE))
    } else {
      NA_real_
    },
    
    message = if (!is.null(best)) NA_character_ else paste("All starts failed for", spec$model)
  )
}

# ============================================================
# 7. Extract summaries
# ============================================================

safe_sdr_summary <- function(sdr) {
  if (is.null(sdr)) {
    return(tibble())
  }
  
  ss <- tryCatch(
    as.data.frame(summary(sdr)),
    error = function(e) tibble()
  )
  
  if (nrow(ss) == 0) {
    return(tibble())
  }
  
  ss |>
    tibble::rownames_to_column("parameter") |>
    as_tibble() |>
    rename(
      estimate = Estimate,
      se = `Std. Error`
    ) |>
    mutate(
      z_value = estimate / se
    )
}

extract_sdr_estimates <- function(fit, prefix) {
  if (is.null(fit) || is.null(fit$sdr)) {
    return(tibble())
  }
  
  ss <- safe_sdr_summary(fit$sdr)
  
  if (nrow(ss) == 0) {
    return(tibble())
  }
  
  keep_pars <- c("log_r", "a", "log_sigma_theta")
  
  out <- ss |>
    filter(parameter %in% keep_pars) |>
    select(parameter, estimate, se, z_value) |>
    pivot_wider(
      names_from = parameter,
      values_from = c(estimate, se, z_value),
      names_glue = paste0(prefix, "_{parameter}_sdr_{.value}")
    )
  
  if ("log_r" %in% ss$parameter) {
    est <- ss$estimate[ss$parameter == "log_r"][1]
    se <- ss$se[ss$parameter == "log_r"][1]
    
    out <- bind_cols(
      out,
      tibble(
        !!paste0(prefix, "_r_sdr_estimate") := exp(est),
        !!paste0(prefix, "_r_sdr_se_delta") := exp(est) * se
      )
    )
  }
  
  if ("log_sigma_theta" %in% ss$parameter) {
    est <- ss$estimate[ss$parameter == "log_sigma_theta"][1]
    se <- ss$se[ss$parameter == "log_sigma_theta"][1]
    
    out <- bind_cols(
      out,
      tibble(
        !!paste0(prefix, "_sigma_theta_sdr_estimate") := exp(est),
        !!paste0(prefix, "_sigma_theta_sdr_se_delta") := exp(est) * se
      )
    )
  }
  
  out
}

extract_fit_row <- function(fit, prefix) {
  if (is.null(fit)) {
    return(tibble(
      !!paste0(prefix, "_ok") := FALSE,
      !!paste0(prefix, "_logLik") := NA_real_,
      !!paste0(prefix, "_AIC") := NA_real_,
      !!paste0(prefix, "_random_ok") := FALSE,
      !!paste0(prefix, "_random_message") := NA_character_
    ))
  }
  
  base <- tibble(
    !!paste0(prefix, "_ok") := fit$ok,
    !!paste0(prefix, "_logLik") := fit$logLik,
    !!paste0(prefix, "_nll") := fit$nll,
    !!paste0(prefix, "_AIC") := fit$AIC,
    !!paste0(prefix, "_k_aic") := fit$k_aic,
    
    !!paste0(prefix, "_convergence") := fit$convergence,
    !!paste0(prefix, "_message") := fit$message,
    
    !!paste0(prefix, "_max_grad") := fit$max_grad,
    !!paste0(prefix, "_grad_log_r") := fit$grad_log_r,
    !!paste0(prefix, "_grad_a") := fit$grad_a,
    !!paste0(prefix, "_grad_log_sigma_theta") := fit$grad_log_sigma_theta,
    
    !!paste0(prefix, "_hess_ok") := fit$hess_ok,
    !!paste0(prefix, "_hess_message") := fit$hess_message,
    !!paste0(prefix, "_hess_n_par") := fit$hess_n_par,
    
    !!paste0(prefix, "_hess_min_eig") := fit$hess_min_eig,
    !!paste0(prefix, "_hess_max_eig") := fit$hess_max_eig,
    !!paste0(prefix, "_hess_min_abs_eig") := fit$hess_min_abs_eig,
    !!paste0(prefix, "_hess_max_abs_eig") := fit$hess_max_abs_eig,
    
    !!paste0(prefix, "_hess_abs_cond") := fit$hess_abs_cond,
    !!paste0(prefix, "_hess_signed_cond_ratio") := fit$hess_signed_cond_ratio,
    !!paste0(prefix, "_hess_log10_abs_cond") := fit$hess_log10_abs_cond,
    
    !!paste0(prefix, "_hess_n_negative_eig") := fit$hess_n_negative_eig,
    !!paste0(prefix, "_hess_n_near_zero_eig") := fit$hess_n_near_zero_eig,
    !!paste0(prefix, "_hess_positive_definite") := fit$hess_positive_definite,
    !!paste0(prefix, "_hess_near_singular") := fit$hess_near_singular,
    !!paste0(prefix, "_hess_indefinite") := fit$hess_indefinite,
    !!paste0(prefix, "_hess_status") := fit$hess_status,
    
    !!paste0(prefix, "_sdr_ok") := fit$sdr_ok,
    !!paste0(prefix, "_sdr_message") := fit$sdr_message,
    
    !!paste0(prefix, "_random_ok") := fit$random_ok,
    !!paste0(prefix, "_random_message") := fit$random_message,
    
    !!paste0(prefix, "_log_r_est") := fit$log_r_hat,
    !!paste0(prefix, "_r_est") := fit$r_hat,
    
    !!paste0(prefix, "_a_est") := fit$a_hat,
    !!paste0(prefix, "_a_on_boundary") := fit$a_on_boundary,
    
    !!paste0(prefix, "_log_sigma_theta_est") := fit$log_sigma_theta_hat,
    !!paste0(prefix, "_sigma_theta_est") := fit$sigma_theta_hat,
    
    !!paste0(prefix, "_best_a_start") := fit$a_start,
    !!paste0(prefix, "_best_log_sigma_theta_start") := fit$log_sigma_theta_start
  )
  
  bind_cols(base, extract_sdr_estimates(fit, prefix))
}

extract_latent_positions <- function(fit, dataset_name, model_name) {
  if (is.null(fit) || is.null(fit$x_hat)) {
    return(tibble())
  }
  
  x <- fit$x_hat
  
  as_tibble(x, .name_repair = ~ paste0("x", seq_along(.x))) |>
    mutate(
      individual = row_number(),
      dataset_name = dataset_name,
      model = model_name,
      .before = 1
    )
}

extract_dyad_effects <- function(fit, dataset_name, model_name) {
  if (is.null(fit) || is.null(fit$dyad_data)) {
    return(tibble())
  }
  
  fit$dyad_data |>
    select(
      i, j, z, s, p_sat, eta_hat, p_hat,
      any_of(c("theta_raw_hat", "theta_hat"))
    ) |>
    mutate(
      dataset_name = dataset_name,
      model = model_name,
      .before = 1
    )
}

# ============================================================
# 8. Fit one dataset
# ============================================================
results <- results |>
  mutate(
    full_hess_status_clean = case_when(
      is.na(full_ok) | !full_ok ~ "Fit failed",
      is.na(full_hess_ok) ~ "Hessian status missing",
      full_hess_ok == FALSE & grepl("non-finite", full_hess_message, ignore.case = TRUE) ~ "Hessian non-finite",
      full_hess_ok == FALSE ~ "Hessian failed",
      is.na(full_hess_min_eig) ~ "Eigenvalue missing",
      full_hess_indefinite %in% TRUE ~ "Indefinite",
      full_hess_near_singular %in% TRUE ~ "Near-singular",
      full_hess_positive_definite %in% TRUE ~ "Positive definite",
      TRUE ~ "Borderline"
    ), 
    full_theta_hess_status_clean = case_when(
      is.na(full_theta_ok) | !full_theta_ok ~ "Fit failed",
      is.na(full_theta_hess_ok) ~ "Hessian status missing",
      full_theta_hess_ok == FALSE & grepl("non-finite", full_hess_message, ignore.case = TRUE) ~ "Hessian non-finite",
      full_theta_hess_ok == FALSE ~ "Hessian failed",
      is.na(full_theta_hess_min_eig) ~ "Eigenvalue missing",
      full_theta_hess_indefinite %in% TRUE ~ "Indefinite",
      full_theta_hess_near_singular %in% TRUE ~ "Near-singular",
      full_theta_hess_positive_definite %in% TRUE ~ "Positive definite",
      TRUE ~ "Borderline"
  ))
results |>
  count(full_hess_status_clean)
lrt <- function(logLik0, logLik1) {
  ifelse(
    is.finite(logLik0) & is.finite(logLik1),
    pmax(2 * (logLik1 - logLik0), 0),
    NA_real_
  )
}
results%>%
  filter(is.na(full_theta_logLik))
fit_dataset <- function(name, dataset,
                        specs = model_specs,
                        control = list(eval.max = 1000, iter.max = 1000)) {
  M <- dataset$matrix
  
  if (!is.matrix(M)) {
    message("Skipping ", name, ": missing matrix")
    return(NULL)
  }
  
  data <- tryCatch(
    make_pair_data(M, drop_zero = TRUE),
    error = function(e) {
      message("Skipping ", name, ": ", conditionMessage(e))
      NULL
    }
  )
  
  if (is.null(data)) return(NULL)
  
  message("Fitting: ", name)
  
  fits_tbl <- specs |>
    mutate(fit = pmap(
      list(model, func, random, constrain_a, k_aic, start_grid),
      function(model, func, random, constrain_a, k_aic, start_grid) {
        spec <- tibble(
          model = model,
          func = list(func),
          random = list(random),
          constrain_a = constrain_a,
          k_aic = k_aic,
          start_grid = list(start_grid)
        )
        
        fit_model_multistart(
          data = data,
          spec = spec,
          control = control
        )
      }
    ))
  
  fits <- set_names(fits_tbl$fit, fits_tbl$model)
  best_fits <- map(fits, "best")
  
  meta <- as_tibble(dataset$metadata)
  
  model_rows <- imap(
    best_fits,
    ~ extract_fit_row(.x, .y)
  ) |>
    bind_cols()
  
  ll <- map_dbl(best_fits, ~ if (is.null(.x)) NA_real_ else .x$logLik)
  aa <- map_dbl(best_fits, ~ if (is.null(.x)) NA_real_ else .x$AIC)
  
  comparison_row <- tibble(
    dataset_name = name,
    n_agents = data$n,
    n_observed_dyads = data$n_dyads_observed,
    n_interactions = data$n_interactions,
    interactions_per_observed_dyad = data$n_interactions / data$n_dyads_observed,
    
    lrt_null_vs_full = lrt(ll["null"], ll["full"]),
    lrt_null_vs_only_theta = lrt(ll["null"], ll["only_theta"]),
    lrt_null_vs_full_theta = lrt(ll["null"], ll["full_theta"]),
    lrt_full_vs_full_theta = lrt(ll["full"], ll["full_theta"]),
    lrt_only_theta_vs_full_theta = lrt(ll["only_theta"], ll["full_theta"]),
    
    delta_aic_full_vs_null = aa["full"] - aa["null"],
    delta_aic_only_theta_vs_null = aa["only_theta"] - aa["null"],
    delta_aic_full_theta_vs_null = aa["full_theta"] - aa["null"],
    delta_aic_full_theta_vs_full = aa["full_theta"] - aa["full"],
    delta_aic_full_theta_vs_only_theta = aa["full_theta"] - aa["only_theta"],
    delta_aic_full_vs_only_theta = aa["full"] - aa["only_theta"],
    
    best_aic_model = if (all(is.na(aa))) NA_character_ else names(which.min(aa)),
    best_aic = if (all(is.na(aa))) NA_real_ else min(aa, na.rm = TRUE)
  )
  
  multistart_row <- imap_dfr(fits, \(fit, model_name) {
    tibble(
      model = model_name,
      n_starts = fit$n_starts,
      n_ok = fit$n_ok,
      logLik_range = fit$logLik_range
    )
  }) |>
    pivot_wider(
      names_from = model,
      values_from = c(n_starts, n_ok, logLik_range),
      names_glue = "{model}_{.value}"
    )
  
  result_row <- bind_cols(meta, comparison_row, model_rows, multistart_row) |>
    mutate(across(where(is.numeric), ~ ifelse(is.infinite(.x), NA_real_, .x)))
  
  start_diag <- map_dfr(fits, \(fit) {
    tibble(
      model = fit$model,
      ok = fit$ok,
      logLik = fit$logLik,
      nll = fit$nll,
      AIC = fit$AIC,
      convergence = fit$convergence,
      message = fit$message,
      
      max_grad = fit$max_grad,
      grad_log_r = fit$grad_log_r,
      grad_a = fit$grad_a,
      grad_log_sigma_theta = fit$grad_log_sigma_theta,
      
      hess_ok = fit$hess_ok,
      hess_message = fit$hess_message,
      hess_n_par = fit$hess_n_par,
      hess_min_eig = fit$hess_min_eig,
      hess_max_eig = fit$hess_max_eig,
      hess_min_abs_eig = fit$hess_min_abs_eig,
      hess_max_abs_eig = fit$hess_max_abs_eig,
      hess_abs_cond = fit$hess_abs_cond,
      hess_signed_cond_ratio = fit$hess_signed_cond_ratio,
      hess_log10_abs_cond = fit$hess_log10_abs_cond,
      hess_n_negative_eig = fit$hess_n_negative_eig,
      hess_n_near_zero_eig = fit$hess_n_near_zero_eig,
      hess_positive_definite = fit$hess_positive_definite,
      hess_near_singular = fit$hess_near_singular,
      hess_indefinite = fit$hess_indefinite,
      hess_status = fit$hess_status,
      
      sdr_ok = fit$sdr_ok,
      sdr_message = fit$sdr_message,
      
      random_ok = fit$random_ok,
      random_message = fit$random_message,
      
      a_start = fit$a_start,
      log_sigma_theta_start = fit$log_sigma_theta_start,
      
      log_r_hat = fit$log_r_hat,
      r_hat = fit$r_hat,
      a_hat = fit$a_hat,
      a_on_boundary = fit$a_on_boundary,
      log_sigma_theta_hat = fit$log_sigma_theta_hat,
      sigma_theta_hat = fit$sigma_theta_hat
    )
  }) |>
    arrange(desc(logLik))
  
  fitted_probs <- imap_dfr(best_fits, \(fit, model_name) {
    if (is.null(fit) || is.null(fit$dyad_data)) {
      return(tibble())
    }
    
    fit$dyad_data |>
      mutate(
        dataset_name = name,
        model = model_name,
        .before = 1
      )
  })
  latent_positions <- imap_dfr(best_fits, \(fit, model_name) {
    extract_latent_positions(
      fit = fit,
      dataset_name = name,
      model_name = model_name
    )
  })
  
  dyad_effects <- imap_dfr(best_fits, \(fit, model_name) {
    extract_dyad_effects(
      fit = fit,
      dataset_name = name,
      model_name = model_name
    )
  })
  
  list(
    result = result_row,
    start_diag = start_diag,
    fitted_probs = fitted_probs,
    latent_positions = latent_positions,
    dyad_effects = dyad_effects,
    pair_data = data
  )
}
# ============================================================
# 9. Run all datasets
# ============================================================
# ============================================================
#   a. Paths
# ============================================================

out_dir <- "DomArchiveResults"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

tag <- "four_models_comparison_filtered_betterdiag"

paths <- list(
  results          = file.path(out_dir, paste0(tag, "_results.rds")),
  start_diag       = file.path(out_dir, paste0(tag, "_start_diagnostics.rds")),
  fitted_probs     = file.path(out_dir, paste0(tag, "_fitted_probs.rds")),
  latent_positions = file.path(out_dir, paste0(tag, "_latent_positions.rds")),
  dyad_effects     = file.path(out_dir, paste0(tag, "_dyad_effects.rds")),
  agreement        = file.path(out_dir, paste0(tag, "_probability_agreement.rds")),
  outputs_light    = file.path(out_dir, paste0(tag, "_outputs_light.rds"))
)
# ============================================================
#   b. Filter DomArchive datasets
# ============================================================

dom.data_filt <- dom.metadata |>
  filter(
    countbinary == "Count",
    matrix_edgelist == "Matrix"
  )

dom.data_filtered <- dom.data[names(dom.data) %in% dom.data_filt$fileid]
length(dom.data_filtered)
# ============================================================
#   c. Fit all datasets
# ============================================================

all_outputs <- imap(
  dom.data_filtered,
  \(dataset, dataset_name) {
    fit_dataset(
      name = dataset_name,
      dataset = dataset,
      specs = model_specs,
      control = list(eval.max = 1000, iter.max = 1000)
    )
  }
) |>
  compact()
# ============================================================
#   d. Extract result tables and fit objects
# ============================================================

results <- map_dfr(all_outputs, "result")

start_diag <- map_dfr(all_outputs, "start_diag")

fitted_probs <- map_dfr(all_outputs, "fitted_probs")

latent_positions <- map_dfr(all_outputs, "latent_positions")

dyad_effects <- map_dfr(all_outputs, "dyad_effects")
# ============================================================
#   e. Probability agreement metrics
# ============================================================

calc_probability_agreement <- function(fitted_probs) {
  fitted_probs |>
    group_by(dataset_name, model) |>
    summarise(
      n_dyads = n(),
      total_interactions = sum(s, na.rm = TRUE),
      
      D_WMAD = sum(
        s * abs(p_hat - p_sat),
        na.rm = TRUE
      ) / sum(s, na.rm = TRUE),
      
      weighted_frobenius_norm = sqrt(
        sum(
          s * (p_hat - p_sat)^2 +
            s * ((1 - p_hat) - (1 - p_sat))^2,
          na.rm = TRUE
        ) / sum(2 * s, na.rm = TRUE)
      ),
      
      weighted_frobenius_norm_raw = sqrt(
        sum(
          s * (p_hat - p_sat)^2 +
            s * ((1 - p_hat) - (1 - p_sat))^2,
          na.rm = TRUE
        )
      ),
      
      .groups = "drop"
    )
}

probability_agreement <- calc_probability_agreement(fitted_probs)
agreement_wide <- probability_agreement |>
  pivot_wider(
    names_from = model,
    values_from = c(
      D_WMAD,
      weighted_frobenius_norm,
      weighted_frobenius_norm_raw
    ),
    names_glue = "{model}_{.value}"
  )

results <- results |>
  left_join(agreement_wide, by = "dataset_name")
# ============================================================
#   f. Save outputs
# ============================================================

saveRDS(results, paths$results)
saveRDS(start_diag, paths$start_diag)
saveRDS(fitted_probs, paths$fitted_probs)
saveRDS(latent_positions, paths$latent_positions)
saveRDS(dyad_effects, paths$dyad_effects)
saveRDS(probability_agreement, paths$agreement)
saveRDS(all_outputs, paths$outputs_light)


highlight_fileids <- c(
  "Alados_1992b",
  "Adcock_2015a",
  "Blatrix_2004c",
  "Poisbleau_2005c",
  "Shimoji_2014c", 
  "Cote_2000d",
  "Correa_2013a",
  "Kolodziejczyk_2005",
  "Cui_2014",
  "Prieto_1978",
  "Mwamende_2009a",
  "ScottLockhard_1999b"
)

selected_colours <- c(
  "Other" = "grey75", 
  "Correa_2013a" = "brown",
  "Alados_1992b" = "red", 
  "Adcock_2015a" = "orange", 
  "Blatrix_2004c" = "yellow", 
  "Poisbleau_2005c" = "green", 
  "Shimoji_2014c" = "blue", 
  "Kolodziejczyk_2005" = "purple", 
  "Cote_2000d" = "violet", 
  "Cui_2014" = "darkblue", 
  "Prieto_1978" = "black",
  "Mwamende_2009a" = "magenta",
  "ScottLockhard_1999b" = "cyan"
)
# =============================================================
# 10. Tables 
# =============================================================

escape_latex <- function(x) {
  x %>%
    str_replace_all("\\\\", "\\\\textbackslash{}") %>%
    str_replace_all("_", "\\\\_") %>%
    str_replace_all("&", "\\\\&") %>%
    str_replace_all("%", "\\\\%") %>%
    str_replace_all("#", "\\\\#")
}

fmt_num <- function(x, digits = 2) {
  ifelse(
    is.na(x) | is.nan(x) | is.infinite(x),
    "--",
    formatC(x, digits = digits, format = "f", drop0trailing = TRUE)
  )
}

lrt_aic_table <- results %>%
  filter(fileid %in% highlight_fileids) %>%
  arrange(fileid) %>%
  transmute(
    Dataset = escape_latex(fileid),
    
    `LRT: Null vs Full` = fmt_num(lrt_null_vs_full),
    `LRT: Null vs Only theta` = fmt_num(lrt_null_vs_only_theta),
    `LRT: Full vs Full theta` = fmt_num(lrt_full_vs_full_theta),
    `LRT: Only theta vs Full theta` = fmt_num(lrt_only_theta_vs_full_theta),
    
    `AIC Null` = fmt_num(null_AIC),
    `AIC Full` = fmt_num(full_AIC),
    `AIC Only theta` = fmt_num(only_theta_AIC),
    `AIC Full theta` = fmt_num(full_theta_AIC)
  )

latex_rows <- lrt_aic_table %>%
  mutate(
    row = paste(
      Dataset,
      `LRT: Null vs Full`,
      `LRT: Null vs Only theta`,
      `LRT: Full vs Full theta`,
      `LRT: Only theta vs Full theta`,
      `AIC Null`,
      `AIC Full`,
      `AIC Only theta`,
      `AIC Full theta`,
      sep = " & "
    ),
    row = paste0(row, " \\\\")
  ) %>%
  pull(row)

latex_code <- paste(
  "\\begin{table}[htbp]",
  "\\centering",
  "\\scriptsize",
  "\\setlength{\\tabcolsep}{3pt}",
  "\\caption{Likelihood ratio statistics and AIC scores for the selected datasets.}",
  "\\label{tab:highlight_lrt_aic}",
  "\\resizebox{\\textwidth}{!}{%",
  "\\begin{tabular}{lcccccccc}",
  "\\toprule",
  "Dataset & $\\Lambda_{N,F}$ & $\\Lambda_{N,O}$ & $\\Lambda_{F,F_\\theta}$ & $\\Lambda_{O,F_\\theta}$ & AIC$_N$ & AIC$_F$ & AIC$_O$ & AIC$_{F_\\theta}$ \\\\",
  "\\midrule",
  paste(latex_rows, collapse = "\n"),
  "\\bottomrule",
  "\\end{tabular}%",
  "}",
  "\\end{table}",
  sep = "\n"
)

hessian_diag_table <- results %>%
  filter(fileid %in% highlight_fileids) %>%
  arrange(fileid) %>%
  select(
    fileid,
    
    null_hess_min_eig,
    null_hess_abs_cond,
    null_max_grad,
    null_hess_status,
    
    full_hess_min_eig,
    full_hess_abs_cond,
    full_max_grad,
    full_hess_status,
    
    only_theta_hess_min_eig,
    only_theta_hess_abs_cond,
    only_theta_max_grad,
    only_theta_hess_status,
    
    full_theta_hess_min_eig,
    full_theta_hess_abs_cond,
    full_theta_max_grad,
    full_theta_hess_status
  ) %>%
  pivot_longer(
    cols = -fileid,
    names_to = c("model", ".value"),
    names_pattern = "^(null|full|only_theta|full_theta)_(hess_min_eig|hess_abs_cond|max_grad|hess_status)$"
  ) %>%
  mutate(
    model = factor(model, levels = model_order),
    Dataset = escape_latex(fileid),
    Model = unname(model_labels[as.character(model)]),
    min_eig = fmt_num(hess_min_eig),
    cond_num = fmt_num(hess_abs_cond),
    max_grad = fmt_num(max_grad),
    hess_status = ifelse(
      is.na(hess_status),
      "--",
      escape_latex(as.character(hess_status))
    )
  ) %>%
  arrange(Dataset, model) %>%
  select(
    Dataset,
    Model,
    min_eig,
    cond_num,
    max_grad,
    hess_status
  )

latex_rows <- hessian_diag_table %>%
  mutate(
    row = paste(
      Dataset,
      Model,
      min_eig,
      cond_num,
      max_grad,
      hess_status,
      sep = " & "
    ),
    row = paste0(row, " \\\\")
  ) %>%
  pull(row)

latex_code <- paste(
  "\\begin{table}[htbp]",
  "\\centering",
  "\\scriptsize",
  "\\setlength{\\tabcolsep}{4pt}",
  "\\caption{Hessian diagnostics for the selected datasets and fitted models. The table reports the minimum eigenvalue of the fixed-parameter Hessian, the absolute condition number, the maximum absolute gradient component, and the Hessian status.}",
  "\\label{tab:highlight_hessian_diagnostics}",
  "\\resizebox{\\textwidth}{!}{%",
  "\\begin{tabular}{llrrrr}",
  "\\toprule",
  "Dataset & Model & $\\lambda_{\\min}(H)$ & $\\kappa(H)$ & $\\|\\nabla f\\|_\\infty$ & Hessian status \\\\",
  "\\midrule",
  paste(latex_rows, collapse = "\n"),
  "\\bottomrule",
  "\\end{tabular}%",
  "}",
  "\\end{table}",
  sep = "\n"
)

summary_table <- results %>%
  select(
    n_agents,
    interactions_per_observed_dyad,
    proportion_unknown,
    dci,
    ttri,
    modified_landaus_h
  ) %>%
  pivot_longer(
    cols = everything(),
    names_to = "metric",
    values_to = "value"
  ) %>%
  group_by(metric) %>%
  summarise(
    Mean = mean(value, na.rm = TRUE),
    Median = median(value, na.rm = TRUE),
    SD = sd(value, na.rm = TRUE),
    IQR = IQR(value, na.rm = TRUE),
    Min = min(value, na.rm = TRUE),
    Max = max(value, na.rm = TRUE),
    .groups = "drop"
  )
summary_table


# ============================================================
# 11. Figures
# ============================================================

data_explanalysis_ploting <- results %>%
  mutate(
    highlighted = if_else(fileid %in% highlight_fileids, fileid, "Other")
  )
cor_fun <- function(data, mapping, ...) {
  x <- eval_data_col(data, mapping$x)
  y <- eval_data_col(data, mapping$y)
  
  test <- cor.test(x, y, use = "complete.obs")
  r <- test$estimate
  p <- test$p.value
  
  stars <- case_when(
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE      ~ ""
  )
  
  ggplot() +
    annotate(
      "text",
      x = 0.5,
      y = 0.5,
      label = paste0("Corr: ", round(r, 3), stars),
      size = 3.5
    ) +
    xlim(0, 1) +
    ylim(0, 1) +
    theme_void()
}
lower_fun <- function(data, mapping, ...) {
  x_var <- as_label(mapping$x)
  y_var <- as_label(mapping$y)
  
  ggplot(data, aes(x = .data[[x_var]], y = .data[[y_var]])) +
    geom_point(
      data = data %>% filter(highlighted == "Other"),
      color = "grey75",
      alpha = 0.25,
      size = 0.8
    ) +
    geom_point(
      data = data %>% filter(highlighted != "Other"),
      aes(color = highlighted),
      alpha = 1,
      size = 1.8
    )
}

data_explanalysis_ploting %>%
  select(fileid, highlighted, all_of(metrics)) %>%
  ggpairs(
    columns = metrics,
    columnLabels = pretty_labs,
    aes(color = highlighted, alpha = highlighted),
    upper = list(continuous = cor_fun),
    lower = list(continuous = lower_fun
                 ),
    diag = list(
      continuous = wrap("densityDiag", alpha = 0.4)
    )
  ) +
  scale_color_manual(
        values = selected_colours
  ) +
  theme_bw()
p_a <- ggplot(
  plot_df,
  aes(x = full_a_est + eps, y = full_theta_a_est + eps, color = fileid_highlight)
) +
  geom_abline(
    intercept = 0,
    slope = 1,
    linewidth = 0.6,
    linetype = "dashed"
  ) +
  geom_point(
    data = filter(plot_df, fileid_highlight == "Other"),
    colour = "grey75",
    alpha = 0.6,
    size = 2.0
  ) +
  geom_point(
    data = filter(plot_df, fileid_highlight != "Other"),
    aes(colour = fileid_highlight),
    alpha = 1,
    size = 5.2
  ) +
  scale_x_log10() +
  scale_y_log10() +
  coord_equal() +
  scale_colour_manual(
    values = selected_colours,
    guide = "none"
  ) +
  labs(
    x = expression(hat(a)[LID]),
    y = expression(hat(a)[DALID]),
    title = ""
  ) +
  theme_bw(base_size = 28) +
  theme(
    axis.title = element_text(face = "bold"),
    axis.text = element_text(size = 22),
    panel.grid.minor = element_blank()
  )
p_sigma <- ggplot(
  plot_df,
  aes(x = only_theta_sigma_theta_est + eps, y = full_theta_sigma_theta_est + eps, color = fileid_highlight)
) +
  geom_abline(
    intercept = 0,
    slope = 1,
    linewidth = 0.6,
    linetype = "dashed"
  ) +
  geom_point(
    data = filter(plot_df, fileid_highlight == "Other"),
    colour = "grey75",
    alpha = 0.6,
    size = 2.0
  ) +
  geom_point(
    data = filter(plot_df, fileid_highlight != "Other"),
    aes(colour = fileid_highlight),
    alpha = 1,
    size = 5.2
  ) +
  scale_x_log10() +
  scale_y_log10() +
  scale_colour_manual(
    values = selected_colours,
    guide = "none"
  ) +
  labs(
    x = expression(hat(sigma)[theta*","*DABT]),
    y = expression(hat(sigma)[theta*","*DALID]),
    title = ""
  ) +
  theme_bw(base_size = 28) +
  theme(
    axis.title = element_text(face = "bold"),
    axis.text = element_text(size = 22),
    panel.grid.minor = element_blank(),
    aspect.ratio = 1
  )

  # ============================================================
  # Likelihood diagnostics
  # ============================================================

eps <- 1e-8
p_lrt_bltlid <- res |>
  filter(!is.na(full_a_est), !is.na(lrt_null_vs_full)) |>
  mutate(a_plot = full_a_est + eps) |>
  ggplot(aes(x = a_plot, y = lrt_null_vs_full)) +
  geom_point(alpha = 0.55) +
  scale_x_log10() +
  scale_y_continuous(trans = "sqrt") +
  labs(
    x = expression(hat(a)[LID] + epsilon),
    y = expression(Lambda[LBT-LID]),
    title = "Estimated intransitivity coefficient and model improvement LBT vs LID"
  ) +
  theme_bw()

p_lrt_liddalid <- res |>
  filter(!is.na(full_theta_sigma_theta_est), !is.na(lrt_full_vs_full_theta)) |>
  mutate(a_plot = full_theta_sigma_theta_est + eps) |>
  ggplot(aes(x = a_plot, y = lrt_full_vs_full_theta)) +
  geom_point(alpha = 0.55) +
  scale_x_log10() +
  scale_y_continuous(trans = "sqrt") +
  labs(
    x = expression(hat(sigma)[DALID] + epsilon),
    y = expression(Lambda[LID-DALID]),
    title = "Estimated random effect scale and model improvement LID vs DALID"
  ) +
  theme_bw()

p_lrt_lbtdabt <- res |>
  mutate(a_plot = only_theta_sigma_theta_est + eps) |>
  ggplot(aes(x = a_plot, y = lrt_null_vs_only_theta)) +
  geom_point(alpha = 0.55) +
  scale_x_log10() +
  scale_y_continuous(trans = "sqrt") +
  labs(
    x = expression(hat(sigma)[DABT] + epsilon),
    y = expression(Lambda[LBT-DABT]),
    title = "Estimated random effect scale and model improvement LBT vs DABT"
  ) +
  theme_bw()

p_lrt_dabtdalid <- res |>
  mutate(a_plot = full_theta_a_est + eps) |>
  ggplot(aes(x = a_plot, y = lrt_only_theta_vs_full_theta)) +
  geom_point(alpha = 0.55) +
  scale_x_log10() +
  scale_y_continuous(trans = "sqrt") +
  labs(
    x = expression(hat(a)[DALID] + epsilon),
    y = expression(Lambda[DABT-DALID]),
    title = "Estimated intransitivity coefficient and model improvement DABT vs DALID"
  ) +
  theme_bw()

crit_mix_005 <- qchisq(0.90, df = 1)  # 2.705543

lrt_summary <- results %>%
  select(
    lrt_null_vs_full,
    lrt_null_vs_only_theta,
    lrt_full_vs_full_theta,
    lrt_only_theta_vs_full_theta
  ) %>%
  pivot_longer(
    cols = everything(),
    names_to = "comparison",
    values_to = "lrt"
  ) %>%
  mutate(
    comparison = recode(
      comparison,
      lrt_null_vs_full = "LBT vs LID",
      lrt_null_vs_only_theta = "LBT vs DABT",
      lrt_full_vs_full_theta = "LID vs DALID",
      lrt_only_theta_vs_full_theta = "DABT vs DALID"
    ),
    category = case_when(
      is.na(lrt) ~ "Missing",
      lrt < 1e-8 ~ "< 1e-8",
      lrt >= 1e-8 & lrt <= crit_mix_005 ~ "1e-8 to 2.7055",
      lrt > crit_mix_005 ~ "> 2.7055"
    )
  ) %>%
  count(comparison, category) %>%
  pivot_wider(
    names_from = category,
    values_from = n,
    values_fill = 0
  ) 
lrt_summary
lrt_summary_gt <- lrt_summary %>%
  gt() %>%
  tab_caption(
    md("Summary of likelihood ratio statistics for the four model comparisons.")
  ) %>%
  cols_label(
    comparison = "Model comparison",
    `< 1e-8` = md("LRT $< 10^{-8}$"),
    `1e-8 to 2.7055` = md("$10^{-8} \\leq$ LRT $\\leq 2.7055$"),
    `> 2.7055` = md("LRT $> 2.7055$")
  ) %>%
  tab_options(
    table.font.size = px(12),
    column_labels.font.weight = "bold"
  )
gtsave(
  lrt_summary_gt,
  filename = "lrt_summary_table.tex"
)

# ============================================================
# Parameters against LRT
# ============================================================

make_estimate_lrt_plot <- function(data,
                                   estimate_col,
                                   lrt_col,
                                   x_lab,
                                   y_lab = "LRT statistic",
                                   title_lab = NULL,
                                   log_x = FALSE,
                                   log_y = FALSE,
                                   eps = 1e-6,
                                   highlight_fileids,
                                   highlight_values) {
  crit_mix_005 <- qchisq(0.90, df = 1)
  plot_data <- data %>%
    mutate(
      highlighted = if_else(fileid %in% highlight_fileids, fileid, "Other"),
      estimate = .data[[estimate_col]],
      lrt = .data[[lrt_col]],
      estimate_plot = if_else(rep(log_x, 410), estimate + eps, estimate),
      lrt_plot = if_else(rep(log_y, 410), lrt + eps, lrt)
    ) %>%
    filter(
      is.finite(estimate_plot),
      is.finite(lrt_plot),
      !is.na(estimate_plot),
      !is.na(lrt_plot)
    )
  
  if (log_x) {
    plot_data <- plot_data %>%
      filter(estimate_plot > 0)
  }
  
  if (log_y) {
    plot_data <- plot_data %>%
      filter(lrt_plot > 0)
  }
  
  p <- ggplot(plot_data, aes(x = estimate_plot, y = lrt_plot)) +
    geom_point(
      data = \(d) d %>% filter(highlighted == "Other"),
      colour = "grey75",
      alpha = 0.6,
      size = 2.0
    ) +
    geom_hline(
      yintercept = crit_mix_005,
      linetype = "dashed",
      linewidth = 0.9,
      colour = "black"
    ) +
    geom_point(
      data = \(d) d %>% filter(highlighted != "Other"),
      aes(colour = highlighted),
      alpha = 1,
      size = 5.2
    ) +
    scale_colour_manual(
      values = highlight_values,
      guide = "none"
    ) +
    labs(
      x = x_lab,
      y = y_lab,
      title = title_lab
    ) +
    theme_bw(base_size = 28) +
    theme(
      plot.title = element_text(face = "bold", size = 22),
      axis.title = element_text(face = "bold", size = 22),
      axis.text = element_text(size = 22, colour = "black"),
      panel.grid.minor = element_blank(),
      legend.position = "none",
      aspect.ratio = 1
    )
  
  if (log_x) {
    p <- p + scale_x_log10()
  }
  
  if (log_y) {
    p <- p + scale_y_log10()
  }
  
  p
}

p_a_lid_lrt <- make_estimate_lrt_plot(
  data = results,
  estimate_col = "full_a_est",
  lrt_col = "lrt_null_vs_full",
  x_lab = expression(hat(a)[LID] ),
  y_lab = expression(LRT[LBT~vs.~LID]),
  title_lab = "",
  log_x = TRUE,
  log_y = TRUE,
  eps = 1e-8,
  highlight_fileids = highlight_fileids,
  highlight_values = selected_colours
)
p_a_lid_lrt
p_sigma_dalid_lrt <- make_estimate_lrt_plot(
  data = results,
  estimate_col = "full_theta_sigma_theta_est",
  lrt_col = "lrt_full_vs_full_theta",
  x_lab = expression(hat(sigma)[theta*","*DALID]),
  y_lab = expression(LRT[LID~vs.~DALID]),
  title_lab = "",
  log_x = TRUE,
  log_y = TRUE,
  eps = 1e-8,
  highlight_fileids = highlight_fileids,
  highlight_values = selected_colours
)
p_sigma_dalid_lrt
p_sigma_dabt_lrt <- make_estimate_lrt_plot(
  data = results,
  estimate_col = "only_theta_sigma_theta_est",
  lrt_col = "lrt_null_vs_only_theta",
  x_lab = expression(hat(sigma)[theta*","*DABT]),
  y_lab = expression(LRT[LBT~vs.~DABT]),
  title_lab = "",
  log_x = TRUE,
  log_y = TRUE,
  eps = 1e-8,
  highlight_fileids = highlight_fileids,
  highlight_values = selected_colours
)
p_sigma_dabt_lrt
p_sigma_dalid_lrt <- make_estimate_lrt_plot(
  data = results,
  estimate_col = "full_theta_a_est",
  lrt_col = "lrt_only_theta_vs_full_theta",
  x_lab = expression(hat(a)[DALID]),
  y_lab = expression(LRT[DABT~vs.~DALID]),
  title_lab = "",
  log_x = TRUE,
  log_y = TRUE,
  eps = 1e-8,
  highlight_fileids = highlight_fileids,
  highlight_values = selected_colours
)
p_sigma_dalid_lrt


# ============================================================
# Hessian diagnostics
# ============================================================

highlighted_datasets <- setdiff(names(highlight_values), "Other")

signed_log10 <- function(x) {
  sign(x) * log10(1 + abs(x))
}

hess_diag_long <- results |>
  select(
    dataset_name,
    matches("^(null|full|only_theta|full_theta)_(max_grad|hess_)")
  ) |>
  pivot_longer(
    cols = -dataset_name,
    names_to = c("model", ".value"),
    names_pattern = "^(full_theta|only_theta|full|null)_(.*)$"
  ) |>
  mutate(
    model = recode(
      model,
      null = "LBT",
      full = "LID",
      only_theta = "DABT",
      full_theta = "DALID"
    ),
    model = factor(model, levels = c("LBT", "LID", "DABT", "DALID")),
    
    hess_min_eig_slog10 = signed_log10(hess_min_eig),
    hess_abs_cond_log10 = if_else(
      hess_abs_cond > 0 & is.finite(hess_abs_cond),
      log10(hess_abs_cond),
      NA_real_
    ),
    
    highlight_group = if_else(
      dataset_name %in% highlighted_datasets,
      dataset_name,
      "Other"
    ),
    highlight_group = factor(highlight_group, levels = names(highlight_values))
  )
glimpse(hess_diag_long)

make_hess_boxplot <- function(data, yvar, ylab, title, log_scale = FALSE) {
  
  plot_data <- data |>
    filter(
      !is.na(.data[[yvar]]),
      is.finite(.data[[yvar]])
    )
  
  p <- ggplot(plot_data, aes(x = model, y = .data[[yvar]])) +
    geom_boxplot(
      outlier.shape = NA,
      width = 0.7
    ) +
    geom_jitter(
      data = plot_data |> filter(highlight_group == "Other"),
      color = unname(highlight_values["Other"]),
      width = 0.18,
      height = 0,
      size = 1.5,
      alpha = 0.5
    ) +
    geom_jitter(
      data = plot_data |> filter(highlight_group != "Other"),
      aes(color = highlight_group),
      width = 0.18,
      height = 0,
      size = 5.2,
      alpha = 1.0
    ) +
    scale_color_manual(values = highlight_values, drop = FALSE) +
    labs(
      x = "Model",
      y = ylab,
      title = title
    ) +
    guides(color = "none") +
    theme_bw(base_size = 24) +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold"),
      axis.text.x = element_text(angle = 0, hjust = 0.5)
    )
  
  if (log_scale) {
    p <- p + scale_y_log10()
  }
  
  p
}

compress_expand_trans <- function(neg_factor = 0.25, pos_factor = 2) {
  trans_new(
    name = "compress_expand",
    transform = function(x) {
      ifelse(x < 0, x * neg_factor, x * pos_factor)
    },
    inverse = function(x) {
      ifelse(x < 0, x / neg_factor, x / pos_factor)
    }
  )
}

p_min_eig <- hess_diag_long |>
  make_hess_boxplot(
    yvar = "hess_min_eig_slog10",
    ylab = expression(sign(lambda[min]) %.% log[10](1 + abs(lambda[min]))),
    title = "",
    log_scale = FALSE
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.4
  ) +
  scale_y_continuous(
    trans = compress_expand_trans(
      neg_factor = 0.22,   # komprimer negativ del
      pos_factor = 2.0     # strekk positiv del
    ),
    breaks = c(-7, -5, -2.5, 0, 0.5, 1, 1.5),
    labels = c("-7", "-5", "-2.5", "0", "0.5", "1", "1.5")
  )
p_min_eig

loglog10_trans <- scales::trans_new(
  name = "loglog10",
  transform = function(x) log10(1 + log10(x)),
  inverse = function(y) 10^(10^y - 1)
)
p_abs_cond <- hess_diag_long |>
  filter(hess_abs_cond > 0) |>
  make_hess_boxplot(
    yvar = "hess_abs_cond",
    ylab = "Absolute Hessian condition number",
    title = "",
    log_scale = FALSE
  ) +
  scale_y_continuous(
    trans = loglog10_trans,
    breaks = c(1e0, 1e2, 1e5, 1e8, 1e12, 1e18, 1e26),
    labels = scales::scientific
  )
p_abs_cond

hess_status_pattern_table <- results |>
  transmute(
    LBT   = null_hess_status_clean,
    LID   = full_hess_status_clean,
    DABT  = only_theta_hess_status_clean,
    DALID = full_theta_hess_status_clean
  ) |>
  mutate(
    across(
      everything(),
      ~ if_else(is.na(.x), "Missing status", .x)
    )
  ) |>
  count(LBT, LID, DABT, DALID, name = "n_datasets") |>
  mutate(
    proportion = n_datasets / sum(n_datasets),
    
    n_positive_definite = rowSums(
      across(c(LBT, LID, DABT, DALID), ~ .x == "Positive definite")
    ),
    n_near_singular = rowSums(
      across(c(LBT, LID, DABT, DALID), ~ .x == "Near-singular")
    ),
    n_indefinite = rowSums(
      across(c(LBT, LID, DABT, DALID), ~ .x == "Indefinite")
    ),
    n_failed_or_nonfinite = rowSums(
      across(
        c(LBT, LID, DABT, DALID),
        ~ .x %in% c(
          "Fit failed",
          "Hessian failed",
          "Hessian non-finite",
          "Eigen failed",
          "Eigenvalue missing",
          "Hessian status missing",
          "Missing status"
        )
      )
    ),
    
    pattern_type = case_when(
      n_positive_definite == 4 ~
        "All positive definite",
      
      n_failed_or_nonfinite > 0 ~
        "At least one failed/non-finite diagnostic",
      
      n_indefinite > 0 ~
        "At least one indefinite Hessian",
      
      n_near_singular > 0 ~
        "At least one near-singular Hessian",
      
      TRUE ~
        "Other mixed pattern"
    ),
    
    pattern_type = factor(pattern_type, levels = pattern_levels)
  ) |>
  left_join(pattern_codes, by = "pattern_type") |>
  rowwise() |>
  mutate(
    problematic_models = paste(
      c("LBT", "LID", "DABT", "DALID")[
        c_across(c(LBT, LID, DABT, DALID)) %in% problem_statuses
      ],
      collapse = ", "
    ),
    problematic_models = if_else(
      problematic_models == "",
      "None",
      problematic_models
    )
  ) |>
  ungroup() |>
  arrange(
    pattern,
    desc(n_datasets),
    LBT, LID, DABT, DALID
  )

hess_status_pattern_table

hess_status_pattern_compact <- hess_status_pattern_table |>
  mutate(
    across(
      c(LBT, LID, DABT, DALID),
      ~ recode(.x, !!!status_abbrev)
    )
  ) |>
  select(
    pattern,
    problematic_models,
    LBT, LID, DABT, DALID,
    n_datasets,
    proportion
  )

hess_status_pattern_compact

highlight_hessian_patterns <- results |>
  filter(dataset_name %in% highlighted_datasets) |>
  transmute(
    dataset_name,
    LBT   = null_hess_status_clean,
    LID   = full_hess_status_clean,
    DABT  = only_theta_hess_status_clean,
    DALID = full_theta_hess_status_clean
  ) |>
  mutate(
    across(
      c(LBT, LID, DABT, DALID),
      ~ if_else(is.na(.x), "Missing status", .x)
    ),
    
    n_positive_definite = rowSums(
      across(c(LBT, LID, DABT, DALID), ~ .x == "Positive definite")
    ),
    n_near_singular = rowSums(
      across(c(LBT, LID, DABT, DALID), ~ .x == "Near-singular")
    ),
    n_indefinite = rowSums(
      across(c(LBT, LID, DABT, DALID), ~ .x == "Indefinite")
    ),
    n_failed_or_nonfinite = rowSums(
      across(
        c(LBT, LID, DABT, DALID),
        ~ .x %in% c(
          "Fit failed",
          "Hessian failed",
          "Hessian non-finite",
          "Eigen failed",
          "Eigenvalue missing",
          "Hessian status missing",
          "Missing status"
        )
      )
    ),
    
    pattern_type = case_when(
      n_positive_definite == 4 ~
        "All positive definite",
      
      n_failed_or_nonfinite > 0 ~
        "At least one failed/non-finite diagnostic",
      
      n_indefinite > 0 ~
        "At least one indefinite Hessian",
      
      n_near_singular > 0 ~
        "At least one near-singular Hessian",
      
      TRUE ~
        "Other mixed pattern"
    ),
    
    pattern_type = factor(pattern_type, levels = pattern_levels)
  ) |>
  left_join(pattern_codes, by = "pattern_type") |>
  rowwise() |>
  mutate(
    problematic_models = paste(
      c("LBT", "LID", "DABT", "DALID")[
        c_across(c(LBT, LID, DABT, DALID)) %in% problem_statuses
      ],
      collapse = ", "
    ),
    problematic_models = if_else(
      problematic_models == "",
      "None",
      problematic_models
    )
  ) |>
  ungroup() |>
  mutate(
    across(
      c(LBT, LID, DABT, DALID),
      ~ recode(.x, !!!status_abbrev)
    )
  ) |>
  select(
    dataset_name,
    pattern,
    problematic_models,
    LBT, LID, DABT, DALID
  ) |>
  arrange(pattern, dataset_name)

highlight_hessian_patterns

# ============================================================
# Log-likelihood surface visualizations
# ============================================================

get_matrix_from_fileid <- function(fileid, dom.data) {
  M <- dom.data[[fileid]]$matrix
  
  if (is.null(M)) {
    stop("Fant ikke matrix for fileid = ", fileid)
  }
  
  M
}

make_profile_obj <- function(fileid,
                             model = c("LID", "DABT"),
                             results,
                             dom.data,
                             make_pair_data,
                             f_full,
                             f_only_theta) {
  
  model <- match.arg(model)
  
  row <- results |> 
    filter(.data$fileid == !!fileid)
  
  if (nrow(row) != 1) {
    stop("Forventet n??yaktig ??n rad i results for fileid = ", fileid)
  }
  
  M <- get_matrix_from_fileid(fileid, dom.data)
  dat <- make_pair_data(M, drop_zero = TRUE)
  n_dyads <- length(dat$z)
  n <- nrow(M)
  
  if (model == "LID") {
    
    parameters <- list(
      x = matrix(0, nrow = n, ncol = 2),
      log_r = row$full_log_r_est,
      a = row$full_a_est
    )
    
    obj <- RTMB::MakeADFun(
      cmb(f_full, dat),
      parameters = parameters,
      random = "x",
      silent = TRUE
    )
    
  } else if (model == "DABT") {
    
    parameters <- list(
      x = matrix(0, nrow = n, ncol = 2),
      theta_raw = rep(0, n_dyads),
      log_r = row$only_theta_log_r_est,
      log_sigma_theta = row$only_theta_log_sigma_theta_est
    )
    
    obj <- RTMB::MakeADFun(
      cmb(f_only_theta, dat),
      parameters = parameters,
      random = c("x", "theta_raw"),
      silent = TRUE
    )
  }
  
  obj
}

profile_loglik_grid <- function(fileid,
                                model = c("LID", "DABT"),
                                results,
                                dom.data,
                                make_pair_data,
                                f_full,
                                f_only_theta,
                                n_grid = 50,
                                beta_factor = c(0.01, 640),
                                a_width = NULL,
                                sigma_factor = c(0.01, 8),
                                min_a = 0,
                                min_sigma = 1e-8,
                                min_beta = 1e-16,
                                fallback_beta = 2,
                                fallback_a = 0.4,
                                fallback_sigma = 0.5,
                                verbose = TRUE) {
  
  model <- match.arg(model)
  
  is_bad <- function(x) {
    length(x) != 1 || is.na(x) || is.nan(x) || is.infinite(x)
  }
  
  safe_positive <- function(x, fallback, lower = 1e-16) {
    if (is_bad(x) || x <= lower) {
      fallback
    } else {
      x
    }
  }
  
  safe_real <- function(x, fallback) {
    if (is_bad(x)) {
      fallback
    } else {
      x
    }
  }
  
  row <- results |> 
    dplyr::filter(.data$fileid == !!fileid)
  
  if (nrow(row) != 1) {
    stop("Forventet n??yaktig ??n rad i results for fileid = ", fileid)
  }
  
  obj <- make_profile_obj(
    fileid = fileid,
    model = model,
    results = results,
    dom.data = dom.data,
    make_pair_data = make_pair_data,
    f_full = f_full,
    f_only_theta = f_only_theta
  )
  
  if (model == "LID") {
    
    raw_log_r <- row$full_log_r_est
    raw_a     <- row$full_a_est
    
    beta_hat <- if (is_bad(raw_log_r)) {
      fallback_beta
    } else {
      exp(raw_log_r)
    }
    
    beta_hat <- safe_positive(
      beta_hat,
      fallback = fallback_beta,
      lower = min_beta
    )
    
    a_hat <- safe_real(
      raw_a,
      fallback = fallback_a
    )
    
    a_hat <- max(min_a, a_hat)
    
    used_fallback <- is_bad(raw_log_r) || is_bad(raw_a) || beta_hat <= min_beta
    
    if (verbose && used_fallback) {
      message(
        "Using fallback/moderate parameter values for fileid = ",
        fileid,
        ", model = LID. Values used: beta = ",
        signif(beta_hat, 4),
        ", a = ",
        signif(a_hat, 4)
      )
    }
    
    if (is.null(a_width)) {
      a_width <- max(0.5, 1.5 * abs(a_hat), 0.25 * beta_hat)
    }
    
    #beta_grid <- seq(
      #1e-8,
      #max(min_beta, beta_factor[2] * beta_hat),
      #length.out = n_grid
    #)
    beta_grid <- seq(
      1e-8,
      2,
      length.out = n_grid
    )
    
    a_grid <- seq(
      1e-8,
      max(min_a, a_hat + a_width),
      length.out = n_grid
    )
    
    grid <- tidyr::expand_grid(
      beta = beta_grid,
      a = a_grid
    )
    
    grid_eval <- grid |>
      dplyr::mutate(
        log_r = log(beta),
        nll = purrr::map2_dbl(log_r, a, \(lr, aa) {
          par <- obj$par
          par["log_r"] <- lr
          par["a"] <- aa
          
          val <- tryCatch(
            obj$fn(par),
            error = function(e) NA_real_
          )
          
          val
        }),
        logLik = -nll,
        rel_logLik = logLik - max(logLik, na.rm = TRUE),
        model = "LID",
        fileid = fileid,
        grid_center_beta = beta_hat,
        grid_center_a = a_hat,
        used_fallback = used_fallback
      )
    
  } else if (model == "DABT") {
    
    raw_log_r       <- row$only_theta_log_r_est
    raw_sigma_theta <- row$only_theta_sigma_theta_est
    
    beta_hat <- if (is_bad(raw_log_r)) {
      fallback_beta
    } else {
      exp(raw_log_r)
    }
    
    beta_hat <- safe_positive(
      beta_hat,
      fallback = fallback_beta,
      lower = min_beta
    )
    
    sigma_hat <- safe_positive(
      raw_sigma_theta,
      fallback = fallback_sigma,
      lower = min_sigma
    )
    
    used_fallback <- is_bad(raw_log_r) || 
      is_bad(raw_sigma_theta) || 
      beta_hat <= min_beta || 
      sigma_hat <= min_sigma
    
    if (verbose && used_fallback) {
      message(
        "Using fallback/moderate parameter values for fileid = ",
        fileid,
        ", model = DABT. Values used: beta = ",
        signif(beta_hat, 4),
        ", sigma_theta = ",
        signif(sigma_hat, 4)
      )
    }
    
    beta_grid <- seq(
      max(min_beta, beta_factor[1] * beta_hat),
      max(min_beta, beta_factor[2] * beta_hat),
      length.out = n_grid
    )
    
    sigma_grid <- exp(seq(
      log(max(min_sigma, sigma_factor[1] * sigma_hat)),
      log(max(min_sigma, sigma_factor[2] * sigma_hat)),
      length.out = n_grid
    ))
    
    grid <- tidyr::expand_grid(
      beta = beta_grid,
      sigma_theta = sigma_grid
    )
    
    grid_eval <- grid |>
      dplyr::mutate(
        log_r = log(beta),
        log_sigma_theta = log(sigma_theta),
        nll = purrr::map2_dbl(log_r, log_sigma_theta, \(lr, lsig) {
          par <- obj$par
          par["log_r"] <- lr
          par["log_sigma_theta"] <- lsig
          
          val <- tryCatch(
            obj$fn(par),
            error = function(e) NA_real_
          )
          
          val
        }),
        logLik = -nll,
        rel_logLik = logLik - max(logLik, na.rm = TRUE),
        model = "DABT",
        fileid = fileid,
        grid_center_beta = beta_hat,
        grid_center_sigma_theta = sigma_hat,
        used_fallback = used_fallback
      )
  }
  
  grid_eval
}

plot_profile_contour_2d <- function(grid_eval,
                                    use_relative = FALSE,
                                    bins = 80,
                                    legend_n = 6,
                                    legend_breaks = NULL) {
  
  model <- unique(grid_eval$model)
  
  if (length(model) != 1) {
    stop("grid_eval m?? komme fra ??n modell.")
  }
  
  # Remove invalid log-likelihood values
  grid_eval <- grid_eval |>
    dplyr::filter(is.finite(logLik))
  
  if (nrow(grid_eval) == 0) {
    stop("Ingen endelige logLik-verdier ?? plotte.")
  }
  
  # Recompute relative log-likelihood using only finite values
  grid_eval <- grid_eval |>
    dplyr::mutate(
      rel_logLik = logLik - max(logLik, na.rm = TRUE)
    )
  
  z_var <- if (use_relative) "rel_logLik" else "logLik"
  
  z_range <- range(grid_eval[[z_var]], na.rm = TRUE)
  
  if (!all(is.finite(z_range))) {
    stop("z_range er ikke endelig: ", paste(z_range, collapse = ", "))
  }
  
  if (diff(z_range) == 0) {
    stop("Alle gyldige ", z_var, "-verdier er like. Kan ikke lage contour-plot.")
  }
  
  contour_breaks <- seq(
    z_range[1],
    z_range[2],
    length.out = bins + 1
  )
  
  if (is.null(legend_breaks)) {
    legend_breaks <- pretty(z_range, n = legend_n)
    legend_breaks <- legend_breaks[
      legend_breaks >= z_range[1] & legend_breaks <= z_range[2]
    ]
  }
  
  fill_lab <- if (use_relative) {
    expression(ell - max(ell))
  } else {
    "Log-likelihood"
  }
  
  fill_scale <- scale_fill_gradientn(
    colours = hcl.colors(256, "viridis"),
    breaks = legend_breaks,
    limits = z_range,
    name = fill_lab,
    guide = guide_colorbar(
      ticks = TRUE,
      barheight = unit(5, "cm")
    )
  )
  
  mle_point <- grid_eval |>
    dplyr::slice_max(logLik, n = 1, with_ties = FALSE)
  
  if (model == "LID") {
    
    p <- ggplot(grid_eval, aes(x = a, y = beta, z = .data[[z_var]])) +
      geom_contour_filled(
        aes(fill = after_stat(level_mid)),
        breaks = contour_breaks
      ) +
      geom_point(
        data = mle_point,
        aes(x = a, y = beta),
        inherit.aes = FALSE,
        size = 3.6
      ) +
      fill_scale +
      labs(
        x = expression(a),
        y = expression(beta)
      ) +
      theme_bw(base_size = 22)
    
  } else if (model == "DABT") {
    
    p <- ggplot(grid_eval, aes(x = sigma_theta, y = beta, z = .data[[z_var]])) +
      geom_contour_filled(
        aes(fill = after_stat(level_mid)),
        breaks = contour_breaks
      ) +
      geom_point(
        data = mle_point,
        aes(x = sigma_theta, y = beta),
        inherit.aes = FALSE,
        size = 3.6
      ) +
      scale_y_log10() +
      fill_scale +
      labs(
        x = expression(sigma[theta]),
        y = expression(beta)
      ) +
      theme_bw(base_size = 22)
  }
  
  p
}
#dom.data_filtered$ScottLockhard_1999b$matrix

plot_profile_contour_3d <- function(grid_eval,
                                   use_relative = TRUE) {
  
  z_var <- if (use_relative) "rel_logLik" else "logLik"
  model <- unique(grid_eval$model)
  
  if (length(model) != 1) {
    stop("grid_eval m?? komme fra ??n modell.")
  }
  
  if (model == "LID") {
    
    zmat <- grid_eval |>
      select(beta, a, z = all_of(z_var)) |>
      tidyr::pivot_wider(names_from = a, values_from = z) |>
      arrange(beta)
    
    x <- zmat$beta
    y <- as.numeric(names(zmat)[-1])
    z <- as.matrix(zmat[, -1])
    
    fig <- plot_ly(
      x = x,
      y = y,
      z = z,
      type = "surface",
      contours = list(
        z = list(
          show = TRUE,
          usecolormap = TRUE,
          highlightcolor = "#ff0000",
          project = list(z = TRUE)
        )
      )
    ) |>
      layout(
        title = paste("LID profile log-likelihood:", unique(grid_eval$fileid)),
        scene = list(
          xaxis = list(title = "beta"),
          yaxis = list(title = "a"),
          zaxis = list(title = if (use_relative) "relative logLik" else "logLik")
        )
      )
    
  } else if (model == "DABT") {
    
    zmat <- grid_eval |>
      mutate(log_sigma_theta_plot = log(sigma_theta)) |>
      select(beta, log_sigma_theta_plot, z = all_of(z_var)) |>
      tidyr::pivot_wider(names_from = log_sigma_theta_plot, values_from = z) |>
      arrange(beta)
    
    x <- zmat$beta
    y <- as.numeric(names(zmat)[-1])
    z <- as.matrix(zmat[, -1])
    
    fig <- plot_ly(
      x = x,
      y = y,
      z = z,
      type = "surface",
      contours = list(
        z = list(
          show = TRUE,
          usecolormap = TRUE,
          highlightcolor = "#ff0000",
          project = list(z = TRUE)
        )
      )
    ) |>
      layout(
        title = paste("DABT profile log-likelihood:", unique(grid_eval$fileid)),
        scene = list(
          xaxis = list(title = "beta"),
          yaxis = list(title = "log(sigma_theta)"),
          zaxis = list(title = if (use_relative) "relative logLik" else "logLik")
        )
      )
  }
  
  fig
}

fileid_i <- "Poisbleau_2005c"

# LID: loglikelihood som funksjon av beta og a
grid_lid <- profile_loglik_grid(
  fileid = fileid_i,
  model = "LID",
  results = results,
  dom.data = dom.data_filtered,
  make_pair_data = make_pair_data,
  f_full = f_full,
  f_only_theta = f_only_theta,
  n_grid = 150
)

p_lid_2d <- plot_profile_contour_2d(grid_lid)
p_lid_2d