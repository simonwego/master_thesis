library(DomArchive)
library(dplyr)
library(purrr)
library(tidyr)
library(ggplot2)
library(plotly)
library(RTMB)

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