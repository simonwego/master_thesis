library(DomArchive)
library(RTMB)
library(dplyr)
results <- readRDS("~/DomArchiveResults/four_models_comparison_filtered_betterdiag_results.rds")

invlogit_ad <- function(eta) {
  1 / (1 + exp(-eta))
}
dom.data$Bonanni_2017c

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
    prob = invlogit_ad(eta),
    log = TRUE
  ))
  
  nll
}
cmb <- function(f, d) function(p) f(p, d)

## =========================================================
## Tilpasning
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

## =========================================================
## Hjelpefunksjon: matrix -> dyad data
## =========================================================

matrix_to_dyad_data <- function(M) {
  n <- nrow(M)
  
  ij <- which(upper.tri(M), arr.ind = TRUE)
  i <- ij[, 1]
  j <- ij[, 2]
  
  s <- M[cbind(i, j)] + M[cbind(j, i)]
  keep <- is.finite(s) & s > 0
  
  list(
    n = n,
    i = i[keep],
    j = j[keep],
    z = M[cbind(i[keep], j[keep])],
    s = s[keep]
  )
}


## =========================================================
## Fit valgt modell
## =========================================================

fit_one_model <- function(M,
                          model = c("LBT", "LID", "DABT", "DALID"),
                          init_log_r = 0,
                          init_a = 0.1,
                          init_log_sigma_theta = 0.1) {
  
  model <- match.arg(model)
  data <- matrix_to_dyad_data(M)
  n <- data$n
  m <- length(data$z)
  
  if (model == "LBT") {
    parms <- list(
      x = matrix(0, n, 2),
      log_r = init_log_r
    )
    
    fit <- fit_model(
      data = data,
      func = f_null,
      parms = parms,
      random = "x"
    )
  }
  
  if (model == "LID") {
    parms <- list(
      x = matrix(0, n, 2),
      log_r = init_log_r,
      a = init_a
    )
    
    fit <- fit_model(
      data = data,
      func = f_full,
      parms = parms,
      random = "x",
      lower = c(log_r = -Inf, a = 0),
      upper = c(log_r =  Inf, a = Inf)
    )
  }
  
  if (model == "DABT") {
    parms <- list(
      x = matrix(0, n, 2),
      theta_raw = rep(0, m),
      log_r = init_log_r,
      log_sigma_theta = init_log_sigma_theta
    )
    
    fit <- fit_model(
      data = data,
      func = f_only_theta,
      parms = parms,
      random = c("x", "theta_raw")
    )
  }
  
  if (model == "DALID") {
    parms <- list(
      x = matrix(0, n, 2),
      theta_raw = rep(0, m),
      log_r = init_log_r,
      a = init_a,
      log_sigma_theta = init_log_sigma_theta
    )
    
    fit <- fit_model(
      data = data,
      func = f_full_theta,
      parms = parms,
      random = c("x", "theta_raw"),
      lower = c(log_r = -Inf, a = 0, log_sigma_theta = -Inf),
      upper = c(log_r =  Inf, a = Inf, log_sigma_theta =  Inf)
    )
  }
  
  fit$model <- model
  fit
}


## =========================================================
## Simuler datasett under en fitted nullmodell
## =========================================================

simulate_from_fitted_model <- function(M_template,
                                       fit_null,
                                       null_model = c("LBT", "LID", "DABT", "DALID"),
                                       d = 2) {
  
  null_model <- match.arg(null_model)
  
  if (is.null(fit_null$opt)) {
    stop("fit_null$opt is NULL, cannot simulate from fitted null model.")
  }
  
  par <- fit_null$opt$par
  
  n <- nrow(M_template)
  S <- M_template + t(M_template)
  
  X <- matrix(rnorm(n * d), nrow = n, ncol = d)
  
  M_sim <- matrix(0, nrow = n, ncol = n)
  
  r_hat <- exp(par[["log_r"]])
  
  a_hat <- if ("a" %in% names(par)) {
    par[["a"]]
  } else {
    0
  }
  
  sigma_theta_hat <- if ("log_sigma_theta" %in% names(par)) {
    exp(par[["log_sigma_theta"]])
  } else {
    0
  }
  
  for (i in seq_len(n - 1)) {
    for (j in (i + 1):n) {
      s_ij <- S[i, j]
      
      if (is.finite(s_ij) && s_ij > 0) {
        
        bt_part <- r_hat * (X[i, 1] - X[j, 1])
        
        lid_part <- a_hat * (
          X[i, 2] * X[j, 1] -
            X[i, 1] * X[j, 2]
        )
        
        theta_part <- if (sigma_theta_hat > 0) {
          sigma_theta_hat * rnorm(1)
        } else {
          0
        }
        
        eta_ij <- bt_part + lid_part + theta_part
        p_ij <- plogis(eta_ij)
        
        m_ij <- rbinom(1, size = s_ij, prob = p_ij)
        
        M_sim[i, j] <- m_ij
        M_sim[j, i] <- s_ij - m_ij
      }
    }
  }
  
  diag(M_sim) <- 0
  M_sim
}


## =========================================================
## Bootstrap LRT for modell-sammenligning
## =========================================================

bootstrap_lrt_one_comparison <- function(M,
                                         null_model,
                                         full_model,
                                         B = 500,
                                         d = 2,
                                         seed = NULL,
                                         tol = 1e-7,
                                         verbose = TRUE) {
  
  if (!is.null(seed)) set.seed(seed)
  
  fit_null_obs <- fit_one_model(M, null_model)
  fit_full_obs <- fit_one_model(M, full_model)
  
  lrt_obs <- if (
    isTRUE(fit_null_obs$converged) &&
    isTRUE(fit_full_obs$converged) &&
    is.finite(fit_null_obs$logLik) &&
    is.finite(fit_full_obs$logLik)
  ) {
    2 * (fit_full_obs$logLik - fit_null_obs$logLik)
  } else {
    NA_real_
  }
  
  if (is.finite(lrt_obs) && lrt_obs < tol) {
    lrt_obs <- 0
  }
  
  lrt_boot <- rep(NA_real_, B)
  conv_null <- rep(FALSE, B)
  conv_full <- rep(FALSE, B)
  
  for (b in seq_len(B)) {
    if (verbose && b %% 25 == 0) {
      message(null_model, " vs ", full_model, ": bootstrap ", b, " / ", B)
    }
    
    M_b <- try(
      simulate_from_fitted_model(
        M_template = M,
        fit_null = fit_null_obs,
        null_model = null_model,
        d = d
      ),
      silent = TRUE
    )
    
    if (inherits(M_b, "try-error")) {
      next
    }
    
    fit_null_b <- try(fit_one_model(M_b, null_model), silent = TRUE)
    fit_full_b <- try(fit_one_model(M_b, full_model), silent = TRUE)
    
    if (
      !inherits(fit_null_b, "try-error") &&
      !inherits(fit_full_b, "try-error")
    ) {
      conv_null[b] <- isTRUE(fit_null_b$converged)
      conv_full[b] <- isTRUE(fit_full_b$converged)
      
      if (
        conv_null[b] &&
        conv_full[b] &&
        is.finite(fit_null_b$logLik) &&
        is.finite(fit_full_b$logLik)
      ) {
        lrt_b <- 2 * (fit_full_b$logLik - fit_null_b$logLik)
        
        if (is.finite(lrt_b) && lrt_b < tol) {
          lrt_b <- 0
        }
        
        lrt_boot[b] <- max(lrt_b, 0)
      }
    }
  }
  
  lrt_valid <- lrt_boot[is.finite(lrt_boot)]
  
  p_boot <- if (length(lrt_valid) > 0 && is.finite(lrt_obs)) {
    (1 + sum(lrt_valid >= lrt_obs)) / (length(lrt_valid) + 1)
  } else {
    NA_real_
  }
  
  list(
    comparison = paste(null_model, "vs", full_model),
    null_model = null_model,
    full_model = full_model,
    lrt_obs = lrt_obs,
    lrt_boot = lrt_boot,
    p_boot = p_boot,
    n_success = length(lrt_valid),
    success_rate = length(lrt_valid) / B,
    conv_null_rate = mean(conv_null),
    conv_full_rate = mean(conv_full),
    fit_null_obs = fit_null_obs,
    fit_full_obs = fit_full_obs
  )
}


## =========================================================
## Kjør alle fire sammenligningene
## =========================================================

bootstrap_four_comparisons <- function(M,
                                       B = 500,
                                       d = 2,
                                       seed = NULL,
                                       verbose = TRUE) {
  
  comparisons <- list(
    c("LBT",  "LID"),
    c("LBT",  "DABT"),
    c("LID",  "DALID"),
    c("DABT", "DALID")
  )
  
  out <- vector("list", length(comparisons))
  
  for (k in seq_along(comparisons)) {
    null_model <- comparisons[[k]][1]
    full_model <- comparisons[[k]][2]
    
    out[[k]] <- bootstrap_lrt_one_comparison(
      M = M,
      null_model = null_model,
      full_model = full_model,
      B = B,
      d = d,
      seed = if (!is.null(seed)) seed + k - 1 else NULL,
      verbose = verbose
    )
  }
  
  names(out) <- vapply(
    comparisons,
    function(x) paste(x[1], "vs", x[2]),
    character(1)
  )
  
  out
}


## =========================================================
## Kompakt oppsummeringstabell
## =========================================================

summarise_bootstrap_results <- function(boot_results) {
  do.call(
    rbind,
    lapply(boot_results, function(res) {
      data.frame(
        comparison = res$comparison,
        null_model = res$null_model,
        full_model = res$full_model,
        lrt_obs = res$lrt_obs,
        p_boot = res$p_boot,
        n_success = res$n_success,
        success_rate = res$success_rate,
        conv_null_rate = res$conv_null_rate,
        conv_full_rate = res$conv_full_rate
      )
    })
  )
}

fileid <- "Bonanni_2017c"
M <- dom.data[[fileid]]$matrix

fit_LBT_obs <- fit_one_model(M, "LBT")
fit_LID_obs <- fit_one_model(M, "LID")

c(
  logLik_LBT = fit_LBT_obs$logLik,
  logLik_LID = fit_LID_obs$logLik,
  lrt = 2 * (fit_LID_obs$logLik - fit_LBT_obs$logLik),
  LBT_converged = fit_LBT_obs$converged,
  LID_converged = fit_LID_obs$converged
)
fit_LID_obs$opt$par
results %>%
  filter(fileid == "Bonanni_2017c") %>%
  select(fileid, null_logLik, full_logLik, lrt_null_vs_full, full_a_est)

dataset_names <- c("Alados_1992b, Adcock_2015a, Bonanni_2017c") #High lrt, low lrt, border lrt

fileid <- "Bonanni_2017c"
fileid
M <- dom.data[[fileid]]$matrix
M
a <- results%>%
  filter(fileid =="Bonanni_2017c")
a$lrt_null_vs_full
boot_res <- bootstrap_four_comparisons(
  M = M,
  B = 1000,
  seed = 123,
  verbose = TRUE
)

boot_summary <- summarise_bootstrap_results(boot_res)
print(boot_summary)

bootstrap_dir <- file.path(
  "DomArchiveResults",
  "bootstrap_lrt"
)
dir.create(bootstrap_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(
  boot_res,
  file = file.path(
    bootstrap_dir,
    paste0("bootstrap_lrt_", fileid, "_s123.rds")
  )
)

saveRDS(
  boot_summary,
  file = file.path(
    bootstrap_dir,
    paste0("bootstrap_lrt_summary_", fileid, "_s123.rds")
  )
)

dir.create(bootstrap_dir, recursive = TRUE, showWarnings = FALSE)

all_boot <- lapply(dataset_names, function(id) {
  message("Running bootstrap for ", id)
  
  M <- dom.data[[id]]$matrix
  
  res <- try(
    bootstrap_four_comparisons(
      M = M,
      B = 500,
      seed = 123,
      verbose = FALSE
    ),
    silent = TRUE
  )
  
  if (inherits(res, "try-error")) {
    return(NULL)
  }
  
  summary <- summarise_bootstrap_results(res)
  summary$fileid <- id
  summary
})

all_boot_summary <- do.call(rbind, all_boot)
all_boot_summary

## =========================================================
# no refit, Parameter bootstrap 
## =========================================================

get_col_value <- function(row, possible_names, default = NA_real_) {
  nm <- intersect(possible_names, names(row))
  
  if (length(nm) == 0) {
    return(default)
  }
  
  val <- row[[nm[1]]]
  
  if (length(val) == 0 || is.null(val)) {
    return(default)
  }
  
  as.numeric(val[1])
}

get_fitted_params_from_results <- function(results_row,
                                           model = c("LBT", "LID", "DABT", "DALID")) {
  
  model <- match.arg(model)
  
  if (model == "LBT") {
    return(list(
      model = "LBT",
      log_r = get_col_value(
        results_row,
        c("null_log_r_est", "log_r_null", "null_log_r", "LBT_log_r")
      ),
      a = 0,
      log_sigma_theta = NA_real_
    ))
  }
  
  if (model == "LID") {
    return(list(
      model = "LID",
      log_r = get_col_value(
        results_row,
        c("full_log_r_est", "log_r_full", "full_log_r", "LID_log_r")
      ),
      a = get_col_value(
        results_row,
        c("full_a_est", "a_full", "full_a", "LID_a")
      ),
      log_sigma_theta = NA_real_
    ))
  }
  
  if (model == "DABT") {
    sigma_theta <- get_col_value(
      results_row,
      c("only_theta_sigma_theta_est")
    )
    
    log_sigma_theta <- get_col_value(
      results_row,
      c("only_theta_log_sigma_theta_est")
    )
    
    if (!is.finite(sigma_theta) && is.finite(log_sigma_theta)) {
      sigma_theta <- exp(log_sigma_theta)
    }
    
    return(list(
      model = "DABT",
      log_r = get_col_value(
        results_row,
        c("only_theta_log_r_est")
      ),
      a = 0,
      log_sigma_theta = log_sigma_theta
    ))
  }
  
  if (model == "DALID") {
    sigma_theta <- get_col_value(
      results_row,
      c("full_theta_sigma_theta_est")
    )
    
    log_sigma_theta <- get_col_value(
      results_row,
      c("full_theta_log_sigma_theta_est")
    )
    
    if (!is.finite(sigma_theta) && is.finite(log_sigma_theta)) {
      sigma_theta <- exp(log_sigma_theta)
    }
    
    return(list(
      model = "DALID",
      log_r = get_col_value(
        results_row,
        c("full_theta_log_r_est", "log_r_full_theta",
          "full_theta_log_r", "DALID_log_r")
      ),
      a = get_col_value(
        results_row,
        c("full_theta_a_est", "a_full_theta",
          "full_theta_a", "DALID_a")
      ),
      log_sigma_theta = log_sigma_theta
    ))
  }
}
get_fitted_params_from_fit <- function(fit, model) {
  
  par <- fit$opt$par
  
  tibble::tibble(
    model = model,
    log_r = unname(par["log_r"]),
    a = unname(par["a"]),
    log_sigma_theta = unname(par["log_sigma_theta"]))
}
simulate_from_results_model <- function(M_template,
                                        results_row,
                                        null_model = c("LBT", "LID", "DABT", "DALID"),
                                        d = 2) {
  
  null_model <- match.arg(null_model)
  pars <- get_fitted_params_from_results(results_row, null_model)
  
  if (!is.finite(pars$log_r)) {
    stop("Missing or non-finite log_r estimate for model ", null_model)
  }
  
  if (!is.finite(pars$a)) {
    pars$a <- 0
  }
  
  if (!is.finite(pars$log_sigma_theta)) {
    pars$log_sigma_theta <- 0
  }
  
  n <- nrow(M_template)
  S <- M_template + t(M_template)
  
  X <- matrix(rnorm(n * d), nrow = n, ncol = d)
  
  r_hat <- exp(pars$log_r)
  a_hat <- pars$a
  sigma_theta_hat <- exp(pars$log_sigma_theta)
  
  M_sim <- matrix(0, nrow = n, ncol = n)
  
  for (i in seq_len(n - 1)) {
    for (j in (i + 1):n) {
      s_ij <- S[i, j]
      
      if (is.finite(s_ij) && s_ij > 0) {
        
        bt_part <- r_hat * (X[i, 1] - X[j, 1])
        
        lid_part <- a_hat * (
          X[i, 2] * X[j, 1] -
            X[i, 1] * X[j, 2]
        )
        
        theta_part <- if (sigma_theta_hat > 0) {
          sigma_theta_hat * rnorm(1)
        } else {
          0
        }
        
        eta_ij <- bt_part + lid_part + theta_part
        p_ij <- plogis(eta_ij)
        
        m_ij <- rbinom(1, size = s_ij, prob = p_ij)
        
        M_sim[i, j] <- m_ij
        M_sim[j, i] <- s_ij - m_ij
      }
    }
  }
  
  diag(M_sim) <- 0
  M_sim
}

bootstrap_parameter_sample <- function(results, dom.data, dataset_name, model_name, B = 1000) {
  
  row <- results %>%
    filter(.data$fileid == .env$dataset_name)
  row
  if (nrow(row) != 1) {
    stop("Expected exactly one row in results for dataset_name = ", dataset_name,
         ", but found ", nrow(row), ".")
  }
  
  if (!dataset_name %in% names(dom.data)) {
    stop("dataset_name not found in dom.data: ", dataset_name)
  }
  
  params <- get_fitted_params_from_results(row, model_name)
  bootstrap_params <- vector("list", B)
  
  for (b in seq_len(B)) {
    
    bootstrap_params[[b]] <- tryCatch({
      
      M_b <- simulate_from_results_model(
        dom.data[[dataset_name]]$matrix,
        row,
        model_name
      )
      
      fit_b <- fit_one_model(
        M_b,
        params$model,
        init_log_r = -1,
        init_a = 0.185,
        init_log_sigma_theta = params$log_sigma_theta
      )
      
      get_fitted_params_from_fit(fit_b, model_name) %>%
        mutate(
          bootstrap_id = b,
          status = ifelse(fit_b$converged, "ok", "not_converged"),
          convergence = if (!is.null(fit_b$opt)) fit_b$opt$convergence else NA_integer_,
          logLik = fit_b$logLik
        )
      
    }, error = function(e) {
      
      tibble::tibble(
        bootstrap_id = b,
        model = model_name,
        log_r = NA_real_,
        a = NA_real_,
        log_sigma_theta = NA_real_,
        status = "failed",
        error_message = conditionMessage(e)
      )
    })
  }
  
  bind_rows(bootstrap_params)
}

highlighted_datasets <- c(
  "Alados_1992b",
  "Adcock_2015a",
  "Blatrix_2004c",
  "Correa_2013a",
  "Cote_2000d",
  "Cui_2014",
  "Kolodziejczyk_2005",
  "Mwamende_2009a",
  "Poisbleau_2005c",
  "Prieto_1978",
  "ScottLockhard_1999b",
  "Shimoji_2014c"
)
models <- c("LBT", "LID", "DABT", "DALID")
B <- 1000

all_bootstrap_params <- tidyr::crossing(
  dataset_name = highlighted_datasets,
  model_name = models
) %>%
  mutate(
    boot = purrr::map2(
      dataset_name,
      model_name,
      ~ {
        message("Running bootstrap: ", .x, " | ", .y)
        
        bootstrap_parameter_sample(
          results = results,
          dom.data = dom.data,
          dataset_name = .x,
          model_name = .y,
          B = B
        )
      }
    )
  ) %>%
  tidyr::unnest(boot) %>%
  select(-model)
all_bootstrap_params %>%
  filter(dataset_name == "Kolodziejczyk_2005", model_name == "DALID") %>%
  count(status)
B_used_table <- all_bootstrap_params %>%
  filter(status == "ok") %>%
  group_by(dataset_name, model_name) %>%
  summarise(
    B_used = n_distinct(bootstrap_id),
    .groups = "drop"
  )

bootstrap_param_summary_wide <- all_bootstrap_params %>%
  filter(status == "ok") %>%
  select(
    dataset_name,
    model_name,
    bootstrap_id,
    any_of(c("log_r", "a", "log_sigma_theta"))
  ) %>%
  pivot_longer(
    cols = any_of(c("log_r", "a", "log_sigma_theta")),
    names_to = "parameter",
    values_to = "estimate"
  ) %>%
  filter(is.finite(estimate)) %>%
  group_by(dataset_name, model_name, parameter) %>%
  summarise(
    bootstrap_mean = mean(estimate),
    bootstrap_sd = sd(estimate),
    .groups = "drop"
  ) %>%
  mutate(
    parameter = factor(
      parameter,
      levels = c("log_r", "a", "log_sigma_theta")
    ),
    value = sprintf("%.3f (%.3f)", bootstrap_mean, bootstrap_sd)
  ) %>%
  select(
    dataset_name,
    model_name,
    parameter,
    value
  ) %>%
  pivot_wider(
    names_from = parameter,
    values_from = value
  ) %>%
  left_join(
    B_used_table,
    by = c("dataset_name", "model_name")
  ) %>%
  mutate(
    model_name = factor(
      model_name,
      levels = c("LBT", "LID", "DABT", "DALID")
    )
  ) %>%
  arrange(dataset_name, model_name) %>%
  select(
    dataset_name,
    model_name,
    B_used,
    log_r,
    a,
    log_sigma_theta
  )
bootstrap_param_summary_wide <- bootstrap_param_summary_wide %>%
  mutate(
    across(
      c(log_r, log_sigma_theta, a),
      ~ replace_na(.x, "--")
    )
  )
bootstrap_param_summary_wide_latex <- bootstrap_param_summary_wide %>%
  rename(
    Dataset = dataset_name,
    Model = model_name,
    `$B$` = B_used,
    `$\\log r$` = log_r,
    `$\\log \\sigma_\\theta$` = log_sigma_theta,
    `$a$` = a
  )
latex_boot_table <- bootstrap_param_summary_wide_latex %>%
  knitr::kable(
    format = "latex",
    booktabs = TRUE,
    escape = FALSE,
    caption = "Bootstrap means with bootstrap standard deviations in parentheses.",
    label = "tab:bootstrap_parameter_summary"
  ) %>%
  kableExtra::kable_styling(
    latex_options = c("hold_position", "scale_down")
  )
latex_boot_table
writeLines(latex_boot_table, "Tables/bootstrap_parameter_summary.tex")
bootstrap_param_summary_long <- all_bootstrap_params %>%
  filter(status == "ok") %>%
  select(
    dataset_name,
    model_name,
    bootstrap_id,
    any_of(c("log_r", "r", "a", "log_sigma_theta", "sigma_theta"))
  ) %>%
  pivot_longer(
    cols = any_of(c("log_r", "r", "a", "log_sigma_theta", "sigma_theta")),
    names_to = "parameter",
    values_to = "estimate"
  ) %>%
  filter(is.finite(estimate)) %>%
  group_by(dataset_name, model_name, parameter) %>%
  summarise(
    B_ok = n(),
    bootstrap_mean = mean(estimate),
    bootstrap_sd = sd(estimate),
    
    mcse_mean = bootstrap_sd / sqrt(B_ok),
    relative_mcse_mean = mcse_mean / bootstrap_sd,
    
    mcse_sd_approx = bootstrap_sd / sqrt(2 * (B_ok - 1)),
    relative_mcse_sd_approx = mcse_sd_approx / bootstrap_sd,
    
    q025 = quantile(estimate, 0.025),
    q975 = quantile(estimate, 0.975),
    .groups = "drop"
  )
bootstrap_param_summary_long
bootstrap_param_summary_long %>%
  filter(dataset_name == "ScottLockhard_1999b")
bootstrap_table <- bootstrap_param_summary_long %>%
  mutate(
    model_name = factor(model_name, levels = c("LBT", "LID", "DABT", "DALID")),
    parameter = factor(parameter, levels = c("log_r", "a", "log_sigma_theta")),
    value = sprintf("%.3f (%.3f)", bootstrap_mean, bootstrap_sd)) %>%
  select(dataset_name, model_name, parameter, value) %>%
  pivot_wider(names_from = parameter,
              values_from = value) %>%
  select(dataset_name, model_name, log_r, a, log_sigma_theta) %>%
  arrange(dataset_name, model_name)%>%
  mutate(
      across(c(log_r, a, log_sigma_theta), ~replace_na(.x, "--"))
    ) %>%
  rename(
    Dataset = dataset_name,
    Model = model_name,
    '$\\log r$' = log_r,
    '$\\log \\sigma_\\theta$' = log_sigma_theta,
    '$a$' = a
    )
latex_boot_table <- bootstrap_table %>%
  knitr::kable(
    format = "latex",
    booktabs = TRUE,
    escape = FALSE,
    caption = "Bootstrap means with bootstrap standard deviations in parentheses.",
    label = "tab:bootstrap_parameter_summary"
  ) %>%
  kableExtra::kable_styling(
    latex_options = c("hold_position", "scale_down")
  )
latex_boot_table
all_bootstrap_params %>%
  count(dataset_name, model_name, status) %>%
  group_by(dataset_name, model_name) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()
bootstrap_dir <- file.path(
  "DomArchiveResults",
  "bootstrap_parameters"
)
all_bootstrap_params
dir.create(bootstrap_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(
  all_bootstrap_params,
  file = file.path(
    bootstrap_dir,
    paste0("all_bootstrap_param_samples.rds")
  )
)

##############################
# LRT
##############################
get_lrt_obs_from_results <- function(results_row,
                                     null_model,
                                     full_model) {
  
  comparison <- paste(null_model, full_model, sep = "_vs_")
  
  if (comparison == "LBT_vs_LID") {
    return(get_col_value(
      results_row,
      c("lrt_null_vs_full", "lrt_LBT_vs_LID")
    ))
  }
  
  if (comparison == "LBT_vs_DABT") {
    return(get_col_value(
      results_row,
      c("lrt_null_vs_theta", "lrt_null_vs_only_theta", "lrt_LBT_vs_DABT")
    ))
  }
  
  if (comparison == "LID_vs_DALID") {
    return(get_col_value(
      results_row,
      c("lrt_full_vs_full_theta", "lrt_LID_vs_DALID")
    ))
  }
  
  if (comparison == "DABT_vs_DALID") {
    return(get_col_value(
      results_row,
      c("lrt_theta_vs_full_theta", "lrt_only_theta_vs_full_theta",
        "lrt_DABT_vs_DALID")
    ))
  }
  
  NA_real_
}

bootstrap_lrt_one_comparison_from_results <- function(M,
                                                      results_row,
                                                      null_model,
                                                      full_model,
                                                      B = 500,
                                                      d = 2,
                                                      seed = NULL,
                                                      tol = 1e-8,
                                                      verbose = TRUE) {
  
  if (!is.null(seed)) set.seed(seed)
  
  lrt_obs <- get_lrt_obs_from_results(
    results_row = results_row,
    null_model = null_model,
    full_model = full_model
  )
  
  if (!is.finite(lrt_obs)) {
    warning(
      "Could not find observed LRT in results for ",
      null_model, " vs ", full_model,
      ". Refitting observed data instead."
    )
    
    fit_null_obs <- fit_one_model(M, null_model)
    fit_full_obs <- fit_one_model(M, full_model)
    
    if (
      is.finite(fit_null_obs$logLik) &&
      is.finite(fit_full_obs$logLik)
    ) {
      lrt_obs <- 2 * (fit_full_obs$logLik - fit_null_obs$logLik)
    } else {
      lrt_obs <- NA_real_
    }
  }
  
  if (is.finite(lrt_obs) && lrt_obs < tol) {
    lrt_obs <- 0
  }
  
  lrt_boot <- rep(NA_real_, B)
  
  conv_null <- rep(FALSE, B)
  conv_full <- rep(FALSE, B)
  
  finite_null_logLik <- rep(FALSE, B)
  finite_full_logLik <- rep(FALSE, B)
  finite_lrt <- rep(FALSE, B)
  
  for (b in seq_len(B)) {
    if (verbose && b %% 25 == 0) {
      message(null_model, " vs ", full_model, ": bootstrap ", b, " / ", B)
    }
    
    M_b <- try(
      simulate_from_results_model(
        M_template = M,
        results_row = results_row,
        null_model = null_model,
        d = d
      ),
      silent = TRUE
    )
    
    if (inherits(M_b, "try-error")) {
      next
    }
    
    fit_null_b <- try(fit_one_model(M_b, null_model), silent = TRUE)
    fit_full_b <- try(fit_one_model(M_b, full_model), silent = TRUE)
    
    if (
      !inherits(fit_null_b, "try-error") &&
      !inherits(fit_full_b, "try-error")
    ) {
      conv_null[b] <- isTRUE(fit_null_b$converged)
      conv_full[b] <- isTRUE(fit_full_b$converged)
      
      finite_null_logLik[b] <- is.finite(fit_null_b$logLik)
      finite_full_logLik[b] <- is.finite(fit_full_b$logLik)
      
      if (
        finite_null_logLik[b] &&
        finite_full_logLik[b]
      ) {
        lrt_b <- 2 * (fit_full_b$logLik - fit_null_b$logLik)
        
        if (is.finite(lrt_b)) {
          if (lrt_b < tol) {
            lrt_b <- 0
          }
          
          lrt_boot[b] <- max(lrt_b, 0)
          finite_lrt[b] <- TRUE
        }
      }
    }
  }
  
  lrt_valid <- lrt_boot[is.finite(lrt_boot)]
  
  p_boot <- if (length(lrt_valid) > 0 && is.finite(lrt_obs)) {
    (1 + sum(lrt_valid >= lrt_obs)) / (length(lrt_valid) + 1)
  } else {
    NA_real_
  }
  
  list(
    comparison = paste(null_model, "vs", full_model),
    null_model = null_model,
    full_model = full_model,
    
    lrt_obs = lrt_obs,
    lrt_boot = lrt_boot,
    p_boot = p_boot,
    
    n_finite_lrt = length(lrt_valid),
    finite_lrt_rate = length(lrt_valid) / B,
    
    n_success = length(lrt_valid),
    success_rate = length(lrt_valid) / B,
    
    conv_null_rate = mean(conv_null),
    conv_full_rate = mean(conv_full),
    
    finite_null_logLik_rate = mean(finite_null_logLik),
    finite_full_logLik_rate = mean(finite_full_logLik),
    
    conv_null = conv_null,
    conv_full = conv_full,
    finite_null_logLik = finite_null_logLik,
    finite_full_logLik = finite_full_logLik,
    finite_lrt = finite_lrt
  )
}
bootstrap_four_comparisons_from_results <- function(fileid,
                                                    dom_data,
                                                    results,
                                                    B = 500,
                                                    d = 2,
                                                    seed = NULL,
                                                    verbose = TRUE) {
  
  M <- dom_data[[fileid]]$matrix
  
  results_row <- results %>%
    dplyr::filter(.data$fileid == !!fileid)
  
  if (nrow(results_row) != 1) {
    stop("Expected exactly one row in results for fileid = ", fileid)
  }
  
  comparisons <- list(
    c("LBT",  "LID"),
    c("LBT",  "DABT"),
    c("LID",  "DALID"),
    c("DABT", "DALID")
  )
  
  out <- vector("list", length(comparisons))
  
  for (k in seq_along(comparisons)) {
    null_model <- comparisons[[k]][1]
    full_model <- comparisons[[k]][2]
    
    out[[k]] <- bootstrap_lrt_one_comparison_from_results(
      M = M,
      results_row = results_row,
      null_model = null_model,
      full_model = full_model,
      B = B,
      d = d,
      seed = if (!is.null(seed)) seed + k - 1 else NULL,
      verbose = verbose
    )
  }
  
  names(out) <- vapply(
    comparisons,
    function(x) paste(x[1], "vs", x[2]),
    character(1)
  )
  
  out
}

library(dplyr)
library(purrr)

fileids <- c(
  "Correa_2013a",
  "Alados_1992b",
  "Adcock_2015a",
  "Blatrix_2004c",
  "Poisbleau_2005c",
  "Shimoji_2014c",
  "Kolodziejczyk_2005",
  "Cote_2000d",
  "Cui_2014",
  "Prieto_1978",
  "Mwamende_2009a",
  "Bonanni_2017c",
  "ScottLockhard_1999b"
)

B <- 1000
seed <- 123

bootstrap_dir <- file.path(
  "DomArchiveResults",
  "bootstrap_lrt"
)

dir.create(bootstrap_dir, recursive = TRUE, showWarnings = FALSE)

all_boot_summaries <- list()
failed_bootstraps <- list()

for (fileid in fileids) {
  
  message("\n==================================================")
  message("Running bootstrap for: ", fileid)
  message("==================================================")
  
  res <- try(
    {
      boot_res <- bootstrap_four_comparisons_from_results(
        fileid = fileid,
        dom_data = dom.data,
        results = results,
        B = B,
        seed = seed,
        verbose = TRUE
      )
      
      boot_summary <- summarise_bootstrap_results(boot_res)
      
      boot_summary <- boot_summary %>%
        mutate(
          fileid = fileid,
          B = B,
          seed = seed,
          .before = 1
        )
      
      saveRDS(
        boot_res,
        file = file.path(
          bootstrap_dir,
          paste0("bootstrap_lrt_", fileid, "_B", B, "_s", seed, ".rds")
        )
      )
      
      saveRDS(
        boot_summary,
        file = file.path(
          bootstrap_dir,
          paste0("bootstrap_lrt_summary_", fileid, "_B", B, "_s", seed, ".rds")
        )
      )
      
      boot_summary
    },
    silent = TRUE
  )
  
  if (inherits(res, "try-error")) {
    warning("Bootstrap failed for ", fileid)
    failed_bootstraps[[fileid]] <- res
  } else {
    all_boot_summaries[[fileid]] <- res
  }
}

all_boot_summary <- bind_rows(all_boot_summaries)

saveRDS(
  all_boot_summary,
  file = file.path(
    bootstrap_dir,
    paste0("bootstrap_lrt_summary_all_datasets_B", B, "_s", seed, ".rds")
  )
)

write.csv(
  all_boot_summary,
  file = file.path(
    bootstrap_dir,
    paste0("bootstrap_lrt_summary_all_datasets_B", B, "_s", seed, ".csv")
  ),
  row.names = FALSE
)

if (length(failed_bootstraps) > 0) {
  saveRDS(
    failed_bootstraps,
    file = file.path(
      bootstrap_dir,
      paste0("bootstrap_lrt_failed_datasets_B", B, "_s", seed, ".rds")
    )
  )
}
adcock_lrt_bootstrap <- readRDS("~/DomArchiveResults/bootstrap_lrt/bootstrap_lrt_Adcock_2015a_B1000_s123.rds")
names(adcock_lrt_bootstrap$`LBT vs LID`)
diagnostic_summary <- all_boot_summary %>%
  group_by(comparison) %>%
  summarise(
    n_datasets = n(),
    mean_success_rate = mean(success_rate, na.rm = TRUE),
    mean_conv_null_rate = mean(conv_null_rate, na.rm = TRUE),
    mean_conv_full_rate = mean(conv_full_rate, na.rm = TRUE),
    min_success_rate = min(success_rate, na.rm = TRUE),
    min_conv_null_rate = min(conv_null_rate, na.rm = TRUE),
    min_conv_full_rate = min(conv_full_rate, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    across(
      where(is.numeric) & !matches("n_datasets"),
      ~ sprintf("%.3f", .x)
    )
  )

diagnostic_summary
latex_code <- diagnostic_summary %>%
  kable(
    format = "latex",
    booktabs = TRUE,
    escape = FALSE,
    caption = "Aggregated bootstrap refit diagnostics by model comparison.",
    col.names = c(
      "Comparison",
      "Datasets",
      "Mean finite LRT rate",
      "Mean null convergence rate",
      "Mean full convergence rate",
      "Minimum finite LRT rate",
      "Minimum null convergence rate",
      "Minimum full convergence rate"
    )
  ) %>%
  kable_styling(
    latex_options = c("hold_position", "scale_down")
  )
asdf <- results%>%
  filter(fileid == "ScottLockhard_1999b")
asdf$
cat(latex_code)
library(knitr)
install.packages("kableExtra")
library(kableExtra)

latex_boot_table <- all_boot_summary %>%
  mutate(
    comparison_clean = case_when(
      comparison == "LBT vs LID" ~ "LBT_LID",
      comparison == "LBT vs DABT" ~ "LBT_DABT",
      comparison == "LID vs DALID" ~ "LID_DALID",
      comparison == "DABT vs DALID" ~ "DABT_DALID",
      TRUE ~ comparison
    )
  ) %>%
  select(
    fileid,
    comparison_clean,
    lrt_obs,
    p_boot
  ) %>%
  pivot_wider(
    names_from = comparison_clean,
    values_from = c(lrt_obs, p_boot),
    names_glue = "{comparison_clean}_{.value}"
  ) %>%
  arrange(fileid) %>%
  mutate(
    across(
      starts_with("LBT_") | starts_with("LID_") | starts_with("DABT_"),
      ~ round(.x, 3)
    )
  )

latex_boot_table
latex_code <- latex_boot_table %>%
  kable(
    format = "latex",
    booktabs = TRUE,
    escape = FALSE,
    digits = 3,
    caption = "Parametric bootstrap likelihood ratio test results for selected datasets.",
    col.names = c(
      "Dataset",
      "\\lambda_{\\mathrm{obs}}$ LBT--LID",
      "$\\lambda_{\\mathrm{obs}}$ LBT--DABT",
      "$\\lambda_{\\mathrm{obs}}$ LID--DALID",
      "$\\lambda_{\\mathrm{obs}}$ DABT--DALID",
      "$p$ LBT--LID",
      "$p$ LBT--DABT",
      "$p$ LID--DALID",
      "$p$ DABT--DALID"
    )
  ) %>%
  kable_styling(
    latex_options = c("hold_position", "scale_down")
  )

cat(latex_code)

summarise_one_lrt_comparison <- function(x, dataset_name) {
  
  lrt_boot <- unlist(x$lrt_boot)
  conv_null <- unlist(x$conv_null)
  conv_full <- unlist(x$conv_full)
  finite_null_logLik <- unlist(x$finite_null_logLik)
  finite_full_logLik <- unlist(x$finite_full_logLik)
  finite_lrt <- unlist(x$finite_lrt)
  
  valid <- conv_null &
    finite_null_logLik &
    finite_full_logLik &
    finite_lrt &
    is.finite(lrt_boot)
  
  lrt_boot_valid <- lrt_boot[valid]
  
  B_total <- length(lrt_boot)
  B_valid <- length(lrt_boot_valid)
  
  tibble(
    dataset_name = dataset_name,
    comparison = x$comparison,
    null_model = x$null_model,
    full_model = x$full_model,
    lrt_obs = x$lrt_obs,
    
    B_total = B_total,
    B_valid = B_valid,
    valid_rate = B_valid / B_total,
    
    conv_null_rate_strict = mean(conv_null, na.rm = TRUE),
    conv_full_rate_strict = mean(conv_full, na.rm = TRUE),
    finite_null_logLik_rate = mean(finite_null_logLik, na.rm = TRUE),
    finite_full_logLik_rate = mean(finite_full_logLik, na.rm = TRUE),
    finite_lrt_rate = mean(finite_lrt, na.rm = TRUE),
    
    p_boot_original = x$p_boot,
    
    p_boot_strict = if_else(
      B_valid > 0,
      mean(lrt_boot_valid >= x$lrt_obs),
      NA_real_
    ),
    
    p_boot_strict_plus1 = if_else(
      B_valid > 0,
      (1 + sum(lrt_boot_valid >= x$lrt_obs)) / (B_valid + 1),
      NA_real_
    )
  )
}
highlighted_datasets <- c(
  "Alados_1992b",
  "Adcock_2015a",
  "Blatrix_2004c",
  "Bonanni_2017c",
  "Correa_2013a",
  "Cote_2000d",
  "Cui_2014",
  "Kolodziejczyk_2005",
  "Mwamende_2009a",
  "Poisbleau_2005c",
  "Prieto_1978",
  "ScottLockhard_1999b",
  "Shimoji_2014c"
)

summarise_one_lrt_file <- function(dataset_name,
                                   base_dir = "~/DomArchiveResults/bootstrap_lrt",
                                   B = 1000,
                                   seed = 123) {
  
  file_path <- file.path(
    base_dir,
    paste0("bootstrap_lrt_", dataset_name, "_B", B, "_s", seed, ".rds")
  )
  
  lrt_bootstrap <- readRDS(file_path)
  
  map_dfr(
    lrt_bootstrap,
    ~ summarise_one_lrt_comparison(.x, dataset_name = dataset_name)
  )
}
lrt_bootstrap_summary_strict <- map_dfr(
  highlighted_datasets,
  summarise_one_lrt_file
)
lrt_bootstrap_summary_strict <- lrt_bootstrap_summary_strict %>%
  mutate(
    comparison = factor(
      comparison,
      levels = c(
        "LBT vs LID",
        "LBT vs DABT",
        "LID vs DALID",
        "DABT vs DALID"
      )
    )
  ) %>%
  arrange(dataset_name, comparison)
p_mixture_chisq <- function(lrt) {
  case_when(
    is.na(lrt) ~ NA_real_,
    lrt <= 0 ~ 1,
    TRUE ~ 0.5 * pchisq(lrt, df = 1, lower.tail = FALSE)
  )
}
lrt_bootstrap_summary_strict <- lrt_bootstrap_summary_strict %>%
  mutate(
    lrt_obs_nonneg = pmax(lrt_obs, 0),
    p_asymp_mix = p_mixture_chisq(lrt_obs_nonneg)
  )
lrt_bootstrap_summary_strict %>%
  select(
    dataset_name,
    comparison,
    B_valid,
    p_asymp_mix,
    p_boot_strict_plus1
  )
lrt_bootstrap_table_latex <- lrt_bootstrap_summary_strict %>%
  transmute(
    Dataset = dataset_name,
    Comparison = comparison,
    `$B_{\\mathrm{valid}}$` = B_valid,
    `$p_{\\mathrm{boot}, +1}$` = sprintf("%.3f", p_boot_strict_plus1),
    `$p_{\\mathrm{mix}}$` = sprintf("%.3f", p_asymp_mix)
  )
lrt_bootstrap_table_latex
latex_lrt_bootstrap_table <- lrt_bootstrap_table_latex %>%
  knitr::kable(
    format = "latex",
    booktabs = TRUE,
    escape = FALSE,
    caption = paste(
      "Likelihood-ratio test results.",
      "Bootstrap p-values are computed only from replicates where both the null",
      "and full model converged and both log-likelihood values were finite.",
      "Asymptotic p-values are computed using the mixture reference distribution",
      "$0.5\\chi^2_0 + 0.5\\chi^2_1$."
    ),
    label = "tab:lrt_bootstrap_strict"
  ) %>%
  kableExtra::kable_styling(
    latex_options = c("hold_position", "scale_down")
  )
latex_lrt_bootstrap_table
## =========================================================
## Theoretical CDF and quantile functions for LRT
## =========================================================
library(ggplot2)
theoretical_lrt_cdf <- function(x, distribution = c("chisq1", "mixture_0_1")) {
  distribution <- match.arg(distribution)
  
  if (distribution == "chisq1") {
    return(pchisq(x, df = 1))
  }
  
  if (distribution == "mixture_0_1") {
    return(ifelse(
      x < 0,
      0,
      0.5 + 0.5 * pchisq(x, df = 1)
    ))
  }
}

theoretical_lrt_quantile <- function(p, distribution = c("chisq1", "mixture_0_1")) {
  distribution <- match.arg(distribution)
  
  if (distribution == "chisq1") {
    return(qchisq(p, df = 1))
  }
  
  if (distribution == "mixture_0_1") {
    return(ifelse(
      p <= 0.5,
      0,
      qchisq(2 * p - 1, df = 1)
    ))
  }
}

theoretical_lrt_tail <- function(x, distribution = c("chisq1", "mixture_0_1")) {
  distribution <- match.arg(distribution)
  
  if (!is.finite(x)) {
    return(NA_real_)
  }
  
  if (distribution == "chisq1") {
    return(1 - pchisq(x, df = 1))
  }
  
  if (distribution == "mixture_0_1") {
    if (x <= 0) {
      return(1)
    } else {
      return(0.5 * (1 - pchisq(x, df = 1)))
    }
  }
}
compare_empirical_to_asymptotic <- function(lrt_boot,
                                            lrt_obs = NA_real_,
                                            distribution = c("mixture_0_1", "chisq1"),
                                            trim_p = 1e-6) {
  
  distribution <- match.arg(distribution)
  
  x <- sort(lrt_boot[is.finite(lrt_boot)])
  n <- length(x)
  
  if (n == 0) {
    stop("No finite bootstrap LRT values.")
  }
  
  ## Empirical plotting positions
  p_emp <- (seq_len(n) - 0.5) / n
  p_emp <- pmin(pmax(p_emp, trim_p), 1 - trim_p)
  
  ## Theoretical CDF evaluated at bootstrap values
  F_theory <- theoretical_lrt_cdf(x, distribution)
  
  ## KS distance: max vertical difference between CDFs
  ks_distance <- max(abs(p_emp - F_theory), na.rm = TRUE)
  
  ## Cram??r--von Mises-type distance: average squared CDF difference
  cvm_distance <- mean((p_emp - F_theory)^2, na.rm = TRUE)
  
  ## Quantile RMSE: horizontal quantile discrepancy
  q_theory <- theoretical_lrt_quantile(p_emp, distribution)
  quantile_rmse <- sqrt(mean((x - q_theory)^2, na.rm = TRUE))
  quantile_mae <- mean(abs(x - q_theory), na.rm = TRUE)
  
  ## Bootstrap p-value and asymptotic p-value at observed statistic
  p_boot <- if (is.finite(lrt_obs)) {
    (1 + sum(x >= lrt_obs)) / (n + 1)
  } else {
    NA_real_
  }
  
  p_asymp <- if (is.finite(lrt_obs)) {
    theoretical_lrt_tail(lrt_obs, distribution)
  } else {
    NA_real_
  }
  
  tail_p_diff <- abs(p_boot - p_asymp)
  
  tibble(
    distribution = distribution,
    n_finite_lrt = n,
    ks_distance = ks_distance,
    cvm_distance = cvm_distance,
    quantile_rmse = quantile_rmse,
    quantile_mae = quantile_mae,
    lrt_obs = lrt_obs,
    p_boot = p_boot,
    p_asymp = p_asymp,
    tail_p_diff = tail_p_diff
  )
}
distance_summaries <- list()

for (fileid in fileids) {
  
  boot_res <- readRDS(
    file.path(
      "DomArchiveResults",
      "bootstrap_lrt",
      paste0("bootstrap_lrt_", fileid, "_B1000_s123.rds")
    )
  )
  
  comparison <- "LBT vs LID"
  
  distance_summary <- compare_empirical_to_asymptotic(
    lrt_boot = boot_res[[comparison]]$lrt_boot,
    lrt_obs = boot_res[[comparison]]$lrt_obs,
    distribution = "mixture_0_1"
  ) %>%
    mutate(
      fileid = fileid,
      comparison = comparison,
      .before = 1
    )
  
  distance_summaries[[fileid]] <- distance_summary
}

distance_summaries <- bind_rows(distance_summaries)
distance_summaries
fileid <- "Bonanni_2017c"

boot_res <- readRDS(
  file.path(
    "DomArchiveResults",
    "bootstrap_lrt",
    paste0("bootstrap_lrt_", fileid, "_B1000_s123.rds")
  )
)

comparison <- "LBT vs LID"

distance_summary <- compare_empirical_to_asymptotic(
  lrt_boot = boot_res[[comparison]]$lrt_boot,
  lrt_obs = boot_res[[comparison]]$lrt_obs,
  distribution = "mixture_0_1"
)
make_lrt_qq_pp_plots <- function(boot_res,
                                 comparison = "LBT vs LID",
                                 distribution = c("mixture_0_1", "chisq1"),
                                 fileid = NULL,
                                 trim_p = 1e-6) {
  
  distribution <- match.arg(distribution)
  
  if (!comparison %in% names(boot_res)) {
    stop("comparison must be one of: ", paste(names(boot_res), collapse = ", "))
  }
  
  res <- boot_res[[comparison]]
  
  lrt_boot <- res$lrt_boot
  lrt_boot <- lrt_boot[is.finite(lrt_boot)]
  lrt_boot <- sort(lrt_boot)
  
  n <- length(lrt_boot)
  
  if (n == 0) {
    stop("No finite bootstrap LRT values found.")
  }
  
  ## Plotting positions
  p_emp <- (seq_len(n) - 0.5) / n
  
  ## Avoid numerical issues at exactly 0 or 1
  p_emp <- pmin(pmax(p_emp, trim_p), 1 - trim_p)
  
  q_theory <- theoretical_lrt_quantile(
    p = p_emp,
    distribution = distribution
  )
  
  F_theory_at_boot <- theoretical_lrt_cdf(
    x = lrt_boot,
    distribution = distribution
  )
  
  plot_df <- tibble(
    p_emp = p_emp,
    lrt_boot = lrt_boot,
    q_theory = q_theory,
    F_theory = F_theory_at_boot
  )
  
  title_prefix <- if (is.null(fileid)) {
    comparison
  } else {
    paste(fileid, "-", comparison)
  }
  
  qq_plot <- ggplot(plot_df, aes(x = q_theory, y = lrt_boot)) +
    geom_point(alpha = 0.65, size = 3.6) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    labs(
      x = "Theoretical LRT quantiles",
      y = "Bootstrap LRT quantiles"
    ) +
    theme_bw(base_size = 20)
  
  pp_plot <- ggplot(plot_df, aes(x = F_theory, y = p_emp)) +
    geom_point(alpha = 0.65, size = 3.6) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    labs(
      x = "Theoretical cumulative probability",
      y = "Empirical cumulative probability"
    ) +
    theme_bw(base_size = 20)
  
  list(
    comparison = comparison,
    distribution = distribution,
    fileid = fileid,
    n_finite_lrt = n,
    plot_data = plot_df,
    qq_plot = qq_plot,
    pp_plot = pp_plot
  )
}

fileid <- "Shimoji_2014c"

boot_res <- readRDS(
  file.path(
    "DomArchiveResults",
    "bootstrap_lrt",
    paste0("bootstrap_lrt_", fileid, "_B1000_s123.rds")
  )
)

plots <- make_lrt_qq_pp_plots(
  boot_res = boot_res,
  comparison = "LBT vs LID",
  distribution = "mixture_0_1",
  fileid = fileid
)
plots$qq_plot
plots$pp_plot

