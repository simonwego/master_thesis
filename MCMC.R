library(dplyr)
library(purrr)
library(tidyr)
library(ggplot2)
library(DomArchive)
# ------------------------------------------------------------
# Helper: which models have dyad-level random effects?
# ------------------------------------------------------------

has_theta_random_effects <- function(model) {
  model %in% c("only_theta", "full_theta")
}

# ------------------------------------------------------------
# 1. Unpack joint latent vector b into x and possibly theta_raw
# ------------------------------------------------------------

unpack_b <- function(b_vec, data, model, d = 2) {
  n <- data$n
  n_x <- n * d
  
  x_vec <- b_vec[seq_len(n_x)]
  x <- matrix(x_vec, nrow = n, ncol = d)
  
  if (has_theta_random_effects(model)) {
    q <- length(data$s)
    theta_start <- n_x + 1
    theta_end <- n_x + q
    
    if (length(b_vec) != n_x + q) {
      stop(
        "Wrong length of b_vec. Expected ", n_x + q,
        " but got ", length(b_vec), "."
      )
    }
    
    theta_raw <- b_vec[theta_start:theta_end]
    
    list(
      x = x,
      theta_raw = theta_raw
    )
  } else {
    if (length(b_vec) != n_x) {
      stop(
        "Wrong length of b_vec. Expected ", n_x,
        " but got ", length(b_vec), "."
      )
    }
    
    list(
      x = x
    )
  }
}

# ------------------------------------------------------------
# 2. Pack latent variables into one vector b
# ------------------------------------------------------------

pack_b <- function(x, theta_raw = NULL) {
  if (is.null(theta_raw)) {
    as.vector(x)
  } else {
    c(as.vector(x), as.vector(theta_raw))
  }
}

# ------------------------------------------------------------
# 3. Build parameter list for model function
# ------------------------------------------------------------

make_par_from_b <- function(b_vec, data, theta_hat, model, d = 2) {
  b <- unpack_b(
    b_vec = b_vec,
    data = data,
    model = model,
    d = d
  )
  
  par <- list(
    x = b$x,
    log_r = unname(theta_hat$log_r)
  )
  
  if (model %in% c("full", "full_theta")) {
    par$a <- unname(theta_hat$a)
  }
  
  if (model %in% c("only_theta", "full_theta")) {
    par$log_sigma_theta <- unname(theta_hat$log_sigma_theta)
    par$theta_raw <- b$theta_raw
  }
  
  par
}

# ------------------------------------------------------------
# 4. Conditional NLL for b | y, theta_hat
# ------------------------------------------------------------

make_conditional_b_nll <- function(data, theta_hat, model, fns, d = 2) {
  n <- data$n
  
  if (is.null(n)) {
    stop("data$n is missing.")
  }
  
  if (is.null(data$z) || is.null(data$s) || is.null(data$i) || is.null(data$j)) {
    stop("data must contain z, s, i, j and n. Use make_pair_data() first.")
  }
  
  if (!is.numeric(data$z)) stop("data$z must be numeric.")
  if (!is.numeric(data$s)) stop("data$s must be numeric.")
  
  fn <- fns[[model]]
  
  if (is.null(fn)) {
    stop("No objective function supplied for model = ", model)
  }
  
  force(data)
  force(theta_hat)
  force(model)
  force(fn)
  force(d)
  
  function(b_vec) {
    par <- make_par_from_b(
      b_vec = b_vec,
      data = data,
      theta_hat = theta_hat,
      model = model,
      d = d
    )
    
    val <- fn(parms = par, data = data)
    
    as.numeric(val)
  }
}

# ------------------------------------------------------------
# 5. Random-walk Metropolis-Hastings for b | y, theta_hat
# ------------------------------------------------------------

sample_b_conditional_mh <- function(data,
                                    theta_hat,
                                    model,
                                    fns,
                                    d = 2,
                                    n_iter = 50000,
                                    burn = 5000,
                                    thin = 10,
                                    proposal_sd_x = 0.05,
                                    proposal_sd_theta = 0.05,
                                    adapt = TRUE,
                                    seed = NULL,
                                    verbose = TRUE) {
  if (!is.null(seed)) set.seed(seed)
  
  n <- data$n
  q <- length(data$s)
  n_x <- n * d
  
  target_nll <- make_conditional_b_nll(
    data = data,
    theta_hat = theta_hat,
    model = model,
    fns = fns,
    d = d
  )
  
  # Initial values
  x0 <- matrix(rnorm(n * d, mean = 0, sd = 0.1), nrow = n, ncol = d)
  
  if (has_theta_random_effects(model)) {
    theta_raw0 <- rnorm(q, mean = 0, sd = 0.1)
    b0 <- pack_b(x0, theta_raw0)
  } else {
    theta_raw0 <- NULL
    b0 <- pack_b(x0)
  }
  
  # Find conditional mode of b
  opt <- nlminb(
    start = b0,
    objective = target_nll,
    control = list(eval.max = 1000, iter.max = 1000)
  )
  
  b_current <- opt$par
  nll_current <- target_nll(b_current)
  
  if (!is.finite(nll_current)) {
    stop("Initial NLL is not finite. Check data, parameter estimates and objective function.")
  }
  
  keep_id <- seq(from = burn + 1, to = n_iter, by = thin)
  n_keep <- length(keep_id)
  
  x_samples <- array(NA_real_, dim = c(n_keep, n, d))
  
  if (has_theta_random_effects(model)) {
    theta_raw_samples <- matrix(NA_real_, nrow = n_keep, ncol = q)
  } else {
    theta_raw_samples <- NULL
  }
  
  nll_trace <- rep(NA_real_, n_iter)
  accept_trace <- rep(FALSE, n_iter)
  
  keep_counter <- 0
  
  # Proposal SD vector: possibly different scale for x and theta_raw
  if (has_theta_random_effects(model)) {
    proposal_sd_vec <- c(
      rep(proposal_sd_x, n_x),
      rep(proposal_sd_theta, q)
    )
  } else {
    proposal_sd_vec <- rep(proposal_sd_x, n_x)
  }
  
  for (iter in seq_len(n_iter)) {
    b_prop <- b_current + rnorm(length(b_current), mean = 0, sd = proposal_sd_vec)
    
    nll_prop <- target_nll(b_prop)
    
    log_alpha <- -nll_prop + nll_current
    
    if (is.finite(log_alpha) && log(runif(1)) < log_alpha) {
      b_current <- b_prop
      nll_current <- nll_prop
      accept_trace[iter] <- TRUE
    }
    
    nll_trace[iter] <- nll_current
    
    if (iter %in% keep_id) {
      keep_counter <- keep_counter + 1
      
      b_unpacked <- unpack_b(
        b_vec = b_current,
        data = data,
        model = model,
        d = d
      )
      
      x_samples[keep_counter, , ] <- b_unpacked$x
      
      if (has_theta_random_effects(model)) {
        theta_raw_samples[keep_counter, ] <- b_unpacked$theta_raw
      }
    }
    
    # Simple adaptation during burn-in
    if (adapt && iter <= burn && iter %% 200 == 0) {
      acc_rate_window <- mean(accept_trace[(iter - 199):iter])
      
      if (acc_rate_window < 0.15) {
        proposal_sd_vec <- proposal_sd_vec * 0.8
      }
      
      if (acc_rate_window > 0.35) {
        proposal_sd_vec <- proposal_sd_vec * 1.2
      }
    }
    
    if (verbose && iter %% 5000 == 0) {
      message(
        "iter = ", iter,
        ", acceptance = ", round(mean(accept_trace[1:iter]), 3),
        ", proposal_sd median = ", signif(median(proposal_sd_vec), 3),
        ", nll = ", signif(nll_current, 5)
      )
    }
  }
  
  mode_unpacked <- unpack_b(
    b_vec = opt$par,
    data = data,
    model = model,
    d = d
  )
  
  list(
    model = model,
    samples = x_samples,
    x_samples = x_samples,
    theta_raw_samples = theta_raw_samples,
    nll_trace = nll_trace,
    accept_trace = accept_trace,
    acceptance_rate = mean(accept_trace),
    proposal_sd_vec_final = proposal_sd_vec,
    proposal_sd_x_final = median(proposal_sd_vec[seq_len(n_x)]),
    proposal_sd_theta_final = if (has_theta_random_effects(model)) {
      median(proposal_sd_vec[(n_x + 1):length(proposal_sd_vec)])
    } else {
      NA_real_
    },
    mode = mode_unpacked$x,
    x_mode = mode_unpacked$x,
    theta_raw_mode = if (has_theta_random_effects(model)) mode_unpacked$theta_raw else NULL,
    b_mode = opt$par,
    mode_nll = opt$objective,
    opt = opt
  )
}


# ============================================================
# Model components
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
# Linear predictors
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
# Models
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

fns <- list(
  null = f_null,
  full = f_full,
  only_theta = f_only_theta,
  full_theta = f_full_theta
)

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

mcmc_x_samples_to_df <- function(mcmc, fileid, model) {
  x_samples <- mcmc$x_samples
  
  if (is.null(x_samples)) {
    x_samples <- mcmc$samples
  }
  
  if (length(dim(x_samples)) != 3) {
    stop("x_samples must be a 3D array: draw x agent x dimension.")
  }
  
  dims <- dim(x_samples)
  
  expand.grid(
    draw = seq_len(dims[1]),
    agent = seq_len(dims[2]),
    dimension = seq_len(dims[3])
  ) |>
    as_tibble() |>
    mutate(
      x = as.vector(x_samples),
      fileid = fileid,
      model = model
    ) |>
    select(fileid, model, draw, agent, dimension, x)
}

mcmc_theta_samples_to_df <- function(mcmc, data, fileid, model) {
  if (is.null(mcmc$theta_raw_samples)) {
    return(tibble())
  }
  
  theta_samples <- mcmc$theta_raw_samples
  q <- ncol(theta_samples)
  
  dyad_lookup <- tibble(
    theta_index = seq_len(q),
    i = data$i,
    j = data$j,
    s = data$s,
    z = data$z
  )
  
  expand.grid(
    draw = seq_len(nrow(theta_samples)),
    theta_index = seq_len(q)
  ) |>
    as_tibble() |>
    mutate(
      theta_raw = as.vector(theta_samples),
      fileid = fileid,
      model = model
    ) |>
    left_join(dyad_lookup, by = "theta_index") |>
    select(fileid, model, draw, theta_index, i, j, z, s, theta_raw)
}

make_latent_mcmc_diagnostics <- function(mcmc, fileid, model, data_i, results, burn) {
  res_row <- results |> filter(.data$fileid == !!fileid)
  
  n_iter <- length(mcmc$nll_trace)
  post_idx <- (burn + 1):n_iter
  
  tibble(
    fileid = fileid,
    model = model,
    
    sampled_latent_state = ifelse(
      has_theta_random_effects(model),
      "x + theta_raw",
      "x"
    ),
    
    n_agents = data_i$n,
    n_observed_dyads = data_i$n_dyads_observed,
    n_interactions = data_i$n_interactions,
    
    n_iter = n_iter,
    burn = burn,
    n_saved = dim(mcmc$x_samples)[1],
    n_x_dim = prod(dim(mcmc$x_samples)[2:3]),
    n_theta_raw_dim = ifelse(
      is.null(mcmc$theta_raw_samples),
      0L,
      ncol(mcmc$theta_raw_samples)
    ),
    n_total_latent_dim = length(mcmc$b_mode),
    
    acceptance_rate = mcmc$acceptance_rate,
    proposal_sd_x_final = mcmc$proposal_sd_x_final,
    proposal_sd_theta_final = mcmc$proposal_sd_theta_final,
    
    mode_nll = mcmc$mode_nll,
    final_nll = tail(mcmc$nll_trace, 1),
    mean_nll_after_burn = mean(mcmc$nll_trace[post_idx]),
    sd_nll_after_burn = sd(mcmc$nll_trace[post_idx]),
    
    opt_convergence = mcmc$opt$convergence,
    opt_message = mcmc$opt$message,
    
    model_logLik = res_row[[paste0(model, "_logLik")]][1],
    model_AIC = res_row[[paste0(model, "_AIC")]][1],
    model_max_grad = res_row[[paste0(model, "_max_grad")]][1],
    model_hess_min_eig = res_row[[paste0(model, "_hess_min_eig")]][1],
    model_hess_abs_cond = res_row[[paste0(model, "_hess_abs_cond")]][1],
    model_hess_status = res_row[[paste0(model, "_hess_status")]][1],
    
    log_r_hat = res_row[[paste0(model, "_log_r_est")]][1],
    r_hat = res_row[[paste0(model, "_r_est")]][1],
    a_hat = if (model %in% c("full", "full_theta")) {
      res_row[[paste0(model, "_a_est")]][1]
    } else {
      NA_real_
    },
    log_sigma_theta_hat = if (model %in% c("only_theta", "full_theta")) {
      res_row[[paste0(model, "_log_sigma_theta_est")]][1]
    } else {
      NA_real_
    },
    sigma_theta_hat = if (model %in% c("only_theta", "full_theta")) {
      res_row[[paste0(model, "_sigma_theta_est")]][1]
    } else {
      NA_real_
    }
  )
}

run_one_latent_mcmc <- function(fileid,
                                model,
                                dom.data,
                                results,
                                fns,
                                settings,
                                out_dir,
                                seed_base = 1000) {
  message("\n==================================================")
  message("Running latent MCMC for fileid = ", fileid)
  message("Model = ", model)
  message("==================================================")
  
  tryCatch({
    M_i <- dom.data[[fileid]]$matrix
    
    if (is.null(M_i)) {
      stop("No matrix found for fileid = ", fileid)
    }
    
    data_i <- make_pair_data(M_i)
    
    theta_hat_i <- extract_theta_hat(
      results = results,
      fileid = fileid,
      model = model
    )
    
    seed_i <- seed_base + match(fileid, settings$fileids)
    
    mcmc_i <- sample_b_conditional_mh(
      data = data_i,
      theta_hat = theta_hat_i,
      model = model,
      fns = fns,
      d = settings$d,
      n_iter = settings$n_iter,
      burn = settings$burn,
      thin = settings$thin,
      proposal_sd_x = settings$proposal_sd_x,
      proposal_sd_theta = settings$proposal_sd_theta,
      adapt = settings$adapt,
      seed = seed_i,
      verbose = settings$verbose)
    
   
    x_samples_i <- mcmc_x_samples_to_df(
      mcmc = mcmc_i,
      fileid = fileid,
      model = model
    )
    
    theta_samples_i <- mcmc_theta_samples_to_df(
      mcmc = mcmc_i,
      data = data_i,
      fileid = fileid,
      model = model
    )
    
    diagnostics_i <- make_latent_mcmc_diagnostics(
      mcmc = mcmc_i,
      fileid = fileid,
      model = model,
      data_i = data_i,
      results = results,
      burn = settings$burn
    ) |>
      mutate(
        status = "ok",
        error_message = NA_character_
      )
    
    saveRDS(
      object = list(
        fileid = fileid,
        model = model,
        data = data_i,
        theta_hat = theta_hat_i,
        mcmc = mcmc_i,
        x_samples = x_samples_i,
        theta_samples = theta_samples_i,
        diagnostics = diagnostics_i
      ),
      file = file.path(out_dir, paste0("latent_mcmc_", model, "_", fileid, ".rds"))
    )
    
    list(
      fileid = fileid,
      model = model,
      status = "ok",
      data = data_i,
      theta_hat = theta_hat_i,
      mcmc = mcmc_i,
      x_samples = x_samples_i,
      theta_samples = theta_samples_i,
      diagnostics = diagnostics_i,
      error = NULL
    )
    
  }, error = function(e) {
    message("ERROR for fileid = ", fileid, ", model = ", model, ": ", conditionMessage(e))
    
    diagnostics_i <- tibble(
      fileid = fileid,
      model = model,
      sampled_latent_state = ifelse(
        has_theta_random_effects(model),
        "x + theta_raw",
        "x"
      ),
      n_agents = NA_integer_,
      n_observed_dyads = NA_integer_,
      n_interactions = NA_real_,
      n_iter = settings$n_iter,
      burn = settings$burn,
      n_saved = NA_integer_,
      n_x_dim = NA_integer_,
      n_theta_raw_dim = NA_integer_,
      n_total_latent_dim = NA_integer_,
      acceptance_rate = NA_real_,
      proposal_sd_x_final = NA_real_,
      proposal_sd_theta_final = NA_real_,
      mode_nll = NA_real_,
      final_nll = NA_real_,
      mean_nll_after_burn = NA_real_,
      sd_nll_after_burn = NA_real_,
      opt_convergence = NA_integer_,
      opt_message = NA_character_,
      model_logLik = NA_real_,
      model_AIC = NA_real_,
      model_max_grad = NA_real_,
      model_hess_min_eig = NA_real_,
      model_hess_abs_cond = NA_real_,
      model_hess_status = NA_character_,
      log_r_hat = NA_real_,
      r_hat = NA_real_,
      a_hat = NA_real_,
      log_sigma_theta_hat = NA_real_,
      sigma_theta_hat = NA_real_,
      status = "error",
      error_message = conditionMessage(e)
    )
    
    list(
      fileid = fileid,
      model = model,
      status = "error",
      data = NULL,
      theta_hat = NULL,
      mcmc = NULL,
      x_samples = tibble(),
      theta_samples = tibble(),
      diagnostics = diagnostics_i,
      error = e
    )
  })
}
###################################
# Running MCMC
###################################

fileid <- "Blatrix_2004c"
model <- "full"

M_i <- dom.data[[fileid]]$matrix
data_i <- make_pair_data(M_i)

theta_hat_i <- extract_theta_hat(
  results = results,
  fileid = fileid,
  model = model
)

mcmc_blatrix_full <- sample_b_conditional_mh(
  data = data_i,
  theta_hat = theta_hat_i,
  model = model,
  fns = fns,
  seed = 1
)

dim(mcmc_blatrix_full$x_samples)
mcmc_blatrix_full$theta_raw_samples
table <- mcmc_x_samples_to_df(mcmc_blatrix_full, fileid, model)
make_latent_mcmc_diagnostics

fileid <- "Blatrix_2004c"
model <- "full_theta"

M_i <- dom.data[[fileid]]$matrix
data_i <- make_pair_data(M_i)

theta_hat_i <- extract_theta_hat(
  results = results,
  fileid = fileid,
  model = model
)


mcmc_blatrix_full_theta <- sample_b_conditional_mh(
  data = data_i,
  theta_hat = theta_hat_i,
  model = model,
  fns = fns,
  seed = 1
)

dim(mcmc_blatrix_full_theta$x_samples)
dim(mcmc_blatrix_full_theta$theta_raw_samples)
mcmc_blatrix_full_theta$acceptance_rate
asdf <- results%>% 
  filter(null_r_est < 0.1)
asdf
mcmc_settings <- list(
  fileids = c(
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
    "ScottLockhard_1999b"
  ),
  d = 2,
  n_iter = 100000,
  burn = 10000,
  thin = 1,
  proposal_sd_x = 0.04,
  proposal_sd_theta = 0.04,
  adapt = TRUE,
  verbose = TRUE
)

out_dir <- "latent_mcmc_results_no_thinning"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

models_to_run <- c("null", "full", "only_theta", "full_theta")

all_latent_mcmc_runs <- list()

for (model_i in models_to_run) {
  
  runs_i <- vector("list", length(mcmc_settings$fileids))
  names(runs_i) <- mcmc_settings$fileids
  
  for (fileid_i in mcmc_settings$fileids) {
    runs_i[[fileid_i]] <- run_one_latent_mcmc(
      fileid = fileid_i,
      model = model_i,
      dom.data = dom.data,
      results = results,
      fns = fns,
      settings = mcmc_settings,
      out_dir = out_dir,
      seed_base = 1000
    )
  }
  
  all_latent_mcmc_runs[[model_i]] <- runs_i
}

#############################################
# Effective sample size og IAT
#############################################

safe_iat <- function(x, max_lag = 1000) {
  x <- x[is.finite(x)]
  n <- length(x)
  
  if (n < 10) return(NA_real_)
  if (!is.finite(sd(x)) || sd(x) == 0) return(NA_real_)
  
  max_lag <- min(max_lag, n - 1)
  
  acf_vals <- as.numeric(stats::acf(
    x,
    lag.max = max_lag,
    plot = FALSE,
    na.action = na.pass
  )$acf)[-1]
  
  # Initial positive sequence truncation:
  # stop summing when autocorrelation first becomes non-positive.
  first_nonpositive <- which(acf_vals <= 0)[1]
  
  if (!is.na(first_nonpositive)) {
    if (first_nonpositive == 1) return(1)
    acf_used <- acf_vals[seq_len(first_nonpositive - 1)]
  } else {
    acf_used <- acf_vals
  }
  
  iat <- 1 + 2 * sum(acf_used)
  
  if (!is.finite(iat) || iat < 1) return(1)
  
  iat
}

safe_ess <- function(x, max_lag = 1000) {
  x <- x[is.finite(x)]
  iat <- safe_iat(x, max_lag = max_lag)
  
  if (!is.finite(iat) || iat <= 0) return(NA_real_)
  
  length(x) / iat
}

# ------------------------------------------------------------
# 2. ESS / IAT for x_samples in one run
# ------------------------------------------------------------

compute_x_ess_iat_one_run <- function(run_i, max_lag = 1000) {
  if (is.null(run_i$mcmc)) {
    return(tibble())
  }
  
  x_samples <- run_i$mcmc$x_samples
  
  if (is.null(x_samples)) {
    x_samples <- run_i$mcmc$samples
  }
  
  if (is.null(x_samples)) {
    return(tibble())
  }
  
  if (length(dim(x_samples)) != 3) {
    stop("x_samples must be a 3D array: draw x agent x dimension.")
  }
  
  dims <- dim(x_samples)
  n_draws <- dims[1]
  n_agents <- dims[2]
  d <- dims[3]
  
  map_dfr(seq_len(n_agents), function(agent_i) {
    map_dfr(seq_len(d), function(dim_i) {
      chain <- x_samples[, agent_i, dim_i]
      
      tibble(
        fileid = run_i$fileid,
        model = run_i$model,
        latent_type = "x",
        agent = agent_i,
        dimension = dim_i,
        theta_index = NA_integer_,
        i = NA_integer_,
        j = NA_integer_,
        n_draws = n_draws,
        ess = safe_ess(chain, max_lag = max_lag),
        iat = safe_iat(chain, max_lag = max_lag),
        mean = mean(chain, na.rm = TRUE),
        sd = sd(chain, na.rm = TRUE)
      )
    })
  })
}

# ------------------------------------------------------------
# 3. ESS / IAT for theta_raw_samples in one run
# ------------------------------------------------------------

compute_theta_ess_iat_one_run <- function(run_i, max_lag = 1000) {
  if (is.null(run_i$mcmc) || is.null(run_i$mcmc$theta_raw_samples)) {
    return(tibble())
  }
  
  theta_samples <- run_i$mcmc$theta_raw_samples
  
  if (!is.matrix(theta_samples)) {
    stop("theta_raw_samples must be a matrix: draw x theta_index.")
  }
  
  n_draws <- nrow(theta_samples)
  q <- ncol(theta_samples)
  
  # Optional dyad lookup, if data is present
  if (!is.null(run_i$data)) {
    i_vec <- run_i$data$i
    j_vec <- run_i$data$j
  } else {
    i_vec <- rep(NA_integer_, q)
    j_vec <- rep(NA_integer_, q)
  }
  
  map_dfr(seq_len(q), function(theta_i) {
    chain <- theta_samples[, theta_i]
    
    tibble(
      fileid = run_i$fileid,
      model = run_i$model,
      latent_type = "theta_raw",
      agent = NA_integer_,
      dimension = NA_integer_,
      theta_index = theta_i,
      i = i_vec[theta_i],
      j = j_vec[theta_i],
      n_draws = n_draws,
      ess = safe_ess(chain, max_lag = max_lag),
      iat = safe_iat(chain, max_lag = max_lag),
      mean = mean(chain, na.rm = TRUE),
      sd = sd(chain, na.rm = TRUE)
    )
  })
}

# ------------------------------------------------------------
# 4. ESS / IAT for all runs
# ------------------------------------------------------------

compute_latent_ess_iat <- function(all_latent_mcmc_runs,
                                   max_lag = 1000,
                                   include_x = TRUE,
                                   include_theta = TRUE) {
  map_dfr(all_latent_mcmc_runs, function(model_runs) {
    map_dfr(model_runs, function(run_i) {
      out <- list()
      
      if (include_x) {
        out$x <- compute_x_ess_iat_one_run(
          run_i = run_i,
          max_lag = max_lag
        )
      }
      
      if (include_theta) {
        out$theta <- compute_theta_ess_iat_one_run(
          run_i = run_i,
          max_lag = max_lag
        )
      }
      
      bind_rows(out)
    })
  }, .id = "model_group")
}

latent_ess_iat <- compute_latent_ess_iat(
  all_latent_mcmc_runs = all_latent_mcmc_runs,
  max_lag = 1000,
  include_x = TRUE,
  include_theta = TRUE
)

latent_ess_iat

# ------------------------------------------------------------
# 5. Summary table by fileid, model and latent type
# ------------------------------------------------------------

latent_ess_iat_summary <- latent_ess_iat |>
  group_by(fileid, model, latent_type) |>
  summarise(
    n_coordinates = n(),
    n_draws = first(n_draws),
    
    median_ess = median(ess, na.rm = TRUE),
    min_ess = min(ess, na.rm = TRUE),
    q10_ess = quantile(ess, 0.10, na.rm = TRUE),
    q90_ess = quantile(ess, 0.90, na.rm = TRUE),
    
    median_iat = median(iat, na.rm = TRUE),
    max_iat = max(iat, na.rm = TRUE),
    q90_iat = quantile(iat, 0.90, na.rm = TRUE),
    
    median_sd = median(sd, na.rm = TRUE),
    max_sd = max(sd, na.rm = TRUE),
    
    prop_ess_below_50 = mean(ess < 50, na.rm = TRUE),
    prop_ess_below_100 = mean(ess < 100, na.rm = TRUE),
    prop_ess_below_200 = mean(ess < 200, na.rm = TRUE),
    
    .groups = "drop"
  ) |>
  arrange(model, fileid, latent_type)

latent_ess_iat_summary

latent_ess_iat_summary_wide <- latent_ess_iat_summary |>
  select(
    fileid,
    model,
    latent_type,
    n_coordinates,
    median_ess,
    min_ess,
    median_iat,
    max_iat,
    prop_ess_below_100
  ) |>
  tidyr::pivot_wider(
    names_from = latent_type,
    values_from = c(
      n_coordinates,
      median_ess,
      min_ess,
      median_iat,
      max_iat,
      prop_ess_below_100
    )
  ) |>
  arrange(model, fileid)

latent_ess_iat_summary_wide

saveRDS(
  latent_ess_iat,
  file.path(out_dir, "latent_ess_iat_by_coordinate.rds")
)

saveRDS(
  latent_ess_iat_summary,
  file.path(out_dir, "latent_ess_iat_summary.rds")
)

saveRDS(
  latent_ess_iat_summary_wide,
  file.path(out_dir, "latent_ess_iat_summary_wide.rds")
)

latent_ess_iat_latex_table <- latent_ess_iat_summary_wide |>
  mutate(
    Model = recode(
      model,
      "null" = "LBT",
      "full" = "LID",
      "only_theta" = "DABT",
      "full_theta" = "DALID",
      .default = model
    )
  ) |>
  transmute(
    Dataset = fileid,
    Model,
    `$n_x$` = n_coordinates_x,
    `Med. ESS$_x$` = median_ess_x,
    `Min. ESS$_x$` = min_ess_x,
    `Med. IAT$_x$` = median_iat_x,
    `$p(\\mathrm{ESS}_x<100)$` = prop_ess_below_100_x,
    `$n_\\theta$` = n_coordinates_theta_raw,
    `Med. ESS$_\\theta$` = median_ess_theta_raw,
    `Min. ESS$_\\theta$` = min_ess_theta_raw,
    `Med. IAT$_\\theta$` = median_iat_theta_raw,
    `$p(\\mathrm{ESS}_\\theta<100)$` = prop_ess_below_100_theta_raw
  ) |>
  mutate(
    across(
      where(is.numeric),
      ~ round(.x, 2)
    ),
    across(
      everything(),
      ~ ifelse(is.na(.x), "--", as.character(.x))
    )
  ) |>
  arrange(Model, Dataset)

latent_ess_iat_latex_code <- latent_ess_iat_latex_table |>
  kable(
    format = "latex",
    booktabs = TRUE,
    escape = FALSE,
    caption = paste(
      "Effective sample size and integrated autocorrelation time diagnostics",
      "for the conditional MCMC samples. The table reports summaries for",
      "the latent position coordinates $x$ and, for the random-effect models,",
      "the dyad-level random effects $\\theta^{\\mathrm{raw}}$."
    ),
    label = "tab:latent_ess_iat_summary"
  ) |>
  kable_styling(
    latex_options = c("hold_position", "scale_down")
  )

latent_ess_iat_latex_code

# ------------------------------------------------------------
# 6. Safe helpers
# ------------------------------------------------------------

safe_skewness <- function(x) {
  x <- x[is.finite(x)]
  n <- length(x)
  
  if (n < 3) return(NA_real_)
  
  m <- mean(x)
  xc <- x - m
  m2 <- mean(xc^2)
  m3 <- mean(xc^3)
  
  if (!is.finite(m2) || m2 <= 0) return(NA_real_)
  
  m3 / m2^(3 / 2)
}

safe_kurtosis <- function(x) {
  x <- x[is.finite(x)]
  n <- length(x)
  
  if (n < 4) return(NA_real_)
  
  m <- mean(x)
  xc <- x - m
  m2 <- mean(xc^2)
  m4 <- mean(xc^4)
  
  if (!is.finite(m2) || m2 <= 0) return(NA_real_)
  
  m4 / m2^2
}

safe_excess_kurtosis <- function(x) {
  k <- safe_kurtosis(x)
  if (!is.finite(k)) return(NA_real_)
  k - 3
}

safe_median <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  median(x)
}

compute_x_shape_one_run <- function(run_i) {
  if (is.null(run_i$mcmc)) {
    return(tibble())
  }
  
  x_samples <- run_i$mcmc$x_samples
  
  if (is.null(x_samples)) {
    x_samples <- run_i$mcmc$samples
  }
  
  if (is.null(x_samples)) {
    return(tibble())
  }
  
  if (length(dim(x_samples)) != 3) {
    stop("x_samples must be a 3D array: draw x agent x dimension.")
  }
  
  dims <- dim(x_samples)
  n_agents <- dims[2]
  d <- dims[3]
  
  map_dfr(seq_len(n_agents), function(agent_i) {
    map_dfr(seq_len(d), function(dim_i) {
      chain <- x_samples[, agent_i, dim_i]
      
      tibble(
        fileid = run_i$fileid,
        model = run_i$model,
        latent_type = "x",
        agent = agent_i,
        dimension = dim_i,
        skew = safe_skewness(chain),
        excess_kurtosis = safe_excess_kurtosis(chain)
      )
    })
  })
}

compute_theta_shape_one_run <- function(run_i) {
  if (is.null(run_i$mcmc) || is.null(run_i$mcmc$theta_raw_samples)) {
    return(tibble())
  }
  
  theta_samples <- run_i$mcmc$theta_raw_samples
  
  if (!is.matrix(theta_samples)) {
    stop("theta_raw_samples must be a matrix: draw x theta_index.")
  }
  
  q <- ncol(theta_samples)
  
  map_dfr(seq_len(q), function(theta_i) {
    chain <- theta_samples[, theta_i]
    
    tibble(
      fileid = run_i$fileid,
      model = run_i$model,
      latent_type = "theta_raw",
      theta_index = theta_i,
      skew = safe_skewness(chain),
      excess_kurtosis = safe_excess_kurtosis(chain)
    )
  })
}

compute_shape_from_all_runs <- function(all_latent_mcmc_runs,
                                        include_x = TRUE,
                                        include_theta = TRUE) {
  map_dfr(all_latent_mcmc_runs, function(model_runs) {
    map_dfr(model_runs, function(run_i) {
      out <- list()
      
      if (include_x) {
        out$x <- compute_x_shape_one_run(run_i)
      }
      
      if (include_theta) {
        out$theta <- compute_theta_shape_one_run(run_i)
      }
      
      bind_rows(out)
    })
  }, .id = "model_group")
}

latent_shape_by_coordinate <- compute_shape_from_all_runs(
  all_latent_mcmc_runs = all_latent_mcmc_runs,
  include_x = TRUE,
  include_theta = TRUE
)

latent_shape_summary_wide <- latent_shape_by_coordinate |>
  group_by(fileid, model, latent_type) |>
  summarise(
    n_shape_coordinates = n(),
    median_abs_skew = safe_median(abs(skew)),
    median_abs_excess_kurtosis = safe_median(abs(excess_kurtosis)),
    max_abs_skew = {
      x <- abs(skew)
      x <- x[is.finite(x)]
      if (length(x) == 0) NA_real_ else max(x)
    },
    max_abs_excess_kurtosis = {
      x <- abs(excess_kurtosis)
      x <- x[is.finite(x)]
      if (length(x) == 0) NA_real_ else max(x)
    },
    .groups = "drop"
  ) |>
  pivot_wider(
    names_from = latent_type,
    values_from = c(
      n_shape_coordinates,
      median_abs_skew,
      median_abs_excess_kurtosis,
      max_abs_skew,
      max_abs_excess_kurtosis
    )
  )

mcmc_final_summary <- all_diagnostics |>
  select(
    fileid,
    model,
    status,
    sampled_latent_state,
    n_total_latent_dim,
    acceptance_rate
  ) |>
  left_join(
    latent_shape_summary_wide,
    by = c("fileid", "model")
  ) |>
  left_join(
    latent_ess_iat_summary_wide |>
      select(
        fileid,
        model,
        n_coordinates_x,
        n_coordinates_theta_raw,
        median_ess_x,
        median_ess_theta_raw,
        median_iat_x,
        median_iat_theta_raw
      ),
    by = c("fileid", "model")
  ) |>
  mutate(
    Model = recode(
      model,
      "null" = "LBT",
      "full" = "LID",
      "only_theta" = "DABT",
      "full_theta" = "DALID",
      .default = model
    )
  ) |>
  transmute(
    Dataset = fileid,
    Model,
    Status = status,
    `Latent state` = sampled_latent_state,
    `Latent dim.` = n_total_latent_dim,
    `Acc.` = acceptance_rate,
    
    `$n_x$` = n_coordinates_x,
    `Med. |skew_x|` = median_abs_skew_x,
    `Med. |ex.kurt._x|` = median_abs_excess_kurtosis_x,
    `Med. ESS_x` = median_ess_x,
    `Med. IAT_x` = median_iat_x,
    
    `$n_theta$` = n_coordinates_theta_raw,
    `Med. |skew_theta|` = median_abs_skew_theta_raw,
    `Med. |ex.kurt._theta|` = median_abs_excess_kurtosis_theta_raw,
    `Med. ESS_theta` = median_ess_theta_raw,
    `Med. IAT_theta` = median_iat_theta_raw
  ) |>
  arrange(Model, Dataset)

mcmc_final_summary

mcmc_final_summary_compact <- mcmc_final_summary |>
  select(
    Dataset,
    Model,
    `Latent dim.`,
    `Acc.`,
    `Med. |skew_x|`,
    `Med. |ex.kurt._x|`,
    `Med. ESS_x`,
    `Med. IAT_x`,
    `Med. |skew_theta|`,
    `Med. |ex.kurt._theta|`,
    `Med. ESS_theta`,
    `Med. IAT_theta`
  )

mcmc_final_summary_compact

mcmc_final_summary_latex <- mcmc_final_summary_compact |>
  mutate(
    across(where(is.numeric), ~ round(.x, 3)),
    across(everything(), ~ ifelse(is.na(.x), "--", as.character(.x)))
  ) |>
  rename(
    `Latent dim.` = `Latent dim.`,
    `Acc.` = `Acc.`,
    `Med. $|S_x|$` = `Med. |skew_x|`,
    `Med. $|K_x-3|$` = `Med. |ex.kurt._x|`,
    `Med. ESS$_x$` = `Med. ESS_x`,
    `Med. IAT$_x$` = `Med. IAT_x`,
    `Med. $|S_\\theta|$` = `Med. |skew_theta|`,
    `Med. $|K_\\theta-3|$` = `Med. |ex.kurt._theta|`,
    `Med. ESS$_\\theta$` = `Med. ESS_theta`,
    `Med. IAT$_\\theta$` = `Med. IAT_theta`
  )

mcmc_final_summary_latex_code <- mcmc_final_summary_latex |>
  kable(
    format = "latex",
    booktabs = TRUE,
    escape = FALSE,
    caption = paste(
      "Summary of conditional MCMC diagnostics for the selected datasets and models.",
      "The table reports the latent dimension, acceptance rate, median absolute",
      "skewness, median absolute excess kurtosis, median effective sample size",
      "and median integrated autocorrelation time. Shape diagnostics are computed",
      "directly from the full retained MCMC arrays."
    ),
    label = "tab:mcmc_final_summary"
  ) |>
  kable_styling(
    latex_options = c("hold_position", "scale_down")
  )

mcmc_final_summary_latex_code

# ------------------------------------------------------------
# Plotting marginal 
# ------------------------------------------------------------

get_mcmc_run <- function(all_latent_mcmc_runs, fileid, model) {
  run_i <- all_latent_mcmc_runs[[model]][[fileid]]
  
  if (is.null(run_i)) {
    stop("Could not find run for model = ", model, ", fileid = ", fileid)
  }
  if (is.null(run_i$mcmc)) {
    stop("Run exists, but run_i$mcmc is NULL for model = ", model, ", fileid = ", fileid)
  }
  
  run_i
}

extract_x_marginal_samples <- function(run_i, components = NULL) {
  x_samples <- run_i$mcmc$x_samples
  
  if (is.null(x_samples)) {
    x_samples <- run_i$mcmc$samples
  }
  
  if (is.null(x_samples) || length(dim(x_samples)) != 3) {
    stop("x_samples must be a 3D array: draw x agent x dimension.")
  }
  
  n_draws <- dim(x_samples)[1]
  n_agents <- dim(x_samples)[2]
  d <- dim(x_samples)[3]
  
  if (is.null(components)) {
    components <- expand.grid(
      agent = seq_len(n_agents),
      dimension = seq_len(d)
    ) |>
      as_tibble()
  }
  
  purrr::map_dfr(seq_len(nrow(components)), function(k) {
    a <- components$agent[k]
    dim_i <- components$dimension[k]
    chain <- x_samples[, a, dim_i]
    
    tibble(
      draw = seq_len(n_draws),
      value = chain,
      agent = a,
      dimension = dim_i,
      component_label = paste0("x[", a, ",", dim_i, "]")
    )
  })
}

extract_theta_marginal_samples <- function(run_i, theta_indices = NULL) {
  theta_samples <- run_i$mcmc$theta_raw_samples
  
  if (is.null(theta_samples)) {
    stop("No theta_raw_samples found in this run.")
  }
  
  n_draws <- nrow(theta_samples)
  q <- ncol(theta_samples)
  
  if (is.null(theta_indices)) {
    theta_indices <- seq_len(q)
  }
  
  if (!is.null(run_i$data)) {
    i_vec <- run_i$data$i
    j_vec <- run_i$data$j
  } else {
    i_vec <- rep(NA_integer_, q)
    j_vec <- rep(NA_integer_, q)
  }
  
  purrr::map_dfr(theta_indices, function(idx) {
    chain <- theta_samples[, idx]
    
    tibble(
      draw = seq_len(n_draws),
      value = chain,
      theta_index = idx,
      i = i_vec[idx],
      j = j_vec[idx],
      component_label = ifelse(
        is.na(i_vec[idx]) || is.na(j_vec[idx]),
        paste0("theta[", idx, "]"),
        paste0("theta[", idx, "] (", i_vec[idx], ",", j_vec[idx], ")")
      )
    )
  })
}

select_top_x_components_by_skew <- function(run_i, n_select = 6) {
  x_samples <- run_i$mcmc$x_samples
  
  if (is.null(x_samples)) {
    x_samples <- run_i$mcmc$samples
  }
  
  n_agents <- dim(x_samples)[2]
  d <- dim(x_samples)[3]
  
  comps <- expand.grid(
    agent = seq_len(n_agents),
    dimension = seq_len(d)
  ) |>
    as_tibble()
  
  comps |>
    rowwise() |>
    mutate(
      skew = safe_skewness(x_samples[, agent, dimension]),
      abs_skew = abs(skew),
      excess_kurtosis = safe_excess_kurtosis(x_samples[, agent, dimension]),
      abs_excess_kurtosis = abs(excess_kurtosis)
    ) |>
    ungroup() |>
    arrange(desc(abs_skew), desc(abs_excess_kurtosis)) |>
    slice_head(n = n_select)
}

select_top_theta_components_by_skew <- function(run_i, n_select = 6) {
  theta_samples <- run_i$mcmc$theta_raw_samples
  
  if (is.null(theta_samples)) {
    stop("No theta_raw_samples found in this run.")
  }
  
  q <- ncol(theta_samples)
  
  tibble(theta_index = seq_len(q)) |>
    rowwise() |>
    mutate(
      skew = safe_skewness(theta_samples[, theta_index]),
      abs_skew = abs(skew),
      excess_kurtosis = safe_excess_kurtosis(theta_samples[, theta_index]),
      abs_excess_kurtosis = abs(excess_kurtosis)
    ) |>
    ungroup() |>
    arrange(desc(abs_skew), desc(abs_excess_kurtosis)) |>
    slice_head(n = n_select)
}

plot_x_marginal_distributions <- function(run_i, components = NULL, n_select = 6, bins = 30) {
  if (is.null(components)) {
    components <- select_top_x_components_by_skew(run_i, n_select = n_select) |>
      select(agent, dimension)
  }
  
  df <- extract_x_marginal_samples(run_i, components = components)
  
  ggplot(df, aes(x = value)) +
    geom_histogram(aes(y = after_stat(density)), bins = bins, alpha = 0.5) +
    geom_density(linewidth = 0.7) +
    facet_wrap(~ component_label, scales = "free") +
    labs(
      #title = paste0("Marginal distributions of selected x-components: ",
       #              run_i$fileid, " (", run_i$model, ")"),
      x = "Sample value",
      y = "Density"
    ) +
    theme_bw(base_size = 20)
}

plot_theta_marginal_distributions <- function(run_i, theta_indices = NULL, n_select = 6, bins = 30) {
  if (is.null(theta_indices)) {
    theta_indices <- select_top_theta_components_by_skew(run_i, n_select = n_select)$theta_index
  }
  
  df <- extract_theta_marginal_samples(run_i, theta_indices = theta_indices)
  
  ggplot(df, aes(x = value)) +
    geom_histogram(aes(y = after_stat(density)), bins = bins, alpha = 0.5) +
    geom_density(linewidth = 0.7) +
    facet_wrap(~ component_label, scales = "free") +
    labs(
      #title = paste0("Marginal distributions of selected theta-components: ",
       #              run_i$fileid, " (", run_i$model, ")"),
      x = "Sample value",
      y = "Density"
    ) +
    theme_bw(base_size = 20)
}

run_scott_lid <- get_mcmc_run(
  all_latent_mcmc_runs,
  fileid = "ScottLockhard_1999b",
  model = "full"
)

plot_x_marginal_distributions(
  run_i = run_scott_lid,
  n_select = 6,
  bins = 50
)


run_pois_lid <- get_mcmc_run(
  all_latent_mcmc_runs,
  fileid = "Poisbleau_2005c",
  model = "full"
)

plot_x_marginal_distributions(
  run_i = run_pois_lid,
  n_select = 6,
  bins = 50
)

run_pois_dalid <- get_mcmc_run(
  all_latent_mcmc_runs,
  fileid = "Poisbleau_2005c",
  model = "full_theta"
)

plot_x_marginal_distributions(
  run_i = run_pois_dalid,
  n_select = 6,
  bins = 50
)
plot_theta_marginal_distributions(
  run_i = run_pois_dalid,
  n_select = 6,
  bins = 50
)

run_cui_lid <- get_mcmc_run(
  all_latent_mcmc_runs,
  fileid = "Cui_2014",
  model = "full"
)

plot_x_marginal_distributions(
  run_i = run_cui_dalid,
  n_select = 6,
  bins = 50
)

plot_theta_marginal_distributions(
  run_i = run_cui_dalid,
  n_select = 6,
  bins = 50
)

run_cui_dalid <- get_mcmc_run(
  all_latent_mcmc_runs,
  fileid = "Cui_2014",
  model = "full_theta"
)

plot_x_marginal_distributions(
  run_i = run_cui_dalid,
  n_select = 6,
  bins = 50
)

plot_theta_marginal_distributions(
  run_i = run_cui_dalid,
  n_select = 6,
  bins = 50
)

run_correa_dabt <- get_mcmc_run(
  all_latent_mcmc_runs,
  fileid = "Correa_2013a",
  model = "only_theta"
)

plot_x_marginal_distributions(
  run_i = run_correa_dabt,
  n_select = 6,
  bins = 50
)

plot_theta_marginal_distributions(
  run_i = run_correa_dabt,
  n_select = 6,
  bins = 50
)

run_adcock_lid <- get_mcmc_run(
  all_latent_mcmc_runs,
  fileid = "Adcock_2015a",
  model = "full"
)

plot_x_marginal_distributions(
  run_i = run_adcock_lid,
  n_select = 6,
  bins = 50
)

run_adcock_dalid <- get_mcmc_run(
  all_latent_mcmc_runs,
  fileid = "Adcock_2015a",
  model = "full_theta"
)

plot_x_marginal_distributions(
  run_i = run_adcock_dalid,
  n_select = 6,
  bins = 50
)

plot_theta_marginal_distributions(
  run_i = run_adcock_dalid,
  n_select = 6,
  bins = 50
)

run_general <- get_mcmc_run(
  all_latent_mcmc_runs,
  fileid = "Mwamende_2009a",
  model = "full"
)

plot_x_marginal_distributions(
  run_i = run_general,
  n_select = 6,
  bins = 50
)

plot_theta_marginal_distributions(
  run_i = run_general,
  n_select = 6,
  bins = 50
)

#######################################
# Standardizing marginals
#######################################

make_standardized_x_df <- function(all_latent_mcmc_runs,
                                   fileid,
                                   model,
                                   thin_by = 1) {
  run_i <- get_mcmc_run(
    all_latent_mcmc_runs = all_latent_mcmc_runs,
    fileid = fileid,
    model = model
  )
  
  x_samples <- run_i$mcmc$x_samples
  
  if (is.null(x_samples)) {
    x_samples <- run_i$mcmc$samples
  }
  
  if (is.null(x_samples)) {
    stop("No x_samples found in run.")
  }
  
  if (length(dim(x_samples)) != 3) {
    stop("x_samples must be a 3D array: draw x agent x dimension.")
  }
  
  dims <- dim(x_samples)
  n_draws <- dims[1]
  n_agents <- dims[2]
  d <- dims[3]
  
  keep_draws <- seq(1, n_draws, by = thin_by)
  
  map_dfr(seq_len(n_agents), function(agent_i) {
    map_dfr(seq_len(d), function(dim_i) {
      chain <- x_samples[keep_draws, agent_i, dim_i]
      
      chain_mean <- mean(chain, na.rm = TRUE)
      chain_sd <- sd(chain, na.rm = TRUE)
      
      if (!is.finite(chain_sd) || chain_sd <= 0) {
        z <- rep(NA_real_, length(chain))
      } else {
        z <- (chain - chain_mean) / chain_sd
      }
      
      tibble(
        fileid = fileid,
        model = model,
        draw_original = keep_draws,
        draw = seq_along(keep_draws),
        agent = agent_i,
        dimension = dim_i,
        x = chain,
        x_mean = chain_mean,
        x_sd = chain_sd,
        z = z
      )
    })
  }) |>
    filter(is.finite(z))
}

plot_standardized_x_marginals <- function(all_latent_mcmc_runs,
                                          fileid,
                                          model,
                                          thin_by = 1,
                                          by_dimension = TRUE,
                                          bins = 60,
                                          show_histogram = TRUE) {
  plot_df <- make_standardized_x_df(
    all_latent_mcmc_runs = all_latent_mcmc_runs,
    fileid = fileid,
    model = model,
    thin_by = thin_by
  )
  
  model_label <- recode(
    model,
    "null" = "LBT",
    "full" = "LID",
    "only_theta" = "DABT",
    "full_theta" = "DALID",
    .default = model
  )
  
  p <- ggplot(plot_df, aes(x = z))
  
  if (show_histogram) {
    p <- p +
      geom_histogram(
        aes(y = after_stat(density)),
        bins = bins,
        alpha = 0.35
      )
  }
  
  p <- p +
    geom_density(linewidth = 0.8) +
    stat_function(
      fun = dnorm,
      linetype = "dashed",
      linewidth = 0.7
    ) +
    labs(
      title = paste0(
        "Aggregated standardized marginal distributions of x: ",
        fileid, " (", model_label, ")"
      ),
      subtitle = paste0(
        "Each marginal coordinate x[i,k] is centered and scaled before aggregation. ",
        "Dashed curve is N(0,1)."
      ),
      x = expression((x[ik] - bar(x)[ik]) / widehat(sd)(x[ik])),
      y = "Density"
    ) +
    theme_minimal()
  
  if (by_dimension) {
    p <- p + facet_wrap(~ dimension, scales = "free_y")
  }
  
  p
}

plot_standardized_x_marginals(
  all_latent_mcmc_runs = all_latent_mcmc_runs,
  fileid = "ScottLockhard_1999b",
  model = "full",
  thin_by = 5
)

plot_standardized_x_marginals(
  all_latent_mcmc_runs = all_latent_mcmc_runs,
  fileid = "Poisbleau_2005c",
  model = "full",
  thin_by = 5
)

plot_standardized_x_marginals(
  all_latent_mcmc_runs = all_latent_mcmc_runs,
  fileid = "Cui_2014",
  model = "full_theta",
  thin_by = 5
)
plot_standardized_x_marginals(
  all_latent_mcmc_runs = all_latent_mcmc_runs,
  fileid = "Alados_1992b",
  model = "full",
  thin_by = 5
)

mcmc_x_samples_to_df_thin <- function(mcmc, fileid, model, thin_by = 5) {
  x_samples <- mcmc$x_samples
  
  if (is.null(x_samples)) {
    x_samples <- mcmc$samples
  }
  
  if (length(dim(x_samples)) != 3) {
    stop("x_samples must be a 3D array: draw x agent x dimension.")
  }
  
  dims <- dim(x_samples)
  
  keep_draws <- seq(1, dims[1], by = thin_by)
  x_samples_thin <- x_samples[keep_draws, , , drop = FALSE]
  dims_thin <- dim(x_samples_thin)
  
  expand.grid(
    draw_thinned = seq_len(dims_thin[1]),
    agent = seq_len(dims_thin[2]),
    dimension = seq_len(dims_thin[3])
  ) |>
    dplyr::as_tibble() |>
    dplyr::mutate(
      draw_original = rep(keep_draws, each = dims_thin[2] * dims_thin[3]),
      x = as.vector(x_samples_thin),
      fileid = fileid,
      model = model
    ) |>
    dplyr::select(
      fileid,
      model,
      draw_original,
      draw_thinned,
      agent,
      dimension,
      x
    )
}

mcmc_theta_samples_to_df_thin <- function(mcmc, data, fileid, model, thin_by = 5) {
  theta_samples <- mcmc$theta_raw_samples
  
  if (is.null(theta_samples)) {
    return(tibble::tibble())
  }
  
  keep_draws <- seq(1, nrow(theta_samples), by = thin_by)
  theta_samples_thin <- theta_samples[keep_draws, , drop = FALSE]
  
  q <- ncol(theta_samples_thin)
  
  dyad_lookup <- tibble::tibble(
    theta_index = seq_len(q),
    i = data$i,
    j = data$j,
    s = data$s,
    z = data$z
  )
  
  expand.grid(
    draw_thinned = seq_len(nrow(theta_samples_thin)),
    theta_index = seq_len(q)
  ) |>
    dplyr::as_tibble() |>
    dplyr::mutate(
      draw_original = rep(keep_draws, each = q),
      theta_raw = as.vector(theta_samples_thin),
      fileid = fileid,
      model = model
    ) |>
    dplyr::left_join(dyad_lookup, by = "theta_index") |>
    dplyr::select(
      fileid,
      model,
      draw_original,
      draw_thinned,
      theta_index,
      i,
      j,
      z,
      s,
      theta_raw
    )
}

all_diagnostics <- purrr::map_dfr(
  all_latent_mcmc_runs,
  ~ purrr::map_dfr(.x, "diagnostics"),
  .id = "model_group"
)
all_diagnostics
thin_by <- 5

all_x_samples <- purrr::map_dfr(
  all_latent_mcmc_runs,
  function(model_runs) {
    purrr::map_dfr(model_runs, function(run_i) {
      if (is.null(run_i$mcmc)) {
        return(tibble::tibble())
      }
      
      mcmc_x_samples_to_df_thin(
        mcmc = run_i$mcmc,
        fileid = run_i$fileid,
        model = run_i$model,
        thin_by = thin_by
      )
    })
  },
  .id = "model_group"
)

all_theta_samples <- purrr::map_dfr(
  all_latent_mcmc_runs,
  function(model_runs) {
    purrr::map_dfr(model_runs, function(run_i) {
      if (is.null(run_i$mcmc) || is.null(run_i$data)) {
        return(tibble::tibble())
      }
      
      mcmc_theta_samples_to_df_thin(
        mcmc = run_i$mcmc,
        data = run_i$data,
        fileid = run_i$fileid,
        model = run_i$model,
        thin_by = thin_by
      )
    })
  },
  .id = "model_group"
)

saveRDS(all_latent_mcmc_runs, file.path(out_dir, "all_latent_mcmc_runs.rds"))
saveRDS(all_diagnostics, file.path(out_dir, "all_latent_mcmc_diagnostics.rds"))
saveRDS(all_x_samples, file.path(out_dir, "all_latent_mcmc_x_samples.rds"))
saveRDS(all_theta_samples, file.path(out_dir, "all_latent_mcmc_theta_samples.rds"))

readr::write_csv(all_diagnostics, file.path(out_dir, "all_latent_mcmc_diagnostics.csv"))

#######################################
# Results 
#######################################

library(stringr)
library(readr)
library(gt)

# ------------------------------------------------------------
# Table 1: Overall run status by model
# ------------------------------------------------------------

table_mcmc_status <- all_diagnostics |>
  group_by(model) |>
  summarise(
    n_runs = n(),
    n_ok = sum(status == "ok", na.rm = TRUE),
    n_error = sum(status == "error", na.rm = TRUE),
    prop_ok = n_ok / n_runs,
    median_acceptance = median(acceptance_rate, na.rm = TRUE),
    min_acceptance = min(acceptance_rate, na.rm = TRUE),
    max_acceptance = max(acceptance_rate, na.rm = TRUE),
    median_latent_dim = median(n_total_latent_dim, na.rm = TRUE),
    max_latent_dim = max(n_total_latent_dim, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(
    across(
      c(prop_ok, median_acceptance, min_acceptance, max_acceptance),
      ~ round(.x, 3)
    )
  )

table_mcmc_status
library(stringr)
library(readr)
library(gt)

# ------------------------------------------------------------
# Table 1: Overall run status by model
# ------------------------------------------------------------

table_mcmc_status <- all_diagnostics |>
  group_by(model) |>
  summarise(
    n_runs = n(),
    n_ok = sum(status == "ok", na.rm = TRUE),
    n_error = sum(status == "error", na.rm = TRUE),
    prop_ok = n_ok / n_runs,
    median_acceptance = median(acceptance_rate, na.rm = TRUE),
    min_acceptance = min(acceptance_rate, na.rm = TRUE),
    max_acceptance = max(acceptance_rate, na.rm = TRUE),
    median_latent_dim = median(n_total_latent_dim, na.rm = TRUE),
    max_latent_dim = max(n_total_latent_dim, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(
    across(
      c(prop_ok, median_acceptance, min_acceptance, max_acceptance),
      ~ round(.x, 3)
    )
  )

table_mcmc_status
gt_table_mcmc_status <- table_mcmc_status |>
  gt() |>
  tab_header(
    title = "Overall MCMC run status by model"
  ) |>
  fmt_number(
    columns = c(prop_ok, median_acceptance, min_acceptance, max_acceptance),
    decimals = 3
  )

gt_table_mcmc_status

# ------------------------------------------------------------
# Table 2: Detailed diagnostics per dataset and model
# ------------------------------------------------------------

table_mcmc_diagnostics <- all_diagnostics |>
  mutate(
    acceptance_flag = case_when(
      is.na(acceptance_rate) ~ "missing",
      acceptance_rate < 0.05 ~ "very low",
      acceptance_rate < 0.15 ~ "low",
      acceptance_rate <= 0.35 ~ "reasonable",
      acceptance_rate <= 0.60 ~ "high",
      TRUE ~ "very high"
    ),
    opt_flag = case_when(
      is.na(opt_convergence) ~ "missing",
      opt_convergence == 0 ~ "ok",
      TRUE ~ "check"
    ),
    hessian_flag = case_when(
      is.na(model_hess_status) ~ "missing",
      model_hess_status %in% c("positive definite", "pd", "PD") ~ "ok",
      TRUE ~ as.character(model_hess_status)
    )
  ) |>
  select(
    fileid,
    model,
    sampled_latent_state,
    status,
    n_agents,
    n_observed_dyads,
    n_interactions,
    n_total_latent_dim,
    n_x_dim,
    n_theta_raw_dim,
    acceptance_rate,
    acceptance_flag,
    proposal_sd_x_final,
    proposal_sd_theta_final,
    mode_nll,
    final_nll,
    mean_nll_after_burn,
    sd_nll_after_burn,
    opt_convergence,
    opt_flag,
    model_hess_min_eig,
    model_hess_abs_cond,
    model_hess_status,
    log_r_hat,
    a_hat,
    sigma_theta_hat,
    error_message
  ) |>
  arrange(model, fileid)

table_mcmc_diagnostics
# ------------------------------------------------------------
# Table 3: Compact diagnostics table for thesis text
# ------------------------------------------------------------

table_mcmc_compact <- all_diagnostics |>
  transmute(
    Dataset = fileid,
    Model = model,
    `Latent state` = sampled_latent_state,
    Status = status,
    `n agents` = n_agents,
    `Observed dyads` = n_observed_dyads,
    `Latent dim.` = n_total_latent_dim,
    `Acceptance rate` = acceptance_rate,
    `Final proposal sd, x` = proposal_sd_x_final,
    `Final proposal sd, theta` = proposal_sd_theta_final,
    `Mean NLL after burn-in` = mean_nll_after_burn,
    `SD NLL after burn-in` = sd_nll_after_burn,
    `Hessian min. eig.` = model_hess_min_eig,
    `Hessian cond.` = model_hess_abs_cond,
    `Hessian status` = model_hess_status
  ) |>
  arrange(Model, Dataset)

gt_table_mcmc_compact <- table_mcmc_compact |>
  gt() |>
  tab_header(
    title = "Conditional MCMC diagnostics for latent variables"
  ) |>
  fmt_number(
    columns = c(
      `Acceptance rate`,
      `Final proposal sd, x`,
      `Final proposal sd, theta`,
      `Mean NLL after burn-in`,
      `SD NLL after burn-in`,
      `Hessian min. eig.`,
      `Hessian cond.`
    ),
    decimals = 3
  )

gt_table_mcmc_compact

# ------------------------------------------------------------
# Table 4: Aggregated x summaries by dataset/model/dimension
# ------------------------------------------------------------

table_x_by_dataset_dimension <- all_x_samples |>
  group_by(fileid, model, dimension) |>
  summarise(
    n_samples_total = n(),
    n_agents = n_distinct(agent),
    mean_x = mean(x, na.rm = TRUE),
    sd_x = sd(x, na.rm = TRUE),
    median_x = median(x, na.rm = TRUE),
    q025_x = quantile(x, 0.025, na.rm = TRUE),
    q975_x = quantile(x, 0.975, na.rm = TRUE),
    min_x = min(x, na.rm = TRUE),
    max_x = max(x, na.rm = TRUE),
    .groups = "drop"
  ) |>
  arrange(model, fileid, dimension)

table_x_by_dataset_dimension

# ------------------------------------------------------------
# Table 5: x summaries by agent and dimension
# ------------------------------------------------------------

table_x_by_agent <- all_x_samples |>
  group_by(fileid, model, agent, dimension) |>
  summarise(
    n_draws = n(),
    mean_x = mean(x, na.rm = TRUE),
    sd_x = sd(x, na.rm = TRUE),
    median_x = median(x, na.rm = TRUE),
    q025_x = quantile(x, 0.025, na.rm = TRUE),
    q975_x = quantile(x, 0.975, na.rm = TRUE),
    ci_width_x = q975_x - q025_x,
    .groups = "drop"
  ) |>
  arrange(model, fileid, dimension, desc(sd_x))

table_x_by_agent

table_x_most_uncertain <- table_x_by_agent |>
  group_by(model) |>
  slice_max(sd_x, n = 20, with_ties = FALSE) |>
  ungroup() |>
  arrange(model, desc(sd_x))

table_x_most_uncertain

# ------------------------------------------------------------
# Helper functions for skewness and kurtosis
# ------------------------------------------------------------

sample_skewness <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 3) return(NA_real_)
  m <- mean(x)
  s <- sd(x)
  if (s == 0) return(NA_real_)
  mean((x - m)^3) / s^3
}

sample_kurtosis <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 4) return(NA_real_)
  m <- mean(x)
  s <- sd(x)
  if (s == 0) return(NA_real_)
  mean((x - m)^4) / s^4
}

# ------------------------------------------------------------
# Table 6: Shape diagnostics for x
# ------------------------------------------------------------

table_x_shape_by_agent <- all_x_samples |>
  group_by(fileid, model, agent, dimension) |>
  summarise(
    n_draws = n(),
    mean_x = mean(x, na.rm = TRUE),
    sd_x = sd(x, na.rm = TRUE),
    skew_x = sample_skewness(x),
    kurtosis_x = sample_kurtosis(x),
    excess_kurtosis_x = kurtosis_x - 3,
    q025_x = quantile(x, 0.025, na.rm = TRUE),
    q975_x = quantile(x, 0.975, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(
    nonnormality_score =
      abs(skew_x) +
      abs(excess_kurtosis_x) / 3
  ) |>
  arrange(model, fileid, desc(nonnormality_score))

table_x_shape_by_agent

table_x_most_nonnormal <- table_x_shape_by_agent |>
  group_by(model) |>
  slice_max(nonnormality_score, n = 20, with_ties = FALSE) |>
  ungroup() |>
  arrange(model, desc(nonnormality_score))

table_x_most_nonnormal

# ------------------------------------------------------------
# Table 7: Aggregated x normality diagnostics by dataset/model
# ------------------------------------------------------------

table_x_shape_by_dataset <- table_x_shape_by_agent |>
  group_by(fileid, model) |>
  summarise(
    n_coordinates = n(),
    
    median_abs_skew = median(abs(skew_x), na.rm = TRUE),
    max_abs_skew = max(abs(skew_x), na.rm = TRUE),
    
    median_abs_excess_kurtosis = median(abs(excess_kurtosis_x), na.rm = TRUE),
    max_abs_excess_kurtosis = max(abs(excess_kurtosis_x), na.rm = TRUE),
    
    prop_abs_skew_gt_1 = mean(abs(skew_x) > 1, na.rm = TRUE),
    prop_abs_excess_kurtosis_gt_2 = mean(abs(excess_kurtosis_x) > 2, na.rm = TRUE),
    
    median_sd_x = median(sd_x, na.rm = TRUE),
    max_sd_x = max(sd_x, na.rm = TRUE),
    
    median_ci_width_x = median(q975_x - q025_x, na.rm = TRUE),
    max_ci_width_x = max(q975_x - q025_x, na.rm = TRUE),
    
    .groups = "drop"
  ) |>
  arrange(desc(median_abs_skew + median_abs_excess_kurtosis))

table_x_shape_by_dataset

# ------------------------------------------------------------
# Table 8: theta_raw summaries by dyad
# ------------------------------------------------------------

table_theta_by_dyad <- all_theta_samples |>
  group_by(fileid, model, theta_index, i, j) |>
  summarise(
    n_draws = n(),
    z = first(z),
    s = first(s),
    mean_theta_raw = mean(theta_raw, na.rm = TRUE),
    sd_theta_raw = sd(theta_raw, na.rm = TRUE),
    median_theta_raw = median(theta_raw, na.rm = TRUE),
    q025_theta_raw = quantile(theta_raw, 0.025, na.rm = TRUE),
    q975_theta_raw = quantile(theta_raw, 0.975, na.rm = TRUE),
    skew_theta_raw = sample_skewness(theta_raw),
    kurtosis_theta_raw = sample_kurtosis(theta_raw),
    excess_kurtosis_theta_raw = kurtosis_theta_raw - 3,
    .groups = "drop"
  ) |>
  mutate(
    theta_raw_ci_width = q975_theta_raw - q025_theta_raw,
    theta_raw_abs_mean = abs(mean_theta_raw),
    theta_raw_nonnormality_score =
      abs(skew_theta_raw) +
      abs(excess_kurtosis_theta_raw) / 3
  ) |>
  arrange(model, fileid, desc(theta_raw_abs_mean))

table_theta_by_dyad

table_theta_largest_effects <- table_theta_by_dyad |>
  group_by(model) |>
  slice_max(theta_raw_abs_mean, n = 30, with_ties = FALSE) |>
  ungroup() |>
  arrange(desc(theta_raw_abs_mean))

table_theta_largest_effects

table_theta_most_uncertain <- table_theta_by_dyad |>
  group_by(model) |>
  slice_max(sd_theta_raw, n = 30, with_ties = FALSE) |>
  ungroup() |>
  arrange(model, desc(sd_theta_raw))

table_theta_most_uncertain

# ------------------------------------------------------------
# Table 9: theta_raw summaries by dataset/model
# ------------------------------------------------------------

table_theta_by_dataset <- table_theta_by_dyad |>
  group_by(fileid, model) |>
  summarise(
    n_theta_effects = n(),
    
    median_abs_mean_theta_raw = median(abs(mean_theta_raw), na.rm = TRUE),
    max_abs_mean_theta_raw = max(abs(mean_theta_raw), na.rm = TRUE),
    
    median_sd_theta_raw = median(sd_theta_raw, na.rm = TRUE),
    max_sd_theta_raw = max(sd_theta_raw, na.rm = TRUE),
    
    median_abs_skew_theta_raw = median(abs(skew_theta_raw), na.rm = TRUE),
    max_abs_skew_theta_raw = max(abs(skew_theta_raw), na.rm = TRUE),
    
    median_abs_excess_kurtosis_theta_raw =
      median(abs(excess_kurtosis_theta_raw), na.rm = TRUE),
    max_abs_excess_kurtosis_theta_raw =
      max(abs(excess_kurtosis_theta_raw), na.rm = TRUE),
    
    prop_abs_skew_theta_gt_1 =
      mean(abs(skew_theta_raw) > 1, na.rm = TRUE),
    prop_abs_excess_kurtosis_theta_gt_2 =
      mean(abs(excess_kurtosis_theta_raw) > 2, na.rm = TRUE),
    
    .groups = "drop"
  ) |>
  arrange(model, desc(median_abs_mean_theta_raw))

table_theta_by_dataset

# ------------------------------------------------------------
# Table 10: Combined model diagnostics and MCMC shape diagnostics
# ------------------------------------------------------------

table_combined_diagnostics <- all_diagnostics |>
  select(
    fileid,
    model,
    status,
    sampled_latent_state,
    n_agents,
    n_observed_dyads,
    n_interactions,
    n_total_latent_dim,
    acceptance_rate,
    proposal_sd_x_final,
    proposal_sd_theta_final,
    mean_nll_after_burn,
    sd_nll_after_burn,
    model_hess_min_eig,
    model_hess_abs_cond,
    model_hess_status,
    log_r_hat,
    a_hat,
    sigma_theta_hat
  ) |>
  left_join(
    table_x_shape_by_dataset |>
      select(
        fileid,
        model,
        n_coordinates,
        median_abs_skew_x = median_abs_skew,
        max_abs_skew_x = max_abs_skew,
        median_abs_excess_kurtosis_x = median_abs_excess_kurtosis,
        max_abs_excess_kurtosis_x = max_abs_excess_kurtosis,
        prop_abs_skew_gt_1,
        prop_abs_excess_kurtosis_gt_2,
        median_sd_x,
        max_sd_x
      ),
    by = c("fileid", "model")
  ) |>
  left_join(
    table_theta_by_dataset |>
      select(
        fileid,
        model,
        n_theta_effects,
        median_abs_mean_theta_raw,
        max_abs_mean_theta_raw,
        median_sd_theta_raw,
        max_sd_theta_raw,
        median_abs_skew_theta_raw,
        median_abs_excess_kurtosis_theta_raw,
        prop_abs_skew_theta_gt_1,
        prop_abs_excess_kurtosis_theta_gt_2
      ),
    by = c("fileid", "model")
  ) |>
  arrange(model, fileid)

table_combined_diagnostics

# ------------------------------------------------------------
# Table 11: Most problematic MCMC runs
# ------------------------------------------------------------

table_problem_ranking <- table_combined_diagnostics |>
  mutate(
    acceptance_penalty = case_when(
      is.na(acceptance_rate) ~ 5,
      acceptance_rate < 0.05 ~ 4,
      acceptance_rate < 0.15 ~ 2,
      acceptance_rate > 0.60 ~ 2,
      acceptance_rate > 0.35 ~ 1,
      TRUE ~ 0
    ),
    hessian_penalty = case_when(
      is.na(model_hess_min_eig) ~ 1,
      model_hess_min_eig <= 0 ~ 3,
      TRUE ~ 0
    ),
    x_shape_penalty =
      coalesce(median_abs_skew_x, 0) +
      coalesce(median_abs_excess_kurtosis_x, 0) / 3 +
      2 * coalesce(prop_abs_skew_gt_1, 0) +
      2 * coalesce(prop_abs_excess_kurtosis_gt_2, 0),
    theta_shape_penalty =
      coalesce(median_abs_skew_theta_raw, 0) +
      coalesce(median_abs_excess_kurtosis_theta_raw, 0) / 3 +
      2 * coalesce(prop_abs_skew_theta_gt_1, 0) +
      2 * coalesce(prop_abs_excess_kurtosis_theta_gt_2, 0),
    problem_score =
      acceptance_penalty +
      hessian_penalty +
      x_shape_penalty +
      theta_shape_penalty
  ) |>
  arrange(model, desc(problem_score)) |>
  select(
    fileid,
    model,
    status,
    sampled_latent_state,
    n_total_latent_dim,
    acceptance_rate,
    model_hess_min_eig,
    model_hess_abs_cond,
    median_abs_skew_x,
    median_abs_excess_kurtosis_x,
    prop_abs_skew_gt_1,
    prop_abs_excess_kurtosis_gt_2,
    median_abs_skew_theta_raw,
    median_abs_excess_kurtosis_theta_raw,
    problem_score
  )

table_problem_ranking

library(knitr)

latex_mcmc_table <- table_problem_ranking |>
  mutate(
    acceptance_rate = round(acceptance_rate, 3),
    model_hess_min_eig = signif(model_hess_min_eig, 3),
    model_hess_abs_cond = signif(model_hess_abs_cond, 3),
    median_abs_skew_x = round(median_abs_skew_x, 3),
    median_abs_excess_kurtosis_x = round(median_abs_excess_kurtosis_x, 3),
    problem_score = round(problem_score, 3)
  ) |>
  select(
    fileid,
    model,
    n_total_latent_dim,
    acceptance_rate,
    median_abs_skew_x,
    median_abs_excess_kurtosis_x,
    problem_score
  )

knitr::kable(
  latex_mcmc_table,
  format = "latex",
  booktabs = TRUE,
  caption = "Summary of conditional MCMC diagnostics for selected fitted models.",
  digits = 3
)


# ------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------

sample_skewness <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 3) return(NA_real_)
  m <- mean(x)
  s <- sd(x)
  if (!is.finite(s) || s == 0) return(NA_real_)
  mean((x - m)^3) / s^3
}

sample_kurtosis <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 4) return(NA_real_)
  m <- mean(x)
  s <- sd(x)
  if (!is.finite(s) || s == 0) return(NA_real_)
  mean((x - m)^4) / s^4
}

safe_iat <- function(x, max_lag = 1000) {
  x <- x[is.finite(x)]
  n <- length(x)
  
  if (n < 10) return(NA_real_)
  if (sd(x) == 0) return(NA_real_)
  
  max_lag <- min(max_lag, n - 1)
  
  acf_vals <- as.numeric(stats::acf(
    x,
    lag.max = max_lag,
    plot = FALSE,
    na.action = na.pass
  )$acf)[-1]
  
  # Geyer-like simple positive sequence truncation
  positive_acf <- acf_vals[acf_vals > 0]
  
  if (length(positive_acf) == 0) {
    return(1)
  }
  
  1 + 2 * sum(positive_acf)
}

safe_ess <- function(x, max_lag = 1000) {
  iat <- safe_iat(x, max_lag = max_lag)
  if (!is.finite(iat) || iat <= 0) return(NA_real_)
  length(x[is.finite(x)]) / iat
}

# ------------------------------------------------------------
# 1. Shape and trace diagnostics for x by coordinate
# ------------------------------------------------------------

x_coordinate_diagnostics <- all_x_samples |>
  group_by(fileid, model, agent, dimension) |>
  summarise(
    n_draws = n(),
    mean_x = mean(x, na.rm = TRUE),
    sd_x = sd(x, na.rm = TRUE),
    skew_x = sample_skewness(x),
    excess_kurtosis_x = sample_kurtosis(x) - 3,
    ess_x = safe_ess(x),
    iat_x = safe_iat(x),
    .groups = "drop"
  ) |>
  mutate(
    abs_skew_x = abs(skew_x),
    abs_excess_kurtosis_x = abs(excess_kurtosis_x)
  )

# ------------------------------------------------------------
# 2. Aggregate x diagnostics by dataset and model
# ------------------------------------------------------------

x_summary_by_run <- x_coordinate_diagnostics |>
  group_by(fileid, model) |>
  summarise(
    median_sd_x = median(sd_x, na.rm = TRUE),
    max_sd_x = max(sd_x, na.rm = TRUE),
    
    median_abs_skew_x = median(abs_skew_x, na.rm = TRUE),
    max_abs_skew_x = max(abs_skew_x, na.rm = TRUE),
    
    median_abs_excess_kurtosis_x = median(abs_excess_kurtosis_x, na.rm = TRUE),
    max_abs_excess_kurtosis_x = max(abs_excess_kurtosis_x, na.rm = TRUE),
    
    median_ess_x = median(ess_x, na.rm = TRUE),
    min_ess_x = min(ess_x, na.rm = TRUE),
    
    median_iat_x = median(iat_x, na.rm = TRUE),
    max_iat_x = max(iat_x, na.rm = TRUE),
    
    .groups = "drop"
  )

# ------------------------------------------------------------
# 3. MCMC run-level diagnostics from all_diagnostics
# ------------------------------------------------------------

run_summary <- all_diagnostics |>
  group_by(fileid, model) |>
  summarise(
    status = first(status),
    sampled_latent_state = first(sampled_latent_state),
    
    n_agents = first(n_agents),
    n_observed_dyads = first(n_observed_dyads),
    n_interactions = first(n_interactions),
    
    n_total_latent_dim = first(n_total_latent_dim),
    n_x_dim = first(n_x_dim),
    n_theta_raw_dim = first(n_theta_raw_dim),
    
    acceptance_rate = mean(acceptance_rate, na.rm = TRUE),
    proposal_sd_x_final = first(proposal_sd_x_final),
    proposal_sd_theta_final = first(proposal_sd_theta_final),
    
    mean_nll_after_burn = first(mean_nll_after_burn),
    sd_nll_after_burn = first(sd_nll_after_burn),
    final_nll = first(final_nll),
    mode_nll = first(mode_nll),
    
    opt_convergence = first(opt_convergence),
    opt_message = first(opt_message),
    
    model_hess_min_eig = first(model_hess_min_eig),
    model_hess_abs_cond = first(model_hess_abs_cond),
    model_hess_status = first(model_hess_status),
    
    error_message = first(error_message),
    .groups = "drop"
  )

# ------------------------------------------------------------
# 4. Main summary table: x diagnostics
# ------------------------------------------------------------

mcmc_summary_table <- run_summary |>
  left_join(x_summary_by_run, by = c("fileid", "model")) |>
  mutate(
    acceptance_flag = case_when(
      is.na(acceptance_rate) ~ "missing",
      acceptance_rate < 0.05 ~ "very low",
      acceptance_rate < 0.15 ~ "low",
      acceptance_rate <= 0.35 ~ "reasonable",
      acceptance_rate <= 0.60 ~ "high",
      TRUE ~ "very high"
    ),
    trace_flag = case_when(
      is.na(median_ess_x) ~ "missing",
      median_ess_x < 50 ~ "poor",
      median_ess_x < 200 ~ "moderate",
      TRUE ~ "ok"
    ),
    shape_flag_x = case_when(
      is.na(median_abs_skew_x) | is.na(median_abs_excess_kurtosis_x) ~ "missing",
      median_abs_skew_x > 1 | median_abs_excess_kurtosis_x > 2 ~ "strong deviation",
      median_abs_skew_x > 0.5 | median_abs_excess_kurtosis_x > 1 ~ "moderate deviation",
      TRUE ~ "mild"
    )
  ) |>
  select(
    fileid,
    model,
    status,
    sampled_latent_state,
    n_agents,
    n_observed_dyads,
    n_total_latent_dim,
    n_x_dim,
    n_theta_raw_dim,
    acceptance_rate,
    acceptance_flag,
    median_abs_skew_x,
    median_abs_excess_kurtosis_x,
    median_sd_x,
    max_sd_x,
    median_ess_x,
    min_ess_x,
    median_iat_x,
    max_iat_x,
    trace_flag,
    shape_flag_x,
    mean_nll_after_burn,
    sd_nll_after_burn,
    model_hess_min_eig,
    model_hess_abs_cond,
    model_hess_status,
    error_message
  ) |>
  arrange(model, fileid)

mcmc_summary_table

mcmc_summary_compact_x <- mcmc_summary_table |>
  mutate(
    Model = recode(
      model,
      "null" = "LBT",
      "full" = "LID",
      "only_theta" = "DABT",
      "full_theta" = "DALID",
      .default = model
    ),
    flag = case_when(
      status == "error" ~ "error",
      acceptance_rate < 0.05 ~ "very low acc.",
      acceptance_rate < 0.15 ~ "low acc.",
      acceptance_rate > 0.60 ~ "very high acc.",
      acceptance_rate > 0.35 ~ "high acc.",
      median_ess_x < 50 ~ "low ESS",
      median_abs_skew_x > 1 ~ "high skew",
      median_abs_excess_kurtosis_x > 2 ~ "heavy tails",
      TRUE ~ "ok"
    )
  ) |>
  transmute(
    Dataset = fileid,
    Model,
    `Latent dim.` = n_total_latent_dim,
    `Acc.` = acceptance_rate,
    `Med. |skew|` = median_abs_skew_x,
    `Med. |ex. kurt.|` = median_abs_excess_kurtosis_x,
    `Med. sd` = median_sd_x,
    `Med. ESS` = median_ess_x,
    `Med. IAT` = median_iat_x,
    `SD NLL` = sd_nll_after_burn,
    Flag = flag
  ) |>
  arrange(Model, Dataset)

mcmc_summary_compact_x

mcmc_summary_compact_x_latex <- mcmc_summary_compact_x |>
  mutate(
    `Acc.` = round(`Acc.`, 3),
    `Med. |skew|` = round(`Med. |skew|`, 3),
    `Med. |ex. kurt.|` = round(`Med. |ex. kurt.|`, 3),
    `Med. sd` = round(`Med. sd`, 3),
    `Med. ESS` = round(`Med. ESS`, 1),
    `Med. IAT` = round(`Med. IAT`, 1),
    `SD NLL` = round(`SD NLL`, 3)
  )

knitr::kable(
  mcmc_summary_compact_x_latex,
  format = "latex",
  booktabs = TRUE,
  caption = "Compact summary of conditional MCMC diagnostics for the latent coordinates.",
  label = "tab:mcmc_summary_compact_x"
)


#################################
# Normality test D???Agostino???Pearson-test
#################################

dagostino_pearson_test <- function(x) {
  x <- x[is.finite(x)]
  n <- length(x)
  
  if (n < 20) {
    return(tibble::tibble(
      n = n,
      skew = NA_real_,
      kurtosis = NA_real_,
      excess_kurtosis = NA_real_,
      z_skew = NA_real_,
      z_kurtosis = NA_real_,
      k2 = NA_real_,
      p_value = NA_real_
    ))
  }
  
  m <- mean(x)
  xc <- x - m
  
  m2 <- mean(xc^2)
  m3 <- mean(xc^3)
  m4 <- mean(xc^4)
  
  if (!is.finite(m2) || m2 <= 0) {
    return(tibble::tibble(
      n = n,
      skew = NA_real_,
      kurtosis = NA_real_,
      excess_kurtosis = NA_real_,
      z_skew = NA_real_,
      z_kurtosis = NA_real_,
      k2 = NA_real_,
      p_value = NA_real_
    ))
  }
  
  # Moment skewness and Pearson kurtosis
  skew <- m3 / m2^(3 / 2)
  kurtosis <- m4 / m2^2
  excess_kurtosis <- kurtosis - 3
  
  # ----------------------------------------------------------
  # D'Agostino skewness transformation
  # ----------------------------------------------------------
  
  y <- skew * sqrt(((n + 1) * (n + 3)) / (6 * (n - 2)))
  
  beta2 <- (3 * (n^2 + 27 * n - 70) * (n + 1) * (n + 3)) /
    ((n - 2) * (n + 5) * (n + 7) * (n + 9))
  
  w2 <- -1 + sqrt(2 * (beta2 - 1))
  delta <- 1 / sqrt(log(sqrt(w2)))
  alpha <- sqrt(2 / (w2 - 1))
  
  z_skew <- delta * log(y / alpha + sqrt((y / alpha)^2 + 1))
  
  # ----------------------------------------------------------
  # Anscombe-Glynn kurtosis transformation
  # ----------------------------------------------------------
  
  expected_kurtosis <- 3 * (n - 1) / (n + 1)
  
  var_kurtosis <- (24 * n * (n - 2) * (n - 3)) /
    ((n + 1)^2 * (n + 3) * (n + 5))
  
  x_kurt <- (kurtosis - expected_kurtosis) / sqrt(var_kurtosis)
  
  sqrt_beta1 <- (
    6 * (n^2 - 5 * n + 2) /
      ((n + 7) * (n + 9))
  ) *
    sqrt(
      (6 * (n + 3) * (n + 5)) /
        (n * (n - 2) * (n - 3))
    )
  
  A <- 6 + (8 / sqrt_beta1) *
    (2 / sqrt_beta1 + sqrt(1 + 4 / sqrt_beta1^2))
  
  term1 <- 1 - 2 / (9 * A)
  term2 <- (1 - 2 / A) / (1 + x_kurt * sqrt(2 / (A - 4)))
  
  z_kurtosis <- (term1 - sign(term2) * abs(term2)^(1 / 3)) /
    sqrt(2 / (9 * A))
  
  k2 <- z_skew^2 + z_kurtosis^2
  p_value <- stats::pchisq(k2, df = 2, lower.tail = FALSE)
  
  tibble::tibble(
    n = n,
    skew = skew,
    kurtosis = kurtosis,
    excess_kurtosis = excess_kurtosis,
    z_skew = z_skew,
    z_kurtosis = z_kurtosis,
    k2 = k2,
    p_value = p_value
  )
}

dagostino_x_by_coordinate <- all_x_samples |>
  group_by(fileid, model, agent, dimension) |>
  summarise(
    dagostino_pearson_test(x),
    .groups = "drop"
  ) |>
  group_by(fileid, model) |>
  mutate(
    p_adj_BH = p.adjust(p_value, method = "BH"),
    reject_005 = p_value < 0.05,
    reject_BH_005 = p_adj_BH < 0.05
  ) |>
  ungroup()

dagostino_x_summary <- dagostino_x_by_coordinate |>
  group_by(fileid, model) |>
  summarise(
    n_coordinates = n(),
    n_tested = sum(is.finite(p_value)),
    
    prop_reject_005 = mean(reject_005, na.rm = TRUE),
    prop_reject_BH_005 = mean(reject_BH_005, na.rm = TRUE),
    
    median_abs_skew = median(abs(skew), na.rm = TRUE),
    max_abs_skew = max(abs(skew), na.rm = TRUE),
    
    median_abs_excess_kurtosis =
      median(abs(excess_kurtosis), na.rm = TRUE),
    max_abs_excess_kurtosis =
      max(abs(excess_kurtosis), na.rm = TRUE),
    
    median_k2 = median(k2, na.rm = TRUE),
    max_k2 = max(k2, na.rm = TRUE),
    
    .groups = "drop"
  ) |>
  arrange(model, desc(prop_reject_BH_005))

dagostino_x_summary

dagostino_theta_by_coordinate <- all_theta_samples |>
  group_by(fileid, model, theta_index, i, j) |>
  summarise(
    dagostino_pearson_test(theta_raw),
    .groups = "drop"
  ) |>
  group_by(fileid, model) |>
  mutate(
    p_adj_BH = p.adjust(p_value, method = "BH"),
    reject_005 = p_value < 0.05,
    reject_BH_005 = p_adj_BH < 0.05
  ) |>
  ungroup()

dagostino_theta_summary <- dagostino_theta_by_coordinate |>
  group_by(fileid, model) |>
  summarise(
    n_theta_effects = n(),
    n_tested = sum(is.finite(p_value)),
    
    prop_reject_005 = mean(reject_005, na.rm = TRUE),
    prop_reject_BH_005 = mean(reject_BH_005, na.rm = TRUE),
    
    median_abs_skew_theta = median(abs(skew), na.rm = TRUE),
    max_abs_skew_theta = max(abs(skew), na.rm = TRUE),
    
    median_abs_excess_kurtosis_theta =
      median(abs(excess_kurtosis), na.rm = TRUE),
    max_abs_excess_kurtosis_theta =
      max(abs(excess_kurtosis), na.rm = TRUE),
    
    median_k2_theta = median(k2, na.rm = TRUE),
    max_k2_theta = max(k2, na.rm = TRUE),
    
    .groups = "drop"
  ) |>
  arrange(model, desc(prop_reject_BH_005))

dagostino_theta_summary

dagostino_combined_summary <- dagostino_x_summary |>
  left_join(
    dagostino_theta_summary,
    by = c("fileid", "model")
  ) |>
  arrange(model, fileid)

dagostino_combined_summary

dagostino_compact_table <- dagostino_combined_summary |>
  transmute(
    Dataset = fileid,
    Model = recode(
      model,
      "null" = "LBT",
      "full" = "LID",
      "only_theta" = "DABT",
      "full_theta" = "DALID",
      .default = model
    ),
    
    `x coords.` = n_coordinates,
    `Prop. reject x` = prop_reject_BH_005.x,
    `Med. |skew(x)|` = median_abs_skew,
    `Med. |ex. kurt.(x)|` = median_abs_excess_kurtosis,
    
    `theta coords.` = n_theta_effects,
    `Prop. reject theta` = prop_reject_BH_005.y,
    `Med. |skew(theta)|` = median_abs_skew_theta,
    `Med. |ex. kurt.(theta)|` = median_abs_excess_kurtosis_theta
  ) |>
  arrange(Model, Dataset)

dagostino_compact_table

dagostino_latex_table <- dagostino_compact_table |>
  mutate(
    across(
      where(is.numeric),
      ~ round(.x, 3)
    ),
    across(
      everything(),
      ~ ifelse(is.na(.x), "--", as.character(.x))
    )
  ) |>
  rename(
    Dataset = Dataset,
    Model = Model,
    `$n_x$` = `x coords.`,
    `$p_{\\mathrm{BH},x}$` = `Prop. reject x`,
    `Med. $|\\mathrm{skew}_x|$` = `Med. |skew(x)|`,
    `Med. $|\\mathrm{ex.kurt.}_x|$` = `Med. |ex. kurt.(x)|`,
    `$n_{\\theta}$` = `theta coords.`,
    `$p_{\\mathrm{BH},\\theta}$` = `Prop. reject theta`,
    `Med. $|\\mathrm{skew}_{\\theta}|$` = `Med. |skew(theta)|`,
    `Med. $|\\mathrm{ex.kurt.}_{\\theta}|$` = `Med. |ex. kurt.(theta)|`
  )

dagostino_latex_code <- dagostino_latex_table |>
  kable(
    format = "latex",
    booktabs = TRUE,
    escape = FALSE,
    caption = paste(
      "D'Agostino--Pearson omnibus normality diagnostics for the conditional",
      "MCMC samples. The table reports the proportion of marginal latent",
      "coordinates rejected after Benjamini--Hochberg correction, together",
      "with median absolute skewness and median absolute excess kurtosis."
    ),
    label = "tab:dagostino_compact"
  ) |>
  kable_styling(
    latex_options = c("hold_position", "scale_down")
  )

dagostino_latex_code

# ------------------------------------------------------------
# Effective sample size / IAT
# ------------------------------------------------------------

safe_iat <- function(x, max_lag = 1000) {
  x <- x[is.finite(x)]
  n <- length(x)
  
  if (n < 10) return(NA_real_)
  if (sd(x) == 0) return(NA_real_)
  
  max_lag <- min(max_lag, n - 1)
  
  acf_vals <- as.numeric(stats::acf(
    x,
    lag.max = max_lag,
    plot = FALSE,
    na.action = na.pass
  )$acf)[-1]
  
  # Simple positive-sequence truncation
  positive_acf <- acf_vals[acf_vals > 0]
  
  if (length(positive_acf) == 0) {
    return(1)
  }
  
  1 + 2 * sum(positive_acf)
}

safe_ess <- function(x, max_lag = 1000) {
  x <- x[is.finite(x)]
  iat <- safe_iat(x, max_lag = max_lag)
  
  if (!is.finite(iat) || iat <= 0) return(NA_real_)
  
  length(x) / iat
}

# ------------------------------------------------------------
# ESS-adjusted D'Agostino-Pearson omnibus diagnostic
# ------------------------------------------------------------

dagostino_pearson_test_ess <- function(x,
                                       ess = NULL,
                                       max_lag = 1000,
                                       min_ess = 20) {
  x <- x[is.finite(x)]
  n_raw <- length(x)
  
  if (n_raw < 20) {
    return(tibble::tibble(
      n_raw = n_raw,
      ess = NA_real_,
      skew = NA_real_,
      kurtosis = NA_real_,
      excess_kurtosis = NA_real_,
      z_skew = NA_real_,
      z_kurtosis = NA_real_,
      k2 = NA_real_,
      p_value = NA_real_
    ))
  }
  
  if (is.null(ess)) {
    ess <- safe_ess(x, max_lag = max_lag)
  }
  
  if (!is.finite(ess) || ess < min_ess) {
    return(tibble::tibble(
      n_raw = n_raw,
      ess = ess,
      skew = NA_real_,
      kurtosis = NA_real_,
      excess_kurtosis = NA_real_,
      z_skew = NA_real_,
      z_kurtosis = NA_real_,
      k2 = NA_real_,
      p_value = NA_real_
    ))
  }
  
  # D'Agostino formulas use sample size. Here we use ESS as a heuristic.
  n <- ess
  
  m <- mean(x)
  xc <- x - m
  
  m2 <- mean(xc^2)
  m3 <- mean(xc^3)
  m4 <- mean(xc^4)
  
  if (!is.finite(m2) || m2 <= 0) {
    return(tibble::tibble(
      n_raw = n_raw,
      ess = ess,
      skew = NA_real_,
      kurtosis = NA_real_,
      excess_kurtosis = NA_real_,
      z_skew = NA_real_,
      z_kurtosis = NA_real_,
      k2 = NA_real_,
      p_value = NA_real_
    ))
  }
  
  skew <- m3 / m2^(3 / 2)
  kurtosis <- m4 / m2^2
  excess_kurtosis <- kurtosis - 3
  
  # ----------------------------------------------------------
  # Skewness transformation
  # ----------------------------------------------------------
  
  y <- skew * sqrt(((n + 1) * (n + 3)) / (6 * (n - 2)))
  
  beta2 <- (3 * (n^2 + 27 * n - 70) * (n + 1) * (n + 3)) /
    ((n - 2) * (n + 5) * (n + 7) * (n + 9))
  
  w2 <- -1 + sqrt(2 * (beta2 - 1))
  delta <- 1 / sqrt(log(sqrt(w2)))
  alpha <- sqrt(2 / (w2 - 1))
  
  z_skew <- delta * log(y / alpha + sqrt((y / alpha)^2 + 1))
  
  # ----------------------------------------------------------
  # Kurtosis transformation
  # ----------------------------------------------------------
  
  expected_kurtosis <- 3 * (n - 1) / (n + 1)
  
  var_kurtosis <- (24 * n * (n - 2) * (n - 3)) /
    ((n + 1)^2 * (n + 3) * (n + 5))
  
  x_kurt <- (kurtosis - expected_kurtosis) / sqrt(var_kurtosis)
  
  sqrt_beta1 <- (
    6 * (n^2 - 5 * n + 2) /
      ((n + 7) * (n + 9))
  ) *
    sqrt(
      (6 * (n + 3) * (n + 5)) /
        (n * (n - 2) * (n - 3))
    )
  
  A <- 6 + (8 / sqrt_beta1) *
    (2 / sqrt_beta1 + sqrt(1 + 4 / sqrt_beta1^2))
  
  term1 <- 1 - 2 / (9 * A)
  term2 <- (1 - 2 / A) / (1 + x_kurt * sqrt(2 / (A - 4)))
  
  z_kurtosis <- (term1 - sign(term2) * abs(term2)^(1 / 3)) /
    sqrt(2 / (9 * A))
  
  k2 <- z_skew^2 + z_kurtosis^2
  p_value <- stats::pchisq(k2, df = 2, lower.tail = FALSE)
  
  tibble::tibble(
    n_raw = n_raw,
    ess = ess,
    skew = skew,
    kurtosis = kurtosis,
    excess_kurtosis = excess_kurtosis,
    z_skew = z_skew,
    z_kurtosis = z_kurtosis,
    k2 = k2,
    p_value = p_value
  )
}

dagostino_x_by_coordinate_ess <- all_x_samples |>
  group_by(fileid, model, agent, dimension) |>
  summarise(
    dagostino_pearson_test_ess(x),
    .groups = "drop"
  ) |>
  group_by(fileid, model) |>
  mutate(
    p_adj_BH = p.adjust(p_value, method = "BH"),
    reject_005 = p_value < 0.05,
    reject_BH_005 = p_adj_BH < 0.05
  ) |>
  ungroup()

safe_max <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  max(x)
}

safe_median <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  median(x)
}

safe_mean <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  mean(x)
}

dagostino_x_summary_ess <- dagostino_x_by_coordinate_ess |>
  group_by(fileid, model) |>
  summarise(
    n_coordinates = n(),
    n_tested = sum(is.finite(p_value)),
    
    median_ess = safe_median(ess),
    min_ess = {
      ess_finite <- ess[is.finite(ess)]
      if (length(ess_finite) == 0) NA_real_ else min(ess_finite)
    },
    median_n_raw = safe_median(n_raw),
    
    prop_reject_005 = ifelse(
      n_tested > 0,
      mean(reject_005, na.rm = TRUE),
      NA_real_
    ),
    prop_reject_BH_005 = ifelse(
      n_tested > 0,
      mean(reject_BH_005, na.rm = TRUE),
      NA_real_
    ),
    
    median_abs_skew = safe_median(abs(skew)),
    max_abs_skew = safe_max(abs(skew)),
    
    median_abs_excess_kurtosis =
      safe_median(abs(excess_kurtosis)),
    max_abs_excess_kurtosis =
      safe_max(abs(excess_kurtosis)),
    
    median_k2 = safe_median(k2),
    max_k2 = safe_max(k2),
    
    .groups = "drop"
  ) |>
  arrange(model, desc(prop_reject_BH_005))

dagostino_x_summary_ess

dagostino_theta_by_coordinate_ess <- all_theta_samples |>
  group_by(fileid, model, theta_index, i, j) |>
  summarise(
    dagostino_pearson_test_ess(theta_raw),
    .groups = "drop"
  ) |>
  group_by(fileid, model) |>
  mutate(
    p_adj_BH = p.adjust(p_value, method = "BH"),
    reject_005 = p_value < 0.05,
    reject_BH_005 = p_adj_BH < 0.05
  ) |>
  ungroup()

dagostino_theta_summary_ess <- dagostino_theta_by_coordinate_ess |>
  group_by(fileid, model) |>
  summarise(
    n_theta_effects = n(),
    n_tested = sum(is.finite(p_value)),
    
    median_ess_theta = safe_median(ess),
    min_ess_theta = {
      ess_finite <- ess[is.finite(ess)]
      if (length(ess_finite) == 0) NA_real_ else min(ess_finite)
    },
    
    prop_reject_005_theta = ifelse(
      n_tested > 0,
      mean(reject_005, na.rm = TRUE),
      NA_real_
    ),
    prop_reject_BH_005_theta = ifelse(
      n_tested > 0,
      mean(reject_BH_005, na.rm = TRUE),
      NA_real_
    ),
    
    median_abs_skew_theta = safe_median(abs(skew)),
    max_abs_skew_theta = safe_max(abs(skew)),
    
    median_abs_excess_kurtosis_theta =
      safe_median(abs(excess_kurtosis)),
    max_abs_excess_kurtosis_theta =
      safe_max(abs(excess_kurtosis)),
    
    median_k2_theta = safe_median(k2),
    max_k2_theta = safe_max(k2),
    
    .groups = "drop"
  ) |>
  arrange(model, desc(prop_reject_BH_005_theta))

dagostino_theta_summary_ess

dagostino_x_comparison <- dagostino_x_summary |>
  select(
    fileid,
    model,
    prop_reject_BH_005_raw_n = prop_reject_BH_005
  ) |>
  left_join(
    dagostino_x_summary_ess |>
      select(
        fileid,
        model,
        prop_reject_BH_005_ess = prop_reject_BH_005,
        median_ess,
        min_ess
      ),
    by = c("fileid", "model")
  ) |>
  mutate(
    reduction =
      prop_reject_BH_005_raw_n - prop_reject_BH_005_ess
  ) |>
  arrange(model, desc(reduction))

dagostino_x_comparison
