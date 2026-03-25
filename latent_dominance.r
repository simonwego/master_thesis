repos <- "https://cloud.r-project.org"
pkgs  <- c("TMB", "glmmTMB", "lme4")
to_install <- setdiff(pkgs, rownames(installed.packages()))
if (length(to_install)) {
  install.packages(to_install, repos = repos, dependencies = TRUE)
}
library(TMB)

inv_logit <- function(eta) 1 / (1 + exp(-eta))

simulate_dominance_data <- function(
  s = 100,
  params = list(a = 1.0, b = 0.5, c = 0.8),
  lam = 10,
  seed = NULL,
  return_matrix = TRUE
) {
  if (!is.null(seed)) set.seed(seed)

  a <- params$a; b <- params$b; c <- params$c

  # Latente og observerte individkovariater
  z <- rnorm(s)
  y <- rnorm(s)

  # Bygg p-matrise vektorisert (raskt og kompakt)
  dz <- outer(z, z, "-")                       # z_i - z_j
  dy <- outer(y, y, "-")                       # y_i - y_j
  cross <- outer(z, y, "*") - outer(y, z, "*") # z_i y_j - z_j y_i
  eta <- a * dz + b * dy + c * cross
  p <- inv_logit(eta)

  # Sample kun C8vre trekant (i<j)
  iu <- which(upper.tri(p), arr.ind = TRUE)  # matrise med kolonner: row, col
  K <- nrow(iu)

  n_ij <- rpois(K, lambda = lam)
  x_ij <- rbinom(K, size = n_ij, prob = p[iu])

  # Edge list (0- eller 1-indeksering avhenger av C++-koden din)
  # TMB/R er 1-indeksert, men i C++ bruker du ofte 0-indekser:
  i <- iu[, 1]
  j <- iu[, 2]

  # Valgfritt: full X og N matrise for debugging/plotting
  X <- N <- NULL
  if (return_matrix) {
    X <- matrix(0L, s, s)
    N <- matrix(0L, s, s)

    X[cbind(i, j)] <- x_ij
    N[cbind(i, j)] <- n_ij

    X[cbind(j, i)] <- n_ij - x_ij
    N[cbind(j, i)] <- n_ij

    diag(X) <- 0L
    diag(N) <- 0L
  }

  # Det du typisk sender til TMB: edge list + y (og evt. df etc.)
  list(
    y = y,
    # z er latent i modellen din; i simulering kan du returnere den for evaluering
    z_true = z,
    params_true = params,

    # Edge list data (for likelihood-sum i C++)
    i = i,
    j = j,
    x = x_ij,
    n = n_ij,

    # Optional debug objects
    X = X,
    N = N,
    p = p
  )
}

# simulere data
data <- simulate_dominance_data(s = 100, lam = 0.5, seed = 1)
# endre til nullindeksering for C++
data$i <- data$i - 1L
data$j <- data$j - 1L
data$use_c <- 1
str(data[c("y","i","j","x","n","use_c")])
data$x
data$n

parameters <- list(
  log_a = 0,
  b = 0,
  c = 0,
  z = rep(0, length(dat$y))
)

#Getting TMB template
tmb_dir  <- tmb_dir <- "C:/Users/simon/OneDrive/Skrivebord/Masteroppgave/master_thesis"
cpp      <- file.path(tmb_dir, "latent_dominance.cpp")
dll_base <- file.path(tmb_dir, "latent_dominance")
# Unload hvis allerede lastet
try(dyn.unload(TMB::dynlib(dll_base)), silent=TRUE)
# Alltid recompile
TMB::compile(cpp)
dyn.load(TMB::dynlib(dll_base))


obj <- TMB::MakeADFun(data, parameters, random = "z", DLL = "latent_dominance")
opt <- nlminb(obj$par, obj$fn, obj$gr)
rep <- sdreport(obj)
summary(rep)[c("a","b","c"), ]
rep$cov
  