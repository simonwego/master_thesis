library(RTMB)
set.seed(3)
n <- 50
x <- rnorm(n)
y <- rnorm(n)
y[n] <- 0 # Just for visualization purposes in final plot
a <- 0
b <- 0
c <- 1
tmp <- combn(1:n,2)
i <- tmp[1,]
tmp
j <- tmp[2,]
eta <- a*(x[i] - x[j]) + b*(y[i] - y[j]) + c*(x[i]*y[j] - x[j]*y[i])
eta
size <- 2
z <- rbinom(length(i), size=size, prob=pnorm(eta))
z
plot(x,y)
data <- data.frame(i,j,size,z)
data
dim(data)

f <- function(parms, data) {
  getAll(data, parms, warn=FALSE)
  z <- OBS(z)
  nll <- 0.0
  nll <- nll - sum(dnorm(x, log=TRUE))
  nll <- nll - sum(dnorm(y, log=TRUE))
  r2 <- x[1]^2 + y[1]^2
  r <- sqrt(r2) + 1e-12
  ct <- x[1] / r
  st <- y[1] / r
  xr <- ct*x - st*y
  yr <- st*x + ct*y
  yr[1] <- 0
  eta <- a*(xr[i] - xr[j]) + b*(yr[i] - yr[j]) + c*(xr[i]*yr[j] - xr[j]*yr[i])
  p <- pnorm(eta)
  ADREPORT(p)
  nll <- nll - sum(dbinom(z, size, prob=p, log=TRUE))
  nll
}

parameters <- list(
  u_r=0.0,
  u_theta=0.0,
  a=1.0,
  b=1.0,
  c=1,
  xf=rep(0,n-1),
  yf=rep(0,n-1)
)

f_constrained <- function(parms, data) {
  getAll(data, parms, warn=FALSE)
  z <- OBS(z)
  
  #Type <- typeof(a)  # RTMB pleier ikke trenge dette; bare for lesbarhet
  
  nll <- 0
  
  # 1) u_r ~ N(0,1)
  nll <- nll - dnorm(u, 0, 1, log=TRUE)
  
  # 2) r^2 ~ chisq_2  <=> Gamma(shape=1, scale=2)
  v  <- pnorm(u, 0, 1)
  r2 <- qgamma(v, shape=1, scale=2)
  r  <- sqrt(r2) + 1e-12
  
  # 3) Priors for free coords (x[2:n], y[2:n])
  nll <- nll - sum(dnorm(xf, 0, 1, log=TRUE))
  nll <- nll - sum(dnorm(yf, 0, 1, log=TRUE))
  
  # 4) Build full latent vectors with gauge-fix
  x <- c(r, xf)       # length n
  y <- c(0.0, yf)     # length n
  
  #a <- exp(log_a)
  #b <- exp(log_b)
  # 5) Likelihood
  eta <- a*(x[i] - x[j]) + b*(y[i] - y[j]) + c*(x[i]*y[j] - x[j]*y[i])
  p <- pnorm(eta)
  
  ADREPORT(p)
  ADREPORT(r)
  
  nll <- nll - sum(dbinom(z, size, prob=p, log=TRUE))
  nll
}

f_gauge <- function(parms, data) {
  getAll(data, parms, warn=FALSE)
  z <- OBS(z)
  nll <- 0
  
  ## --- 1) r^2 ~ chisq_2 via transform-trikset
  nll <- nll - dnorm(u_r, 0, 1, log=TRUE)
  v_r  <- pnorm(u_r, 0, 1)                 # U(0,1)
  r2   <- qgamma(v_r, shape=1, scale=2)    # chisq_2
  r    <- sqrt(r2) + 1e-12                 # chi_2 / Rayleigh
  
  ## --- 2) theta ~ Uniform(-pi, pi)
  nll <- nll - dnorm(u_theta, 0, 1, log=TRUE)
  v_th   <- pnorm(u_theta, 0, 1)           # U(0,1)
  theta  <- (2*pi) * v_th - pi
  
  ## --- 3) priors for the remaining points (i=2..n)
  nll <- nll - sum(dnorm(xf, 0, 1, log=TRUE))
  nll <- nll - sum(dnorm(yf, 0, 1, log=TRUE))
  
  ## --- 4) build "world-frame" points (before gauge-rotation)
  x_raw <- c(r*cos(theta), xf)   # length n
  y_raw <- c(r*sin(theta), yf)
  
  ## --- 5) rotate ALL points by -theta so that y1' = 0
  ct <- cos(theta); st <- sin(theta)
  x <-  ct*x_raw + st*y_raw          # rotation by -theta
  y <- -st*x_raw + ct*y_raw
  
  # Now: x[1] == r, y[1] == 0 (up to numerical eps)
  
  ## --- 6) your model
  eta <- a*(x[i] - x[j]) + b*(y[i] - y[j]) + c*(x[i]*y[j] - x[j]*y[i])
  p <- pnorm(eta)
  
  ADREPORT(r)
  ADREPORT(theta)
  ADREPORT(p)
  
  nll <- nll - sum(dbinom(z, size, prob=p, log=TRUE))
  nll
}


cmb <- function(f, d) function(p) f(p, d)
obj <- MakeADFun(cmb(f_gauge, data), parameters, random=c("xf","yf"))
obj$fn(obj$par)
opt <- nlminb(obj$par, obj$fn, obj$gr)
sdr <- sdreport(obj)
exp(sdr$par.fixed['log_a'])
sdr
plot(x,y)
xhat <- c(sdr$par.random[1:n])
yhat <- c(sdr$par.random[(n+1):(2*n-1)],0)
arrows(x,y,xhat,yhat)
abline(h=0, lty=2)

res <- replicate(1e+4, {
  # Simulate two bivariate standard normal vectors x1 and x2
  x1 <- rnorm(2)
  x2 <- rnorm(2)
  # Rotate both vectors such that x1 aligns with the first axis
  theta <- -atan2(x1[2],x1[1])
  R <- rbind(c(cos(theta), -sin(theta)), c(sin(theta),cos(theta)))
  x1p <- R %*% x1
  x2r <- R %*% x2
  # Return the resulting non-zero components
  c(x11=x1p[1], x2=x2r)
})
# Plot the results
dchi <- function(x, df) {
  dchisq(x^2, df=df)*2*x
}
par(mfcol=c(3,1))
MASS::truehist(res["x11",]); curve(dchi(x, df=2), add=TRUE)
MASS::truehist(res["x21",]); curve(dnorm, add=TRUE)
MASS::truehist(res["x22",]); curve(dnorm, add=TRUE)
par(mfrow=c(1,1))
pairs(t(res), pch=".")
t(res)
