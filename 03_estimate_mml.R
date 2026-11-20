# ==============================================================================
# MML Estimation
#
# Conventions
#   i = 1, ..., N : students (rows)
#   j = 1, ..., J : questions (columns)
#   r = 1, ..., R : simulation draws
#   X             : N x J matrix of simulated responses (the gradebook)
#
# True values (from the simulation, used only for comparison)
#   b_true        : difficulties
#   sigma_true    : SD of abilities, theta ~ N(0, sigma_true^2)
#
# Estimated parameters
#   b[j]          : difficulty of question j
#   sigma         : SD of abilities in the population
#
# Simulation of the integral over ability
#   draws[i, r]   : fixed standard normal draw r for student i
#   theta[i, r]   : simulated ability = sigma * draws[i, r]
#   P[i, r, j]    : probability of a correct answer at that ability
#                   = logistic(theta[i, r] - b[j])
# ==============================================================================

source("01_simulate_data.R")  # creates N, J, sigma_true, theta_true, b_true, X

# ---- Simulation draws --------------------------------------------------------

# Fixed standard normal draws, R per student, created once outside the
# likelihood so it stays smooth across optimizer iterations.

set.seed(2) 
R <- 100
draws <- matrix(rnorm(N * R, 0, 1), nrow = N, ncol = R) # N * R

# ---- Marginal log-likelihood -------------------------------------------------

# logsumexp helper

marginal_log_likelihood <- function(params, X, draws) {
  N <- nrow(X)
  J <- ncol(X)
  R <- ncol(draws)
  
  b <- params[1:J]
  log_sigma <- params[J + 1] # optimizer searches over any real number
  sigma <- exp(log_sigma) # so sigma stays positive
  
  theta <- sigma * draws # N x R, simulated abilities
  gap <- outer(theta, b, "-") # N x R x J
  
  log_p <- plogis(gap, log.p = TRUE) # log(P)
  log_1mp <- plogis(gap, log.p = TRUE, lower.tail = FALSE) # log(1 - P)
  
  # Log-probability of each student's answer pattern at each draw
  ll_ir <- sapply(1:R, function(r) {
    rowSums(X * log_p[, r, ] + (1 - X) * log_1mp[, r, ])
  }) # N * R
  
  # Average over draws for each student (log-sum-exp needed for long tests)
  ll_i <- log(rowMeans(exp(ll_ir))) # length N
  
  sum(ll_i)
}

# ---- Optimization ------------------------------------------------------------

start <- rep(0, J + 1) # all b = 0; log_sigma = 0, i.e. sigma = 1
names(start) <- c(paste0("b_", 1:J), "log_sigma")

fit_mml <- optim(
  par     = start,
  fn      = marginal_log_likelihood,
  X       = X,
  draws   = draws,
  method  = "BFGS",
  hessian = TRUE,
  control = list(fnscale = -1, maxit = 1000)
)

fit_mml$convergence # 0 = converged

# ---- Estimates ---------------------------------------------------------------

b_mml     <- fit_mml$par[1:J]
sigma_mml <- unname(exp(fit_mml$par["log_sigma"]))

# ---- Standard errors ---------------------------------------------------------

V_mml    <- solve(-fit_mml$hessian)  # variance matrix of (b, log_sigma)
b_se_mml <- sqrt(diag(V_mml)[1:J])

# sigma was estimated as log(sigma), so its SE comes from the delta method:
# SE(sigma) = sigma * SE(log sigma)
sigma_se_mml <- sigma_mml * sqrt(V_mml[J + 1, J + 1])

# ---- Results -----------------------------------------------------------------

# MML sets zero at the population's mean ability, while the true b's sum to 0,
# so the difficulties are centered before comparing them with the truth.
round(cbind(true = b_true, mml = b_mml - mean(b_mml)), 2)

c(sigma_true = sigma_true, sigma_mml = sigma_mml)

# ---- Check against mirt ------------------------------------------------------

# mirt fits the same model (Rasch, mean theta = 0, sigma estimated), but
# integrates by quadrature instead of simulation draws. Both use the same
# zero, so the estimates and SEs are compared without centering.

library(mirt)

fit_mirt <- mirt(as.data.frame(X), model = 1, itemtype = "Rasch",
                 SE = TRUE, verbose = FALSE)
cf <- coef(fit_mirt, printSE = TRUE)  # one matrix per item, plus GroupPars

# mirt writes each item as logistic(theta + d), so b = -d and SE(b) = SE(d)
b_mirt    <- -sapply(cf[1:J], function(m) m["par", "d"])
b_se_mirt <-  sapply(cf[1:J], function(m) m["SE",  "d"])

# mirt reports the variance of theta; SE(sigma) = SE(variance) / (2 * sigma)
var_mirt      <- cf$GroupPars["par", "COV_11"]
sigma_mirt    <- sqrt(var_mirt)
sigma_se_mirt <- cf$GroupPars["SE", "COV_11"] / (2 * sigma_mirt)

round(cbind(mml     = b_mml,  se_mml  = b_se_mml,
            mirt    = b_mirt, se_mirt = b_se_mirt), 3)

round(c(sigma_mml  = sigma_mml,  se_mml  = sigma_se_mml,
        sigma_mirt = sigma_mirt, se_mirt = sigma_se_mirt), 3)