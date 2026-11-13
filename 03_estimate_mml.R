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

source("01_simulate_data.R") # creates N, J, theta, b, X


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

fit <- optim(
  par     = start,
  fn      = marginal_log_likelihood,
  X       = X,
  draws   = draws,
  method  = "BFGS",
  hessian = TRUE,
  control = list(fnscale = -1, maxit = 1000)
)

fit$convergence # 0 = converged

# ---- Results -----------------------------------------------------------------

b_mml <- fit$par[1:J]
round(cbind(true = b, mml = b_mml - mean(b_mml)), 2) # centered difficulties
exp(fit$par["log_sigma"]) # estimated sigma
