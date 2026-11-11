# ==============================================================================
# JML Estimation with bias correction
#
# Conventions
#   i = 1, ..., N : students (rows)
#   j = 1, ..., J : questions (columns)
#   theta[i]      : ability of student i
#   b[j]          : difficulty of question j (true values, from simulation)
#   P[i, j]       : probability that student i answers question j correctly
#                   = logistic(theta[i] - b[j])
#   X             : N x J matrix of simulated responses (the gradebook)
#   X_fit         : X without students who scored 0 or J
# ==============================================================================

# -------------------------------------------------------------------------



source("01_simulate_data.R")  # creates N, J, theta, b, X

# ---- log-likelihood -------------------------------------

joint_log_likelihood <- function(params, X) {
  N <- nrow(X)
  J <- ncol(X)
  
  # N + J - 1 free parameters; the last difficulty is fixed by sum(b) = 0
  theta  <- params[1:N]
  b_free <- params[ (N+1) : (N+J-1)]
  b      <- c(b_free,-sum(b_free))  #last question difficulty is normalized
  
  gap <- outer(theta,b,"-")  # gap per student - difficulty pair
  P   <- plogis(gap)         # theoretical probability of being correct
  
  sum(X * log(P)  + (1-X)  * log((1-P)))
}

# ---- Data for estimation -----------------------------------------------------

# Perfect and zero scores have no finite JML estimate, so drop them
scores <- rowSums(X)
keep <- scores > 0 & scores < J
X_fit <- X[keep,]
sum(!keep)  # number of students dropped

N_fit <- nrow(X_fit) 
idx_b <- (N_fit + 1):(N_fit + J - 1)   # positions of the free difficulties

# ---- Optimization ------------------------------------------------------------

# BFGS follows the slope, which suits this smooth, single-peaked log-likelihood.
# Nelder-Mead uses function values only and was far too slow with ~500 params.
start <- rep(0, N_fit + J - 1)  # start all parameters at 0

fit <- optim(
  par = start,
  fn =joint_log_likelihood, 
  X=X_fit, 
  hessian = TRUE,
  method = "BFGS", 
  control = list(fnscale = -1, maxit = 1000)
) 

fit$convergence  # 0 = converged

# ---- Estimates and bias correction -------------------------------------------

b_free     <- fit$par[idx_b]
b_hat      <- c(b_free, -sum(b_free))
b_hat_corr <- b_hat * (J - 1) / J  # approximate JML bias correction

# ---- Standard errors ---------------------------------------------------------

V       <- solve(-fit$hessian)  # variance matrix of all estimates
V_b     <- V[idx_b, idx_b]
se_free <- sqrt(diag(V_b))
se_last <- sqrt(sum(V_b))       # Var(-sum(b_free)) = sum of all entries of V_b

b_se      <- c(se_free, se_last)
b_se_corr <- b_se * (J - 1) / J

# ---- Results -----------------------------------------------------------------

round(cbind(true      = b,
            estimate  = b_hat,
            corrected = b_hat_corr,
            se_corr   = b_se_corr), 2)
