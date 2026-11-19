# ==============================================================================
# Simulating data from a 1PL model
#
# Conventions
#   i = 1, ..., N : students (rows)
#   j = 1, ..., J : questions (columns)
#   sigma_true    : SD of the true abilities, theta_true ~ N(0, sigma_true^2)
#   theta_true[i] : true ability of student i
#   b_true[j]     : true difficulty of question j
#   P[i, j]       : probability that student i answers question j correctly
#                   = logistic(theta_true[i] - b_true[j])
#   X             : N x J matrix of simulated responses (the gradebook)
# ==============================================================================

set.seed(1)

# ---- Settings ----------------------------------------------------------------
N <- 500  # number of test takers
J <- 10   # number of test questions

# ---- True parameters ---------------------------------------------------------
sigma_true <- 2 # standard deviation of the true ability distribution
theta_true <- rnorm(n = N, mean = 0, sd = sigma_true)  # abilities, length N
b_true     <- seq(-2, 2, length.out = J)      # difficulties, length J

# ---- Response probabilities --------------------------------------------------
# gap[i, j] = theta[i] - b[j]; rows = students, columns = questions
gap <- outer(theta_true, b_true, "-")  # N x J 
P   <- plogis(gap)           # Pr(correct) 

# ---- Simulated responses  ----------------------------------------------------
# rbinom reads results as a long vector, column by column
# matrix() refills the results in the same order.
x_vec <- rbinom(N * J, 1, P) # N*J results, question by question
X     <- matrix(x_vec, nrow = N, ncol = J)
 
stopifnot(dim(X) == c(N, J))

# ---- Sanity check ------------------------------------------
# Share correct per question; should fall as b rises
colMeans(X)

