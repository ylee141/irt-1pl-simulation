# ==============================================================================
# Simulating data from a 1PL model
#
# Conventions
#   i = 1, ..., N : students (rows)
#   j = 1, ..., J : questions (columns)
#   theta[i]      : ability of student i
#   b[j]          : difficulty of question j
#   P[i, j]       : probability of student i getting question j correctly
# ==============================================================================

set.seed(1)

# ---- Settings ----------------------------------------------------------------
N <- 500  # number of test takers
J <- 10   # number of test questions
sigma_true <- 2 # standard deviation of the true ability distribution

# ---- True parameters ---------------------------------------------------------
theta <- rnorm(n = N, mean = 0, sd = sigma_true)  # abilities, length N
b     <- seq(-2, 2, length.out = J)      # difficulties, length J

# ---- Response probabilities --------------------------------------------------
# gap[i, j] = theta[i] - b[j]; rows = students, columns = questions
gap <- outer(theta, b, "-")  # N x J 
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

