# Examples to see how the calipeRmatching package works.

###############################################################################
# Load package
###############################################################################
library("calipeRmatching")
cm_set_number_of_threads(1) # set number of parallel threads

###############################################################################
# Generate data
###############################################################################
# --- x
# generate n-by-k covariate matrix.
generate_x <- function(n, k){
  x <- array(0, dim=c(n, k))
  for (i in 1:n){
    for (j in 1:k){
      x[i, j] <- runif(1, min=-10, max=10) # Uniform[-10,10] covariates
    }
  }
  return(x)
}
# --- propensity score
# computes propensity score based on logistic model
compute_propscore <- function(x, theta, modeltype){
  stopifnot(ncol(x)==2) # k=2 is assumed
  stopifnot(length(theta)== (2 + 1)) # k=2 is assumed
  stopifnot(modeltype == "logit") # logit model is assumed
  n <- nrow(x)
  k <- ncol(x)
  score <- vector(mode="numeric", n) # X*theta
  propscore <- vector(mode="numeric", n) # g(X*theta)
  for (i in 1:n){
    score_i = theta[1] # intercept
    for (j in 2:(k+1)){
      score_i <- score_i + theta[j] * x[i, j-1]
    }
    score[i] <- score_i
    propscore[i] <- plogis(score_i) # g(x_i^T*theta) propensity score
  }
  return(propscore)
}
# --- treatment
# generates treatment vector
generate_d <- function(propscore){
  n <- length(propscore)
  d <- vector(mode="numeric", n)
  for (i in 1:n){
    d[i] <- rbinom(1, 1, propscore[i]) # D_i | X_i ~ Bernoulli(g(X_i^T*theta))
  }
  return(d)
}
# --- outcome
# generates outcome y
generate_y <- function(x, d){
  n <- nrow(x)
  stopifnot(ncol(x)==2) # k=2 is assumed
  y <- vector(mode="numeric", n)
  mu0 <- 0.0
  mu1 <- 0.0 # imply ATE = mu1 - mu0 if E[X] = 0_k
  for (i in 1:n){
    y0i <- mu0 + 0.8 * x[i, 1] + 1.5 * x[i, 2] + 0.9 * rnorm(n=1)
    y1i <- mu1 + 3.1 * sin(x[i, 1]) - 1.8 * x[i, 2] + 1.1 * rnorm(n=1)
    y[i] <- y0i + d[i] * (y1i - y0i)
  }
  return(y) 
}

# setings
n <- 5000
k <- 2 # if changed, theta and the way that `y` is generated must  be changed too.
modeltype <- "logit" # if changed, the way that `propscore` is generated must be changed too.
theta <- c(0.004, 0.5, 0.2) # propensity score parameters; first entry is intercept

# generate data
x <- generate_x(n, k)
propscore <- compute_propscore(x, theta, modeltype)
d <- generate_d(propscore)
y <- generate_y(x, d)
# print(x)
# print(propscore)
# print(d)
# print(y)

###############################################################################
# Estimation (known propensity score)
###############################################################################

# test library
test_cm()

# estimation
delta <- 0.0 # caliper; set to zero for default value
estimate_variance <- TRUE # do estimate variance
alpha <- 0.0
beta <- 0.0
kappa_a <- 0.0
kappa_gamma <- 0.0
# for categorical treatment indicator use # d <- as.numeric(d=="string_for_treated") 
help("cm_cm_known_propscore")
result_known_propscore <- cm_cm_known_propscore(y, d, propscore, delta, estimate_variance, beta, alpha, kappa_a, kappa_gamma)
print(result_known_propscore)
ate_hat <- result_known_propscore$point_estimates$ate_hat
var_hat_ate <- result_known_propscore$point_estimates$var_hat_ate
ate_ci95_low <- ate_hat - qnorm(0.975) * sqrt(var_hat_ate/length(y))
ate_ci95_high <- ate_hat + qnorm(0.975) * sqrt(var_hat_ate/length(y))
print(sprintf("95%% asymptotically valid confidence interval for ATE: [%f, %f]", ate_ci95_low, ate_ci95_high))

###############################################################################
# Estimation (estimated propensity score)
###############################################################################

# split sample
# (For real-wrold data use random partition of indices; here splitting at half is suitable)
ps_estimation_index <- 1:(n/2) # use first half of sample for propensity score estimation
cm_estimation_index <- ((n/2)+1):n # use second half of sample for matching estimation

# propensity score estimation
x_ps_df <- data.frame(x[ps_estimation_index,])
x_ps_model <- model.matrix(~ X1 + X2, x_ps_df) # covariate matrix with intercept
d_ps_model <- d[ps_estimation_index]
propscore_model_estimated <- glm.fit(x=x_ps_model, y=d_ps_model, family=binomial(link=modeltype))
theta_hat <- propscore_model_estimated$coefficients

# caliper matching
x_cm_df <- data.frame(x[cm_estimation_index,])
x_cm_model <- model.matrix(~ X1 + X2, x_cm_df)
x_cm_model <- x_cm_model[,-(1)] # caliper matching doesn't take intercept as input, so drop it
d_cm_model <- d[cm_estimation_index]
y_cm_model <- y[cm_estimation_index]

kappa_gamma_derivative <- 0.0
help("cm_cm_estimated_propscore")
result_estimated_propscore <- cm_cm_estimated_propscore(y_cm_model, d_cm_model, x_cm_model, modeltype, theta_hat, delta, estimate_variance, beta, alpha, kappa_a, kappa_gamma, kappa_gamma_derivative)
print(result_estimated_propscore)
ate_hat <- result_estimated_propscore$point_estimates$ate_hat
var_hat_ate <- result_estimated_propscore$point_estimates$var_hat_ate
ate_ci95_low <- ate_hat - qnorm(0.975) * sqrt(var_hat_ate/length(y_cm_model)) # NOTE the use of y_cm_model not y!
ate_ci95_high <- ate_hat + qnorm(0.975) * sqrt(var_hat_ate/length(y_cm_model)) # NOTE the use of y_cm_model not y!
print(sprintf("95%% asymptotically valid confidence interval for ATE: [%f, %f]", ate_ci95_low, ate_ci95_high))






