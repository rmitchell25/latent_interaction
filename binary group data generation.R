# remotes::update_packages('rblimp')

library(Matrix)
library(mvtnorm)
library(rblimp)
options(scipen = 999)

################################################################################
# Simulation Conditions
################################################################################

corr_Xs <- .20
rsq_baseline <- .13
rsq_prod <- c(0,0.03,0.07)
group_probs3 <- rbind(c(.34, .33, .33),
                      c(.40, .40, .20),
                      c(.60, .20, .20))
group_probs2 <- rbind(c(.50,.50),
                      c(.60,.40),
                      c(.80,.20))
sample_size <- c(seq(100,400, by = 50),500,1000)
loading_size <- c(.5,.8)
num_loadings <- c(6,12)

# selected test conditions
N <- sample_size[4]
loading <- loading_size[2]
probs <- group_probs2[2,]
rsq_product <- rsq_prod[2]
num_indicators <- num_loadings[1]

################################################################################
# Step 1: Define Group-Specific Means of X to Induce Correlation
################################################################################

solve_factor_means_2group <- function(probs, target_corr) {
  p1 <- probs[1] 
  p2 <- probs[2]
  
  stopifnot(abs(p1 + p2 - 1) < 1e-6)
  
  # Define objective function to minimize the deviation from target correlation
  objective <- function(mu2) {
    mu1 <- 0  # fix group 1 mean to zero
    mu_bar <- p1 * mu1 + p2 * mu2  # expected value of X (should be 0)
    
    # Between-group variance
    var_between <- p1 * (mu1 - mu_bar)^2 + p2 * (mu2 - mu_bar)^2
    var_within <- 1 - var_between  # to keep total variance of X at 1
    var_X <- var_between + var_within  # should be ~1
    
    # Covariance between X and G2 dummy code
    cov_X_G2 <- p2 * (mu2 - mu_bar)
    
    # Variance of G2 dummy variable
    var_G2 <- p2 * (1 - p2)
    
    # Compute actual correlation
    corr_G2 <- cov_X_G2 / sqrt(var_X * var_G2)
    
    # Return squared error from target
    return((corr_G2 - target_corr)^2)
  }
  
  # Optimize over a reasonable range
  opt <- optimize(objective, interval = c(-3, 3))
  
  mu1 <- 0
  mu2 <- opt$minimum
  
  return(c(mu1 = mu1, mu2 = mu2))
}

factor_means <- solve_factor_means_2group(probs = probs, target_corr = corr_Xs)

################################################################################
# Step 2: Solve for Within-Group Variance of X to Achieve Var(X) = 1
################################################################################

solve_factor_variance_2group <- function(mu_X, probs, target_var = 1) {
  stopifnot(length(mu_X) == 2, length(probs) == 2, abs(sum(probs) - 1) < 1e-6)
  
  # Pooled (overall) mean of X
  mu_bar <- sum(probs * mu_X)
  
  # Between-group variance component
  var_between <- sum(probs * (mu_X - mu_bar)^2)
  
  # Solve for within-group variance to achieve total variance = target_var
  var_within <- target_var - var_between
  
  # All group variances assumed equal
  var_group <- rep(var_within, 2)
  names(var_group) <- paste0("sigma_sq_g", 1:2)
  
  return(list(
    mu_bar = mu_bar,
    var_between = var_between,
    var_within = var_within,
    group_variances = var_group
  ))
}

factor_vars <- solve_factor_variance_2group(mu_X = factor_means, 
                                     probs = probs)$group_variances

################################################################################
# Step 3: Solve for Group-Specific Slopes and Intercepts (Structural Parameters)
################################################################################

solve_group_parameters_2group <- function(mu_X, var_X, probs, rsq_baseline, rsq_product) {
  
  stopifnot(length(mu_X) == 2, length(var_X) == 2, length(probs) == 2, abs(sum(probs) - 1) < 1e-6)
  
  mu_X <- as.numeric(mu_X)
  var_X <- as.numeric(var_X)
  probs <- as.numeric(probs)
  
  # 1. Fix baseline slope (assumes full sample Var(X) = 1)
  beta_baseline <- sqrt(rsq_baseline)
  
  # 2. Optimize delta2 (only one delta needed in two-group case)
  objective <- function(delta2) {
    deltas <- c(0, delta2)  # Group 1 is reference
    betas <- beta_baseline + deltas
    
    # Compute group-level explained variance
    y_means <- betas * mu_X
    mean_y <- sum(probs * y_means)
    
    explained_var <- sum(probs * (betas^2 * var_X + (y_means - mean_y)^2))
    achieved_rsq <- explained_var
    
    return((achieved_rsq - (rsq_baseline + rsq_product))^2)
  }
  
  delta2_opt <- optimize(objective, interval = c(0, 5))$minimum
  deltas <- c(0, delta2_opt)
  betas <- beta_baseline + deltas
  
  # 3. Shared intercept to ensure E[Y] = 0
  alpha_shared <- -sum(probs * betas * mu_X)
  alphas <- rep(alpha_shared, 2)
  
  # 4. Residual variance to ensure Var(Y) = 1
  y_means <- betas * mu_X + alphas
  var_Y <- sum(probs * (betas^2 * var_X + (y_means - sum(probs * y_means))^2))
  residual_var <- 1 - var_Y
  sigmas_sq <- rep(residual_var, 2)
  
  names(betas) <- names(alphas) <- names(sigmas_sq) <- paste0("G", 1:2)
  
  return(list(
    intercepts = alphas,
    slopes = betas,
    residual_variances = sigmas_sq
  ))
}


group_params <- solve_group_parameters_2group(mu_X = factor_means, 
                                       var_X = factor_vars, 
                                       probs = probs,
                                       rsq_baseline = rsq_baseline,
                                       rsq_product = rsq_product)
print(group_params)

################################################################################
# Step 4: Convert to Pooled Moderated Regression Parameters
################################################################################

group_to_moderated_2group <- function(mu_X, var_X, probs, group_params) {
  
  stopifnot(length(mu_X) == 2, length(var_X) == 2, length(probs) == 2)
  
  # Extract components from group_params
  intercepts <- group_params$intercepts
  slopes <- group_params$slopes
  resid_vars <- group_params$residual_variances
  
  # 1. Reference group: Group 1
  beta_0 <- intercepts[1]
  beta_X <- slopes[1]
  
  # 2. Dummy-coded differences for Group 2
  beta_G2 <- intercepts[2] - intercepts[1]
  beta_XG2 <- slopes[2] - slopes[1]
  
  # 3. Pooled residual variance
  residual_variance <- sum(probs * resid_vars)
  
  return(list(
    beta_0 = beta_0,
    beta_G2 = beta_G2,
    beta_X = beta_X,
    beta_XG2 = beta_XG2,
    residual_variance = residual_variance
  ))
}

moderated_params <- group_to_moderated_2group(
  mu_X = factor_means,
  var_X = factor_vars,
  probs = probs,
  group_params = group_params
)

print(moderated_params)

################################################################################
# Check Structural Parameters With a Large Data Set
################################################################################

set.seed(90291)
N_large <- 10000000
ng <- N_large * probs

X1 <- rnorm(ng[1],factor_means[1],sqrt(factor_vars[1]))
X2 <- rnorm(ng[2],factor_means[2],sqrt(factor_vars[2]))

X1s <- cbind(0,X1,0*X1)
X2s <- cbind(1,X2,1*X2)
Xs <- rbind(X1s,X2s)

Y1 <- rnorm(ng[1], group_params$intercepts[1] + X1*group_params$slopes[1], sqrt(group_params$residual_variances[1]))
Y2 <- rnorm(ng[2], group_params$intercepts[2] + X2*group_params$slopes[2], sqrt(group_params$residual_variances[2]))
Y <- c(Y1,Y2)

dat <- as.data.frame(cbind(Y,Xs))
names(dat) <- c('Y','G2','X','G2X')

summary(lm(Y ~ G2 + X + G2X, data = dat)) # ALL GOOD


################################################################################
# Step 5: Measurement Model for Latent Factors
################################################################################

model_implied_moments <- function(mu_X, var_X,
                                  alpha, beta, residual_var_Y,
                                  n_X = num_indicators, n_Y = num_indicators,
                                  stanload = loading) {
  
  # Measurement parameters
  lambda_X <- stanload
  theta_X  <- 1 - lambda_X^2  
  
  lambda_Y <- 1.00
  theta_Y  <- 1 / stanload^2 - 1 # assumes that all raw Y loadings = 1
  
  # Latent means
  mu_Y <- alpha + beta * mu_X
  mu_indicators <- c(rep(lambda_X * mu_X, n_X),
                     rep(lambda_Y * mu_Y, n_Y))
  names(mu_indicators) <- c(paste0("X",1:n_X), paste0("Y",1:n_Y))
  
  # Latent covariance matrix (2×2)
  # Var(X) = var_X; Cov(X,Y) = beta * var_X; Var(Y) = beta^2 var_X + residual_var_Y
  Sigma_latent <- matrix(c(var_X, beta * var_X, beta * var_X, beta^2 * var_X + residual_var_Y),nrow = 2, byrow = TRUE)
  
  # Loading matrix
  Lambda_X <- matrix(lambda_X, n_X, 1)
  Lambda_Y <- matrix(lambda_Y, n_Y, 1)
  Lambda   <- rbind(cbind(Lambda_X, matrix(0, n_X, 1)),
                    cbind(matrix(0, n_Y, 1),  Lambda_Y))
  
  # Measurement error covariance
  Theta_X <- diag(rep(theta_X, n_X))
  Theta_Y <- diag(rep(theta_Y, n_Y))
  Theta   <- as.matrix(bdiag(Theta_X, Theta_Y))
  
  # Model-implied observed covariance
  Sigma_obs <- Lambda %*% Sigma_latent %*% t(Lambda) + Theta
  rownames(Sigma_obs) <- colnames(Sigma_obs) <- names(mu_indicators)
  
  return(list(mean = mu_indicators, covariance = Sigma_obs))
}

# Group 1 mean vector and covariance matrix
moments_G1 <- model_implied_moments(
  mu_X           = factor_means["mu1"],
  var_X          = factor_vars["sigma_sq_g1"],
  alpha          = group_params$intercepts["G1"],
  beta           = group_params$slopes      ["G1"],
  residual_var_Y = group_params$residual_variances["G1"]
)

# Group 2 mean vector and covariance matrix
moments_G2 <- model_implied_moments(
  mu_X           = factor_means["mu2"],
  var_X          = factor_vars["sigma_sq_g2"],
  alpha          = group_params$intercepts["G2"],
  beta           = group_params$slopes      ["G2"],
  residual_var_Y = group_params$residual_variances["G2"]
)

################################################################################
# Check Latent Interaction Model With a Large Data Set
################################################################################

set.seed(90291)
N_large <- 100000
ng <- N_large * probs

g1dat <- cbind(1,rmvnorm(ng[1], as.vector(moments_G1$mean), as.matrix(moments_G1$covariance)))
g2dat <- cbind(2,rmvnorm(ng[2], as.vector(moments_G2$mean), as.matrix(moments_G2$covariance)))
dat <- as.data.frame(rbind(g1dat,g2dat))
names(dat) <- c('G',paste0('X',1:num_indicators),paste0('Y',1:num_indicators))

syntax <- list(
  structural_model = c(
    'Y ~ 1 G X G*X'
  ),
  X_measurement = c(
    'X ~~ X@1',
    paste0('X -> X1@lx1 X2:X', num_indicators, collapse = '')
  ),
  Y_measurement = c(
    paste0('Y -> Y1:Y', num_indicators, collapse = ''),
    'Y1 ~ 1@0'
  ),
  predictors = c(
    'X ~ 1@0',
    'G ~ X'
  )
)

mymodel <- rblimp(
  data = dat,
  latent = 'X Y',
  nominal = 'G',
  model = syntax,
  simple = 'X | G',
  seed = 90291,
  burn = 5000,
  iter = 5000
)
output(mymodel)
simple_plot(Y ~ X | G.2, mymodel)

################################################################################
# Conduct Small Simulation With N Test Value to Verify Results
################################################################################

N <- 250
ng <- N * probs

results <- NULL
for(r in 1:25){
  
  g1dat <- cbind(1,rmvnorm(ng[1], as.vector(moments_G1$mean), as.matrix(moments_G1$covariance)))
  g2dat <- cbind(2,rmvnorm(ng[2], as.vector(moments_G2$mean), as.matrix(moments_G2$covariance)))
  dat <- as.data.frame(rbind(g1dat,g2dat))
  names(dat) <- c('G',paste0('X',1:num_indicators),paste0('Y',1:num_indicators))
  
  syntax <- list(
    structural_model = c(
      'Y ~ 1@beta0 G.2@beta_g2 X@beta_X G.2*X@beta_XG2',
      'Y ~~ Y@resvar'
    ),
    X_measurement = c(
      'X ~~ X@1',
      paste0('X -> X1@lx1 X2:X', num_indicators, collapse = '')
    ),
    Y_measurement = c(
      paste0('Y -> Y1:Y', num_indicators, collapse = ''),
      'Y1 ~ 1@0'
    ),
    predictors = c(
      'X ~ 1@0',
      'G ~ X'
    )
  )
  
  mymodel <- rblimp(
    data = dat,
    latent = 'X Y',
    nominal = 'G',
    model = syntax,
    parameters = c(
      paste0('beta0_bi = beta0 - ', round(moderated_params$beta_0,6), collapse = ''),
      paste0('beta_g2_bi = beta_g2 - ', round(moderated_params$beta_G2,6), collapse = ''),
      paste0('beta_X_bi = beta_X - ', round(moderated_params$beta_X,6), collapse = ''),
      paste0('beta_XG2_bi = beta_XG2 - ', round(moderated_params$beta_XG2,6), collapse = ''),
      paste0('resvar_bi = resvar - ', round(moderated_params$residual_variance,6), collapse = '')),
    seed = 90291,
    burn = 5000,
    iter = 5000
  )
  output(mymodel)
  
  results_tmp <- cbind(r, 1:5,mymodel@estimates[(nrow(mymodel@estimates)-4):nrow(mymodel@estimates),1:2])
  results <- rbind(results,results_tmp)
  
}

colnames(results) <- c('rep','param','bias','sd')
results <- as.data.frame(results)
agg_results <- aggregate(cbind(bias, sd) ~ param, data = results, FUN = mean, na.rm = TRUE)
agg_results$bias_sd <- agg_results$bias / agg_results$sd
print(agg_results)

