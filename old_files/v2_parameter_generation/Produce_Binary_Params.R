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


################################################################################
# Derive parameters for each set of conditions
################################################################################

binary_parameters <- NULL

for (p in 1:3) {
  probs <- group_probs2[p,]
  for (rsq_product in rsq_prod) {

    factor_means <- solve_factor_means_2group(probs = probs, target_corr = corr_Xs)
    
    factor_vars <- solve_factor_variance_2group(mu_X = factor_means, 
                                                probs = probs)$group_variances
    
    group_params <- solve_group_parameters_2group(mu_X = factor_means, 
                                                  var_X = factor_vars, 
                                                  probs = probs,
                                                  rsq_baseline = rsq_baseline,
                                                  rsq_product = rsq_product)
    
    moderated_params <- group_to_moderated_2group(
      mu_X = factor_means,
      var_X = factor_vars,
      probs = probs,
      group_params = group_params
    )
    
    id <- c(paste0("b",0:3),"res_var")
    params <- cbind.data.frame(id,as.numeric(moderated_params), p, rsq_product)
    
    binary_parameters <- rbind(binary_parameters, params)
  }
}

colnames(binary_parameters) <- c("param_id","param","prob","rsq_p")
binary_parameters <- as.data.frame(binary_parameters)


save(binary_parameters, file = "binary_structural_parameters.rda")
