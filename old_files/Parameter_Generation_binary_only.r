# remotes::update_packages('rblimp')

library(Matrix)
library(mvtnorm)
library(rblimp)

# Simulation Conditions ----
corr_Xs <- .20
rsq_baseline <- .13
rsq_prod <- c(0,0.03,0.07)
# group_probs3 <- rbind(c(.34, .33, .33), 
#                       c(.46, .27, .27),
#                       c(.74, .13, .13))
# group_probs2 <- rbind(c(.50,.50),
#                       c(.60,.40),
#                       c(.80,.20))
group_probs <- rbind(c(.50,.50),
                      c(.80,.20))
# sample_size <- c(seq(100,400, by = 50),500,1000)
loading_size <- c(.5,.8)
num_loadings <- c(6,12)


# Step 1: Define Group-Specific Means of X to Induce Correlation ----

solve_factor_means <- function(probs, target_corr = corr_Xs) {
  
  p1 <- probs[1]
  p2 <- probs[2]
  
  stopifnot(abs(p1 + p2 - 1) < 1e-6)
  
  objective <- function(mu2) {
    mu1         <- -(p2 * mu2) / p1
    mu_bar      <- 0
    var_between <- p1 * mu1^2 + p2 * mu2^2
    var_within  <- 1 - var_between
    var_X       <- var_between + var_within
    cov_X_G2    <- p2 * mu2
    var_G2      <- p2 * (1 - p2)
    corr_G2     <- cov_X_G2 / sqrt(var_X * var_G2)
    (corr_G2 - target_corr)^2
  }
  
  mu2 <- optimize(objective, interval = c(-3, 3))$minimum
  mu1 <- -(p2 * mu2) / p1
  
  c(mu1 = mu1, mu2 = mu2)
}


# Step 2: Solve for Within-Group Variance of X to Achieve Var(X) = 1 ----

solve_factor_variance <- function(mu_X, probs, target_var = 1) {
  
  stopifnot(length(mu_X) == 2,
            length(probs) == 2,
            abs(sum(probs) - 1) < 1e-6)
  
  mu_bar      <- sum(probs * mu_X)
  var_between <- sum(probs * (mu_X - mu_bar)^2)
  var_within  <- target_var - var_between
  
  var_group <- rep(var_within, 2)
  names(var_group) <- paste0("sigma_sq_g", 1:2)
  
  list(
    mu_bar          = mu_bar,
    var_between     = var_between,
    var_within      = var_within,
    group_variances = var_group
  )
}


# Step 3: Solve for Group-Specific Slopes and Intercepts (Structural Parameters) ----

solve_group_parameters <- function(mu_X, var_X, probs,
                                   rsq_baseline, rsq_product) {
  
  stopifnot(length(mu_X) == 2,
            length(var_X) == 2,
            length(probs) == 2,
            abs(sum(probs) - 1) < 1e-6)
  
  mu_X  <- as.numeric(mu_X)
  var_X <- as.numeric(var_X)
  probs <- as.numeric(probs)
  
  # --- baseline slope ---
  beta_baseline <- sqrt(rsq_baseline)
  
  # --- match TOTAL explained variance in Y (same as cluster code) ---
  objective <- function(delta2) {
    
    deltas <- c(0, delta2)   # asymmetric (matches cluster version)
    betas  <- beta_baseline + deltas
    
    # implied Y means
    y_means <- betas * mu_X
    mean_y  <- sum(probs * y_means)
    
    # total explained variance in Y
    explained_var <- sum(
      probs * (betas^2 * var_X + (y_means - mean_y)^2)
    )
    
    (explained_var - (rsq_baseline + rsq_product))^2
  }
  
  # --- solve for interaction ---
  delta2_opt <- optimize(objective, interval = c(0, 5))$minimum
  
  deltas <- c(0, delta2_opt)
  betas  <- beta_baseline + deltas
  
  # --- intercepts (same as cluster version) ---
  alpha_shared <- -sum(probs * betas * mu_X)
  alphas       <- rep(alpha_shared, 2)
  
  # --- compute TRUE residual variance (not fixed) ---
  y_means <- betas * mu_X + alphas
  var_Y   <- sum(
    probs * (betas^2 * var_X + (y_means - sum(probs * y_means))^2)
  )
  
  residual_var <- 1 - var_Y
  sigmas_sq    <- rep(residual_var, 2)
  
  names(betas) <- names(alphas) <- names(sigmas_sq) <- paste0("G", 1:2)
  
  list(
    intercepts         = alphas,
    slopes             = betas,
    residual_variances = sigmas_sq
  )
}


# Step 4: Convert to Pooled Moderated Regression Parameters ----

group_to_moderated <- function(mu_X, var_X, probs, group_params) {
  
  intercepts <- group_params$intercepts
  slopes     <- group_params$slopes
  resid_vars <- group_params$residual_variances
  
  beta_0  <- intercepts[1]
  beta_X  <- slopes[1]
  beta_G  <- intercepts[2] - intercepts[1]
  beta_XG <- slopes[2]     - slopes[1]
  
  residual_variance <- sum(probs * resid_vars)
  
  list(
    beta_0            = beta_0,
    beta_G            = beta_G,
    beta_X            = beta_X,
    beta_XG           = beta_XG,
    residual_variance = residual_variance
  )
}

# Step 5: Measurement Model for Latent Factors ----


model_implied_moments <- function(mu_X, var_X, alpha, beta, residual_var_Y, 
                                  n_X, n_Y, stanload, binary = F) {
  
  if(binary == T){
    
    # Measurement parameters
    lambda_X <- stanload
    theta_X  <- 1 - lambda_X^2
    lambda_Y <- 1.00
    theta_Y  <- 1 / stanload^2 - 1
    
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
    
  } else{
    
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
  
}


# Step 6: Create sum score moderation parameters ----

generate_sumscore_parameters <- function(
    group_params,
    factor_means,
    factor_vars,
    probs,
    n_indicators = num_indicators,
    loading = loading
) {
  
  # --- Extract group-specific values ---
  beta  <- group_params$slopes
  alpha <- group_params$intercepts
  resid <- group_params$residual_variances
  
  mu_X  <- factor_means
  var_X <- factor_vars
  
  # --- Measurement model constants ---
  lambda_X <- loading
  lambda_Y <- 1
  
  # Sum score scaling
  k <- n_indicators
  
  # --- Expected sum score means ---
  mu_X_sum <- k * lambda_X * mu_X
  mu_Y_sum <- k * lambda_Y * (alpha + beta * mu_X)
  
  # --- Variances of sum scores ---
  var_X_sum <- k^2 * lambda_X^2 * var_X + k * (1 - lambda_X^2)
  
  var_Y_sum <- k^2 * lambda_Y^2 * (beta^2 * var_X + resid) +
    k * (1 / lambda_X^2 - 1)
  
  # --- Covariance ---
  cov_XY_sum <- k^2 * lambda_X * lambda_Y * beta * var_X
  
  # --- Convert to regression within each group ---
  b_X_g <- cov_XY_sum / var_X_sum
  b_0_g <- mu_Y_sum - b_X_g * mu_X_sum
  
  # --- Convert to moderated regression ---
  beta_0  <- b_0_g[1]
  beta_X  <- b_X_g[1]
  beta_G  <- b_0_g[2] - b_0_g[1]
  beta_XG <- b_X_g[2] - b_X_g[1]
  
  # --- Residual variance (pooled) ---
  resid_g <- var_Y_sum - (b_X_g^2 * var_X_sum)
  resid_pooled <- sum(probs * resid_g)
  
  return(list(
    beta_0            = beta_0,
    beta_X            = beta_X,
    beta_G            = beta_G,
    beta_XG           = beta_XG,
    residual_variance = resid_pooled
  ))
}


# Derive parameters for each set of conditions ----
parameter_values <- list()

for (p in 1:2) {
  
  probs <- group_probs[p,]
  
  for (rsq_product in rsq_prod) {
    
    # --- Step 1 ---
    factor_means <- solve_factor_means(probs = probs)
    
    # --- Step 2 ---
    factor_vars <- solve_factor_variance(
      mu_X = factor_means,
      probs = probs
    )$group_variances
    
    # --- Step 3 ---
    group_params <- solve_group_parameters(
      mu_X         = factor_means,
      var_X        = factor_vars,
      probs        = probs,
      rsq_baseline = rsq_baseline,
      rsq_product  = rsq_product
    )
    
    # --- Step 4 ---
    moderated_params <- group_to_moderated(
      mu_X         = factor_means,
      var_X        = factor_vars,
      probs        = probs,
      group_params = group_params
    )
    
    for (items in num_loadings) {
      for (load in loading_size) {
        
        # --- Step 5: Moments ---
        moments_G1 <- model_implied_moments(
          mu_X           = factor_means["mu1"],
          var_X          = factor_vars["sigma_sq_g1"],
          alpha          = group_params$intercepts["G1"],
          beta           = group_params$slopes["G1"],
          residual_var_Y = group_params$residual_variances["G1"],
          n_X            = items,
          n_Y            = items,
          stanload       = load
        )
        
        moments_G2 <- model_implied_moments(
          mu_X           = factor_means["mu2"],
          var_X          = factor_vars["sigma_sq_g2"],
          alpha          = group_params$intercepts["G2"],
          beta           = group_params$slopes["G2"],
          residual_var_Y = group_params$residual_variances["G2"],
          n_X            = items,
          n_Y            = items,
          stanload       = load
        )
        
        # --- Step 6: UPDATED SUM SCORE PARAMETERS ---
        sum_mod <- generate_sumscore_parameters(
          group_params  = group_params,
          factor_means  = factor_means,
          factor_vars   = factor_vars,
          probs         = probs,
          n_indicators  = items,
          loading       = load
        )
        
        # --- Naming ---
        name <- paste0(
          "prob", p,
          "_rsq", rsq_product,
          "_items", items,
          "_loading", load
        )
        
        # --- Store ---
        parameter_values[[name]] <- list(
          G1 = moments_G1,
          G2 = moments_G2,
          group_parameters     = group_params,
          mod_parameters       = moderated_params,
          sum_score_parameters = sum_mod
        )
      }
    }
  }
}

save(parameter_values, file = "parameter_values.rda")

