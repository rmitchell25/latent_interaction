# remotes::update_packages('rblimp')

library(Matrix)
library(mvtnorm)
library(rblimp)

# Simulation Conditions ----
corr_Xs <- .20
rsq_baseline <- .13
rsq_prod <- c(0,0.03,0.07)
group_probs3 <- rbind(c(.34, .33, .33), 
                      c(.46, .27, .27),
                      c(.74, .13, .13))
group_probs2 <- rbind(c(.50,.50),
                      c(.60,.40),
                      c(.80,.20))
sample_size <- c(seq(100,400, by = 50),500,1000)
loading_size <- c(.5,.8)
num_loadings <- c(6,12)


# Step 1: Define Group-Specific Means of X to Induce Correlation ----

solve_factor_means <- function(probs, target_corr = corr_Xs, binary = F) {
  
  stopifnot(length(probs) ==2 & binary==T | length(probs) ==3 & binary==F)
  stopifnot(sum(probs)==1)
  
  if(binary == T){
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
      return((corr_G2 - target_corr)^2)
    }
    
    mu2 <- optimize(objective, interval = c(-3, 3))$minimum
    mu1 <- -(p2 * mu2) / p1
    
    # validation
    var_between <- p1 * mu1^2 + p2 * mu2^2
    if (var_between >= 1) {
      stop("Between-group variance >= 1 → invalid (negative within variance).")
    }
    
    var_X    <- 1
    cov_X_G2 <- p2 * mu2
    var_G2   <- p2 * (1 - p2)
    corr_G2  <- cov_X_G2 / sqrt(var_X * var_G2)
    
    if (abs(corr_G2 - target_corr) > 1e-3) {
      warning(sprintf("Achieved corr (%.4f) differs from target (%.4f)", corr_G2, target_corr))
    }
    
    return(c(mu1 = mu1, mu2 = mu2))
    
  } else{
    p1 <- probs[1]
    p2 <- probs[2]
    p3 <- probs[3]
    
    stopifnot(abs(p1 + p2 + p3 - 1) < 1e-6)
    
    # Define objective function for optimizer
    objective <- function(mu2) {
      mu1 <- 0
      mu3 <- -(p1 * mu1 + p2 * mu2) / p3  # ensures E[X] = 0
      mu_bar <- 0 
      
      # Between-group variance
      var_between <- p1 * (mu1 - mu_bar)^2 + p2 * (mu2 - mu_bar)^2 + p3 * (mu3 - mu_bar)^2
      var_within <- 1 - var_between
      var_X <- var_between + var_within  # should = 1
      
      # Covariances
      cov_X_G2 <- p2 * (mu2 - mu_bar)
      cov_X_G3 <- p3 * (mu3 - mu_bar)
      
      # Variances of dummy codes
      var_G2 <- p2 * (1 - p2)
      var_G3 <- p3 * (1 - p3)
      
      # Correlations
      corr_G2 <- cov_X_G2 / sqrt(var_X * var_G2)
      corr_G3 <- cov_X_G3 / sqrt(var_X * var_G3)
      
      # Minimize squared deviation from +target and -target
      return((corr_G2 - target_corr)^2 + (corr_G3 + target_corr)^2)
    }
    
    # Optimize over a wider interval to avoid getting stuck at zero
    opt <- optimize(objective, interval = c(0, 3))
    
    mu1 <- 0
    mu2 <- opt$minimum
    mu3 <- -(p1 * mu1 + p2 * mu2) / p3
    
    # validation
    var_between <- p1 * mu1^2 + p2 * mu2^2 + p3 * mu3^2
    if (var_between >= 1) {
      stop("Between-group variance >= 1 → invalid (negative within variance).")
    }
    
    var_X <- 1
    
    cov_X_G2 <- p2 * mu2
    cov_X_G3 <- p3 * mu3
    
    var_G2 <- p2 * (1 - p2)
    var_G3 <- p3 * (1 - p3)
    
    corr_G2 <- cov_X_G2 / sqrt(var_X * var_G2)
    corr_G3 <- cov_X_G3 / sqrt(var_X * var_G3)
    
    if (abs(corr_G2 - target_corr) > 1e-3 ||
        abs(corr_G3 + target_corr) > 1e-3) {
      warning(sprintf(
        "Achieved corr: G2=%.4f (target %.4f), G3=%.4f (target %.4f)",
        corr_G2, target_corr, corr_G3, -target_corr))
    }
    
    return(c(mu1 = mu1, mu2 = mu2, mu3 = mu3))
  }
  
}


# Step 2: Solve for Within-Group Variance of X to Achieve Var(X) = 1 ----

solve_factor_variance <- function(mu_X, probs, target_var = 1, binary = F) {
  
  if(binary == T){
    stopifnot(length(mu_X) == 2,
              length(probs) == 2,
              abs(sum(probs) - 1) < 1e-6)
    
    mu_bar <- sum(probs * mu_X)
    var_between <- sum(probs * (mu_X - mu_bar)^2)
    var_within  <- target_var - var_between
    
    if(var_within <= 0) {
      stop("Cannot achieve Var(X)=1 with these group means. var_within = ", var_within)
    }
    
    var_group <- rep(var_within, 2)
    names(var_group) <- paste0("sigma_sq_g", 1:2)
    
    return(list(mu_bar = mu_bar, var_between = var_between,var_within = var_within,
                group_variances = var_group))
    
  } else{
    
    stopifnot(length(mu_X) == 3, length(probs) == 3, abs(sum(probs) - 1) < 1e-6)
    
    # Pooled mean of X
    mu_bar <- sum(probs * mu_X)
    
    # Between-group variance component
    var_between <- sum(probs * (mu_X - mu_bar)^2)
    
    # Solve for within-group variance to ensure total variance = target_var
    var_within <- target_var - var_between
    
    if(var_within <= 0) {
      stop("Cannot achieve Var(X)=1 with these group means. var_within = ", var_within)
    }
    
    # All group variances set equal to this value
    var_group <- rep(var_within, 3)
    names(var_group) <- paste0("sigma_sq_g", 1:3)
    
    return(list(
      mu_bar = mu_bar,
      var_between = var_between,
      var_within = var_within,
      group_variances = var_group
    ))
  }
}



# Step 3: Solve for Group-Specific Slopes and Intercepts (Structural Parameters) ----

solve_group_parameters <- function(mu_X, var_X, probs, rsq_baseline, rsq_product,
                                   binary = F) {
  
  if(binary == T){
    stopifnot(length(mu_X) == 2, length(var_X) == 2, length(probs) == 2,
              abs(sum(probs) - 1) < 1e-6)
    stopifnot(all(var_X > 0), rsq_baseline >= 0, rsq_product >= 0, 
              rsq_baseline + rsq_product <= 1)
    
    mu_X  <- as.numeric(mu_X)
    var_X <- as.numeric(var_X)
    probs <- as.numeric(probs)
    
    # baseline slope
    beta_baseline <- sqrt(rsq_baseline)
    
    # match TOTAL explained variance in Y (same as cluster code) ---
    objective <- function(delta2) {
      
      deltas <- c(0, delta2)   # asymmetric (matches cluster version)
      betas  <- beta_baseline + deltas
      
      # implied Y means
      y_means <- betas * mu_X
      mean_y  <- sum(probs * y_means)
      
      # total explained variance in Y
      explained_var <- sum(probs * (betas^2 * var_X + (y_means - mean_y)^2))
      
      (explained_var - (rsq_baseline + rsq_product))^2
    }
    
    # solve for interaction
    delta2_opt <- optimize(objective, interval = c(0, 5))$minimum
    
    deltas <- c(0, delta2_opt)
    betas  <- beta_baseline + deltas
    
    # intercepts (same as cluster version)
    alpha_shared <- -sum(probs * betas * mu_X)
    alphas       <- rep(alpha_shared, 2)
    
    # compute TRUE residual variance (not fixed) 
    y_means <- betas * mu_X + alphas
    var_predicted <- sum(probs * (betas^2 * var_X + (y_means - sum(probs * y_means))^2))
    
    if(var_predicted >= 1) {
      stop("Predicted variance >= 1. Cannot achieve R² = ", 
           rsq_baseline + rsq_product, " with Var(Y)=1")
    }
    
    if (abs(var_predicted - (rsq_baseline + rsq_product)) > 1e-3) {
      stop(sprintf("Optimizer failed: var_predicted=%.6f, target=%.6f",
                   var_predicted, rsq_baseline + rsq_product))
    }
    
    residual_var <- 1 - var_predicted
    sigmas_sq    <- rep(residual_var, 2)
    
    names(betas) <- names(alphas) <- names(sigmas_sq) <- paste0("G", 1:2)
    
    return(list(
      intercepts = alphas,
      slopes = betas,
      residual_variances = sigmas_sq,
      predicted_variance = var_predicted,
      target_r_squared = rsq_baseline + rsq_product))
    
  } else{
    
    stopifnot(length(mu_X) == 3, length(var_X) == 3, length(probs) == 3, abs(sum(probs) - 1) < 1e-6)
    stopifnot(all(var_X > 0), rsq_baseline >= 0, rsq_product >= 0, 
              rsq_baseline + rsq_product <= 1)
    
    mu_X <- as.numeric(mu_X)
    var_X <- as.numeric(var_X)
    probs <- as.numeric(probs)
    
    # 1. Solve for baseline slope that gives exact rsq_baseline
    beta_baseline <- sqrt(rsq_baseline)
    
    
    # 2. Optimize delta2 under constraint delta3 = 2 * delta2
    objective <- function(delta2) {
      delta3 <- 2 * delta2
      deltas <- c(0, delta2, delta3)
      betas <- beta_baseline + deltas
      
      # Compute shared intercept to ensure E[Y] = 0
      alpha <- -sum(probs * betas * mu_X)
      
      # Compute predicted means per group
      pred_means <- alpha + betas * mu_X
      grand_mean_y <- sum(probs * pred_means)
      
      # Var(Yhat) = E[Var(Yhat|G)] + Var(E[Yhat|G])
      within_component <- sum(probs * betas^2 * var_X)
      between_component <- sum(probs * (pred_means - grand_mean_y)^2)
      var_predicted <- within_component + between_component
      
      target_rsq <- rsq_baseline + rsq_product
      return((var_predicted - target_rsq)^2)
    }
    
    delta2_opt <- optimize(objective, interval = c(0, 5))$minimum
    delta3_opt <- 2 * delta2_opt
    deltas <- c(0, delta2_opt, delta3_opt)
    betas <- beta_baseline + deltas
    
    # 3. Shared intercept to ensure E[Y] = 0
    alpha_shared <- -sum(probs * betas * mu_X)
    alphas <- rep(alpha_shared, 3)
    
    # 4. Residual variance to ensure Var(Y) = 1
    pred_means <- alpha_shared + betas * mu_X
    grand_mean_y <- sum(probs * pred_means)
    within_component <- sum(probs * betas^2 * var_X)
    between_component <- sum(probs * (pred_means - grand_mean_y)^2)
    var_predicted <- within_component + between_component
    
    if(var_predicted >= 1) {
      stop("Predicted variance >= 1. Cannot achieve R² = ", 
           rsq_baseline + rsq_product, " with Var(Y)=1")
    }
    
    if (abs(var_predicted - (rsq_baseline + rsq_product)) > 1e-3) {
      stop(sprintf("Optimizer failed: var_predicted=%.6f, target=%.6f",
                   var_predicted, rsq_baseline + rsq_product))
    }
    
    residual_var <- 1 - var_predicted
    sigmas_sq <- rep(residual_var, 3)
    
    names(betas) <- names(alphas) <- names(sigmas_sq) <- paste0("G", 1:3)
    
    return(list(
      intercepts = alphas,
      slopes = betas,
      residual_variances = sigmas_sq,
      predicted_variance = var_predicted,
      target_r_squared = rsq_baseline + rsq_product))
  }
}




# Step 4: Convert to Pooled Moderated Regression Parameters ----

group_to_moderated <- function(mu_X, probs, group_params, binary = F) {
  
  if(binary == T){
    
    stopifnot(length(mu_X) == 2, length(probs) == 2, length(group_params$intercepts) == 2)
    
    intercepts <- group_params$intercepts
    slopes     <- group_params$slopes
    resid_vars <- group_params$residual_variances
    
    beta_0  <- intercepts[1]
    beta_X  <- slopes[1]
    beta_G  <- intercepts[2] - intercepts[1]
    beta_XG <- slopes[2] - slopes[1]
    
    residual_variance <- sum(probs * resid_vars)
    
    return(list(
      beta_0 = beta_0,
      beta_G = beta_G,
      beta_X = beta_X,
      beta_XG = beta_XG,
      residual_variance = residual_variance))
    
  } else{
    
    stopifnot(length(mu_X) == 3, length(probs) == 3, length(group_params$intercepts) == 3)
    
    # Extract components from group_params
    intercepts <- group_params$intercepts
    slopes <- group_params$slopes
    resid_vars <- group_params$residual_variances
    
    # 1. Reference group: Group 1
    beta_0 <- intercepts[1]
    beta_X <- slopes[1]
    
    # 2. Dummy-coded differences
    beta_G2 <- intercepts[2] - intercepts[1]
    beta_G3 <- intercepts[3] - intercepts[1]
    
    beta_XG2 <- slopes[2] - slopes[1]
    beta_XG3 <- slopes[3] - slopes[1]
    
    # 3. Pooled residual variance
    residual_variance <- sum(probs * resid_vars)
    
    return(list(
      beta_0 = beta_0,
      beta_G2 = beta_G2,
      beta_G3 = beta_G3,
      beta_X = beta_X,
      beta_XG2 = beta_XG2,
      beta_XG3 = beta_XG3,
      residual_variance = residual_variance))
  }
}



# Step 5: Measurement Model for Latent Factors ----


model_implied_moments <- function(mu_X, var_X, alpha, beta, residual_var_Y, 
                                  n_X, n_Y, stanload) {
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
}


# Step 6: Create sum score moderation parameters ----

latent_to_sumscore_moderation <- function(beta_0, beta_G, beta_X, beta_XG, 
                                          res.var, mu_X, probs, var_X,
                                          n_items, lambda, binary) {
  
  var_X <- unname(var_X[1]) 
  
  # Within-group X_sum variance
  VxsW <- (n_items*lambda)^2 * var_X + n_items*(1 - lambda^2)
  
  # Within-group reliability of X_sum = rho_w
  one_mr <- n_items*(1 - lambda^2) / VxsW  # 1 - rho_w
  
  # Measurement error for Y
  theta_Y <- 1/lambda^2 - 1
  
  # Factor for adjusting slope parameters
  slope_factor <- (n_items^2 * lambda * var_X) / VxsW
  
  
  if (binary) {
    # per-group latent slopes and intercepts
    beta_g <- c(beta_X, beta_X + beta_XG)
    alpha_g <- c(beta_0, beta_0 + beta_G)
    # dilution term in group main effects
    dG <- n_items * one_mr * (beta_g[2]*mu_X[2] - beta_g[1]*mu_X[1])
    # residual variance
    Eb2 <- sum(probs * beta_g^2)
    res.var_sum <- n_items^2 * res.var + n_items * theta_Y + n_items^2 * var_X * one_mr * Eb2
   
     out <- list(
      beta_0  = n_items * alpha_g[1] + n_items * one_mr * beta_g[1] * mu_X[1],
      beta_G  = n_items * beta_G + dG,
      beta_X  = slope_factor * beta_X,
      beta_XG = slope_factor * beta_XG,
      res.var = res.var_sum
    )
  } else {
    beta_g <- c(beta_X, beta_X + beta_XG[1], beta_X + beta_XG[2])
    alpha_g <- c(beta_0, beta_0 + beta_G[1], beta_0 + beta_G[2])
    dG2 <- n_items * one_mr * (beta_g[2]*mu_X[2] - beta_g[1]*mu_X[1])
    dG3 <- n_items * one_mr * (beta_g[3]*mu_X[3] - beta_g[1]*mu_X[1])
    Eb2 <- sum(probs * beta_g^2)
    res.var_sum <- n_items^2 * res.var + n_items * theta_Y + n_items^2 * var_X * one_mr * Eb2
    
    out <- list(
      beta_0   = n_items * alpha_g[1] + n_items * one_mr * beta_g[1] * mu_X[1],
      beta_G2  = n_items * beta_G[1] + dG2,
      beta_G3  = n_items * beta_G[2] + dG3,
      beta_X   = slope_factor * beta_X,
      beta_XG2 = slope_factor * beta_XG[1],
      beta_XG3 = slope_factor * beta_XG[2],
      res.var  = res.var_sum)
  }
  return(out)
}



# Derive parameters for each set of conditions ----

parameter_values <- list()

for (categories in 2:3) {
  for (p in 1:3) {
    
    if(categories == 2){
      probs <- group_probs2[p,]
      bin <- T
    } else {
      probs <- group_probs3[p,]
      bin <- F
    }
    
    for (rsq_product in rsq_prod) {
      
      factor_means <- solve_factor_means(probs = probs, binary = bin)
      
      factor_vars <- solve_factor_variance(mu_X = factor_means, binary = bin,
                                           probs = probs)$group_variances
      
      group_params <- solve_group_parameters(mu_X = factor_means, 
                                             var_X = factor_vars, 
                                             probs = probs,
                                             rsq_baseline = rsq_baseline,
                                             rsq_product = rsq_product,
                                             binary = bin)
      
      moderated_params <- group_to_moderated(mu_X = factor_means,probs = probs,
                                             group_params = group_params,
                                             binary = bin)
      
      
      for (items in num_loadings) {
        for (load in loading_size) {
          
          # Group 1 mean vector and covariance matrix
          moments_G1 <- model_implied_moments(
            mu_X           = factor_means["mu1"],
            var_X          = factor_vars["sigma_sq_g1"],
            alpha          = group_params$intercepts["G1"],
            beta           = group_params$slopes      ["G1"],
            residual_var_Y = group_params$residual_variances["G1"],
            n_X            = items,
            n_Y            = items,
            stanload       = load)
          
          # Group 2 mean vector and covariance matrix
          moments_G2 <- model_implied_moments(
            mu_X           = factor_means["mu2"],
            var_X          = factor_vars["sigma_sq_g2"],
            alpha          = group_params$intercepts["G2"],
            beta           = group_params$slopes      ["G2"],
            residual_var_Y = group_params$residual_variances["G2"],
            n_X            = items,
            n_Y            = items,
            stanload       = load)
          
          if (bin == F){
            # Group 3 mean vector and covariance matrix
            moments_G3 <- model_implied_moments(
              mu_X           = factor_means["mu3"],
              var_X          = factor_vars["sigma_sq_g3"],
              alpha          = group_params$intercepts["G3"],
              beta           = group_params$slopes      ["G3"],
              residual_var_Y = group_params$residual_variances["G3"],
              n_X            = items,
              n_Y            = items,
              stanload       = load)
          } 
          
          
          
          if(bin == T) {
            slopes <- moderated_params$beta_G
            int_slopes <- moderated_params$beta_XG
          } else{
            slopes <- c(moderated_params$beta_G2, moderated_params$beta_G3)
            int_slopes <- c(moderated_params$beta_XG2, moderated_params$beta_XG3)
          }
          sum_mod <- latent_to_sumscore_moderation(beta_0 = moderated_params$beta_0,
                                                   beta_G = slopes,   # vector of group effects (b_G2, b_G3)
                                                   beta_X = moderated_params$beta_X,
                                                   beta_XG = int_slopes,  # vector of interaction effects
                                                   res.var = moderated_params$residual_variance,
                                                   mu_X = factor_means,
                                                   probs = probs,
                                                   var_X = factor_vars,
                                                   n_items = items, 
                                                   lambda = load,
                                                   binary = bin)
          
          
          name <- paste0("cat",categories,"_prob",p,"_rsq",rsq_product,"_items",
                         items,"_loading",load)
          
          parameter_values[[name]] <- list(
            G1 = moments_G1,
            G2 = moments_G2
          )
          if (bin == FALSE) {
            parameter_values[[name]]$G3 <- moments_G3
          }
          
          parameter_values[[name]]$group_parameters <- group_params
          
          parameter_values[[name]]$mod_parameters <- moderated_params
          
          parameter_values[[name]]$sum_score_parameters <- sum_mod
          
        }
      }
    }
  }
}


save(parameter_values, file = "parameter_values.rda")
