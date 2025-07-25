library(mvtnorm)
library(ggplot2)
library(gridExtra)
library(grid)
library(psych)
library(dplyr)
library(lavaan)
library(rblimp)

options(scipen = 999)
set.seed(90291)

# pre-sets
muY <- 0
varY <- 1
corr_Xs <- .20
rsq_baseline <- .13
XtoD_split <- 1
simp_slopes <- 'diffuse' # interaction pattern is concentrated vs. diffuse


# conditions
rsq_prod <- c(0,0.03,0.07)
group_probs <- rbind(c(.34, .33, .33),
                     c(.40, .40, .20),
                     c(.60, .20, .20))

binary_probs <- rbind(c(.50,.50),
                      c(.60,.40),
                      c(.80,.20))


############################################################
# function to solve for latent response variable means
############################################################

solve_mu <- function(mu, group_probs) {
  
  sigma <- matrix(c(1,.5,.5,1), nrow = 2)
  
  # compute probability c = 1 (reference) as area below 0
  lower <- rep(-Inf, (length(mu) + 1) - 1); upper <- rep(0, (length(mu) + 1) - 1)
  prob1 <- pmvnorm(lower, upper, mu, corr = NULL, sigma)
  
  # transform mean vector and covariance matrix
  trans_matrix <- matrix(0, ncol = (length(mu) + 1) - 1, nrow = (length(mu) + 1) - 1)
  trans_matrix[,1] <- -1; trans_matrix[2,2] <- 1
  mu_trans <- as.vector(trans_matrix %*% (mu))
  sigma_trans <- trans_matrix %*% sigma %*% t(trans_matrix)
  
  # compute probability c = 2 as area below 0 in the transformed space
  lower <- rep(-Inf, (length(mu) + 1) - 1); upper <- rep(0, (length(mu) + 1) - 1)
  prob2 <- pmvnorm(lower, upper, mu_trans, corr = NULL, sigma_trans)
  
  # transformation matrix, mean vector, covariance matrix for c = 3
  trans_matrix <- matrix(0, ncol = (length(mu) + 1) - 1, nrow = (length(mu) + 1) - 1)
  trans_matrix[,2] <- -1; trans_matrix[2,1] <- 1
  mu_trans <- as.vector(trans_matrix %*% (mu))
  sigma_trans <- trans_matrix %*% sigma %*% t(trans_matrix)
  
  # compute probability c = 3 as area below 0 in the transformed space
  lower <- rep(-Inf, (length(mu) + 1) - 1); upper <- rep(0, (length(mu) + 1) - 1)
  prob3 <- pmvnorm(lower, upper, mu_trans, corr = NULL, sigma_trans)
  
  # loss function: sum of square differences between target and current values
  return(sum((c(prob1, prob2, prob3) - group_probs)^2))
}

############################################################
# # function to obtain predictor covariance matrix
############################################################

solve_cov <- function(muX){
  
  sigma <- matrix(corr_Xs, nrow = 3, ncol = 3)
  diag(sigma) <- 1; sigma[2,1] <- sigma[1,2] <- .5
  
  coding <- 'dummy'
  N <- 25000000
  Xs <- rmvnorm(N,muX,sigma)
  nomX <- ifelse(Xs[,1] < 0 & Xs[,2] < 0, 1, 999)
  nomX <- ifelse(Xs[,1] > 0 & Xs[,1] > Xs[,2], 2, nomX)
  nomX <- ifelse(Xs[,2] > 0 & Xs[,2] > Xs[,1], 3, nomX)
  if(coding == 'dummy'){
    C2 <- ifelse(nomX == '2', 1, 0) 
    C3 <- ifelse(nomX == '3', 1, 0)
  }
  if(coding == 'effect'){
    C2 <- ifelse(nomX == "1", -1, ifelse(nomX == "2", 1, 0))
    C3 <- ifelse(nomX == "1", -1, ifelse(nomX == "3", 1, 0))
  }
  
  code_data <- cbind(C2,C3,Xs[,3],C2*Xs[,3],C3*Xs[,3])
  stats_CXs <- rbind(colMeans(code_data),cov(code_data))
  
  return(stats_CXs)
}

############################################################
# function to solve for slopes
############################################################

solve_slopes <- function(slopes, statsX, rsq_baseline, rsq_prod){
  
  # slopes <- c(0.13073646,0.13208590,0.18246931,0.06743022,0.13484283)
  # statsX = stats_CXs[[1]]
  # muX = statsX[1,]
  # covX = statsX[-1,]
  # rsq_baseline = rsq_baseline[1]
  # rsq_prod = rsq_prod[1]
  
  # check whether slopes are positive
  if(any(slopes <= 0)) return(1e10)
  
  # extract predictor means and covariance matrix
  muX <- statsX[1, ]
  covX <- statsX[-1, ]
  num_vars <- nrow(covX) + 1
  
  # construct beta and psi matrices
  beta_mat <- psi_mat <- matrix(0, nrow = num_vars, ncol = num_vars)
  beta_mat[num_vars, 1:(num_vars - 1)] <- slopes
  psi_mat[1:(num_vars - 1), 1:(num_vars - 1)] <- covX
  
  # compute the residual variance and make sure it is positive
  resvar <- varY - (slopes %*% covX %*% slopes)
  if(!is.finite(resvar) || resvar <= 0){
    return(1e10)
  }
  psi_mat[num_vars, num_vars] <- resvar
  
  # solve for the overall mean vector and covariance matrix
  inv_mat <- solve(diag(num_vars) - beta_mat)
  mu_all <- inv_mat %*% c(muX, muY)
  cov_all <- inv_mat %*% psi_mat %*% t(inv_mat)
  
  # assign row/column names
  row.names(mu_all) <- c(paste0("predictor ", 1:(num_vars-1)), "outcome")
  colnames(cov_all) <- row.names(cov_all) <- c(paste0("predictor ", 1:(num_vars-1)), "outcome")
  
  # target r-square values
  rsq_tot_target <- rsq_baseline + rsq_prod
  rsq_X_target <- XtoD_split * rsq_baseline
  rsq_codes_target <- rsq_baseline - rsq_X_target
  rsq_cha_target <- rsq_prod
  
  # current R-square values from covariance matrix
  rsq_tot_current <- (cov_all[1:(num_vars-1), num_vars] %*% 
                        solve(cov_all[1:(num_vars-1), 1:(num_vars-1)]) %*% 
                        cov_all[num_vars, 1:(num_vars-1)]) / varY
  rsq_codes_X_current <- (cov_all[1:3, num_vars] %*% 
                            solve(cov_all[1:3, 1:3]) %*% 
                            cov_all[num_vars, 1:3]) / varY
  rsq_codes_current <- (cov_all[1:2, num_vars] %*% 
                          solve(cov_all[1:2, 1:2]) %*% 
                          cov_all[num_vars, 1:2]) / varY
  rsq_X_current <- rsq_codes_X_current - rsq_codes_current
  rsq_cha_current <- rsq_tot_current - rsq_codes_X_current
  
  # interaction slope difference
  int_ratio <- slopes[5] / slopes[4]
  if(simp_slopes == "concentrated"){pop_ratio <- 1}
  if(simp_slopes == "diffuse"){pop_ratio <- 2}
  
  # form target and current vectors
  target_vec <- c(rsq_tot_target, rsq_X_target, rsq_codes_target, rsq_cha_target, pop_ratio)
  current_vec <- c(rsq_tot_current, rsq_X_current, rsq_codes_current, rsq_cha_current, int_ratio)
  
  # loss function: sum of squared differences.
  loss <- sum((current_vec - target_vec)^2)
  return(loss)
}

############################################################
# solve for latent predictor means
# get predictor means covariance matrices
############################################################

lat_means <- list()
stats_CXs <- list()

# solve predictor means and covariance matrices for each group size condition
for(g in 1:3){
  print(paste0('computing large-N predictor covariance matrix for group ', g))
  # latent response variable means
  lat_means[[g]] <- c(optim(c(0,0), solve_mu, group_probs = group_probs[g,], method = "BFGS")$par,0)
  # code variable means and covariance matrix
  stats_CXs[[g]] <- solve_cov(lat_means[[g]]) # cbind(C2,C3,Xs[,3],C2*Xs[,3],C3*Xs[,3])
}

############################################################
# solve for population parameters
############################################################

optim_results <- NULL
starting <-  0.115

for(g in 1:3){
  for(p in 1:length(rsq_prod)){
    condition_name <- paste0('group probs = (', paste(group_probs[g,], collapse = ", "),'), product rsq = ',rsq_prod[p], " varY = ", varY, " , start = ",starting)
    
    # solve for slopes
    result <- optim(rep(starting, 5), solve_slopes, 
                    statsX = stats_CXs[[g]],
                    rsq_baseline = rsq_baseline, 
                    rsq_prod = rsq_prod[p],
                    method = "BFGS",
                    control = list(maxit = 5000, reltol = 1e-8))
    optim_results <- rbind(optim_results, c(g,rsq_baseline,rsq_prod[p],varY,result$convergence,result$value,round(result$par,3)))
    
  }
}
colnames(optim_results) <- c('group','rsq_base','rsq_prod','varY','conv_status','loss_fun','b1','b2','b3','b4','b5')
optim_results <- as.data.frame(optim_results)


############################################################
# produce population parameters across conditions - 3 group
############################################################

structural_params <- NULL

for (group in 1:3) {
  for (prod in rsq_prod) {
    # population slopes
    slopes_opt <- optim_results[optim_results$group == group & optim_results$rsq_prod == prod, 7:11]
    slopes_opt <- as.numeric(slopes_opt)
    
    # compute population residual variance
    sigma_pop <- stats_CXs[[group]][2:nrow(stats_CXs[[group]]), ]
    explained_var <- t(matrix(slopes_opt, ncol = 1)) %*% sigma_pop %*% matrix(slopes_opt, ncol = 1)
    resid_var <- 1 - explained_var
    
    # population intercept
    X_means_pop <- stats_CXs[[group]][1,]
    Y_mu <- 0
    B0 <- Y_mu - slopes_opt %*% X_means_pop
    
    params <- c(group,prod,B0,slopes_opt,resid_var)
    
    structural_params <- rbind(structural_params,params)
  }
}
colnames(structural_params) <- c("group","prod_rsq","b0", "b1","b2","b3","b4","b5","res_var")
structural_params <- as.data.frame(structural_params)


save(structural_params, file = "structural_params.rda")


latent_means <- matrix(data = NA, nrow = 3, ncol = 3)
latent_means[1,] <- lat_means[[1]]
latent_means[2,] <- lat_means[[2]]
latent_means[3,] <- lat_means[[3]]
save(latent_means, file ="latent_means.rda")


############################################################
# produce population parameters across conditions - 2 group
############################################################

solve_mu_binary <- function(mu, group_probs) {
  # Define 2x2 covariance matrix
  sigma <- matrix(c(1, .5, .5, 1), nrow = 2)
  
  # Probability of being in group 1 (e.g., reference group)
  # This corresponds to both latent responses being below zero
  lower <- rep(-Inf, 2)
  upper <- rep(0, 2)
  prob1 <- pmvnorm(lower, upper, mean = mu, sigma = sigma)
  
  # Probability of group 2 is just 1 - prob1
  prob2 <- 1 - prob1
  
  # Loss function: sum of squared differences between target and current values
  return(sum((c(prob1, prob2) - group_probs)^2))
}

solve_cov_binary <- function(muX) {
  
  # 2 latent classification variables: C* and X
  sigma <- matrix(c(1, .5, .5, 1), nrow = 2)
  
  N <- 25000000
  Xs <- rmvnorm(N, muX, sigma)
  
  # Binary classification: 1 if both latent variables below 0, else 2
  nomX <- ifelse(Xs[,1] < 0, 1, 2)
  
  # Choose coding scheme
  coding <- 'dummy'  # or 'effect'
  
  if (coding == 'dummy') {
    C <- ifelse(nomX == 2, 1, 0)
  } else if (coding == 'effect') {
    C <- ifelse(nomX == 1, -1, 1)
  }
  
  # Construct variables for interaction
  X <- Xs[,2]
  CX <- C * X
  
  # Create data matrix: [C, X, C*X]
  code_data <- cbind(C, X, CX)
  
  # Return means and covariance matrix
  stats_CXs <- rbind(colMeans(code_data), cov(code_data))
  return(stats_CXs)
}

solve_slopes_binary <- function(slopes, statsX, rsq_baseline, rsq_prod) {
  
  # Check for non-positive slopes
  if (any(slopes <= 0)) return(1e10)
  
  # Extract predictor means and covariance matrix
  muX <- statsX[1, ]
  covX <- statsX[-1, ]
  num_vars <- nrow(covX) + 1
  
  # Construct beta and psi matrices
  beta_mat <- psi_mat <- matrix(0, nrow = num_vars, ncol = num_vars)
  beta_mat[num_vars, 1:(num_vars - 1)] <- slopes
  psi_mat[1:(num_vars - 1), 1:(num_vars - 1)] <- covX
  
  # Compute residual variance and check positivity
  resvar <- varY - (slopes %*% covX %*% slopes)
  if (!is.finite(resvar) || resvar <= 0) return(1e10)
  psi_mat[num_vars, num_vars] <- resvar
  
  # Solve for overall means and covariance matrix
  inv_mat <- solve(diag(num_vars) - beta_mat)
  mu_all <- inv_mat %*% c(muX, muY)
  cov_all <- inv_mat %*% psi_mat %*% t(inv_mat)
  
  # Assign row/col names
  row.names(mu_all) <- c("C", "X", "C*X", "Y")
  colnames(cov_all) <- row.names(cov_all) <- c("C", "X", "C*X", "Y")
  
  # Target R² values
  rsq_tot_target <- rsq_baseline + rsq_prod
  rsq_X_target <- XtoD_split * rsq_baseline
  rsq_codes_target <- rsq_baseline - rsq_X_target
  rsq_cha_target <- rsq_prod
  
  # Current R² values based on covariance matrix
  # Total R²
  rsq_tot_current <- (cov_all[1:3, 4] %*% 
                        solve(cov_all[1:3, 1:3]) %*% 
                        cov_all[4, 1:3]) / varY
  # R² for C and X
  rsq_codes_X_current <- (cov_all[1:2, 4] %*% 
                            solve(cov_all[1:2, 1:2]) %*% 
                            cov_all[4, 1:2]) / varY
  # R² for C alone
  rsq_codes_current <- (cov_all[1, 4] %*% 
                          solve(cov_all[1, 1]) %*% 
                          cov_all[4, 1]) / varY
  # Unique R² for X
  rsq_X_current <- rsq_codes_X_current - rsq_codes_current
  # Unique R² for interaction
  rsq_cha_current <- rsq_tot_current - rsq_codes_X_current
  
  # Interaction slope ratio: beta_int / beta_X
  int_ratio <- slopes[3] / slopes[2]
  if (simp_slopes == "concentrated") {
    pop_ratio <- 1
  } else if (simp_slopes == "diffuse") {
    pop_ratio <- 2
  }
  
  # Target and current vectors
  target_vec <- c(rsq_tot_target, rsq_X_target, rsq_codes_target, rsq_cha_target, pop_ratio)
  current_vec <- c(rsq_tot_current, rsq_X_current, rsq_codes_current, rsq_cha_current, int_ratio)
  
  # Return loss (squared error)
  loss <- sum((current_vec - target_vec)^2)
  return(loss)
}



lat_bin_means <- list()
stats_CXs_bin <- list()
for(g in 1:3){
  print(paste0('computing large-N predictor covariance matrix for group ', g))
  # latent response variable means
  lat_bin_means[[g]] <- c(optim(0, solve_mu_binary, group_probs = binary_probs[g,], method = "BFGS")$par,0)
  # code variable means and covariance matrix
  stats_CXs_bin[[g]] <- solve_cov_binary(lat_bin_means[[g]]) # cbind(C2,C3,Xs[,3],C2*Xs[,3],C3*Xs[,3])
}


optim_results_bin <- NULL
starting <-  0.115

for(g in 1:3){
  for(p in 1:length(rsq_prod)){
    # solve for slopes
    result <- optim(rep(starting, 3), solve_slopes_binary, 
                    statsX = stats_CXs_bin[[g]],
                    rsq_baseline = rsq_baseline, 
                    rsq_prod = rsq_prod[p],
                    method = "BFGS",
                    control = list(maxit = 5000, reltol = 1e-8))
    
    optim_results_bin <- rbind(optim_results_bin, c(g,rsq_baseline,rsq_prod[p],varY,result$convergence,result$value,round(result$par,3)))
    
  }
}
colnames(optim_results_bin) <- c('group','rsq_base','rsq_prod','varY','conv_status','loss_fun','b1','b2','b3')
optim_results_bin <- as.data.frame(optim_results_bin)


structural_params_bin <- NULL

for (group in 1:3) {
  for (prod in rsq_prod) {
    # population slopes
    slopes_opt <- optim_results_bin[optim_results_bin$group == group & optim_results_bin$rsq_prod == prod, 7:9]
    slopes_opt <- as.numeric(slopes_opt)
    
    # compute population residual variance
    sigma_pop <- stats_CXs_bin[[group]][2:nrow(stats_CXs_bin[[group]]), ]
    explained_var <- t(matrix(slopes_opt, ncol = 1)) %*% sigma_pop %*% matrix(slopes_opt, ncol = 1)
    resid_var <- 1 - explained_var
    
    # population intercept
    X_means_pop <- stats_CXs_bin[[group]][1,]
    Y_mu <- 0
    B0 <- Y_mu - slopes_opt %*% X_means_pop
    
    params <- c(group,prod,B0,slopes_opt,resid_var)
    
    structural_params_bin <- rbind(structural_params_bin,params)
  }
}
colnames(structural_params_bin) <- c("group","prod_rsq","b0", "b1","b2","b3","res_var")
structural_params_bin <- as.data.frame(structural_params_bin)


save(structural_params_bin, file = "structural_params_bin.rda")


lat_bin_means <- matrix(data = NA, nrow = 3, ncol = 2)
lat_bin_means[1,] <- lat_bin_means[[1]]
lat_bin_means[2,] <- lat_bin_means[[2]]
lat_bin_means[3,] <- lat_bin_means[[3]]
save(lat_bin_means, file ="lat_bin_means.rda")
