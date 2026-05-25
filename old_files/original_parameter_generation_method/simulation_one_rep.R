library(mvtnorm)
library(gridExtra)
library(grid)
library(psych)
library(dplyr)
library(lavaan)
library(rblimp)
options(scipen=999)

runfromshell <- F


# import arguments from unix shell script or specify manually in else{}
if(runfromshell){
  runvars <- commandArgs(T)
  runoncluster <- as.integer(runvars[1])
  dirname <- (runvars[2])
  filename <- (runvars[3])
  
  group_probs <- as.numeric(runvars[4])
  rsq_prod <- as.numeric(runvars[5])
  N <- as.numeric(runvars[6])
  n_items <- as.numeric(runvars[7])
  loading <- as.numeric(runvars[8])
}else{
  runoncluster <- 0
  dirname <- "latent_interaction"
  filename<- "g1prod03N350load5item6"
  
  group_prob <- 1  # 1:3
  rsq_prod <- 0.03    # 0, 0.03, 0.07
  N <- 350      # seq(100,400, by = 50), 500, 1000
  loading <- .5   # .5 or .8
  n_items <- 6    # 6 or 12
}

# set paths
if(runoncluster == 1){setwd(dirname)} else if(runoncluster == 0){setwd(paste0(dirname))}


# pre-set values
muY <- 0
varY <- 1
corr_Xs <- .20
rsq_baseline <- .13
communality <- .75^2
error_var_Y <- (1/communality)-1


############################################################
# Produce population parameters
############################################################

# Pull structural parameters
load(file = "structural_params.rda")
params <- structural_params[structural_params$group == group_prob & structural_params$prod_rsq == rsq_prod, ]
colnames(params) <- rownames(params) <- NULL
params <- as.numeric(params[3:9])  # b0:b5, resvar

# Create loading vector of X and find error_var
lambdaX <- rep(loading, n_items)
error_var_X <- 1 - loading^2

# Create loading vector of X and find error_var  
lambdaY <- rep(1, n_items)
communality <- .75^2
error_var_Y <- (1/communality)-1

# covariance matrix of latent difference scores and latent X
load(file = "latent_means.rda")
muX <- latent_means[group_prob,]
sigma <- matrix(corr_Xs, nrow = 3, ncol = 3)
diag(sigma) <- 1; sigma[2,1] <- sigma[1,2] <- .5


############################################################
# Generate data 
############################################################

# generate Xs
Xs <- rmvnorm(N,muX,sigma)
nomX <- ifelse(Xs[,1] < 0 & Xs[,2] < 0, 1, 999)
nomX <- ifelse(Xs[,1] > 0 & Xs[,1] > Xs[,2], 2, nomX)
nomX <- ifelse(Xs[,2] > 0 & Xs[,2] > Xs[,1], 3, nomX)
# summarytools::freq(nomX)

# create dummy codes for Y generation
C2 <- ifelse(nomX == '2', 1, 0) 
C3 <- ifelse(nomX == '3', 1, 0)

# generate Y
code_data <- cbind(1,C2,C3,Xs[,3],C2*Xs[,3],C3*Xs[,3])
Y <- code_data %*% params[1:6] + rnorm(N, 0, sqrt(params[7]))

# generate each observed indicator of X = lambda * eta_X + error
Xindicator_matrix <- sapply(1:n_items, function(i) {
  lambdaX[i] * Xs[,3] + rnorm(N, mean = 0, sd = sqrt(error_var_X))
})

# generate each observed indicator of Y = lambda * eta_Y + error
Yindicator_matrix <- sapply(1:n_items, function(i) {
  lambdaY[i] * Y + rnorm(N, mean = 0, sd = sqrt(error_var_Y))
})


data <- cbind.data.frame(Y,nomX,Xindicator_matrix,Xs[,3],Yindicator_matrix)
colnames(data) <- c('Y','nomX', paste0("X", 1:n_items),'eta',paste0("Y", 1:n_items))

############################################################
# Fit in Blimp
############################################################

blimp_model <- rblimp(
  data = data,
  burn = 2000,
  iter = 5000,
  seed = 91030,
  nominal = 'nomX',
  latent = 'X_eta',
  model = '
  structural:
  Y ~ nomX X_eta nomX*X_eta;
  measurement:
  X_eta ~~ X_eta@1;
  X_eta -> X1@lo1 X2:X6;
  Y -> Y1@1 Y2:Y6;
  predictors:
  nomX ~ X_eta',
  parameters = 'lo1 ~ trunc(0,inf)'
)


############################################################
# Fit multigroup version in Laavan 
############################################################



############################################################
# Fit version with sum score
############################################################



############################################################
# Save results
############################################################







