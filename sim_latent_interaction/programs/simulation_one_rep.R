library(mvtnorm)
library(gridExtra)
library(grid)
library(psych)
library(dplyr)
library(lavaan)
library(rblimp)
library(MplusAutomation)
options(scipen=999)

runfromshell <- F


# import arguments from unix shell script or specify manually in else{}
if(runfromshell){
  runvars <- commandArgs(T)
  runoncluster <- as.integer(runvars[1])
  dirname <- (runvars[2])
  filename <- (runvars[3])
  
  cat <- as.numeric(runvars[4])
  group_prob <- as.numeric(runvars[5])
  rsq_prod <- as.numeric(runvars[6])
  N <- as.numeric(runvars[7])
  loading <- as.numeric(runvars[8])
  n_items <- as.numeric(runvars[9])
  rep <- as.numeric(runvars[10])
  seed <- as.numeric(runvars[11])
}else{
  runoncluster <- 0
  dirname <- "C:/Users/remus/OneDrive/Documents/GitHub/latent_interaction/sim_latent_interaction"
  # dirname <- "~/Documents/GitHub/latent_interaction/sim_latent_interaction"
  filename<- "g1prod03N350load5item6"
  
  cat <- 2
  group_prob <- 1  # 1:3
  rsq_prod <- 0.03    # 0, 0.03, 0.07
  N <- 350      # seq(100,400, by = 50), 500, 1000
  loading <- .5   # .5 or .8
  n_items <- 6    # 6 or 12
  rep <- 1
  seed <- 75080
}

# set paths
if(runoncluster == 1){setwd(dirname)} else if(runoncluster == 0){setwd(paste0(dirname))}


############################################################
# Select group probabilities based on category condition
############################################################

group_probs3 <- rbind(c(.34, .33, .33),
                      c(.40, .40, .20),
                      c(.60, .20, .20))
group_probs2 <- rbind(c(.50,.50),
                      c(.60,.40),
                      c(.80,.20))
if(cat == 2){
  probs <- group_probs2[group_prob,]
  bin <- T
} else {
  probs <- group_probs3[group_prob,]
  bin <- F
}


############################################################
# Pull population parameters
############################################################

# Load parameters file
load(file = paste0(dirname,"/misc/parameter_values.rda"))

name <- paste0("cat",cat,"_prob",group_prob,"_rsq",rsq_prod,"_items",
               n_items,"_loading",loading)

# Moderation parameters
mod_parameters <- parameter_values[[name]]$mod_parameters 

# Multigroup parameters
group_parameters <- parameter_values[[name]]$group_parameters 

# Covariance matrix and mean vector per group
moments_G1 <- parameter_values[[name]]$G1
moments_G2 <- parameter_values[[name]]$G2
if(bin == F){moments_G3 <- parameter_values[[name]]$G3}

# sum score parameters
sum_score_params <- parameter_values[[name]]$sum_score_parameters 


############################################################
# Generate data
############################################################

ng <- N * probs

g1dat <- cbind(1,rmvnorm(ng[1], as.vector(moments_G1$mean), as.matrix(moments_G1$covariance)))
g2dat <- cbind(2,rmvnorm(ng[2], as.vector(moments_G2$mean), as.matrix(moments_G2$covariance)))

if (bin == F){
  g3dat <- cbind(3,rmvnorm(ng[3], as.vector(moments_G3$mean), as.matrix(moments_G3$covariance)))
  dat <- as.data.frame(rbind(g1dat,g2dat,g3dat))
  names(dat) <- c('G',paste0('X',1:n_items),paste0('Y',1:n_items))
} else {
  dat <- as.data.frame(rbind(g1dat,g2dat))
  names(dat) <- c('G',paste0('X',1:n_items),paste0('Y',1:n_items))
}


############################################################
# Fit in Blimp and save results
############################################################

blimp_model <- rblimp(
  data = dat,
  burn = 2000,
  iter = 5000,
  seed = 91030,
  nominal = 'G',
  latent = 'X_eta Y',
  model = '
  structural:
  Y ~ 1 G X_eta G*X_eta;
  measurement:
  X_eta ~~ X_eta@1;
  X_eta -> X1@lo1 X2:X6;
  Y -> Y1@1 Y2:Y6;
  predictors:
  G ~ X_eta',
  parameters = 'lo1 ~ trunc(0,inf)'
)


############################################################
# Fit multigroup version with MPLus Automation 
############################################################

# Just includes general structure... working on Mplus script

if (bin == T & n_items == 6) {
  model_syntax <- "
  MODEL:
    Y BY y1 y2 y3 y4 y5 y6;
    X BY x1* x2 x3 x4 x5 x6;

    Y ON X;

  model group1:
  
      [Y@0];   
      Y;     
      [X@0];   
      X@1;     
  
  model group2:
  
      [Y@0];
      Y;  
      [X];
      X@1;"
  mplus_model <- mplusObject(
    TITLE = "Multigroup Model;",
    VARIABLE = c('G',paste0('X',1:n_items),paste0('Y',1:n_items)),
    GROUPING = "G (1 = group1 2 = group2);",
    MODEL = model_syntax,
    rdata = dat,
    usevariables = c('G',paste0('X',1:n_items),paste0('Y',1:n_items)),
    OUTPUT = "SAMPSTAT STANDARDIZED TECH1 TECH4;"
  )
} else {
  
}

fit <- mplusModeler(
  mplus_model,
  modelout = "multigroup.inp",
  run = 1L  # set to 0 to write the input but not run
)


############################################################
# Fit version with sum score
############################################################

# add sum scores to data set
X_sum <- dat %>% select(starts_with("X")) %>% rowSums()
Y_sum <- dat %>% select(starts_with("Y")) %>% rowSums()
dat <- cbind(dat, X_sum, Y_sum)

# run analysis
sum_model <- summary(lm(Y_sum ~ G + X_sum + G*X_sum, data = dat))


############################################################
# Results: Power, Type 2 Error, Relative Bias, MSE, CI Coverage
############################################################

if (bin == F){
  results <- as.data.frame(matrix(999, nrow = 7, ncol = 15))
  param.id <- c("beta_0","beta_G2","beta_G3","beta_X","beta_XG2","beta_XG3","res.var")
  mod.est <- as.numeric(c(blimp_model@estimates[2:7,1],blimp_model@estimates[1,1]))
  
} else {
  results <- as.data.frame(matrix(999, nrow = 5, ncol = 15))
  param.id <- c("beta_0","beta_G2","beta_X","beta_XG2","res.var")
  mod.est <- as.numeric(c(blimp_model@estimates[2:5,1],blimp_model@estimates[1,1]))
  
}


# Moderation Significance and CI Coverage
pvalues1 <- as.numeric(c(blimp_model@estimates[2:nrow(results),6],blimp_model@estimates[1,6]))
sig_mod<- rep(0,nrow(results))
sig_mod[pvalues1 < .05]<-1

CI_lower <- as.numeric(c(blimp_model@estimates[2:nrow(results),3],blimp_model@estimates[1,3]))
CI_upper <- as.numeric(c(blimp_model@estimates[2:nrow(results),4],blimp_model@estimates[1,4]))
CI_cov_mod<- rep(0,nrow(results))
CI_cov_mod[CI_lower < as.numeric(mod_parameters) & CI_upper > as.numeric(mod_parameters)] <-1


# Save results
results[,1] <- cat
results[,2] <- group_prob
results[,3] <- rsq_prod
results[,4] <- N
results[,5] <- loading
results[,6] <- n_items
results[,7] <- 1            # 1 = moderation parameters
results[,8] <- param.id
results[,9] <- mod.est
results[,10] <- as.numeric(mod_parameters)
results[,11] <- results[,9] - results[,10]
results[,12] <- results[,11]/results[,10]
results[,13] <- results[,11]^2
results[,14] <- sig_mod
results[,15] <- CI_cov_mod
  
colnames(results) <- c("categories", "group_prob", "rsq_prod", "N", "loading", 
                       "n_items", "model_type", "param.id", "mod_est", "true_mod", 
                       "bias", "rel.bias", "squared_bias", "sig", "ci.cov")
write.table(results,paste0(dirname,'/results/',filename,'.dat'),row.names = F,col.names = F)

