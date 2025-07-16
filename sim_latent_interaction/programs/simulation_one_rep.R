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
  
  categories <- as.numeric(runvars[4])
  group_prob <- as.numeric(runvars[5])
  rsq_prod <- as.numeric(runvars[6])
  N <- as.numeric(runvars[7])
  loading <- as.numeric(runvars[8])
  n_items <- as.numeric(runvars[9])
  rep <- as.numeric(runvars[10])
  seed <- as.numeric(runvars[11])
}else{
  runoncluster <- 0
  dirname <- "~/Documents/GitHub/latent_interaction/sim_latent_interaction"
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


############################################################
# Generate data
############################################################

ng <- N * probs

g1dat <- cbind(1,rmvnorm(ng[1], as.vector(moments_G1$mean), as.matrix(moments_G1$covariance)))
g2dat <- cbind(2,rmvnorm(ng[2], as.vector(moments_G2$mean), as.matrix(moments_G2$covariance)))

if (bin == F){
  g3dat <- cbind(3,rmvnorm(ng[3], as.vector(moments_G3$mean), as.matrix(moments_G3$covariance)))
  dat <- as.data.frame(rbind(g1dat,g2dat,g3dat))
  names(dat) <- c('G',paste0('X',1:num_indicators),paste0('Y',1:num_indicators))
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
  Y ~ G X_eta G*X_eta;
  measurement:
  X_eta ~~ X_eta@1;
  X_eta -> X1@lo1 X2:X6;
  Y -> Y1@1 Y2:Y6;
  predictors:
  G ~ X_eta',
  parameters = 'lo1 ~ trunc(0,inf)'
)

blimp_est <- blimp_model@estimates[1:9,]

############################################################
# Fit multigroup version with MPLus Automation 
############################################################



############################################################
# Fit version with sum score
############################################################



############################################################
# Save results
############################################################

results <- matrix(c(group_prob,rsq_prod,N,loading,n_items, rep, seed), nrow = 1)
colnames(results) <- c("group_prob","rsq_prod","N","loading","n_items", "rep", "seed")

write.table(results,paste0(dirname,'/results/',filename,'.dat'),row.names = F,col.names = F)


