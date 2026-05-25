# Load packages ----

library(mvtnorm) 
library(dplyr) 
library(lavaan)
library(rblimp)
library(car)
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
  # dirname <- "C:/Users/remus/OneDrive/Documents/GitHub/latent_interaction/sim_latent_interaction"
  dirname <- "~/Documents/GitHub/latent_interaction/sim_latent_interaction"
  filename<- "g1prod03N350load5item6"
  
  cat <- 2
  group_prob <- 1  
  rsq_prod <- 0.07  
  N <- 10000
  loading <- .8 
  n_items <- 12  
  rep <- 1
  seed <- 92834191
}

# Select group probabilities based on category condition ----

group_probs3 <- rbind(c(.34, .33, .33),
                      c(.40, .40, .20),
                      c(.60, .20, .20))
group_probs2 <- rbind(c(.50,.50),
                      # c(.60,.40),
                      c(.80,.20))
if(cat == 2){
  probs <- group_probs2[group_prob,]
  bin <- T
  numparams <- 5
} else {
  probs <- group_probs3[group_prob,]
  bin <- F
  numparams <- 7
}



# Pull population parameters ----

# Load parameters file
load(file = paste0(dirname,"/misc/parameter_values.rda"))

# TESTS
# print(paste("Loading from:", paste0(dirname,"/misc/parameter_values.rda")))
# print(file.exists(paste0(dirname,"/misc/parameter_values.rda")))
# 
# print("Objects in environment after load:")
# print(ls())
# 
# print("parameter_values exists?")
# print(exists("parameter_values"))
# 
# print("Names inside parameter_values:")
# print(names(parameter_values)[1:5])
# print(paste0("dirname = ", dirname))

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



# Generate data ----

ng <- round(N * probs)

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

 dat$G1 <- ifelse(dat$G == 2, 1, 0)
# dat$G2 <- ifelse(dat$G == 3, 1, 0)
 
setwd("/Users/remus/Documents/GitHub/latent_interaction/mplus_test")

write.table(dat,'mplus_test_data.txt', row.names = F,
            col.names = F)
