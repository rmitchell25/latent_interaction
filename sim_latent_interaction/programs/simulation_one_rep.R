library(mvtnorm)
library(gridExtra)
library(grid)
library(psych)
library(dplyr)
library(lavaan)
library(rblimp)
options(scipen=999)

runfromshell <- T


# import arguments from unix shell script or specify manually in else{}
if(runfromshell){
  runvars <- commandArgs(T)
  runoncluster <- as.integer(runvars[1])
  dirname <- (runvars[2])
  filename <- (runvars[3])
  
  group_prob <- as.numeric(runvars[4])
  rsq_prod <- as.numeric(runvars[5])
  N <- as.numeric(runvars[6])
  loading <- as.numeric(runvars[7])
  n_items <- as.numeric(runvars[8])
  rep <- as.numeric(runvars[9])
  seed <- as.numeric(runvars[10])
}else{
  runoncluster <- 0
  dirname <- "~/Desktop/sim_latent_interaction"
  filename<- "g1prod03N350load5item6"
  
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


results <- matrix(c(group_prob,rsq_prod,N,loading,n_items, rep, seed), nrow = 1)
colnames(results) <- c("group_prob","rsq_prod","N","loading","n_items", "rep", "seed")

write.table(results,paste0(dirname,'/results/',filename,'.dat'),row.names = F,col.names = F)
