library(mvtnorm)
library(gridExtra)
library(grid)
library(psych)
library(dplyr)
library(lavaan)
library(rblimp)
library(MplusAutomation)
options(scipen=999)
set.seed(75080)


dirname <- "C:/Users/remus/OneDrive/Documents/GitHub/latent_interaction/sim_latent_interaction"
load(file = paste0(dirname,"/misc/parameter_values.rda"))


corr_Xs <- .20
rsq_baseline <- .13
rsq_prod <- c(0,0.03,0.07)
group_probs3 <- rbind(c(.34, .33, .33),
                      c(.40, .40, .20),
                      c(.60, .20, .20))
group_probs2 <- rbind(c(.50,.50),
                      c(.60,.40),
                      c(.80,.20))
# sample_size <- c(seq(100,400, by = 50),500,1000)
loading_size <- c(.5,.8)
num_loadings <- c(6,12)
sample_size <- c(seq(100,400, by = 50), 500, 1000)


resultsfull <- NULL
counter <- 0

for (cat in 2:3) {
  for (group_prob in 1:3) {
    for (rsq_product in rsq_prod) {
      for (n_items in num_loadings) {
        for (load in loading_size) {
          for (N in sample_size) {
            
            # Load parameters file
            load(file = paste0(dirname,"/misc/parameter_values.rda"))
            
            name <- paste0("cat",cat,"_prob",group_prob,"_rsq",rsq_product,"_items",
                           n_items,"_loading",load)
            
            if(cat == 2){
              probs <- group_probs2[group_prob,]
              bin <- T
            } else {
              probs <- group_probs3[group_prob,]
              bin <- F
            }
            
            # Covariance matrix and mean vector per group
            moments_G1 <- parameter_values[[name]]$G1
            moments_G2 <- parameter_values[[name]]$G2
            if(bin == F){moments_G3 <- parameter_values[[name]]$G3}
            
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
            
            
            if (exists("blimp_model")){
              if(bin == T){
                psr_check <- max(blimp_model@psr[20,], na.rm = T)
                neff_check <- min(blimp_model@estimates[1:5,7], na.rm = T)
              } else {
                psr_check <- max(blimp_model@psr[20,], na.rm = T)
                neff_check <- min(blimp_model@estimates[1:7,7], na.rm = T)
              }
              
              results <- matrix(999,nrow = 1,ncol=6)
              results[,1]<- cat
              results[,2]<-group_prob
              results[,3]<- rsq_product
              results[,4]<-n_items
              results[,5]<- load 
              results[,6]<- N 
              
              results<-cbind(results,psr_check,neff_check)
              results<-as.data.frame(results)
            }else{
              results[,1:8]<- -999
            }
            
            colnames(results)<- c("categories","group_prob","rsq_product","n_items", 
                                  "loading","N","psr_max", "n_eff_greater")
            resultsfull<- rbind(resultsfull, results)
            
            counter <- counter + 1
            print(counter)
            
          }
        }
      }
    }
  }
}
