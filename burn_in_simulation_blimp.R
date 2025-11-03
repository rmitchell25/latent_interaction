library(mvtnorm) 
library(dplyr) 
library(lavaan)
library(rblimp)
library(car)
options(scipen=999)


# dirname <- "~/Documents/GitHub/latent_interaction/sim_latent_interaction"
dirname <- "C:/Users/remus/OneDrive/Documents/GitHub/latent_interaction/sim_latent_interaction"

cats <- c(2,3)
group_probs <- seq(1:3)
rsq_prod <- 0.03    # c(0, 0.03, 0.07)
Ns <- c(100,1000)      # seq(100,400, by = 50), 500, 1000
loadings <- c(.5,.8)   # .5 or .8
n_item <- c(6,12)    # 6 or 12
reps <- 100
seed <- 91030

# burn <- seq(5000,30000, by = 5000)
# iteration <- seq(5000,30000, by = 5000)

burnin <- 15000
iter <- 15000

resultsfull<- NULL

#for (burnin in burn) {
 # for (iter in iteration) {
    for (rep in 1:reps) {
      for (cat in cats) {
        for (group_prob in group_probs){
          for (N in Ns) {
            for (loading in loadings) {
              for (n_items in n_item) {
                # Select group probabilities based on category condition
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
                
                
                # Pull population parameters
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
                
                
                # Generate data
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
                
                
                # Fit in Blimp and save results
                syntax <- list(
                  structural_model = c('Y ~ 1 G X G*X@int'),
                  X_measurement = c('X ~~ X@1',
                                    paste0('X -> X1@lx1 X2:X', n_items)),
                  Y_measurement = c(paste0('Y -> Y1:Y', n_items),
                                    'Y1 ~ 1@0'),
                  predictors = c('X ~ 1@0',
                                 'G ~ X')
                )
                
                no_cvg_blimp <- 0 
                
                blimp_model <- tryCatch({
                  rblimp(
                    data = dat,
                    latent = 'X Y',
                    nominal = 'G',
                    model = syntax,
                    simple = 'X | G',
                    seed = seed,
                    burn = burnin,
                    iter = iter,
                    output = "default wald pvalue",
                    waldtest = list('int = 0')
                  )
                }, error = function(e) {
                  message("Blimp model failed to converge: ", e$message)
                  no_cvg_blimp <<- 1
                  return(NULL)
                })
                
                
                results <- matrix(NA, nrow = 1, ncol = 10)
                
                results[,1] <- cat
                results[,2] <- group_prob
                results[,3] <- N
                results[,4] <- loading
                results[,5] <- n_items
                results[,6] <- rep
                results[,7] <- burnin
                results[,8] <- iter
                
                if (no_cvg_blimp == 0){
                  psr_check <- max(blimp_model@psr[20,], na.rm = T)
                  neff_check <- min(blimp_model@estimates[,7], na.rm = T)
                  results[,9] <- psr_check
                  results[,10] <- neff_check
                } else {
                  results[,9] <- NA
                  results[,10] <- NA
                }
                
                resultsfull<- rbind(resultsfull,results)
                
                
              }
            }
          }
        }
      }
    }
#  }
#}

resultsfull <- as.data.frame(resultsfull)
colnames(resultsfull) <- c("cat","group_prob","N","loading","n_items","rep",
                           "psr_check","neff_check")



