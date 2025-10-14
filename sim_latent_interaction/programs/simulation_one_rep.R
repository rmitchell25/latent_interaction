# Load packages ----

library(mvtnorm) 
library(dplyr) 
library(lavaan)
library(rblimp)
library(car)
options(scipen=999)

runfromshell <- T


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
  
  cat <- 3
  group_prob <- 1  # 1:3
  rsq_prod <- 0.07    # 0, 0.03, 0.07
  N <- 100      # seq(100,400, by = 50), 500, 1000
  loading <- .5   # .5 or .8
  n_items <- 6    # 6 or 12
  rep <- 1
  seed <- 75080
}

# set paths
if(runoncluster == 1){setwd(dirname)} else if(runoncluster == 0){setwd(paste0(dirname))}



# Select group probabilities based on category condition ----

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



# Pull population parameters ----

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



# Fit in Blimp and save results ----

syntax <- list(
  structural_model = c('Y ~ 1 G X G*X@int'),
  X_measurement = c('X ~~ X@1',
                   paste0('X -> X1@lx1 X2:X', n_items)),
  Y_measurement = c(paste0('Y -> Y1:Y', n_items),
                  'Y1 ~ 1@0'),
  predictors = c('X ~ 1@0',
                  'G ~ X')
)

blimp_model <- rblimp(
  data = dat,
  latent = 'X Y',
  nominal = 'G',
  model = syntax,
  simple = 'X | G',
  seed = seed,
  burn = 10000,
  iter = 10000,
  output = "default wald pvalue",
  waldtest = list('int = 0')
)
# output(blimp_model)




# Fit multigroup version with Lavaan ----

if(n_items == 6){
  if(bin == T){
    model <- '
        # structural X
        X ~ c(0,intX2)*1      # fix the first latent mean to 0 and free the other
        X ~~ 1*X              # set within-group variance to 1 (equal across groups)

        # structural Y
        Y ~ c(intY1,intY2)*1 + c(slope1, slope2)*X  # group-specific intercepts map to dummy code effects
        Y ~~ c(resvarY,resvarY)*Y                   # pooled residual variance

        # measurement: loadings & intercepts equal across groups
        X =~ NA*X1 + X2 + X3 + X4 + X5 + X6   # free first loading because var(Xw) = 1
        X1 ~ 0*1                                     # fix first intercept to 0
        Y =~ Y1 + Y2 + Y3 + Y4 + Y5 + Y6      # fix first loading because var(Y) is estimated
        Y1 ~ 0*1                              # fix first intercept to 0

        # conversions to MR parameters
        B0 := intY1
        G_slope := intY2 - intY1
        X_slope := slope1
        GbyX_slope := slope2 - slope1
        Y_resvar := resvarY'

  } else {
    model <- '
        # structural X
        X ~ c(0,intX2,intX3)*1      # fix the first latent meam to 0 and free the others
        X ~~ 1*X                    # set within-group variance to 1

        # structural Y
        Y ~ c(intY1,intY2,intY3)*1 + c(slope1, slope2, slope3)*X    # group-specific intercepts map to the dummy code effects
        Y ~~ c(resvarY,resvarY,resvarY)*Y                           # pooled residual variance

        # measurement models: all intercepts and residual variances estimated with equality constraints across groups
        X =~ NA*X1 + X2 + X3 + X4 + X5 + X6          # free first loading because var(Xw) = 1
        X1 ~ 0*1                                     # fix first intercept to 0
        Y =~ Y1 + Y2 + Y3 + Y4 + Y5 + Y6             # fix first loading because var(Y) is estimated
        Y1 ~ 0*1                                     # fix first intercept to 0

        B0 := intY1
        G2_slope := intY2 - intY1
        G3_slope := intY3 - intY1
        X_slope := slope1
        G2byX_slope := slope2 - slope1
        G3byX_slope := slope3 - slope1
        Y_resvar := resvarY'
  }

} else if(n_items == 12){
  if(bin == T){
    model <- '
        # structural X
        X ~ c(0,intX2)*1
        X ~~ 1*X

        # structural Y
        Y ~ c(intY1,intY2)*1 + c(slope1, slope2)*X
        Y ~~ c(resvarY,resvarY)*Y

        # measurement: loadings & intercepts equal across groups
        X =~ NA*X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + X11 + X12
        X1 ~ 0*1
        Y =~ Y1 + Y2 + Y3 + Y4 + Y5 + Y6 + Y7 + Y8 + Y9 + Y10 + Y11 + Y12
        Y1 ~ 0*1

        # conversions to MR parameters
        B0 := intY1
        G_slope := intY2 - intY1
        X_slope := slope1
        GbyX_slope := slope2 - slope1
        Y_resvar := resvarY'
  } else {
    model <- '
        # structural X
        X ~ c(0,intX2,intX3)*1
        X ~~ 1*X

        # structural Y
        Y ~ c(intY1,intY2,intY3)*1 + c(slope1, slope2, slope3)*X
        Y ~~ c(resvarY,resvarY,resvarY)*Y

        # measurement models: all intercepts and residual variances estimated with equality constraints across groups
        X =~ NA*X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + X11 + X12
        X1 ~ 0*1
        Y =~ Y1 + Y2 + Y3 + Y4 + Y5 + Y6 + Y7 + Y8 + Y9 + Y10 + Y11 + Y12
        Y1 ~ 0*1

        B0 := intY1
        G2_slope := intY2 - intY1
        G3_slope := intY3 - intY1
        X_slope := slope1
        G2byX_slope := slope2 - slope1
        G3byX_slope := slope3 - slope1
        Y_resvar := resvarY'
  }
}

# fit multiple group model
fit <- sem(model, data = dat, group = "G",
           group.equal = c("loadings", "residuals","intercepts"))

# summarize
mg_model <-  summary(fit, standardized = T)
mg_est <- parameterEstimates(fit, ci = TRUE, level = 0.95)

if(!bin){
  lavJoint <- lavTestWald(fit, constraints = "G2byX_slope == G3byX_slope == 0")
}






# Fit version with sum score ----

# Calculate SD of sum scores
Sigma1 <- moments_G1$covariance
Sigma2 <- moments_G2$covariance
if(!bin){
  Sigma3 <- moments_G3$covariance
  Sigma_pooled <- ((ng[1] - 1) * Sigma1 + (ng[2] - 1) * Sigma2 + (ng[3] - 1) * Sigma3) / 
    (ng[1] + ng[2] + ng[3] - 3)
} else{
  Sigma_pooled <- ((ng[1] - 1) * Sigma1 + (ng[2] - 1) * Sigma2) / (ng[1] + ng[2] - 2)
}

# Compute pooled covariance matrix

ones_X <- c(rep(1, n_items), rep(0, n_items))  
ones_Y <- c(rep(0, n_items), rep(1, n_items)) 

var_X <- t(ones_X) %*% Sigma_pooled %*% ones_X
var_Y <- t(ones_Y) %*% Sigma_pooled %*% ones_Y

sd_X <- as.numeric(sqrt(var_X))
sd_Y <- as.numeric(sqrt(var_Y))


# add sum scores to data set after scaling them
X_sum <- (dat %>% select(starts_with("X")) %>% rowSums())
Y_sum <- (dat %>% select(starts_with("Y")) %>% rowSums())
X_sum <- X_sum/sd_X
Y_sum <- Y_sum*((mod_parameters[["residual_variance"]]/(1-rsq_prod))/(sd_Y))

dat <- cbind(dat, X_sum, Y_sum)

# run analysis
sum_score_model <- lm(Y_sum ~ factor(G) + X_sum + factor(G)*X_sum, data = dat)
sum_model <- summary(lm(Y_sum ~ factor(G) + X_sum + factor(G)*X_sum, data = dat))

if (!bin){
  sumJoint <- linearHypothesis(sum_score_model, c("factor(G)2:X_sum = 0", "factor(G)3:X_sum  = 0"))
}





# Results: Power, Type 2 Error, Relative Bias, MSE, CI Coverage ----

# Set up data frame differently based on number of categories
if (bin == F){
  results <- as.data.frame(matrix(999, nrow = 21, ncol = 22))
  param.id <- c("beta_0","beta_G2","beta_G3","beta_X","beta_XG2","beta_XG3","res.var")
} else {
  results <- as.data.frame(matrix(999, nrow = 15, ncol = 22))
  param.id <- c("beta_0","beta_G2","beta_X","beta_XG2","res.var")
}
colnames(results) <- c("categories", "group_prob", "rsq_prod", "N", "loading",
                       "n_items", "model_type", "param.id", "est", "se","true",
                       "bias", "rel.bias", "squared_bias","pval", "sig", "ci.cov",
                       "ci.zero", "psr", "neff","joint.stat","joint.p")

# Save info on set of conditions
results[,1] <- cat
results[,2] <- group_prob
results[,3] <- rsq_prod
results[,4] <- N
results[,5] <- loading
results[,6] <- n_items
results[,7] <- c(rep(1,(nrow(results)/3)),rep(2,(nrow(results)/3)),rep(3,(nrow(results)/3)))
results[,8] <- rep(param.id,3)


# Save estimates, sd, true parameters, and then calculate bias
blimp.est <- as.numeric(c(blimp_model@estimates[2:(nrow(results)/3),1],blimp_model@estimates[1,1]))
mg.est <- as.numeric(mg_model[["pe"]][["est"]][(((length(mg_model[["pe"]][["est"]])+1)-(nrow(results)/3))):length(mg_model[["pe"]][["est"]])])
sum.est <- as.numeric(c(sum_model[["coefficients"]][,1], sum_model[["sigma"]]^2))
results[,9] <- c(blimp.est, mg.est, sum.est)
results[,10] <- c(as.numeric(c(blimp_model@estimates[2:(nrow(results)/3),2],blimp_model@estimates[1,2])),
                  as.numeric(mg_model[["pe"]][["se"]][(((length(mg_model[["pe"]][["se"]])+1)-(nrow(results)/3))):length(mg_model[["pe"]][["se"]])]),
                  as.numeric(sum_model[["coefficients"]][,2]),NA)

results[,10] <- c(rep(as.numeric(mod_parameters),3))

results[,12] <- results[,9] - results[,10]
results[,13] <- results[,11]/results[,10]
results[,14] <- results[,11]^2


# P-values and Significance
pvalues_blimp <- as.numeric(c(blimp_model@estimates[2:(nrow(results)/3),8], blimp_model@estimates[1,8]))
pvalues_mg <- as.numeric(mg_model[["pe"]][["pvalue"]][(((length(mg_model[["pe"]][["pvalue"]])+1)-(nrow(results)/3))):length(mg_model[["pe"]][["pvalue"]])])
pvalues_sum <- as.numeric(c(sum_model[["coefficients"]][,4], NA))
pvals <- c(pvalues_blimp, pvalues_mg, pvalues_sum)

sig_mod<- rep(0,nrow(results))
sig_mod[pvals < .05]<-1

results[,15] <- pvals
results[,16] <- sig_mod


# CI Coverage - True value within CI
CIL_blimp <- as.numeric(c(blimp_model@estimates[2:(nrow(results)/3),3],blimp_model@estimates[1,3]))
CIU_blimp <- as.numeric(c(blimp_model@estimates[2:(nrow(results)/3),4],blimp_model@estimates[1,4]))
CIL_mg <- as.numeric(mg_est[(1 + nrow(mg_est) - (nrow(results)/3)):nrow(mg_est),11])
CIU_mg <- as.numeric(mg_est[(1 + nrow(mg_est) - (nrow(results)/3)):nrow(mg_est),12])
CIL_sum <- as.numeric(c(confint(sum_score_model, level = 0.95)[,1],NA))
CIU_sum <- as.numeric(c(confint(sum_score_model, level = 0.95)[,2],NA))

CI_lower <- c(CIL_blimp ,CIL_mg, CIL_sum)
CI_upper <- c(CIU_blimp, CIU_mg, CIU_sum)

CI_cov <- rep(0,nrow(results))
CI_cov[CI_lower < as.numeric(results[,9]) & CI_upper > as.numeric(results[,9])] <-1

results[,17] <- CI_cov


# Zero within CI
CI_zero <- rep(0,nrow(results))
CI_zero[CI_lower < 0 & CI_upper > 0] <- 1
results[,18] <- CI_zero


# Save PSR and N-effective values from Blimp model
results[1:(nrow(results)/3),19] <- as.numeric(c(blimp_model@psr[20,2:(nrow(results)/3)],blimp_model@psr[20,1]))
results[(1+nrow(results)/3):nrow(results),19] <- results[(1+nrow(results)/3):nrow(results),20] <- NA
results[1:(nrow(results)/3),20] <- as.numeric(c(blimp_model@estimates[2:(nrow(results)/3),6], blimp_model@estimates[1,6]))


# Save results from joint tests if cat = 3
if(!bin){
  results[,21] <- c(rep(blimp_model@waldtest[["statistic"]],(nrow(results)/3)),
                    rep(lavJoint[["stat"]],((nrow(results)/3))),
                    rep(sumJoint$F[2],((nrow(results)/3))))
  results[,22] <- c(rep(blimp_model@waldtest[["probability"]],(nrow(results)/3)),
                    rep(lavJoint[["p.value"]],(nrow(results)/3)),
                    rep(sumJoint$`Pr(>F)`[2],(nrow(results)/3)))
} else {
  results[,21] <- results[,22] <- NA
}




# Test for shell script
# results <- c(cat, group_prob, rsq_prod, N, loading, n_items, rep, seed)

write.table(results,paste0(dirname,'/results/',filename,'.dat'),row.names = F,col.names = F)

