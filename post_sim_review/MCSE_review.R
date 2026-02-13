options(scipen = 999)
library(dplyr)
library(ggplot2)

# Load data ----
colnames(rep_level) <- c("categories", "group_prob", "rsq_prod", "sample_size", "loading",
                         "n_items", "model_type", "param.id", "est", "se","true",
                         "bias", "rel.bias", "squared_bias","pval", "sig", "ci.cov",
                         "ci.zero","joint.stat","joint.p","no_cvg_check",
                         "psr.max","neff.min","ci.low","ci.up")


# POWER/TYPE I ERROR ----
power_mcse <- rep_binary %>%
  group_by(categories, group_prob, rsq_prod, sample_size, loading, n_items, 
           model_type, param.id) %>%
  summarise(n_sim = n(),
            power_est = mean(sig),
            mcse_power = sqrt(power_est * (1 - power_est) / n_sim),
            .groups = "drop")
summary(power_mcse$mcse_power)


# RELATIVE BIAS ----
bias_mcse <- rep_binary %>%
  group_by(categories, group_prob, rsq_prod, sample_size, loading,
           n_items, model_type, param.id) %>%
  summarise(n_sim = n(),
            S2_est = var(est),
            rel_bias = mean(rel.bias),
            mcse_rel_bias = sqrt(S2_est / (true^2 * n_sim)),
            .groups = "drop")
summary(bias_mcse$mcse_rel_bias)


# MEAN SQUARE ERROR ----
MSE_mcse <- rep_binary %>%
  group_by(categories, group_prob, rsq_prod, sample_size, loading, n_items, 
           model_type, param.id) %>%
  summarise(n_sim = n(), 
            S2 = var(squared_bias),  # sample variance of squared errors
            mcse_MSE = sqrt(S2 / n_sim), # MCSE
            .groups = "drop")
summary(MSE_mcse$mcse_MSE)


# CI COVERAGE ----
coverage_mcse <- rep_binary %>%
  group_by(categories, group_prob, rsq_prod, sample_size, loading, n_items, 
           model_type, param.id) %>%
  summarise(n_sim = n(),
            coverage = mean(ci.cov),
            mcse_coverage = sqrt(coverage * (1 - coverage) / n_sim),
            .groups = "drop")

summary(coverage_mcse$mcse_coverage)

