library(dplyr)
library(tidyr)
options(scipen = 999)

# Load data ----
setwd("/Users/remus/Documents/GitHub/latent_interaction/post_sim_review/binary_AERA")
load(file = "agg_results.rda")

# thresholds for meaningful differences
###     relative bias - 5%
###     power/TIE -  1% 
###     MSE ratios - .05
###     CI coverage - based on MCSE


# CI Coverage ----
# Incorporate MCSEs for coverage into aggregated dataset
coverage_mcse <- agg_results %>%
  group_by(categories, group_prob, rsq_prod, sample_size, loading, n_items, 
           model_type, param.id) %>%
  summarise(n_sim = n(), coverage = mean(CI_cov),
    mcse_coverage = sqrt(coverage * (1 - coverage) / n_sim), .groups = "drop")
agg_results <- agg_results %>%
  left_join(coverage_mcse, by = c("categories","group_prob","rsq_prod","sample_size",
                                  "loading","n_items","model_type","param.id"))
hist(agg_results$mcse_coverage)


# get coverage and coverage MCSEs in wide format for comparison
coverage_wide <- agg_results %>%
  select(categories, group_prob, rsq_prod, sample_size, loading, n_items, 
         param.id, model_type, coverage, mcse_coverage) %>%
  pivot_wider(names_from = model_type, 
              values_from = c(coverage, mcse_coverage))

# save unique model types and empty list for storage
model_types <- unique(agg_results$model_type)
cov_results_list <- list()

# loop through all comparisons for specific set of conditions
for (i in 1:length(model_types)) {
  for (j in (i+1):length(model_types)) {
    if (j > length(model_types)) next
    
    model1 <- model_types[i]
    model2 <- model_types[j]
    
    cov1_col <- paste("coverage", model1, sep = "_")
    cov2_col <- paste("coverage", model2, sep = "_")
    mcse1_col <- paste("mcse_coverage", model1, sep = "_")
    mcse2_col <- paste("mcse_coverage", model2, sep = "_")
    
    temp_df <- coverage_wide %>%
      mutate(
        diff = get(cov1_col) - get(cov2_col),
        diff_se = sqrt(get(mcse1_col)^2 + get(mcse2_col)^2),
        z_score = diff / diff_se,
        p_value = 2 * pnorm(-abs(z_score)),
        comparison = paste(model1, "vs", model2)
      ) %>%
      select(categories, group_prob, rsq_prod, sample_size, loading, n_items,
             param.id, diff, diff_se, z_score, p_value, comparison)
    
    cov_results_list[[paste(model1, model2, sep = "_")]] <- temp_df
  }
}

pairwise_results <- bind_rows(cov_results_list)
significant_diffs <- pairwise_results %>%
  filter(p_value < 0.10)
print(significant_diffs)


# Relative Bias ----
# get relative bias in wide format for comparison
relbias_wide <- agg_results %>%
  select(categories, group_prob, rsq_prod, sample_size, loading, n_items, 
         param.id, model_type, rel.bias) %>%
  pivot_wider(names_from = model_type, values_from = c(rel.bias))
colnames(relbias_wide) <- c("categories", "group_prob", "rsq_prod", "sample_size", 
                            "loading", "n_items", "param.id", "relbias_1",
                            "relbias_2", "relbias_3")

# save empty list for storage
rel_results_list <- list()

# loop through all comparisons for specific set of conditions
for (i in 1:length(model_types)) {
  for (j in (i+1):length(model_types)) {
    if (j > length(model_types)) next
    
    model1 <- model_types[i]
    model2 <- model_types[j]
    
    # Get column names for relative bias
    relbias1_col <- paste("relbias", model1, sep = "_")
    relbias2_col <- paste("relbias", model2, sep = "_")
    
    # Just flag >10% differences
    temp_df <- relbias_wide %>%
      mutate(
        diff = abs(get(relbias1_col) - get(relbias2_col)),
        diff_gt_5pct = diff > 0.5,  # Direct check if absolute difference > 0.1 (10%)
        comparison = paste(model1, "vs", model2)
      ) %>%
      select(categories, group_prob, rsq_prod, sample_size, loading, n_items,
             param.id, diff, diff_gt_5pct, comparison)
    
    
    rel_results_list[[paste(model1, model2, sep = "_")]] <- temp_df
  }
}

# Combine all results into a single dataframe
rel_results_df <- bind_rows(rel_results_list)

# Count how many comparisons have difference > 5%
summary_stats_rel <- rel_results_df %>%
  group_by(comparison) %>%
  summarise(
    total_comparisons = n(),
    gt_5pct_count = sum(diff_gt_5pct, na.rm = TRUE),
    gt_5pct_percent = mean(diff_gt_5pct, na.rm = TRUE) * 100,
    .groups = "drop")

# Print summary
print(summary_stats_rel)

# Filter for cases where difference > 5%
significant_diffs_rel <- rel_results_df %>%
  filter(diff_gt_5pct == TRUE); significant_diffs_rel





# Power/TIE ----
# get rejection rate in wide format for comparison
reject_wide <- agg_results %>%
  select(categories, group_prob, rsq_prod, sample_size, loading, n_items, 
         param.id, model_type, sig) %>%
  pivot_wider(names_from = model_type, values_from = c(sig))
colnames(reject_wide) <- c("categories", "group_prob", "rsq_prod", "sample_size", 
                            "loading", "n_items", "param.id", "reject_1",
                            "reject_2", "reject_3")

# save empty list for storage
reject_results_list <- list()

# loop through all comparisons for specific set of conditions
for (i in 1:length(model_types)) {
  for (j in (i+1):length(model_types)) {
    if (j > length(model_types)) next
    
    model1 <- model_types[i]
    model2 <- model_types[j]
    
    # Get column names for relative bias
    reject1_col <- paste("reject", model1, sep = "_")
    reject2_col <- paste("reject", model2, sep = "_")
    
    # Just flag >1% differences
    temp_df <- reject_wide %>%
      mutate(
        diff = abs(get(reject1_col) - get(reject2_col)),
        diff_gt_1pct = diff > 0.01,  # Direct check if difference > 0.01
        comparison = paste(model1, "vs", model2)
      ) %>%
      select(categories, group_prob, rsq_prod, sample_size, loading, n_items,
             param.id, diff, diff_gt_1pct, comparison)
    
    
    reject_results_list[[paste(model1, model2, sep = "_")]] <- temp_df
  }
}


# Combine all results into a single dataframe
reject_results_df <- bind_rows(reject_results_list)

# Count how many comparisons have difference > 10%
summary_stats_rej <- reject_results_df %>%
  group_by(comparison) %>%
  summarise(
    total_comparisons = n(),
    gt_1pct_count = sum(diff_gt_1pct, na.rm = TRUE),
    gt_1pct_percent = mean(diff_gt_1pct, na.rm = TRUE) * 100,
    .groups = "drop")

# Print summary
print(summary_stats_rej)

# Filter for cases where difference > 10%
significant_diffs_rej <- reject_results_df %>%
  filter(diff_gt_1pct == TRUE); significant_diffs_rej





# MSE ----
# Calculate MSE ratios

agg_results <- agg_results %>%
  group_by(categories, group_prob, rsq_prod, sample_size, loading, n_items) %>%
  mutate(
    reference_value = if (any(model_type == 2)) squared_bias[model_type == 2] else NA_real_,
    mse_ratio = squared_bias / reference_value) %>%
  ungroup()



# save empty list for storage
mse_results_list <- list()

# loop through all comparisons for specific set of conditions
for (i in 1:nrow(agg_results)) {
  if (!is.na(agg_results$mse_ratio[i])) {
    if (agg_results$mse_ratio[i] < 0.95 | agg_results$mse_ratio[i] > 1.05) {
      mse_results_list[[length(mse_results_list) + 1]] <- data.frame(
        row        = i,
        model_type  = agg_results$model_type[i],
        categories  = agg_results$categories[i],
        group_prob  = agg_results$group_prob[i],
        rsq_prod    = agg_results$rsq_prod[i],
        sample_size = agg_results$sample_size[i],
        loading     = agg_results$loading[i],
        n_items     = agg_results$n_items[i],
        mse_ratio   = round(agg_results$mse_ratio[i], 4)
      )
    }
  }
}

# Combine all list elements into a single dataframe
mse_results_df <- do.call(rbind, mse_results_list)
print(mse_results_df)



