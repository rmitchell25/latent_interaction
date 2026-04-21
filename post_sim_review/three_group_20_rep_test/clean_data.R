options(scipen = 999)
library(dplyr)

# Load data ----
rep_level <- read.table("~/Desktop/rep_test.dat", 
                        quote="\"", comment.char="")
colnames(rep_level) <- c("categories", "group_prob", "rsq_prod", "sample_size", 
                         "loading", "n_items", "model_type", "param.id", "est", 
                         "se","true", "bias", "rel.bias", "squared_bias","pval", 
                         "sig", "ci.cov", "ci.zero","joint.stat","joint.p",
                         "no_cvg_check", "psr.max","neff.min","ci.low","ci.up",
                         "rep","seed")


# Check number of replications ----
# Ensure the grouping columns exist
group_cols <- c("categories", "group_prob", "rsq_prod", "sample_size",
                "loading", "n_items", "model_type", "param.id")

# Use aggregate to count rows for each combination of grouping levels
agg_counts <- aggregate(
  x = rep(1, nrow(rep_level)),           # a vector to count each row
  by = rep_level[, group_cols, drop = FALSE],
  FUN = sum)

names(agg_counts)[names(agg_counts) == "x"] <- "row_count"
agg_counts

any(agg_counts$row_count > 20)        # Checked 




# Check Convergence Issues ----
conv_diagnostics <- rep_level %>%
  mutate(
    no_cvg_check = as.logical(no_cvg_check),
    inadmiss1 = psr.max > 1.10 & neff.min < 100,
    inadmiss2 = psr.max < 1.10 & neff.min < 100,
    
    total_cut = case_when(
      no_cvg_check ~ TRUE,
      model_type == 1 & inadmiss1 ~ TRUE,
      model_type == 1 & inadmiss2 ~ TRUE,
      TRUE ~ FALSE)) %>%
  group_by(categories, group_prob, rsq_prod, sample_size, loading, n_items, model_type) %>%
  summarise(
    n_rep = n(),
    n_fail = sum(no_cvg_check, na.rm = TRUE),
    n_inadmiss1 = sum(inadmiss1, na.rm = TRUE),
    n_inadmiss2 = sum(inadmiss2, na.rm = TRUE),
    total_cut = sum(total_cut, na.rm = TRUE),
    prop_total_cut = total_cut/n_rep,
    .groups = "drop")

nrow(conv_diagnostics[conv_diagnostics$prop_total_cut > .50,])/nrow(conv_diagnostics)

rep_level <- rep_level %>%
  left_join(conv_diagnostics,
    by = c("categories", "group_prob", "rsq_prod","sample_size", "loading", 
           "n_items", "model_type")) %>%
  filter(prop_total_cut <= 0.50 | is.na(prop_total_cut))




# Clean up CI Coverage variable ----
# rep_binary$CI_cov <- ifelse(rep_binary$true > rep_binary$ci.low &
#                             rep_binary$true < rep_binary$ci.up, 1, 0)
# mean(rep_binary$CI_cov == rep_binary$ci.cov, na.rm = TRUE)
# rep_binary[rep_binary$ci.cov != rep_binary$CI_cov, ]
# 
# 
# rep_binary$CI_zero <- ifelse(0 > rep_binary$ci.low & 0 < rep_binary$ci.up, 1, 0)
# mean(rep_binary$CI_zero == rep_binary$ci.zero, na.rm = TRUE)
# rep_binary[rep_binary$ci.zero != rep_binary$CI_zero, ]





# Flag Bias ----
flag_sim_bias <- function(data, z_threshold, rel_threshold) {
  
  # Variance across replications
  variance_summary <- sapply(data["bias"], var, na.rm = TRUE)
  
  # Standardize bias
  z_bias <- as.data.frame(scale(data["bias"], center = TRUE, scale = TRUE))
  
  # z > 5  indicators 
  z_indicators <- as.integer(abs(z_bias) > z_threshold)
  
  # rel.bias > 500
  rel_indicators <- as.integer(abs(data["rel.bias"]) > 500)
  
  # Bind to data
  data_full <- cbind(data, z_bias, z_indicators, rel_indicators)
  
  # 6. Trimmed dataset
  data_trimmed <- data_full[data_full$z_indicators == 0, ]
  
  # 7. Summary of cuts
  cut_summary <- list(
    total_rows = nrow(data_full),
    rows_cut = sum(data_full$z_indicators, na.rm = T),
    percent_cut = mean(data_full$z_indicators, na.rm = T) * 100,
    rows_removed_indices = which(data_full$z_indicators == 1))
  
  return(list(full_data = data_full, trimmed_data = data_trimmed,
              cut_summary = cut_summary))
}


result <- flag_sim_bias(rep_level, z_threshold = 5, rel_threshold = 500)

rep_level <- result$trimmed_data
cut_info <- result$cut_summary
rep_level_full <- result$full_data


setwd("/Users/remus/Documents/GitHub/latent_interaction/post_sim_review/three_group_20_rep_test")
save(rep_level, file = "rep_level_20rep.rda")


