options(scipen = 999)
library(dplyr)

# Load data ----
colnames(rep_level) <- c("categories", "group_prob", "rsq_prod", "sample_size", 
                         "loading", "n_items", "model_type", "param.id", "est", 
                         "se","true", "bias", "rel.bias", "squared_bias","pval", 
                         "sig", "ci.cov", "ci.zero","joint.stat","joint.p",
                         "no_cvg_check", "psr.max","neff.min","ci.low","ci.up")




# Check Convergence Issues ----
conv_summary_rep <- rep_level %>%
  group_by(categories, group_prob, rsq_prod, sample_size, loading, n_items, model_type) %>%
  summarise(no_cvg = max(no_cvg_check, na.rm = TRUE), .groups = "drop")

# Filter out models that did not converge at all 
rep_clean <- rep_level[rep_level$no_cvg_check == 0, ]

# Check for inadmissable convergence
inadmiss_summary_rep <- rep_clean %>%
  filter(!is.na(psr.max), !is.na(neff.min)) %>%
  group_by(categories, group_prob, rsq_prod, sample_size, loading, n_items, model_type) %>%
  reframe(inadmiss_check = (psr.max > 1.05 & neff.min < 100), .groups = "drop")
summary(inadmiss_summary_rep$inadmiss_check)

# Filter out Blimp models that idadmissably converged
rep_clean <- rep_clean[rep_clean$model_type != 1 |(rep_clean$model_type == 1 &
                                                     rep_clean$neff.min >= 100 &
                                                     rep_clean$psr.max <= 1.10),]



# Filter out three-category cases and non-interaction params ----
rep_binary <- rep_clean[rep_clean$categories == 2, ]
rep_binary <- rep_binary[rep_binary$param.id == "beta_XG2",]


# Clean up CI Coverage variable ----
rep_binary$test <- ifelse(rep_binary$true > rep_binary$ci.low &
                            rep_binary$true < rep_binary$ci.up, 1, 0)
mean(rep_binary$test == rep_binary$ci.cov, na.rm = TRUE)



# Standardize Performance Measures (z > 5 get cut) ----
standardize_check_sim <- function(data, vars, cut_threshold = 5) {
  
  # 1. Variance across replications
  variance_summary <- sapply(data[vars], var, na.rm = TRUE)
  
  # 2. Standardize variables
  z_scores <- as.data.frame(scale(data[vars], center = TRUE, scale = TRUE))
  names(z_scores) <- paste0(vars, "_z")
  
  # 3. Cut indicators
  cut_indicators <- as.data.frame(abs(z_scores) > cut_threshold)
  names(cut_indicators) <- paste0(vars, "_cut")
  cut_indicators[] <- lapply(cut_indicators, as.integer)
  
  # 4. Bind back to data
  data_full <- cbind(data, z_scores, cut_indicators)
  
  # 5. Overall cut flag
  data_full$any_cut <- as.integer(rowSums(cut_indicators, na.rm = TRUE) > 0)
  
  # 6. Trimmed dataset
  data_trimmed <- data_full[data_full$any_cut == 0, ]
  
  # 7. Summary of cuts
  cut_summary <- list(
    total_rows = nrow(data_full),
    rows_cut = sum(data_full$any_cut),
    percent_cut = mean(data_full$any_cut) * 100,
    cuts_by_variable = colSums(cut_indicators, na.rm = TRUE),
    rows_removed_indices = which(data_full$any_cut == 1))
  
  return(list(full_data = data_full, trimmed_data = data_trimmed,
              variance = variance_summary, cut_summary = cut_summary))
}

vars_to_check <- c("sig", "ci.cov", "ci.zero", "bias")

result <- standardize_check_sim(rep_binary, vars_to_check, cut_threshold = 5)

rep_binary <- result$trimmed_data
variance_summary   <- result$variance
cut_info <- result$cut_summary


save(rep_binary, file = "rep_binary.rda")


# AGGREGATE

