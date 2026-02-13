# Purpose: use ANOVA to parse out most influential effects
options(scipen = 999)
load(file = "rep_binary.rda")

# Set Variables as factors for ANOVAs ----
rep_binary$group_prob <- as.factor(rep_binary$group_prob)
rep_binary$rsq_prod <- as.factor(rep_binary$rsq_prod)
rep_binary$sample_size <- as.factor(rep_binary$sample_size)
rep_binary$loading <- as.factor(rep_binary$loading)
rep_binary$n_items <- as.factor(rep_binary$n_items)
rep_binary$model_type <- as.factor(rep_binary$model_type)


# Function to sort effect sizes ----
sort_effects <- function(data, outcome, predictors, order = 1, 
                         effect_cutoff = 0.01){
  # Build formula
  if (order == 1) {
    formula_str <- paste(outcome, "~", paste(predictors, collapse = " + "))
  } else {
    formula_str <- paste(outcome,"~ (",paste(predictors, collapse = " + "),")^", order)
  }
  
  formula_obj <- as.formula(formula_str)
  
  # Fit ANOVA
  model <- aov(formula_obj, data = data)
  model_summary <- summary(model)[[1]]
  
  # Remove residual row
  model_summary <- model_summary[rownames(model_summary) != "Residuals", ]
  
  # Compute proportion of variance explained
  total_ss <- sum(model_summary[["Sum Sq"]])
  effects_prop <- model_summary[["Sum Sq"]] / total_ss
  
  # Build output table
  results <- data.frame(Effect = rownames(model_summary), 
                        Sum_Sq = model_summary[["Sum Sq"]],
                        Prop_SS = effects_prop)
  
  # Filter + sort
  results <- results[results$Prop_SS > effect_cutoff, ]
  results <- results[order(-results$Prop_SS), ]
  
  rownames(results) <- NULL
  return(results)
}



# STANDARDIZED BIAS ----
# One-way 
sort_effects(rep_binary, "bias_z", 
             c("group_prob","rsq_prod","sample_size","loading","n_items","model_type"), 
             order = 1, effect_cutoff = 0.01)

# Two-way 
sort_effects(rep_binary, "bias_z", 
             c("group_prob","rsq_prod","sample_size","loading","n_items","model_type"), 
             order = 2, effect_cutoff = 0.01)

# Three-way 
sort_effects(rep_binary, "bias_z", 
             c("group_prob","rsq_prod","sample_size","loading","n_items","model_type"), 
             order = 3, effect_cutoff = 0.01)


# PERCENT BIAS ----
# One-way 
sort_effects(rep_binary, "rel.bias", 
             c("group_prob","rsq_prod","sample_size","loading","n_items","model_type"), 
             order = 1, effect_cutoff = 0.01)

# Two-way 
sort_effects(rep_binary, "rel.bias", 
             c("group_prob","rsq_prod","sample_size","loading","n_items","model_type"), 
             order = 2, effect_cutoff = 0.01)

# Three-way 
sort_effects(rep_binary, "rel.bias", 
             c("group_prob","rsq_prod","sample_size","loading","n_items","model_type"), 
             order = 3, effect_cutoff = 0.01)


# CI COVERAGE ----
# One-way 
sort_effects(rep_binary, "ci.cov", 
             c("group_prob","rsq_prod","sample_size","loading","n_items","model_type"), 
             order = 1, effect_cutoff = 0.01)

# Two-way 
sort_effects(rep_binary, "ci.cov", 
             c("group_prob","rsq_prod","sample_size","loading","n_items","model_type"), 
             order = 2, effect_cutoff = 0.01)

# Three-way 
sort_effects(rep_binary, "ci.cov", 
             c("group_prob","rsq_prod","sample_size","loading","n_items","model_type"), 
             order = 3, effect_cutoff = 0.01)


# MSE ----
# One-way 
sort_effects(rep_binary, "squared_bias", 
             c("group_prob","rsq_prod","sample_size","loading","n_items","model_type"), 
             order = 1, effect_cutoff = 0.01)

# Two-way 
sort_effects(rep_binary, "squared_bias", 
             c("group_prob","rsq_prod","sample_size","loading","n_items","model_type"), 
             order = 2, effect_cutoff = 0.01)

# Three-way 
sort_effects(rep_binary, "squared_bias", 
             c("group_prob","rsq_prod","sample_size","loading","n_items","model_type"), 
             order = 3, effect_cutoff = 0.01)

