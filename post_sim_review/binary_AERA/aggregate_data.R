
load(file = "rep_binary.rda")

# Save results to plug into AI ----
test <- rep_binary[,c("categories","group_prob","rsq_prod","sample_size","loading",
                      "n_items","model_type","param.id","est","se","true","bias",
                      "rel.bias","squared_bias","pval","sig", "ci.low","ci.up",
                      "CI_cov","CI_zero","bias_z")]

agg_results <- aggregate(
  cbind(est, se, bias, rel.bias, squared_bias, sig, CI_cov, CI_zero, bias_z) ~
    categories + group_prob + rsq_prod + sample_size + loading + n_items + 
    model_type + param.id,
  data = test, FUN = function(x) mean(x, na.rm = TRUE))

# write.csv(agg_results, file = "results_AI_test.csv")


# Aggregate Entire Data set ----

agg_results <- aggregate(
  cbind(est, se, true, bias, rel.bias, squared_bias, pval, sig, ci.low, ci.up, 
        CI_cov, CI_zero, bias_z) ~ categories + group_prob + rsq_prod + 
    sample_size + loading + n_items + model_type + param.id,
  data = rep_binary, FUN = function(x) mean(x, na.rm = TRUE))

save(agg_results, file = "agg_results.rda")







