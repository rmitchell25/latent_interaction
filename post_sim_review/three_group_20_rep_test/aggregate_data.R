
load(file = "rep_level_20rep.rda")


# Aggregate Entire Data set ----
colnames(rep_level) <- c("categories","group_prob" ,"rsq_prod","sample_size",
                          "loading","n_items" ,"model_type","param.id","est","se",
                          "true" ,"bias" ,"rel.bias","squared_bias","pval","sig",
                          "ci.cov","ci.zero", "joint.stat" ,"joint.p","no_cvg_check",
                          "psr.max" ,"neff.min","ci.low", "ci.up","rep","seed", "n_rep", 
                          "n_fail", "n_inadmiss1", "n_inadmiss2", "total_cut", 
                          "prop_total_cut","bias_z","z_indicators",
                          "rel_indicators")

agg_results <- aggregate(
  cbind(est, se, true, bias, rel.bias, squared_bias, pval, sig, ci.low, ci.up, 
        ci.cov, ci.zero, bias_z) ~ categories + group_prob + rsq_prod + 
    sample_size + loading + n_items + model_type + param.id,
  data = rep_level, FUN = function(x) mean(x, na.rm = TRUE))

save(agg_results, file = "agg_results.rda")







