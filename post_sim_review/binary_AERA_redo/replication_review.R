library(dplyr)

summarise_condition <- function(data, loading_val = 0.5, group_prob_val = 2) {
  
  dat <- data %>%
    filter(
      loading == loading_val,
      group_prob == group_prob_val,
      rsq_prod == 0.03
    ) %>%
    mutate(
      crashed = is.na(est),
      exclude_psr  = psr.max > 1.05,
      exclude_neff = neff.min < 100,
      exclude_any  = exclude_psr | exclude_neff
    )
  
  total_reps <- length(unique(dat$rep))
  cat("Found", total_reps, "replications\n")
  cat("Found", nrow(dat), "rows of results\n")
  
  for (n in sort(unique(dat$sample_size))) {
    
    cat("\n########################################\n")
    cat("N =", n, "\n")
    cat("########################################\n")
    
    d <- dat %>% filter(sample_size == n)
    
    attempted <- nrow(d)
    crashed <- sum(d$crashed)
    successful <- sum(!d$crashed)
    
    cat(sprintf("Crashed reps: %d (%.1f%%)\n", crashed, 100 * crashed / attempted))
    cat(sprintf("Successful reps: %d (%.1f%%)\n", successful, 100 * successful / attempted))
    
    # ---------------------------
    # ALL successful reps
    # ---------------------------
    cat("\n-- Parameter Recovery (all successful reps) --\n")
    
    d_all <- d %>% filter(!crashed)
    
    emp_sd_df <- d_all %>%
      group_by(param.id) %>%
      summarise(emp_sd = sd(est, na.rm = TRUE), .groups = "drop")
    
    summ_all <- d_all %>%
      group_by(param.id) %>%
      summarise(
        true = mean(true, na.rm = TRUE),
        est = mean(est, na.rm = TRUE),
        avg_post_sd = mean(se, na.rm = TRUE),
        pct_sig = mean(sig, na.rm = TRUE),
        .groups = "drop")
    
    summ_all <- summ_all %>%
      left_join(emp_sd_df, by = "param.id") %>%
      mutate(
        bias_sd = (est - true) / emp_sd,
        z = (est - true) / avg_post_sd)
    
    p1 <-summ_all %>%
      mutate(across(where(is.numeric), ~ round(.x, 4))) 
    
    print(as.data.frame(p1))
    
    # ---------------------------
    # FILTERED reps
    # ---------------------------
    cat("\n-- Parameter Recovery (excluding high PSR or low N_eff) --\n")
    
    d_filt <- d_all %>% filter(!exclude_any)
    
    emp_sd_filt <- d_filt %>%
      group_by(param.id) %>%
      summarise(emp_sd = sd(est, na.rm = TRUE), .groups = "drop")
    
    summ_filt <- d_filt %>%
      group_by(param.id) %>%
      summarise(
        true = mean(true, na.rm = TRUE),
        est = mean(est, na.rm = TRUE),
        avg_post_sd = mean(se, na.rm = TRUE),
        pct_sig = mean(sig, na.rm = TRUE),
        .groups = "drop")
    
    summ_filt <- summ_filt %>%
      left_join(emp_sd_filt, by = "param.id") %>%
      mutate(
        bias_sd = (est - true) / emp_sd,
        z = (est - true) / avg_post_sd)
    
    
    p2<- summ_filt %>%
            mutate(across(where(is.numeric), ~ round(.x, 4)))
    
    print(as.data.frame(p2))
    
    # ---------------------------
    # Convergence summary
    # ---------------------------
    cat("\n-- Convergence / N_eff Summary (all successful reps) --\n")
    
    conv <- d_all %>%
      summarise(
        psr_min = min(psr.max, na.rm = TRUE),
        psr_q1 = quantile(psr.max, .25, na.rm = TRUE),
        psr_med = median(psr.max, na.rm = TRUE),
        psr_mean = mean(psr.max, na.rm = TRUE),
        psr_q3 = quantile(psr.max, .75, na.rm = TRUE),
        psr_max = max(psr.max, na.rm = TRUE),
        
        neff_min = min(neff.min, na.rm = TRUE),
        neff_q1 = quantile(neff.min, .25, na.rm = TRUE),
        neff_med = median(neff.min, na.rm = TRUE),
        neff_mean = mean(neff.min, na.rm = TRUE),
        neff_q3 = quantile(neff.min, .75, na.rm = TRUE),
        neff_max = max(neff.min, na.rm = TRUE)
      )
    
    print(round(conv, 3))
  }
}

summarise_condition(rep_binary, loading_val = 0.5, group_prob_val = 1)
summarise_condition(rep_binary, loading_val = 0.8, group_prob_val = 1)
summarise_condition(rep_binary, loading_val = 0.5, group_prob_val = 2)
summarise_condition(rep_binary, loading_val = 0.8, group_prob_val = 2)
