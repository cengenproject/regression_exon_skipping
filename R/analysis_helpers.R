# helper functions for the analysis


read_one_res <- function(run_name, dir){
  filename <- paste0(run_name, ".csv")
  read_csv(file.path(dir, filename),
           col_types = cols(
             penalty = col_double(),
             fold = col_double(),
             permutation = col_double(),
             Rsquared = col_double(),
             mean_FEV = col_double(),
             loss_frobenius = col_double(),
             loss_quadratic = col_double(),
             bias_loss_frobenius = col_double(),
             bias_loss_quadratic = col_double(),
             literature_TPR = col_double(),
             literature_FPR = col_double(),
             power_law = col_double(),
             sparsity = col_double()
           ))
}






summarize_metrics_by_pen <- function(tib){
  tib |>
    filter(permutation == 0) |> select(-permutation) |>
    mutate(`TPR/FPR` = literature_TPR/literature_FPR) |>
    summarize(across(-c(fold),
                     list(mean = partial(mean, na.rm = TRUE),
                          sd = partial(sd, na.rm = TRUE))),
              .by = penalty ) |>
    pivot_longer(-penalty,
                 names_to = c("metric", "type"),
                 names_pattern = "(.+)_(mean|sd|pval)$",
                 values_to = "value") |>
    pivot_wider(names_from = "type",
                values_from = "value")
}

# returns wide format, need to fuse with above
summarize_metrics_by_pen2 <- function(tib){
  tib |>
    filter(permutation == 0) |> select(-permutation) |>
    mutate(`TPR/FPR` = literature_TPR/literature_FPR) |>
    summarize(across(-c(fold),
                     list(mean = partial(mean, na.rm = TRUE),
                          sd = partial(sd, na.rm = TRUE))),
              .by = penalty )
}


summarize_metrics_by_spars <- function(tib){
  tib |>
    filter(permutation == 0) |> select(-permutation) |>
    mutate(`TPR/FPR` = literature_TPR/literature_FPR) |>
    summarize(across(-c(fold),
                     list(mean = partial(mean, na.rm = TRUE),
                          sd = partial(sd, na.rm = TRUE))),
              .by = penalty ) |>
    select(- sparsity_sd) |>
    pivot_longer(-c(penalty, sparsity_mean),
                 names_to = c("metric", "type"),
                 names_pattern = "(.+)_(mean|sd|pval)$",
                 values_to = "value") |>
    pivot_wider(names_from = "type",
                values_from = "value")
}




permutation_pval <- function(statistic, permutation){
  # count prop of times the null is as big as the observed statistic
  
  permuted <- abs(statistic[ as.logical(permutation) ])
  observed <- abs(statistic[ !permutation ])
  
  mean(permuted >= observed)
}



get_pvals <- function(tib){
  tib |>
    summarize(across(-c(permutation),
                     list(pval = ~ permutation_pval(.x, permutation))),
              .by = c(penalty, fold)) |>
    summarize(across(-c(fold), partial(mean, na.rm = TRUE)), .by = penalty )
}




