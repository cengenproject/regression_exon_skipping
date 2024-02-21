# Process results from running graph_power4 on cluster
# note this is an interactive draft


# From cluster ----
library(tidyverse) |> suppressPackageStartupMessages()





tib_quic <- read_csv("data/intermediates/231109_permutations_cv/231107_quic_50perm_nosep.csv")
tib_quic <- read_csv("data/intermediates/231109_permutations_cv/240116_6penlti_50perm.csv")
tib_quic <- read_csv("data/intermediates/231109_permutations_cv/240118_quic_7penlt_noperm.csv")
tib_quic <- read_csv("data/intermediates/231109_permutations_cv/240126_params_noperm_7penalties.csv")


# use permutations
perm_pval <- function(statistic, permutation){
  # count number of times the null is as big as the observed statistic
  mean(abs(statistic[ as.logical(permutation) ]) >= abs(statistic[ !permutation ]))
}

# precompute p-values
pvals_quic <- tib_quic |>
  select(penalty, fold, permutation, Rsquared, sum_abs_residuals,
         mean_FEV, loss_frobenius, loss_quadratic,
         prop_non_zero_coefs_litt, prop_non_zero_coefs_nonlitt) |>
  summarize(across(-c(permutation),
                   list(pval = ~perm_pval(.x, permutation))),
            .by = c(penalty, fold)) |>
  summarize(across(-c(fold), partial(mean, na.rm = TRUE)), .by = penalty )



summary_metrics <- tib_quic |>
  filter(permutation == 0) |> select(-permutation) |>
  summarize(across(-c(fold),
                   list(mean = partial(mean, na.rm = TRUE),
                        sd = partial(sd, na.rm = TRUE))),
            .by = penalty ) |>
  left_join(pvals_quic, by = "penalty") |>
  pivot_longer(-penalty,
               names_to = c("metric", "type"),
               names_pattern = "(.+)_(mean|sd|pval)$",
               values_to = "value") |>
  pivot_wider(names_from = "type",
              values_from = "value") |>
  mutate(
    metric = case_when(
      metric == "prop_non_zero_coefs_litt" ~ "literature_TPR",
      metric == "prop_non_zero_coefs_nonlitt" ~ "literature_FPR",
      .default = metric) |>
      fct_inorder()
  )

clipr::write_clip(summary_metrics)
summary_metrics |>
  pivot_wider(id_cols = metric, names_from = penalty, values_from = mean) |>
  clipr::write_clip()

summary_metrics |>
  ggplot(aes(x = penalty, y = mean, ymin = mean - sd, ymax = mean + sd)) +
  theme_classic() +
  facet_wrap(~metric, scales = "free_y") +
  geom_line() +
  geom_errorbar(width = .1) +
  geom_point() +
  # geom_point(aes(color = pval < .05, shape = pval < .05), size = 2) +
  scale_x_log10()




# Compare version ----
tib_quic <- read_csv("data/graph_power4/outputs/240118_tests.csv")
tib_quic <- read_csv("data/graph_power4/outputs/240118_tests_object.csv")
tib_quic <- read_csv("data/graph_power4/outputs/240118_tests_man_npn.csv")
tib_quic <- read_csv("data/graph_power4/outputs/240126_tests_parmaterized_npn.csv")
tib_quic <- read_csv("data/graph_power4/outputs/240126_tests_use_parameters_npn.csv")
tib_quic <- read_csv("data/graph_power4/outputs/240126_tests_params_nosepPSI.csv")
tib_quic <- read_csv("data/graph_power4/outputs/240130_tests_revert_psi.csv")

tib_quic <- read_csv("data/graph_power4/outputs/240131_revertpsi_nosep_noperm_7penalties.csv")
tib_quic <- read_csv("data/graph_power4/outputs/240202_recompandrevertpsi_noperm_7penalties.csv")

# 11 penalties (on cluster), PSI vs counts
tib_quic <- read_csv("data/graph_power4/outputs/240208_recompandrevertpsi_noperm_11penalties.csv")
tib_quic <- read_csv("data/graph_power4/outputs/240208_revertpsi_nosep_noperm_11penalties.csv")


# No imputation on validation PSI
tib_quic <- read_csv("data/graph_power4/outputs/240212_npntrunc_imput_psi_noperm_7penalties.csv")
tib_quic <- read_csv("data/graph_power4/outputs/240220_npnshrink_imptrain_noperm_7penalties.csv")

# added power law
tib_quic <- read_csv("data/graph_power4/outputs/240220_npnshrink_imptrain_noperm_7penalties_powerlaw.csv")



summary_metrics <- tib_quic |>
  filter(permutation == 0) |> select(-permutation) |>
  summarize(across(-c(fold),
                   list(mean = partial(mean, na.rm = TRUE),
                        sd = partial(sd, na.rm = TRUE))),
            .by = penalty ) |>
  pivot_longer(-penalty,
               names_to = c("metric", "type"),
               names_pattern = "(.+)_(mean|sd|pval)$",
               values_to = "value") |>
  pivot_wider(names_from = "type",
              values_from = "value") |>
  mutate(
    metric = case_when(
      metric == "prop_non_zero_coefs_litt" ~ "literature_TPR",
      metric == "prop_non_zero_coefs_nonlitt" ~ "literature_FPR",
      .default = metric) |>
      fct_inorder()
  )

summary_metrics |>
  ggplot(aes(x = penalty, y = mean, ymin = mean - sd, ymax = mean + sd)) +
  theme_classic() +
  facet_wrap(~metric, scales = "free_y") +
  geom_line() +
  geom_errorbar(width = .1) +
  geom_point() +
  scale_x_log10()


summary_metrics |>
  filter(grepl("loss", metric)) |>
  ggplot(aes(x = penalty, y = mean, ymin = mean - sd, ymax = mean + sd)) +
  theme_classic() +
  facet_wrap(~metric, scales = "free_y") +
  geom_line() +
  geom_errorbar(width = .1) +
  geom_point() +
  scale_x_log10()




tib_quic1 <- read_csv("data/graph_power4/outputs/240212_npnshrink_imp_noperm_7penalties.csv")
tib_quic2 <- read_csv("data/graph_power4/outputs/240220_npnshrink_imptrain_noperm_7penalties.csv")

tib_quic <- bind_rows(
  tib_quic1 |>
    add_column(run = "both imputed"),
  tib_quic2 |>
    add_column(run = "training imputed")
)


summary_metrics <- tib_quic |>
  filter(permutation == 0) |> select(-permutation) |>
  summarize(across(-c(fold),
                   list(mean = partial(mean, na.rm = TRUE),
                        sd = partial(sd, na.rm = TRUE))),
            .by = c(penalty, run) ) |>
  pivot_longer(-c(penalty, run),
               names_to = c("metric", "type"),
               names_pattern = "(.+)_(mean|sd|pval)$",
               values_to = "value") |>
  pivot_wider(names_from = "type",
              values_from = "value") |>
  mutate(
    metric = case_when(
      metric == "prop_non_zero_coefs_litt" ~ "literature_TPR",
      metric == "prop_non_zero_coefs_nonlitt" ~ "literature_FPR",
      .default = metric) |>
      fct_inorder()
  )


summary_metrics |>
  ggplot(aes(x = penalty, y = mean, ymin = mean - sd, ymax = mean + sd, color = run)) +
  theme_classic() +
  facet_wrap(~metric, scales = "free_y") +
  geom_line() +
  geom_errorbar(width = .1) +
  geom_point() +
  scale_x_log10()



# Compare sep and PSI ----

# 11 penalties (on cluster), PSI vs counts
tib_quic_psi <- read_csv("data/graph_power4/outputs/240208_recompandrevertpsi_noperm_11penalties.csv")
tib_quic_counts <- read_csv("data/graph_power4/outputs/240208_revertpsi_nosep_noperm_11penalties.csv")

tib_quic <- bind_rows(tib_quic_psi |>
                        add_column(run = "psi"),
                      tib_quic_counts |>
                        add_column(run = "counts"))

summary_metrics <- tib_quic |>
  filter(permutation == 0) |> select(-permutation) |>
  summarize(across(-c(fold),
                   list(mean = partial(mean, na.rm = TRUE),
                        sd = partial(sd, na.rm = TRUE))),
            .by = c(penalty, run) ) |>
  pivot_longer(-c(penalty, run),
               names_to = c("metric", "type"),
               names_pattern = "(.+)_(mean|sd|pval)$",
               values_to = "value") |>
  pivot_wider(names_from = "type",
              values_from = "value") |>
  mutate(
    metric = case_when(
      metric == "prop_non_zero_coefs_litt" ~ "literature_TPR",
      metric == "prop_non_zero_coefs_nonlitt" ~ "literature_FPR",
      .default = metric) |>
      fct_inorder()
  )

summary_metrics |>
  ggplot(aes(x = penalty, y = mean, ymin = mean - sd, ymax = mean + sd, color = run)) +
  theme_classic() +
  facet_wrap(~metric, scales = "free_y") +
  geom_line() +
  geom_errorbar(width = .1) +
  geom_point() +
  scale_x_log10()




# Compare normalizations ----



# 7 penalties (on PC), NPN shrunken vs Zscore
tib_quic_npn <- read_csv("data/graph_power4/outputs/240202_recompandrevertpsi_noperm_7penalties.csv")
tib_quic_zscore <- read_csv("data/graph_power4/outputs/240208_zscore_counts_noperm_7penalties.csv")

# NPN shrinage vs truncation
tib_quic_shrink <- read_csv("data/graph_power4/outputs/240208_psi_noperm_7penalties.csv")
tib_quic_trunc <- read_csv("data/graph_power4/outputs/240209_npntrunc_psi_noperm_7penalties.csv")

# with imputation: shrinkage vs truncation
tib_quic_shrink <- read_csv("data/graph_power4/outputs/240212_npnshrink_imp_noperm_7penalties.csv")
tib_quic_trunc <- read_csv("data/graph_power4/outputs/240212_npntrunc_imput_psi_noperm_7penalties.csv")
tib_quic_zscore <- read_csv("data/graph_power4/outputs/240212_zscore_imput_psi_noperm_7penalties.csv")


tib_quic <- bind_rows(tib_quic_shrink |>
                        add_column(run = "NPN (shrunken)"),
                      tib_quic_trunc |>
                        add_column(run = "NPN (truncation)"),
                      tib_quic_zscore |>
                        add_column(run = "Zscore"))

summary_metrics <- tib_quic |>
  filter(permutation == 0) |> select(-permutation) |>
  summarize(across(-c(fold),
                   list(mean = partial(mean, na.rm = TRUE),
                        sd = partial(sd, na.rm = TRUE))),
            .by = c(penalty, run) ) |>
  pivot_longer(-c(penalty, run),
               names_to = c("metric", "type"),
               names_pattern = "(.+)_(mean|sd|pval)$",
               values_to = "value") |>
  pivot_wider(names_from = "type",
              values_from = "value") |>
  mutate(
    metric = case_when(
      metric == "prop_non_zero_coefs_litt" ~ "literature_TPR",
      metric == "prop_non_zero_coefs_nonlitt" ~ "literature_FPR",
      .default = metric) |>
      fct_inorder()
  )

summary_metrics |>
  ggplot(aes(x = penalty, y = mean, ymin = mean - sd, ymax = mean + sd, color = run)) +
  theme_classic() +
  facet_wrap(~metric, scales = "free_y") +
  geom_line() +
  geom_errorbar(width = .1) +
  geom_point() +
  scale_x_log10()




# Compare imputations ----


# with imputation: shrinkage vs truncation
tib_quic_knn <- read_csv("data/graph_power4/outputs/240220_npnshrink_imptrain_noperm_7penalties_powerlaw.csv")
tib_quic_median <- read_csv("data/graph_power4/outputs/240220_npnshrink_impmedian_noperm_7penalties_powerlaw.csv")



tib_quic <- bind_rows(tib_quic_knn |>
                        add_column(run = "KNN"),
                      tib_quic_median |>
                        add_column(run = "median"))

summary_metrics <- tib_quic |>
  filter(permutation == 0) |> select(-permutation) |>
  summarize(across(-c(fold),
                   list(mean = partial(mean, na.rm = TRUE),
                        sd = partial(sd, na.rm = TRUE))),
            .by = c(penalty, run) ) |>
  pivot_longer(-c(penalty, run),
               names_to = c("metric", "type"),
               names_pattern = "(.+)_(mean|sd|pval)$",
               values_to = "value") |>
  pivot_wider(names_from = "type",
              values_from = "value") |>
  mutate(
    metric = case_when(
      metric == "prop_non_zero_coefs_litt" ~ "literature_TPR",
      metric == "prop_non_zero_coefs_nonlitt" ~ "literature_FPR",
      .default = metric) |>
      fct_inorder()
  )

summary_metrics |>
  ggplot(aes(x = penalty, y = mean, ymin = mean - sd, ymax = mean + sd, color = run)) +
  theme_classic() +
  facet_wrap(~metric, scales = "free_y") +
  geom_line() +
  geom_errorbar(width = .1) +
  geom_point() +
  scale_x_log10()



# ROC curve
summary_metrics |>
  filter(startsWith(as.character(metric),"literature")) |>
  mutate(metric = str_split_i(metric, "_", 2)) |>
  pivot_wider(names_from = "metric",
              values_from = c("mean", "sd")) |>
  ggplot(aes(x = mean_FPR, y = mean_TPR,
             ymin = mean_TPR - sd_TPR,
             ymax = mean_TPR + sd_TPR,
             xmin = mean_FPR - sd_FPR,
             xmax = mean_FPR + sd_FPR,
             color = run)) +
  theme_classic() +
  # facet_wrap(~loss, scales = "free_y") +
  geom_line() +
  geom_errorbar(width = .01) +
  geom_errorbarh(height = .005) +
  geom_point() +
  geom_abline(linetype = "dashed", color = "grey") +
  scale_x_continuous(limits = c(0,1)) +
  scale_y_continuous(limits = c(0,1))





# Bias loss ----
tib_quic <- read_csv("data/graph_power4/outputs/240220_npnshrink_impmedian_noperm_7penalties_powerlaw.csv")



summary_metrics <- tib_quic |>
  filter(permutation == 0) |> select(-permutation) |>
  summarize(across(-c(fold),
                   list(mean = partial(mean, na.rm = TRUE),
                        sd = partial(sd, na.rm = TRUE))),
            .by = penalty ) |>
  pivot_longer(-penalty,
               names_to = c("metric", "type"),
               names_pattern = "(.+)_(mean|sd|pval)$",
               values_to = "value") |>
  pivot_wider(names_from = "type",
              values_from = "value") |>
  mutate(
    metric = case_when(
      metric == "prop_non_zero_coefs_litt" ~ "literature_TPR",
      metric == "prop_non_zero_coefs_nonlitt" ~ "literature_FPR",
      .default = metric) |>
      fct_inorder()
  )

summary_metrics |>
  ggplot(aes(x = penalty, y = mean, ymin = mean - sd, ymax = mean + sd)) +
  theme_classic() +
  facet_wrap(~metric, scales = "free_y") +
  geom_line() +
  geom_errorbar(width = .1) +
  geom_point() +
  scale_x_log10()


# Bias loss
summary_metrics |>
  filter(grepl("loss", metric)) |>
  mutate(type = if_else(startsWith(as.character(metric), "bias_"),
                        "Training bias",
                        "Validation sampling"),
         loss = str_match(metric, "loss_(.*)$")[,2]) |>
  ggplot(aes(x = penalty, y = mean, ymin = mean - sd, ymax = mean + sd, color = type)) +
  theme_classic() +
  facet_wrap(~loss, scales = "free_y") +
  geom_line() +
  geom_errorbar(width = .1) +
  geom_point() +
  scale_x_log10()


# ROC curve
summary_metrics |>
  filter(startsWith(as.character(metric),"literature")) |>
  mutate(metric = str_split_i(metric, "_", 2)) |>
  pivot_wider(names_from = "metric",
              values_from = c("mean", "sd")) |>
  ggplot(aes(x = mean_FPR, y = mean_TPR,
             ymin = mean_TPR - sd_TPR,
             ymax = mean_TPR + sd_TPR,
             xmin = mean_FPR - sd_FPR,
             xmax = mean_FPR + sd_FPR)) +
  theme_classic() +
  # facet_wrap(~loss, scales = "free_y") +
  geom_line() +
  geom_errorbar(width = .01) +
  geom_errorbarh(height = .005) +
  geom_point() +
  geom_abline(linetype = "dashed", color = "grey") +
  scale_x_continuous(limits = c(0,1)) +
  scale_y_continuous(limits = c(0,1))

