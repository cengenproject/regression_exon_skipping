# Process results from running graph_power4 on cluster
# note this is an interactive draft


# From cluster ----
library(tidyverse) |> suppressPackageStartupMessages()

design <- "
 CD#
 AB#
 EFG
"

export_dir <- "presentations/240308_figures"



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
tib_quic <- read_csv("data/graph_power4/outputs/240314_npnshrink_median_21_5_QUIC.csv")
tib_quic <- read_csv("data/graph_power4/outputs/240314b_npnshrink_median_21_5_QUIC.csv")
tib_quic <- read_csv("data/graph_power4/outputs/240314c_npnshrink_median_21_5_QUIC.csv")
tib_quic <- read_csv("data/graph_power4/outputs/240314d_QUIC_PSI_npnshrink_median_21_5.csv")
tib_quic <- read_csv("data/graph_power4/outputs/240314e_QUIC_counts_npnshrink_median_21_5.csv")
tib_quic <- read_csv("data/graph_power4/outputs/240314f_QUIC_counts_npnshrink_median_2_6.csv")
tib_quic <- read_csv("data/graph_power4/outputs/240314g_glasso_PSI_npnshrink_median_2_4.csv")
tib_quic <- read_csv("data/graph_power4/outputs/240314g_CLIME_PSI_npnshrink_median_1_4.csv")
tib_quic <- read_csv("data/graph_power4/outputs/240314g_SCIO_PSI_npnshrink_median_1_4.csv")


# plot as a function of sparsity


summary_metrics <- tib_quic |>
  filter(permutation == 0) |> select(-permutation) |>
  mutate(`TPR/FPR` = literature_TPR/literature_FPR) |>
  summarize(across(-c(fold),
                   list(mean = partial(mean, na.rm = TRUE),
                        sd = partial(sd, na.rm = TRUE))),
            .by = c(penalty, sparsity) ) |>
  pivot_longer(-c(penalty, sparsity),
               names_to = c("metric", "type"),
               names_pattern = "(.+)_(mean|sd|pval)$",
               values_to = "value") |>
  pivot_wider(names_from = "type",
              values_from = "value")



design <- "
 EF
 AB
 IG
 DC
 KH
"

summary_metrics |>
  ggplot(aes(x = sparsity, y = mean, ymin = mean - sd, ymax = mean + sd)) +
  theme_classic() +
  # facet_wrap(~metric) +
  ggh4x::facet_manual(~ metric, scales = "free_y",design = design) +
  geom_line() +
  geom_errorbar(width = .1) +
  geom_point() +
  scale_x_log10() +
  ylab(NULL) + xlab("Sparsity")








# Previous version: as a function of penalty
summary_metrics <- tib_quic |>
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



design <- "
 EF#
 AB#
 IG#
 DCL
 KH#
"

summary_metrics |>
  ggplot(aes(x = penalty, y = mean, ymin = mean - sd, ymax = mean + sd)) +
  theme_classic() +
  # facet_wrap(~metric) +
  ggh4x::facet_manual(~ metric, scales = "free_y",design = design) +
  geom_line() +
  geom_errorbar(width = .1) +
  geom_point() +
  scale_x_log10() +
  ylab(NULL) + xlab("Penalty (log)")

# ggsave("metrics.png", path = export_dir,
#        width = 18, height = 14, units = "cm")

summary_metrics |>
  filter(! metric == "sum_abs_residuals") |>
  ggplot(aes(x = penalty, y = mean, ymin = mean - sd, ymax = mean + sd)) +
  theme_classic() +
  ggh4x::facet_manual(~ metric, scales = "free_y",design = design) +
  geom_line() +
  geom_errorbar(width = .1) +
  geom_point() +
  scale_x_log10() +
  ylab(NULL) + xlab("Penalty (log)") +
  geom_vline(aes(xintercept = .2), linewidth = 10, alpha = .2, color = 'grey')

# ggsave("metrics_with_line.png", path = export_dir,
#        width = 18, height = 14, units = "cm")





#~ compare on same plot ----
tib_quic_old <- read_csv("data/graph_power4/outputs/240314g_glasso_PSI_npnshrink_median_2_4.csv")
tib_quic_new <- read_csv("data/graph_power4/outputs/240314g_SCIO_PSI_npnshrink_median_1_4.csv")

design <- "
 EF
 AB
 IG
 DC
 JH
"

tib_quic <- bind_rows(
  tib_quic_old |>
    add_column(run = "glasso"),
  tib_quic_new |>
    add_column(run = "SCIO")
)


summary_metrics <- tib_quic |>
  filter(permutation == 0) |> select(-permutation) |>
  mutate(`TPR/FPR` = literature_TPR/literature_FPR) |>
  summarize(across(-c(fold),
                   list(mean = partial(mean, na.rm = TRUE),
                        sd = partial(sd, na.rm = TRUE))),
            .by = c(penalty, run) ) |>
  select(- sparsity_sd) |>
  pivot_longer(-c(penalty, run, sparsity_mean),
               names_to = c("metric", "type"),
               names_pattern = "(.+)_(mean|sd|pval)$",
               values_to = "value") |>
  pivot_wider(names_from = "type",
              values_from = "value")



summary_metrics |>
  ggplot(aes(x = penalty, y = mean, ymin = mean - sd, ymax = mean + sd, color = run)) +
  theme_classic() +
  ggh4x::facet_manual(~ metric, scales = "free_y",design = design) +
  # facet_wrap(~metric) +
  geom_line() +
  geom_errorbar(width = .1) +
  geom_point() +
  scale_x_log10() +
  ylab(NULL) + xlab("Penalty (log)")

summary_metrics |>
  ggplot(aes(x = 100*(1-sparsity_mean), y = mean, ymin = mean - sd, ymax = mean + sd, color = run)) +
  theme_classic() +
  ggh4x::facet_manual(~ metric, scales = "free_y",design = design) +
  # facet_wrap(~metric) +
  geom_line() +
  geom_errorbar(width = .2) +
  geom_point() +
  ylab(NULL) + xlab("Sparsity (%)")









# ---



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
tib_quic_psi <- read_csv("data/graph_power4/outputs/240208_revertpsi_nosep_noperm_11penalties.csv")
tib_quic_counts <- read_csv("data/graph_power4/outputs/240208_recompandrevertpsi_noperm_11penalties.csv")

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






summary_metrics |>
  filter(! metric == "sum_abs_residuals") |>
  ggplot(aes(x = penalty, y = mean, ymin = mean - sd, ymax = mean + sd, color = run)) +
  theme_classic() +
  ggh4x::facet_manual(~ metric, scales = "free_y",design = design) +
  geom_line() +
  geom_errorbar(width = .1) +
  geom_point() +
  scale_x_log10() +
  ylab(NULL) + xlab("Penalty (log)")

# ggsave("metrics_counts_PSI.png", path = export_dir,
#        width = 18, height = 14, units = "cm")




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




summary_metrics |>
  filter(! metric == "sum_abs_residuals") |>
  ggplot(aes(x = penalty, y = mean, ymin = mean - sd, ymax = mean + sd, color = run)) +
  theme_classic() +
  ggh4x::facet_manual(~ metric, scales = "free_y",design = design) +
  geom_line() +
  geom_errorbar(width = .1) +
  geom_point() +
  scale_x_log10() +
  ylab(NULL) + xlab("Penalty (log)")

# ggsave("metrics_counts_normalization.png", path = export_dir,
#        width = 18, height = 14, units = "cm")






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




summary_metrics |>
  filter(! metric == "sum_abs_residuals") |>
  ggplot(aes(x = penalty, y = mean, ymin = mean - sd, ymax = mean + sd, color = run)) +
  theme_classic() +
  ggh4x::facet_manual(~ metric, scales = "free_y",design = design) +
  geom_line() +
  geom_errorbar(width = .1) +
  geom_point() +
  scale_x_log10() +
  ylab(NULL) + xlab("Penalty (log)")

# ggsave("metrics_counts_imputation.png", path = export_dir,
#        width = 18, height = 14, units = "cm")






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









# Visualize PSI and counts ----

res_quic_psi <- qs::qread("data/graph_power4/outputs/240131_revertpsi_nosep_noperm_7penalties.qs")
res_quic_counts <- qs::qread("data/graph_power4/outputs/240202_recompandrevertpsi_noperm_7penalties.qs")




# Recheck results
tib_quic_psi <- res_quic_psi |>
  select(penalty, fold, permutation, Rsquared, sum_abs_residuals,
         mean_FEV, loss_frobenius, loss_quadratic,
         prop_non_zero_coefs_litt, prop_non_zero_coefs_nonlitt)
tib_quic_counts <- res_quic_counts |>
  select(penalty, fold, permutation, Rsquared, sum_abs_residuals,
         mean_FEV, loss_frobenius, loss_quadratic,
         prop_non_zero_coefs_litt, prop_non_zero_coefs_nonlitt)

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



#~ Plot raw PSI and counts ----

res_quic_psi |> names()


plot(
res_quic_psi$psi_valid_u[[1]],
res_quic_psi$psi_valid_hat_u[[1]]
)

psi_valid <- pmap_dfr(list(res_quic_psi$penalty,
                           res_quic_psi$fold,
                            res_quic_psi$psi_valid_u,
                            res_quic_psi$psi_valid_hat_u),
                      \(.penalty, .fold, .psi_valid_u, .psi_valid_hat_u){
                        tibble(penalty = .penalty,
                               fold = .fold,
                               PSI_measured = as.numeric(.psi_valid_u),
                               PSI_hat = as.numeric(.psi_valid_hat_u))
                      })

psi_valid |>
  ggplot() +
  theme_classic() +
  geom_point(aes(x = PSI_measured, y = PSI_hat),
             alpha = .2) +
  facet_grid(rows = vars(fold), cols = vars(penalty))

psi_valid |>
  filter(fold == 1,
         penalty == .05 | penalty == 10) |>
  ggplot(aes(x = PSI_measured, y = PSI_hat)) +
  theme_classic() +
  geom_point(alpha = .2) +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~penalty) +
  xlab(expression(PSI[test] ~ (measured))) +
  ylab(expression(widehat(PSI)[test] ~ (estimated)))

# ggsave("PSI_fit.png", path = export_dir,
#        width = 15, height = 7, units = "cm")


counts_valid <- pmap_dfr(list(res_quic_counts$penalty,
                              res_quic_counts$fold,
                              res_quic_counts$psi_valid_u,
                              res_quic_counts$psi_valid_hat_u),
                         \(.penalty, .fold, .psi_valid_u, .psi_valid_hat_u){
                           tibble(penalty = .penalty,
                                  fold = .fold,
                                  PSI_measured = as.numeric(.psi_valid_u),
                                  PSI_hat = as.numeric(.psi_valid_hat_u))
                         })

counts_valid |>
  ggplot() +
  theme_classic() +
  geom_point(aes(x = PSI_measured, y = PSI_hat),
             alpha = .2) +
  facet_grid(rows = vars(fold), cols = vars(penalty)) +
  scale_x_log10() + scale_y_log10()


## Only first fold

psi_valid |>
  filter(fold == 1) |>
  ggplot() +
  theme_classic() +
  geom_point(aes(x = PSI_measured, y = PSI_hat),
             alpha = .2) +
  facet_wrap( ~penalty) +
  ggtitle("PSI")

counts_valid |>
  filter(fold == 1) |>
  ggplot() +
  theme_classic() +
  geom_point(aes(x = PSI_measured, y = PSI_hat),
             alpha = .2) +
  facet_wrap( ~penalty) +
  scale_x_log10() + scale_y_log10() +
  ggtitle("counts")




# Number of SF

res_quic_psi |> select(1:2)

xx <- res_quic_psi$adj[[6]]

dim(xx)

table(xx$coefficient == 0)


xx$sf_id[xx$coefficient != 0 ] |> unique() |> length()

xx$sf_id |> unique() |> length()



xx[xx$coefficient != 0, ] |> head()

xx |>
  summarize(nb_nonzero = sum(coefficient != 0),
            nb_total = n(),
            .by = event_id) |>
  pull(nb_nonzero) |>
  hist()








