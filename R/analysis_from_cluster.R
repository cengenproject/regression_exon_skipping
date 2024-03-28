# Partly copied from `graph_power4_postanalysis.R``




# From cluster ----
library(tidyverse) |> suppressPackageStartupMessages()

source("R/analysis_helpers.R")

resdir <- "data/graph_power4/from_cluster/240326/"

fl <- list.files(resdir)



res <- tibble(run_name = str_remove(fl, "\\.csv$")) |>
  mutate(results = map(run_name,
                       read_one_res, resdir)) |>
  unnest(cols = results) |>
  mutate(`TPR/FPR` = (literature_TPR/literature_FPR) |> replace_na(1)) |>
  separate_wider_regex(run_name,
                       patterns = c("^(?:24[0-9]{4}[a-z]{1,2}_)?",
                                    run_algo = "[QUICSOglasso]+", "_",
                                    run_exonsInput = "PSI", "_",
                                    run_transformation = "(?:npnshrink|npntrunc)", "_",
                                    run_imputation = "median", "_",
                                    "[0-9]{1,2}", "_",
                                    run_nb_penalties = "[0-9]{1,2}")
                       ) |>
  nest(.by = starts_with("run_"), .key = "results")




res_noperm <- res |>
  mutate(summary_by_penalty = map(results, summarize_metrics_by_penalty),
         summary_by_sparsity = map(results, summarize_metrics_by_sparsity))



# No permutations ----


#~ Plotted against penalty ----

design <- "
 EF#
 AB#
 IG#
 DCK
 JH#
"

res_noperm |>
  select(-c(results, summary_by_sparsity)) |>
  unnest(cols = summary_by_penalty) |>
  ggplot(aes(x = penalty, y = mean, ymin = mean - sd, ymax = mean + sd, color = run_algo)) +
  theme_classic() +
  # facet_wrap(~metric) +
  ggh4x::facet_manual(~ metric, scales = "free_y",design = design) +
  geom_line() +
  geom_errorbar(width = .1) +
  geom_point() +
  scale_x_log10() +
  ylab(NULL) + xlab("Penalty (log)")

#~ Plotted against sparsity ----

design <- "
 EF
 AB
 IG
 DC
 JH
"


res_noperm |>
  select(-c(results, summary_by_penalty)) |>
  unnest(cols = summary_by_sparsity) |>
  ggplot(aes(x = 100*(1-sparsity_mean),
             y = mean,
             ymin = mean - sd,
             ymax = mean + sd,
             color = run_algo)) +
  theme_classic() +
  ggh4x::facet_manual(~ metric, scales = "free_y",design = design) +
  # facet_wrap(~metric) +
  geom_line() +
  geom_errorbar(width = .2) +
  geom_point() +
  ylab(NULL) + xlab("Sparsity (%)")






# with Permutations ----



# # Diagnostics

# expect 5 folds * 11 penalties = 55, for 201 permutations
sapply(res$results,
       \(.x) table(.x[["permutation"]]) |> table())

# setdiff(0:200,
#         res$results[[4]]$permutation |> unique() |> sort())

#
# table(res$results[[3]]$permutation) |> table()
# res2 <- tibble(run_name = str_remove(fl, "\\.csv$")) |>
#   mutate(results = map(run_name,
#                        read_one_res, resdir)) |>
#   unnest(cols = results) |>
#   mutate(`TPR/FPR` = (literature_TPR/literature_FPR) |> replace_na(1)) |>
#   separate_wider_regex(run_name,
#                        patterns = c("^(?:24[0-9]{4}[1-z]{1,2}_)?",
#                                     run_algo = "[QUICSOglasso]+", "_",
#                                     run_exonsInput = "PSI", "_",
#                                     run_transformation = "npnshrink", "_",
#                                     run_imputation = "median", "_",
#                                     "[0-9]{1,2}", "_",
#                                     run_nb_penalties = "[0-9]{1,2}"),
#                        cols_remove = FALSE
#   )
# 
# res2 |> filter(run_algo == "glasso",
#                permutation == 54,
#                penalty == 10) |> pull(run_name)


res_perm <- res |>
  mutate(summary_by_penalty = map(results, summarize_metrics_by_penalty),
         summary_by_sparsity = map(results, summarize_metrics_by_sparsity),
         pvals_res = map(results, get_pvals)) |>
  mutate(summary_by_penalty = map2(summary_by_penalty, pvals_res,
                    ~left_join(.x, .y, by = c("penalty", "metric"))),
         summary_by_sparsity = map2(summary_by_sparsity, pvals_res,
                                    ~left_join(.x, .y, by = c("penalty", "metric")))) |>
  select(-c(results, pvals_res))



#~ Plotted against penalty ----
design <- "
 EF#
 AB#
 IG#
 DCK
 JH#
"

res_perm |>
  unnest(summary_by_penalty) |>
  select(-c(run_exonsInput, run_transformation,
            run_imputation, run_nb_penalties,
            summary_by_sparsity)) |>
  ggplot(aes(x = penalty,
             y = mean,
             ymin = mean - sd, ymax = mean + sd,
             color = run_algo)) +
  theme_classic() +
  # facet_wrap(~metric) +
  ggh4x::facet_manual(~ metric, scales = "free_y",design = design) +
  geom_line() +
  geom_errorbar(width = .1) +
  geom_point(aes(shape = pval < .05, size = pval < .05)) +
  scale_x_log10() +
  ylab(NULL) + xlab("Penalty (log)") +
  scale_size_manual(values = c(.5,2.5))


#~ Plotted against sparsity ----
design <- "
 EF
 AB
 IG
 DC
 JH
"
res_perm |>
  unnest(summary_by_sparsity) |>
  select(-c(run_exonsInput, run_transformation,
            run_imputation, run_nb_penalties,
            summary_by_penalty)) |>
  ggplot(aes(x = 100*(1-sparsity_mean),
             y = mean,
             ymin = mean - sd, ymax = mean + sd,
             color = run_algo)) +
  theme_classic() +
  # facet_wrap(~metric) +
  ggh4x::facet_manual(~ metric, scales = "free_y",design = design) +
  geom_line() +
  geom_errorbar(width = .001) +
  geom_point(aes(shape = pval < .05, size = pval < .05)) +
  scale_x_log10() +
  ylab(NULL) + xlab("Sparsity (%)") +
  scale_size_manual(values = c(.5,2.5))









#~ Selected algo ----



res_sparsity_perm |>
  filter(run_algo == "glasso") |>
  select(-c(run_exonsInput, run_transformation,
            run_imputation, run_nb_penalties)) |>
  ggplot(aes(x = 100*(1-(sparsity_mean)),
             y = mean,
             ymin = mean - sd, ymax = mean + sd)) +
  theme_classic() +
  # facet_wrap(~metric) +
  ggh4x::facet_manual(~ metric, scales = "free_y",design = design) +
  geom_line() +
  geom_errorbar(width = .001) +
  geom_point(aes(shape = pval < .05, color = pval < .05)) +
  # scale_x_log10() +
  scale_x_continuous(trans = \(x){1-x/100})
  ylab(NULL) + xlab("Sparsity (%)")







# just glasso at specific penalty vs permutations

res$results[[1]]$penalty |> unique()

res$results[[1]] |>
  filter(penalty == .2) |>
  select(permutation, fold, `TPR/FPR`) |>
  ggplot() +
  ggbeeswarm::geom_quasirandom(aes(y = `TPR/FPR`, x = fold, color = permutation == 0))

# power law
res$results[[1]] |>
  filter(penalty == .2) |>
  select(permutation, fold, power_law) |>
  ggplot() +
  ggbeeswarm::geom_quasirandom(aes(y = power_law, x = fold, color = permutation == 0))


# SCIO FEV
res$results[[3]] |>
  filter(penalty == .2) |>
  select(permutation, fold, mean_FEV) |>
  ggplot() +
  ggbeeswarm::geom_quasirandom(aes(y = mean_FEV, x = fold, color = permutation == 0))








