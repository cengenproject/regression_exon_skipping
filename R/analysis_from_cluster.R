# Partly copied from `graph_power4_postanalysis.R``




# From cluster ----
library(tidyverse) |> suppressPackageStartupMessages()

source("R/analysis_helpers.R")

resdir <- "data/graph_power4/from_cluster/240315/"

fl <- list.files(resdir)





res <- tibble(run_name = str_remove(fl, "\\.csv$")) |>
  mutate(results = map(run_name,
                       read_one_res, resdir)) |>
  mutate(run_name = str_remove(run_name, "^240315[a-z]_")) |>
  unnest(cols = results) |>
  nest(.by = run_name, .key = "results")

all.equal(res$results[[2]] |> select(-permutation), res$results[[3]] |> select(-permutation))


res_noperm <- res |>
  mutate(summary_by_penalty = map(results, summarize_metrics_by_pen),
         summary_by_sparsity = map(results, summarize_metrics_by_spars))



# Ignoring permutations ----


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
  ggplot(aes(x = penalty, y = mean, ymin = mean - sd, ymax = mean + sd, color = run_name)) +
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
  ggplot(aes(x = 100*(1-sparsity_mean), y = mean, ymin = mean - sd, ymax = mean + sd, color = run_name)) +
  theme_classic() +
  ggh4x::facet_manual(~ metric, scales = "free_y",design = design) +
  # facet_wrap(~metric) +
  geom_line() +
  geom_errorbar(width = .2) +
  geom_point() +
  ylab(NULL) + xlab("Sparsity (%)")




# Permutation test ----
a <- res |>
  mutate(summary_by_penalty = map(results, summarize_metrics_by_pen2),
         pvals_res = map(results, get_pvals)) |>
  mutate(both = map2(summary_by_penalty, pvals_res,
                    ~left_join(.x, .y, by = "penalty"))) |>
  select(-c(results, summary_by_penalty, pvals_res)) |>
  unnest(cols = both) |>
  pivot_longer(-c(penalty, run_name),
               names_to = c("metric", "type"),
               names_pattern = "(.+)_(mean|sd|pval)$",
               values_to = "value") |>
  pivot_wider(names_from = "type",
              values_from = "value")

z <- a$summary_by_penalty[[1]]
s <- a$pvals_res[[1]]
left_join(z, s, by = "penalty")



#####
a |>
  separate_wider_delim(cols = run_name,
                       delim = "_",
                       names = c("date", "algo", "exonsInput", "transformation",
                                 "imputation", "nb_permutations", "nb_penalties"),
                       cols_remove = FALSE)





res |>
  unnest(cols = results)









