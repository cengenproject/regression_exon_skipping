# Partly copied from `graph_power4_postanalysis.R``




# From cluster ----
library(tidyverse) |> suppressPackageStartupMessages()

source("R/analysis_helpers.R")

resdir <- "data/graph_power4/from_cluster/noperm//"

fl <- list.files(resdir)



res <- tibble(run_name = str_remove(fl, "\\.csv$")) |>
  mutate(results = map(run_name,
                       read_one_res, resdir)) |>
  unnest(cols = results) |>
  mutate(`TPR/FPR` = (literature_TPR/literature_FPR) |> replace_na(1)) |>
  separate_wider_regex(run_name,
                       patterns = c("^(?:24[0-9]{4}[a-z]{0,2}_)?",
                                    run_algo = "[QUICSOglassoLME]+", "_",
                                    run_exonsInput = "PSI|counts", "_",
                                    run_transformation = "npnshrink|npntrunc|zscore", "_",
                                    run_imputation = "median|knn", "_",
                                    "[0-9]{1,2}", "_",
                                    run_nb_penalties = "[0-9]{1,2}"),
                       cols_remove = FALSE
                       ) |>
  nest(.by = starts_with("run_"), .key = "results")




res_noperm <- res |>
  mutate(summary_by_penalty = map(results, summarize_metrics_by_penalty),
         summary_by_sparsity = map(results, summarize_metrics_by_sparsity))


metr_levels <- c("loss_frobenius", "loss_quadratic",
                 "bias_loss_frobenius", "bias_loss_quadratic",
                 "Rsquared", "mean_FEV",
                 "literature_TPR", "literature_FPR",
                 "TPR/FPR", "power_law", "sparsity")




# No permutations ----


#~ Plotted against penalty ----

design <- "
 AB#
 CD#
 EF#
 GHI
 JK#
"

res_noperm |>
  select(-c(results, summary_by_sparsity)) |>
  unnest(cols = summary_by_penalty) |>
  mutate(metric = factor(metric, levels = metr_levels)) |>
  ggplot(aes(x = penalty, y = mean, ymin = mean - sd, ymax = mean + sd,
             color = run_algo, group = run_name)) +
  theme_classic() +
  # facet_wrap(~metric) +
  ggh4x::facet_manual(~ metric, scales = "free_y",design = design) +
  geom_line( alpha = .2) +
  # geom_errorbar(width = .1) +
  geom_point(alpha = .2) +
  scale_x_log10() +
  ylab(NULL) + xlab("Penalty (log)") +
  theme(legend.position = c(.9,.9))




#~~ manual selection of "best" penalty for each condition ----

design <- "
 AB
 CD
 EF
 G#
"


by_run <- res_noperm |>
  select(-c(results, summary_by_sparsity)) |>
  mutate(gg = map(summary_by_penalty,
                  \(tib){
                  
                    tib |>
                      filter(! startsWith(metric, "literature_"),
                             ! startsWith(metric, "bias_")) |>
                      mutate(metric = factor(metric, levels = metr_levels)) |>
                      ggplot(aes(x = penalty, y = mean, ymin = mean - sd, ymax = mean + sd)) +
                      theme_classic() +
                      ggh4x::facet_manual(~ metric, scales = "free_y",design = design) +
                      geom_line() +
                      # geom_errorbar(width = .1) +
                      geom_point() +
                      scale_x_log10() +
                      ylab(NULL) + xlab("Penalty (log)") +
                      theme(legend.position = 'none')
                  }) |>
           set_names(run_name))
  


# by_run <- by_run |> filter(run_algo == "SCIO")
for(i in seq_len(nrow(by_run))){

  ggplot <- by_run[["gg"]][[i]] +
    ggtitle(i, by_run[["run_name"]][[i]]) +
    geom_vline(aes(xintercept = .2), linewidth = 10,
               alpha = .2, color = "grey")+
    geom_errorbar(width = .05)

  ggsave(filename = paste0(i,"_",by_run[["run_name"]][[i]],".png"),
         plot = ggplot,
         path = "data/intermediates/240401_best_penalty_noperm/penalties_by_run",
         width = 18, height = 20, units = "cm")
}

manual_best_penalty <- readxl::read_excel("data/intermediates/240401_best_penalty_noperm/best_penalty.xlsx")

# same penalty for each algo
manual_best_penalty <- tibble(run_id = 1:48,
                              penalty = c(
                                rep(.2, 12),
                                rep(.2, 12),
                                rep(.2, 12),
                                rep(.2, 12)
                              ))

best_by_run <- res_noperm |>
  select(-c(results, summary_by_sparsity)) |>
  mutate(manual_penalty = manual_best_penalty$penalty) |>
  mutate(summary_by_penalty = map2(summary_by_penalty, manual_penalty,
                                    ~ {.x |> filter(penalty == .y)})) |>
  unnest(cols = summary_by_penalty)

best_by_run |>
  # mutate(mean = (mean-mean(mean))/sd(mean),
  #           .by = metric) |>
  ggplot() +
  theme_classic() +
  geom_tile(aes(x = metric, y = run_name, fill = mean)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

mat_best_by_run <- best_by_run |>
  pivot_wider(id_cols = run_name,
              names_from = "metric",
              values_from = "mean") |>
  column_to_rownames("run_name") |>
  as.matrix()
annot_best_by_run <- best_by_run |>
  select(starts_with("run_")) |>
  distinct() |>
  column_to_rownames("run_name")

mat_best_by_run <- scale(mat_best_by_run)

# ensure higher is better
mat_best_by_run[,"loss_frobenius"] <- - mat_best_by_run[,"loss_frobenius"]
mat_best_by_run[,"loss_quadratic"] <- - mat_best_by_run[,"loss_quadratic"]
mat_best_by_run[,"bias_loss_frobenius"] <- - mat_best_by_run[,"bias_loss_frobenius"]
mat_best_by_run[,"bias_loss_quadratic"] <- - mat_best_by_run[,"bias_loss_quadratic"]
mat_best_by_run[,"literature_FPR"] <- - mat_best_by_run[,"literature_FPR"]


pheatmap::pheatmap((mat_best_by_run),
                   scale = "none",
                   annotation_row = annot_best_by_run,
                   show_rownames = FALSE,
                   main = "Zscore, higher is better")


for(.metric in unique(best_by_run$metric)){
  gg <- best_by_run |>
    filter(metric == .metric) |>
    ggplot() +
    theme_classic() +
    facet_grid(rows = vars(run_transformation),
               cols = vars(run_imputation)) +
    geom_tile(aes(x = run_algo, y = run_exonsInput, fill = mean)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    ggtitle(.metric)
  
  ggsave(paste0(.metric, ".png"),
         plot = gg,
         path = "data/intermediates/240401_best_penalty_noperm/runs_by_metric",
         width = 15, height = 15, units = "cm")
}


best_by_run |>
  ggplot() +
  theme_bw() +
  facet_grid(rows = vars(metric), scales = "free_y") +
  geom_point(aes(x = interaction(run_algo, run_exonsInput, run_transformation, run_imputation), y = mean)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(axis.text.x = element_blank()) +
  scale_y_log10()

ggsave("runs_by_metric.png",
       path = "data/intermediates/240401_best_penalty_noperm/runs_by_metric/",
       width = 18, height = 35, units = "cm")



#> Proceed by elimination: First pass ----
#> 
#> First, loss Frobenius: Zscore does terrible.
#> 


res_noperm |>
  select(-c(results, summary_by_sparsity)) |>
  unnest(cols = summary_by_penalty) |>
  mutate(metric = factor(metric, levels = metr_levels)) |>
  mutate(ymin = pmax(0, mean - sd),
         ymax = mean + sd) |>
  filter(metric == "loss_frobenius") |>
  ggplot(aes(x = penalty, y = mean, ymin = ymin, ymax = ymax,
             color = run_algo)) +
  theme_classic() +
  facet_grid(rows = vars(run_exonsInput, run_imputation),
             cols = vars(run_transformation)) +
  geom_line() +
  geom_errorbar(width = .1) +
  geom_point() +
  scale_x_log10() +
  scale_y_log10(oob = scales::oob_keep) +
  ylab("Loss Frobenius (log)") + xlab("Penalty (log)")


res_noperm |>
  select(-c(results, summary_by_sparsity)) |>
  unnest(cols = summary_by_penalty) |>
  mutate(metric = factor(metric, levels = metr_levels)) |>
  mutate(ymin = pmax(0, mean - sd),
         ymax = mean + sd) |>
  filter(metric == "loss_quadratic") |>
  ggplot(aes(x = penalty, y = mean, ymin = ymin, ymax = ymax,
             color = run_algo)) +
  theme_classic() +
  facet_grid(rows = vars(run_exonsInput, run_imputation),
             cols = vars(run_transformation)) +
  geom_line() +
  geom_errorbar(width = .1) +
  geom_point() +
  scale_x_log10() +
  scale_y_log10(oob = scales::oob_keep) +
  ylab("Loss Quadratic (log)") + xlab("Penalty (log)")


res_noperm |>
  select(-c(results, summary_by_sparsity)) |>
  unnest(cols = summary_by_penalty) |>
  mutate(metric = factor(metric, levels = metr_levels)) |>
  mutate(ymin = pmax(0, mean - sd),
         ymax = mean + sd) |>
  filter(metric == "Rsquared") |>
  ggplot(aes(x = penalty, y = mean, ymin = ymin, ymax = ymax,
             color = run_algo)) +
  theme_classic() +
  facet_grid(rows = vars(run_exonsInput, run_imputation),
             cols = vars(run_transformation)) +
  geom_line() +
  geom_errorbar(width = .1) +
  geom_point() +
  scale_x_log10() +
  # scale_y_log10(oob = scales::oob_keep) +
  ylab("Rsquare") + xlab("Penalty (log)") +
  theme(legend.position = "none")

res_noperm |>
  select(-c(results, summary_by_sparsity)) |>
  unnest(cols = summary_by_penalty) |>
  mutate(metric = factor(metric, levels = metr_levels)) |>
  mutate(ymin = pmax(0, mean - sd),
         ymax = mean + sd) |>
  filter(metric == "mean_FEV") |>
  ggplot(aes(x = penalty, y = mean, ymin = ymin, ymax = ymax,
             color = run_algo)) +
  theme_classic() +
  facet_grid(rows = vars(run_exonsInput, run_imputation),
             cols = vars(run_transformation)) +
  geom_line() +
  geom_errorbar(width = .1) +
  geom_point() +
  scale_x_log10() +
  # scale_y_log10(oob = scales::oob_keep) +
  ylab("FEV") + xlab("Penalty (log)") +
  theme(legend.position = "none")



res_noperm |>
  select(-c(results, summary_by_sparsity)) |>
  unnest(cols = summary_by_penalty) |>
  mutate(metric = factor(metric, levels = metr_levels)) |>
  mutate(ymin = pmax(0, mean - sd),
         ymax = mean + sd) |>
  filter(metric == "TPR/FPR") |>
  ggplot(aes(x = penalty, y = mean, ymin = ymin, ymax = ymax,
             color = run_algo)) +
  theme_bw() +
  facet_grid(rows = vars(run_exonsInput, run_imputation),
             cols = vars(run_transformation),
             scales = "free_y") +
  geom_line() +
  geom_errorbar(width = .1) +
  geom_point() +
  scale_x_log10(limits = c(.04, .9)) +
  # scale_y_log10(oob = scales::oob_keep) +
  ylab("TPR/FPR") + xlab("Penalty (log)") +
  theme(legend.position = "none")

res_noperm |>
  select(-c(results, summary_by_sparsity)) |>
  unnest(cols = summary_by_penalty) |>
  mutate(metric = factor(metric, levels = metr_levels)) |>
  mutate(ymin = pmax(0, mean - sd),
         ymax = mean + sd) |>
  filter(metric == "power_law") |>
  ggplot(aes(x = penalty, y = mean, ymin = ymin, ymax = ymax,
             color = run_algo)) +
  theme_bw() +
  facet_grid(rows = vars(run_exonsInput, run_imputation),
             cols = vars(run_transformation)) +
  geom_line() +
  geom_errorbar(width = .1) +
  geom_point() +
  scale_x_log10(limits = c(.04, .9)) +
  # scale_y_log10(oob = scales::oob_keep) +
  ylab("Power law") + xlab("Penalty (log)") +
  theme(legend.position = "none")



#> Second pass ----
#> Eliminate SCIO and Zscore
#> Focus on imputation method


res_noperm |>
  select(-c(results, summary_by_sparsity)) |>
  unnest(cols = summary_by_penalty) |>
  mutate(metric = factor(metric, levels = metr_levels)) |>
  filter(run_transformation != "zscore", run_algo %in% c("glasso", "CLIME")) |>
  mutate(ymin = pmax(0, mean - sd),
         ymax = mean + sd) |>
  filter(metric == "loss_frobenius") |>
  ggplot(aes(x = penalty, y = mean, ymin = ymin, ymax = ymax,
             shape = run_imputation, color = run_algo)) +
  theme_bw() +
  facet_grid(rows = vars(run_exonsInput),
             cols = vars(run_transformation)) +
  geom_line() +
  geom_errorbar(width = .1) +
  geom_point() +
  scale_x_log10() +
  # scale_y_log10(oob = scales::oob_keep) +
  ylab("Loss Frobenius") + xlab("Penalty (log)") +
  theme(legend.position = "none")

res_noperm |>
  select(-c(results, summary_by_sparsity)) |>
  unnest(cols = summary_by_penalty) |>
  mutate(metric = factor(metric, levels = metr_levels)) |>
  filter(run_transformation != "zscore", run_algo %in% c("glasso", "CLIME")) |>
  mutate(ymin = pmax(0, mean - sd),
         ymax = mean + sd) |>
  filter(metric == "loss_quadratic") |>
  ggplot(aes(x = penalty, y = mean, ymin = ymin, ymax = ymax,
             shape = run_imputation, color = run_algo)) +
  theme_bw() +
  facet_grid(rows = vars(run_exonsInput),
             cols = vars(run_transformation)) +
  geom_line() +
  geom_errorbar(width = .1) +
  geom_point() +
  scale_x_log10() +
  scale_y_log10(oob = scales::oob_keep) +
  ylab("Loss Quadratic (log)") + xlab("Penalty (log)") +
  theme(legend.position = "none")


res_noperm |>
  select(-c(results, summary_by_sparsity)) |>
  unnest(cols = summary_by_penalty) |>
  mutate(metric = factor(metric, levels = metr_levels)) |>
  filter(run_transformation != "zscore", run_algo %in% c("glasso", "CLIME")) |>
  mutate(ymin = pmax(0, mean - sd),
         ymax = mean + sd) |>
  filter(metric == "Rsquared") |>
  ggplot(aes(x = penalty, y = mean, ymin = ymin, ymax = ymax,
             shape = run_imputation, color = run_algo)) +
  theme_bw() +
  facet_grid(rows = vars(run_exonsInput),
             cols = vars(run_transformation)) +
  geom_line() +
  geom_errorbar(width = .1) +
  geom_point() +
  scale_x_log10() +
  # scale_y_log10(oob = scales::oob_keep) +
  ylab("R squared") + xlab("Penalty (log)") +
  theme(legend.position = "none")

res_noperm |>
  select(-c(results, summary_by_sparsity)) |>
  unnest(cols = summary_by_penalty) |>
  mutate(metric = factor(metric, levels = metr_levels)) |>
  filter(run_transformation != "zscore", run_algo %in% c("glasso", "CLIME")) |>
  mutate(ymin = pmax(0, mean - sd),
         ymax = mean + sd) |>
  filter(metric == "mean_FEV") |>
  ggplot(aes(x = penalty, y = mean, ymin = ymin, ymax = ymax,
             shape = run_imputation, color = run_algo)) +
  theme_bw() +
  facet_grid(rows = vars(run_exonsInput),
             cols = vars(run_transformation)) +
  geom_line() +
  geom_errorbar(width = .1) +
  geom_point() +
  scale_x_log10() +
  # scale_y_log10(oob = scales::oob_keep) +
  ylab("FEV") + xlab("Penalty (log)") +
  theme(legend.position = "none")

res_noperm |>
  select(-c(results, summary_by_sparsity)) |>
  unnest(cols = summary_by_penalty) |>
  mutate(metric = factor(metric, levels = metr_levels)) |>
  filter(run_transformation != "zscore", run_algo %in% c("glasso", "CLIME")) |>
  mutate(ymin = pmax(0, mean - sd),
         ymax = mean + sd) |>
  filter(metric == "TPR/FPR") |>
  ggplot(aes(x = penalty, y = mean, ymin = ymin, ymax = ymax,
             shape = run_imputation, color = run_algo)) +
  theme_bw() +
  facet_grid(rows = vars(run_exonsInput),
             cols = vars(run_transformation)) +
  geom_line() +
  geom_errorbar(width = .05) +
  geom_point() +
  scale_x_log10(limits = c(.04, .9)) +
  scale_y_continuous(limits = c(0,10), oob = scales::oob_keep) +
  # scale_y_log10(oob = scales::oob_keep) +
  ylab("TPR/FPR") + xlab("Penalty (log)") +
  theme(legend.position = "none")

res_noperm |>
  select(-c(results, summary_by_sparsity)) |>
  unnest(cols = summary_by_penalty) |>
  mutate(metric = factor(metric, levels = metr_levels)) |>
  filter(run_transformation != "zscore", run_algo %in% c("glasso", "CLIME")) |>
  mutate(ymin = pmax(0, mean - sd),
         ymax = mean + sd) |>
  filter(metric == "power_law") |>
  ggplot(aes(x = penalty, y = mean, ymin = ymin, ymax = ymax,
             shape = run_imputation, color = run_algo)) +
  theme_bw() +
  facet_grid(rows = vars(run_exonsInput),
             cols = vars(run_transformation)) +
  geom_line() +
  geom_errorbar(width = .1) +
  geom_point() +
  scale_x_log10(limits = c(.04, .9)) +
  # scale_y_log10(oob = scales::oob_keep) +
  ylab("Power law") + xlab("Penalty (log)") +
  theme(legend.position = "none")






#> Third pass ----
#> Focus on exonsInput method


res_noperm |>
  select(-c(results, summary_by_sparsity)) |>
  unnest(cols = summary_by_penalty) |>
  mutate(metric = factor(metric, levels = metr_levels)) |>
  filter(run_transformation != "zscore", run_algo %in% c("glasso", "CLIME")) |>
  mutate(ymin = pmax(0, mean - sd),
         ymax = mean + sd) |>
  filter(metric == "loss_frobenius") |>
  ggplot(aes(x = penalty, y = mean, ymin = ymin, ymax = ymax,
             shape = run_exonsInput, color = run_algo)) +
  theme_bw() +
  facet_grid(rows = vars(run_imputation),
             cols = vars(run_transformation)) +
  geom_line() +
  geom_errorbar(width = .1) +
  geom_point() +
  scale_x_log10() +
  # scale_y_log10(oob = scales::oob_keep) +
  ylab("Loss Frobenius") + xlab("Penalty (log)") +
  theme(legend.position = "none")

res_noperm |>
  select(-c(results, summary_by_sparsity)) |>
  unnest(cols = summary_by_penalty) |>
  mutate(metric = factor(metric, levels = metr_levels)) |>
  filter(run_transformation != "zscore", run_algo %in% c("glasso", "CLIME")) |>
  mutate(ymin = pmax(0, mean - sd),
         ymax = mean + sd) |>
  filter(metric == "loss_quadratic") |>
  ggplot(aes(x = penalty, y = mean, ymin = ymin, ymax = ymax,
             shape = run_exonsInput, color = run_algo)) +
  theme_bw() +
  facet_grid(rows = vars(run_imputation),
             cols = vars(run_transformation)) +
  geom_line() +
  geom_errorbar(width = .1) +
  geom_point() +
  scale_x_log10() +
  scale_y_log10(oob = scales::oob_keep) +
  ylab("Loss Quadratic (log)") + xlab("Penalty (log)") +
  theme(legend.position = "none")


res_noperm |>
  select(-c(results, summary_by_sparsity)) |>
  unnest(cols = summary_by_penalty) |>
  mutate(metric = factor(metric, levels = metr_levels)) |>
  filter(run_transformation != "zscore", run_algo %in% c("glasso", "CLIME")) |>
  mutate(ymin = pmax(0, mean - sd),
         ymax = mean + sd) |>
  filter(metric == "Rsquared") |>
  ggplot(aes(x = penalty, y = mean, ymin = ymin, ymax = ymax,
             shape = run_exonsInput, color = run_algo)) +
  theme_bw() +
  facet_grid(rows = vars(run_imputation),
             cols = vars(run_transformation)) +
  geom_line() +
  geom_errorbar(width = .1) +
  geom_point() +
  scale_x_log10() +
  # scale_y_log10(oob = scales::oob_keep) +
  ylab("R squared") + xlab("Penalty (log)") +
  theme(legend.position = "none")

res_noperm |>
  select(-c(results, summary_by_sparsity)) |>
  unnest(cols = summary_by_penalty) |>
  mutate(metric = factor(metric, levels = metr_levels)) |>
  filter(run_transformation != "zscore", run_algo %in% c("glasso", "CLIME")) |>
  mutate(ymin = pmax(0, mean - sd),
         ymax = mean + sd) |>
  filter(metric == "mean_FEV") |>
  ggplot(aes(x = penalty, y = mean, ymin = ymin, ymax = ymax,
             shape = run_exonsInput, color = run_algo)) +
  theme_bw() +
  facet_grid(rows = vars(run_imputation),
             cols = vars(run_transformation)) +
  geom_line() +
  geom_errorbar(width = .1) +
  geom_point() +
  scale_x_log10() +
  # scale_y_log10(oob = scales::oob_keep) +
  ylab("FEV") + xlab("Penalty (log)") +
  theme(legend.position = "none")

res_noperm |>
  select(-c(results, summary_by_sparsity)) |>
  unnest(cols = summary_by_penalty) |>
  mutate(metric = factor(metric, levels = metr_levels)) |>
  filter(run_transformation != "zscore", run_algo %in% c("glasso", "CLIME")) |>
  mutate(ymin = pmax(0, mean - sd),
         ymax = mean + sd) |>
  filter(metric == "TPR/FPR") |>
  ggplot(aes(x = penalty, y = mean, ymin = ymin, ymax = ymax,
             shape = run_exonsInput, color = run_algo)) +
  theme_bw() +
  facet_grid(rows = vars(run_imputation),
             cols = vars(run_transformation)) +
  geom_line() +
  geom_errorbar(width = .03) +
  geom_point() +
  scale_x_log10(limits = c(.04, .9)) +
  scale_y_continuous(limits = c(0,10), oob = scales::oob_keep) +
  # scale_y_log10(oob = scales::oob_keep) +
  ylab("TPR/FPR") + xlab("Penalty (log)") +
  theme(legend.position = "none")

res_noperm |>
  select(-c(results, summary_by_sparsity)) |>
  unnest(cols = summary_by_penalty) |>
  mutate(metric = factor(metric, levels = metr_levels)) |>
  filter(run_transformation != "zscore", run_algo %in% c("glasso", "CLIME")) |>
  mutate(ymin = pmax(0, mean - sd),
         ymax = mean + sd) |>
  filter(metric == "power_law") |>
  ggplot(aes(x = penalty, y = mean, ymin = ymin, ymax = ymax,
             shape = run_exonsInput, color = run_algo)) +
  theme_bw() +
  facet_grid(rows = vars(run_imputation),
             cols = vars(run_transformation)) +
  geom_line() +
  geom_errorbar(width = .1) +
  geom_point() +
  scale_x_log10(limits = c(.04, .9)) +
  # scale_y_log10(oob = scales::oob_keep) +
  ylab("Power law") + xlab("Penalty (log)") +
  theme(legend.position = "none")


#> Fourth pass ----
#> Focus on transformation method


res_noperm |>
  select(-c(results, summary_by_sparsity)) |>
  unnest(cols = summary_by_penalty) |>
  mutate(metric = factor(metric, levels = metr_levels)) |>
  filter(run_transformation != "zscore", run_algo == "glasso") |>
  mutate(ymin = pmax(0, mean - sd),
         ymax = mean + sd) |>
  filter(metric == "loss_frobenius") |>
  ggplot(aes(x = penalty, y = mean, ymin = ymin, ymax = ymax,
             color = run_transformation)) +
  theme_bw() +
  facet_grid(rows = vars(run_imputation),
             cols = vars(run_exonsInput)) +
  geom_line() +
  geom_errorbar(width = .1) +
  geom_point() +
  scale_x_log10() +
  # scale_y_log10(oob = scales::oob_keep) +
  ylab("Loss Frobenius") + xlab("Penalty (log)") +
  theme(legend.position = "none")

res_noperm |>
  select(-c(results, summary_by_sparsity)) |>
  unnest(cols = summary_by_penalty) |>
  mutate(metric = factor(metric, levels = metr_levels)) |>
  filter(run_transformation != "zscore", run_algo == "glasso") |>
  mutate(ymin = pmax(0, mean - sd),
         ymax = mean + sd) |>
  filter(metric == "loss_quadratic") |>
  ggplot(aes(x = penalty, y = mean, ymin = ymin, ymax = ymax,
             color = run_transformation)) +
  theme_bw() +
  facet_grid(rows = vars(run_imputation),
             cols = vars(run_exonsInput)) +
  geom_line() +
  geom_errorbar(width = .1) +
  geom_point() +
  scale_x_log10() +
  # scale_y_log10(oob = scales::oob_keep) +
  ylab("Loss Quadratic") + xlab("Penalty (log)") +
  theme(legend.position = "none")


res_noperm |>
  select(-c(results, summary_by_sparsity)) |>
  unnest(cols = summary_by_penalty) |>
  mutate(metric = factor(metric, levels = metr_levels)) |>
  filter(run_transformation != "zscore", run_algo == "glasso") |>
  mutate(ymin = pmax(0, mean - sd),
         ymax = mean + sd) |>
  filter(metric == "Rsquared") |>
  ggplot(aes(x = penalty, y = mean, ymin = ymin, ymax = ymax,
             color = run_transformation)) +
  theme_bw() +
  facet_grid(rows = vars(run_imputation),
             cols = vars(run_exonsInput)) +
  geom_line() +
  geom_errorbar(width = .1) +
  geom_point() +
  scale_x_log10() +
  # scale_y_log10(oob = scales::oob_keep) +
  ylab("R squared") + xlab("Penalty (log)") +
  theme(legend.position = "none")

res_noperm |>
  select(-c(results, summary_by_sparsity)) |>
  unnest(cols = summary_by_penalty) |>
  mutate(metric = factor(metric, levels = metr_levels)) |>
  filter(run_transformation != "zscore", run_algo == "glasso") |>
  mutate(ymin = pmax(0, mean - sd),
         ymax = mean + sd) |>
  filter(metric == "mean_FEV") |>
  ggplot(aes(x = penalty, y = mean, ymin = ymin, ymax = ymax,
             color = run_transformation)) +
  theme_bw() +
  facet_grid(rows = vars(run_imputation),
             cols = vars(run_exonsInput)) +
  geom_line() +
  geom_errorbar(width = .1) +
  geom_point() +
  scale_x_log10() +
  # scale_y_log10(oob = scales::oob_keep) +
  ylab("FEV") + xlab("Penalty (log)") +
  theme(legend.position = "none")

res_noperm |>
  select(-c(results, summary_by_sparsity)) |>
  unnest(cols = summary_by_penalty) |>
  mutate(metric = factor(metric, levels = metr_levels)) |>
  filter(run_transformation != "zscore", run_algo == "glasso") |>
  mutate(ymin = pmax(0, mean - sd),
         ymax = mean + sd) |>
  filter(metric == "TPR/FPR") |>
  ggplot(aes(x = penalty, y = mean, ymin = ymin, ymax = ymax,
             color = run_transformation)) +
  theme_bw() +
  facet_grid(rows = vars(run_imputation),
             cols = vars(run_exonsInput)) +
  geom_line() +
  geom_errorbar(width = .03) +
  geom_point() +
  scale_x_log10(limits = c(.04, .9)) +
  scale_y_continuous(limits = c(0,4), oob = scales::oob_keep) +
  # scale_y_log10(oob = scales::oob_keep) +
  ylab("TPR/FPR") + xlab("Penalty (log)") +
  theme(legend.position = "none")

res_noperm |>
  select(-c(results, summary_by_sparsity)) |>
  unnest(cols = summary_by_penalty) |>
  mutate(metric = factor(metric, levels = metr_levels)) |>
  filter(run_transformation != "zscore", run_algo == "glasso") |>
  mutate(ymin = pmax(0, mean - sd),
         ymax = mean + sd) |>
  filter(metric == "power_law") |>
  ggplot(aes(x = penalty, y = mean, ymin = ymin, ymax = ymax,
             color = run_transformation)) +
  theme_bw() +
  facet_grid(rows = vars(run_imputation),
             cols = vars(run_exonsInput)) +
  geom_line() +
  geom_errorbar(width = .1) +
  geom_point() +
  scale_x_log10(limits = c(.04, .9)) +
  # scale_y_log10(oob = scales::oob_keep) +
  ylab("Power law") + xlab("Penalty (log)") +
  theme(legend.position = "none")















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








