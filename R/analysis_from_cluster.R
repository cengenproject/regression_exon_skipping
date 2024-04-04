# Partly copied from `graph_power4_postanalysis.R``

# Inits ----

library(tidyverse) |> suppressPackageStartupMessages()

source("R/analysis_helpers.R")



# various conditions, no permutation ----


resdir <- "data/graph_power4/from_cluster/240326/"

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



# knn k ----


resdirknn <- "data/graph_power4/from_cluster/240403_reknn/"

flknn <- list.files(resdirknn)



resknn <- tibble(run_name = str_remove(flknn, "\\.csv$")) |>
  mutate(results = map(run_name,
                       read_one_res, resdirknn)) |>
  unnest(cols = results) |>
  mutate(`TPR/FPR` = (literature_TPR/literature_FPR) |> replace_na(1)) |>
  separate_wider_regex(run_name,
                       patterns = c("^(?:24[0-9]{4}[a-z]{0,2}_)?",
                                    run_algo = "[QUICSOglassoLME]+", "_",
                                    run_exonsInput = "PSI|counts", "_",
                                    run_transformation = "npnshrink|npntrunc|zscore", "_",
                                    run_imputation = "median|knn", "_k",
                                    run_k = "[0-9]+", "_",
                                    "[0-9]{1,2}", "_",
                                    run_nb_penalties = "[0-9]{1,2}"),
                       cols_remove = FALSE
  ) |>
  nest(.by = starts_with("run_"), .key = "results")




res_knn <- resknn |>
  mutate(summary_by_penalty = map(results, summarize_metrics_by_penalty),
         summary_by_sparsity = map(results, summarize_metrics_by_sparsity))




#> Fifth pass ----
#> 
#> Selected at fixed penalty

res_knn |>
  select(-c(results, summary_by_sparsity)) |>
  unnest(cols = summary_by_penalty) |>
  mutate(run_k = as.integer(run_k),
         penalty = as.factor(penalty)) |>
  mutate(ymin = pmax(0, mean - sd),
         ymax = mean + sd) |>
  filter(metric == "loss_frobenius",
         penalty == .2 | penalty == .3) |>
  ggplot(aes(x = run_k, y = mean, ymin = ymin, ymax = ymax,
             color = penalty, shape = penalty)) +
  theme_classic() +
  facet_grid(rows = vars(run_transformation)) +
  geom_errorbar(width = .1) +
  geom_point() +
  geom_line() +
  # scale_x_log10() +
  # scale_y_log10(oob = scales::oob_keep) +
  ylab("Loss Frobenius") + xlab("k")


res_knn |>
  select(-c(results, summary_by_sparsity)) |>
  unnest(cols = summary_by_penalty) |>
  mutate(run_k = as.integer(run_k),
         penalty = as.factor(penalty)) |>
  mutate(ymin = pmax(0, mean - sd),
         ymax = mean + sd) |>
  filter(metric == "loss_quadratic",
         penalty == .2 | penalty == .3) |>
  ggplot(aes(x = run_k, y = mean, ymin = ymin, ymax = ymax,
             color = penalty, shape = penalty)) +
  theme_classic() +
  facet_grid(rows = vars(run_transformation)) +
  geom_line() +
  geom_errorbar(width = .1) +
  geom_point() +
  # scale_x_log10() +
  # scale_y_log10(oob = scales::oob_keep) +
  ylab("Loss Quadratic") + xlab("k") +
  theme(legend.position = "none")


res_knn |>
  select(-c(results, summary_by_sparsity)) |>
  unnest(cols = summary_by_penalty) |>
  mutate(run_k = as.integer(run_k),
         penalty = as.factor(penalty)) |>
  mutate(ymin = pmax(0, mean - sd),
         ymax = mean + sd) |>
  filter(metric == "Rsquared",
         penalty == .2 | penalty == .3) |>
  ggplot(aes(x = run_k, y = mean, ymin = ymin, ymax = ymax,
             color = penalty, shape = penalty)) +
  theme_classic() +
  facet_grid(rows = vars(run_transformation)) +
  geom_line() +
  geom_errorbar(width = .1) +
  geom_point() +
  # scale_x_log10() +
  # scale_y_log10(oob = scales::oob_keep) +
  ylab("Rsquare") + xlab("k") +
  theme(legend.position = "none")

res_knn |>
  select(-c(results, summary_by_sparsity)) |>
  unnest(cols = summary_by_penalty) |>
  mutate(run_k = as.integer(run_k),
         penalty = as.factor(penalty)) |>
  mutate(ymin = pmax(0, mean - sd),
         ymax = mean + sd) |>
  filter(metric == "mean_FEV",
         penalty == .2 | penalty == .3) |>
  ggplot(aes(x = run_k, y = mean, ymin = ymin, ymax = ymax,
             color = penalty, shape = penalty)) +
  theme_classic() +
  facet_grid(rows = vars(run_transformation)) +
  geom_line() +
  geom_errorbar(width = .1) +
  geom_point() +
  # scale_x_log10() +
  # scale_y_log10(oob = scales::oob_keep) +
  ylab("FEV") + xlab("k") +
  theme(legend.position = "none")



res_knn |>
  select(-c(results, summary_by_sparsity)) |>
  unnest(cols = summary_by_penalty) |>
  mutate(run_k = as.integer(run_k),
         penalty = as.factor(penalty)) |>
  mutate(ymin = pmax(0, mean - sd),
         ymax = mean + sd) |>
  filter(metric == "TPR/FPR",
         penalty == .2 | penalty == .3) |>
  ggplot(aes(x = run_k, y = mean, ymin = ymin, ymax = ymax,
             color = penalty, shape = penalty)) +
  theme_bw() +
  facet_grid(rows = vars(run_transformation)) +
  geom_line() +
  geom_errorbar(width = .1) +
  geom_point() +
  # scale_x_log10(limits = c(.04, .9)) +
  # scale_y_log10(oob = scales::oob_keep) +
  ylab("TPR/FPR") + xlab("k") +
  theme(legend.position = "none")

res_knn |>
  select(-c(results, summary_by_sparsity)) |>
  unnest(cols = summary_by_penalty) |>
  mutate(run_k = as.integer(run_k),
         penalty = as.factor(penalty)) |>
  mutate(ymin = pmax(0, mean - sd),
         ymax = mean + sd) |>
  filter(metric == "power_law",
         penalty == .2 | penalty == .3) |>
  ggplot(aes(x = run_k, y = mean, ymin = ymin, ymax = ymax,
             color = penalty, shape = penalty)) +
  theme_bw() +
  facet_grid(rows = vars(run_transformation)) +
  geom_line() +
  geom_errorbar(width = .1) +
  geom_point() +
  # scale_x_log10(limits = c(.04, .9)) +
  # scale_y_log10(oob = scales::oob_keep) +
  ylab("Power law") + xlab("k") +
  theme(legend.position = "none")



#> Sixth pass ----
#> 
#> Selected at fixed penalty



res_knn |>
  select(-c(results, summary_by_sparsity)) |>
  unnest(cols = summary_by_penalty) |>
  filter(penalty == .3,
         run_k == 10,
         ! startsWith(metric, "literature"),
         metric != "sparsity") |>
  ggplot(aes(x = run_transformation, y = mean, ymin = mean - sd, ymax = mean + sd,
             color = run_transformation)) +
  theme_classic() +
  facet_wrap(~metric, scales = "free_y") +
  # ggh4x::facet_manual(~ metric, scales = "free_y",design = design) +
  # geom_line() +
  geom_errorbar(width = .1) +
  geom_point() +
  # scale_x_log10() +
  ylab(NULL) + xlab("Metric")










#~ Plotted against sparsity ----

design <- "
 EF
 AB
 IG
 DC
 JH
"


res_knn |>
  select(-c(results, summary_by_penalty)) |>
  unnest(cols = summary_by_sparsity) |>
  filter(run_transformation == "npnshrink") |>
  ggplot(aes(x = 100*(1-sparsity_mean),
             y = mean,
             ymin = mean - sd,
             ymax = mean + sd,
             color = run_k)) +
  theme_classic() +
  ggh4x::facet_manual(~ metric, scales = "free_y",design = design) +
  # facet_wrap(~metric) +
  geom_line() +
  geom_errorbar(width = .2) +
  geom_point() +
  ylab(NULL) + xlab("Sparsity (%)")

res_knn |>
  select(-c(results, summary_by_penalty)) |>
  unnest(cols = summary_by_sparsity) |>
  filter(run_transformation == "npntrunc") |>
  ggplot(aes(x = 100*(1-sparsity_mean),
             y = mean,
             ymin = mean - sd,
             ymax = mean + sd,
             color = run_k)) +
  theme_classic() +
  ggh4x::facet_manual(~ metric, scales = "free_y",design = design) +
  # facet_wrap(~metric) +
  geom_line() +
  geom_errorbar(width = .2) +
  geom_point() +
  ylab(NULL) + xlab("Sparsity (%)")






# Final with Permutations ----


resdirfin <- "data/graph_power4/from_cluster/final_perm/"

flfin <- list.files(resdirfin)



resfin <- tibble(run_name = str_remove(flfin, "\\.csv$")) |>
  mutate(results = map(run_name,
                       read_one_res, resdirfin)) |>
  unnest(cols = results) |>
  mutate(`TPR/FPR` = (literature_TPR/literature_FPR) |> replace_na(1)) |>
  separate_wider_regex(run_name,
                       patterns = c("^(?:2404[0-9]{2}p[0-9]{0,4}_)?",
                                    run_algo = "[QUICSOglassoLME]+", "_",
                                    run_exonsInput = "PSI|counts", "_",
                                    run_transformation = "npnshrink|npntrunc|zscore", "_",
                                    run_imputation = "median|knn", "_k",
                                    run_k = "[0-9]+", "_",
                                    "[0-9]{1,2}", "_",
                                    run_nb_penalties = "[0-9]{1,2}"),
                       cols_remove = TRUE
  ) |>
  nest(.by = starts_with("run_"), .key = "results")



# # Diagnostics

# expect 1 fold * 16 penalties = 16 values, for 1000+1 permutations
sapply(resfin$results,
       \(.x) table(.x[["permutation"]]) |> table())

resfin |>
  # filter(run_algo == "glasso", run_exonsInput == "PSI",
  #        run_transformation == "npnshrink", run_imputation=="knn") |>
  unnest(results) |>
  pull(permutation) |>
  table() |>
  table()

# Check that permutations not identical
res_by_run <- tibble(run_name = str_remove(flfin, "\\.csv$")) |>
  mutate(results = map(run_name,
                       read_one_res, resdirfin))

if(identical(
  res_by_run$results[[2]] |> select(-permutation),
  res_by_run$results[[3]] |> select(-permutation)
)){
  stop("the second and third run are identical, permutations were not random")
}



# check p-val distribution
resfin |>
  filter(run_algo == "glasso", run_exonsInput == "PSI",
         run_transformation == "npnshrink", run_imputation=="knn") |>
  unnest(results) |>
  # filter(penalty == 0.3, fold == 2) |>
  nest(.by = starts_with("run_")) |>
  mutate(summary_by_penalty = map(data, get_pvals)) |>
  select(-data) |>
  unnest(summary_by_penalty) |>
  pull(pval) |>
  table()


res_perm <- resfin |>
  mutate(summary_by_penalty = map(results, summarize_metrics_by_penalty),
         summary_by_sparsity = map(results, summarize_metrics_by_sparsity),
         pvals_res = map(results, get_pvals)) |>
  mutate(summary_by_penalty = map2(summary_by_penalty, pvals_res,
                                   ~left_join(.x, .y, by = c("penalty", "metric"))),
         summary_by_sparsity = map2(summary_by_sparsity, pvals_res,
                                    ~left_join(.x, .y, by = c("penalty", "metric")))) |>
  select(-c(results, pvals_res))


# Plot individual metrics permutations ----

for(metr in unique(res_perm$summary_by_penalty[[1]]$metric)){
  gg <- resfin |>
    unnest(results) |>
    select(-c(starts_with("run_"))) |>
    ggplot() +
    ggbeeswarm::geom_quasirandom(aes(y = .data[[metr]], x = penalty,
                                     color = permutation == 0,
                                     alpha = permutation == 0,
                                     size = permutation == 0,
                                     shape = permutation == 0)) +
    scale_x_log10() +
    # ggtitle(metr) +
    scale_alpha_manual(values = c(0.1,1)) +
    scale_size_manual(values = c(0.5,1.5)) +
    theme(legend.position = "none")
  
  ggsave(paste0("pval_", metr, ".png"), plot = gg,
         path = "data/intermediates/240404_perm/",
         width = 12, height = 7, units = "cm")
  
}

gg




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
  select(-starts_with("run_")) |>
  ggplot(aes(x = penalty,
             y = mean,
             ymin = mean - sd, ymax = mean + sd)) +
  theme_classic() +
  # facet_wrap(~metric) +
  ggh4x::facet_manual(~ metric, scales = "free_y",design = design) +
  geom_line() +
  geom_errorbar(width = .1) +
  geom_point(aes(shape = pval_adj < .05, color = pval_adj < .05)) +
  scale_x_log10() +
  ylab(NULL) + xlab("Penalty (log)") +
  scale_color_manual(values = c("grey", 'red4')) +
  theme(legend.position = c(.9,.9))


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
  select(-starts_with("run_")) |>
  ggplot(aes(x = 100*(1-sparsity_mean),
             y = mean,
             ymin = mean - sd, ymax = mean + sd)) +
  theme_classic() +
  # facet_wrap(~metric) +
  ggh4x::facet_manual(~ metric, scales = "free_y",design = design) +
  geom_line() +
  geom_errorbar(width = .001) +
  geom_point(aes(shape = pval_adj < .05, color = pval_adj < .05)) +
  scale_x_log10() +
  ylab(NULL) + xlab("Sparsity (%)") +
  scale_size_manual(values = c(.5,2.5)) +
  scale_color_manual(values = c("grey", 'red4'))



design <- "
AC
ED
FB
"

res_perm |>
  unnest(summary_by_sparsity) |>
  select(-starts_with("run_")) |>
  filter(! startsWith(metric, "bias_"),
         ! startsWith(metric, "literature")) |>
  mutate(metric = case_match(
    metric,
    "loss_frobenius" ~ "Frobenius loss",
    "loss_quadratic" ~ "Quadratic loss",
    "Rsquared" ~ "Reconstruction R squared",
    "mean_FEV" ~ "Reconstruction Fraction explained variance",
    "TPR/FPR" ~ "TPR/FPR",
    "power_law" ~ "Network Power law",
  )) |>
  ggplot(aes(x = 100*(1-sparsity_mean),
             y = mean,
             ymin = mean - sd, ymax = mean + sd)) +
  theme_classic() +
  # facet_wrap(~metric) +
  ggh4x::facet_manual(~ metric, scales = "free_y",design = design) +
  geom_line() +
  geom_errorbar(width = .01) +
  geom_point(aes(shape = pval_adj < .05,
                 color = pval_adj < .05),
             size = 2) +
  # scale_x_log10() +
  ylab(NULL) + xlab("Sparsity (%)") +
  scale_color_manual(values = c("grey", 'red4')) +
  theme(legend.position = 'none')



res_perm |>
  unnest(summary_by_sparsity) |>
  select(-starts_with("run_")) |>
  filter(! startsWith(metric, "bias_"),
         ! startsWith(metric, "literature")) |>
  mutate(category = case_match(
    metric,
    "loss_frobenius" ~ "Global structure",
    "loss_quadratic" ~ "Global structure",
    "Rsquared" ~ "PSI reconstruction",
    "mean_FEV" ~ "PSI reconstruction",
    "TPR/FPR" ~ "Network",
    "power_law" ~ "Network",
  ) |>
    factor(levels = c("Global structure", "PSI reconstruction", "Network")),
  metric = case_match(
    metric,
    "loss_frobenius" ~ "Frobenius loss",
    "loss_quadratic" ~ "Quadratic loss",
    "Rsquared" ~ "R squared",
    "mean_FEV" ~ "Fraction explained variance",
    "TPR/FPR" ~ "TPR/FPR",
    "power_law" ~ "Power law",
  ) |>
    factor(levels = c("Frobenius loss", "Quadratic loss",
                      "R squared", "Fraction explained variance",
                      "TPR/FPR", "Power law"))) |>
  ggplot(aes(x = 100*(1-sparsity_mean),
             y = mean,
             ymin = mean - sd, ymax = mean + sd)) +
  theme_bw() +
  facet_grid(cols = vars(metric), rows = vars(category),
             scales = "free_y") +
  # ggh4x::facet_manual(~ metric, scales = "free_y",design = design) +
  ggh4x::facet_nested_wrap(vars(category, metric),
                           ncol = 2, nrow = 3,
                           scales = "free_y") +
  geom_line() +
  geom_errorbar(width = .01) +
  geom_point(aes(shape = pval_adj < .05,
                 color = pval_adj < .05),
             size = 2) +
  geom_vline(aes(xintercept = 97.9), linewidth = 5, alpha = .2) +
  # scale_x_log10() +
  ylab(NULL) + xlab("Sparsity (%)") +
  scale_color_manual(values = c("grey", 'red4')) +
  theme(legend.position = 'none')




#~ Selected case ----

res <- qs::qread("data/graph_power4/from_cluster/final_save/240404_glasso_PSI_npnshrink_knn_k10_2_16.qs")

final_tib <- read_one_res("240404_glasso_PSI_npnshrink_knn_k10_2_16",
                          "data/graph_power4/from_cluster/final_save/")

# One with and one without permutations
final_tib$permutation |> table()
final_tib$fold |> table()

design <- "
 EF#
 AB#
 IG#
 DCK
 JH#
"

# re-check non-permuted (same as above but no pvals)
final_tib |>
  mutate(`TPR/FPR` = (literature_TPR/literature_FPR) |> replace_na(1)) |>
  nest() |>
  mutate(data = map(data, summarize_metrics_by_penalty)) |>
  unnest(data) |>
  ggplot(aes(x = penalty,
             y = mean,
             ymin = mean - sd, ymax = mean + sd)) +
  theme_classic() +
  # facet_wrap(~metric) +
  ggh4x::facet_manual(~ metric, scales = "free_y",design = design) +
  geom_line() +
  geom_errorbar(width = .1) +
  geom_point() +
  scale_x_log10() +
  ylab(NULL) + xlab("Penalty (log)") +
  scale_color_manual(values = c("grey", 'red4')) +
  theme(legend.position = c(.9,.9))

# look at (double check) the permuted case
final_tib |>
  mutate(`TPR/FPR` = (literature_TPR/literature_FPR) |> replace_na(1)) |>
  filter(permutation ==1) |> mutate(permutation = 0) |>
  nest() |>
  mutate(data = map(data, summarize_metrics_by_penalty)) |>
  unnest(data) |>
  ggplot(aes(x = penalty,
             y = mean,
             ymin = mean - sd, ymax = mean + sd)) +
  theme_classic() +
  # facet_wrap(~metric) +
  ggh4x::facet_manual(~ metric, scales = "free_y",design = design) +
  geom_line() +
  geom_errorbar(width = .1) +
  geom_point() +
  scale_x_log10() +
  ylab(NULL) + xlab("Penalty (log)") +
  scale_color_manual(values = c("grey", 'red4')) +
  theme(legend.position = c(.9,.9))














