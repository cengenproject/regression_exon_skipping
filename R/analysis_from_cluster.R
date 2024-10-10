# Partly copied from `graph_power4_postanalysis.R``

# Inits ----

library(tidyverse) |> suppressPackageStartupMessages()

source("R/analysis_helpers.R")



# various conditions, no permutation ----


resdir <- "data/graph_power4/from_cluster/240429_noperm/"

fl <- list.files(resdir)



res <- tibble(run_name = str_remove(fl, "\\.csv$")) |>
  mutate(results = map(run_name,
                       read_one_res, resdir)) |>
  unnest(cols = results) |>
  mutate(`TPR/FPR` = (literature_TPR/literature_FPR) |> replace_na(1),
         run_name = str_replace(run_name, "knn_k10", "knn")) |>
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
                 "loss_frobenius_adj", "bias_loss_frobenius_adj",
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
         path = "data/intermediates/240429_best_penalty_noperm/penalties_by_run",
         width = 18, height = 20, units = "cm")
}


# same penalty for each algo
# manual_best_penalty <- tibble(run_id = 1:48,
#                               penalty = c(
#                                 rep(.2, 12),
#                                 rep(.2, 12),
#                                 rep(.2, 12),
#                                 rep(.2, 12)
#                               ))

manual_best_penalty <- tibble(run_id = 1:36,
                              penalty = c(
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
  # filter(run_algo != "QUIC") |>
  pivot_wider(id_cols = run_name,
              names_from = "metric",
              values_from = "mean") |>
  column_to_rownames("run_name") |>
  as.matrix()
annot_best_by_run <- best_by_run |>
  # filter(run_algo != "QUIC") |>
  select(starts_with("run_")) |>
  distinct() |>
  column_to_rownames("run_name") |>
  mutate(run_transformation = recode(run_transformation,
                                     npnshrink = "NPN (shrunken)",
                                     npntrunc = "NPN (truncated)",
                                     zscore = "Z-score"))

mat_best_by_run <- scale(mat_best_by_run)

# ensure higher is better
mat_best_by_run[,"loss_frobenius"] <- - mat_best_by_run[,"loss_frobenius"]
mat_best_by_run[,"loss_frobenius_adj"] <- - mat_best_by_run[,"loss_frobenius_adj"]
mat_best_by_run[,"loss_quadratic"] <- - mat_best_by_run[,"loss_quadratic"]
mat_best_by_run[,"bias_loss_frobenius"] <- - mat_best_by_run[,"bias_loss_frobenius"]
mat_best_by_run[,"bias_loss_quadratic"] <- - mat_best_by_run[,"bias_loss_quadratic"]
mat_best_by_run[,"literature_FPR"] <- - mat_best_by_run[,"literature_FPR"]



# pheatmap::pheatmap((mat_best_by_run[,more_useful_metrics]),
#                    scale = "none",
#                    annotation_row = annot_best_by_run,
#                    show_rownames = FALSE,
#                    main = "Zscore, higher is better")
# 
# 
# pheatmap::pheatmap((mat_best_by_run[,more_useful_metrics]),
#                    scale = "none",
#                    annotation_row = annot_best_by_run |>
#                      select(imputation = run_imputation,
#                             transformation = run_transformation,
#                             input = run_exonsInput,
#                             algorithm = run_algo),
#                    show_rownames = FALSE,
#                    main = "Zscore, higher is better")
# 
# 
# ComplexHeatmap::pheatmap(mat_best_by_run[,more_useful_metrics],
#                          annotation_row = annot_best_by_run |>
#                            select(imputation = run_imputation,
#                                   transformation = run_transformation,
#                                   input = run_exonsInput,
#                                   algorithm = run_algo),
#                          cluster_cols = FALSE,
#                          show_rownames = FALSE)


row_annot <- ComplexHeatmap::HeatmapAnnotation(
  df = annot_best_by_run |>
    select(algorithm = run_algo,
           `input format` = run_exonsInput,
           imputation = run_imputation,
           transformation = run_transformation),
  which = "row",
  col =  list(algorithm = c(CLIME = '#C2AD4B', glasso = '#418C82',
                            QUIC = '#00BCBD', SCIO = '#C55E2D'),
              `input format` = c(counts = "#7fc97f", PSI = "#beaed4"),
              imputation = c(knn = "#b3cde3", median = "#ccebc5"),
              transformation = c(`NPN (shrunken)` = "#dfc27d",
                                 `NPN (truncated)` = "#a6611a",
                                 `Z-score` = "#80cdc1"))
)


more_useful_metrics <- c("loss_frobenius_adj",
                         "mean_FEV",
                         "TPR/FPR", "power_law")

ComplexHeatmap::Heatmap(mat_best_by_run[,more_useful_metrics],
                        left_annotation = row_annot,
                        show_row_names = FALSE,
                        rect_gp = grid::gpar(col="grey"),
                        cluster_columns = FALSE,
                        column_labels = c("Frobenius loss\n(inverted)",
                                          "Fraction\nExplained Variance",
                                          "TPR / FPR", "Power Law"),
                        name = "Z-score")

pdf("presentations/240308_figures/heatmap.pdf",
    width = 8, height = 8)
dev.off()

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
  
  # ggsave(paste0(.metric, ".png"),
  #        plot = gg,
  #        path = "data/intermediates/240415_best_penalty_noperm/runs_by_metric",
  #        width = 15, height = 15, units = "cm")
}


best_by_run |>
  ggplot() +
  theme_bw() +
  facet_grid(rows = vars(metric), scales = "free_y") +
  geom_point(aes(x = interaction(run_algo, run_exonsInput, run_transformation, run_imputation), y = mean)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(axis.text.x = element_blank()) +
  scale_y_log10()

# ggsave("runs_by_metric.png",
#        path = "data/intermediates/240401_best_penalty_noperm/runs_by_metric/",
#        width = 18, height = 35, units = "cm")



# Old version: keeping because can be more readable ----

# #> Proceed by elimination: First pass ---
# #> 
# #> First, loss Frobenius: Zscore does terrible.
# #> 
# 
# 
# res_noperm |>
#   select(-c(results, summary_by_sparsity)) |>
#   unnest(cols = summary_by_penalty) |>
#   mutate(metric = factor(metric, levels = metr_levels)) |>
#   mutate(ymin = pmax(0, mean - sd),
#          ymax = mean + sd) |>
#   filter(metric == "loss_frobenius") |>
#   ggplot(aes(x = penalty, y = mean, ymin = ymin, ymax = ymax,
#              color = run_algo)) +
#   theme_classic() +
#   facet_grid(rows = vars(run_exonsInput, run_imputation),
#              cols = vars(run_transformation)) +
#   geom_line() +
#   geom_errorbar(width = .1) +
#   geom_point() +
#   scale_x_log10() +
#   scale_y_log10(oob = scales::oob_keep) +
#   ylab("Loss Frobenius (log)") + xlab("Penalty (log)")
# 
# 
# res_noperm |>
#   select(-c(results, summary_by_sparsity)) |>
#   unnest(cols = summary_by_penalty) |>
#   mutate(metric = factor(metric, levels = metr_levels)) |>
#   mutate(ymin = pmax(0, mean - sd),
#          ymax = mean + sd) |>
#   filter(metric == "loss_quadratic") |>
#   ggplot(aes(x = penalty, y = mean, ymin = ymin, ymax = ymax,
#              color = run_algo)) +
#   theme_classic() +
#   facet_grid(rows = vars(run_exonsInput, run_imputation),
#              cols = vars(run_transformation)) +
#   geom_line() +
#   geom_errorbar(width = .1) +
#   geom_point() +
#   scale_x_log10() +
#   scale_y_log10(oob = scales::oob_keep) +
#   ylab("Loss Quadratic (log)") + xlab("Penalty (log)")
# 
# 
# 
# res_noperm |>
#   select(-c(results, summary_by_sparsity)) |>
#   unnest(cols = summary_by_penalty) |>
#   mutate(metric = factor(metric, levels = metr_levels)) |>
#   mutate(ymin = pmax(0, mean - sd),
#          ymax = mean + sd) |>
#   filter(metric == "loss_frobenius_adj") |>
#   ggplot(aes(x = penalty, y = mean, ymin = ymin, ymax = ymax,
#              color = run_algo)) +
#   theme_classic() +
#   facet_grid(rows = vars(run_exonsInput, run_imputation),
#              cols = vars(run_transformation)) +
#   geom_line() +
#   geom_errorbar(width = .1) +
#   geom_point() +
#   scale_x_log10() +
#   # scale_y_log10(oob = scales::oob_keep) +
#   ylab("Frob adjacency") + xlab("Penalty (log)") +
#   theme(legend.position = "none")

# 
# res_noperm |>
#   select(-c(results, summary_by_sparsity)) |>
#   unnest(cols = summary_by_penalty) |>
#   mutate(metric = factor(metric, levels = metr_levels)) |>
#   mutate(ymin = pmax(0, mean - sd),
#          ymax = mean + sd) |>
#   filter(metric == "Rsquared") |>
#   ggplot(aes(x = penalty, y = mean, ymin = ymin, ymax = ymax,
#              color = run_algo)) +
#   theme_classic() +
#   facet_grid(rows = vars(run_exonsInput, run_imputation),
#              cols = vars(run_transformation)) +
#   geom_line() +
#   geom_errorbar(width = .1) +
#   geom_point() +
#   scale_x_log10() +
#   # scale_y_log10(oob = scales::oob_keep) +
#   ylab("Rsquare") + xlab("Penalty (log)") +
#   theme(legend.position = "none")
# 
# res_noperm |>
#   select(-c(results, summary_by_sparsity)) |>
#   unnest(cols = summary_by_penalty) |>
#   mutate(metric = factor(metric, levels = metr_levels)) |>
#   mutate(ymin = pmax(0, mean - sd),
#          ymax = mean + sd) |>
#   filter(metric == "mean_FEV") |>
#   ggplot(aes(x = penalty, y = mean, ymin = ymin, ymax = ymax,
#              color = run_algo)) +
#   theme_classic() +
#   facet_grid(rows = vars(run_exonsInput, run_imputation),
#              cols = vars(run_transformation)) +
#   geom_line() +
#   geom_errorbar(width = .1) +
#   geom_point() +
#   scale_x_log10() +
#   # scale_y_log10(oob = scales::oob_keep) +
#   ylab("FEV") + xlab("Penalty (log)") +
#   theme(legend.position = "none")
# 
# 
# 
# res_noperm |>
#   select(-c(results, summary_by_sparsity)) |>
#   unnest(cols = summary_by_penalty) |>
#   mutate(metric = factor(metric, levels = metr_levels)) |>
#   mutate(ymin = pmax(0, mean - sd),
#          ymax = mean + sd) |>
#   filter(metric == "TPR/FPR") |>
#   ggplot(aes(x = penalty, y = mean, ymin = ymin, ymax = ymax,
#              color = run_algo)) +
#   theme_bw() +
#   facet_grid(rows = vars(run_exonsInput, run_imputation),
#              cols = vars(run_transformation),
#              scales = "free_y") +
#   geom_line() +
#   geom_errorbar(width = .1) +
#   geom_point() +
#   scale_x_log10(limits = c(.04, .9)) +
#   # scale_y_log10(oob = scales::oob_keep) +
#   ylab("TPR/FPR") + xlab("Penalty (log)") +
#   theme(legend.position = "none")
# 
# res_noperm |>
#   select(-c(results, summary_by_sparsity)) |>
#   unnest(cols = summary_by_penalty) |>
#   mutate(metric = factor(metric, levels = metr_levels)) |>
#   mutate(ymin = pmax(0, mean - sd),
#          ymax = mean + sd) |>
#   filter(metric == "power_law") |>
#   ggplot(aes(x = penalty, y = mean, ymin = ymin, ymax = ymax,
#              color = run_algo)) +
#   theme_bw() +
#   facet_grid(rows = vars(run_exonsInput, run_imputation),
#              cols = vars(run_transformation)) +
#   geom_line() +
#   geom_errorbar(width = .1) +
#   geom_point() +
#   scale_x_log10(limits = c(.04, .9)) +
#   # scale_y_log10(oob = scales::oob_keep) +
#   ylab("Power law") + xlab("Penalty (log)") +
#   theme(legend.position = "none")
# 
# 
# 
# #> Second pass ---
# #> Eliminate SCIO and Zscore
# #> Focus on imputation method
# 
# 
# res_noperm |>
#   select(-c(results, summary_by_sparsity)) |>
#   unnest(cols = summary_by_penalty) |>
#   mutate(metric = factor(metric, levels = metr_levels)) |>
#   filter(run_transformation != "zscore", run_algo %in% c("glasso")) |>
#   mutate(ymin = pmax(0, mean - sd),
#          ymax = mean + sd) |>
#   filter(metric == "loss_frobenius") |>
#   ggplot(aes(x = penalty, y = mean, ymin = ymin, ymax = ymax,
#              shape = run_imputation, color = run_algo)) +
#   theme_bw() +
#   facet_grid(rows = vars(run_exonsInput),
#              cols = vars(run_transformation)) +
#   geom_line() +
#   geom_errorbar(width = .1) +
#   geom_point() +
#   scale_x_log10() +
#   # scale_y_log10(oob = scales::oob_keep) +
#   ylab("Loss Frobenius") + xlab("Penalty (log)") +
#   theme(legend.position = "none")
# 
# res_noperm |>
#   select(-c(results, summary_by_sparsity)) |>
#   unnest(cols = summary_by_penalty) |>
#   mutate(metric = factor(metric, levels = metr_levels)) |>
#   filter(run_transformation != "zscore", run_algo %in% c("glasso", "CLIME")) |>
#   mutate(ymin = pmax(0, mean - sd),
#          ymax = mean + sd) |>
#   filter(metric == "loss_quadratic") |>
#   ggplot(aes(x = penalty, y = mean, ymin = ymin, ymax = ymax,
#              shape = run_imputation, color = run_algo)) +
#   theme_bw() +
#   facet_grid(rows = vars(run_exonsInput),
#              cols = vars(run_transformation)) +
#   geom_line() +
#   geom_errorbar(width = .1) +
#   geom_point() +
#   scale_x_log10() +
#   scale_y_log10(oob = scales::oob_keep) +
#   ylab("Loss Quadratic (log)") + xlab("Penalty (log)") +
#   theme(legend.position = "none")
# 
# 
# res_noperm |>
#   select(-c(results, summary_by_sparsity)) |>
#   unnest(cols = summary_by_penalty) |>
#   mutate(metric = factor(metric, levels = metr_levels)) |>
#   filter(run_transformation != "zscore", run_algo %in% c("glasso", "CLIME")) |>
#   mutate(ymin = pmax(0, mean - sd),
#          ymax = mean + sd) |>
#   filter(metric == "loss_frobenius_adj") |>
#   ggplot(aes(x = penalty, y = mean, ymin = ymin, ymax = ymax,
#              shape = run_imputation, color = run_algo)) +
#   theme_classic() +
#   facet_grid(rows = vars(run_exonsInput),
#              cols = vars(run_transformation)) +
#   geom_line() +
#   geom_errorbar(width = .1) +
#   geom_point() +
#   scale_x_log10() +
#   # scale_y_log10(oob = scales::oob_keep) +
#   ylab("Frob adjacency") + xlab("Penalty (log)") +
#   theme(legend.position = "none")
# 
# 
# res_noperm |>
#   select(-c(results, summary_by_sparsity)) |>
#   unnest(cols = summary_by_penalty) |>
#   mutate(metric = factor(metric, levels = metr_levels)) |>
#   filter(run_transformation != "zscore", run_algo %in% c("glasso", "CLIME")) |>
#   mutate(ymin = pmax(0, mean - sd),
#          ymax = mean + sd) |>
#   filter(metric == "Rsquared") |>
#   ggplot(aes(x = penalty, y = mean, ymin = ymin, ymax = ymax,
#              shape = run_imputation, color = run_algo)) +
#   theme_bw() +
#   facet_grid(rows = vars(run_exonsInput),
#              cols = vars(run_transformation)) +
#   geom_line() +
#   geom_errorbar(width = .1) +
#   geom_point() +
#   scale_x_log10() +
#   # scale_y_log10(oob = scales::oob_keep) +
#   ylab("R squared") + xlab("Penalty (log)") +
#   theme(legend.position = "none")
# 
# res_noperm |>
#   select(-c(results, summary_by_sparsity)) |>
#   unnest(cols = summary_by_penalty) |>
#   mutate(metric = factor(metric, levels = metr_levels)) |>
#   filter(run_transformation != "zscore", run_algo %in% c("glasso", "CLIME")) |>
#   mutate(ymin = pmax(0, mean - sd),
#          ymax = mean + sd) |>
#   filter(metric == "mean_FEV") |>
#   ggplot(aes(x = penalty, y = mean, ymin = ymin, ymax = ymax,
#              shape = run_imputation, color = run_algo)) +
#   theme_bw() +
#   facet_grid(rows = vars(run_exonsInput),
#              cols = vars(run_transformation)) +
#   geom_line() +
#   geom_errorbar(width = .1) +
#   geom_point() +
#   scale_x_log10() +
#   # scale_y_log10(oob = scales::oob_keep) +
#   ylab("FEV") + xlab("Penalty (log)") +
#   theme(legend.position = "none")
# 
# res_noperm |>
#   select(-c(results, summary_by_sparsity)) |>
#   unnest(cols = summary_by_penalty) |>
#   mutate(metric = factor(metric, levels = metr_levels)) |>
#   filter(run_transformation != "zscore", run_algo %in% c("glasso", "CLIME")) |>
#   mutate(ymin = pmax(0, mean - sd),
#          ymax = mean + sd) |>
#   filter(metric == "TPR/FPR") |>
#   ggplot(aes(x = penalty, y = mean, ymin = ymin, ymax = ymax,
#              shape = run_imputation, color = run_algo)) +
#   theme_bw() +
#   facet_grid(rows = vars(run_exonsInput),
#              cols = vars(run_transformation)) +
#   geom_line() +
#   geom_errorbar(width = .05) +
#   geom_point() +
#   scale_x_log10(limits = c(.04, .9)) +
#   scale_y_continuous(limits = c(0,10), oob = scales::oob_keep) +
#   # scale_y_log10(oob = scales::oob_keep) +
#   ylab("TPR/FPR") + xlab("Penalty (log)") +
#   theme(legend.position = "none")
# 
# res_noperm |>
#   select(-c(results, summary_by_sparsity)) |>
#   unnest(cols = summary_by_penalty) |>
#   mutate(metric = factor(metric, levels = metr_levels)) |>
#   filter(run_transformation != "zscore", run_algo %in% c("glasso", "CLIME")) |>
#   mutate(ymin = pmax(0, mean - sd),
#          ymax = mean + sd) |>
#   filter(metric == "power_law") |>
#   ggplot(aes(x = penalty, y = mean, ymin = ymin, ymax = ymax,
#              shape = run_imputation, color = run_algo)) +
#   theme_bw() +
#   facet_grid(rows = vars(run_exonsInput),
#              cols = vars(run_transformation)) +
#   geom_line() +
#   geom_errorbar(width = .1) +
#   geom_point() +
#   scale_x_log10(limits = c(.04, .9)) +
#   # scale_y_log10(oob = scales::oob_keep) +
#   ylab("Power law") + xlab("Penalty (log)") +
#   theme(legend.position = "none")
# 
# 
# 
# 
# 
# 
# #> Third pass ---
# #> Focus on exonsInput method
# 
# 
# res_noperm |>
#   select(-c(results, summary_by_sparsity)) |>
#   unnest(cols = summary_by_penalty) |>
#   mutate(metric = factor(metric, levels = metr_levels)) |>
#   filter(run_transformation != "zscore", run_algo %in% c("glasso")) |>
#   mutate(ymin = pmax(0, mean - sd),
#          ymax = mean + sd) |>
#   filter(metric == "loss_frobenius") |>
#   ggplot(aes(x = penalty, y = mean, ymin = ymin, ymax = ymax,
#              shape = run_exonsInput, color = run_algo)) +
#   theme_bw() +
#   facet_grid(rows = vars(run_imputation),
#              cols = vars(run_transformation)) +
#   geom_line() +
#   geom_errorbar(width = .1) +
#   geom_point() +
#   scale_x_log10() +
#   # scale_y_log10(oob = scales::oob_keep) +
#   ylab("Loss Frobenius") + xlab("Penalty (log)") +
#   theme(legend.position = "none")
# 
# res_noperm |>
#   select(-c(results, summary_by_sparsity)) |>
#   unnest(cols = summary_by_penalty) |>
#   mutate(metric = factor(metric, levels = metr_levels)) |>
#   filter(run_transformation != "zscore", run_algo %in% c("glasso", "CLIME")) |>
#   mutate(ymin = pmax(0, mean - sd),
#          ymax = mean + sd) |>
#   filter(metric == "loss_quadratic") |>
#   ggplot(aes(x = penalty, y = mean, ymin = ymin, ymax = ymax,
#              shape = run_exonsInput, color = run_algo)) +
#   theme_bw() +
#   facet_grid(rows = vars(run_imputation),
#              cols = vars(run_transformation)) +
#   geom_line() +
#   geom_errorbar(width = .1) +
#   geom_point() +
#   scale_x_log10() +
#   scale_y_log10(oob = scales::oob_keep) +
#   ylab("Loss Quadratic (log)") + xlab("Penalty (log)") +
#   theme(legend.position = "none")
# 
# 
# 
# res_noperm |>
#   select(-c(results, summary_by_sparsity)) |>
#   unnest(cols = summary_by_penalty) |>
#   mutate(metric = factor(metric, levels = metr_levels)) |>
#   filter(run_transformation != "zscore", run_algo %in% c("glasso", "CLIME")) |>
#   mutate(ymin = pmax(0, mean - sd),
#          ymax = mean + sd) |>
#   filter(metric == "loss_frobenius_adj") |>
#   ggplot(aes(x = penalty, y = mean, ymin = ymin, ymax = ymax,
#              shape = run_exonsInput, color = run_algo)) +
#   theme_classic() +
#   facet_grid(rows = vars(run_imputation),
#              cols = vars(run_transformation)) +
#   geom_line() +
#   geom_errorbar(width = .1) +
#   geom_point() +
#   scale_x_log10() +
#   # scale_y_log10(oob = scales::oob_keep) +
#   ylab("Frob adjacency") + xlab("Penalty (log)") +
#   theme(legend.position = "none")
# 
# res_noperm |>
#   select(-c(results, summary_by_sparsity)) |>
#   unnest(cols = summary_by_penalty) |>
#   mutate(metric = factor(metric, levels = metr_levels)) |>
#   filter(run_transformation != "zscore", run_algo %in% c("glasso", "CLIME")) |>
#   mutate(ymin = pmax(0, mean - sd),
#          ymax = mean + sd) |>
#   filter(metric == "Rsquared") |>
#   ggplot(aes(x = penalty, y = mean, ymin = ymin, ymax = ymax,
#              shape = run_exonsInput, color = run_algo)) +
#   theme_bw() +
#   facet_grid(rows = vars(run_imputation),
#              cols = vars(run_transformation)) +
#   geom_line() +
#   geom_errorbar(width = .1) +
#   geom_point() +
#   scale_x_log10() +
#   # scale_y_log10(oob = scales::oob_keep) +
#   ylab("R squared") + xlab("Penalty (log)") +
#   theme(legend.position = "none")
# 
# res_noperm |>
#   select(-c(results, summary_by_sparsity)) |>
#   unnest(cols = summary_by_penalty) |>
#   mutate(metric = factor(metric, levels = metr_levels)) |>
#   filter(run_transformation != "zscore", run_algo %in% c("glasso", "CLIME")) |>
#   mutate(ymin = pmax(0, mean - sd),
#          ymax = mean + sd) |>
#   filter(metric == "mean_FEV") |>
#   ggplot(aes(x = penalty, y = mean, ymin = ymin, ymax = ymax,
#              shape = run_exonsInput, color = run_algo)) +
#   theme_bw() +
#   facet_grid(rows = vars(run_imputation),
#              cols = vars(run_transformation)) +
#   geom_line() +
#   geom_errorbar(width = .1) +
#   geom_point() +
#   scale_x_log10() +
#   # scale_y_log10(oob = scales::oob_keep) +
#   ylab("FEV") + xlab("Penalty (log)") +
#   theme(legend.position = "none")
# 
# res_noperm |>
#   select(-c(results, summary_by_sparsity)) |>
#   unnest(cols = summary_by_penalty) |>
#   mutate(metric = factor(metric, levels = metr_levels)) |>
#   filter(run_transformation != "zscore", run_algo %in% c("glasso", "CLIME")) |>
#   mutate(ymin = pmax(0, mean - sd),
#          ymax = mean + sd) |>
#   filter(metric == "TPR/FPR") |>
#   ggplot(aes(x = penalty, y = mean, ymin = ymin, ymax = ymax,
#              shape = run_exonsInput, color = run_algo)) +
#   theme_bw() +
#   facet_grid(rows = vars(run_imputation),
#              cols = vars(run_transformation)) +
#   geom_line() +
#   geom_errorbar(width = .03) +
#   geom_point() +
#   scale_x_log10(limits = c(.04, .9)) +
#   scale_y_continuous(limits = c(0,10), oob = scales::oob_keep) +
#   # scale_y_log10(oob = scales::oob_keep) +
#   ylab("TPR/FPR") + xlab("Penalty (log)") +
#   theme(legend.position = "none")
# 
# res_noperm |>
#   select(-c(results, summary_by_sparsity)) |>
#   unnest(cols = summary_by_penalty) |>
#   mutate(metric = factor(metric, levels = metr_levels)) |>
#   filter(run_transformation != "zscore", run_algo %in% c("glasso", "CLIME")) |>
#   mutate(ymin = pmax(0, mean - sd),
#          ymax = mean + sd) |>
#   filter(metric == "power_law") |>
#   ggplot(aes(x = penalty, y = mean, ymin = ymin, ymax = ymax,
#              shape = run_exonsInput, color = run_algo)) +
#   theme_bw() +
#   facet_grid(rows = vars(run_imputation),
#              cols = vars(run_transformation)) +
#   geom_line() +
#   geom_errorbar(width = .1) +
#   geom_point() +
#   scale_x_log10(limits = c(.04, .9)) +
#   # scale_y_log10(oob = scales::oob_keep) +
#   ylab("Power law") + xlab("Penalty (log)") +
#   theme(legend.position = "none")
# 
# 
# 
# #> Fourth pass ---
# #> Focus on transformation method
# 
# 
# res_noperm |>
#   select(-c(results, summary_by_sparsity)) |>
#   unnest(cols = summary_by_penalty) |>
#   mutate(metric = factor(metric, levels = metr_levels)) |>
#   filter(run_transformation != "zscore", run_algo == "glasso") |>
#   mutate(ymin = pmax(0, mean - sd),
#          ymax = mean + sd) |>
#   filter(metric == "loss_frobenius") |>
#   ggplot(aes(x = penalty, y = mean, ymin = ymin, ymax = ymax,
#              color = run_transformation)) +
#   theme_bw() +
#   facet_grid(rows = vars(run_imputation),
#              cols = vars(run_exonsInput)) +
#   geom_line() +
#   geom_errorbar(width = .1) +
#   geom_point() +
#   scale_x_log10() +
#   # scale_y_log10(oob = scales::oob_keep) +
#   ylab("Loss Frobenius") + xlab("Penalty (log)") +
#   theme(legend.position = "none")
# 
# res_noperm |>
#   select(-c(results, summary_by_sparsity)) |>
#   unnest(cols = summary_by_penalty) |>
#   mutate(metric = factor(metric, levels = metr_levels)) |>
#   filter(run_transformation != "zscore", run_algo == "glasso") |>
#   mutate(ymin = pmax(0, mean - sd),
#          ymax = mean + sd) |>
#   filter(metric == "loss_quadratic") |>
#   ggplot(aes(x = penalty, y = mean, ymin = ymin, ymax = ymax,
#              color = run_transformation)) +
#   theme_bw() +
#   facet_grid(rows = vars(run_imputation),
#              cols = vars(run_exonsInput)) +
#   geom_line() +
#   geom_errorbar(width = .1) +
#   geom_point() +
#   scale_x_log10() +
#   # scale_y_log10(oob = scales::oob_keep) +
#   ylab("Loss Quadratic") + xlab("Penalty (log)") +
#   theme(legend.position = "none")
# 
# res_noperm |>
#   select(-c(results, summary_by_sparsity)) |>
#   unnest(cols = summary_by_penalty) |>
#   mutate(metric = factor(metric, levels = metr_levels)) |>
#   filter(run_transformation != "zscore", run_algo == "glasso") |>
#   mutate(ymin = pmax(0, mean - sd),
#          ymax = mean + sd) |>
#   filter(metric == "loss_frobenius_adj") |>
#   ggplot(aes(x = penalty, y = mean, ymin = ymin, ymax = ymax,
#              color = run_transformation)) +
#   theme_bw() +
#   facet_grid(rows = vars(run_imputation),
#              cols = vars(run_exonsInput)) +
#   geom_line() +
#   geom_errorbar(width = .1) +
#   geom_point() +
#   scale_x_log10() +
#   # scale_y_log10(oob = scales::oob_keep) +
#   ylab("Loss Quadratic") + xlab("Penalty (log)") +
#   theme(legend.position = "none")
# 
# 
# res_noperm |>
#   select(-c(results, summary_by_sparsity)) |>
#   unnest(cols = summary_by_penalty) |>
#   mutate(metric = factor(metric, levels = metr_levels)) |>
#   filter(run_transformation != "zscore", run_algo == "glasso") |>
#   mutate(ymin = pmax(0, mean - sd),
#          ymax = mean + sd) |>
#   filter(metric == "Rsquared") |>
#   ggplot(aes(x = penalty, y = mean, ymin = ymin, ymax = ymax,
#              color = run_transformation)) +
#   theme_bw() +
#   facet_grid(rows = vars(run_imputation),
#              cols = vars(run_exonsInput)) +
#   geom_line() +
#   geom_errorbar(width = .1) +
#   geom_point() +
#   scale_x_log10() +
#   # scale_y_log10(oob = scales::oob_keep) +
#   ylab("R squared") + xlab("Penalty (log)") +
#   theme(legend.position = "none")
# 
# res_noperm |>
#   select(-c(results, summary_by_sparsity)) |>
#   unnest(cols = summary_by_penalty) |>
#   mutate(metric = factor(metric, levels = metr_levels)) |>
#   filter(run_transformation != "zscore", run_algo == "glasso") |>
#   mutate(ymin = pmax(0, mean - sd),
#          ymax = mean + sd) |>
#   filter(metric == "mean_FEV") |>
#   ggplot(aes(x = penalty, y = mean, ymin = ymin, ymax = ymax,
#              color = run_transformation)) +
#   theme_bw() +
#   facet_grid(rows = vars(run_imputation),
#              cols = vars(run_exonsInput)) +
#   geom_line() +
#   geom_errorbar(width = .1) +
#   geom_point() +
#   scale_x_log10() +
#   # scale_y_log10(oob = scales::oob_keep) +
#   ylab("FEV") + xlab("Penalty (log)") +
#   theme(legend.position = "none")
# 
# res_noperm |>
#   select(-c(results, summary_by_sparsity)) |>
#   unnest(cols = summary_by_penalty) |>
#   mutate(metric = factor(metric, levels = metr_levels)) |>
#   filter(run_transformation != "zscore", run_algo == "glasso") |>
#   mutate(ymin = pmax(0, mean - sd),
#          ymax = mean + sd) |>
#   filter(metric == "TPR/FPR") |>
#   ggplot(aes(x = penalty, y = mean, ymin = ymin, ymax = ymax,
#              color = run_transformation)) +
#   theme_bw() +
#   facet_grid(rows = vars(run_imputation),
#              cols = vars(run_exonsInput)) +
#   geom_line() +
#   geom_errorbar(width = .03) +
#   geom_point() +
#   scale_x_log10(limits = c(.04, .9)) +
#   scale_y_continuous(limits = c(0,4), oob = scales::oob_keep) +
#   # scale_y_log10(oob = scales::oob_keep) +
#   ylab("TPR/FPR") + xlab("Penalty (log)") +
#   theme(legend.position = "none")
# 
# res_noperm |>
#   select(-c(results, summary_by_sparsity)) |>
#   unnest(cols = summary_by_penalty) |>
#   mutate(metric = factor(metric, levels = metr_levels)) |>
#   filter(run_transformation != "zscore", run_algo == "glasso") |>
#   mutate(ymin = pmax(0, mean - sd),
#          ymax = mean + sd) |>
#   filter(metric == "power_law") |>
#   ggplot(aes(x = penalty, y = mean, ymin = ymin, ymax = ymax,
#              color = run_transformation)) +
#   theme_bw() +
#   facet_grid(rows = vars(run_imputation),
#              cols = vars(run_exonsInput)) +
#   geom_line() +
#   geom_errorbar(width = .1) +
#   geom_point() +
#   scale_x_log10(limits = c(.04, .9)) +
#   # scale_y_log10(oob = scales::oob_keep) +
#   ylab("Power law") + xlab("Penalty (log)") +
#   theme(legend.position = "none")
# 





# knn k ----


resdirknn <- "data/graph_power4/from_cluster/240430_k/"

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
  filter(metric == "loss_frobenius_adj",
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



#~ plot knn vs k ----

res_knn |>
  select(-c(results, summary_by_sparsity)) |>
  unnest(cols = summary_by_penalty) |>
  mutate(run_k = as.integer(run_k),
         penalty = as.factor(penalty)) |>
  mutate(ymin = pmax(0, mean - sd),
         ymax = mean + sd) |>
  filter(metric %in% more_useful_metrics,
         penalty == .2) |>
  mutate(metric = factor(metric, levels = more_useful_metrics),
         metric = recode_factor(metric,
                    "loss_frobenius_adj" = "Frobenius loss",
                    "mean_FEV" = "Fraction\nExplained Variance",
                    "TPR/FPR" = "TPR / FPR",
                    "power_law" = "Power Law")) |>
  ggplot(aes(x = run_k, y = mean, ymin = ymin, ymax = ymax)) +
  theme_bw() +
  facet_grid(rows = vars(metric), scales = "free_y") +
  geom_line() +
  geom_errorbar(width = .1) +
  geom_point() +
  geom_vline(aes(xintercept = 4), linewidth = 5, alpha = .2) +
  ylab(NULL) + xlab("k") +
  theme(legend.position = "none")

# ggsave("knn.pdf", path = "presentations/240308_figures/",
#        width = 12, height = 17, units = "cm")



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


resdirfin <- "data/graph_power4/from_cluster/240426_perm/"

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
  # filter(run_algo == "glasso", run_exonsInput == "PSI",
  #        run_transformation == "npntr", run_imputation=="knn") |>
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
  
  # ggsave(paste0("pval_", metr, ".png"), plot = gg,
  #        path = "data/intermediates/240426_perm/",
  #        width = 12, height = 7, units = "cm")
  
}

# For fig Suppl: against sparsity, show all permutations ----

design <- "
A#
B#
CD
"

toplot <- resfin |>
  unnest(results) |>
  select(-c(starts_with("run_"), fold)) |>
  pivot_longer(-c(penalty, permutation),
               names_to = "metric",
               values_to = "value") |>
  filter(metric %in% more_useful_metrics) |>
  mutate(permuted = if_else(permutation == 0,
                              "non permuted",
                              "permutation") |>
           factor(levels = c("permutation", "non permuted"))) |>
  mutate(metric = factor(metric, levels = more_useful_metrics),
         metric = recode_factor(metric,
                                "loss_frobenius_adj" = "Frobenius loss",
                                "mean_FEV" = "Fraction\nExplained Variance",
                                "TPR/FPR" = "TPR / FPR",
                                "power_law" = "Power Law"))

ggplot(toplot |> filter(permuted == "permutation")) +
  theme_bw() +
  ggbeeswarm::geom_quasirandom(aes(x = penalty,
                                   y = value,
                                   color = permuted,
                                   # alpha = permuted,
                                   size = permuted,
                                   shape = permuted)) +
  ggbeeswarm::geom_quasirandom(aes(x = penalty,
                                   y = value,
                                   color = permuted,
                                   # alpha = permuted,
                                   size = permuted,
                                   shape = permuted),
                               data = toplot |> filter(permuted == "non permuted")) +
  scale_x_log10() +
  ggh4x::facet_manual(~ metric, scales = "free_y",design = design) +
  # scale_alpha_manual(values = c(0.05,1)) +
  scale_size_manual(values = c(0.5,1.5)) +
  scale_color_manual(values = c("grey80", "darkred")) +
  theme(legend.position = "inside",
        legend.position.inside = c(.8,.8)) +
  ylab(NULL) +
  aes(group=rev(permuted))


# ggsave("presentations/240308_figures/permutations_points.pdf",
#        width = 15, height = 20, units = "cm")








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
  filter(metric %in% more_useful_metrics) |>
  mutate(metric = factor(metric, levels = more_useful_metrics),
         metric = recode_factor(metric,
                                "loss_frobenius_adj" = "Frobenius loss",
                                "mean_FEV" = "Fraction\nExplained Variance",
                                "TPR/FPR" = "TPR / FPR",
                                "power_law" = "Power Law")) |>
  ggplot(aes(x = 100*(1-sparsity_mean),
             y = mean,
             ymin = mean - sd, ymax = mean + sd)) +
  theme_bw() +
  facet_grid(rows = vars(metric),
             scales = "free_y") +
  # ggh4x::facet_manual(~ metric, scales = "free_y",design = design) +
  geom_line() +
  geom_errorbar(width = .01) +
  geom_point(aes(shape = pval_adj < .05,
                 color = pval_adj < .05),
             size = 2) +
  geom_vline(aes(xintercept = 96.97), linewidth = 5, alpha = .2) +
  # scale_x_log10() +
  ylab(NULL) + xlab("Sparsity (%)") +
  scale_color_manual(values = c("grey", 'red4')) +
  theme(legend.position = 'none')

ggsave("presentations/240308_figures/permutations.pdf",
       width = 10, height = 14, units = "cm")





#~ Selected case ----

#~ Check plots ----
final_tib <- read_one_res("240426_final_glasso_PSI_npntrunc_knn_k4_2_16",
                          "data/graph_power4/from_cluster/240426_final")

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




# Check PSI reconstruction ----
all_qs <- qs::qread("data/graph_power4/from_cluster/240415_final/240415_final_glasso_PSI_npntrunc_knn_k10_2_16.qs")

main <- all_qs |>
  filter(penalty == 0.3, permutation == 0)

permuted <- all_qs |>
  filter(penalty == 0.3, permutation == 1)


PSI_measured <- main$se_valid_u[[1]]
PSI_estimated <- main$psi_valid_hat_u[[1]]

plot(PSI_measured, PSI_estimated)


lm(estimated ~ measured,
   data = tibble(measured = as.numeric(PSI_measured),
                 estimated = as.numeric(PSI_estimated))) |>
  summary()



# Compare to permuted 

PSI_measured_perm <- permuted$se_valid_u[[1]]
PSI_estimated_perm <- permuted$psi_valid_hat_u[[1]]

plot(PSI_measured_perm, PSI_estimated_perm)


lm(estimated ~ measured,
   data = tibble(measured = as.numeric(PSI_measured_perm),
                 estimated = as.numeric(PSI_estimated_perm))) |>
  summary()



# plot both
tibble(
  measured = c(
    as.numeric(main$se_valid_u[[1]]),
    as.numeric(permuted$se_valid_u[[1]])
  ),
  estimated = c(
    as.numeric(main$psi_valid_hat_u[[1]]),
    as.numeric(permuted$psi_valid_hat_u[[1]])
  ),
  permuted = rep(c("no", "yes"),
                 times = c(
                   length(as.numeric(main$se_valid_u[[1]])),
                   length(as.numeric(permuted$se_valid_u[[1]]))
                 ))
) |>
  ggplot() +
  theme_classic() +
  geom_point(aes(x = measured, y = estimated), alpha = .1) +
  facet_wrap(~permuted) +
  xlab("PSI measured") + ylab("PSI estimated")



#~~ For small set of events ----
set.seed(123)
selected_se <- sample(colnames(PSI_estimated), 10)

gg_selected_noperm <- tibble(measured = as.numeric(main$se_valid_u[[1]][,selected_se]),
                             estimated = as.numeric(main$psi_valid_hat_u[[1]][,selected_se]),
                             se_id = rep(selected_se, each = nrow(main$se_valid_u[[1]]))) |>
  ggplot() +
  theme_classic() +
  geom_point(aes(x = measured, y = estimated, color = se_id)) +
  xlab("PSI measured") + ylab("PSI estimated") +
  theme(legend.position = "none")


gg_selected_perm <- tibble(measured = permuted$se_valid_u[[1]][,selected_se] |> as.numeric(),
                           estimated = permuted$psi_valid_hat_u[[1]][,selected_se] |> as.numeric(),
                           se_id = rep(selected_se, each = nrow(permuted$se_valid_u[[1]]))) |>
  ggplot() +
  theme_classic() +
  geom_point(aes(x = measured, y = estimated, color = se_id)) +
  xlab("PSI measured") + ylab("PSI estimated (from permutation)")


# patchwork::wrap_plots(gg_selected_noperm, gg_selected_perm)
cowplot::plot_grid(gg_selected_noperm, gg_selected_perm)




#~ Residuals ----

tib_resid <-tibble(se_id = colnames(PSI_estimated) |> rep(each = nrow(PSI_estimated)),
                   sample_id = rownames(PSI_estimated) |> rep(times = ncol(PSI_estimated)),
                   measured = as.numeric(PSI_measured),
                   estimated = as.numeric(PSI_estimated)) |>
  mutate(residual = estimated - measured)

tib_resid |>
  ggplot() +
  theme_classic() +
  geom_point(aes(x = measured, y = residual), alpha = .1)

#~~ break down (ordered by mean resid) ----
tib_resid |>
  mutate(mean_resid = mean(residual, na.rm = TRUE),
         .by = se_id) |>
  arrange(desc(mean_resid)) |>
  mutate(se_id = fct_inorder(se_id)) |>
  ggplot() +
  theme_classic() +
  geom_point(aes(x = se_id, y = residual)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 5))

tib_resid |>
  mutate(mean_resid = mean(residual, na.rm = TRUE),
         .by = sample_id) |>
  arrange(desc(mean_resid)) |>
  mutate(sample_id = fct_inorder(sample_id)) |>
  ggplot() +
  theme_classic() +
  geom_point(aes(x = sample_id, y = residual)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 6))


#~~ ordered by mean abs resid ----
tib_resid |>
  mutate(mean_resid = mean(abs(residual), na.rm = TRUE),
         .by = se_id) |>
  arrange(mean_resid) |>
  mutate(se_id = fct_inorder(se_id)) |>
  ggplot() +
  theme_classic() +
  geom_point(aes(x = se_id, y = residual), alpha = .1) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 5))

#~~ check extremes ----
extr_se_ids <- tib_resid |>
  mutate(mean_resid = mean(abs(residual), na.rm = TRUE),
         .by = se_id) |>
  arrange(mean_resid) |>
  mutate(se_id = fct_inorder(se_id)) |>
  pull(se_id) |> levels() |>
  (\(x) {n <- length(x); c(x[1:3], x[(n-2):n])})()


tib_resid |>
  filter(se_id %in% extr_se_ids) |>
  mutate(se_id = factor(se_id, levels = extr_se_ids)) |>
  ggplot() +
  theme_classic() +
  ggbeeswarm::geom_quasirandom(aes(x = se_id, y = residual,
                                   color = se_id, shape = se_id), alpha = .7) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_color_brewer(type = "qual", palette = 2) +
  scale_shape_manual(values = rep(16:17, each = 3)) +
  theme(legend.position = 'none')

tib_resid |>
  filter(se_id %in% extr_se_ids) |>
  mutate(se_id = factor(se_id, levels = extr_se_ids)) |>
  ggplot() +
  theme_classic() +
  geom_point(aes(x = measured, y = residual,
                 color = se_id, shape = se_id),
             size = 2, alpha = .8) +
  scale_color_brewer(type = "qual", palette = 2) +
  scale_shape_manual(values = rep(16:17, each = 3))

# same by sample
tib_resid |>
  mutate(mean_resid = mean(abs(residual), na.rm = TRUE),
         .by = sample_id) |>
  arrange(mean_resid) |>
  mutate(sample_id = fct_inorder(sample_id)) |>
  ggplot() +
  theme_classic() +
  geom_point(aes(x = sample_id, y = residual)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 6))




#~ Source of error for each metric ----

#~~ Loss Frobenius ----

Sts <- main$S_valid_t[[1]]
Str <- main$S_train_hat_t[[1]]

Sts[1:3,1:3]
Str[1:3,1:3]

printMat::matimage(Sts)

printMat::matimage(permuted$S_valid_t[[1]])

opar <- par()

par(mfrow = c(1,2))
printMat::matimage(main$S_valid_t[[1]], main = "S test")
printMat::matimage(permuted$S_valid_t[[1]], main = "S test (permuted)")
par(opar)


# loss Frobenius on adjacency

is_se <- startsWith(colnames(all_qs$S_valid_t[[1]]), "SE_")
is_sf <- !startsWith(colnames(all_qs$S_valid_t[[1]]), "SE_")



loss_frob_adj <- function(Sts, Str){
  adj_test <- Sts[is_se, is_sf]
  adj_train <- Str[is_se, is_sf]
  
  sqrt(sum( (adj_test - adj_train)^2 ))
}

all_qs$loss_frobenius_adj <- map2_dbl(all_qs$S_valid_t,
                                    all_qs$S_train_hat_t,
                                    loss_frob_adj)

plot(all_qs$penalty, all_qs$loss_frobenius_adj, col = c("red4","green4")[1+all_qs[["permutation"]]])
plot(all_qs$penalty, all_qs$loss_frobenius, col = c("red4","green4")[1+all_qs[["permutation"]]])



#~~ Loss quadratic ----

sum(diag( ( Sts %*na% OMtr - diag(nrow(Sts)) )^2 ))

OMtr <- main$OM_train[[1]]

Sts[1:3,1:3]
OMtr[1:3,1:3]

printMat::matimage(OMtr)
table(OMtr)
OMtr_nodiag <- OMtr
diag(OMtr_nodiag) <- 0
printMat::matimage(OMtr_nodiag)

printMat::matimage(Sts %*na% OMtr)

annot_df <- data.frame(type = if_else(str_detect(colnames(Sts), "^SE_[0-9]{1,4}$"),
                                      "SE",
                                      "SF"),
                       row.names = colnames(Sts))

pheatmap::pheatmap(Sts %*na% OMtr,
                   cluster_rows = FALSE,
                   cluster_cols = FALSE,
                   show_rownames = FALSE,
                   show_colnames = FALSE,
                   annotation_row = annot_df,
                   annotation_col = annot_df)











