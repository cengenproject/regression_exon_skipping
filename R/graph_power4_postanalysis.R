# Process results from running graph_power4 on cluster
# note this is an interactive draft


# From cluster ----






tib_quic <- read_csv("data/intermediates/231109_permutations_cv/231107_quic_50perm_nosep.csv")
tib_quic <- read_csv("data/intermediates/231109_permutations_cv/240116_6penlti_50perm.csv")
tib_quic <- read_csv("data/intermediates/231109_permutations_cv/240118_quic_7penlt_noperm.csv")


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











# Are some PSI easier to predict ----
tib_quic




huge::huge.npn()


# Compare sep and PSI ----

quantifs_filtered |>
  mutate(Nincl = round(PSI * nb_reads),
         Nexcl = round((1-PSI) * nb_reads),
         rPSI = Nincl/(Nincl+Nexcl)) |>
  ggplot() + theme_classic() +
  geom_point(aes(x = PSI, y = rPSI), alpha = .2)




tib_quic <- qs::qread("data/intermediates/231012_cv/231013_quic.qs")

xx <- tib_quic$psi_train[[which(tib_quic$fold == 1 & tib_quic$penalty == 0.5 & tib_quic$permutation == 0)]]

dim(xx)

xx[1:3,1:3]

xx2 <- xx |>
  as.data.frame() |>
  rownames_to_column("sample_id") |>
  pivot_longer(-sample_id,
               names_to = c("event_id", "type"),
               names_sep = "\\.",
               values_to = "count_normalized") |>
  pivot_wider(names_from = type,
              values_from = count_normalized)



# rescale a centered normal into a uniform
# see https://math.stackexchange.com/questions/1063865/transforming-a-normal-distribution-to-a-uniform-one
# and https://math.stackexchange.com/questions/2343952/how-to-transform-gaussiannormal-distribution-to-uniform-distribution
rescale_distr <- function(x){
  2*pnorm(x/sd(x))
}

xx2 |>
  mutate(rPSI = Nincl/(Nincl+Nexcl))


quantifs_filtered |>
  mutate(Nincl = round(PSI * nb_reads),
         Nexcl = round((1-PSI) * nb_reads))

mat_psi[1:3,1:3]
mat_psi_npn[1:3,1:3]

mat_psi_npn <- huge::huge.npn(mat_psi)

# Invert npn ----

xx <- QUIC::QUIC(mat_psi_npn, rho = .1)







mat_psi_unif <- apply(mat_psi_npn, 2, rescale_distr)

# opar <- par(no.readonly = TRUE)

i <- 2
# par(mfrow = c(3,1))
plot(mat_psi[,i], mat_psi_npn[,i])
plot(mat_psi_npn[,i], mat_psi_unif[,i])
plot(mat_psi[,i], mat_psi_unif[,i])


par(opar)




make_npn <- function(x){
  x = qnorm(apply(x, 2, rank)/(nrow(x) + 1))
  x/sd(x[, 1])
}
mat_psi_npn2 <- make_npn(mat_psi)


all.equal(mat_psi_npn2,mat_psi_npn)
xx <- qnorm(apply(x, 2, rank)/(nrow(x) + 1))

x2 <- pnorm(xx)
x2 <- x2 * (nrow(x2)+1)




plot(sort(xx[,1]))
sd(xx[,1])
apply(xx, 1, sd) |> hist(breaks = 100)
x2 <- mat_psi_npn2*sd(mat_psi_npn2[,1])

xx <- apply(x, 2, rank)/(nrow(x) + 1)
x2 <- pnorm(qnorm(xx))
all.equal(x2, xx)




# ~~~~~ -----



xx <- tib_quic |> filter(fold == 1, penalty == .1, permutation == 0)


xx$S_train[[1]][1:3,1:3]


hist(diag(xx$S_train[[1]]))


xx$psi_estimated[[1]][1:3,1:3]

aa <- xx$psi_valid[[1]]


dim(aa)

aa_npn <- huge::huge.npn(aa)
aa[1:3,1:3]


aa_npn_1row <- huge::huge.npn(aa[1,])




