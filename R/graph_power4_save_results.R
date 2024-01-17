# Convert a tib_quic into results only
# removing the columns with actual data
# Note, the full tib_quic is 25-30GB, can't be loaded interactively

#tib_quic <- qs::qread("data/graph_power4/outputs/231107_quic_50perm_nosep.qs")
#
#tib_quic |>
#  dplyr::select(penalty, fold, permutation, Rsquared, sum_abs_residuals,
#         mean_FEV, loss_frobenius, loss_quadratic,
#         prop_non_zero_coefs_litt, prop_non_zero_coefs_nonlitt) |>
#  readr::write_csv("data/graph_power4/outputs/231107_quic_50perm.csv")


tib_quic <- qs::qread("data/graph_power4/outputs/240116_6penlti_50perm.qs")

tib_quic |>
  dplyr::select(penalty, fold, permutation, Rsquared, sum_abs_residuals,
                mean_FEV, loss_frobenius, loss_quadratic,
                prop_non_zero_coefs_litt, prop_non_zero_coefs_nonlitt) |>
  readr::write_csv("data/graph_power4/outputs/240116_6penlti_50perm.csv")


