# analyze permutation tests from Ruddle

# Initializations ----
#~ Packages ----
library(tidyverse)
library(wbData)


tx2g <- wb_load_tx2gene(281)
gids <- wb_load_gene_ids(281) |>
  add_row(X="0", gene_id = "(Intercept)", symbol ="(Intercept)",
          sequence = "(Intercept)", status="Live",biotype="none",name="(Intercept)")



root_breaks <- function(n = 10, exponent, signed){
  n_default <- n
  function(x, n = n_default){
    if(signed){
      min <- min(abs(x))
      max <- max(abs(x))
    } else{
      min <- min(x)
      max <- max(x)
    }
    
    by <- (max - min)/n
    breaks <- seq(min, 1.1*max^(1/exponent), by = by)^exponent
    breaks <- round(breaks, digits = 2)
    
    if(signed){
      return(c(-breaks, 0, breaks))
    } else{
      return(c(0,breaks))
    }
  }
}

root_trans <- function(exponent = 2, signed = FALSE){
  if(signed){
    scales::trans_new("root_trans",
                      \(x) sign(x)*abs(x)^(1/exponent),
                      \(x) sign(x)*abs(x)^exponent,
                      breaks = root_breaks(exponent = exponent, signed = signed),
                      domain = c(-Inf,Inf))
  } else{
    scales::trans_new("root_trans",
                      \(x) x^(1/exponent),
                      \(x) x^exponent,
                      breaks = root_breaks(exponent = exponent, signed = signed),
                      domain = c(0,Inf))
  }
}





#~ Read general data ----

putative_splice_factors <- wormDatasets::worm_putative_splice_factors |>
  filter(keep == 1 | keep == 2) |>
  pull(gene_id)

sf_expression <- qs::qread("data/intermediates/230117_tx_filtered.qs") |>
  filter(gene_id %in% putative_splice_factors) |>
  mutate(gene_name = i2s(gene_id, gids, warn_missing = TRUE))

sf_tx2g <- sf_expression |>
  select(transcript_id, gene_id, gene_name) |>
  add_row(transcript_id = "(Intercept)",
          gene_id =  "(Intercept)",
          gene_name = "(Intercept)") |>
  distinct()

convert_sf_tx2g <- function(tx_names, warn_missing = TRUE){
  res <- sf_tx2g$gene_id[match(tx_names, sf_tx2g$transcript_id, incomparables = NA)]
  if (warn_missing && any(is.na(res))) {
    warning("tx2g: ", sum(is.na(res)), " tx names could not be converted. NA are returned.")
  }
  res
}


events_coordinates <- read_tsv("data/export_for_arman/221111_events_coordinates.tsv")

convert_event2_gene_id <- function(event_ids, warn_missing = TRUE){
  res <- events_coordinates$gene_id[match(event_ids, events_coordinates$event_id, incomparables = NA)]
  if (warn_missing && any(is.na(res))) {
    warning("converts: ", sum(is.na(res)), " event names could not be converted. NA are returned.")
  }
  res
}



#~ Read simulation results ----

# # With lambda min
# permutations_psi_a <- qs::qread("data/intermediates/230203_regression_permutations_psi_a.qs") |>
#   mutate(sf_id = convert_sf_tx2g(transcript_id),
#          target_id = convert_event2_gene_id(event_id),
#          sf_name = i2s(sf_id, gids, warn_missing = TRUE),
#          target_name = i2s(target_id, gids, warn_missing = TRUE))
# permutations_psi_b <- qs::qread("data/intermediates/230203_regression_permutations_psi_b.qs") |>
#   mutate(sf_id = convert_sf_tx2g(transcript_id),
#          target_id = convert_event2_gene_id(event_id),
#          sf_name = i2s(sf_id, gids, warn_missing = TRUE),
#          target_name = i2s(target_id, gids, warn_missing = TRUE))
# permutations_dpsi_a <- qs::qread("data/intermediates/230203_regression_permutations_dpsi_a.qs") |>
#   mutate(sf_id = convert_sf_tx2g(transcript_id),
#          target_id = convert_event2_gene_id(event_id),
#          sf_name = i2s(sf_id, gids, warn_missing = TRUE),
#          target_name = i2s(target_id, gids, warn_missing = TRUE))
# permutations_dpsi_b <- qs::qread("data/intermediates/230203_regression_permutations_dpsi_b.qs") |>
#   mutate(sf_id = convert_sf_tx2g(transcript_id),
#          target_id = convert_event2_gene_id(event_id),
#          sf_name = i2s(sf_id, gids, warn_missing = TRUE),
#          target_name = i2s(target_id, gids, warn_missing = TRUE))


# with lambda 1se
permutations_psi_a <- qs::qread("data/intermediates/simultation/perm_tests/lambda_1se/230203_regression_permutations_psi_1se_a.qs") |>
  mutate(sf_id = convert_sf_tx2g(transcript_id),
         target_id = convert_event2_gene_id(event_id),
         sf_name = i2s(sf_id, gids, warn_missing = TRUE),
         target_name = i2s(target_id, gids, warn_missing = TRUE))
permutations_psi_b <- qs::qread("data/intermediates/simultation/perm_tests/lambda_1se/230203_regression_permutations_psi_1se_b.qs") |>
  mutate(sf_id = convert_sf_tx2g(transcript_id),
         target_id = convert_event2_gene_id(event_id),
         sf_name = i2s(sf_id, gids, warn_missing = TRUE),
         target_name = i2s(target_id, gids, warn_missing = TRUE))




# Plots ----

permutations_psi_a |>
  ggplot() +
  theme_classic() +
  geom_point(aes(x= coef_effect_size, y = -log10(p_adj)), alpha = .2) +
  # scale_x_continuous(limits = c(-.05,.05)) +
  NULL

# On log scale
permutations_psi_a |>
  ggplot(aes(x= log10(abs(coef_effect_size)), y = -log10(p_adj))) +
  theme_classic() +
  geom_point(aes(color = as.factor(sign(coef_effect_size))))
# ggrepel::geom_text_repel(aes(label = paste0(sf_name,"_",target_name)),
#                          data = perm_res |> filter(p_adj < 1))



permutations_psi_a |>
  filter(p_adj < 1) |>
  ggplot() +
  theme_classic() +
  geom_histogram(aes(x = coef_effect_size), bins = 100, color = 'white') +
  scale_x_continuous(trans = root_trans(4, signed = TRUE))


# Volcano Plot
permutations_psi_a |>#slice_sample(n=100) |> 
  mutate(selected = abs(coef_effect_size) >= 0.005 & p_adj <= 0.1) |>
  ggplot(aes(x= coef_effect_size, y = -log10(p_adj), color = selected)) +
  theme_classic() +
  geom_point(alpha = .2) +
  scale_x_continuous(trans = root_trans(4, signed = TRUE)) +
  xlab(expression(sqrt("effect size", 4))) +
  theme(axis.text = element_text(size = 7)) +
  geom_vline(aes(xintercept = .005)) +
  geom_vline(aes(xintercept = -.005)) +
  geom_hline(aes(yintercept = -log10(.1))) +
  scale_color_manual(values = c("black", "red")) +
  ggtitle("PSI, replicate a")

permutations_psi_b |>#slice_sample(n=100) |> 
  mutate(selected = abs(coef_effect_size) >= 0.005 & p_adj <= 0.1) |>
  ggplot(aes(x= coef_effect_size, y = -log10(p_adj), color = selected)) +
  theme_classic() +
  geom_point(alpha = .2) +
  scale_x_continuous(trans = root_trans(4, signed = TRUE)) +
  xlab(expression(sqrt("effect size", 4))) +
  theme(axis.text = element_text(size = 7)) +
  geom_vline(aes(xintercept = .005)) +
  geom_vline(aes(xintercept = -.005)) +
  geom_hline(aes(yintercept = -log10(.1))) +
  scale_color_manual(values = c("black", "red")) +
  ggtitle("PSI, replicate b")





# Extract significant interactions



interactions_psi_a <- permutations_psi_a |>
  mutate(selected = abs(coef_effect_size) >= 0.005 & p_adj <= 0.1) |>
  filter(selected) |>
  mutate(interaction = paste0(sf_name,"_",target_name)) |>
  pull(interaction)

interactions_psi_b <- permutations_psi_b |>
  mutate(selected = abs(coef_effect_size) >= 0.005 & p_adj <= 0.1) |>
  filter(selected) |>
  mutate(interaction = paste0(sf_name,"_",target_name)) |>
  pull(interaction)




eulerr::euler(list(psi_a = unique(interactions_psi_a),
                   psi_b = unique(interactions_psi_b))) |>
  plot(quantities = TRUE, main = "PSI")




UpSetR::fromList(list(psi_a = unique(interactions_psi_a),
                      psi_b = unique(interactions_psi_b),
                      dpsi_a = unique(interactions_dpsi_a),
                      dpsi_b = unique(interactions_dpsi_b))) |>
  UpSetR::upset()


