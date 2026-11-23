# -------------------------------------------------------------------------------------------------
#
#                                    HEYMANS VARIABLE SELECTION
#
# -------------------------------------------------------------------------------------------------

library(tidyverse)
library(mice)
library(rsample)

seed <- 63

df_base <- readRDS("donnees/mimic_iv/derived/bc_candida.rds")

imp <- readRDS("donnees/mimic_iv/derived/bc_candida_imp.rds")

m <- imp$m

# -------------------------------------------------------------------------------------------------
# 1/ Bootstrapper 200 échantillons sur "m" jeux
# -------------------------------------------------------------------------------------------------
n_boot <- 200

boot_un_jeu <- function(df, n_boot, seed) {
  set.seed(seed)
  if ("iep" %in% names(df)) {
    group_bootstraps(df, group = iep, times = n_boot)
  } else {
    bootstraps(df, times = n_boot)
  }
}

boot_mids <- tibble(.imp = seq_len(imp$m)) |>
  mutate(
    data = map(.imp, \(i) complete(imp, action = i)),
    boot = map(data, \(d) boot_un_jeu(d, n_boot = n_boot, seed = seed))
  ) |>
  select(-data)

# Format long : une ligne par (imputation x bootstrap) -> m * n_boot lignes
boot_long <- boot_mids |>
  unnest(boot) |>
  rename(.boot = id)

message(sprintf(
  "[Heymans] %d imputations x %d bootstraps = %d échantillons",
  imp$m,
  n_boot,
  nrow(boot_long)
))

# -------------------------------------------------------------------------------------------------
# 2/ Sélection ascendante stepwise (AIC) sur chaque échantillon bootstrap
# -------------------------------------------------------------------------------------------------
outcome <- "candida"
candidats <- c("surgery", "surgery_abdo", "transfusion", "creat_max", "eer")

stepwise_un_echantillon <- function(split) {
  d <- analysis(split)
  d$y <- as.integer(d[[outcome]] == "Positive")

  mod_vide <- glm(y ~ 1, data = d, family = binomial)
  mod_final <- step(
    mod_vide,
    scope = reformulate(candidats, response = "y"),
    direction = "forward",
    trace = 0
  )

  attr(terms(mod_final), "term.labels")
}

boot_long <- boot_long |>
  mutate(vars_retenues = map(splits, stepwise_un_echantillon))

# Fréquence de sélection de chaque variable sur les m * n_boot échantillons
freq_selection <- boot_long |>
  select(.imp, .boot, vars_retenues) |>
  unnest(vars_retenues) |>
  count(vars_retenues, name = "n") |>
  mutate(freq = n / nrow(boot_long)) |>
  arrange(desc(freq))

freq_selection
