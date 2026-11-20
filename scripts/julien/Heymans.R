# -------------------------------------------------------------------------------------------------
#
#                                    HEYMANS VARIABLE SELECTION
#
# -------------------------------------------------------------------------------------------------

library(tidyverse)
library(mice)
library(rsample)

df_base <- source("donnees/df_clean.rds")
df_mids <- readRDS("donnees/df_imputes.rds")
m <- df_mids$m

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
    boot = map(data, \(d) boot_un_jeu(d, n_boot = n_boot, seed = seed_boot))
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
