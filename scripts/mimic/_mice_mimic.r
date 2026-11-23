library(tidyverse)
library(mice)

df_base <- readRDS("donnees/mimic_iv/derived/bc_candida.rds")

df_base <- df_base[, 8:14]

df_base$candida[1:15] <- "Positive"

imp <- mice(
  df_base,
  m = 5,
  maxit = 20
)


saveRDS(imp, "donnees/mimic_iv/derived/bc_candida_imp.rds")
