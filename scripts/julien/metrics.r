# ================================================================================================= #
#
#                                           METRICS
#
# =================================================================================================
# A paramétrer :
# - dataset_reference
# - model
# - truth_vector = la colonne de la variable à prédire
# - prediction_class = seuil de prédiction à définir

library(tidyverse)
library(yardstick)

dataset_reference <- readRDS("mimic_iv/derived/bc_candida.rds")
model <- readRDS("")
truth_vector <- 2

df_roc <- data.frame(
  truth = as.factor(truth_vector),
  prediction = predict(model, dataset_reference, type = "response")
) |>
  mutate(prediction_class = as.factor(ifelse(prediction >= 0.5, "1", "0")))

# ----
# Metrics
# ----

conf_mat(df_roc, truth = truth, estimate = prediction_class) |>
  summary()
