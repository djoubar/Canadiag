# =================================================================================================
#
#                                           ROC CURVE
#
# =================================================================================================

library(tidyverse)
library(pROC)
library(yardstick)

dataset_reference <- readRDS("")
model <- readRDS("")
truth_vector <- dataset_reference$df_roc <- data.frame(
  truth = as.factor(truth_vector),
  prediction = predict(model, dataset_reference, type = "response")
) |>
  mutate(prediction_class = as.factor(ifelse(prediction >= 0.5, "1", "0")))

# ------------------
# Figure
# ------------------
fig_roc <- df_roc |>
  roc_curve(truth, prediction, event_level = "second") |>
  ggplot(aes(x = 1 - specificity, y = sensitivity)) +
  geom_line(alpha = 0.7) +
  geom_text(
    data = df_roc |> roc_auc(truth, prediction, event_level = "second"),
    mapping = aes(label = sprintf("AUC = %.2f", .estimate)),
    x = 0.2,
    y = 0.85,
    size = 5,
    inherit.aes = FALSE
  ) +
  geom_path(linewidth = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey60") +
  coord_equal() +
  labs(x = "Non spécificité", y = "Sensibilité") +
  theme_classic() +
  annotation_custom(
    grob = df_roc |>
      conf_mat(prediction_class, truth) |>
      autoplot(type = "heatmap") |>
      ggplotGrob(),
    xmin = 0.65,
    xmax = 1.05,
    ymin = 0.05,
    ymax = 0.45
  )
