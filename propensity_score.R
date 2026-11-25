# ---------------
#
# ---------------

library(tidyverse)
library(patchwork)
library(MatchIt)
library(cobalt)

dataset_reference <- readRDS("mimic_iv/derived/bc_candida.rds") |>
  select(-candida_org) |>
  na.omit()


model <- glm(
  surgery ~ transfusion + creat_max + surgery_abdo + eer,
  data = dataset_reference,
  family = "binomial"
)

df_propensity <- data.frame(
  ps_score = predict(model, type = "response"),
  surgery = dataset_reference$surgery
)


head(df_propensity)

# --------
#
# --------

ggplot(df_propensity, aes(ps_score, fill = factor(surgery))) +
  geom_histogram(alpha = .5, position = "identity")
# ------
#
# ------

mod_match <- matchit(
  surgery ~ transfusion + creat_max + surgery_abdo + eer,
  data = dataset_reference,
  method = "nearest",
  caliper = 0.002
)

df_matched <- match.data(mod_match)

ggplot(df_matched, aes(ps_score, fill = factor(surgery))) +
  geom_histogram(alpha = .5, position = "identity")

plot(mod_match, type = "hist")
plot(summary(mod_match), abs = TRUE, threshold = 0.1)

love.plot(
  mod_match,
  addl = dataset_reference |> select(-c(surgery, distance)),
  stats = "mean.diffs",
  abs = TRUE,
  binary = "std",
  thresholds = c(m = 0.1),
  var.order = "unadjusted",
  sample.names = c("Avant appariement", "Après appariement"),
  stars = "raw"
)

#
#
#
