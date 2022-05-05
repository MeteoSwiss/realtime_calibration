
library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(magrittr)
library(padr)
library(kableExtra)
library(purrr)
library(ggplot2)
library(ggpubr)
library(here)
library(lubridate)
library(caret)
library(psych)
library(scales)
library(nparcomp)

library(conflicted)
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")

devtools::load_all()



# Main colors for the analysis
theme_set(theme_minimal(base_size = 14))
col_names <- c("measurement", "baseline", "calibration")
col_hex <- c("#3b3a3a", "#e98962", "#66C2A5")
names(col_hex) <- col_names



load(paste0(here(), "/data/other/species.RData"))
load(paste0(here(), "/data/other/stations.RData"))
load(paste0(here(), "/data/cosmo/data_cosmo.RData"))
load(paste0(here(), "/data/dwh/data_dwh.RData"))



data_combined <- data_dwh %>%
  pivot_wider(names_from = type) %>%
  inner_join(data_cosmo %>%
    pivot_wider(names_from = type), by = c("taxon", "station", "date")) %>%
  select(taxon, station, datetime = datetime.x, date, hour = hour.x, measurement, baseline, calibration) %>%
  pivot_longer(measurement:calibration, names_to = "type", values_to = "value") %>%
  select(taxon, station, type, value, datetime, date, hour)

above10_dwh <- data_dwh %>%
  filter(value >= 10) %>%
  select(date, taxon, type, station)

data_above10 <- data_combined %>%
  inner_join(above10_dwh, by = c("date", "taxon", "station")) %>%
  select(taxon, station, type = type.x, value, datetime, date, hour)

data_log10 <- data_above10 %>%
  mutate(value = if_else(value == 0, log10(value + 1), log10(value)))


data_impact_categories <- data_above10 %>%
  pivot_wider(names_from = type) %>%
  mutate(
    categories_measurement = case_when(
      taxon == "Alnus" & measurement < 1 ~ "nothing",
      taxon == "Alnus" & measurement >= 1 & measurement < 10 ~ "weak",
      taxon == "Alnus" & measurement >= 10 & measurement < 70 ~ "medium",
      taxon == "Alnus" & measurement >= 70 & measurement < 250 ~ "strong",
      taxon == "Alnus" & measurement >= 250 ~ "verystrong",
      taxon == "Corylus" & measurement < 1 ~ "nothing",
      taxon == "Corylus" & measurement >= 1 & measurement < 10 ~ "weak",
      taxon == "Corylus" & measurement >= 10 & measurement < 70 ~ "medium",
      taxon == "Corylus" & measurement >= 70 & measurement < 250 ~ "strong",
      taxon == "Corylus" & measurement >= 250 ~ "verystrong",
      taxon == "Betula" & measurement < 1 ~ "nothing",
      taxon == "Betula" & measurement >= 1 & measurement < 10 ~ "weak",
      taxon == "Betula" & measurement >= 10 & measurement < 70 ~ "medium",
      taxon == "Betula" & measurement >= 70 & measurement < 300 ~ "strong",
      taxon == "Betula" & measurement >= 300 ~ "verystrong",
      taxon == "Poaceae" & measurement < 1 ~ "nothing",
      taxon == "Poaceae" & measurement >= 1 & measurement < 20 ~ "weak",
      taxon == "Poaceae" & measurement >= 20 & measurement < 50 ~ "medium",
      taxon == "Poaceae" & measurement >= 50 & measurement < 150 ~ "strong",
      taxon == "Poaceae" & measurement >= 150 ~ "verystrong"
    ),
    categories_baseline = case_when(
      taxon == "Alnus" & baseline < 1 ~ "nothing",
      taxon == "Alnus" & baseline >= 1 & baseline < 10 ~ "weak",
      taxon == "Alnus" & baseline >= 10 & baseline < 70 ~ "medium",
      taxon == "Alnus" & baseline >= 70 & baseline < 250 ~ "strong",
      taxon == "Alnus" & baseline >= 250 ~ "verystrong",
      taxon == "Corylus" & baseline < 1 ~ "nothing",
      taxon == "Corylus" & baseline >= 1 & baseline < 10 ~ "weak",
      taxon == "Corylus" & baseline >= 10 & baseline < 70 ~ "medium",
      taxon == "Corylus" & baseline >= 70 & baseline < 250 ~ "strong",
      taxon == "Corylus" & baseline >= 250 ~ "verystrong",
      taxon == "Betula" & baseline < 1 ~ "nothing",
      taxon == "Betula" & baseline >= 1 & baseline < 10 ~ "weak",
      taxon == "Betula" & baseline >= 10 & baseline < 70 ~ "medium",
      taxon == "Betula" & baseline >= 70 & baseline < 300 ~ "strong",
      taxon == "Betula" & baseline >= 300 ~ "verystrong",
      taxon == "Poaceae" & baseline < 1 ~ "nothing",
      taxon == "Poaceae" & baseline >= 1 & baseline < 20 ~ "weak",
      taxon == "Poaceae" & baseline >= 20 & baseline < 50 ~ "medium",
      taxon == "Poaceae" & baseline >= 50 & baseline < 150 ~ "strong",
      taxon == "Poaceae" & baseline >= 150 ~ "verystrong"
    ),
    categories_calibration = case_when(
      taxon == "Alnus" & calibration < 1 ~ "nothing",
      taxon == "Alnus" & calibration >= 1 & calibration < 10 ~ "weak",
      taxon == "Alnus" & calibration >= 10 & calibration < 70 ~ "medium",
      taxon == "Alnus" & calibration >= 70 & calibration < 250 ~ "strong",
      taxon == "Alnus" & calibration >= 250 ~ "verystrong",
      taxon == "Corylus" & calibration < 1 ~ "nothing",
      taxon == "Corylus" & calibration >= 1 & calibration < 10 ~ "weak",
      taxon == "Corylus" & calibration >= 10 & calibration < 70 ~ "medium",
      taxon == "Corylus" & calibration >= 70 & calibration < 250 ~ "strong",
      taxon == "Corylus" & calibration >= 250 ~ "verystrong",
      taxon == "Betula" & calibration < 1 ~ "nothing",
      taxon == "Betula" & calibration >= 1 & calibration < 10 ~ "weak",
      taxon == "Betula" & calibration >= 10 & calibration < 70 ~ "medium",
      taxon == "Betula" & calibration >= 70 & calibration < 300 ~ "strong",
      taxon == "Betula" & calibration >= 300 ~ "verystrong",
      taxon == "Poaceae" & calibration < 1 ~ "nothing",
      taxon == "Poaceae" & calibration >= 1 & calibration < 20 ~ "weak",
      taxon == "Poaceae" & calibration >= 20 & calibration < 50 ~ "medium",
      taxon == "Poaceae" & calibration >= 50 & calibration < 150 ~ "strong",
      taxon == "Poaceae" & calibration >= 150 ~ "verystrong"
    )
  ) %>%
  mutate_at(
    vars(categories_measurement, categories_baseline, categories_calibration),
    ~ factor(., levels = c("nothing", "weak", "medium", "strong", "verystrong"))
  ) %>%
  filter(categories_measurement != "nothing")


data_impact_categories_mean <- data_above10 %>%
  pivot_wider(names_from = type, values_from = value) %>%
  mutate(mean = (measurement + baseline + calibration) / 3,
    categories_mean = case_when(
      taxon == "Alnus" & mean < 1 ~ "nothing",
      taxon == "Alnus" & mean >= 1 & mean < 10 ~ "weak",
      taxon == "Alnus" & mean >= 10 & mean < 70 ~ "medium",
      taxon == "Alnus" & mean >= 70 & mean < 250 ~ "strong",
      taxon == "Alnus" & mean >= 250 ~ "verystrong",
      taxon == "Corylus" & mean < 1 ~ "nothing",
      taxon == "Corylus" & mean >= 1 & mean < 10 ~ "weak",
      taxon == "Corylus" & mean >= 10 & mean < 70 ~ "medium",
      taxon == "Corylus" & mean >= 70 & mean < 250 ~ "strong",
      taxon == "Corylus" & mean >= 250 ~ "verystrong",
      taxon == "Betula" & mean < 1 ~ "nothing",
      taxon == "Betula" & mean >= 1 & mean < 10 ~ "weak",
      taxon == "Betula" & mean >= 10 & mean < 70 ~ "medium",
      taxon == "Betula" & mean >= 70 & mean < 300 ~ "strong",
      taxon == "Betula" & mean >= 300 ~ "verystrong",
      taxon == "Poaceae" & mean < 1 ~ "nothing",
      taxon == "Poaceae" & mean >= 1 & mean < 20 ~ "weak",
      taxon == "Poaceae" & mean >= 20 & mean < 50 ~ "medium",
      taxon == "Poaceae" & mean >= 50 & mean < 150 ~ "strong",
      taxon == "Poaceae" & mean >= 150 ~ "verystrong"
  )) %>%
  pivot_longer(measurement:calibration, names_to = "type", values_to = "value") %>%
  pivot_wider(names_from = categories_mean, values_from = value)





data_timeseries <-
  data_combined %>%
  filter(
    between(date, as.Date("2020-01-15"), as.Date("2020-04-01")),
    station == "Zürich",
    taxon == "Alnus"
  )

gg1 <- data_timeseries %>%
  ggplot(aes(x = datetime)) +
  geom_line(aes(y = value, col = type), alpha = 0.8) +
  labs(y = "Mean Conc. [Pollen/m³]", x = "") +
  theme(legend.position = "none") +
  scale_color_manual(values = col_hex)



gg2 <- data_timeseries %>%
  ggplot() +
  geom_boxplot(aes(y = log10(value), fill = type)) +
  theme(
    legend.position = "none",
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank()
  ) +
  labs(y = "Log Mean Conc. [Pollen/m³]", x = "") +
  scale_fill_manual(values = col_hex)




sd_hirst <- data_timeseries %>%
  group_by(type) %>%
  summarise(sd = sd(value, na.rm = TRUE))

(gg3 <- data_timeseries %>%
  ggplot() +
  geom_histogram(aes(
    y = log10(value),
    fill = type
  ),
  binwidth = 0.1
  ) +
  geom_label(
    data = sd_hirst,
    aes(
      label = paste("Standard Deviation:\n", round(sd), "Pollen / m³"),
      x = 8,
      y = 2,
      group = type
    ),
    size = 3
  ) +
  facet_wrap(vars(type), ncol = 1) +
  theme(legend.position = "bottom") +
  scale_fill_manual(values = col_hex) +
  coord_flip() +
  labs(x = "Occurence of Pollen Concentrations", y = "Log Mean Conc. [Pollen/m³]"))



## ----echo=FALSE, fig.height = 8, fig.width = 13, fig.dpi=300, out.width="100%"----------------------------------------------------------------------------------------------------------
ggarrange(ggarrange(gg1, gg2, nrow = 2), gg3) %>%
  annotate_figure(top = "Comparison of Models for one Station in one Year")




methods <- c("pearson", "spearman", "kendall")

# For the robust methods (spearman, kendall)
# it doesn't matter whether the transformed data is used or the original

data_corr <- data_log10 %>%
  pivot_wider(names_from = type) %>%
  select(datetime, measurement, baseline, calibration)

corr_matrix <- map(methods, ~ corr.test(
  data_corr %>% select(-datetime),
  use = "complete",
  method = .x,
  adjust = "holm",
  alpha = .05,
  ci = TRUE,
  minlength = 5
))



ci <- map(corr_matrix, ~ .x %>%
  pluck(10)) %>%
  bind_rows() %>%
  round(2) %>%
  mutate(
    method = rep(methods, each = 3),
    metric = rep(c("R-", "rho-", "tau-"), each = 3),
    label = tools::toTitleCase(paste0(metric, method, ": ", lower, " - ", upper)),
    comparison = rep(c("mb", "mc", "bc"), times = 3),
    x = rep(c(0.6, 0.64, 0.55), each = 3),
    y = rep(c(3.2, 3.4, 3.6), each = 3)
  )




gg_corr1 <- data_corr %>%
  ggplot(aes(x = measurement, y = baseline)) +
  geom_point(alpha = 0.1, col = "#222225") +
  geom_smooth(data = data_corr, alpha = 0.3, col = "#3081b8", fill = "#74cbee") +
  geom_abline(slope = 1, intercept = 0, col = "#cc2d2d") +
  geom_label(data = ci %>% filter(comparison == "mb"), aes(label = label, x = x, y = y), parse = TRUE) +
  coord_cartesian(ylim = c(0, 4), xlim = c(0, 4))



gg_corr2 <- data_corr %>%
  ggplot(aes(x = measurement, y = calibration)) +
  geom_point(alpha = 0.1, col = "#222225") +
  geom_smooth(data = data_corr, alpha = 0.3, col = "#3081b8", fill = "#74cbee") +
  geom_abline(slope = 1, intercept = 0, col = "#cc2d2d") +
  geom_label(data = ci %>% filter(comparison == "mc"), aes(label = label, x = x, y = y), parse = TRUE) +
  coord_cartesian(ylim = c(0, 4), xlim = c(0, 4))


## ----echo=FALSE, fig.height = 8, fig.width = 13, fig.dpi=300, out.width="100%"----------------------------------------------------------------------------------------------------------
(gg_corr <- ggarrange(gg_corr1, gg_corr2, ncol = 2) %>%
  annotate_figure(
    top = "Pairwise correlation between Models",
    bottom = text_grob("Pairwise correlation between Models; blue line shows the Loess smoother; the red line shows a theroratical perfect correlation of 1. \n In the text box one can see the 95% confidence intervals of the R-values (adjusted for multiple comparison) as obtained by Pearson and two robust methods.",
      face = "italic", size = 10
    )
  ))




data_altman <- data_above10 %>%
  pivot_wider(names_from = type) %>%
  select(datetime, measurement, baseline, calibration) %>%
  mutate(
    mean_baseline = (measurement + baseline) / 2,
    mean_calibration = (measurement + calibration) / 2,
    diff_baseline = baseline - measurement,
    diff_calibration = calibration - measurement
  )


sd_diff <- data_altman %>%
  select(starts_with("diff")) %>%
  summarise_all(~ sd(.)) %>%
  pivot_longer(1:2, values_to = "sd", names_to = "dummy") %>%
  pull(sd)

gg_ab1 <- data_altman %>%
  ggplot(aes(x = mean_baseline, y = diff_baseline)) +
  geom_point(alpha = 0.1, col = "#222225") +
  coord_cartesian(
    xlim = c(0, 1000),
    ylim = c(
      -sd_diff[1] * 4,
      sd_diff[1] * 4
    )
  ) +
  geom_abline(
    slope = 0, intercept = 0, alpha = 0.8, col = "#cc2d2d"
  ) +
  geom_abline(
    slope = 0, col = "#cc2d2d",
    intercept = sd_diff[1] * 2, alpha = 0.8, linetype = 3
  ) +
  geom_abline(
    slope = 0, col = "#cc2d2d",
    intercept = sd_diff[1] * (-2), alpha = 0.8, linetype = 3
  ) +
  geom_smooth(alpha = 0.3, col = "#3081b8", fill = "#74cbee") +
  labs(y = "Difference(Baseline - Measurement)", x = "Mean(Baseline, Measurement)")




gg_ab2 <- data_altman %>%
  ggplot(aes(x = mean_calibration, y = diff_calibration)) +
  geom_point(alpha = 0.1, col = "#222225") +
  coord_cartesian(
    xlim = c(0, 1000),
    ylim = c(
      -sd_diff[2] * 4,
      sd_diff[2] * 4
    )
  ) +
  geom_abline(
    slope = 0, intercept = 0, alpha = 0.8, col = "#cc2d2d"
  ) +
  geom_abline(
    slope = 0, col = "#cc2d2d",
    intercept = sd_diff[2] * 2, alpha = 0.8, linetype = 3
  ) +
  geom_abline(
    slope = 0, col = "#cc2d2d",
    intercept = sd_diff[2] * (-2), alpha = 0.8, linetype = 3
  ) +
  geom_smooth(alpha = 0.3, col = "#3081b8", fill = "#74cbee") +
  labs(y = "Difference(Calibration - Measurement)", x = "Mean(Calibration, Measurement)")


## ----echo=FALSE, fig.height = 8, fig.width = 13, fig.dpi=300, out.width="100%"----------------------------------------------------------------------------------------------------------
(gg_ab <- ggarrange(gg_ab1, gg_ab2, ncol = 2) %>%
  annotate_figure(top = "Altman-Bland Plots", bottom = text_grob("Pairwise comparison of Models; blue line shows the Loess smother; the red line shows a theroratical perfect
  agreement between Model and Measurements of zero. The dashed red line shows the 2 * sd of the differences, where we expect the points to lie within.",
    face = "italic", size = 12
  )))




categs <- c("weak", "medium", "strong", "verystrong")

gg_conc_dens <- list()

labels_y <- list(0.1, 0.015, 0.004, 0.002)
labels_y_hist <- list(15, 13, 10, 10, 7.5, 3)
names(labels_y) <- categs

for (j in categs) {
  if (j %in% names(data_impact_categories_mean)) {
    obs <- data_impact_categories_mean %>%
      filter(!is.na(!!sym(j))) %>%
      summarise(n() / 2) %>%
      pull()
    obs <- paste("# of Observations:", obs)

    gg_conc_dens[[j]] <- data_impact_categories_mean %>%
      filter(!is.na(!!sym(j))) %>%
      ggplot() +
      # The area under that whole curve should be 1.
      # To get an estimate of the probability of certain values,
      # you'd have to integrate over an interval on your 'y' axis,
      # and that value should never be greater than 1.
      geom_density(aes(x = !!sym(j), col = type, fill = type), alpha = 0.15) +
      geom_label(label = obs, aes(x = max(!!sym(j)) * 0.7), y = labels_y[[j]]) +
      scale_colour_manual("", values = col_hex) +
      scale_fill_manual(values = col_hex) +
      coord_cartesian(xlim = c(0, NA)) +
      guides(fill = "none")
  }
}



## ----echo=FALSE, fig.height = 8, fig.width = 13, fig.dpi=300, out.width="100%"----------------------------------------------------------------------------------------------------------

(gg_dens_conc <- ggarrange(plotlist = gg_conc_dens) %>%
  annotate_figure(
    top = paste(
      "Comparison of Measurements and Model Predictions",
      "for All Stations and Different Concentration Groups."
    ),
    bottom = text_grob(paste0("We are looking at Density Kernel Estimators ",
      "for all three traps to compare the measurements between them. ",
      "\n The area under each curve adds up to 1 and makes it possible ",
      "to vizualise the (dis-)similarities of measurements from the ",
      "three traps. \n It is basically a smoothed histogram. ",
      "The buckets are based on the mean concentrations of measurements and model."),
      face = "italic",
      size = 10
    )
  ))



## ----echo=FALSE-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
metrics_baseline <- data_above10 %>%
  pivot_wider(names_from = type) %>%
  filter(!is.na(baseline)) %>%
  mutate(
    error = baseline - measurement,
  ) %>%
  summarise(
    R2 = cor(baseline, measurement, method = "spearman", use = "complete.obs")^2,
    Bias = mean(error),
    SD = sd(error),
    MAE = mean(abs(error)),
    RMSE = sqrt(mean((error)^2)),
    MSLE = mean((log(1 + baseline) - log(1 + measurement))^2, na.rm = TRUE),
    RMSLE = sqrt(MSLE)
  ) %>%
  mutate(type = "baseline")

metrics_calibration <- data_above10 %>%
  pivot_wider(names_from = type) %>%
  filter(!is.na(calibration)) %>%
  mutate(
    error = calibration - measurement,
  ) %>%
  summarise(
    R2 = cor(calibration, measurement, method = "spearman", use = "complete.obs")^2,
    Bias = mean(error),
    SD = sd(error),
    MAE = mean(abs(error)),
    RMSE = sqrt(mean((error)^2)),
    MSLE = mean((log(1 + calibration) - log(1 + measurement))^2, na.rm = TRUE),
    RMSLE = sqrt(MSLE)
  ) %>%
  mutate(type = "calibration")

metrics_baseline %>%
  bind_rows(metrics_calibration) %>%
  kable() %>%
  kable_styling("striped", full_width = FALSE)


## ----echo=FALSE-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
matrix_baseline <- confusionMatrix(
  data_impact_categories$categories_baseline,
  data_impact_categories$categories_measurement
)
kappa_baseline <- matrix_baseline$overall[1:2] %>%
  tibble() %>%
  mutate(
    type = "baseline",
    metric = c("Accuracy", "Kappa")
  )

matrix_calibration <- confusionMatrix(
  data_impact_categories$categories_calibration,
  data_impact_categories$categories_measurement
)
kappa_calibration <- matrix_calibration$overall[1:2] %>%
  tibble() %>%
  mutate(
    type = "calibration",
    metric = c("Accuracy", "Kappa")
  )

kappa_baseline %>%
  bind_rows(kappa_calibration) %>%
  pivot_wider(names_from = metric, values_from = ".") %>%
  kable() %>%
  kable_styling("striped", full_width = FALSE)


## ----echo=FALSE-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
matrix_baseline$byClass %>%
  as_tibble() %>%
  mutate(
    class = rownames(matrix_baseline$byClass),
    type = "baseline"
  ) %>%
  filter(class != "Class: nothing") %>%
  bind_rows(
    matrix_calibration$byClass %>%
      as_tibble() %>%
      mutate(
        class = rownames(matrix_calibration$byClass),
        type = "calibration"
      ) %>%
      filter(class != "Class: nothing")
  ) %>%
  mutate(across(c(-class, -type), ~ round(., digits = 3))) %>%
  kable() %>%
  kable_styling("striped", full_width = FALSE)


## ----echo=FALSE-------------------------------------------------------------------------------------------------------------------------------------------------------------------------

npar_contr <- map(species$taxon[-3], ~
  nparcomp(
    value ~ type,
    data = data_log10 %>%
      filter(taxon == .x),
    # slice(1:3870),
    conf.level = 0.95,
    alternative = "two.sided",
    type = "Dunnet",
    control = "measurement"
  ))

title <- paste(
  "Robust Contrasts and Confidence Intervals for Corylus Measurements"
)
myheader <- c(title = 6)
names(myheader) <- title

npar_contr[[1]]$Analysis %>%
  bind_rows(
    npar_contr[[2]]$Analysis,
    npar_contr[[2]]$Analysis
  ) %>%
  mutate(taxon = rep(species$taxon[-3], each = 2)) %>%
  mutate(across(where(is.numeric), round, 3)) %>%
  select(Taxon = taxon, Comparison, Estimator, Lower, Upper, pValue = p.Value) %>%
  mutate(pValue = if_else(pValue < 0.05,
    cell_spec(pValue, color = "red"),
    cell_spec(pValue)
  )) %>%
  kable(escape = FALSE) %>%
  kable_styling("striped", full_width = FALSE) %>%
  add_header_above(myheader)

