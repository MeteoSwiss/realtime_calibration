#' ---
#' title: "A real-time calibration method for the numerical pollen forecast model COSMO-ART"
#' author: "Simon Adamov & Andreas Pauling"
#' date: "`r format(Sys.Date(), '%B %d, %Y')`"
#' always_allow_html: TRUE
#' output:
#'   html_document:
#'     df_print: paged 
#'   pdf_document: default
#'   word_document: default
#' ---
#' 
## ----------------------------------------------------------------------------------------------------------------------------------------------------------
# THE USER CAN SELECT A SPECIES HERE!
# Alnus, Betula, Corylus or Poaceae (Alder, Birch, Hazel or Grasses)

species_sel <- "Alnus"

#' 
#' In this project we want to evaluate the new pollen forecast module that uses "realtime data"
#' for calibration. As new automatic pollen monitors are being deployed, realtime calibration
#' of pollen forecast from weather models becomes more and more relevant.
#' 
#' This study uses old measurements and weather model reforcasts to investigate various
#' implementations and possibilities in terms of forecast improvements.
#' 
#' For very detailed information about the final implementation in COSMO-1E please refer to this 
#' documentation page (ask meteoswiss for access): https://service.meteoswiss.ch/confluence/x/dYQYBQ
#' 
## ----include=FALSE-----------------------------------------------------------------------------------------------------------------------------------------
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

#' 
## ----include=FALSE-----------------------------------------------------------------------------------------------------------------------------------------
# Main colors for the analysis
theme_set(theme_minimal(base_size = 14))
col_names <- c("measurement", "baseline", "calibration")
col_hex <- c("#3b3a3a", "#e98962", "#66C2A5")
names(col_hex) <- col_names

#' 
#' 
## ----include=FALSE-----------------------------------------------------------------------------------------------------------------------------------------
load(paste0(here(), "/data/other/species.RData"))
load(paste0(here(), "/data/other/stations.RData"))
load(paste0(here(), "/data/cosmo/data_cosmo.RData"))
load(paste0(here(), "/data/dwh/data_dwh.RData"))

#' 
#' ## The Data
#' 
#' Three datasets are available for that endeavour:
#' 
#' - _measurement_ referring to daily pollen concentration measurements as collected by a Hirst trap and available in the Data-Warehouse at meteoswiss
#' - _baseline_ referring to daily pollen concentration forecasts from the COSMO-1E model *without* calibraton module
#' - _calibration_ referring to daily pollen concentration forecasts from the COSMO-1E model *with* calibraton module
#' 
#' We have four species available for this analysis:
#' 
#' - Alnus from for the pollen seasons 2020-2021
#' - Betula from for the pollen seasons 2020-2021
#' - Corylus from for the pollen seasons 2020-2021
#' - Poaceae from for the pollen seasons 2019-2020
#' 
#' The analysis will be carried out for each species individually, but all for all stations and years combined.
#' 
#' The observations are provided at 14 different stations, whereof one is excluded from the analysis:
#' Davos is high up in the mountains and pollen measurements are almost always zero.
#' 
#' The following settings are crucial and should always be remembered when running the chunks below:
#' 
#' - The temporal resolution for this analysis is daily averages - pollen verification can be sensitive to temporal resolution
#' - The species of pollen can be assessed individually or combined 
#' - The threshold below which Pollen measurements are excluded from the data is set to 10 currently, as they become unreliant
#' 
#' For all statistical analyses the values <10 are removed from all three timeseries.
#' For some plots these values are not removed. This will be denoted in the comments and
#' should allow the reader to get a full picture of the data.
#' 
## ----include=FALSE-----------------------------------------------------------------------------------------------------------------------------------------
data_combined <- data_dwh %>%
  pivot_wider(names_from = type) %>%
  inner_join(data_cosmo %>%
    pivot_wider(names_from = type), by = c("taxon", "station", "date")) %>%
  select(taxon, station, datetime = datetime.x, date, hour = hour.x, measurement, baseline, calibration) %>%
  pivot_longer(measurement:calibration, names_to = "type", values_to = "value") %>%
  select(taxon, station, type, value, datetime, date, hour) %>%
  filter(taxon == species_sel,
         station != "Davos")

data_above10 <- data_combined %>%
  pivot_wider(names_from = type) %>%
  filter(measurement > 10, baseline > 10, calibration > 10) %>%
  pivot_longer(measurement:calibration, names_to = "type", values_to = "value")


data_log10 <- data_above10 %>%
  mutate(value = log10(value + 0.0001))

data_altman <- data_above10 %>%
  pivot_wider(names_from = type) %>%
  select(datetime, measurement, baseline, calibration) %>%
  mutate(
    mean_baseline = (measurement + baseline) / 2,
    mean_calibration = (measurement + calibration) / 2,
    diff_baseline = baseline - measurement,
    diff_calibration = calibration - measurement
  )

data_corr <- data_log10 %>%
  pivot_wider(names_from = type) %>%
  select(datetime, measurement, baseline, calibration)


data_impact_categories <- data_combined %>%
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
  )


#' 
#' ## Visual Assessment
#' 
#' ### Basic Plots
#' 
#' General overview of the daily concentration values as represented in the three timeseries.
#' First as a boxplot and second as a histogram. 
#' 
## ----echo=FALSE--------------------------------------------------------------------------------------------------------------------------------------------
(gg1 <- data_log10 %>%
  ggplot() +
  geom_boxplot(aes(y = value, fill = type)) +
  theme(
    legend.position = "bottom",
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank()
  ) +
  labs(y = "Log Mean Conc. [Pollen/m³]", x = "") +
  scale_fill_manual(values = col_hex)) +
  ggtitle(paste0("Boxplot of Daily ", species_sel, " Concentrations"))

#' 
## ----echo=FALSE--------------------------------------------------------------------------------------------------------------------------------------------
sd_hirst <- data_combined %>%
  group_by(type) %>%
  summarise(sd = sd(value, na.rm = TRUE))

(gg2 <- data_log10 %>%
  ggplot() +
  geom_histogram(aes(
    y = value,
    fill = type
  ),
  binwidth = 0.1
  ) +
  geom_label(
    data = sd_hirst,
    aes(
      label = paste("Standard Deviation:\n", round(sd), "Pollen / m³"),
      x = 30,
      y = 3,
      group = type
    ),
    size = 3
  ) +
  # xlim(0, 100) +
  facet_wrap(vars(type), ncol = 1) +
  theme(legend.position = "bottom") +
  scale_fill_manual(values = col_hex) +
  coord_flip() +
  labs(x = "Occurence of Pollen Concentrations", y = "Log Mean Conc. [Pollen/m³]")) +
  ggtitle(paste0("Histogram of Daily ", species_sel, " Concentrations"))


#' 
#' ### Correlation Plots
#' 
#' The correlation between the Model and Measurements can be calculated easily and then the CI and p-values must be adjusted for multiple comparison.
#' The corr-test function from the psych handily offers this functionality.
#' 
#' Careful the correlation coefficients method have some serious shortcomings:
#' 
#' The correlation coefficient measures linear agreement--whether the measurements go up-and-down together. 
#' Certainly, we want the measures to go up-and-down together, but the correlation coefficient itself is 
#' deficient in at least three ways as a measure of agreement. (http://www.jerrydallal.com/LHSP/compare.htm)
#' 
#' - The correlation coefficient can be close to 1 (or equal to 1!) even when there is considerable bias between the two methods. For example, if one method gives measurements that are always 10 units higher than the other method, the correlation will be 1 exactly, but the measurements will always be 10 units apart.
#' - The magnitude of the correlation coefficient is affected by the range of subjects/units studied. The correlation coefficient can be made smaller by measuring samples that are similar to each other and larger by measuring samples that are very different from each other. The magnitude of the correlation says nothing about the magnitude of the differences between the paired measurements which, when you get right down to it, is all that really matters.
#' - The usual significance test involving a correlation coefficient-- whether the population value is 0--is irrelevant to the comparability problem. What is important is not merely that the correlation coefficient be different from 0. Rather, it should be close to (ideally, equal to) 1! 
#' 
#' A good summary of the methods and their shortcomings can be found here: https://www.statisticssolutions.com/correlation-Pearson-Kendall-spearman/
#' 
## ----include=FALSE-----------------------------------------------------------------------------------------------------------------------------------------

methods <- c("pearson", "spearman", "kendall")

# For the robust methods (spearman, kendall)
# it doesn't matter whether the transformed data is used or the original

corr_matrix <- map(methods, ~ corr.test(
  data_corr %>% select(-datetime),
  use = "complete",
  method = .x,
  adjust = "holm",
  alpha = .05,
  ci = TRUE,
  minlength = 5
))

#' 
## ----include=FALSE-----------------------------------------------------------------------------------------------------------------------------------------
ci <- map(corr_matrix, ~ .x %>%
  pluck(10)) %>%
  bind_rows() %>%
  round(2) %>%
  mutate(
    method = rep(methods, each = 3),
    metric = rep(c("R-", "rho-", "tau-"), each = 3),
    label = tools::toTitleCase(paste0(metric, method, ": ", lower, " - ", upper)),
    comparison = rep(c("mb", "mc", "bc"), times = 3),
    x = rep(c(0.4, 0.4, 0.41), each = 3),
    y = rep(c(0.01, 0.03, 0.05), each = 3)
  )


#' 
## ----include=FALSE-----------------------------------------------------------------------------------------------------------------------------------------
gg_corr1 <- data_corr %>%
  ggplot(aes(x = log10(measurement), y = log10(baseline))) +
  geom_point(alpha = 0.3, col = "#222225") +
  xlim(0, 0.5) +
  ylim(0, 0.5) +
  geom_smooth(data = data_corr, alpha = 0.2, col = "#3081b8", fill = "#74cbee") +
  geom_abline(slope = 1, intercept = 0, col = "#cc2d2d") +
  geom_label(data = ci %>% filter(comparison == "mb"), aes(label = label, x = x, y = y), parse = TRUE) 

#' 
## ----include=FALSE-----------------------------------------------------------------------------------------------------------------------------------------
gg_corr2 <- data_corr %>%
  ggplot(aes(x = log10(measurement), y = log10(calibration))) +
  geom_point(alpha = 0.3, col = "#222225") +
  xlim(0, 0.5) +
  ylim(0, 0.5) +
  geom_smooth(data = data_corr, alpha = 0.2, col = "#3081b8", fill = "#74cbee") +
  geom_abline(slope = 1, intercept = 0, col = "#cc2d2d") +
  geom_label(data = ci %>% filter(comparison == "mc"), aes(label = label, x = x, y = y), parse = TRUE)

#' 
#' 
## ----echo=FALSE, fig.height = 8, fig.width = 13, fig.dpi=300, out.width="100%"-----------------------------------------------------------------------------
(gg_corr <- ggarrange(gg_corr1, gg_corr2, ncol = 2) %>%
  annotate_figure(
    top = paste0("Pairwise Correlation Between Model and Measurements for ", species_sel, " Pollen"),
    bottom = text_grob("Pairwise correlation between Models; blue line shows the Loess smoother; 
    the red line shows a theroratical perfect correlation of 1. In the text box one can see the 95% 
    confidence intervals of the R-values (adjusted for multiple comparison) as obtained by Pearson and two robust methods.",
      face = "italic", size = 10
    )
  ))

#' 
#' ## Altman Bland Plots
#' 
#' The well established AB-method for clinical trials can be used here as well to compare the means and differences between datasets. 
#' If the points lie within the two SD-line for the differences the datasets can be assumed to be strongly associated with each other. 
#' 
## ----include=FALSE-----------------------------------------------------------------------------------------------------------------------------------------
sd_diff <- data_altman %>%
  select(starts_with("diff")) %>%
  summarise_all(~ sd(.)) %>%
  pivot_longer(1:2, values_to = "sd", names_to = "dummy") %>%
  pull(sd)

gg_ab1 <- data_altman %>%
  ggplot(aes(x = mean_baseline, y = diff_baseline)) +
  geom_point(alpha = 0.2, col = "#222225") +
  xlim(0, 1000) +
  ylim(-sd_diff[2] * 3, sd_diff[2] * 3) +
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

#' 
## ----include=FALSE-----------------------------------------------------------------------------------------------------------------------------------------

gg_ab2 <- data_altman %>%
  ggplot(aes(x = mean_calibration, y = diff_calibration)) +
  geom_point(alpha = 0.2, col = "#222225") +
  xlim(0, 1000) +
  ylim(-sd_diff[2] * 3, sd_diff[2] * 3) +
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

#' 
## ----echo=FALSE, fig.height = 8, fig.width = 13, fig.dpi=300, out.width="100%"-----------------------------------------------------------------------------
(gg_ab <- ggarrange(gg_ab1, gg_ab2, ncol = 2) %>%
  annotate_figure(top = paste0("Altman-Bland Plots for ", species_sel), bottom = text_grob("Pairwise comparison of Models; blue line shows the Loess smother; the red line shows a theroratical perfect
  agreement between Model and Measurements of zero. The dashed red line shows the 2 * sd of the differences, where we expect the points to lie within.",
    face = "italic", size = 12
  )))


#' 
#' ## Density Plots
#' 
#' These plots allow to observe the error for different concentration categories.
#' 
## ----include=FALSE-----------------------------------------------------------------------------------------------------------------------------------------
categs <- c("weak", "medium", "strong", "verystrong")

gg_conc_dens <- list()
gg_conc_box <- list()

labels_y <- list(0.1, 0.015, 0.004, 0.002)
labels_y_hist <- list(15, 13, 10, 10, 7.5, 3)
xlims <- list(50, 500, 1000, 2000)
names(labels_y) <- categs
names(xlims) <- categs

for (j in categs) {
    obs <- data_impact_categories %>%
      filter(categories_measurement == j) %>%
      summarise(n()) %>%
      pull()
    obs <- paste("# of Observations:", obs)

    gg_conc_dens[[j]] <- data_impact_categories %>%
      filter(categories_measurement == j) %>%
      pivot_longer(measurement:calibration, names_to = "type", values_to = "Concentration") %>%
      ggplot() +
      # The area under that whole curve should be 1.
      # To get an estimate of the probability of certain values,
      # you'd have to integrate over an interval on your 'y' axis,
      # and that value should never be greater than 1.
      geom_density(aes(x = Concentration, col = type, fill = type), alpha = 0.15) +
      geom_label(label = obs, aes(x = max(Concentration) * 0.7), y = labels_y[[j]]) +
      scale_colour_manual("", values = col_hex) +
      scale_fill_manual(values = col_hex) +
      coord_cartesian(xlim = c(0, xlims[[j]])) +
      guides(fill = "none") +
      ggtitle(j)

   gg_conc_box[[j]] <- data_impact_categories %>%
      filter(categories_measurement == j) %>%
      pivot_longer(measurement:calibration, names_to = "type", values_to = "Concentration") %>%
      ggplot() +
      # The area under that whole curve should be 1.
      # To get an estimate of the probability of certain values,
      # you'd have to integrate over an interval on your 'y' axis,
      # and that value should never be greater than 1.
      geom_boxplot(aes(x = Concentration, col = type, fill = type), alpha = 0.15) +
      scale_colour_manual("", values = col_hex) +
      scale_fill_manual(values = col_hex) +
      coord_cartesian(xlim = c(0, xlims[[j]])) +
      guides(fill = "none") +
      ggtitle(j)
      
  }

#' 
## ----echo=FALSE, fig.height = 8, fig.width = 13, fig.dpi=300, out.width="100%"-----------------------------------------------------------------------------

(gg_dens_conc <- ggarrange(plotlist = gg_conc_dens) %>%
  annotate_figure(
    top = paste0("Comparison of Measurements and Model Predictions for ", species_sel,
      " for All Stations and Different Concentration Groups."),
    bottom = text_grob(paste0("We are looking at Density Kernel Estimators ",
      "for all three timeseries to compare the measurements between them. ",
      "\n The area under each curve adds up to 1 and makes it possible ",
      "to vizualise the (dis-)similarities of measurements from the ",
      "three timeseries. \n It is basically a smoothed histogram. ",
      "The buckets are based on the mean concentrations of measurements and model."),
      face = "italic",
      size = 10
    )
  ))


#' 
## ----echo=FALSE, fig.height = 8, fig.width = 13, fig.dpi=300, out.width="100%"-----------------------------------------------------------------------------
(gg_dens_box <- ggarrange(plotlist = gg_conc_box, common.legend = TRUE, legend = "bottom") %>%
  annotate_figure(
    top = paste(
      "Comparison of Measurements and Model Predictions for ", species_sel,
      "for All Stations and Different Concentration Groups."
    ),
    bottom = text_grob(paste0("We are looking at Boxplots for different health impact groups to
    investigate the differences between modeled and measured concentrations"),
      face = "italic",
      size = 10
    )
  ))

#' 
#' 
## ----echo=FALSE--------------------------------------------------------------------------------------------------------------------------------------------
data_impact_categories %>%
  mutate(diff_baseline = baseline - measurement,
         diff_calibration = calibration - measurement) %>%
         pivot_longer(diff_baseline:diff_calibration) %>%
         mutate(name = str_remove(name, "diff_")) %>%
         ggplot() +
         geom_boxplot(aes(x = categories_measurement, y = value, fill = name)) +
         scale_fill_manual("", values = col_hex[2:3]) +
        #  scale_y_continuous(trans='log10') + 
         coord_cartesian(ylim = c(-500, 1000)) +
         ggtitle(paste("Differences Between Modelled and Measured", species_sel, "Concentrations")) +
         xlab("") +
         ylab("Concentration [1/m³]") +
         theme(legend.position = "bottom")

#' 
#' ## Statistical Assessment
#' 
#' First, various metrics are compared where the pollen concentrations are considered a continuous numerical variable.
#' 
## ----echo=FALSE--------------------------------------------------------------------------------------------------------------------------------------------
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

my_header <- c(seciesl_sel = 8)
names(my_header) <- species_sel

metrics_baseline %>%
  bind_rows(metrics_calibration) %>%
  kable() %>%
  kable_styling("striped", full_width = FALSE) %>%
  add_header_above(my_header)

#' Second, the values will be converted into health impact based buckets.
#' The impact classes have been defined https://service.meteoswiss.ch/confluence/x/1ZG4
#' Now we can investigate various metrics that are typically used for categoric variables.
#' The Kappa metric is explained here and was chosen as the most meaningful metric for this analysis:
#' https://towardsdatascience.com/multi-class-metrics-made-simple-the-kappa-score-aka-cohens-kappa-coefficient-bdea137af09c
#' 
## ----echo=FALSE--------------------------------------------------------------------------------------------------------------------------------------------
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

my_header <- c(seciesl_sel = 3)
names(my_header) <- species_sel

kappa_baseline %>%
  bind_rows(kappa_calibration) %>%
  pivot_wider(names_from = metric, values_from = ".") %>%
  kable() %>%
  kable_styling("striped", full_width = FALSE) %>%
  add_header_above(my_header)

#' 
#' To takes into account that the health impact levels are ordered. 
#' We can treat them as numerical values from 1:weak - 4: very strong.
#' 
#' 
## ----echo=FALSE--------------------------------------------------------------------------------------------------------------------------------------------

my_header <- c(seciesl_sel = 2)
names(my_header) <- species_sel

data_impact_categories %>%
  select(taxon, station, datetime, categories_measurement, categories_baseline, categories_calibration) %>%
  mutate(across(categories_measurement:categories_calibration, ~as.numeric(.x)),
  error_baseline = categories_baseline - categories_measurement,
  error_calibration = categories_calibration - categories_measurement) %>%
  filter(if_any(categories_measurement:categories_calibration, ~ . > 1)) %>%
  summarize(
  MAE_baseline = mean(abs(error_baseline)),
  MAE_calibration = mean(abs(error_calibration))
  ) %>%
  kable() %>%
  kable_styling("striped", full_width = FALSE) %>%
  add_header_above(my_header)

#' 
#' 
#' 
#' The following table could be used in the appendix.
#' 
#' Reference Event No Event
#' Predicted 
#' Event     A        B
#' No Event  C        D
#' The formulas used here are:
#' 
#' - Sensitivity = A/(A+C)
#' - Specificity = D/(B+D)
#' - Prevalence = (A+C)/(A+B+C+D)
#' - PPV = (sensitivity * prevalence)/((sensitivity*prevalence) + ((1-specificity)*(1-prevalence)))
#' - NPV = (specificity * (1-prevalence))/(((1-sensitivity)*prevalence) + ((specificity)*(1-prevalence)))
#' - Detection Rate = A/(A+B+C+D)
#' - Detection Prevalence = (A+B)/(A+B+C+D)
#' - Balanced Accuracy = (sensitivity+specificity)/2
#' - Precision = A/(A+B)
#' - Recall = A/(A+C)
#' - F1 = (1+beta^2)*precision*recall/((beta^2 * precision)+recall)
#' 
## ----echo=FALSE--------------------------------------------------------------------------------------------------------------------------------------------

my_header <- c(seciesl_sel = 13)
names(my_header) <- species_sel

matrix_baseline$byClass %>%
  as_tibble() %>%
  mutate(
    class = rownames(matrix_baseline$byClass),
    type = "baseline"
  ) %>%
  filter(
    class != "Class: nothing",
    class != "Class: weak"
  ) %>%
  bind_rows(
    matrix_calibration$byClass %>%
      as_tibble() %>%
      mutate(
        class = rownames(matrix_calibration$byClass),
        type = "calibration"
      )
      %>% filter(
        class != "Class: nothing",
        class != "Class: weak"
      )
  ) %>%
  mutate(across(c(-class, -type), ~ round(., digits = 3))) %>%
  kable(escape = FALSE) %>%
  kable_styling("striped", full_width = FALSE) %>%
  add_header_above(my_header)

#' 
#' ## Robust Contrasts with Confidence Intervals
#' 
#' https://www.researchgate.net/publication/282206980_nparcomp_An_R_Software_Package_for_Nonparametric_Multiple_Comparisons_and_Simultaneous_Confidence_Intervals 
#' The R package nparcomp implements a broad range of rank-based nonparametric methods for multiple comparisons. 
#' The single step procedures provide local test decisions in terms of multiplicity adjusted p-values and simultaneous conﬁdence intervals. 
#' The null hypothesis H0: p = 1/2 is significantly rejected at 5% level of significance for many pairwise comparisons.
#' Whenever the p-Value is < than 5% = the confidence interval contains 0.5 -> the effect from the factor trap is not statistically meaningful.
#' The Estimator can also be interpreted as a proxy for the relative difference in median between Model and Measurements.
#' If the Estimator is > 0.5 then the second trap tends to have larger measurements.
#' 
## ----echo=FALSE--------------------------------------------------------------------------------------------------------------------------------------------

npar_contr <- 
  nparcomp(
    value ~ type,
    data = data_log10,
    # slice(1:3870),
    conf.level = 0.95,
    alternative = "two.sided",
    type = "Dunnet",
    control = "measurement"
  )

title <- paste(
  "Robust Contrasts and Confidence Intervals for", species_sel, "Measurements"
)
myheader <- c(title = 8)
names(myheader) <- title

npar_contr$Analysis %>%
  mutate(Taxon = species_sel) %>%
  select(Taxon, Comparison, Estimator, Lower, Upper, pValue = p.Value) %>%
  mutate(pValue = if_else(pValue < 0.05,
    cell_spec(pValue, color = "red"),
    cell_spec(pValue)
  )) %>%
  kable(escape = FALSE) %>%
  kable_styling("striped", full_width = FALSE) %>%
  add_header_above(myheader)

