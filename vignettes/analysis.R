# Copyright (c) 2022 MeteoSwiss, contributors listed in AUTHORS

# Distributed under the terms of the BSD 3-Clause License.

# SPDX-License-Identifier: BSD-3-Clause

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
## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# THE USER CAN SELECT A SPECIES HERE!
# Alnus, Betula, Corylus or Poaceae (Alder, Birch, Hazel or Grasses)

species_sel <- "Alnus"

#'
#' In this project we want to evaluate the new pollen forecast module that uses "realtime data"
#' for calibration. As new automatic pollen monitors are being deployed, realtime calibration
#' of pollen forecast from weather models becomes more and more relevant.
#'
#' This study uses old Hirst measurements and weather model reforcasts to investigate various
#' implementations and possibilities in terms of forecast improvements.
#'
#' For detailed information about the final implementation in COSMO-1E please refer to this
#' documentation page (ask meteoswiss for access): https://service.meteoswiss.ch/confluence/x/dYQYBQ
#'
## ----include=FALSE---------------------------------------------------------------------------------------------------------------------------------------------------------------
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
library(gridExtra)
library(ggtext)
library(here)
library(lubridate)
library(caret)
library(psych)
library(scales)
library(nparcomp)

library(conflicted)
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("date", "lubridate")

devtools::load_all()

#'
## ----include=FALSE---------------------------------------------------------------------------------------------------------------------------------------------------------------
# Main colors for the analysis
theme_set(theme_minimal(base_size = 12))
col_names <- c("measurement", "baseline", "calibration")
col_hex <- c("#3b3a3a", "#f16863", "#58c9dd")
names(col_hex) <- col_names

#'
#'
## ----include=FALSE---------------------------------------------------------------------------------------------------------------------------------------------------------------
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
#' We have four species available for this analysis.
#' The months were selected with our pollen experts and data was retrieved from the model runs.
#'
#' - Alnus from for the pollen seasons 2020-2021 (January - March)
#' - Betula from for the pollen seasons 2020-2021 (March - June)
#' - Corylus from for the pollen seasons 2020-2021 (January - March)
#' - Poaceae from for the pollen seasons 2019-2020 (March - October)
#'
#' The analysis will be carried out for each species individually, but all for all stations and both years combined.
#'
#' The observations are provided at 14 different stations, whereof one is excluded from the analysis:
#' Davos/Wolfgang is high up in the mountains and pollen measurements are almost always zero.
#'
#' The following settings are crucial and should always be remembered when running the chunks below:
#'
#' - The temporal resolution for this analysis is daily averages - pollen verification can be sensitive to temporal resolution
#' - The species of pollen can be assessed individually or combined
#'
#' The measured and modelled values were artificially increased by 0.001 to enable taking the log and division by zero.
#'
## ----include=FALSE---------------------------------------------------------------------------------------------------------------------------------------------------------------
data_combined <- data_dwh %>%
  pivot_wider(names_from = type) %>%
  inner_join(data_cosmo %>%
    pivot_wider(names_from = type), by = c("taxon", "station", "datetime")) %>%
  select(taxon, station, datetime, measurement, baseline, calibration) %>%
  # Remove timesteps with missing data in observations
  filter(!is.na(measurement) & !is.na(baseline) & !is.na(calibration)) %>%
  pivot_longer(measurement:calibration, names_to = "type", values_to = "value") %>%
  select(taxon, station, type, value, datetime) %>%
  filter(
    taxon == species_sel,
    station != "Wolfgang"
  ) %>%
  mutate(
    type = factor(type, ordered = TRUE, levels = c("measurement", "baseline", "calibration")),
    value = value + 0.001
  )

data_log <- data_combined %>%
  mutate(value = log10(value))

data_altman <- data_combined %>%
  pivot_wider(names_from = type) %>%
  select(datetime, measurement, baseline, calibration) %>%
  mutate(
    mean_baseline = (measurement + baseline) / 2,
    mean_calibration = (measurement + calibration) / 2,
    diff_baseline = baseline - measurement,
    diff_calibration = calibration - measurement
  )

data_corr <- data_log %>%
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
#' ## Residual Analysis
#'
#' First we compare the three data sets with the traditional ANOVA approach. Statistical inference (p-values, confidence intervals, . . . )
#' is only valid if the model assumptions are fulfilled. So far, this means (many paragraphs are quoted from Lukas Meier ETH - Script Applied Statistics ANOVA Course):
#'
#' - are the errors independent?
#' - are the errors normally distributed?
#' - is the error variance constant?
#' - do the errors have mean zero?
#'
#' The first assumption is most crucial (but also most difficult to check). If the independence assumption is
#' violated, statistical inference can be very inaccurate. In the ANOVA setting, the last assumption is typically
#' not as important compared to a regression setting, as we are typically fitting “large” models.
#'
## ----echo=FALSE------------------------------------------------------------------------------------------------------------------------------------------------------------------
op <- options(contrasts = c("contr.sum", "contr.poly"))
fit_anova <- aov(as.formula(paste("value ~ type")), data = data_combined)
fit_anova_log <- aov(as.formula(paste("value ~ type")), data = data_log)

#'
#' ### Are the errors independent?
#'
#' In a QQ-plot we plot the empirical quantiles (“what we see in the data”) vs. the theoretical quantiles (“what
#' we expect from the model”). The plot should show a more or less straight line if the distributional assumption
#' is correct. By default, a standard normal distribution is the theoretical “reference distribution”.
#'
#' They are definitely not and we have to do some adjustments. So for the following plot we logarithmic the data to deal with the right-skewedness.
#' The best results were achieved by first logarithmic the data and then taking the square root.
#'
## ----echo=FALSE------------------------------------------------------------------------------------------------------------------------------------------------------------------
gg_res1 <- tibble(residuals = residuals(fit_anova_log, type = "pearson")) %>%
  ggplot(aes(sample = residuals)) +
  stat_qq(col = "#222225", alpha = 0.1) +
  stat_qq_line(col = "#cc2d2d")

gg_res2 <- tibble(residuals = residuals(fit_anova, type = "pearson")) %>%
  ggplot(aes(sample = residuals), col = col_hex[3]) +
  stat_qq(col = "#222225", alpha = 0.1) +
  stat_qq_line(col = "#cc2d2d")

ggarrange(gg_res1, gg_res2, nrow = 1) %>%
  annotate_figure(top = "QQ-Plot for the ANOVA Residuals With (left) and Without Logarithmizing")

#' ### Do the errors have mean zero? & Is the error variance constant?
#' The Tukey-Anscombe plot plots the residuals vs. the fitted values.
#' It allows us to check whether the residuals have constant variance and whether the residuals have mean zero
#' (i.e. they don’t show any deterministic pattern).
#' We don't plot the smoothing line as loess (and other) algorithms have issues when the same value is repeated a large number of times (jitter did not really help).
#'
## ----echo=FALSE------------------------------------------------------------------------------------------------------------------------------------------------------------------
gg_tukey1 <- tibble(
  resid = residuals(fit_anova_log, type = "pearson"),
  fitted = fit_anova_log$fitted.values
) %>%
  ggplot(aes(x = fitted, y = resid)) +
  geom_point(col = "#222225", alpha = 0.5, position = position_jitter(width = 5, height = 0)) +
  # geom_smooth(method = "loess", alpha = 0.2, col = "#3081b8", fill = "#74cbee") +
  geom_abline(slope = 0, intercept = 0, col = "#cc2d2d")
gg_tukey2 <-
  tibble(resid = residuals(fit_anova, type = "pearson"), fitted = fit_anova$fitted.values) %>%
  ggplot(aes(x = fitted, y = resid)) +
  geom_point(col = "#222225", alpha = 0.5, position = position_jitter(width = 5, height = 0)) +
  # geom_smooth(method = "loess", alpha = 0.2, col = "#3081b8", fill = "#74cbee") +
  geom_abline(slope = 0, intercept = 0, col = "#cc2d2d")

ggarrange(gg_tukey1, gg_tukey2) %>%
  annotate_figure(top = "Tukey Anscombe - Plot for the ANOVA Residuals With (left) and Without Logarithmizing")

#'
#' ### Are the errors normally distributed?
#'
#' If the data has some serial structure (i.e., if observations were recorded in a certain time order), we typically
#' want to check whether residuals close in time are more similar than residuals far apart, as this would be a
#' violation of the independence assumption. We can do so by using a so-called index plot where we plot the
#' residuals against time. For positively dependent residuals we would see time periods where most residuals
#' have the same sign, while for negatively dependent residuals, the residuals would “jump” too often from
#' positive to negative compared to independent residuals.
#'
## ----echo=FALSE------------------------------------------------------------------------------------------------------------------------------------------------------------------
resid <- residuals(fit_anova_log, type = "pearson")
resid_df <- tibble(resid = resid, id = as.numeric(names(resid)))

gg_timeline_log <- tibble(id = seq_len(nrow(data_log)), time = data_log$datetime) %>%
  left_join(resid_df, by = "id") %>%
  ggplot(aes(x = time, y = resid)) +
  geom_point(alpha = 0.3, col = "#222225") +
  geom_line(alpha = 0.3, col = "#222225")

resid <- residuals(fit_anova, type = "pearson")
resid_df <- tibble(resid = resid, id = as.numeric(names(resid)))

gg_timeline <- tibble(id = seq_len(nrow(data_combined)), time = data_combined$datetime) %>%
  left_join(resid_df, by = "id") %>%
  ggplot(aes(x = time, y = resid)) +
  geom_point(alpha = 0.3, col = "#222225") +
  geom_line(alpha = 0.3, col = "#222225")

start_date <- case_when(
  species_sel == "Betula" ~ as.POSIXct("2020-03-01"),
  species_sel == "Poaceae" ~ as.POSIXct("2019-03-01"),
  TRUE ~ as.POSIXct("2020-01-01")
)
end_date <- case_when(
  species_sel == "Betula" ~ as.POSIXct("2020-06-30"),
  species_sel == "Poaceae" ~ as.POSIXct("2019-10-31"),
  TRUE ~ as.POSIXct("2020-03-31")
)

gg_timeline_log1 <- gg_timeline_log + coord_cartesian(x = c(start_date, end_date))
gg_timeline_log2 <- gg_timeline_log + coord_cartesian(x = c(start_date + years(1), end_date + years(1)))
gg_timeline1 <- gg_timeline + coord_cartesian(x = c(start_date, end_date))
gg_timeline2 <- gg_timeline + coord_cartesian(x = c(start_date + years(1), end_date + years(1)))

ggarrange(gg_timeline_log1, gg_timeline1, gg_timeline_log2, gg_timeline2, nrow = 2, ncol = 2) %>%
  annotate_figure(top = "Index-Plot for the ANOVA Residuals With (left) and Without Logarithmizing")

#'
#' Summary: Redidual Analysis shows that the assumptions of "normal" statiscal methods are validated even for log(daily values).
#' It is therefore suggested to continue the analysis with robust and simple metrics.
#'
#' ## Visual Assessment
#'
#' ### Basic Plots
#'
#' General overview of the daily concentration values as represented in the three timeseries.
#' Each plot shows one species in one year. The seperate lines represent individual measurement stations.
#'
## ----echo=FALSE, warning=FALSE---------------------------------------------------------------------------------------------------------------------------------------------------
startdate <- case_when(
  species_sel == "Betula" ~ ymd_hms("2021-03-01 00:00:00"),
  species_sel == "Poaceae" ~ ymd_hms("2020-03-01 00:00:00"),
  TRUE ~ ymd_hms("2021-01-01 00:00:00")
)
enddate <- case_when(
  species_sel == "Betula" ~ ymd_hms("2021-06-30 00:00:00"),
  species_sel == "Poaceae" ~ ymd_hms("2020-09-30 00:00:00"),
  TRUE ~ ymd_hms("2021-03-31 00:00:00")
)

data_otheryear <- data_combined %>%
  filter(between(datetime, startdate - years(1), enddate - years(1))) %>%
  mutate(
    datetime = datetime + years(1),
    station = paste0("other", station)
  )

(gg_timeseries <- data_combined %>%
  filter(
    taxon == species_sel,
    between(datetime, startdate, enddate)
  ) %>%
  bind_rows(data_otheryear) %>%
  ggplot() +
  geom_line(aes(
    x = as.Date(datetime),
    y = log10(value + 1),
    col = type,
    fill = station
  ), alpha = 0.25) +
  facet_wrap(~type, nrow = 3) +
  theme(legend.position = "none") +
  scale_color_manual("", values = col_hex) +
  scale_x_date(date_labels = "%d. %b") +
  xlab("") +
  ylab(paste0("Daily Mean Log. Concentrations [Pollen / m³]")) +
  ggtitle(paste0("Time Series of Daily ", "*", species_sel, "*", " Pollen Concentrations")) +
  theme(plot.title = element_markdown()))

#'
#' First as a boxplot and second as a histogram of the daily differences.
#'
## ----echo=FALSE------------------------------------------------------------------------------------------------------------------------------------------------------------------
(gg_boxplot <- data_combined %>%
  pivot_wider(names_from = type, values_from = value) %>%
  mutate(
    baseline = baseline - measurement,
    calibration = calibration - measurement
  ) %>%
  select(-measurement) %>%
  pivot_longer(baseline:calibration, names_to = "type") %>%
  ggplot() +
  geom_boxplot(aes(y = value, fill = type), alpha = 0.9) +
  theme(
    legend.position = "bottom",
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank()
  ) +
  labs(y = paste0("Modelled - Measured Concentrations [Pollen / m³]", x = "")) +
  coord_cartesian(y = c(-100, 100)) +
  scale_fill_manual("", values = col_hex[2:3]) +
  ggtitle(paste0(
    "Boxplot of Daily ", "*",
    species_sel, "*", " Pollen Concentration Differences"
  )) +
  theme(plot.title = ggtext::element_markdown()))

#'
## ----echo=FALSE------------------------------------------------------------------------------------------------------------------------------------------------------------------
sd_hirst <- data_combined %>%
  group_by(type) %>%
  summarise(sd = sd(value))

(gg_hist <- data_log %>%
  ggplot() +
  geom_histogram(
    aes(
      y = value,
      fill = type
    ),
    binwidth = 0.1
  ) +
  geom_label(
    data = sd_hirst,
    aes(
      label = paste("Standard Deviation:\n", round(sd), "Pollen / m³"),
      x = 400,
      y = 3,
      group = type
    ),
    size = 3
  ) +
  facet_wrap(vars(type), ncol = 1) +
  theme(legend.position = "bottom") +
  scale_fill_manual(values = col_hex) +
  coord_flip() +
  labs(x = "Occurence of Pollen Concentrations", y = "Log Mean Conc. [Pollen / m³]") +
  ggtitle(paste0("Histogram of Daily ", "*", species_sel, "*", " Pollen Concentrations")) +
  theme(plot.title = ggtext::element_markdown()))

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
## ----include=FALSE---------------------------------------------------------------------------------------------------------------------------------------------------------------

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
## ----include=FALSE---------------------------------------------------------------------------------------------------------------------------------------------------------------
ci <- map(corr_matrix, ~ .x %>%
  pluck(10)) %>%
  bind_rows() %>%
  round(2) %>%
  mutate(
    method = tools::toTitleCase(rep(methods, each = 3)),
    metric = rep(c("R", "Rho", "Tau"), each = 3),
    label = paste0(method, "~", metric, ": ", lower, " - ", upper),
    comparison = rep(c("mb", "mc", "bc"), times = 3),
    x = rep(c(2.3, 2.35, 2.3), each = 3),
    y = rep(c(0, 0.19, 0.38), each = 3)
  )

#'
## ----include=FALSE---------------------------------------------------------------------------------------------------------------------------------------------------------------
gg_corr1 <- data_corr %>%
  ggplot(aes(x = measurement, y = baseline)) +
  geom_point(alpha = 0.3, col = "#222225") +
  coord_cartesian(x = c(0, 3), y = c(0, 3.5)) +
  geom_smooth(data = data_corr, alpha = 0.2, col = "#3081b8", fill = "#74cbee") +
  geom_abline(slope = 1, intercept = 0, col = "#cc2d2d") +
  geom_label(
    data = ci %>% filter(comparison == "mb"),
    aes(label = label, x = x, y = y),
    parse = TRUE
  ) +
  xlab("log10(measurement)") +
  ylab("log10(baseline)")

#'
## ----include=FALSE---------------------------------------------------------------------------------------------------------------------------------------------------------------
gg_corr2 <- data_corr %>%
  ggplot(aes(x = measurement, y = calibration)) +
  geom_point(alpha = 0.3, col = "#222225") +
  coord_cartesian(x = c(0, 3), y = c(0, 3.5)) +
  geom_smooth(data = data_corr, alpha = 0.2, col = "#3081b8", fill = "#74cbee") +
  geom_abline(slope = 1, intercept = 0, col = "#cc2d2d") +
  geom_label(
    data = ci %>% filter(comparison == "mc"),
    aes(label = label, x = x, y = y),
    parse = TRUE
  ) +
  xlab("log10(measurement)") +
  ylab("log10(calibration)")

#'
#'
## ----echo=FALSE, warning=FALSE, message=FALSE, fig.height = 8, fig.width = 13, fig.dpi=300, out.width="100%"---------------------------------------------------------------------
(gg_corr <- ggarrange(gg_corr1, gg_corr2, ncol = 2) %>%
  annotate_figure(
    top = arrangeGrob(
      text_grob("Pairwise Correlation Between Model and Measurements for ",
        vjust = 1.2, hjust = 0.15, size = 14
      ),
      text_grob(species_sel, vjust = 1.2, hjust = -3, face = "italic", size = 14),
      text_grob(" Pollen", vjust = 1.2, hjust = 1.85, size = 14),
      nrow = 1, padding = unit(4, "line")
    )
    # bottom = text_grob("Pairwise correlation between Models; blue line shows the Loess smoother;
    # the red line shows a theoretical perfect correlation of 1. In the text box one can see the 95%
    # confidence intervals of the R-values (adjusted for multiple comparison)
    # as obtained by Pearson and two robust methods.",
    # face = "italic", size = 10)
  )
)

#'
#' ## Altman Bland Plots
#'
#' The well established AB-method for clinical trials can be used here as well to compare the means and differences between datasets.
#' If the points lie within the two SD-line for the differences the datasets can be assumed to be strongly associated with each other.
#'
## ----include=FALSE---------------------------------------------------------------------------------------------------------------------------------------------------------------
sd_diff <- data_altman %>%
  select(starts_with("diff")) %>%
  summarise_all(~ sd(.)) %>%
  pivot_longer(1:2, values_to = "sd", names_to = "dummy") %>%
  pull(sd)

gg_ab1 <- data_altman %>%
  ggplot(aes(x = mean_baseline, y = diff_baseline)) +
  geom_point(alpha = 0.2, col = "#222225") +
  coord_cartesian(x = c(0, 300), y = c(-sd_diff[2] * 5, sd_diff[2] * 5)) +
  geom_abline(
    slope = 0, intercept = 0, alpha = 0.8, col = "#cc2d2d"
  ) +
  geom_abline(
    slope = 0, col = "#cc2d2d",
    intercept = sd_diff[1] * 1.96, alpha = 0.8, linetype = 3
  ) +
  geom_abline(
    slope = 0, col = "#cc2d2d",
    intercept = sd_diff[1] * (-1.96), alpha = 0.8, linetype = 3
  ) +
  geom_smooth(alpha = 0.3, col = "#3081b8", fill = "#74cbee") +
  labs(y = "Difference(Baseline - Measurement)", x = "Mean(Baseline, Measurement)") +
  theme(plot.margin = margin(0.6, 0.2, 0.2, 0.2, "cm"))

#'
## ----include=FALSE---------------------------------------------------------------------------------------------------------------------------------------------------------------
gg_ab2 <- data_altman %>%
  ggplot(aes(x = mean_calibration, y = diff_calibration)) +
  coord_cartesian(x = c(0, 300), y = c(-sd_diff[2] * 5, sd_diff[2] * 5)) +
  geom_point(alpha = 0.2, col = "#222225") +
  geom_abline(
    slope = 0, intercept = 0, alpha = 0.8, col = "#cc2d2d"
  ) +
  geom_abline(
    slope = 0, col = "#cc2d2d",
    intercept = sd_diff[2] * 1.96, alpha = 0.8, linetype = 3, size = 1
  ) +
  geom_abline(
    slope = 0, col = "#cc2d2d",
    intercept = sd_diff[2] * (-1.96), alpha = 0.8, linetype = 3, size = 1
  ) +
  geom_smooth(alpha = 0.3, col = "#3081b8", fill = "#74cbee") +
  labs(y = "Difference(Calibration - Measurement)", x = "Mean(Calibration, Measurement)") +
  theme(plot.margin = margin(0.6, 0.2, 0.2, 0.2, "cm"))

#'
## ----echo=FALSE, warning=FALSE, message=FALSE, fig.height = 8, fig.width = 13, fig.dpi=300, out.width="100%"---------------------------------------------------------------------
(gg_ab <- ggarrange(gg_ab1, gg_ab2, ncol = 2) %>%
  annotate_figure(
    top = gridExtra::arrangeGrob(
      text_grob("Bland-Altman Plots for ", vjust = 1.2, hjust = -0.7, size = 14),
      text_grob(species_sel, vjust = 1.2, hjust = -0.83, face = "italic", size = 14),
      text_grob(" Pollen", vjust = 1.2, hjust = 3.7, size = 14),
      nrow = 1, padding = unit(4, "line")
    )
    # bottom = text_grob("Pairwise comparison of Models; blue line shows the Loess smother; the red line shows a theoretical perfect
    # agreement between Model and Measurements of zero. The dashed red line shows the 2 * sd of the differences, where we expect the points to lie within.",
    #   face = "italic", size = 12)
  ))

#'
#' ## Density Plots
#'
#' These plots allow to observe the error for different concentration categories.
#'
## ----include=FALSE---------------------------------------------------------------------------------------------------------------------------------------------------------------
categs <- c("weak", "medium", "strong", "verystrong")

gg_conc_dens <- list()
gg_conc_box <- list()

labels_y <- list(0.2, 0.02, 0.0075, 0.002)
labels_y_hist <- list(15, 13, 10, 10, 7.5, 3)
xlims <- case_when(
  species_sel == "Corylus" ~ list(50, 500, 750, 1000),
  species_sel == "Poaceae" ~ list(50, 200, 400, 750),
  TRUE ~ list(50, 500, 1000, 2000)
)
xlim_box <- 100

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
    annotate("text", x = xlims[[j]] * 0.6, y = labels_y[[j]], label = obs) +
    scale_colour_manual("", values = col_hex) +
    scale_fill_manual(values = col_hex) +
    coord_cartesian(xlim = c(0, xlims[[j]])) +
    guides(fill = "none") +
    ggtitle(j)

  gg_conc_box[[j]] <- data_impact_categories %>%
    mutate(
      baseline = baseline - measurement,
      calibration = calibration - measurement
    ) %>%
    select(-measurement) %>%
    filter(categories_measurement == j) %>%
    pivot_longer(baseline:calibration, names_to = "type", values_to = "Concentration") %>%
    mutate(type = factor(type, levels = c("calibration", "baseline", "measurement"))) %>%
    ggplot() +
    # The area under that whole curve should be 1.
    # To get an estimate of the probability of certain values,
    # you'd have to integrate over an interval on your 'y' axis,
    # and that value should never be greater than 1.
    geom_boxplot(aes(x = Concentration, col = type, fill = type), alpha = 0.15) +
    scale_colour_manual("", values = col_hex[2:3]) +
    scale_fill_manual(values = col_hex[2:3]) +
    annotate("text", x = xlim_box * 0.54, y = 0.4, label = obs, size = 3.5, alpha = 0.7) +
    coord_cartesian(xlim = c(-xlim_box, xlim_box)) +
    guides(fill = "none") +
    xlab("Concentration [Pollen / m³]") +
    ylab("") +
    ggtitle(j) +
    theme(
      axis.text.y = element_blank(),
      axis.text.y.left = element_blank(),
      plot.margin = margin(1, 0, 0, 0, "cm")
    )
  xlim_box <- xlim_box + 150
}

#'
## ----echo=FALSE, fig.height = 8, fig.width = 13, fig.dpi=300, out.width="100%"---------------------------------------------------------------------------------------------------
(gg_dens_conc <- ggarrange(plotlist = gg_conc_dens) %>%
  annotate_figure(
    top = paste0(
      "Comparison of Measurements and Model Predictions for ", species_sel,
      " Pollen for All Stations and Different Concentration Groups."
    ),
    bottom = text_grob(
      paste0(
        "We are looking at Density Kernel Estimators ",
        "for all three timeseries to compare the measurements between them. ",
        "\n The area under each curve adds up to 1 and makes it possible ",
        "to vizualise the (dis-)similarities of measurements from the ",
        "three timeseries. \n It is basically a smoothed histogram. ",
        "The buckets are based on the mean concentrations of measurements and model."
      ),
      face = "italic",
      size = 10
    )
  ))

#'
## ----echo=FALSE, fig.height = 8, fig.width = 13, fig.dpi=300, out.width="100%"---------------------------------------------------------------------------------------------------
(gg_boxplot_conc <- ggarrange(
  plotlist = gg_conc_box,
  common.legend = TRUE,
  legend = "right",
  vjust = 5
) %>%
  annotate_figure(
    top = gridExtra::arrangeGrob(
      text_grob("Model vs. Measurement Boxplot Comparison for Different Health Impact Groups for ",
        vjust = 1.2, hjust = 0.18, size = 12
      ),
      text_grob(species_sel, vjust = 1.2, hjust = -6, face = "italic", size = 12),
      text_grob(" Pollen", vjust = 1.2, hjust = -0.4, size = 12),
      nrow = 1, padding = unit(4, "line")
    )
    # bottom = text_grob(paste0("We are looking at Boxplots for different health impact groups to
    # investigate the differences between modeled and measured concentrations"),
    #   face = "italic",
    #   size = 10)
  ))

#'
#'
## ----echo=FALSE------------------------------------------------------------------------------------------------------------------------------------------------------------------
data_impact_categories %>%
  mutate(
    diff_baseline = baseline - measurement,
    diff_calibration = calibration - measurement
  ) %>%
  pivot_longer(diff_baseline:diff_calibration) %>%
  mutate(name = str_remove(name, "diff_")) %>%
  ggplot() +
  geom_boxplot(aes(x = categories_measurement, y = value, fill = name)) +
  scale_fill_manual("", values = col_hex[2:3]) +
  ggtitle(paste(
    "Differences Between Modelled and Measured",
    species_sel,
    " Pollen Concentrations"
  )) +
  xlab("") +
  ylab("Concentration [Pollen / m³]") +
  theme(legend.position = "bottom")

#'
#' ## Statistical Assessment
#'
#' First, various metrics are compared where the pollen concentrations are considered a continuous numerical variable.
#'
## ----echo=FALSE------------------------------------------------------------------------------------------------------------------------------------------------------------------
metrics_baseline <- data_combined %>%
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
    MSLE = mean((log10(1 + baseline) - log10(1 + measurement))^2, na.rm = TRUE),
    RMSLE = sqrt(MSLE)
  ) %>%
  mutate(type = "baseline")

metrics_calibration <- data_combined %>%
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
    MSLE = mean((log10(1 + calibration) - log10(1 + measurement))^2, na.rm = TRUE),
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
## ----echo=FALSE------------------------------------------------------------------------------------------------------------------------------------------------------------------
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
#' We can treat them as numerical values from 0:nothing - 4: very strong.
#'
#'
## ----echo=FALSE------------------------------------------------------------------------------------------------------------------------------------------------------------------
my_header <- c(seciesl_sel = 2)
names(my_header) <- species_sel

data_impact_categories %>%
  select(
    taxon,
    station,
    datetime,
    categories_measurement,
    categories_baseline,
    categories_calibration
  ) %>%
  mutate(across(categories_measurement:categories_calibration, ~ as.numeric(.x)),
    error_baseline = categories_baseline - categories_measurement,
    error_calibration = categories_calibration - categories_measurement
  ) %>%
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
## ----echo=FALSE------------------------------------------------------------------------------------------------------------------------------------------------------------------
my_header <- c(seciesl_sel = 13)
names(my_header) <- species_sel

matrix_baseline$byClass %>%
  as_tibble() %>%
  mutate(
    class = factor(rownames(matrix_baseline$byClass),
      ordered = TRUE,
      levels = c(
        "Class: nothing",
        "Class: weak",
        "Class: medium",
        "Class: strong",
        "Class: verystrong"
      )
    ),
    type = "baseline"
  ) %>%
  bind_rows(
    matrix_calibration$byClass %>%
      as_tibble() %>%
      mutate(
        class = factor(rownames(matrix_baseline$byClass),
          ordered = TRUE,
          levels = c(
            "Class: nothing",
            "Class: weak",
            "Class: medium",
            "Class: strong",
            "Class: verystrong"
          )
        ),
        type = "calibration"
      )
  ) %>%
  mutate(across(c(-class, -type), ~ round(., digits = 3))) %>%
  arrange(class) %>%
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
#' If the Estimator is > 0.5 then the model tends to have larger measurements.
#'
## ----echo=FALSE------------------------------------------------------------------------------------------------------------------------------------------------------------------
npar_contr <-
  nparcomp_adjusted(
    value ~ type,
    data = data_combined,
    conf.level = 0.95,
    alternative = "two.sided",
    type = "Dunnet",
    control = "measurement"
  )

title <- paste(
  "Robust Contrasts and Confidence Intervals for", species_sel, " Pollen Measurements"
)
myheader <- c(title = 6)
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

#'
## ----paper_plots, include=FALSE--------------------------------------------------------------------------------------------------------------------------------------------------
# Export plots for paper
ggsave("vignettes/gg_timeseries.png", gg_timeseries, width = 9, height = 6, bg = "white", dpi = 300)
ggsave("vignettes/gg_boxplot.png", gg_boxplot, width = 9, height = 7, bg = "white", dpi = 300)
ggsave("vignettes/gg_corr.png", gg_corr, width = 10, height = 6, bg = "white", dpi = 300)
ggsave("vignettes/gg_ab.png", gg_ab, width = 10, height = 6, bg = "white", dpi = 300)
ggsave("vignettes/gg_boxplot_conc.png",
  gg_boxplot_conc,
  width = 9,
  height = 7,
  bg = "white",
  dpi = 300
)

#'
#' Reviewer 1 requested an additional graph displaying the changes of the tune factor over the course
#' of a season. The following chunk loads the data and saves a line-plot:
#'
#' Final variant of the model:
#'
#' - 5 days were regarded for the threshold, the sum of all hourly values had to surpass 720.
#' - In addition the last 24h were regarded for the threshold, the sum of all hourly values had to surpass 240.
#' - The full update to the tuning factor happens once per day -> 24th root
#' - Tuning factor was limited based on climatologically observed minima and maxima
#'
## ----echo=FALSE, message=FALSE---------------------------------------------------------------------------------------------------------------------------------------------------

if (species_sel == "Alnus") {
  tune_20 <- import_data_cosmo(paste0(here(), "/ext-data/cosmo/20_alnu_tune"), type = "tune")
  tune_21 <- import_data_cosmo(paste0(here(), "/ext-data/cosmo/21_alnu_tune"), type = "tune")

  gg_tune <- tune_20 %>%
    aggregate_pollen() %>%
    mutate(
      datetime = datetime + years(1),
      value = value / 0.34,
      season = "season_2020"
    ) %>%
    bind_rows(tune_21 %>%
      aggregate_pollen() %>%
      mutate(season = "season_2021")) %>%
    filter(
      taxon == species_sel,
      station != "Wolfgang",
      between(datetime, startdate, enddate),
    ) %>%
    ggplot() +
    geom_line(aes(
      x = as.Date(datetime),
      y = value,
      col = season
    ), alpha = 0.6) +
    facet_wrap(~station) +
    coord_cartesian(y = c(0, 7)) +
    theme(
      legend.title = element_blank(), legend.position = "bottom",
      axis.text = element_text(size = 7)
    ) +
    scale_x_date(date_labels = "%e.%b") +
    xlab("") +
    ylab(paste0("")) +
    ggtitle(paste0("Daily Mean Tuning Factor Values (ALNUtune) - Final Model Variant"))

  ggsave("vignettes/gg_tune.png", gg_tune, width = 9, height = 6, bg = "white", dpi = 300)

  gg_tune
}

#'
#' Model Variant 1 below had the following settings:
#'
#' - Only last 5 days were regarded for the threshold, the sum of all hourly values had to surpass 720.
#' - The update to the tuning factor happened slowly (conservatively) -> 120th root
#' - Tuning factor was unlimited
#'
## ----echo=FALSE, message=FALSE---------------------------------------------------------------------------------------------------------------------------------------------------
if (species_sel == "Alnus") {
  tune_20_v1 <- import_data_cosmo(paste0(here(), "/ext-data/cosmo/20_alnu_tune_v1"), type = "tune")
  tune_21_v1 <- import_data_cosmo(paste0(here(), "/ext-data/cosmo/21_alnu_tune_v1"), type = "tune")

  gg_tune_v1 <- tune_20_v1 %>%
    aggregate_pollen() %>%
    mutate(
      datetime = datetime + years(1),
      season = "season_2020"
    ) %>%
    bind_rows(tune_21_v1 %>%
      aggregate_pollen() %>%
      mutate(season = "season_2021")) %>%
    mutate(value = value / 0.34) %>%
    filter(
      taxon == species_sel,
      station != "Wolfgang",
      between(datetime, startdate, enddate),
    ) %>%
    ggplot() +
    geom_line(aes(
      x = as.Date(datetime),
      y = value,
      col = season
    ), alpha = 0.6) +
    facet_wrap(~station) +
    coord_cartesian(y = c(0, 7)) +
    theme(
      legend.title = element_blank(), legend.position = "bottom",
      axis.text = element_text(size = 7)
    ) +
    scale_x_date(date_labels = "%e.%b") +
    xlab("") +
    ylab(paste0("")) +
    ggtitle(paste0("Daily Mean Tuning Factor Values (ALNUtune) - Model Variant 1"))

  ggsave("vignettes/gg_tune_v1.png", gg_tune_v1, width = 9, height = 6, bg = "white", dpi = 300)

  gg_tune_v1
}

#'
#' Model variant 2 has the same setup as the final variant except:
#'
#' - The change to the tuning factor happens every hour in full (no root)
#'
## ----echo=FALSE, message=FALSE---------------------------------------------------------------------------------------------------------------------------------------------------
if (species_sel == "Alnus") {
  tune_20_v2 <- import_data_cosmo(paste0(here(), "/ext-data/cosmo/20_alnu_tune_v2"), type = "tune")
  tune_21_v2 <- import_data_cosmo(paste0(here(), "/ext-data/cosmo/21_alnu_tune_v2"), type = "tune")

  gg_tune_v2 <- tune_20_v2 %>%
    aggregate_pollen() %>%
    mutate(
      datetime = datetime + years(1),
      season = "season_2020"
    ) %>%
    bind_rows(tune_21_v2 %>%
      aggregate_pollen() %>%
      mutate(season = "season_2021")) %>%
    mutate(value = value / 0.34) %>%
    filter(
      taxon == species_sel,
      station != "Wolfgang",
      between(datetime, startdate, enddate),
    ) %>%
    ggplot() +
    geom_line(aes(
      x = as.Date(datetime),
      y = value,
      col = season
    ), alpha = 0.6) +
    facet_wrap(~station) +
    coord_cartesian(y = c(0, 7)) +
    theme(
      legend.title = element_blank(), legend.position = "bottom",
      axis.text = element_text(size = 7)
    ) +
    scale_x_date(date_labels = "%e.%b") +
    xlab("") +
    ylab(paste0("")) +
    ggtitle(paste0("Daily Mean Tuning Factor Values (ALNUtune) - Model Variant 2"))

  ggsave("vignettes/gg_tune_v2.png", gg_tune_v2, width = 9, height = 6, bg = "white", dpi = 300)

  gg_tune_v2
}

#'
#' Model variant 6 had the same setup as the final variant, but the changes to the tuning factor were
#' calculated based on the past 24h hours only:
#'
## ----echo=FALSE, message=FALSE---------------------------------------------------------------------------------------------------------------------------------------------------
if (species_sel == "Alnus") {
  tune_20_v3 <- import_data_cosmo(paste0(here(), "/ext-data/cosmo/20_alnu_tune_v3"), type = "tune")

  gg_tune_v3 <- tune_20_v3 %>%
    aggregate_pollen() %>%
    mutate(
      datetime = datetime + years(1),
      season = "season_2020",
      value = value / 0.34
    ) %>%
    filter(
      taxon == species_sel,
      station != "Wolfgang",
      between(datetime, startdate, enddate),
    ) %>%
    ggplot() +
    geom_line(aes(
      x = as.Date(datetime),
      y = value,
      col = season
    ), alpha = 0.6) +
    facet_wrap(~station) +
    coord_cartesian(y = c(0, 7)) +
    theme(
      legend.title = element_blank(), legend.position = "bottom",
      axis.text = element_text(size = 7)
    ) +
    scale_x_date(date_labels = "%e.%b") +
    xlab("") +
    ylab(paste0("")) +
    ggtitle(paste0("Daily Mean Tuning Factor Values (ALNUtune) - Model Variant 3"))

  ggsave("vignettes/gg_tune_v3.png", gg_tune_v2, width = 9, height = 6, bg = "white", dpi = 300)

  gg_tune_v3
}

#'
#' Reviewer 2 requested a similar graph faceted by station but for the differences between modelled and measured concentrations.
#'
## ----echo=FALSE------------------------------------------------------------------------------------------------------------------------------------------------------------------
if (species_sel == "Alnus") {
  gg_timeseries_2020 <- data_combined %>%
    pivot_wider(names_from = type) %>%
    mutate(
      baseline = baseline - measurement,
      calibration = calibration - measurement
    ) %>%
    select(-measurement) %>%
    pivot_longer(baseline:calibration, names_to = "type") %>%
    filter(
      taxon == species_sel,
      year(datetime) == 2020
    ) %>%
    ggplot() +
    geom_line(aes(
      x = as.Date(datetime),
      y = value,
      col = type
    ), alpha = 0.6) +
    facet_wrap(~station) +
    theme(
      legend.title = element_blank(), legend.position = "bottom",
      axis.text = element_text(size = 7)
    ) +
    scale_x_date(date_labels = "%e.%b") +
    xlab("") +
    ylab(paste0("")) +
    ylab(paste0("Modelled - Measured Concentration [Pollen / m³]")) +
    ggtitle(paste0("Differences of Measured vs. Modelled Daily Alnus Pollen Concentrations in 2020")) +
    scale_color_manual("", values = col_hex[2:3]) +
    coord_cartesian(y = c(-500, 1500))


  gg_timeseries_2021 <- data_combined %>%
    pivot_wider(names_from = type) %>%
    mutate(
      baseline = baseline - measurement,
      calibration = calibration - measurement
    ) %>%
    select(-measurement) %>%
    pivot_longer(baseline:calibration, names_to = "type") %>%
    filter(
      taxon == species_sel,
      year(datetime) == 2021
    ) %>%
    ggplot() +
    geom_line(aes(
      x = as.Date(datetime),
      y = value,
      col = type
    ), alpha = 0.6) +
    facet_wrap(~station) +
    theme(
      legend.title = element_blank(), legend.position = "bottom",
      axis.text = element_text(size = 7)
    ) +
    scale_x_date(date_labels = "%e.%b") +
    xlab("") +
    ylab(paste0("")) +
    ylab(paste0("Modelled - Measured Concentration [Pollen / m³]")) +
    ggtitle(paste0("Differences of Measured vs. Modelled Daily Alnus Pollen Concentrations in 2021")) +
    scale_color_manual("", values = col_hex[2:3]) +
    coord_cartesian(y = c(-500, 1500))

  ggsave("vignettes/gg_timeseries_2020.png",
    gg_timeseries_2020,
    width = 9, height = 6, bg = "white", dpi = 300
  )
  ggsave("vignettes/gg_timeseries_2021.png",
    gg_timeseries_2021,
    width = 9, height = 6, bg = "white", dpi = 300
  )

  gg_timeseries_2020
  gg_timeseries_2021
}

#'
## ----paper, include=FALSE--------------------------------------------------------------------------------------------------------------------------------------------------------
# Export plots for paper
ggsave("vignettes/gg_tune.png", gg_tune, width = 10, height = 6.5, bg = "white", dpi = 300)
ggsave("vignettes/gg_tune_v1.png", gg_tune_v1, width = 10, height = 6.5, bg = "white", dpi = 300)
ggsave("vignettes/gg_tune_v2.png", gg_tune_v2, width = 10, height = 6.5, bg = "white", dpi = 300)
ggsave("vignettes/gg_tune_v3.png", gg_tune_v3, width = 10, height = 6.5, bg = "white", dpi = 300)
ggsave("vignettes/gg_timeseries_2020.png",
  gg_timeseries_2020,
  width = 10, height = 6.5, bg = "white", dpi = 300
)
ggsave("vignettes/gg_timeseries_2021.png",
  gg_timeseries_2021,
  width = 10, height = 6.5, bg = "white", dpi = 300
)
