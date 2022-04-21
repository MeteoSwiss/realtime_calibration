# This script pre-processes the raw cosmo data and generates the Rdata object
# in the /data folder.

library(magrittr)
library(dplyr)
library(purrr)
library(tidyr)
library(readr)
library(stringr)
library(lubridate)
library(padr)
library(here)
devtools::load_all()

load(paste0(here(), "/data/other/species.RData"))
load(paste0(here(), "/data/other/stations.RData"))

data_baseline <- list.files(paste0(here(), "/ext-data/cosmo/"), pattern = "^[0-9]{2}.*baseline") %>%
  map(~ import_data_cosmo(paste0(here(), "/ext-data/cosmo/", .x), type = "baseline"))

data_calibration <- list.files(paste0(here(), "/ext-data/cosmo/"), pattern = "^[0-9]{2}.*calibration") %>%
  map(~ import_data_cosmo(paste0(here(), "/ext-data/cosmo/", .x), type = "calibration"))

data_cosmo <- data_baseline %>%
  c(data_calibration) %>%
  bind_rows() %>%
  aggregate_pollen() %>%
  impute_pollen() %>%
  # There has been a pollen explosion during those days with absurdly high values in the baseline
  mutate(value = if_else(
    taxon == "Alnus" &
      between(date, as.Date("2021-01-10"), as.Date("2021-01-12")),
    0,
    value
  )) %>%
  select(taxon, station, type, value, datetime, date, hour)

save(data_cosmo, file = paste0(here(), "/data/cosmo/data_cosmo.RData"))