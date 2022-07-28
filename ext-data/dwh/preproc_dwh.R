# This script pre-processes the raw dwh data and generates the Rdata object
# in the /data folder.

library(magrittr)
library(dplyr)
library(tidyr)
library(readr)
library(lubridate)
library(padr)
library(here)

devtools::load_all()

load(paste0(here(), "/data/other/species.RData"))
load(paste0(here(), "/data/other/stations.RData"))

# The retrieval from the DWH takes roughly 4 hours
# cmd <- paste0("srun -N1 -n1 --time=10:00:00 --partition=postproc --account=s83 sh ",
#   here(), "/ext-data/dwh/dwh_retrieve.sh")
# system(cmd, wait = TRUE, intern = TRUE)

data_dwh <- import_data_dwh(paste0(here(), "/ext-data/dwh/pollen_dwh_hourly.txt")) %>%
  filter(
    taxon == "Alnus" & between(datetime, ymd_hms("2020-01-01 00:00:00"), ymd_hms("2021-12-31 00:00:00")) |
    taxon == "Betula" & between(datetime, ymd_hms("2020-01-01 00:00:00"), ymd_hms("2021-12-31 00:00:00")) |
    taxon == "Corylus" & between(datetime, ymd_hms("2020-01-01 00:00:00"), ymd_hms("2021-12-31 00:00:00")) |
    taxon == "Poaceae" & between(datetime, ymd_hms("2019-01-01 00:00:00"), ymd_hms("2020-12-31 00:00:00"))
  ) %>%
  aggregate_pollen() %>%
  impute_pollen()

save(data_dwh, file = paste0(here(), "/data/dwh/data_dwh.RData"))