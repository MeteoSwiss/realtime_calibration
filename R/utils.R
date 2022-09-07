#' Retrieve Pollendata from txt-file (COSMO output with Fieldextra)
#'
#' @param datapath File path to the txt-file containing the data
#' @param type Hirst, Cosmo or Assim

import_data_cosmo <- function(datapath, type) {
  read_table2(paste(datapath), col_names = TRUE) %>%
    pivot_longer(CHBASE:CHZUER,
      names_to = "station_tag",
      values_to = "value"
    ) %>%
    mutate(
      measurement = case_when(
        str_detect(PARAMETER, pattern = "sdes") ~ "phenology",
        str_detect(PARAMETER, pattern = "saisl") ~ "length",
        str_detect(PARAMETER, pattern = "saisn") ~ "current_day",
        str_detect(PARAMETER, pattern = "tthrs") ~ "tempsum_min",
        str_detect(PARAMETER, pattern = "tthre") ~ "tempsum_max",
        str_detect(PARAMETER, pattern = "ctsum") ~ "tempsum",
        str_detect(PARAMETER, pattern = "tune") ~ "tuning",
        TRUE ~ "concentration"
      ),
      PARAMETER = str_replace_all(PARAMETER, "sdes|saisl|saisn|tthrs|tthre|ctsum|tune", "")
    ) %>%
    inner_join(species, by = c("PARAMETER" = "cosmo_taxon")) %>%
    inner_join(stations, by = c("station_tag" = "cosmo_station")) %>%
    mutate(
      datetime = ymd_h(paste0(
        YYYY, sprintf("%02d", MM),
        sprintf("%02d", DD),
        sprintf("%02d", hh)
      )),
      type = type
    ) %>%
    select(taxon, station, datetime, value, type)
}


#' Retrieve Pollendata from txt-file (DWH)
#'
#' @param datapath File path to the txt-file containing the data
#'
import_data_dwh <- function(datapath) {
  cols <- c("station", "type", "value", "datetime")
  stn_start <- "PLO"
  stn_end <- "PCF"
  meas <- "concentration"
  stn_abbr <- "hirst_station"
  cols <- c("taxon", cols)

  output <- read_delim(paste(datapath), delim = " ", skip = 17) %>%
    pivot_longer({{ stn_start }}:{{ stn_end }}, names_to = "station_short", values_to = "value") %>%
    mutate(
      value = if_else(value == -9999, NA_real_, value),
      type = "measurement",
      measurement = meas,
      datetime = ymd_h(paste0(
        YYYY, sprintf("%02d", MM),
        sprintf("%02d", DD),
        sprintf("%02d", HH)
      ))) %>%
    inner_join(species, by = c("PARAMETER" = "cosmo_taxon")) %>%
    inner_join(stations, by = c("station_short" = stn_abbr)) %>%
    select(cols)
}

#' Aggregate Hourly Data into Daily
#'
#' @param data Data Frame containing hourly concentrations

aggregate_pollen <- function(data, level = "daily") {
  if (level == "daily") {
    data %>%
      mutate(date = lubridate::date(datetime)) %>%
      group_by(taxon, station, date, type) %>%
      summarise(
        nas = sum(is.na(value)),
        value = if_else(nas <= 12, mean(value, na.rm = TRUE), NA_real_)) %>%
      ungroup() %>%
      mutate(
        datetime = ymd_hms(paste0(
          as.character(date),
          "00:00:00"))
      ) %>%
      select(-date, -nas)
  } else if (level == "12h") {
    data %>%
      mutate(
        hour = if_else(hour(datetime) < 12, 12, 0),
        date = if_else(hour == 0, lubridate::date(datetime) + days(1), lubridate::date(datetime))
      ) %>%
      group_by(taxon, station, date, hour, type) %>%
      summarise(
        nas = sum(is.na(value)),
        value = if_else(nas <= 6, mean(value, na.rm = TRUE), NA_real_)) %>%
      ungroup() %>%
      mutate(
        datetime = ymd_hms(paste0(
          as.character(date),
          paste0(sprintf("%02d", hour), ":00:00")
        ))
      ) %>%
      select(-hour, -date, -nas)
  }
}