library(dplyr)
library(here)

# The following types are being modelled:
# Erle - Alder - Aulne - Alnus
# Birke - Birch - Bouleau - Betula
# Gräser - Grasses - Graminées - Poaceae
# Hasel - Hazel - Noiseties - Corylus

species_all <- tibble(
  taxon = c(
    "Castanea",
    "Alnus",
    "Ulmus",
    "Cupressus",
    "Fraxinus",
    "Fagus",
    "Juglans",
    "Plantago",
    "Corylus",
    "Pinus",
    "Quercus",
    "Rumex",
    "Platanus",
    "Populus",
    "Poaceae",
    "Salix",
    "Betula",
    "Carpinus",
    "Urtica",
    "Taxus",
    "Picea",
    "Ambrosia"
  ),
  hirst_taxon = c(
    "kacasth0",
    "kaalnuh0",
    "kaulmuh0",
    "kacuprh0",
    "kafraxh0",
    "kafaguh0",
    "kajuglh0",
    "khplanh0",
    "kacoryh0",
    "kapinuh0",
    "kaquerh0",
    "khrumeh0",
    "kaplath0",
    "kapopuh0",
    "khpoach0",
    "kasalih0",
    "kabetuh0",
    "kacarph0",
    "khurtih0",
    "kataxuh0",
    "kapiceh0",
    "khambrh0"
  ),
  cosmo_taxon = c(
    NA_character_,
    "ALNU",
    NA_character_,
    NA_character_,
    NA_character_,
    NA_character_,
    NA_character_,
    NA_character_,
    "CORY",
    NA_character_,
    NA_character_,
    NA_character_,
    NA_character_,
    NA_character_,
    "POAC",
    NA_character_,
    "BETU",
    NA_character_,
    NA_character_,
    NA_character_,
    NA_character_,
    NA_character_
  ),
  fieldextra_taxon = c(
    NA_character_,
    "ALNU1",
    NA_character_,
    NA_character_,
    NA_character_,
    NA_character_,
    NA_character_,
    NA_character_,
    "CORY1",
    NA_character_,
    NA_character_,
    NA_character_,
    NA_character_,
    NA_character_,
    "POAC1",
    NA_character_,
    "BETU1",
    NA_character_,
    NA_character_,
    NA_character_,
    NA_character_,
    NA_character_
  )
)

species <- species_all %>%
  filter(taxon %in% c("Alnus", "Betula", "Corylus", "Poaceae"))

stations <-
  tibble(
    hirst_station = c(
      "PDS",
      "PBU",
      "PMU",
      "PBS",
      "PZH",
      "PLZ",
      "PBE",
      # "PPY",
      "PNE",
      "PVI",
      "PLS",
      "PGE",
      "PCF",
      "PLO",
      # "BLR",
      "PLU"
    ),
    station = c(
      "Wolfgang",
      "Buchs",
      "Münsterlingen",
      "Basel",
      "Zürich",
      "Luzern",
      "Bern",
      # "Payerne",
      "Neuchâtel",
      "Visp",
      "Lausanne",
      "Genève",
      "La-Chaux-de-Fonds",
      "Locarno",
      # "Balerna",
      "Lugano"
    ),
    cosmo_station = c(
      "CHDAVO",
      "CHBUCH",
      "CHMUEN",
      "CHBASE",
      "CHZUER",
      "CHLUZE",
      "CHBERN",
      # NA_character_,
      "CHNEUC",
      "CHVISP",
      "CHLAUS",
      "CHGENE",
      "CHLACH",
      "CHLOCA",
      # NA_character_,
      "CHLUGA"
    )
  ) %>%
  arrange(hirst_station)

save(species, file = paste0(here(), "/data/other/species.RData"))
save(stations, file = paste0(here(), "/data/other/stations.RData"))

