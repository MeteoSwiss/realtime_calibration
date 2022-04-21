# This retrieves the daily data for ALNUS, BETULA and POACEAE,CORYLUS for the years 2019-2021
# Running this script inside the /ext-data/dwh folder of the project prepares the data

module purge

# # The UTC 0-0 values are only available from 2017 onwards, hence local time are retrieved
dwh_retrieve -s surface_stations -l 45,48,5,11 -i nat_abbr,PBE,PBS,PBU,PCF,PDS,PGE,PLO,PLS,PLU,PLZ,PMU,PNE,PVI,PZH \
  -t 201901010000-202201010000 -f atab -p 1315,1323,1469,1343 -o pollen_dwh_daily.txt

# This retrieves the hourly values for the same three species for the years 2019-2020
dwh_retrieve -s surface_stations -l 45,48,5,11 -i nat_abbr,PBE,PBS,PBU,PCF,PDS,PGE,PLO,PLS,PLU,PLZ,PMU,PNE,PVI,PZH \
  -t 201901010000-202201010000 -f atab -p 5993,6086,5994,6096 -o pollen_dwh_hourly.txt -c 204

# Separator should always be one blank and not multiple
sed -i 's/  */ /g' pollen_dwh_daily.txt
sed -i 's/  */ /g' pollen_dwh_hourly.txt
# Remove trailing blank at end of lines
sed -i 's/[[:blank:]]*$//' pollen_dwh_daily.txt
sed -i 's/[[:blank:]]*$//' pollen_dwh_hourly.txt
