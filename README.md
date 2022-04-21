# Overview

This readme document should give an overview of the project and its main components. For a full documentation of the newly developed pollen code to assimilate real-time pollen date into COSMO/ICON please refer to this confluence page: <https://service.meteoswiss.ch/confluence/x/dYQYBQ>
There is e second git-repo available here: <https://github.com/sadamov/cosmo.git> that contains the newly developed Fortran Subroutines to assimilate the real-time pollen data in COSMO/ICON.

## Setup

The analysis was conducted in R 4.1.3. The vignettes are Rmarkdown notebooks.
The project is set up as a minimal R-package to assure maximum reproducibility (<https://r-pkgs.org/index.html>). 

All necessary libraries and R-packages were installed using conda.
This command installs all required libraries into a new conda environment called "new_env": `conda env create -n new_env -f environment.yml`
After creation of the environment, please activate it and you should be ready to rerun the analysis.

## Branches

There is one branch in this repo:
- Main: Protected main branch containing the latest fully functional version of all vignettes and scripts.


## Data

There is quite a lot of external data required for this analysis and only people working at MeteoSwiss will have access to all of it.
In the folder */ext-data* scripts are stored that will retrieve and preprocess the data from external sources. The ready-made data is then stored in the */data* folder and accessed by various scripts.

- **dwh**: Text files of Pollen-Measurements (Concentrations 1/m^3) averaged daily and hourly. Daily surface temperatures. The data is retrieved with the ruby script dwh_retrieve(). Used in various vignettes.
- **cosmo**: Text file of Pollen Concentrations (1/m^3) predicted by COSMO for one specific hour. The data is retrieved with Fieldextra.
- **other**: Dataframes containing names and abbreviations of Swiss pollen stations and species. Manually typed (information available in various MeteoSwiss documentations) and used in various vignettes.
## Vignettes

- **Analysis.Rmd**: The whole analysis in one location, loading the data from the /data folder and calculating all metrics and plots. Can be knitted into html or pdf format for a nicer reading experience.
