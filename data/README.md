# EPA burn-integrated analysis extracts

The CSV files in this directory are analysis extracts from the U.S. Environmental Protection Agency (EPA) and U.S. Forest Service 2019 and 2021 Missoula Fire Sciences Laboratory burn-integrated sensor-performance dataset.

Dataset: **2019_2021_chamber_study_burn_integrated_data**  
DOI: https://doi.org/10.23719/1532438

Associated publication:

Landis MS, Long R, Krug J, Colón M, Perth J, Urbanski S. Performance evaluation of commercially available non-regulatory instruments and sensors in smoke for PM, CO, CO2, NO2, and SO2. *Journal of the Air & Waste Management Association*. 2026;76(5):383-405. DOI: https://doi.org/10.1080/10962247.2026.2629562.

## Files

- `epa_2019_burn_integrated.csv` — 2019 A10 variables used to parameterize the transport-shift stress scenarios.
- `epa_2021_burn_integrated.csv` — 2021 reference, humidity, raw duplicate-channel, and study-calibrated variables used in the main application and EPA-informed simulations.

## Variables retained

The extracts retain only variables used in the correlation analysis:

- burn identifier;
- chamber relative humidity;
- reference CO and NO2 measurements;
- duplicate Kunak Model A10 CO and NO2 channels;
- duplicate Kunak Air Pro CO and NO2 channels; and
- the study-provided calibrated channels used only for descriptive comparison in the full-sample analysis.

Column names are kept consistent with the analysis code. Measurement units and full metadata should be taken from the EPA DOI record and source documentation.

The original public source data remain available from the DOI above; these extracts are included only to make the reported analysis directly reproducible.
