# eurac---upshift project documentation

**Author:** Marta De Giuli

**Date:** October 2026

**Contact:** marta.degiuli\@eurac.edu

## project overview

This R project contains the primary code and data used to produce datasets, generate figures and tables, and conduct analyses for the *Plant Communications* manuscript entitled:

**“Decoupling the Effect of Reduced Air Pressure and Humidity on Plant Physiology: an ECOTRON Experiment on *Arabidopsis thaliana*”**

The sections below describe the project’s code and data structure and provide context for the analytical workflow. Detailed explanations of individual steps are provided in comments throughout the R scripts. Some code is commented out, including code for quick exploratory assessments of the data. This code has been retained for reference but is not required to reproduce the primary outputs.

Figure references correspond to those used in the published manuscript and its Supporting Information.

## code information

### `plotting.R`

This script generates and saves the main-text figures 1, 2, 3, 4, and supplementary figures 1 and 2 using the datasets explored and analysed in `analysis.R`. It uses a combination of custom functions defined in `functions.R` and the following existing R packages : `ggplot2`, `patchwork`, `grid`, `gridExtra`, `ggpubr` to plot data.

It also generates Table 1 and Supplementary Table 1 - model and data summaries, respectively.

All the data called in this script are found in the `data` directory also found in this project.

### `analysis.R`

This script performs additional analyses, described below:

- **modeling**

- **anova**

### `functions.R`

This script contains the main functions used to simplify exploring and plotting the data in `analysis.R` and `plotting.R`.

The main functions are:

- `base_plot()` creates a baseline scatterplot or boxplot which can be modified to specifications listed in the arguments of the function.

- `plot_interaction_ci()` creates a scatter plot with trend lines according to the specified model.

- `summarySE()` assists in converting complex datasets into plottable data with a specified confidence interval.

## data information

The `data` directory contains data collected over the course of the UPSHIFT experiment conducted in November 2025 by Marta De Giuli, Bouchra El Omari, and Emily Mason in the TerraXcube ecotron chambers in Bolzano, BZ. The datasets included are as follow:

- `beakers.csv` contains the evaporation measurements from our beakers. The columns are explained as follows:

  | Column Name | Explanation                                                   |
  |-----------------|-------------------------------------------------------|
  | beaker      | beaker ID                                                     |
  | treatment   | treatment present in chamber during measurement period        |
  | date_in     | date of start of measurement period                           |
  | time_in     | time of start of measurement period                           |
  | weight_in   | beaker weight at start of measurement period (not tared) in g |
  | date_out    | date of end of measurement period                             |
  | time_out    | time of end of measurement period                             |
  | weight_out  | beaker weight at end of measurement period in g               |
  | elapsed_hrs | elapsed time (hrs) between start and end                      |
  | evap        | evaporation [*start weight - end weight*] in g                |
  | evap_rate   | rate of evaporation [*evaporation / elapsed time*] in g/hr    |

- `co2_treat.csv` contains the real-time information on CO2 inside the chambers during the duration of the experiment

  | Column Name | Explanation                                         |
  |-------------|-----------------------------------------------------|
  | date        | date:time of measurement                            |
  | SEC         | chamber name                                        |
  | co2_ppm     | CO2 concentration in ppm measured in chamber        |
  | treatment   | treatment present in chamber at time of measurement |

- `cn.csv` contains the carbon analyses done on the samples collected at T2 and the d13C values of the air inside each chamber.

  | Column Name | Explanation                                                 |
  |------------------|------------------------------------------------------|
  | Sample_ID   | ID of the sample in the form of treatment_plantid_droughted |
  | N           | percentage of nitrogen in sample                            |
  | sd_N        | standard deviation of N measurement                         |
  | C           | percentage of carbon in sample                              |
  | sd_C        | standard deviation of C measurement                         |
  | d13C        | discrimination of 13C                                       |
  | sd_d13C     | standard deviation of measurement                           |
  | d13Cair     | discrimination of 13C in the chamber air                    |
  | C_discr     | -1\*(d13Cair - d13C)/(1+(d13C/1000))                        |

- `gasexchange.csv` contains the gas exchange measurements at T2 collected with a GFS-3000 from Heinz Walz Gmbh.

  | Column name | Explanation |
  |----|----|
  | Date | date of measurement |
  | Time | time of measurement |
  | Elevation | elevation (treatment) plants grown at |
  | D | pressure normalized VPD for present treatment |
  | Drought | drought status of plant: FALSE = normal, TRUE = drought stressed |
  | Area | area of leaf within cuvette in cm2 |
  | CO2abs | CO2 of ambient air in ppm |
  | H2Oabs | H2O of ambient air in ppm |
  | Flow | flow rate within the cuvette in umol/s |
  | Pamb | ambient pressure in kPa |
  | Tleaf | leaf surface temperature in °C |
  | PARtop | PAR value of the cuvette in umol/m2s |
  | rh | ambient relative humidity |
  | E | measured transpiration |
  | VPD | measured vapour pressure deficit |
  | GH2O | measured stomatal conductance in mmol/m2s |
  | A | net carbon dioxyde assimilation rate in umol/m2s |
  | ci | intercellular CO2 concentration ppm |
  | ca | ambient CO2 concentration ppm |
  | CiCa | internally calculated ci:ca ratio |
  | WUE | water use efficiency calculated as A/E in umol CO2/mmol H2O |

- `pam.csv` contains the chlorophyll fluorescence measurements at T1 and T2 collected from leaf disks through a Closed FluorCam FC 800-C from Photon Systems Instruments.

  | Column Name | Explanation |
  |-----------------------|-------------------------------------------------|
  | sample | sample ID in the form measurement_treatment_time_ID |
  | treatment | treatment sample was gown in |
  | ID | ID of individual sampled |
  | time | T1 or T2, referring to the middle or the end of the experiment |
  | light_intensity | percentage of light power used in actinic light cycle |
  | irr | irradiance of measurement in umol photons /m2\*s |
  | light | identifier of light pulse |
  | Fo | minimum fluorescence in dark-adapted state |
  | Fm | maximum florescence in dark-adapted state (measured after initial saturating light pulse) |
  | Fv | Variable fluorescence (Fm - Fo) |
  | Fv/Fm | Maximum quantum efficiency of PSII photochemistry |
  | QY | Quantum yield of PSII |
  | NPQ | Non-photochemical quenching, proxy for regulated dissipation of absorbed excitation energy as heat |
  | qP | Photochemical quenching coefficient, proxy for proportion of PSII reaction centers that are open |
  | qE | Energy-dependent quenching |
  | ETR | Electron transport rate, estimated from PSII quantum yield and incident PAR (photosynthetically active radiation) calculated as ETR = QY \* 0.5 \* irradiance |
  | QY_max | Maximum quantum yield of PSII for said irradiance |

- `pigments.csv` contains the pigment analyses of leaf disks collected at T1 and T2.

  | Column Name | Explanation |
  |-----------------------|-------------------------------------------------|
  | sample.id | ID of sample (unique within each treatment) |
  | treatment | treatment present in chamber and time of measurement |
  | chlab | measured chlorophyll a : chlorophyll b ratio |
  | chlcar | measured chlorophyll : carotenoid ratio |
  | density | measured chlorophyll density in samples in ug chlorophyll per cm2 |

- `biomass.csv` contains the final aboveground fresh, turgid, and dry biomass of our *Arabidopsis thaliana* individuals.

  | Column Name | Explanation |
  |------------------|------------------------------------------------------|
  | ID | ID of the sample in the form of treatment_plantid_droughted |
  | freshweight | fresh weight of sample in g |
  | turgidweight | turgid weight of sample in g |
  | dryweight | dried weight of sample in g |
  | RWC | relative water content calculated as (freshweight - dryweight)/(turgidweight - dryweight) |

- `biometrics.csv` contains biometric measurements of the individuals taken at T1 and T2.

  | Column Name | Explanation                                                 |
  |------------------|------------------------------------------------------|
  | ID          | ID of the sample in the form of treatment_plantid_droughted |
  | date        | date of measurement                                         |
  | time        | T1 or T2                                                    |
  | height      | maximum height of rosette                                   |
  | diameter    | maximum diameter of rosette                                 |
  | LeafNum     | number of leaves in rosette                                 |
  | LifeLeaf    | number of live leaves in rosette                            |
