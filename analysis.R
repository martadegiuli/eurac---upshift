#' Marta De Giuli
#' Eurac Research
#' UPSHIFT
#' Analysis
#' 
#' 26.10.26
  
# Data description ----
  
## Data Source ----
  
#The data used in this report was taken from an experiment run in the TerraX Cube hypobaric chambers in Bolzano, IT in November 2025 on *Arabidopsis thaliana* individuals. Data collected by Marta De Giuli, Bouchra El Omari, and Emily Mason and compiled by Marta De Giuli in 2026.

## Description of the data ----

# Description of the data and the question(s).
# 
# - 4 chambers for a 2 factorial design :
#   - 2 pairs at the same elevation/pressure
#   - 2 pairs at the same D/pressure corrected VPD
# - In each chamber we have 30 pots each containing 3-4 individuals, kept at 100% soil saturation.
# - We have a total of 120 pots across all chambers.
# - We have temporal data for select measurements. A selection of non-destructive measurements were taken halfway through and at the end of the experiment, while other - namely the more destructive ones were taken only at the end of the two week exposure to experimental conditions.
# - Our main goal is to understand *if* and *how* pressure impacts ecophysiological parameters of *Arabidopsis thaliana* and how pressure interacts with humidity. The goal is to disentangle the effect of pressure and humidity to better understand the stressors high alpine plants are subjected to.
# 
# Covariates:
#   
# - Chlorophyll fluorescence measurements, as taken with a FluorCAM 600 from leaf discs:
#   - PSII (Y), NPQ, 1-qL
# 
# - Photosynthetic measurements, as taken with a WALZ device:
#   - Net carbon assimilation, gas exchange (stomatal conductance), evapotranspiration
#   - and calculated: Ci:Ca
# 
# - Organic content and pigment measurements
#   - C:N ratios, d13C values
#   - chlorophyll content, carotene content, chlorophyll a:b ratio
# 
# - Biometric measurements including:
#   - plant height, rosette diameter and leaf number
#   - fresh weight, turgid weight, dry weight, and calculated RWC

#Background processes ----

##Libraries----

#import libraries
library(tidyverse)
library(ggpubr)
library(slider)

#import functions
source("functions.R")
 

## Import data ----

# growth parameters

biomass <- read_csv("data/biomass.csv") |>
  mutate(treatment = substr(ID, 1, 5),
         D         = substr(ID, 4, 5),
         pressure  = case_when(
           substr(ID, 1, 2) == 26 ~ factor("98", levels = c("61","98")),
           substr(ID, 1, 2) == 40 ~ factor("61", levels = c("61","98"))),
         drought   = ifelse(substr(ID, 9, 9) == "d", "d", "nd")
  )

biometric <- read_csv("data/biometrics.csv") |>
  mutate(treatment = substr(ID, 1, 5),
         D         = substr(ID, 4, 5),
         pressure  = case_when(
           substr(ID, 1, 2) == 26 ~ 98,
           substr(ID, 1, 2) == 40 ~ 61),
         drought   = ifelse(substr(ID, 9, 9) == "d", "d", "nd")
  )

## gas exchange data

walz <- read_csv ("data/gasexchange.csv")|>
  dplyr::select(-CiCa) |>  #we will recalculate cica just in case
  mutate(cica = abs(ci) / abs(ca),
         Tleaf = as.numeric(Tleaf),
         pressure = factor(
           case_when(
             Elevation == 260  ~ 98,
             Elevation == 4000 ~ 61
           ),
           levels = c(98, 61)
         )
         ,
         treatment = paste0(substr(Elevation,1,2), "D",D)
  )

## C:N data
cn <- read_csv("data/cn.csv") |>
  mutate(treatment = substr(Sample_ID, 1,5),
         cn = C/N,
         pressure = case_when(
           substr(treatment, 1, 2) == "26" ~ factor("98", levels = c("61","98")),
           substr(treatment, 1, 2) == "40" ~ factor("61", levels = c("61","98"))),
         D = factor(substr(treatment, 4,5), levels = c("10", "20")),
         drought = case_when(
           nchar(Sample_ID) == 9 ~ substr(Sample_ID, 9,9),
           nchar(Sample_ID) == 8 ~ "nd"
         )
  )

## Pigment data
pigments <- read_csv("data/pigments.csv") |> 
  mutate(pressure = case_when(
    substr(treatment, 1, 2) == "26" ~ factor("98", levels = c("61","98")),
    substr(treatment, 1, 2) == "40" ~ factor("61", levels = c("61","98"))),
    D = factor(substr(treatment, 4,5), levels = c("10", "20")),
    time = (substr(treatment, 7,8)),
    treatment = substr(treatment, 1,5))


## Beakers / Evaporation
beakers <- read_csv("data/beakers.csv",
                    col_types = "ccctdctdddd")|>
  mutate(D = substr(treatment, 4, 5),
         pressure  = case_when(
           substr(treatment, 1, 2) == 26 ~ 98,
           substr(treatment, 1, 2) == 40 ~ 61))


## PAM data

# #Note: PAM parameters are defined as the following:
# F0	                              Minimal chlorophyll fluorescence intensity measured in the dark-adapted state, when all PSII reaction centres are open
# F	                                Average of Ft in 3 s , when a SP is triggered
# Fm	                              Maximal chlorophyll fluorescence intensity measured in the dark-adapted state during the application of a saturating pulse of light
# Fv/Fm	                            Maximum quantum yield of PSII photochemistry measured in the dark-adapted state
# Y(II)= (Fm'− Fs)/Fm'	            Effective quantum yield of photochemical energy conversion in PSII
# NPQ=(Fm − Fm')/Fm'	              Stern-Volmer non-photochemical quenching coefficient
# qP=(Fm'− Fs)/(Fm'− F0′)	          Coefficient of photochemical quenching of variable fluorescence based on the puddle model of PSII
# Y(NO) = 1/(PNQ + 1 + qL(Fm/F0-1))	Quantum yield of non-regulated energy dissipation in PSII

pam <- read_csv("data/pam.csv") |>
  mutate(treatment = substr(sample,5,9),
         D         = as.numeric(substr(treatment, 4, 5)),
         pressure  = case_when(
           substr(treatment, 1, 2) == 26 ~ as.factor(98),
           substr(treatment, 1, 2) == 40 ~ as.factor(61)),
         drought   = case_when(
           time == "T2" & substr(ID,3,3) == "d" ~ "d",
           .default = "nd"
         ))

co2_treat <- read.csv("data/co2_treat.csv") |>
  mutate(date = dmy_hm(date, tz = "Europe/Rome"),
         .keep = "unused")

## Defaults ----

#color defaults for some plots
cols4 <- c("26D10" = "#fde725", 
           "26D20" = "#35b779", 
           "40D10" = "#31688e", 
           "40D20" = "#440154")

cols2 <- c("98" = "#35b779", 
           "61" = "#440154")

cols2a <- c("10" = "#932667",
            "20" = "#fca50a" )

my.comps <- list(c("26D10", "26D20"), c("40D10", "40D20"), c("26D10", "40D10"), c("26D20", "40D20"))

#methods for stat tests:
method1 = "anova"
method2 = "t.test"
method3 = "wilcox.test" #for non-parametric data
 

# Analysis ----

# Before starting with the deeper analyses, it is important to check that our treatments *pressure* x *D* (normalized VPD) worked as expected. To do so, we placed two beakers in each chamber and measured evaporation throughout the experiment.

## Evaporation----

#add 'comparison' column to be able to compare across factors
beakers.comp <- add_comp(beakers, "pressure", "D")

#test the distribution of the data
shapiro.test(beakers$evap_rate)
#p-value = 0.007278, thus < 0.05 and so the data is not normally distributed - use wilcox test for significance

#plot to get an overview
beakers.plot<-
  beakers.comp |>
  ggplot(aes(x = treatment, y = evap_rate, fill = treatment))+
  geom_boxplot(data = beakers, color = "grey",fill = "grey")+ 
  geom_boxplot()+
  geom_jitter(width=0.01)+
  labs(x = NULL,
       y = expression(Evaporation~Rate~"["*mm~hr^-1*"]"))+
  scale_fill_manual(values= c("26D10" = "#fde725", 
                              "26D20" = "#35b779", 
                              "40D10" = "#31688e", 
                              "40D20" = "#440154"))+
  stat_compare_means(method = method3, 
                     label = "p.signif"
                     #label.y = 40,
                     #label.x = c(2, 3, 3.5, 1) #these need to be corrected
  )+
  stat_compare_means(method = method3, 
                     label = "p.format",
                     label.y.npc = "bottom"
  )+
  facet_grid(~comparison, drop = TRUE)+
  #ylim(NA, 45)+
  theme_bw( )+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), 
        legend.position = "bottom")

beakers.plot

#clean up
rm(beakers.comp)

# As expected, we observe a strong statistical difference between humidities at the same pressure (panels P= 61 and P = 98), with chambers at D = 20 showing higher evaporation. Slight differences are also present between elevations at the same humidity, with chambers at lower pressure generally showing a lower evaporation rate (specifically D = 10, where we observe a significant decrease, p = 0.023). Thus we will move forward with our analyses.

## Chamber parameters ----
# Concentration of $CO_2$ in ppm of each chamber, the data is a 30min running mean of the raw (1min resolution) data. Colors indicate treatment present in each chamber, red shaded area indicate periods when measurements were taking place (expect higher ppm $CO_2$ due to human presence), line type indicates the chamber, and red vertical bars indicate rotations of treatments between the chambers.

#create outline of human presence in chambers (i.e. when measurements were taken)
meas <- data.frame(x_start = parse_date_time(c("2025-11-03 15:00:00",
                                               "2025-11-06 8:00:00",
                                               "2025-11-10 8:30:00",
                                               "2025-11-11 8:30:00", 
                                               "2025-11-12 8:30:00" ), "ymd HMS"),
                   x_end = parse_date_time(c("2025-11-03 17:00:00", 
                                             "2025-11-06 16:00:00",
                                             "2025-11-10 12:00:00",
                                             "2025-11-11 16:00:00",
                                             "2025-11-12 16:00:00"), "ymd HMS"))

#create running mean of CO2 measurement
co2_rm <- co2_treat |>
  arrange(SEC, date) |>
  group_by(SEC) |>
  mutate(co2_runmean = slide_dbl(co2_ppm,
                                 mean,
                                 .before = 30,
                                 .after  = 30,
                                 .complete = TRUE))

#create segment to extend running mean to indicate chamber number
segments <- co2_rm |>
  filter(!is.na(co2_runmean)) |>
  group_by(SEC) |>
  filter(date == max(date)) |>
  ungroup() |>
  mutate(
    xend = date + hours(30),
    yend = co2_runmean
  )

#plot
ppm_co2 <- ggplot(co2_rm, aes(x=date, y= co2_runmean, linetype = SEC, color = treatment))+
  geom_rect(data = meas,
            aes(xmin = x_start,
                xmax = x_end,
                fill = " "),
            ymin = -Inf,
            ymax = Inf,
            inherit.aes = F,
            alpha = 0.5
  )+
  geom_line(linewidth = 1)+
  #geom_line(aes(linetype = SEC), color = "black")+
  geom_segment(data = segments, 
               aes(x = date, xend = xend, y = co2_runmean, yend = co2_runmean, linetype = SEC), 
               color = "black")+
  geom_vline(xintercept = parse_date_time(c("2025-11-03 16:00:00", "2025-11-06 16:00:00", "2025-11-10 10:00:00" ), "ymd HMS"), 
             color = "red",
             linewidth = 1)+
  geom_text(data = segments, 
            aes(label = SEC, x = date + hours(24), y = co2_runmean + 10), 
            color = "black")+
  scale_color_manual(values = cols4)+
  scale_fill_manual(values = c(" " = "grey"))+
  coord_cartesian(ylim = c(300,550),
                  xlim = c(parse_date_time("2025-10-30 10:00:00", "ymd HMS"), parse_date_time("2025-11-14 00:00:00", "ymd HMS")))+
  scale_x_datetime(date_breaks = "3 days", 
                   date_labels = "%d %b")+
  labs(y = "CO2 concentration [ppm]",
       fill = "measurements")+
  guides(linetype = "none")+
  theme_bw()+
  theme(legend.position = "bottom")

ppm_co2

# There is variability in the $CO_2$ concentration between the chambers, with SEC2 having a consistently higher concentration of $CO_2$ than the others. Each treatment was exposed to each chamber, and photosynthetic parameters were taken using the in-situ ppm of $CO_2$ (instantaneous reading).

rm(meas, co2_rm, segments)

## Chlorophyll Fluorescence ----

### Sensitivity to treatment factors ----

# In previous iterations, we have observed that humidity has a smaller impact on chlorophyll fluorescence parameters than pressure. To check whether this is true in this dataset, we can compare four nested models:
# response ~ irradiance + pressure
# response ~ irradiance + pressure + time
# response ~ irradiance + pressure + D + time
# response ~ irradiance + pressure + D + pressure * D + time
# 
# This will also enable us to understand whether there are general *interacting effects* of pressure and humidity. Here we look at the case for the yield of PSII.

#use only data from maximum irradiance
subset.pam <- pam |>
  filter(light == "Lss")

psii.model<- lm(QY ~ irr + pressure + D + pressure*D + time,
                      data = subset.pam)
psii.model.b<-lm(QY ~ irr + pressure + D + time,
                     data = subset.pam)
psii.model.c<-lm(QY ~ irr + pressure + time,
                     data = subset.pam)
psii.model.d<-lm(QY ~ irr + pressure,
                     data = subset.pam)

anova(psii.model.d, psii.model.c, psii.model.b, psii.model)

# Since the explanatory power of the models doesn't decrease with the removal of the interaction + humidity terms, it seems that the best model is the one that doesn't include humidity at all. This indicates that we can make our comparisons looking simply at the difference between pressure treatments, aggregating our observations by pressure. Based on these results, it would be statistically sound not to differentiate by measurement time as no *significant* change in explanatory power is seen. This tends towards pseudo-replication though, and would likely become problematic to defend to a potential reviewer.

#clean up
rm(subset.pam, psii.model, psii.model.b, psii.model.c, psii.model.d)

#NOT RUN, for redundancy. But here you see that the same result is found for the other response components.

# for(response in c("NPQ", "ETR", "qE", "qP")){
#   subset.pam <- pam |>
#     filter(light == "Lss")
#   
#   print(paste("ANOVA for", response))
#   
#   model <- lm(get(response) ~ irr + pressure + D + pressure*D + time,
#               data = subset.pam)
#   model.b<-lm(get(response) ~ irr + pressure + D + time,
#               data = subset.pam)
#   model.c<-lm(get(response) ~ irr + pressure + time,
#               data = subset.pam)
#   model.d<-lm(get(response) ~ irr + pressure,
#               data = subset.pam)
#   
#   #print(anova(model, model.b, model.c, model.d))
#   print(anova(model))}
#
# #clean up
# rm(subset.pam, model, model.b, model.c, model.d)


### PSII(Y) ----

#generate aggregated summary table for PSII(Y)
psii.agg <- data.frame(pressure  = factor(),
                       irr       = as.numeric()  ,
                       N         = as.numeric()  ,
                       QY        = as.numeric()  ,
                       sd        = as.numeric()  ,
                       se        = as.numeric()  ,
                       ci        = as.numeric()  ,
                       time      = as.character()
)

#summarize using defined function
for(sampling_time in unique(pam$time)){
  psii.agg <- rbind(psii.agg,
                    summarySE(data = pam |> 
                        filter(light=="Lss", 
                               drought == "nd", 
                               time == sampling_time), 
                      "QY", 
                      groupvars = c("pressure", "irr"))|>
  mutate(time = sampling_time,
         pressure = as.factor(pressure)))
}


##calculate significance statistics - type of test (wilcox for non-parametric option or t test for parametric option is chosen by the function depending on the distribution of the data)

df_pvalues.psii<-  signif(pam|>filter(light == "Lss", drought == "nd"), 
                          comparison_column = "pressure", 
                          variable = "QY")

#take a peak to get a sneak preview of our results
#df_pvalues.psii

#We see that at lower pressure we have a higher quantum efficiency of PSII at both measurement times.

# more general plot
# psii.agg.plot <-  plot_treatment_comparison(psii.agg, 
#                             xvals = "irr", yvals = "QY", 
#                             group = "pressure", 
#                             facets = c("time"), 
#                             type = "line", 
#                             pvalues = df_pvalues.psii, 
#                             y.axis = "PSII (Y)", x.axis = expression(Light~"["*"\u03BC"*E*"]"), 
#                             colscale = "cols2")
# 
# psii.agg.plot


### NPQ ----

#generate aggregated summary table for NPQ
npq.agg <- data.frame(pressure  = factor(),
                       irr       = as.numeric()  ,
                       N         = as.numeric()  ,
                       NPQ        = as.numeric()  ,
                       sd        = as.numeric()  ,
                       se        = as.numeric()  ,
                       ci        = as.numeric()  ,
                       time      = as.character()
)

#summarize using defined function
for(sampling_time in unique(pam$time)){
  npq.agg <- rbind(npq.agg,
                    summarySE(data = pam |> 
                        filter(light=="Lss", 
                               drought == "nd", 
                               time == sampling_time), 
                      "NPQ", 
                      groupvars = c("pressure", "irr"))|>
  mutate(time = sampling_time,
         pressure = as.factor(pressure)))
}


##calculate significance statistics - type of test (wilcox for non-parametric option or t test for parametric option is chosen by the function depending on the distribution of the data)

df_pvalues.npq <-  signif(pam|>filter(light == "Lss", drought == "nd"), comparison_column = "pressure", variable = "NPQ" )
#note the 'ties' mentioned in the warning refer to the 0 irradiance values where NPQ is 0.

#sneak preview
#df_pvalues.npq

#npq values seem to be different beween pressures

#plot
# npq.agg.plot <-  plot_treatment_comparison(npq.agg, xvals = "irr", yvals = "NPQ",
#                            group = "pressure",
#                            facets = c("time"),
#                            type = "line",
#                            pvalues = df_pvalues.npq,
#                            y.axis = "NPQ", x.axis = expression(Light~"["*"\u03BC"*E*"]"),
#                            colscale = "cols2")
# 
# npq.agg.plot

#they're higher at high pressure

### ETR ----

#dataframe to fill for summary data
etr.agg <- data.frame(pressure  = factor(),
                       irr       = as.numeric()  ,
                       N         = as.numeric()  ,
                       ETR        = as.numeric() ,
                       sd        = as.numeric()  ,
                       se        = as.numeric()  ,
                       ci        = as.numeric()  ,
                       time      = as.character()
)

#summarize using defined function
for(sampling_time in unique(pam$time)){
  etr.agg <- rbind(etr.agg,
                    summarySE(data = pam |> 
                        filter(light=="Lss", 
                               drought == "nd", 
                               time == sampling_time), 
                      "ETR", 
                      groupvars = c("pressure", "irr"))|>
  mutate(time = sampling_time,
         pressure = as.factor(pressure)))
}


##calculate significance statistics

df_pvalues.etr<-  signif(pam|>filter(light == "Lss", drought == "nd"), 
                         comparison_column = "pressure", 
                         variable = "ETR" )

#warnings due to expected 0 values across all treatments at irradiance = 0.

#sneakpeak of significances
#df_pvalues.etr

#plotting
# etr.agg.plot <-  plot_treatment_comparison(etr.agg, 
#                            xvals = "irr", yvals = "ETR", 
#                            group = "pressure", 
#                            facets = c("time"), 
#                            type = "line", 
#                            pvalues = df_pvalues.etr, 
#                            y.axis = "ETR", x.axis = expression(Light~"["*"\u03BC"*E*"]"), 
#                            colscale = "cols2")
# 
# etr.agg.plot

#etr is lower at higher pressure

### qE ----

#dataframe to fill for summary data
qE.agg <- data.frame(pressure  = factor(),
                       irr       = as.numeric()  ,
                       N         = as.numeric()  ,
                       qE        = as.numeric()  ,
                       sd        = as.numeric()  ,
                       se        = as.numeric()  ,
                       ci        = as.numeric()  ,
                       time      = as.character()
)

#summarize using defined function
for(sampling_time in unique(pam$time)){
  qE.agg <- rbind(qE.agg,
                    summarySE(data = pam |> 
                        filter(light=="Lss", 
                               drought == "nd", 
                               time == sampling_time), 
                      "qE", 
                      groupvars = c("pressure", "irr"))|>
  mutate(time = sampling_time,
         pressure = as.factor(pressure)))
}


##calculate significance statistics

df_pvalues.qE<-  signif(pam|>filter(light == "Lss", drought == "nd"), 
                        comparison_column = "pressure", 
                        variable = "qE" )
#sneak preview
#df_pvalues.qE

#plot!

# qE.agg.plot <-  plot_treatment_comparison(qE.agg, 
#                           xvals = "irr", yvals = "qE", 
#                           group = "pressure", 
#                           facets = c("time"), 
#                           type = "line", 
#                           pvalues = df_pvalues.qE, 
#                           y.axis = "qE", x.axis = expression(Light~"["*"\u03BC"*E*"]"), 
#                           colscale = "cols2")
# 
# qE.agg.plot
 
#higher values at higher pressure

### 1- qP ----

#dataframe to fill for summary data
qP.agg <- data.frame(pressure  = factor(),
                       irr       = as.numeric()  ,
                       N         = as.numeric()  ,
                       qP        = as.numeric()  ,
                       sd        = as.numeric()  ,
                       se        = as.numeric()  ,
                       ci        = as.numeric()  ,
                       time      = as.character()
)

#summarize using defined function
for(sampling_time in unique(pam$time)){
  qP.agg <- rbind(qP.agg,
                    summarySE(data = pam |> 
                        filter(light=="Lss", 
                               drought == "nd", 
                               time == sampling_time), 
                      "qP", 
                      groupvars = c("pressure", "irr"))|>
  mutate(time = sampling_time,
         pressure = as.factor(pressure)))
}


##calculate sifnificance statistics

df_pvalues.qP<-  signif(pam|>filter(light == "Lss", drought == "nd"), 
                        comparison_column = "pressure",
                        variable = "qP" )

#sneak peak & see the differences
#df_pvalues.qP

#convert to 1 - qP for plotting
qP.agg <- qP.agg |> 
  dplyr::mutate(qP = case_when(
    qP != 0 ~ 1-qP,
    TRUE ~ qP
  ))

#plot
# qP.agg.plot <-  plot_treatment_comparison(qP.agg,
#                           xvals = "irr", yvals = "qP", 
#                           group = "pressure", 
#                           facets = c("time"), 
#                           type = "line", 
#                           pvalues = df_pvalues.qP, 
#                           y.axis = "qP", x.axis = expression(Light~"["*"\u03BC"*E*"]"),
#                           colscale = "cols2")
# 
# qP.agg.plot




# We see aligned trends across all the chlorophyll fluorescence parameters aggregated by pressure: *Arabidopsis thaliana* plants exposed to low pressure had a higher quantum yield of PSII, lower NPQ, a higher electron transport rate and lower qP and qE. These results are consistent in indicating a higher PSII operating efficiency, and a reduced reliance on photoprotective energy dissipation mechanisms in plants grown at low pressure.

## Gas Exchange ----

### Leaf temperature ----

# We know that these measurements depend on temperature as well here we provide the leaf temperature data - which varies from 25.0°C to 24.55°C across treatments. [NOT RUN - but you can for your peace of mind]

#test for normality

#shapiro.test(walz$Tleaf)
#p-value = 0.052, so the data is (barely) normally distributed and we can use parametric tests for our statistics

#base facet plot

# tleaf.plot <- walz|>
#   ggplot(aes(x = treatment, y = Tleaf, fill = treatment))+
#   geom_boxplot()+
#   geom_jitter(width=0.01)+
#   labs(x = NULL,
#        y = "Temperature [°C]")+
#   scale_fill_manual(values= cols4)+
#   stat_compare_means(comparisons = my.comps,
#                      label = "p.signif",
#                      method = method2,
#                      label.y = c(24.9, 24.9, 25,25.07),
#                      paired = F
#   )+
#   stat_compare_means(label.y = 24.6,
#                      label.x.npc = "middle")+
#   theme_bw( )+
#   theme(axis.text.x = element_blank(),
#         axis.ticks.x = element_blank(),
#         legend.position = "bottom")
# 
# tleaf.plot

#clean up
#rm(tleaf.plot)

#We do see a significant difference between some plots but at the scale of 0.1 °C, physiologically insignificant.

### Sensitivity to treatments----

# To better understand the interaction of our manipulated parameters, it would be wise to look at how a similar model scheme compares for our gas exchange parameters.

#### Assimilation Rate ----

# For assimilation rate, we will compare the following three nested models:
# assimilation ~ pressure \\
# assimilation ~ pressure + D\\
# assimilation ~ pressure * D

assim.model.c <- lm(A ~ pressure + D + pressure*D,
                        data = walz)
assim.model.b <- lm(A ~ pressure + D,
                        data = walz)
assim.model.a <- lm(A ~ pressure,
                    data = walz)

anova(assim.model.a, assim.model.b, assim.model.c) #double check comparing all three models

anova(assim.model.c) #most of the variance explained only by pressure
 

# Comparing the three models, we observe no loss in explanatory power with the removal of *D* and the interaction term. Selecting only the best fit model ($assimilation ~ pressure * D$), we can see that a large majority of the variance (\~94.4%) is explained only by pressure. The other two effects (*D* and the interaction between *pressure* and *D*) don't have a significant impact on the model.

# assim.plot <-  plot_treatment_comparison(walz,
#            xvals = "pressure", yvals = "A", 
#            group = "pressure", 
#            colscale = "cols2", 
#            type = "boxplot", 
#            y.axis = expression(Assimilation~Rate~"["*"\u03BC"*mol~m^-2~s^-1*"]"))

# assim.plot

#A is significantly lower at lower pressure

#clean up
rm(assim.model.a, assim.model.b, assim.model.c)

#### gH2O, cica, E ---- 

# Above, we showed that the interaction between *D* and *pressure* and the fixed effect of *D* don't increase the explanatory power of the response of assimilation. For the remaining parameters, namely transpiration, Ci:Ca, and stomatal conductance, we find that the best fit model is $response ~ pressure * D$ .

#[NOT RUN - but you can if it pleases you]
# for(response in c("cica", "E", "GH2O")){
#   
#   cat(paste("\n ANOVA for", response," \n \n"))
#   
#   model.c <- lm(get(response) ~ pressure + D + pressure*D,
#                 data = walz)
#   model.b <- lm(get(response) ~ pressure + D,
#                 data = walz)
#   model.a <- lm(get (response) ~ pressure,
#                 data = walz)
#   print(anova(model.a, model.b, model.c))
#   
#   #print(autoplot(model.c)) # just a check that a gaussian model is correct for this 
#   
#   print(anova(model.c)) # if we want to look more closely
# }
# rm(model.a, model.b, model.c)
# 
# What we observe there is that the interaction between *D* and *pressure* is fundamental in being able to explain our results for the remaining responses. Importantly (seen in the ANOVAs for the individual best-fit model) we see that *D* has no impact by itself, rather it modulates the impact of *pressure*. We can't aggregate for pressure in these cases due to the importance of the interaction term, rather we should try to understand how *D* and *pressure* interact, looking at the slope of their interaction.

### Pressure : D interaction ----

#modify pressure and D in walz dataframe to be able to use them correctly in the lm
walz_f <- walz|>
  mutate(pressure_n = as.numeric(as.character(pressure)),
         D_f = factor(D, levels = c(10,20)))

# Fit interaction model with raw pressure
model.GH2O <- lm(GH2O ~ pressure_n * D_f, data = walz_f)

model.cica <- lm(cica ~ pressure_n * D_f, data = walz_f)

model.E <- lm(E ~ pressure_n * D_f, data = walz_f)

#see direction of interactions

#model.GH2O
#model.cica
#model.E
 
# While the response changes with pressure, the direction of that change depends on D. We observe that pressure increases gas exchange, transpiration, and Ci:Ca under high *D* (low air water content) and decreases the responses of these variables under low *D* (high air water content).

## Pigments (chlorophyll data) ----

### Sensitivity to treatment components ----

for(response in c("chlab", "chlcar", "density")){
  
  cat(paste("\n ANOVA for", response," \n \n"))
  
  model.c <- lm(get(response) ~ pressure + D + pressure*D + time,
                data = pigments)
  model.b <- lm(get(response) ~ pressure + D + time,
                data = pigments)
  model.a <- lm(get (response) ~ pressure + time,
                data = pigments)
  model.a2 <- lm(get (response) ~ D + time,
                 data = pigments)
  model.t <- lm(get (response) ~ time,
                data = pigments)
  
  print(anova(model.t, model.a, model.b, model.c))
  print(anova(model.t, model.a2, model.b, model.c))
  
  #print(autoplot(model.c)) # just a check that a gaussian model is correct for this
  
  #print(anova(model.c)) # if we want to look more closely
}

# There is no statistical evidence of an interaction between *pressure* and *D*, so we can drop all models which include *pressure * D*.

#clean up
rm(model.t, model.a, model.a2, model.b, model.c)

### Chlorophyll A:B ratio ----

#run models to compare sensitivity to treatment elements
model.chlab1 <- lm(chlab ~ pressure + time,
                        data = pigments)
model.chlab2 <- lm(chlab ~ pressure + D + time,
                        data = pigments)

anova(model.chlab1, model.chlab2)

model.chlab3 <- lm(chlab ~ pressure*time,
                   data = pigments)

anova(model.chlab1, model.chlab3)
anova(model.chlab3)


# The chlorophyll A:B ratio seems to be influenced mostly by pressure with a slight, but non significant effect of the interaction of *time*. This seems to indicate that for this variable we could aggregate pressure treatments and time, but returning to our "pseudoreplication" issue, we will maintain our two sampling times separate.

#plot it 
 plot_treatment_comparison(pigments,
           xvals = "pressure", yvals = "chlab", 
           group = "pressure", 
           colscale = "cols2", 
           type = "boxplot", 
           y.axis = "Chlorophyll A:B Ratio")

#clean up
rm(model.chlab1, model.chlab2, model.chlab3) 

#A:B is usually used as an indicator of light availability, with shade-acclimated leaves often showing a lower ratio due to a larger investment in light-harvesting complexes. It can also be interpreted as an attempt to maximize light capture. In our scenario, this seems to indicate that at lower pressure, *Arabidopsis* plants are investing more in light harvesting than processing machinery - a potential response to the increased efficiency of the photosynthetic apparatus.

### Chlorophyll:carotenoid ratio ----

#run models to compare sensitivity to treatment elements
model.chlcar1 <- lm(chlcar ~ D + time,
                        data = pigments)
model.chlcar2 <- lm(chlcar ~ pressure + D + time,
                        data = pigments)

anova(model.chlcar1, model.chlcar2)

model.chlcar3 <- lm(chlcar ~ D*time,
                   data = pigments)

anova(model.chlcar1, model.chlcar3)
anova(model.chlcar1)


#We have a different effect - *pressure* no longer seems to play a role, and rather we have a strong effect of *D* and a minimal effect of *time*, though their interaction does not have a statistically significant impact. It could be interesting to delve into this temporal trend. To do so, let's look at a temporal comparison based on treatments.

ggplot(pigments, aes(x= treatment, y = chlcar, fill = time))+
  geom_boxplot()+
  stat_compare_means(label = "p.signif")+
  labs(y = "Chlorophyll : Carotenoid",
       x = NULL,
       fill = NULL)+
  theme_bw()

#Cool! There's a statistically significant decrease in the ratio at low pressure with time. This could indicate a decrease in stress response as carotenoids are photoprotective pigments and can indicate an investment in stress tolerance. To confirm this interpretation, we have to look at the chlorophyll density trends:

#clean up
rm(model.chlcar1, model.chlcar2, model.chlcar3)
 
### Chlorophyll density ----

#run models to compare sensitivity to treatment elements
model.chlden1 <- lm(density ~ D + time,
                        data = pigments)
model.chlden2 <- lm(density ~ pressure + D + time,
                        data = pigments)

anova(model.chlden1, model.chlden2)

model.chlden3 <- lm(density ~ D*time,
                   data = pigments)

anova(model.chlden1, model.chlden3)
anova(model.chlden1)

 
#The case for chlorophyll density is similar, but cleaner, with the only fixed effect having a significant impact on the variable being *D*. Time doesn't seem to have an impact. But let's double check, to understand our chlorophyll:carotenoid ratios better.


#plot
ggplot(pigments, aes(x= treatment, y = density, fill = time))+
  geom_boxplot()+
  stat_compare_means(label = "p.signif")+
  labs(y = expression(Chlorophyll~Density~"["*"\u03BC"*g~chl~cm^-2*"]"),
       x = NULL,
       fill = NULL)+
  theme_bw()
 

#While we seem to observe a temporal trend - with an increase of chlorophyll at high pressure and a decrease at low pressure over time, it is not statistically significant. This lack of statistical significance is further supported by the anovas above, indicating that the most likely interpretation of the temporal signal in the chlorophyll:carotenoid ratio is a decrease in carotenoids. Thus a decreased stress response - since the decrease in the chlorophyll:carotenoid ratio is only observed in the low pressure treatments, this is potentially an adaptation to pressure.

#clean up
rm(model.chlden1, model.chlden2, model.chlden3)

## Carbon content - C:N ratio and d13C ----

#modifying dataframe for easier faceting
cn.comp <- add_comp(cn, "pressure", "D")
 

### C:N Ratio ----

#run models to compare sensitivity to treatment elements
model.cn1 <- lm(cn ~ pressure + D + pressure*D, data = cn)
model.cn2 <- lm(cn ~ pressure + D , data = cn)
model.cn3 <- lm(cn ~ pressure, data = cn)
model.cn4 <- lm(cn ~ D, data = cn)

anova(model.cn3, model.cn2, model.cn1)
anova(model.cn4, model.cn2, model.cn1)

#deeper look
#anova(model.cn2)

#There is no statistical evidence of an interaction between *pressure* and *D*, but in both anova model comparisons, the model including both *pressure* and *D* performs much better than the one with a single effect. The anova of the single selected model - *cn \~ pressure + D* , indicates that D has a stronger effect than pressure (F = 18, F = 7). This means we cannot aggregate for pressure treatment (or for D treatment).

#check for normality
shapiro.test(cn$cn)
#p-value = .067, so the data is normally distributed and we can use parametric 
#tests for our statistics
# 
# cat(paste("Since the p-value is", 
#       ifelse(shapiro.test(cn$cn)$p.value > 0.05, "non-significant,", "significant,"), 
#       "the data is", 
#       ifelse(shapiro.test(cn$cn)$p.value > 0.05, "normally distributed \n and we will use t-test for comparison going forward.", 
#              "non-normally distributed \n and we will use the Wilcox test for comparison going forward.")))

 plot_treatment_comparison(cn.comp|>filter(drought=="nd",is.na(cn)==F), 
           cn|>filter(drought=="nd",is.na(cn)==F), 
           xvals = "treatment", yvals = "cn", 
           group = "treatment", 
           colscale = "cols4", 
           facets = "comparison", 
           type = "boxplot", 
           y.axis = "C:N" )
 

#clean up
rm(model.cn1, model.cn2, model.cn3, model.cn4)

#There seems to be a modulating signal of *D* on the effect of *pressure* on cn (also observed below in the d13C). Generally, we observe an increase in C:N with a decrease in pressure an increase in D. This signal is significant only in the more stressed plants (i.e. low pressure - P = 61, and high D - D = 20), but is present across all analyses.

### d13C ----

#run models to compare sensitivity to treatment elements
model.d13c1 <- lm(C_discr ~ pressure + D + pressure*D, data = cn|>filter(drought=="nd",is.na(C_discr)==F))
model.d13c2 <- lm(C_discr ~ pressure + D , data = cn|>filter(drought=="nd",is.na(C_discr)==F))
model.d13c3 <- lm(C_discr ~ pressure, data = cn|>filter(drought=="nd",is.na(C_discr)==F))
model.d13c4 <- lm(C_discr ~ D, data = cn|>filter(drought=="nd",is.na(C_discr)==F))

anova(model.d13c3, model.d13c2, model.d13c1)
anova(model.d13c4, model.d13c2, model.d13c1)

anova(model.d13c1)
 

#*Pressure*, *D*, and their interaction have strong effects on $\delta^{13}C$. We sadly then cannot aggregate by either *pressure* or *D*, but can take a look at the interaction term as well as a closer look at how each effect impacted d13C .


#check for normality
shapiro.test(abs(cn$C_discr))
#p-value = 2.841e-05, so the data is not normally distributed and we cannot use parametric tests for our statistics


 plot_treatment_comparison(cn.comp|>filter(drought=="nd",is.na(cn)==F), 
                        cn|>filter(drought=="nd",is.na(cn)==F), 
                        xvals = "treatment", yvals = "C_discr", 
                        group = "treatment", 
                        colscale = "cols4",
                        facets = "comparison",
                        type = "boxplot",
                        y.axis = expression(delta^{13}*C) )

#modifying cn for the function:

 cn_mod <- cn|>
   filter(drought=="nd",
          is.na(C_discr)==F)|>
   mutate(pressure_n = as.numeric(as.character(pressure)),
          D_f = factor(D, levels = c("10", "20")))

model.d13c1 <- lm(C_discr ~ pressure_n*D_f, data = cn_mod) #keep for plotting later

plot_interaction_ci(
  data = cn_mod,
  model = model.d13c1,
  y_var = "C_discr",
  x_var = "pressure_n",
  group_var = "D_f",
  y_label = expression(delta^{13}*C),
  x_label = "Pressure [kPa]",
  legend = "D",
  colscale = "cols2a"
)

#clean up
rm(model.d13c2, model.d13c3, model.d13c4)

## Biometric Data ----

#modifying dataframe for easier faceting

biomass.comp <- add_comp(biomass, "pressure", "D")

biometric.comp <- add_comp(biometric, "pressure", "D")

### Biomass ----

#test for normality
shapiro.test(biomass$dryweight)
#p-value = 0.003224,  so the data is not normally distributed and we must use non parametric tests for our comparisons (i.e. Wilcox i.e. `method 3`)

biomass.plot<- plot_treatment_comparison(biomass.comp |> filter(drought=="nd", !is.na(dryweight)),
                                         biomass |> filter(drought=="nd", !is.na(dryweight)),
                                         xvals = "treatment", yvals = "dryweight",
                                         group = "treatment",
                                         colscale = "cols4",
                                         facets = "comparison",
                                         type = "boxplot",
                                         y.axis = "Aboveground dry biomass [g]")

biomass.plot

biomass.model <- lm(dryweight~ pressure * D, data = biomass)
#anova(biomass.model)

#clean up
rm(biomass.plot)

### Relative Water Content ----

#run models to compare sensitivity to treatment elements
model.rwc1 <- lm(RWC ~ pressure + D + pressure*D, data = biomass|>filter(drought=="nd",is.na(RWC)==F))
model.rwc2 <- lm(RWC ~ pressure + D , data = biomass|>filter(drought=="nd",is.na(RWC)==F))
model.rwc3 <- lm(RWC ~ pressure, data = biomass|>filter(drought=="nd",is.na(RWC)==F))
model.rwc4 <- lm(RWC ~ D, data = biomass|>filter(drought=="nd",is.na(RWC)==F))

anova(model.rwc3, model.rwc2, model.rwc1)
anova(model.rwc4, model.rwc2, model.rwc1)

anova(model.rwc3)
anova(model.rwc1)
 

# The *pressure* signal for relative water content is very strong, indicating that we can aggregate for pressure treatments.

 plot_treatment_comparison(biomass.comp|>filter(drought=="nd",is.na(RWC)==F), 
                       biomass|>filter(drought=="nd",is.na(RWC)==F), 
                       xvals = "pressure", 
                       yvals = "RWC", 
                       group = "pressure", 
                       colscale = "cols2", 
                       type = "boxplot", 
                       y.axis = "Relative Water Content" )


#At lower pressure we have a significantly lower RWC, giving us long-term image of the impact of our treatments. This is biometric support of our interpretation thus far that at low pressure *Arabidopsis thaliana* reduces stomatal conductance.

#clean up
rm(model.rwc1, model.rwc2, model.rwc3, model.rwc4)

### Leaf difference ----

#looking at the *change* in leaf number over time with comparison columns
leaf.diff.comp <- biometric.comp  |> 
  filter(drought != "d")  |>  
  reframe(.by=c(time, ID, comparison), 
          mean = mean(LeafNum), 
          treatment = treatment) |> 
  unique() |> 
  pivot_wider(values_from = mean, names_from = time)  |>  
  mutate(diff = T2-T1)  |>  
  filter(is.na(diff)==F)

#looking at the raw *change* in leaf number
leaf.diff <- biometric |>
  filter(drought != "d") |> 
  reframe(.by=c(time, ID, pressure, D), 
          mean = mean(LeafNum), 
          treatment = treatment) |> 
  pivot_wider(values_from = mean, names_from = time) |> 
  mutate(diff = T2-T1) |> 
  filter(is.na(diff)==F)


#is it normally distributed?
shapiro.test(leaf.diff$diff)
#p-value <<< 0.05, non normal - non parametric tests here we come

#plot
leaf.diff.plot <- plot_treatment_comparison(df = leaf.diff.comp,
                                                 comp = leaf.diff,
                                                 xvals = "treatment",  yvals= "diff",
                                                 group = "treatment",
                                                 colscale = "cols4",
                                                 facets = "comparison",
                                                 title = "Change in leaves per rosette",
                                                 x.axis = NULL, y.axis = "Number of Leaves",
                                                 type = "boxplot")

leaf.diff.plot

#clean up
rm(leaf.diff.plot, leaf.diff.comp)

#model to understand sensitivity to treatments
leaf.diff.model <- lm(diff ~ pressure * D, data = leaf.diff)

#NOT RUN - but show that above model is the one with the best fit
# leaf.diff.model1 <- lm(diff~ pressure, data = leaf.diff)
# leaf.diff.model2 <- lm(diff ~ D, data = leaf.diff)
# leaf.diff.model3 <- lm(diff~ pressure + D, data = leaf.diff)
# 
# anova(leaf.diff.model, leaf.diff.model1)
# anova(leaf.diff.model, leaf.diff.model2, leaf.diff.model3)
#
#clean up
# rm(leaf.diff.model1, leaf.diff.model2, leaf.diff.model3)

#------------- END ---------------------
