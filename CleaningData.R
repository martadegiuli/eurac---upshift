library(readxl)
library(tidyverse)

#clean co2 data
rotations <- read_excel("../RAW Data/rotations.xlsx")

co2 <- read.delim("../RAW Data/CO2.txt", header = TRUE, sep = "\t") |>
  mutate(date = parse_date_time(X, c("%d.%m.%y %H:%M:%S")), 
         .keep = "unused")|>
  pivot_longer(cols = c("SEC1", "SEC2", "SEC3", "SEC4"),
               names_to = "SEC",
               values_to = "co2_ppm")

co2_treat <-  co2 |>
  left_join(
    rotations,
    join_by(
      SEC,
      date >= date        # co2 timestamp >= treatment start time
    )
  ) |>
  group_by(SEC, date.x)|>
  slice_max(date.y, with_ties = FALSE) |>
  ungroup() |>
  dplyr::rename(date = date.x)|>
  dplyr::select(-date.y)

write.csv(co2_treat, "../RAW Data/co2_treat.csv")


#import data
chlab   <- read_excel("../RAW Data/chl_fluorcam_2025.xlsx", sheet = "chl ab")
chlcar  <- read_excel("../RAW Data/chl_fluorcam_2025.xlsx", sheet = "chlcar")
chlden  <- read_excel("../RAW Data/chl_fluorcam_2025.xlsx", sheet = "ug chl per cm2")

chlab<-chlab|>
  rownames_to_column(var = "sample.id")|>
  pivot_longer(2:9,
               names_to = "treatment",
               values_to = "chlab")


chlcar<-chlcar|>
  rownames_to_column(var = "sample.id")|>
  pivot_longer(2:9,
               names_to = "treatment",
               values_to = "chlcar")

chlden<-chlden|>  
  rownames_to_column(var = "sample.id")|>
  pivot_longer(2:9,
               names_to = "treatment",
               values_to = "density")


pigments <- chlab |>
  left_join(chlcar)|>
  left_join(chlden)

write.csv(pigments, "../RAW Data/pigments.csv")


#clean PAM data

#create list that enumerates each file path that needs to be imported:

files <- c("../RAW Data/NPQ 4000 T1 8.xlsx",
           "../RAW Data/NPQ 4000 T1 18.xlsx",
           "../RAW Data/NPQ 4000 T1 30.xlsx",
           "../RAW Data/NPQ 4000 T1 60.xlsx",
           "../RAW Data/NPQ 4000 T1 90.xlsx",
           "../RAW Data/NPQ 4000 T2 8.xlsx",
           "../RAW Data/NPQ 4000 T2 19.xlsx",
           "../RAW Data/NPQ 4000 T2 30.xlsx",
           "../RAW Data/NPQ 4000 T2 60.xlsx",
           "../RAW Data/NPQ 4000 T2 90.xlsx",
           "../RAW Data/NPQ_260_T1 8.xlsx",
           "../RAW Data/NPQ_260_T1 18.xlsx",
           "../RAW Data/NPQ_260_T1 30.xlsx",
           "../RAW Data/NPQ_260_T1 60.xlsx",
           "../RAW Data/NPQ_260_T1 90.xlsx",
           "../RAW Data/NPQ_260_T2 8.xlsx",
           "../RAW Data/NPQ_260_T2 18.xlsx",
           "../RAW Data/NPQ_260_T2 30.xlsx",
           "../RAW Data/NPQ_260_T2 60.xlsx",
           "../RAW Data/NPQ_260_T2 90.xlsx"
           )

# create list to store everything

data.list <- list()

#import data

for(file in files){
  #extract name of file
  name <- substr(file, 17, nchar(file)-5)
  #extrac light intensity of PAM run
  light_intensity <- substr(file, regexpr(" [^ ]*$", file)+1, nchar(file)-5)
  
  #import file
  test <- read_excel(file,
                     range = cell_rows(5:142))
  
    #extract constants to be able to pivot table
  constants <- test |> 
    slice_head(n=8) |>
    pivot_longer(
      cols = starts_with(c("NPQ_", "PAM_")),
      names_to = "sample",
      values_to = "measure"
    )|>
    pivot_wider(
      names_from = `...1`,
      values_from = measure
    ) |>
    dplyr::rename(Fm_const = Fm,
           Fo_const = Fo,
           Fv_const = Fv,
           Fp_const = Fp)
  
  QY_max <- test |>
    filter(`...1` == "QY_max")|>
    pivot_longer(
      cols = starts_with(c("NPQ_", "PAM_")),
      names_to = "sample",
      values_to = "measure"
    )|>
    pivot_wider(
      names_from = `...1`,
      values_from = measure
    )
  
  data <- test |>
    filter(!`...1` == "QY_max")|>
    slice(-(1:8))|>
    separate_wider_delim(
      `...1`,
      names = c("measure", "light"),
      delim = "_"
    ) |>
    pivot_longer(
      cols = starts_with(c("NPQ_", "PAM_")),
      names_to = "sample",
      values_to = "values"
    )|>
    pivot_wider(
      names_from = measure,
      values_from = values
    )|>
    left_join(constants)|>
    left_join(QY_max)|>
    mutate(light_intensity = as.numeric(light_intensity)) |>
    mutate(irr = case_when(
      light_intensity == 8 ~ 150,
      light_intensity == 18 | light_intensity == 19 ~ 300,
      light_intensity == 30 ~ 600,
      light_intensity == 60 ~ 1200,
      light_intensity == 90 ~ 2000)
    ) |>
    mutate(
      NPQ1 = (Fm_const - Fm)/Fm,
      ETR = QY * 0.5 * irr)
    
  qE <- data |>
    filter(light == "Lss"| light == "D4")|>
    dplyr::select(sample, NPQ1, light) |>
    pivot_wider(names_from = light,
                values_from = NPQ1) |>
    dplyr::mutate(qE = Lss - D4, .keep = "unused")

  data2 <- data |>
    left_join(qE)
  
L0_rows <- data2 |>
  distinct(sample, light_intensity) |>
  mutate(light = "L0") |>
  left_join(data2 |>
              slice(1) |> 
              dplyr::select(-sample, -light, -names(constants), QY_max, qE) |>
              dplyr::mutate(QY = QY_max,
                            across(where(is.numeric) &!c(light_intensity, QY, QY_max, irr), ~ 0)
                            ))|>
  left_join(constants)

data3 <- bind_rows(data2, L0_rows)  

  data.list[[name]]<- data3
}

#adds irradiance = 0 values for
irr0 <- bind_rows(data.list)|> 
  dplyr::summarise(.by = sample, QY = mean(QY_max), light = "Lss", `Size [pixels]` = mean(`Size [pixels]`), across(where(is.double) & !c(light, QY, `Size [pixels]`), ~0))

PAM <- bind_rows(data.list)|>
  bind_rows(irr0)|>
  mutate(treatment = substr(sample, 5, 9),
         time = substr(sample, 11,12),
         ID = substr(sample, 14, nchar(sample)))

write.csv(PAM, "../RAW Data/PAM.csv", row.names = F)
