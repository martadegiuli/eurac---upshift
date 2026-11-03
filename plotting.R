#' Marta De Giuli
#' Eurac Research
#' UPSHIFT
#' Plotting
#' 
#' 26.10.26
#' 
#' Purpose: Generate plots and data tables for article.

## Load Libraries + functions ----
library(tidyverse)
library(ggplot2)
library(patchwork)
library(grid) 
library(gridExtra) 
library(ggpubr)

#import functions:
source("functions.R")

## Set presets ----

#color schemes
cols4 <- c("26D10" = "#fde725",  
           "26D20" = "#35b779", 
           "40D10" = "#31688e", 
           "40D20" = "#440154")

cols2 <- c("98" = "#35b779", 
           "61" = "#440154")

cols2a <- c("10" = "#932667",
            "20" = "#fca50a" )

##Plots ----

### Figure 1: Chlorophyll fluorescence responses to pressure----

#psii plot, with no x-axis label because it will be above
psii.agg.plot.paper <- plot_treatment_comparison(psii.agg, 
                                   xvals = "irr", 
                                   yvals = "QY", 
                                   group = "pressure", 
                                   facets = c("time"), 
                                   type = "line", 
                                   pvalues = df_pvalues.psii, 
                                   y.axis = expression(phi~PSII), 
                                   x.axis = NULL, 
                                   colscale = "cols2", 
                                   base_size = 20)

#npq plot, with x-axis label because it will be positioned below
npq.agg.plot.paper <- plot_treatment_comparison(npq.agg, 
                                  xvals = "irr", 
                                  yvals = "NPQ",
                                  group = "pressure", 
                                  facets = c("time"),
                                  type = "line", 
                                  pvalues = df_pvalues.npq,
                                  y.axis = "NPQ",
                                  x.axis = expression(Light~"["*"\u03BC"*E*"]"),
                                  colscale = "cols2", 
                                  base_size = 20)

#join 2 plots above
chlfluor_paper <- ggarrange(psii.agg.plot.paper,
                           npq.agg.plot.paper,
                           common.legend =T, 
                           legend = "bottom", 
                           ncol = 1, 
                           align = "hv")

#save figure
ggsave(
  filename = "figures/light_curve_paper.png",
  plot = chlfluor_paper,
  width = 10,
  height = 8.7,
  units = "in",
  dpi = 300,
  bg = "white"
)

rm(psii.agg.plot.paper, npq.agg.plot.paper)

### Figure 2: Gas exchange responses to pressure and humidity----

#stomatal conductance plot
plot.GH2O <- plot_interaction_ci(
  data = walz_f,
  model = model.GH2O,
  y_var = "GH2O",
  #x_label = "Pressure [kPa]",
  y_label = expression(gH[2]*O),
  legend = "D",
  colscale = "cols2a",
  text_size = 20,
  line_size = 2,
  show_slope = F
)

#ci:ca plot
plot.cica <-plot_interaction_ci(
  data = walz_f,
  model = model.cica,
  y_var = "cica",
  #x_label = "Pressure [kPa]",
  y_label = "Ci:Ca",
  legend = "D",
  colscale = "cols2a",
  text_size = 20,
  line_size = 2,
  show_slope = F
)

#transpiration plot
plot.E <- plot_interaction_ci(
  data = walz_f,
  model = model.E,
  y_var = "E",
  y_label = "Transpiration",
  x_label = "Pressure [kPa]",
  legend = "D",
  colscale = "cols2a",
  text_size = 20,
  line_size = 2,
  show_slope = F
)

gasexchange_paper = ggarrange(plot.GH2O, 
                              plot.cica,
                              plot.E, 
                              common.legend =T, 
                              legend = "right", 
                              ncol = 1, align = "hv")

ggsave(filename = "figures/gasexchange_paper.png",
       plot = gasexchange_paper,
       width = 6, 
       height = 8.7,
       units = "in",
       dpi = 300, 
       bg = "white")

rm(plot.GH2O, plot.cica, plot.E)

### Figure 3: Pigment responses to pressure, humidity, and time ----

#chlorophyll to carotenoid ratio
chlcar_time_paper <- ggplot(pigments, aes(x= treatment, y = chlcar, fill = time))+
  geom_boxplot()+
  stat_compare_means(label = "p.signif")+
  labs(y = "Chlorophyll : Carotenoid",
       x = NULL,
       fill = NULL)+
  theme_classic(
    base_size = 20,
    base_line_size = 1)

#chlorophyll density
chlden_time_paper <- ggplot(pigments, aes(x= treatment, y = density, fill = time))+
  geom_boxplot()+
  stat_compare_means(label = "p.signif")+
  labs(y = expression(Chlorophyll~Density~"["*"\u03BC"*g~chl~cm^-2*"]"),
       x = NULL,
       fill = NULL)+
  theme_classic(
    base_size = 20,
    base_line_size = 1
  )

pigments_paper <- ggarrange(chlcar_time_paper,
                            chlden_time_paper, 
                            nrow = 2, 
                            legend = "right", 
                            common.legend = T, 
                            align = "hv")

#save image
ggsave(filename = "figures/pigments_paper.png",
       plot = pigments_paper,
       width = 8, 
       height = 9.2,
       units = "in",
       dpi = 300, 
       bg = "white")

rm(chlcar_time_paper, chlden_time_paper)

### Figure 4: 13C discrimination normalized for d13C of chamber air----

#box plot of all treatments
d13C_plot_paper <- plot_treatment_comparison(cn.comp|>filter(drought=="nd",is.na(cn)==F), 
                                             cn|>filter(drought=="nd",is.na(cn)==F), 
                                             xvals = "treatment",
                                             yvals = "C_discr", 
                                             group = "treatment",
                                             colscale = "cols4",
                                             facets = "comparison", 
                                             type = "boxplot",
                                             y.axis = expression(delta^{13}*C) ,
                                             base_size = 20,
                                             point_size = 1.5,
                                             add_p_signif = T,
                                             add_p_format = F)

#scatter plot with interaction
d13C_inter_plot_paper <- plot_interaction_ci(
  data = cn_mod,
  model = model.d13c1,
  y_var = "C_discr",
  x_var = "pressure_n",
  group_var = "D_f",
  y_label = expression(delta^{13}*C),
  x_label = "Pressure [kPa]",
  legend = "D",
  colscale = "cols2a",
  text_size = 20,
  line_size = 2,
  show_slope = F
)

#join two plots
d13C_paper <- ggarrange(d13C_plot_paper,
                        d13C_inter_plot_paper, 
                        nrow = 2, 
                        legend = "right", 
                        common.legend = F)

#save image
ggsave(filename = "figures/d13C_paper.png",
       plot = d13C_paper,
       width = 9.5, 
       height = 6,
       units = "in",
       dpi = 300, 
       bg = "white")

rm(d13C_plot_paper, d13C_inter_plot_paper)

### Supplementary Figure 1: Response of 1-qP to air pressure----

qP.agg.plot.paper <- plot_treatment_comparison(qP.agg, 
                                 xvals = "irr", 
                                 yvals = "qP",
                                 group = "pressure", 
                                 facets = c("time"),
                                 type = "line",
                                 pvalues = df_pvalues.qP,
                                 y.axis = "1- qP", 
                                 x.axis = expression(Light~"["*"\u03BC"*E*"]"), 
                                 colscale = "cols2", 
                                 base_size = 20)

ggsave(
  filename = "figures/qP_paper.png",
  plot = qP.agg.plot.paper,
  width = 10,
  height = 5,
  units = "in",
  dpi = 300,
  bg = "white"
)

### Supplementary Figure 2: Response of C:N ratio----
cn_paper <- plot_treatment_comparison(cn.comp|>filter(drought=="nd",is.na(cn)==F),
           cn|>filter(drought=="nd",is.na(cn)==F), 
           xvals = "treatment", 
           yvals = "cn", 
           group = "treatment",
           colscale = "cols4", 
           facets = "comparison", 
           type = "boxplot",
           y.axis = "C:N",
           base_size = 20,
           point_size = 1.5,
           add_p_signif = T,
           add_p_format = F)

ggsave(filename = "figures/cn_paper.png",
       plot = cn_paper,
       width = 9.5, 
       height = 6,
       units = "in",
       dpi = 300, 
       bg = "white")

##Tables----

### Table 1: Summary of statistical analysis of the effect of D, air pressure, their interaction, and time on ecophysiological variables----

#' Format
#'             | D        | pressure    | D:pressure | time  |
#' Plant Trait | F        | F           | F          | F     |
#' A           | ns       | 42.14***    | ns         | --    |
#' Chl. : Car	 |8.387	**	|	ns		      | ns	       |5.99	*|
#' Ci:Ca       | ns       | 9.24**      | 6.734*     | --    |

#initialize table
model.table <- data.frame(trait = as.character(),
                          D.F = as.numeric(),
                          D.p = as.numeric(),
                          "D.p*" = as.character(),
                          P.F = as.numeric(),
                          P.p = as.numeric(),
                          "P.p*" = as.character(),
                          "D*P.F" = as.numeric(),
                          "D*P.p" = as.numeric(),
                          "D*P.p*" = as.character(),
                          time.F = as.numeric(),
                          time.p = as.numeric(),
                          "time.p*" = as.character()
)

#pigments
for(response in c("chlab", "chlcar", "density")){
  
  mod <- lm(get(response) ~ pressure*D + time,
            data = pigments)
  anova_df <- as.data.frame(anova(mod))
  model.table <- 
    rbind(
      model.table,
      data.frame(
        trait = response,
        D.F = anova_df["D", "F value"],
        D.p = anova_df["D", "Pr(>F)"],
        "D.p*" = add_stars(anova_df["D", "Pr(>F)"]),
        P.F = anova_df["pressure", "F value"],
        P.p = anova_df["pressure", "Pr(>F)"],
        "P.p*" = add_stars(anova_df["pressure", "Pr(>F)"]),
        "D*P.F" = anova_df["pressure:D", "F value"],
        "D*P.p" = anova_df["pressure:D", "Pr(>F)"],
        "D*P.p*" = add_stars(anova_df["pressure:D", "Pr(>F)"]),
        time.F = anova_df["time", "F value"],
        time.p = anova_df["time", "Pr(>F)"],
        "time.p*" = add_stars(anova_df["time", "Pr(>F)"])
      )
    )
}

#walz
for(response in c("A", "cica", "GH2O", "E")){
  mod <- lm(get(response) ~ pressure*D,
            data = walz|> filter(Drought == F))
  anova_df <- as.data.frame(anova(mod))
  model.table <- 
    rbind(
      model.table,
      data.frame(
        trait = response,
        D.F = anova_df["D", "F value"],
        D.p = anova_df["D", "Pr(>F)"],
        "D.p*" = add_stars(anova_df["D", "Pr(>F)"]),
        P.F = anova_df["pressure", "F value"],
        P.p = anova_df["pressure", "Pr(>F)"],
        "P.p*" = add_stars(anova_df["pressure", "Pr(>F)"]),
        "D*P.F" = anova_df["pressure:D", "F value"],
        "D*P.p" = anova_df["pressure:D", "Pr(>F)"],
        "D*P.p*" = add_stars(anova_df["pressure:D", "Pr(>F)"]),
        time.F = NA,
        time.p = NA,
        "time.p*" = NA)
    )
  
}

#biomass
mod <- lm(dryweight ~pressure*D,
          data = biomass |> filter(drought == "nd"))
anova_df <- as.data.frame(anova(mod))
model.table <- 
  rbind(
    model.table,
    data.frame(
      trait = "AGB",
      D.F = anova_df["D", "F value"],
      D.p = anova_df["D", "Pr(>F)"],
      "D.p*" = add_stars(anova_df["D", "Pr(>F)"]),
      P.F = anova_df["pressure", "F value"],
      P.p = anova_df["pressure", "Pr(>F)"],
      "P.p*" = add_stars(anova_df["pressure", "Pr(>F)"]),
      "D*P.F" = anova_df["pressure:D", "F value"],
      "D*P.p" = anova_df["pressure:D", "Pr(>F)"],
      "D*P.p*" = add_stars(anova_df["pressure:D", "Pr(>F)"]),
      time.F = NA,
      time.p = NA,
      "time.p*" = NA)
  )

#rwc
mod <- lm(RWC ~ pressure*D,
          data = biomass|> filter(drought == "nd"))
anova_df <- as.data.frame(anova(mod))
model.table <- 
  rbind(
    model.table,
    data.frame(
      trait = "RWC",
      D.F = anova_df["D", "F value"],
      D.p = anova_df["D", "Pr(>F)"],
      "D.p*" = add_stars(anova_df["D", "Pr(>F)"]),
      P.F = anova_df["pressure", "F value"],
      P.p = anova_df["pressure", "Pr(>F)"],
      "P.p*" = add_stars(anova_df["pressure", "Pr(>F)"]),
      "D*P.F" = anova_df["pressure:D", "F value"],
      "D*P.p" = anova_df["pressure:D", "Pr(>F)"],
      "D*P.p*" = add_stars(anova_df["pressure:D", "Pr(>F)"]),
      time.F = NA,
      time.p = NA,
      "time.p*" = NA)
  )

#leafdiff
mod <- lm(diff ~ pressure*D,
          data = leaf.diff)
anova_df <- as.data.frame(anova(mod))
model.table <- 
  rbind(
    model.table,
    data.frame(
      trait = "Leaf Difference",
      D.F = anova_df["D", "F value"],
      D.p = anova_df["D", "Pr(>F)"],
      "D.p*" = add_stars(anova_df["D", "Pr(>F)"]),
      P.F = anova_df["pressure", "F value"],
      P.p = anova_df["pressure", "Pr(>F)"],
      "P.p*" = add_stars(anova_df["pressure", "Pr(>F)"]),
      "D*P.F" = anova_df["pressure:D", "F value"],
      "D*P.p" = anova_df["pressure:D", "Pr(>F)"],
      "D*P.p*" = add_stars(anova_df["pressure:D", "Pr(>F)"]),
      time.F = NA,
      time.p = NA,
      "time.p*" = NA)
  )

#cn
for(response in c("cn", "C_discr", "N", "C")){
  mod <- lm(get(response) ~ pressure*D,
            data = cn|> filter(drought == "nd"))
  anova_df <- as.data.frame(anova(mod))
  model.table <- 
    rbind(
      model.table,
      data.frame(
        trait = response,
        D.F = anova_df["D", "F value"],
        D.p = anova_df["D", "Pr(>F)"],
        "D.p*" = add_stars(anova_df["D", "Pr(>F)"]),
        P.F = anova_df["pressure", "F value"],
        P.p = anova_df["pressure", "Pr(>F)"],
        "P.p*" = add_stars(anova_df["pressure", "Pr(>F)"]),
        "D*P.F" = anova_df["pressure:D", "F value"],
        "D*P.p" = anova_df["pressure:D", "Pr(>F)"],
        "D*P.p*" = add_stars(anova_df["pressure:D", "Pr(>F)"]),
        time.F = NA,
        time.p = NA,
        "time.p*" = NA)
    )
  
}

#chlor fluor
for(response in c("QY", "NPQ", "ETR", "qE", "qP")){
  
  subset.pam <- pam |>
    filter(light == "Lss")
  
  mod <- lm(get(response) ~ irr + pressure*D + time,
            data = subset.pam)
  
  anova_df <- as.data.frame(anova(mod))
  model.table <- 
    rbind(
      model.table,
      data.frame(
        trait = response,
        D.F = anova_df["D", "F value"],
        D.p = anova_df["D", "Pr(>F)"],
        "D.p*" = add_stars(anova_df["D", "Pr(>F)"]),
        P.F = anova_df["pressure", "F value"],
        P.p = anova_df["pressure", "Pr(>F)"],
        "P.p*" = add_stars(anova_df["pressure", "Pr(>F)"]),
        "D*P.F" = anova_df["pressure:D", "F value"],
        "D*P.p" = anova_df["pressure:D", "Pr(>F)"],
        "D*P.p*" = add_stars(anova_df["pressure:D", "Pr(>F)"]),
        time.F = anova_df["time", "F value"],
        time.p = anova_df["time", "Pr(>F)"],
        "time.p*" = add_stars(anova_df["time", "Pr(>F)"])
      )
    )
}

write.csv(model.table |>
            dplyr::mutate(
              dplyr::across(ends_with(c("F","p")), 
                            ~ round(.x, 3))), 
          "models.csv", 
          row.names = F)

### Supplementary Table 1: Summary of variables measured with mean values +/- standard deviations ----

#'Format
#'      |Treatments                                         |
#'Trait | 26D10      | 26D20      | 40D10      | 40D20      |
#' RWC  | mean +/- sd| mean +/- sd| mean +/- sd| mean +/- sd|
#' 
#' For ease of use - first make long table:
#' RWC - treat1
#' RWC - treat2
#' ....
#' 
#' then rotate to wide format


#pigments
pigment.means <- pigments|>
  pivot_longer(c(chlab, chlcar, density), names_to = "trait")|>
  dplyr::summarise(.by = c("trait", "treatment", "time"),
                   mean = mean(value),
                   sd = sd(value),
                   n = n())|>
  arrange(trait, time)

#walz
walz.means <- walz |>
  filter(Drought == F) |>
  pivot_longer(c("A", "cica", "GH2O", "E"), names_to = "trait")|>
  dplyr::summarise(.by = c("trait", "treatment"),
                   mean = mean(value, na.rm = T),
                   sd = sd(value, na.rm = T),
                   n = n())

#agb
biomass.means <- biomass |> 
  filter(drought == "nd") |> 
  pivot_longer(dryweight, names_to = "trait") |> 
  dplyr::summarise(.by =c("trait", "treatment"),
                   mean = mean(value, na.rm = T),
                   sd = sd(value, na.rm = T),
                   n = n())

#rwc
rwc.means <- biomass |>
  filter(drought == "nd") |>
  pivot_longer(RWC, names_to = "trait")|>
  dplyr::summarise(.by = c("trait", "treatment"),
                   mean = mean(value, na.rm = T),
                   sd = sd(value, na.rm = T),
                   n = n())

#leaf.diff
leaf.diff.means <- leaf.diff |>
  pivot_longer(diff, names_to = "trait")|>
  dplyr::summarise(.by = c("trait", "treatment"),
                   mean = mean(value, na.rm = T),
                   sd = sd(value, na.rm = T),
                   n= n())

#cn
cn.means <- cn |>
  filter(drought == "nd")|>
  pivot_longer(c("cn", "C_discr", "N", "C"), names_to = "trait")|>
  dplyr::summarise(.by = c("trait", "treatment"),
                   mean = mean(value, na.rm = T),
                   sd = sd(value, na.rm = T),
                   n= n())

#chlorfluor
pam.means <- pam |>
  filter(light == "Lss",
         drought == "nd")|>
  pivot_longer(c("QY", "NPQ", "ETR", "qE", "qP"), names_to = "trait")|>
  dplyr::summarise(.by = c("trait", "treatment", "time", "irr"),
                   mean = mean(value, na.rm = T),
                   sd = sd(value, na.rm = T),
                   n = n())|>
  arrange(trait, treatment, time, irr)

value.table <- bind_rows(pigment.means,
                         walz.means,
                         biomass.means,
                         rwc.means,
                         leaf.diff.means,
                         cn.means,
                         pam.means)|>
  dplyr::mutate(value = sprintf("%.3f ± %.3f (%i)", mean, sd, n), .keep = "unused")|>
  pivot_wider(names_from = treatment, values_from = value)

#clean up
rm(pigment.means, walz.means, biomass.means, rwc.means, leaf.diff.means, cn.means, pam.means)

write.csv(value.table, "values.csv", row.names = F) 

#------------- END ---------------------