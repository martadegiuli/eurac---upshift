#'function to add comparison column to data tables.
#'agruments:
#'1. `data` - a dataset that needs to be "comparisonized"
#'2. `pressure` - the name of a column of said dataset (can be character/numeric) which contains the pressure
#'3. `D` - the name of a column of said dataset (can be character/numeric) which contains the D value
add_comp = function(data, pressure, D){
  press <-
    data |>
    mutate(comparison = paste0("P = ", as.character(pressure)))
  
  d <-
    data |>
    mutate(comparison = paste0("D = ", as.character(D)))
  
  comp <-
    press|>
    add_row(d)
  
  return(comp)
}

#'function to create base faceted plot.
#'arguments:
#'1. df  - a dataset to be plotted
#'2. comp  - a dataset to use as a background comparison
#'3. xvals - the name of the variable to be plotted on the x axis
#'4. yvals - the name of the variable to be plotted on the y axis
#'5. group  - the name of the variable to fill/color the plot by
#'5. colscale - the name of the color palette to use
#'6. facets- a 2 element list which will be the variables used to facet the plot. If only one facet is wanted use "horizontal" or "vertical" as the first element to indicate the direction of the facet and the name of the variable in the second. Remember the first element is the vertical facet, the second the horizontal one.
#'7. title - a string with the title (if none added, will result in a plot with no title)
#'8. x.axis- a string with the x.axis label (if none added, will result in a plot with no x axis label)
#'9. y.axis- a string with the y.axis label (if none added, will result in a plot with no y axis label)
#'10. type - a string either "boxplot" or "scatter" which will indicate whether to make a boxplot or a scatter plot - with lines (because that's what we are using here)
#'11. pvalues- a df with the pvalues to add to the scatter plot (i.e. the chlorophyll fluorescence plots) - make sure x axis names are the same because they get reused!

base_plot <- function(df, comp, xvals, yvals, group = NULL, colscale , facets, title = NULL, x.axis = NULL, y.axis = NULL, type = NULL, pvalues){
  
  #generate y lims (.lim = limit of plot, .lev = positioning of significance values)
  y.lim = max(df[[yvals]]) * 1.2
  y.lev = max(df[[yvals]]) * 1.1
  
  #generate color palette if colscale is missing:
  
  if (is.null(colscale) && !is.null(group)) {
    n_groups <- length(unique(df[[group]]))
    colscale <- scales::hue_pal()(n_groups)
  }
  
  
  if(shapiro.test(comp[[yvals]])$p.value > 0.05){
    stat_test = "t.test" # if the distribution is normal, we can use regular parametric tests
  }else{
    stat_test = "wilcox.test"  # otherwise we need to use the non-parametric wilcox test.
  }
  
  if(facets[1] != "horizontal" && facets[1] != "vertical"){
    facet_layer =  facet_grid(as.formula(paste(facets[1],"~", facets[2])), drop = TRUE)
  }else{
    if(facets[1] == "horizontal"){
      facet_layer = facet_grid(as.formula(paste("~", facets[2])), drop = TRUE)
    }else{
      if(facets[1] == "vertical"){
        facet_layer = facet_grid(as.formula(paste(facets[2],"~")), drop = TRUE)
      }
    }
  }
  
  if(type == "boxplot"){
  plot <- ggplot(data = df, aes(x = get(xvals), y = get(yvals), fill = get(group)))+
    geom_boxplot(data = comp, color = "grey",fill = "grey")+ 
    geom_boxplot()+
    geom_jitter(width=0.01)+
    labs(title = title,
         x = x.axis,
         y = y.axis)+
    scale_fill_manual(values= get(colscale),
                      name = group)+
    stat_compare_means(method = stat_test, 
                       label = "p.signif",
                       label.y = y.lev,
                       label.x = c(2, 3, 3.5, 1, 2, 3, 3.5, 1)
    )+
    stat_compare_means(method = stat_test, 
                       label = "p.format",
                       label.y.npc = "bottom"
    )+
    facet_layer+
    ylim(NA, y.lim)+
    theme_bw( )+
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(), 
          legend.position = "bottom")
  
  return(plot)
  }else{
    
    #generating y placements for significance values
    #set an offset that makes sense (arbitrarily 1/10 of the max value?)
    
    offset <- max(comp[[yvals]])*.1
    
    ymax <- dplyr::summarise(comp, .by = .data[[xvals]], max = max(.data[[yvals]])+offset)|> pull(max)
    
    plot <-   ggplot(data = df, aes(x = get(xvals), y = get(yvals), colour = get(group)))+
      geom_line()+
      geom_errorbar(data = comp, 
                    color = "grey", 
                    aes(ymin = get(yvals)-ci, ymax = get(yvals)+ci), 
                    width = 0.1)+ 
      geom_point(data = comp, 
                 color = "grey", 
                 aes(x = get(xvals), y = get(yvals)))+
      geom_errorbar(aes(ymin = get(yvals)-ci, ymax = get(yvals)+ci), 
                    width = 0.1)+ 
      geom_point()+
      geom_text(data=pvalues,
                aes(x = get(xvals),
                    y = rep(c(ymax),8),
                    label = sig),
                colour = "black")+
      labs(title = title,
           x = x.axis,
           y = y.axis)+
      scale_color_manual(values= get(colscale),
                         name = group)+
      facet_layer+
      theme_bw( )+
      theme(legend.position = "bottom")
    
    return(plot)
  }
}

#'function to calculate significance
#'arguments:
#'1. data - dataframe from which to calculate
#'2. comparison_column - column from which to extract pairings to make
#'3. variable - value to compare
#'
signif<- function (data, comparison_column, variable){
  
  #initialize dataframe
  df_pvalues <- data.frame(comp1 = as.character(),
                         comp2 = as.character(),
                         irr = as.numeric(),
                         time = as.character(),
                         comparison = as.character(),
                         p.value = as.numeric(),
                         sig = as.character(),
                         n = as.numeric()
)
  
for(sampling_time in unique(data$time)){ # for each time stamp
  subset <- data|>
    filter(time == sampling_time) #select time
  
  for (comp in unique(data[[comparison_column]])){ #for each category
    subset1 <- subset|>
      filter(get(comparison_column) == comp) #select category
    
    for(rad in unique(data$irr)){ # for each level within
      subset2 <- subset1|>
        filter(irr == rad) #select level - here irradiance (hard coded)
      
      group1 = pull(subset2|>
                      filter(treatment == unique(subset2$treatment)[1])|>  #within each level, we have 2 treatments - pulls 1st
                      dplyr::select(variable))
      
      group2 = pull(subset2|>
                      filter(treatment == unique(subset2$treatment)[2])|> #within each level, we have 2 treatments - pulls 2nd
                      dplyr::select(variable))
      
      
      #decide on statistical test to use (shapiro.test() only works for up to 5000 values)
      if(shapiro.test(c(group1,group2))$p.value > 0.05){
        stat_test = t.test(group1, group2) # if the distribution is normal, we can use regular parametric tests
      }else{
        stat_test = wilcox.test(group1,group2)  # otherwise we need to use the non-parametric wilcox test.
      }
      
      #build dataframe
      df_pvalues<-rbind(df_pvalues, 
                        data.frame(comp1 = unique(subset2$treatment)[1],
                                   comp2 = unique(subset2$treatment)[2],
                                   irr = rad,
                                   time = sampling_time,
                                   comparison = comp,
                                   p.value = stat_test$p.value,
                                   sig = case_when(
                                     stat_test$p.value < 0.001 ~ "***",
                                     stat_test$p.value < 0.01  ~ "**",
                                     stat_test$p.value < 0.05  ~ "*",
                                     stat_test$p.value > 0.05  ~ " "
                                   ),
                                   n = (length(group1)+length(group2))/2 #mean sample size
                        )
      )
    }
  }
}
  return(df_pvalues)
}


#'functions to add n to plots:
n_fun_below <- function(x){
return(data.frame(y = -Inf , label = length(x)))
}

n_fun_bottwisk <- function(x){
  return(data.frame(y = median(x) - (1.8 *IQR(x)) , label = length(x)))
}


summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

# Self-starting function for the exponential decay model
SSexpDecay <- selfStart(
  model = ~ Pmin + Pmax * exp(-alpha * I / Pmax),
  initial = function(mCall, data, LHS, ...) {
    # Extract the independent variable and the response variable
    xy <- sortedXyData(mCall[["I"]], LHS, data)
    
    # Ensure there are enough data points
    if (nrow(xy) < 4) stop("bro you don't have enough values")
    
    #extract indep + response variables
    I <- xy[["x"]]
    y <- xy[["y"]]
    
    # Initial guess for Pmax is the maximum observed value of y
    Pmax_init <- max(y)
    Pmin_init <- min(y)
    
    # Fit a linear model to the log-transformed data to estimate alpha
    # Avoid log(0) by replacing 0s with a small positive value
    y_transformed <- ifelse(y == 0, 1e-6, y)
    log_y <- log(y_transformed)
    
    # Linear model: log(y) =log(Pmin)+ log(Pmax) - (alpha / Pmax) * I
    fit <- lm(log_y ~ I)
    alpha_init <- -coef(fit)[2] * Pmax_init  # Extract alpha from the slope
    
    # Return the initial parameter estimates as a named list
    setNames(c(Pmax = Pmax_init, alpha = alpha_init, Pmin = Pmin_init), 
             mCall[c("Pmax", "alpha", "Pmin")])
  },
  parameters = c("Pmax", "alpha", "Pmin")
)

# Function to create compact letter displays
create_cld <- function(subset_data) {
  # Extract group comparisons and p-values
  comparisons <- strsplit(as.character(subset_data$pairs), "-")       #extract pairs
  groups <- unique(unlist(comparisons))                               #extract unique treatments
  p_matrix <- matrix(1, nrow = length(groups), ncol = length(groups)) #create square matrix with for treatment comparisons
  #rename the rows and columns to treatment names
  rownames(p_matrix) <- groups
  colnames(p_matrix) <- groups
  
  #for each comparison (saved in 'comparisons'), save p-value in matrix)
  for (i in 1:nrow(subset_data)) {
    pair <- comparisons[[i]]
    p_matrix[pair[1], pair[2]] <- subset_data$p_value[i]
    p_matrix[pair[2], pair[1]] <- subset_data$p_value[i]
  }
  
  # Create compact letter display
  cld <- multcompLetters(p_matrix)
  return(cld$Letters)
}

# Function to check overlap between two confidence intervals
check_overlap <- function(ci1, ci2) {
  !(ci1[2] < ci2[1] | ci2[2] < ci1[1])
}

# Pairwise comparison of confidence intervals
compare_param <- function(data) {
  comparisons <- data.frame(
    pair = character(),
    overlap = logical(),
    stringsAsFactors = FALSE
  )
  
  for (i in 1:(nrow(data) - 1)) {
    for (j in (i + 1):nrow(data)) {
      pair_name <- paste(data$treatment[i], data$treatment[j], sep = "-")
      overlap <- check_overlap(
        ci1 = c(data$low[i], data$high[i]),
        ci2 = c(data$low[j], data$high[j])
      )
      comparisons <- rbind(comparisons, data.frame(pair = pair_name, overlap = overlap))
    }
  }
  
  return(comparisons)
}

# Compare coefficients of non-nested models as per Clogg et al. (1995) formula
# returns a data.frame with the diff, zdiff, and p.value of each coefficient

compare_coefs <- function(.data1, .data2){
  
  # imports map_dbl() and pluck() functions from purrr library
  import::here(map_dbl, pluck, .from = purrr)  
  
  # extract the relevant data
  b1 = map_dbl(df1[, 1], pluck)     #get reg. coefs for model 1
  se1 = map_dbl(df1[, 2], pluck)    #get std errors for model 1
  
  b2 = map_dbl(df2[, 1], pluck)     #get reg. coefs for model 2
  se2 = map_dbl(df2[, 2], pluck)    #get std. errors for model 2
  
  # Clogg et al. (1995) formula as cited by Ray Paternoster et al. (1998)
  b = b1 - b2
  s1 = se1^2
  s2 = se2^2
  sc = s1 + s2
  v = b / sqrt(sc)
  
  data.frame(diff=b, zdiff=v, `p-value`=format(2*pnorm(-abs(v)), scientific=FALSE))
  
}