#' Function to add a comparison column to a data table, so the same data can be faceted by pressure and by D
#' arguments:
#' 1. data - data frame to be "comparisonized". Must contain the pressure and D columns.
#' 2. pressure - name of the column in data (character or numeric values) which contains the pressure.
#' 3. D - name of the column in data (character or numeric values) which contains the D value.
#'
#' How it works:
#' - Makes a first copy of data with a new column "comparison" built from the pressure values, e.g. "P = 5".
#' - Makes a second copy of data with "comparison" built from the D values, e.g. "D = 20".
#' - Stacks the two copies on top of each other.
#'
#' Returns: a data frame with every row of data appearing twice (so twice as many rows as the input) plus the new column "comparison". The first half is  labelled "P = <pressure>" and the second half "D = <D>". Facetting or grouping by "comparison" therefore shows the same data once organised by pressure and once organised by D.
#'
#' @param data data frame
#' @param pressure string, name of the pressure column
#' @param D string, name of the D column
#' @return data frame with the original columns plus "comparison", with twice as many rows as data
#' @examples
#' add_comp(df, pressure = "pressure", D = "D")
#' @export

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

#' Function to create a base faceted plot (boxplot or line plot with error bars) with an optional grey background comparison layer and significance annotations
#' 
#' arguments:
#' 1. df - data frame to be plotted. For type = "line" it must contain a "ci" column (e.g. the output of summarySE()).
#' 2. comp - data frame to use as a grey background comparison. Can be left NULL (the default) if no comparison layer is wanted. For type = "line" it must also contain a "ci" column.
#' 3. xvals - name (as a string) of the variable plotted on the x axis.
#' 4. yvals - name (as a string) of the variable plotted on the y axis.
#' 5. group - name (as a string) of the variable to fill/colour the plot by. Needed in practice for both plot types.
#' 6. colscale - name (as a string) of an existing vector of colours to use, which is looked up with get(). Can be left NULL (the default), in which case ggplot's default hue palette is used, with one colour per level of group.
#' 7. facets - character vector of the variable(s) to facet by. With two elements, facets are laid out as facets[1] ~ facets[2] (rows ~ columns). With one element, facets are made on that variable only. Can be left NULL (the default) for no faceting.
#' 8. title - string with the plot title. If NULL (the default) there is no title.
#' 9. x.axis - string with the x axis label. If NULL (the default) there is no x axis label.
#' 10. y.axis - string with the y axis label. If NULL (the default) there is no y axis label.
#' 11. type - either "boxplot" or "line". "boxplot" draws boxplots with jittered points; "line" draws lines with points and error bars (mean +/- ci). Defaults to "boxplot".
#' 12. pvalues - data frame with the p-values / significance stars to add to the line plot (e.g. output of signif()). Must contain a "sig" column and the facet and xvals columns, with exactly the same names as in df, because they are used to join the annotations to the plot. Only used when type = "line". Can be left NULL (the default).
#' 13. add_p_signif - boolean, defaults to TRUE. Controls the significance annotation (*, **, ***).
#' 14. add_p_format - boolean, defaults to TRUE. Controls the formatted p-value annotation (p = xyz). Only applies to type = "boxplot".
#' 15. base_size - base font size of the plot (theme_bw), defaults to 11.
#' 16. point_size - size of the points used in the boxplot outliers, jitter and line plot points, defaults to 2.
#'
#' How it works:
#' - Calculates y limits: the top of the y axis is extended 20% beyond the data range to leave room for annotations, and the significance labels are placed at 5% above the maximum.
#' - Sets up the colour scale (custom via colscale, or default hue palette).
#' - Chooses the statistical test: runs a Shapiro-Wilk normality test on the y values (of comp if supplied, otherwise of df). If p > 0.05 a t-test is used, otherwise a Wilcoxon test.
#' - Builds the base plot (theme_bw, legend at the bottom, x/y/title labels).
#' - If comp is supplied, draws it first in grey (boxplots, or points with error bars) so it appears behind the main data.
#' - Draws the main data: for "boxplot", boxplots filled by group plus jittered points; for "line", coloured lines, error bars (y +/- ci) and points.
#' - Adds statistical annotations:
#'     * "boxplot": stat_compare_means() adds the significance stars and/or the formatted p-value, using the test chosen above.
#'     * "line": significance labels from pvalues$sig are placed above the highest value in each facet / x combination.
#' - Facets the plot (facet_grid) if facets are given.
#' - Applies the y limits with coord_cartesian() (zooms, does not drop data).
#' - For boxplots, removes the x axis text and ticks.
#'
#' Returns: a ggplot object, which can be further modified with + or printed.
#'
#'
#' @param df data frame to plot
#' @param comp data frame for the grey background layer (optional)
#' @param xvals string, x variable name
#' @param yvals string, y variable name
#' @param group string, fill/colour variable name
#' @param colscale string, name of a colour vector (optional)
#' @param facets character vector of facet variable(s) (optional)
#' @param title string, plot title (optional)
#' @param x.axis string, x axis label (optional)
#' @param y.axis string, y axis label (optional)
#' @param type "boxplot" or "line"
#' @param pvalues data frame with a "sig" column, e.g. from signif() (optional)
#' @param add_p_signif logical, add significance stars? Defaults to TRUE
#' @param add_p_format logical, add formatted p-value? Defaults to TRUE
#' @param base_size numeric, base font size. Defaults to 11
#' @param point_size numeric, point size. Defaults to 2
#' @return ggplot object
#' @examples
#' plot_treatment_comparison(df, comp = comp_df, xvals = "pressure_n",
#'                           yvals = "Fv_Fm", group = "D_f",
#'                           facets = c("comparison"), type = "boxplot")
#' plot_treatment_comparison(summary_df, xvals = "irr", yvals = "yield",
#'                           group = "treatment", facets = c("time", "species"),
#'                           type = "line", pvalues = signif_df)
#' @export

plot_treatment_comparison <- function(df,
                       comp = NULL,
                       xvals,
                       yvals,
                       group = NULL,
                       colscale = NULL,
                       facets = NULL,
                       title = NULL,
                       x.axis = NULL,
                       y.axis = NULL,
                       type = c("boxplot", "line"),
                       pvalues = NULL,
                       add_p_signif = TRUE,
                       add_p_format = TRUE,
                       base_size = 11,
                       point_size = 2) {
  
  type <- match.arg(type)
  
  ## 1. Common constants
  
  y.lim <- max(df[[yvals]], na.rm = TRUE) +
    (max(df[[yvals]], na.rm = TRUE) -
       min(df[[yvals]], na.rm = TRUE)) * 0.2
  
  y.lev <- max(df[[yvals]], na.rm = TRUE) +
    (max(df[[yvals]], na.rm = TRUE) -
       min(df[[yvals]], na.rm = TRUE)) * 0.05
  
  
  ## 2. Colour scale
  
  if (!is.null(group)) {
    if (is.null(colscale)) {
      colscale <- scales::hue_pal()(
        length(unique(df[[group]]))
      )
    } else {
      colscale <- get(colscale)
    }
  }
  
  
  ## 3. Statistical test
  
  stat_test <- if (!is.null(comp)) {
    if (shapiro.test(comp[[yvals]])$p.value > 0.05) {
      "t.test"
    } else {
      "wilcox.test"
    }
  } else {
    if (shapiro.test(df[[yvals]])$p.value > 0.05) {
      "t.test"
    } else {
      "wilcox.test"
    }
  }
  
  
  ## 4. Base plot
  
  plot <- ggplot(
    df,
    aes(
      x = get(xvals),
      y = get(yvals)
    )
  ) +
    labs(
      title = title,
      x = x.axis,
      y = y.axis
    ) +
    theme_bw(base_size = base_size) +
    theme(
      legend.position = "bottom"
    )
  
  
  ## 5. Reference / comparison layer
  
  if (!is.null(comp)) {
    
    if (type == "boxplot") {
      
      plot <- plot +
        geom_boxplot(
          data = comp,
          aes(group = get(group)),
          outlier.size = point_size,
          color = "grey",
          fill = "grey"
        )
      
    } else {
      
      plot <- plot +
        geom_errorbar(
          data = comp,
          aes(
            ymin = get(yvals) - ci,
            ymax = get(yvals) + ci
          ),
          color = "grey",
          width = 0.1
        ) +
        geom_point(
          data = comp,
          aes(
            x = get(xvals),
            y = get(yvals)
          ),
          color = "grey",
          size = point_size
        )
    }
  }
  
  
  ## 6. Main data
  
  if (type == "boxplot") {
    
    plot <- plot +
      geom_boxplot(
        aes(
          group = get(group),
          fill = get(group)
        ),
        outlier.size = point_size
      ) +
      geom_jitter(
        width = 0.01,
        size = point_size
      )
    
    if (!is.null(group)) {
      plot <- plot +
        scale_fill_manual(
          values = colscale,
          name = group
        )
    }
    
  } else {
    
    plot <- plot +
      geom_line(
        aes(colour = get(group))
      ) +
      geom_errorbar(
        aes(
          ymin = get(yvals) - ci,
          ymax = get(yvals) + ci,
          colour = get(group)
        ),
        width = 0.1
      ) +
      geom_point(
        aes(colour = get(group)),
        size = point_size
      )
    
    if (!is.null(group)) {
      plot <- plot +
        scale_color_manual(
          values = colscale,
          name = group
        )
    }
  }
  
  
  ## 7. Statistical annotation layers
  
  if (type == "boxplot") {
    
    # Add significance annotation (*, **, ***)
    if (add_p_signif) {
      
      if (!is.null(facets)) {
        
        plot <- plot +
          stat_compare_means(
            data = df,
            method = stat_test,
            label = "p.signif",
            label.y = y.lev,
            label.x.npc = "center"
          )
        
      } else {
        
        plot <- plot +
          stat_compare_means(
            data = df,
            method = stat_test,
            label = "p.signif",
            label.y = y.lev
          )
      }
    }
    
    
    # Add formatted p-value
    if (add_p_format) {
      
      plot <- plot +
        stat_compare_means(
          data = df,
          method = stat_test,
          label = "p.format",
          label.y.npc = "bottom"
        )
    }
    
  } else if (!is.null(pvalues) && add_p_signif) {
    
    # Create a dataframe which includes a y position
    # for placing the significance annotations
    
    y_max_df <- df |>
      dplyr::group_by(
        dplyr::across(
          dplyr::all_of(c(facets, xvals))
        )
      ) |>
      dplyr::summarise(
        y_pos = max(.data[[yvals]], na.rm = TRUE) +
          0.1 * max(df[[yvals]], na.rm = TRUE),
        .groups = "drop"
      )
    
    pvalues_annot <- dplyr::left_join(
      pvalues,
      y_max_df,
      by = c(facets, xvals)
    )
    
    plot <- plot +
      geom_text(
        data = pvalues_annot,
        aes(
          x = .data[[xvals]],
          y = y_pos,
          label = sig
        ),
        colour = "black"
      )
  }
  
  
  ## 8. Faceting
  
  if (!is.null(facets)) {
    
    plot <- plot +
      if (length(facets) == 2) {
        facet_grid(
          as.formula(
            paste(facets[1], "~", facets[2])
          )
        )
      } else {
        facet_grid(
          as.formula(
            paste("~", facets)
          )
        )
      }
  }
  
  
  ## 9. Y limits
  
  plot <- plot +
    coord_cartesian(
      ylim = c(NA, y.lim)
    )
  
  
  ## 10. Theme
  
  if (type == "boxplot") {
    plot <- plot +
      theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
      )
  }
  
  
  ## 11. Return plot
  
  return(plot)
}


#' Function to plot the pressure : D interaction in photosynthetic parameters, showing the raw data, the model fit with a 95% confidence ribbon, the slope of each group and the significance of the difference between slopes
#' 
#' arguments:
#' 1. data - data frame from which to take the points to plot. Must contain x_var (numeric), group_var (a factor) and y_var.
#' 2. model - a fitted linear model (lm) used to create prediction points and to extract the slopes and the interaction p-value. It must include the x_var * group_var interaction, using the same column names as in data.
#' 3. y_var - name (as a string) of the column holding the y variable.
#' 4. x_var - name (as a string) of the column holding the x variable. Must be numeric. Defaults to "pressure_n".
#' 5. group_var - name (as a string) of the column holding the grouping variable. Must be a factor. Defaults to "D_f".
#' 6. colscale - name (as a string) of an existing vector of colours to use looked up with get(). Can be left NULL (the default), in which case ggplot's default hue palette is used.
#' 7. title - string with the plot title. Optional, defaults to NULL (no title).
#' 8. y_label - string with the y axis label. Optional, defaults to the y_var column name.
#' 9. x_label - string with the x axis label. Optional, defaults to NULL (no label).
#' 10. legend - string with the legend title. Optional, defaults to the group_var column name.
#' 11. text_size - base font size of the plot (theme_classic), defaults to 12.
#' 12. line_size - base line size of the plot theme (axes etc.), defaults to 1. It does not change the thickness of the fitted lines.
#' 13. show_slope - boolean, defaults to TRUE. Controls whether the slope value of each group is printed on the plot.
#'
#' How it works:
#' - Checks that x_var is numeric and group_var is a factor, and stops with an  error if not.
#' - Builds a prediction grid: 100 evenly spaced x values per group, spanning the range of x_var observed in that group.
#' - Uses predict() on the model to get fitted values and 95% confidence limits (fit, lwr, upr) for the grid.
#' - Extracts the slope of each group from the model coefficients: the x_var coefficient for the reference level, and that plus the interaction coefficient for the second level.
#' - Extracts the p-value of the interaction term (x_var:group_var) and turns it into a label: "ns" if p > 0.05, "p < 0.001" if p < 0.001, otherwise "p = 0.xx" (two decimals).
#' - Works out where to place each slope label at the left edge of the plot: above the confidence ribbon for the group with the higher fitted value, and below the ribbon for the other group, offset by 3% of the y range.
#' - Draws the raw points (semi-transparent), the confidence ribbons, the fitted lines, the colour scale, the labels, the p-value label (centred on the x range, near the bottom) and, if show_slope = TRUE, the slope labels ("Slope = x.xxx").
#'
#' Returns: a ggplot object, which can be further modified with + or printed.

#'
#' @param data data frame
#' @param model fitted lm model with an x_var * group_var interaction
#' @param y_var string, y variable name
#' @param x_var string, x variable name (numeric). Defaults to "pressure_n"
#' @param group_var string, grouping variable name (factor). Defaults to "D_f"
#' @param colscale string, name of a colour vector (optional)
#' @param title string, plot title (optional)
#' @param y_label string, y axis label (optional)
#' @param x_label string, x axis label (optional)
#' @param legend string, legend title (optional)
#' @param text_size numeric, base font size. Defaults to 12
#' @param line_size numeric, base line size of the theme. Defaults to 1
#' @param show_slope logical, print slope labels? Defaults to TRUE
#' @return ggplot object
#' @examples
#' mod <- lm(Fv_Fm ~ pressure_n * D_f, data = df)
#' plot_interaction_ci(df, mod, y_var = "Fv_Fm")
#' plot_interaction_ci(df, mod, y_var = "Fv_Fm", y_label = "Fv/Fm",
#'                     x_label = "Pressure", legend = "D", show_slope = FALSE)
#' @export

plot_interaction_ci <- function(
    data,
    model,
    y_var,
    x_var = "pressure_n",
    group_var = "D_f",
    colscale = NULL,
    title = NULL,
    y_label = NULL,
    x_label = NULL,
    legend = NULL,
    text_size = 12,
    line_size = 1,
    show_slope = TRUE
) {
  
  # Check that x_var is numeric
  if (!is.numeric(data[[x_var]])) {
    stop("x_var must be numeric")
  }
  
  # Check that group_var is a factor
  if (!is.factor(data[[group_var]])) {
    stop("group_var must be a factor")
  }
  
  # Create prediction grid
  newdat <- data |>
    dplyr::summarise(
      !!sym(x_var) := seq(
        min(get(x_var), na.rm = TRUE),
        max(get(x_var), na.rm = TRUE),
        length.out = 100
      ),
      .by = group_var
    )
  
  # Predictions + 95% CI
  pred <- predict(model, newdata = newdat, interval = "confidence")
  newdat <- cbind(newdat, pred)
  
  # Extract slopes from model
  b <- coef(model)
  
  slopes <- data.frame(
    D_f = factor(
      levels(data[[group_var]]),
      levels = levels(data[[group_var]])
    ),
    slope = c(
      b[x_var],
      b[x_var] + b[paste0(x_var, ":", group_var, "20")]
    )
  )
  
  # Extract significance of difference between slopes
  sm <- summary(model)
  
  int_term <- grep(
    paste0("^", x_var, ":"),
    rownames(sm$coefficients),
    value = TRUE
  )
  
  slope_diff_p <- sm$coefficients[int_term, "Pr(>|t|)"]
  
  p_lab <- ifelse(
    slope_diff_p > 0.05,
    "ns",
    ifelse(
      slope_diff_p < 0.001,
      "p < 0.001",
      paste0(
        "p = ",
        formatC(slope_diff_p, digits = 2, format = "f")
      )
    )
  )
  
  # Position slope labels
  slopes <- slopes |>
    dplyr::left_join(
      newdat |>
        dplyr::filter(
          !!sym(x_var) == min(!!sym(x_var), na.rm = TRUE)
        ) |>
        dplyr::summarise(
          max = max(fit),
          upper = max(upr, na.rm = TRUE),
          lower = min(lwr, na.rm = TRUE),
          .by = !!sym(group_var)
        ),
      by = setNames(group_var, "D_f")
    ) |>
    dplyr::mutate(
      y_pos = ifelse(
        max == max(max),
        upper + 0.03 * diff(range(data[[y_var]], na.rm = TRUE)),
        lower - 0.03 * diff(range(data[[y_var]], na.rm = TRUE))
      )
    )
  
  # Build colour palette
  if (!is.null(colscale)) {
    colscale <- get(colscale)
  } else {
    colscale <- scales::hue_pal()(
      length(unique(data[[group_var]]))
    )
  }
  
  # Build base plot
  plot <- ggplot(
    data,
    aes(
      x = get(x_var),
      y = get(y_var),
      color = !!sym(group_var),
      fill = !!sym(group_var)
    )
  ) +
    
    geom_point(alpha = 0.6) +
    
    geom_ribbon(
      data = newdat,
      aes(
        ymin = lwr,
        ymax = upr,
        y = NULL
      ),
      alpha = 0.25,
      color = NA
    ) +
    
    geom_line(
      data = newdat,
      aes(y = fit),
      linewidth = 1.2
    ) +
    
    scale_color_manual(
      values = colscale,
      aesthetics = c("colour", "fill")
    ) +
    
    labs(
      x = x_label,
      y = ifelse(is.null(y_label), y_var, y_label),
      title = title,
      color = ifelse(is.null(legend), group_var, legend),
      fill = ifelse(is.null(legend), group_var, legend)
    ) +
    
    annotate(
      "text",
      x = mean(range(data[[x_var]], na.rm = TRUE)),
      y = min(range(data[[y_var]], na.rm = TRUE)) +
        0.1 * diff(range(data[[y_var]], na.rm = TRUE)),
      label = p_lab,
      hjust = 0.5
    ) +
    
    theme_classic(
      base_size = text_size,
      base_line_size = line_size
    )
  
  # Optional slope annotation
  
  if (show_slope) {
    
    plot <- plot +
      geom_text(
        data = slopes,
        aes(
          x = min(data[[x_var]], na.rm = TRUE),
          y = y_pos,
          label = paste0(
            "Slope = ",
            round(slope, 3)
          )
        ),
        color = "black",
        hjust = 0,
        show.legend = FALSE
      )
  }
  
  
  return(plot)
}


#' Function to calculate significance between pairs of conditions
#' 
#' arguments:
#' 1. data - data frame from which to calculate. Must contain the columns named in comparison_column, grouping_column (if used) and variable, plus a column called "irr" (irradiance level). If time = TRUE, it must also contain a column called "time".
#' 2. grouping_column - name (as a string) of the column whose levels define the separate sets of comparisons to make, e.g. "species". Can be left NULL (the default) if only making one comparison set.
#' 3. comparison_column - name (as a string) of the column in which the final pairs are found, e.g. "treatment". Within each subset it must contain exactly two levels, otherwise that subset is skipped.
#' 4. variable - name (as a string) of the numeric column to compare, e.g. "yield".
#' 5. time - boolean, defaults to TRUE. If TRUE, comparisons are made separately for each unique value in data$time. If FALSE (or NULL), time is ignored and all data is used together.
#'
#' How it works:
#' - Loops through every combination of time (optional), grouping_column level (optional) and irr level, so each comparison is made on a separate subset of the data.
#' - In each subset, finds the levels of comparison_column. If there are not exactly two, the subset is skipped. Otherwise the values of "variable" are extracted for each level as group1 and group2.
#' - Chooses the statistical test automatically:
#'     * "none": if either group has zero variance (only one unique value), no test is run and p = 1 is returned.
#'     * "t": if both groups have more than 2 observations and pass the Shapiro-Wilk normality test (p > 0.05), a Welch two-sample t-test is used.
#'     * "w": otherwise, a Wilcoxon rank-sum (Mann-Whitney) test is used.
#' - Converts the p-value to significance stars: *** (p < 0.001), ** (p < 0.01), * (p < 0.05), . (p < 0.1, marginal), or a blank space (p >= 0.1). A NaN p-value is treated as p = 1 with no stars. These are the same thresholds and labels as add_stars().
#' - Collects one row per comparison and binds them into a single data frame.
#'
#' Returns: a data frame with one row per comparison, containing:
#'   comp1, comp2 - the two levels being compared
#'   irr          - the irradiance level of the subset
#'   p.value      - p-value of the test (1 if no test was possible)
#'   sig          - significance stars
#'   n1, n2       - number of observations in each group
#'   stat         - test used: "t", "w" or "none"
#'   time         - sampling time (only if time = TRUE)
#'   comparison   - level of grouping_column (only if grouping_column is set)
#'
#' @param data data frame
#' @param grouping_column string, name of the grouping column (optional)
#' @param comparison_column string, name of the column containing the pair
#'   to compare
#' @param variable string, name of the numeric column to compare
#' @param time logical, compare separately at each time point? Defaults to TRUE
#' @return data frame with columns: comp1, comp2, irr, p.value, sig, n1, n2,
#'   stat, and optionally time and comparison
#' @examples
#' signif(df, grouping_column = "species", comparison_column = "treatment",
#'        variable = "yield", time = TRUE)
#' signif(df, comparison_column = "treatment", variable = "yield", time = FALSE)
#' @export

signif <- function(data,
                   grouping_column = NULL,
                   comparison_column,
                   variable,
                   time = TRUE) {
  
  rows <- list()
  i <- 1
  
  ## 1. Time Handling
  time_levels <- if (isFALSE(time) || is.null(time)) NA else unique(data$time)
  
  for (sampling_time in time_levels) {
    
    subset_time <- if (isFALSE(time) || is.null(time)) {
      data
    } else {
      dplyr::filter(data, time == sampling_time)
    }
    
    ## 2. Group Handling
    group_levels <- if (is.null(grouping_column)) NA else unique(subset_time[[grouping_column]])
    
    for (comp in group_levels) {
      
      subset_group <- if (is.null(grouping_column)) {
        subset_time
      } else {
        dplyr::filter(subset_time, .data[[grouping_column]] == comp)
      }
      
      for (rad in unique(subset_group$irr)) {
        
        subset2 <- dplyr::filter(subset_group, irr == rad)
        levs <- unique(subset2[[comparison_column]])
        if (length(levs) != 2) next
        
        group1 <- subset2 |>
          dplyr::filter(.data[[comparison_column]] == levs[1]) |>
          dplyr::pull(.data[[variable]])
        
        group2 <- subset2 |>
          dplyr::filter(.data[[comparison_column]] == levs[2]) |>
          dplyr::pull(.data[[variable]])
        
        ## 3. Stat test
        
        #check if groups have 0 variance
        if (length(unique(group1)) <= 1 || length(unique(group2)) <= 1) {
          
          stat_test <- list(p.value = 1)
          stat <- "none"
          
        } else if (length(group1) > 2 && length(group2) > 2 &&
                   shapiro.test(group1)$p.value > 0.05 &&
                   shapiro.test(group2)$p.value > 0.05) {
          
          stat_test <- t.test(group1, group2)
          stat <- "t"
          
        } else {
          
          stat_test <- wilcox.test(group1, group2)
          stat <- "w"
        }
        
        
        ## 4. Conditional row building
        row <- list(
          comp1 = levs[1],
          comp2 = levs[2],
          irr = rad,
          p.value = ifelse(is.nan(stat_test$p.value), 1, stat_test$p.value),
          sig = dplyr::case_when(
            is.nan(stat_test$p.value) ~ " ",
            stat_test$p.value < 0.001 ~ "***",
            stat_test$p.value < 0.01  ~ "**",
            stat_test$p.value < 0.05  ~ "*",
            stat_test$p.value < 0.1   ~ ".",   # NEW: marginal
            TRUE ~ " "
          ),
          n1 = length(group1),
          n2 = length(group2),
          stat = stat
        )
        
        if (!isFALSE(time) && !is.null(time)) {
          row$time <- sampling_time
        }
        
        if (!is.null(grouping_column)) {
          row$comparison <- comp
        }
        
        rows[[i]] <- row
        i <- i + 1
      }
    }
  }
  
  dplyr::bind_rows(rows)
}


#' Function to calculate summary statistics (N, mean, SD, SE, CI) for a variable
#' arguments:
#' 1. data - data frame from which to calculate. Must contain the columns named in 'measurevar' and 'groupvars'.
#' 2. measurevar - name (as a string) of the numeric column to summarize, e.g. "len".
#' 3. groupvars - character vector of column names to group by, e.g. c("supp", "dose"). Can be left NULL (the default) if no grouping is needed, in which case the data is treated as a single group.
#' 4. na.rm - boolean, defaults to FALSE. If TRUE, missing values are excluded from N, mean and SD, so SE and CI are based on the observations actually used.
#' 5. conf.interval - confidence level for the interval, defaults to .95.
#' 6. .drop - boolean, defaults to TRUE. If TRUE, combinations of grouping variables with no observations are dropped from the output.
#'
#' How it works:
#' - Splits the data by groupvars (using plyr::ddply) and computes, per group: N (count, NA-aware via length2), mean, and sample standard deviation.
#' - Renames the "mean" column to the name of measurevar.
#' - Calculates the standard error of the mean: se = sd / sqrt(N).
#' - Calculates the confidence interval half-width: ci = se * t, where t is the t-distribution quantile for the requested confidence level with N - 1 degrees of freedom (computed per group, so smaller groups get wider intervals).
#'
#' Returns: a data frame with one row per group, containing the grouping columns plus N, the mean (named after measurevar), sd, se and ci. The confidence interval for each group is mean +/- ci.
#'
#' @param data data frame
#' @param measurevar string, name of the column to summarize
#' @param groupvars character vector of grouping column names (optional)
#' @param na.rm logical, remove NAs? Defaults to FALSE
#' @param conf.interval numeric, confidence level. Defaults to .95
#' @param .drop logical, drop empty groups? Defaults to TRUE
#' @return data frame with columns: grouping variables, N, mean (named after
#'   measurevar), sd, se, ci
#' @examples
#' summarySE(ToothGrowth, measurevar = "len", groupvars = c("supp", "dose"))
#' @export

summarySE <- function(data=NULL,
                      measurevar, 
                      groupvars=NULL,
                      na.rm=FALSE,
                      conf.interval=.95, 
                      .drop=TRUE) {
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

#' Helper function to convert p-values into significance stars (uses the same thresholds and labels as signif(), so the two are consistent)
#'
#' arguments:
#' 1. p - numeric vector of p-values. Can be a single value or a whole vector (e.g. a column of a results table). NA and NaN values are allowed.
#'
#' How it works:
#' - Goes through the p-values with nested ifelse() calls, so it is vectorised and returns one label per p-value, in the same order.
#' - Assigns the label using these thresholds:
#'     * p < 0.001         -> "***"
#'     * p < 0.01          -> "**"
#'     * p < 0.05          -> "*"
#'     * p < 0.1           -> "."  (marginal / trend)
#'     * p >= 0.1 or NA    -> " " (a single space, not significant)
#'
#' Returns: a character vector the same length as p, containing "***", "**", "*", "." or " ".
#'
#' @param p numeric vector of p-values
#' @return character vector of significance stars, same length as p
#' @examples
#' add_stars(0.0004)                          # "***"
#' add_stars(c(0.03, 0.07, 0.4, NA))          # "*" "." " " " "
#' add_stars(0.05)                            # "."
#' df$sig <- add_stars(df$p.value)
#' @export

add_stars <- function(p) {
  ifelse(is.na(p), " ",
         ifelse(p < 0.001, "***",
                ifelse(p < 0.01, "**",
                       ifelse(p < 0.05, "*",
                              ifelse(p < 0.1, ".", " ")))))
}


#------------- END ---------------------