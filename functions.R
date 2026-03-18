
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

#'function to plot things??
#'function 