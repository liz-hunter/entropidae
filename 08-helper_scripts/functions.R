# function to grab extremes for any summary statistic
get_extremes <- function(df, sum_stat, n = 10, which = c("high", "low")) {
  ToB <- match.arg(which)
  
  df_filtered <- df %>%
    filter(!is.na(.data[[sum_stat]])) %>%
    arrange(.data[[sum_stat]])
  
  extreme_rows <- if (ToB == "low") {
    slice_head(df_filtered, n = n)
  } else {
    slice_tail(df_filtered, n = n)
  }
  
  extreme_rows %>%
    mutate(stat_group = ToB)
}


