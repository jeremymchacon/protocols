# bleedthrough calculation functions

get_start_growth_time = function(x, y, rolling_window = 11){
  sort_order = sort(x, index.return = TRUE)$ix
  x=  x[sort_order]
  y = y[sort_order]
  slopes = numeric()
  r2 = numeric()
  for (k in 1:(length(y)-rolling_window)){
    m1 = lm(y[k:(k+rolling_window)] ~ x[k:(k+rolling_window)])
    slopes[k] = coef(m1)[2]
    r2[k] = summary(m1)$r.squared
  }
  start_time = which(slopes > 0 & r2 > 0.9)[1]
  start_time = x[start_time]
  return(start_time)
}

get_max_growth_time = function(x, y, rolling_window = 21){
  sort_order = sort(x, index.return = TRUE)$ix
  x=  x[sort_order]
  y = y[sort_order]
  y = zoo::rollmean(y, rolling_window)
  max_growth_time = which(diff(y) == max(diff(y)))[1] + 21
  max_growth_time = x[max_growth_time]
  return(max_growth_time)
}

get_coef = function(df, x_col, y_col, R2_thresh = 0){
  x = df[,x_col]
  y = df[,y_col]
  m1 = lm(y ~ x)
  result = coef(m1)[2]
  if (summary(m1)$r.squared < R2_thresh){
    result = 0
  }
  return(result)
}


