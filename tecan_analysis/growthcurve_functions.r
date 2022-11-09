library(tidyverse)
library(zoo)

fit_loglinear = function(x, y, tries = 5){
  # x = time, y = growth, tries = # of times to try if fitting
  # y should be untransformed. 
  #
  # 
  # remove NA time points if there are any
  x = x[!is.na(y)]
  y = y[!is.na(y)]
  x = x[y > 0]
  y = y[y > 0]
  #figure out start lag guess
  ymax = max(y)
  ymin = min(y)
  saturation_time = x[y > (ymin + ((ymax - ymin) * 0.95))][1]
  if(is.na(saturation_time)) saturation_time = max(x)
  lag_time = x[y > (ymin + ((ymax - ymin) * 0.025))][1]
  if(is.na(lag_time)) lag_time = min(x)
  
  y = log(y)
  
  # remove NA time points again if produced from log-transform
  x = x[!is.na(y)]
  y = y[!is.na(y)]
  
  m1 <- NA
  while(is.na(m1[1]) && (tries > 0)){
    r_start = runif(1, 0, 1)
    y0_start = min(y) + rnorm(1, 0, 5)
    lag_start = lag_time + rnorm(1, 0, 8)
    while(lag_start < min(x)){ lag_start = lag_time + rnorm(1, 0, 8)}
    end_start = saturation_time + rnorm(1, 0, 5)
    while(end_start > max(x)){ end_start = saturation_time + rnorm(1, 0, 5)}
    
    m1 <- tryCatch( nls(y ~ log_linear_logunits(x, r, y0, lag, end),
                        start = list(r = r_start, 
                                     y0 = y0_start, 
                                     lag = lag_start, 
                                     end = end_start)),
                    error = function(e) return(NA))
    tries = tries - 1
  }
  if (is.na(m1[1])){
    return(c(NA, NA,NA,NA))
  }
  r = coef(m1)[1]
  y0 = exp(coef(m1)[2])
  lag = coef(m1)[3]
  end = coef(m1)[4]
  return(c(r, y0, lag, end))
}


fit_logistic = function(x, y, tries = 5){
    m1 <- NA
    
    # remove NA time points
    x = x[!is.na(y)]
    y = y[!is.na(y)]
    #figure out start lag guess
    lag_guess = guess_half_max(x,y)
    while(is.na(m1[1]) && (tries > 0)){
      m1 <- tryCatch( nls(y ~ logistic(x, r, lag, K, y0),
                          start = list(r = runif(1, 0.1, 1), 
                                       lag = lag_guess + sample(-20:20, 1), 
                                       K = max(y), 
                                       y0 = min(y))),
                      error = function(e) return(NA))
      tries = tries - 1
    }
    if (is.na(m1[1])){
      return(c(NA, NA,NA,NA))
    }
    r = coef(m1)[1]
    lag = coef(m1)[2]
    K = coef(m1)[3]
    y0 = coef(m1)[4]
    return(c(r, lag, K, y0))
}

fit_baranyi = function(x, y, tries = 100){
  m1 <- NA
  
  # remove NA time points
  x = x[!is.na(y)]
  y = y[!is.na(y)]
  
  # remove timepoints where y <= 0 (because must take log)
  x = x[y > 0]
  y = y[y > 0]
  
  #figure out start lag guess
  lag_guess = guess_half_max(x,y)
  while(is.na(m1[1]) && (tries > 0)){
    m1 <- tryCatch( nls(y ~ baranyi(x, r, lag, ymax, y0),
                        start = list(r = runif(1, 0.1, 0.4), 
                                     lag = lag_guess + sample(-20:20, 1), 
                                     ymax = max(y), 
                                     y0 = min(y)),
                        lower = c(r = 0, lag = 0, ymax = 1e-12, y0 = 1e-12),
                        algorithm = "port"),
                    error = function(e) return(NA))
    tries = tries - 1
  }
  if (is.na(m1[1])){
    return(c(NA, NA,NA,NA))
  }
  r = coef(m1)[1]
  lag = coef(m1)[2]
  ymax = coef(m1)[3]
  y0 = coef(m1)[4]
  return(c(r, lag, ymax, y0))
}

guess_half_max = function(x, y){
  # this function looks at a logistically-increasing time series
  # and guesses when the growth rate is at a maximum (which is also when
  # the trend is at half-max)
  
  # put in order
  y = y[sort(x, index.return = TRUE)$ix]
  x = sort(x)
  
  #find approximate time of max growth rate, using diffs, on smoothed y
  y2 = diff(y)
  y2 = zoo::rollmean(y2, 11, fill = NA, align = "center")
  half_max_idx = which(y2 == max(y2, na.rm = TRUE))[1]
  half_max = x[half_max_idx]
  return(half_max)
}

baranyi <- function(t, r, lag, ymax, y0){
  logy0 = log(y0)
  logymax = log(ymax)
  At = t + (1 / r) * log(exp(-r * t) + exp(-r * lag) - exp(-r * (t + lag)))
  logy = logy0 + r * At - log(1 + (exp(r * At) - 1) / exp(logymax - logy0))
  y = exp(logy)
  return(y)
}

logistic <- function(t, r, lag, K, y0){
  y0 + (K - y0) / (1 + exp(-r * (t-lag)))
}

log_linear <- function(t, r, y0, lag, end){
  # this expects untransformed y0 and gives back untransformed
  y = rep(0, length(t))
  y[t < lag] = y0
  y[t >= lag & t <= end] = y0 *exp(r * (t[t >= lag & t <= end]-lag))
  y[t > end] = y0 * exp(r * (end - lag))
  return(y)
}

log_linear_logunits <- function(t, r, y0, lag, end){
  # this expects that your y0 is log-transformed
  # it also returns log-transformed values
  y = rep(0, length(t))
  y[t < lag] = y0
  y[t >= lag & t <= end] = y0 + (t[t >= lag & t <= end]-lag) * r
  y[t > end] = y0 + (end - lag) * r
  return(y)
}

fit_loglinear_lm = function(x,y, surround = 5){
  # This uses lm, not nls, which is guaranteed not to fail (this is not to say
  # the results will make sense)
  #
  # To do so, it heavily smooths the y data, then find the point where the difference
  # is maximized, then fits a line to log(y) vs. x for the data surrounding
  # the max rate change of the smoothed data. It uses "surround" datapoints on 
  # either size of the max. By defauly this is 5, so in total 11 datapoints are
  # used for fitting.
  y = y[order(x)]
  x = x[order(x)]
  k = floor(length(x)/10)
  ysmooth = rollmean(y,k, fill = NA)
  max_diff_loc = which(diff(log(ysmooth)) == max(diff(log(ysmooth)), na.rm = TRUE))[1]
  startloc = ifelse(max_diff_loc-surround > 0, max_diff_loc-surround, 1)
  endloc = ifelse(max_diff_loc+surround < length(x), max_diff_loc+surround, length(x))
  
  ysub = y[startloc:endloc]
  xsub = x[startloc:endloc]
  logysub = log(ysub)
  m1 = lm(logysub ~ xsub)
  r = coef(m1)[2] %>% unname
  r
}
