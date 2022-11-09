source("./growthcurve_functions.r")

read_tecan_csv = function(path, sep = ",",measure = "OD"){
  OD = read.table(path, sep = sep, as.is = TRUE, row.names = 1)
  # transpose it so wells are columns, then convert to dataframe
  
  OD = t(OD)
  OD = data.frame(OD)
  names(OD)[1:3] = c("cycle","seconds","temp")
  
  # convert to long format by "gathering"
  OD = OD %>%
    pivot_longer(-c(cycle, seconds, temp), names_to = "well", values_to = measure) %>%
    mutate(hour = seconds / 60 / 60) %>%
    mutate(well = factor(well, levels = str_sort(unique(well), numeric = TRUE)))
}

plot_plate = function(OD, measure = "OD", color = "none"){
  measure = ensym(measure)
  if (color == "none"){
    fig =   OD %>%
      mutate(row = str_sub(well, 1, 1),
             col = str_sub(well, 2, -1)) %>%
      mutate(col = factor(col, 
                          levels = str_sort(unique(col),numeric = TRUE))) %>%
      ggplot(aes(x = hour, y = !!measure))+
      geom_line()+
      facet_grid(row~col)
  }else{
    color = ensym(color)
    fig =   OD %>%
      mutate(row = str_sub(well, 1, 1),
             col = str_sub(well, 2, -1)) %>%
      mutate(col = factor(col, 
                          levels = str_sort(unique(col),numeric = TRUE))) %>%
      ggplot(aes(x = hour, y = !!measure, color = !!color))+
      geom_line()+
      facet_grid(row~col)
  }
  return(fig)
}

plot_plate_growthrates = function(growth_rates, 
                                  response = "r", color = "none"){
  response = ensym(response)
  if (color == "none"){
    fig =   growth_rates %>%
      mutate(row = str_sub(well, 1, 1),
             col = str_sub(well, 2, -1)) %>%
      mutate(col = factor(col, 
                          levels = str_sort(unique(col),numeric = TRUE))) %>%
      ggplot(aes(x = 1, y = !!response))+
      geom_bar(stat = "identity")+
      facet_grid(row~col)
  }else{
    color = ensym(color)
    fig =   growth_rates %>%
      mutate(row = str_sub(well, 1, 1),
             col = str_sub(well, 2, -1)) %>%
      mutate(col = factor(col, 
                          levels = str_sort(unique(col),numeric = TRUE))) %>%
      ggplot(aes(x = 1, y = !!response, fill = !!color))+
      geom_bar(stat = "identity")+
      facet_grid(row~col)
  }
  return(fig)
}

fit_all_loglinear = function(OD, measure = "OD", groups = c("well"), tries = 5){
  measure = ensym(measure)
  OD %>%
    group_by_at(groups) %>%
    summarize(model_fit = fit_loglinear(hour, !!measure, tries),
              fit_variable = c("r", "y0", "lag", "end")) %>%
    ungroup() %>%
    pivot_wider(names_from = fit_variable, values_from = model_fit)
}


fit_all_logistic = function(OD, measure = "OD", groups = c("well"), tries = 5){
  measure = ensym(measure)
  OD %>%
    group_by_at(groups) %>%
    summarize(model_fit = fit_logistic(hour, !!measure, tries),
              fit_variable = c("r", "lag", "K", "y0")) %>%
    ungroup() %>%
    pivot_wider(names_from = fit_variable, values_from = model_fit)
}


fit_all_baranyi = function(OD, measure = "OD", groups = c("well"), tries = 5){
  measure = ensym(measure)
  OD %>%
    group_by_at(groups) %>%
    summarize(model_fit = fit_baranyi(hour, !!measure, tries),
              fit_variable = c("r", "lag", "ymax", "y0")) %>%
    ungroup() %>%
    pivot_wider(names_from = fit_variable, values_from = model_fit)
}


fit_all_loglinear_lm = function(OD, measure = "OD", groups = c("well"),
                                surround = 5){
  measure = ensym(measure)
  OD %>%
    group_by_at(groups) %>%
    summarize(model_fit = fit_loglinear_lm(hour, !!measure, surround),
              fit_variable = c("r")) %>%
    ungroup() %>%
    pivot_wider(names_from = fit_variable, values_from = model_fit)
}

get_yields = function(OD, measure = "OD", groups = c("well"), smooth = 1){
  # returns a dataframe grouped by the stated groups, that has the max of the
  # measure as "yield." Optionally, if smooth > 1, it computes a moving average of window
  # "smooth" first
  measure = ensym(measure)
  OD %>%
    group_by(well) %>%
    arrange(cycle) %>%
    mutate(measure = rollmean(!!measure, smooth, fill = NA)) %>%
    group_by_at(groups) %>%
    summarize(yield = max(measure, na.rm = TRUE)) %>%
    ungroup()
}
