

scenario_plot <- function(listdat) {
  ggplot() + 
    geom_density(data = listdat[["df1"]], aes(x = value), fill = "gray70", alpha = 0.6) + 
    geom_density(data = listdat[["df2"]], aes(x = value), fill = "black", alpha = 0.4) + 
    labs(x = "", y = "Density") + 
    theme_classic()
}

create_dat <- function(params){

  df1 <- rpearson(1e5, moments = c(mean = params$mean_g1,  
                               variance = params$variance_g1,  
                               skewness = params$skewness_g1,   
                               kurtosis = params$kurtosis_g1) )
  df2 <- rpearson(1e5, moments = c(mean = params$mean_g2, 
                               variance = params$variance_g2,   
                               skewness = params$skewness_g2,   
                               kurtosis = params$kurtosis_g2) )

  return(list(df1 = data.frame(value = df1), df2 = data.frame(value = df2), params = params))
}


