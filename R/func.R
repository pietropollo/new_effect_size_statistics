

scenario_plot <- function(listdat, print = TRUE, xpos = 6) {
# Compute density for both datasets
  dens1 <- density(listdat[["df1"]]$value)
  dens2 <- density(listdat[["df2"]]$value)
  
  # Get axis limits
  xmin <- min(dens1$x, dens2$x)
  ymax <- max(dens1$y, dens2$y)
  
  # Spread annotation labels more vertically
  n_params <- ncol(listdat$params)
  y_top <- ymax * 0.95
  y_bottom <- ymax * 0.45  # lower bottom for more space
  y_positions <- seq(y_top, y_bottom, length.out = n_params)

 p <-  ggplot() + 
    geom_density(data = listdat[["df1"]], aes(x = value), fill = "gray70", alpha = 0.6) + 
    geom_density(data = listdat[["df2"]], aes(x = value), fill = "black", alpha = 0.4) + 
    labs(x = "", y = "Density") + 
    theme_classic() + theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16))

if(print){
	p  <- p + annotate(
      "text",
      x = xmin+1,
      y = y_positions,
      label = paste(names(listdat$params), "=", listdat$params[1, ]),
      size = 3
    )
}	
return(p)
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

plot_scenario <- function(params, print = TRUE, xpos = 6) {
  dat <- create_dat(params)
  scenario_plot(dat, print = print, xpos = xpos)
}

# Function for visualising the scenarios in a grid layout
multi_plot_scenarios <- function(scenarios, folder = "figs/") {
  plots <- list() # Initialize an empty list to store plots
  for(i in 1:nrow(scenarios)) {
    params <- scenarios[i,]
    plots[[i]] <- plot_scenario(params, print = TRUE, xpos = 3)
  }
# Combine all plots into a grid layout
ggsave(filename = paste0(folder, "scenario_plots_", Sys.time(), ".png"), plot = wrap_plots(plots, ncol = 2), width = 10, height = 10, dpi = 300)
}


