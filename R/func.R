
#' @title scenario_plot
#' @description Plot Overlaid Density Plots with Parameter Annotations
#' @details This function plots overlaid kernel density curves from two datasets and optionally
#' annotates the plot with parameter values (mean, variance, etc.) provided in the input.
#' @param listdat A list containing three named elements:
#'   \describe{
#'     \item{\code{df1}}{A data frame with a numeric column named \code{value}.}
#'     \item{\code{df2}}{A second data frame with a column named \code{value}.}
#'     \item{\code{params}}{A single-row data frame of named parameters (e.g., means, variances).}
#'   }
#' @param print Logical; if \code{TRUE}, annotate the plot with parameter labels. Default is \code{TRUE}.
#' @param xpos Numeric; horizontal position for text annotation. Default is 6.
#' @return A \code{ggplot} object showing overlaid densities and optional annotations.
#' @import ggplot2
#' @export
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

#' @title create_dat
#' @description Generate Random Data from Pearson Distributions
#' @details Generates two synthetic datasets using the Pearson distribution based on provided
#' moments (mean, variance, skewness, kurtosis) for each group.
#' @param params A named data frame with one row and the following columns:
#'   \code{mean_g1}, \code{variance_g1}, \code{skewness_g1}, \code{kurtosis_g1}, 
#'   \code{mean_g2}, \code{variance_g2}, \code{skewness_g2}, \code{kurtosis_g2}.
#' @return A list containing:
#'   \describe{
#'     \item{\code{df1}}{A data frame of 100,000 random values for group 1.}
#'     \item{\code{df2}}{A data frame of 100,000 random values for group 2.}
#'     \item{\code{params}}{The original parameter input.}
#'   }
#' @import PearsonDS
#' @export
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

#' @title plot_scenario
#' @description Wrapper to Create and Plot a Scenario Based on Parameters
#' @details This is a wrapper function that first generates data using \code{create_dat()} and then
#' plots it using \code{scenario_plot()}.
#' @param params A single-row data frame with named parameters used to generate the distributions.
#' @param print Logical; if \code{TRUE}, annotate the plot with parameter labels. Default is \code{TRUE}.
#' @param xpos Numeric; horizontal position for text annotation. Default is 6.
#' @return A \code{ggplot} object.
#' @seealso \code{\link{create_dat}}, \code{\link{scenario_plot}}
#' @export
  plot_scenario <- function(params, print = TRUE, xpos = 6) {
    dat <- create_dat(params)
    scenario_plot(dat, print = print, xpos = xpos)
  }

#' @title multi_plot_scenarios
#' @description Plot Multiple Scenarios in a Grid and Save to File
#' @details Iterates over a set of parameter rows and plots each using \code{plot_scenario()}, combining
#' the results into a grid layout and saving as a PNG file.
#' @param scenarios A data frame where each row contains parameters for one scenario
#'   (columns as in \code{create_dat()}).
#' @param folder Path to the output folder for the plot file. Defaults to \code{"figs/"}.
#' @return Saves a grid of plots as a PNG file and returns \code{NULL} invisibly.
#' @import ggplot2
#' @import patchwork
#' @export
  multi_plot_scenarios <- function(scenarios, folder = "figs/") {
    plots <- list() # Initialize an empty list to store plots
    for(i in 1:nrow(scenarios)) {
      params <- scenarios[i,]
      plots[[i]] <- plot_scenario(params, print = TRUE, xpos = 3)
    }
  # Combine all plots into a grid layout
  ggsave(filename = paste0(folder, "scenario_plots_", Sys.time(), ".png"), plot = wrap_plots(plots, ncol = 2), width = 10, height = 10, dpi = 300)
  }


