## ====================================================================##
##    R Functions for New Effect Size Statistics
##    This file contains functions related to the simulation and analysis
##    of effect size statistics, including correlation, skewness, and kurtosis.
##    It includes functions for simulating data, calculating effect sizes,
##    and performing bootstrap and jackknife resampling methods.
## ====================================================================##

## ================================================================
## 1. delta-method plug-in variance for sample skewness  (g1)
##    – returns a list:  estimate, var, se
## ================================================================
  #' @title Delta-method Variance for Sample Skewness
  #' @description Computes the sample skewness and estimates its variance using the delta method.
  #' @param x A numeric vector.
  #' @return A list with the skewness estimate (`est`), its variance (`var`), and standard error (`se`).
  #' @examples
  #' x <- rgamma(25, shape = 5)
  #' skew_delta(x)
  skew_delta <- function(x) {
    n  <- length(x)
    xc <- x - mean(x)
    
    ## central moments up to order 6
    m2 <- mean(xc^2); m3 <- mean(xc^3); m4 <- mean(xc^4)
    m5 <- mean(xc^5); m6 <- mean(xc^6)
    
    ## standardised moments (λr = μr / σ^r)
    l3 <- m3 / m2^(3/2)
    l4 <- m4 / m2^2
    l5 <- m5 / m2^(5/2)
    l6 <- m6 / m2^3
    
    ## Δ-method variance  (Eq. 1 in the previous reply)
    var_g1 <- ( l6
                - 3 * l3 * l5
                + 2.25 * l3^2 * (l4 - 1)
                - 0.25 * l3^2 ) / n
    
    list(est = l3,
        var = var_g1,
        se  = sqrt(var_g1))
}

## ================================================================
## 2. delta-method plug-in variance for sample *excess* kurtosis  (g2)
##    – returns a list:  estimate, var, se
## ================================================================
  #' @title Delta-method Variance for Sample Excess Kurtosis
  #' @description Computes the sample excess kurtosis and estimates its variance using the delta method.
  #' @param x A numeric vector.
  #' @return A list with the excess kurtosis estimate (`est`), its variance (`var`), and standard error (`se`).
  #' @examples
  #' x <- rgamma(25, shape = 5)
  #' kurt_delta(x)
  kurt_delta <- function(x) {
    n  <- length(x)
    xc <- x - mean(x)
    
    ## central moments up to order 8
    m2 <- mean(xc^2); m4 <- mean(xc^4)
    m6 <- mean(xc^6); m8 <- mean(xc^8)
    
    l4 <- m4 / m2^2        # raw kurtosis
    l6 <- m6 / m2^3
    l8 <- m8 / m2^4
    
    ## Δ-method variance  (Eq. 2 in the previous reply)
    var_g2 <- ( l8
                - 4 * l4 * l6
                + 4 * l4^3
                - l4^2 ) / n
    
    g2 <- l4 - 3            # convert to *excess* kurtosis
    
    list(est = g2,
        var = var_g2,
        se  = sqrt(var_g2))
}

## ================================================================
## 3. Non-parametric bootstrap for skewness
##    – returns a list:  point est, bias-corrected est, se, replicates*
## ================================================================
  #' @title Bootstrap Estimation of Skewness
  #' @description Computes bootstrap estimate and standard error of skewness.
  #' @param x A numeric vector.
  #' @param B Number of bootstrap replicates. Default is 2000.
  #' @param bias.correct Logical, whether to apply bias correction. Default is TRUE.
  #' @param return.replicates Logical, whether to return replicate values. Default is FALSE.
  #' @return A list with estimate (`est`), bias-corrected estimate (`est_bc`), variance, and standard error.
  #' @examples
  #' x <- rgamma(25, shape = 5)
  #' boot_skew(x)
  boot_skew <- function(x, B = 2000, bias.correct = TRUE,
                        return.replicates = FALSE) {
    g1 <- function(z) mean((z - mean(z))^3) /
      mean((z - mean(z))^2)^(3/2)
    
    b.reps <- replicate(B, g1(sample(x, replace = TRUE)))
    
    est    <- g1(x)
    est.bc <- if (bias.correct) 2 * est - mean(b.reps) else est
    out    <- list(est     = est,
                  est_bc  = est.bc,
                  var     = sd(b.reps)^2,
                  se      = sd(b.reps))
    if (return.replicates) out$boot <- b.reps
    out
}

## ================================================================
## 4. Non-parametric bootstrap for excess kurtosis
##    – returns a list:  point est, bias-corrected est, se, replicates*
## ================================================================
  #' @title Bootstrap Estimation of Excess Kurtosis
  #' @description Computes bootstrap estimate and standard error of excess kurtosis.
  #' @param x A numeric vector.
  #' @param B Number of bootstrap replicates. Default is 2000.
  #' @param bias.correct Logical, whether to apply bias correction. Default is TRUE.
  #' @param return.replicates Logical, whether to return replicate values. Default is FALSE.
  #' @return A list with estimate (`est`), bias-corrected estimate (`est_bc`), variance, and standard error.
  #' @examples
  #' x <- rgamma(25, shape = 5)
  #' boot_kurt(x)
  boot_kurt <- function(x, B = 2000, bias.correct = TRUE,
                        return.replicates = FALSE) {
    g2 <- function(z) mean((z - mean(z))^4) /
      mean((z - mean(z))^2)^2 - 3
    
    b.reps <- replicate(B, g2(sample(x, replace = TRUE)))
    
    est    <- g2(x)
    est.bc <- if (bias.correct) 2 * est - mean(b.reps) else est
    out    <- list(est     = est,
                  est_bc  = est.bc,
                  var     = sd(b.reps)^2,
                  se      = sd(b.reps))
    if (return.replicates) out$boot <- b.reps
    out
}

## ================================================================
## 5. Jack-knife SE (and bias-correction) for skewness
## ================================================================
  #' @title Jackknife Estimation of Skewness
  #' @description Computes jackknife estimate and standard error for skewness.
  #' @param x A numeric vector.
  #' @param bias.correct Logical, whether to apply bias correction. Default is TRUE.
  #' @param return.replicates Logical, whether to return jackknife replicates. Default is FALSE.
  #' @return A list with estimate (`est`), bias-corrected estimate (`est_bc`), variance, and standard error.
  #' @examples
  #' x <- rgamma(25, shape = 5)
  #' jack_skew(x)
  jack_skew <- function(x, bias.correct = TRUE,
                        return.replicates = FALSE) {
    
    n   <- length(x)
    g1  <- function(z) mean((z - mean(z))^3) /
      mean((z - mean(z))^2)^(3/2)
    
    ## leave-one-out replicates
    theta_i <- vapply(seq_len(n),
                      function(i) g1(x[-i]),
                      numeric(1))
    
    theta   <- g1(x)                 # full-sample estimate
    theta_bar <- mean(theta_i)
    
    ## jack-knife standard error
    se_jack <- sqrt((n - 1) * mean((theta_i - theta_bar)^2))
    
    ## bias correction
    theta_bc <- if (bias.correct) n * theta - (n - 1) * theta_bar else theta
    
    out <- list(est    = theta,
                est_bc = theta_bc,
                var    = se_jack^2,
                se     = se_jack)
    
    if (return.replicates) out$jack <- theta_i
    out
}

## ================================================================
## 6. Jack-knife SE (and bias-correction) for *excess* kurtosis
## ================================================================
  #' @title Jackknife Estimation of Excess Kurtosis
  #' @description Computes jackknife estimate and standard error for excess kurtosis.
  #' @param x A numeric vector.
  #' @param bias.correct Logical, whether to apply bias correction. Default is TRUE.
  #' @param return.replicates Logical, whether to return jackknife replicates. Default is FALSE.
  #' @return A list with estimate (`est`), bias-corrected estimate (`est_bc`), variance, and standard error.
  #' @examples
  #' x <- rgamma(25, shape = 5)
  #' jack_kurt(x)
  jack_kurt <- function(x, bias.correct = TRUE,
                        return.replicates = FALSE) {
    
    n   <- length(x)
    g2  <- function(z) mean((z - mean(z))^4) /
      mean((z - mean(z))^2)^2 - 3
    
    theta_i <- vapply(seq_len(n),
                      function(i) g2(x[-i]),
                      numeric(1))
    
    theta     <- g2(x)
    theta_bar <- mean(theta_i)
    
    se_jack   <- sqrt((n - 1) * mean((theta_i - theta_bar)^2))
    
    theta_bc  <- if (bias.correct) n * theta - (n - 1) * theta_bar else theta
    
    out <- list(est    = theta,
                est_bc = theta_bc,
                var    = se_jack^2,
                se     = se_jack)
    
    if (return.replicates) out$jack <- theta_i
    out
  }

## ================================================================
## 7. Jack-knife SE (and bias-correction) for correlation
## ================================================================
  #' @title Jackknife Estimation of Correlation
  #' @description Computes jackknife estimate and standard error for the Pearson correlation between two variables.
  #' @param dat A two-column data frame or matrix.
  #' @param bias.correct Logical, whether to apply bias correction. Default is TRUE.
  #' @param return.replicates Logical, whether to return jackknife replicates. Default is FALSE.
  #' @return A list with estimate (`est`), bias-corrected estimate (`est_bc`), variance, and standard error.
  #' @examples
  #' g2 <- MASS::mvrnorm(1000, mu = c(0, 0), Sigma = matrix(c(1, -0.4, -0.4, 1), 2))
  #' jack_cor(g2)
  jack_cor <- function(dat, bias.correct = TRUE,
                        return.replicates = FALSE) {
    
    n   <- nrow(dat)
    g2  <- function(z) r.to.zr(cor(z)[1,2])

    theta_i <- vapply(seq_len(n),
                      function(i) g2(dat[-i, ]),
                      numeric(1))
    
    theta     <- g2(dat)
    theta_bar <- mean(theta_i)
    
    se_jack   <- sqrt((n - 1) * mean((theta_i - theta_bar)^2))
    
    theta_bc  <- if (bias.correct) n * theta - (n - 1) * theta_bar else theta
    
    out <- list(est    = theta,
                est_bc = theta_bc,
                var    = se_jack^2,
                se     = se_jack)
    
    if (return.replicates) out$jack <- theta_i
    out
  }

## ================================================================
## 8. Non-parametric bootstrap for correlation
##    – returns a list:  point est, bias-corrected est, se, replicates*
## ================================================================
  #' @title Bootstrap Estimation of Correlation
  #' @description Computes bootstrap estimate and standard error for the Pearson correlation between two variables.
  #' @param dat A two-column data frame or matrix.
  #' @param B Number of bootstrap replicates. Default is 2000.
  #' @param bias.correct Logical, whether to apply bias correction. Default is TRUE.
  #' @param return.replicates Logical, whether to return replicate values. Default is FALSE.
  #' @return A list with estimate (`est`), bias-corrected estimate (`est_bc`), variance, and standard error.
  #' @examples
  #' g2 <- MASS::mvrnorm(1000, mu = c(0, 0), Sigma = matrix(c(1, -0.4, -0.4, 1), 2))
  #' x = boot_cor(g2, return.replicates = FALSE)
  boot_cor <- function(dat, B = 2000, bias.correct = TRUE,
                      return.replicates = FALSE) {
    g2 <- function(z) r.to.zr(cor(z)[1,2])

    b.reps <- replicate(B, {
      sample_rows <- dat[sample(nrow(dat), replace = TRUE), ]
      res <- g2(sample_rows)
      
      # Replace Inf/-Inf or NA with NA
      if (is.infinite(res) || is.na(res)) {res <- NA}
      res
    })

    est    <- g2(dat)
    est.bc <- if (bias.correct) 2 * est - mean(b.reps, na.rm = TRUE) else est
    out    <- list(est     = est,
                  est_bc  = est.bc,
                  var     = sd(b.reps, na.rm = TRUE)^2,
                  se      = sd(b.reps, na.rm = TRUE))
    if (return.replicates) {out$boot <- b.reps}
    out
  }

## ================================================================
## 9. Plotting Function for Overlaid Density Curves
##    – returns a ggplot object
## ================================================================
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

## ================================================================
## 10. Create Data Function
##    – generates random data from Pearson distributions
##    based on provided moments (mean, variance, skewness, kurtosis)
##    – returns a list of data frames and parameters
## ================================================================
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

## ================================================================
## 11. Plotting Function for Overlaid Density Curves
##    – returns a ggplot object
## ================================================================
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

## ================================================================
## 12. Plotting Function for Multiple Scenarios
##    – iterates over scenarios and saves a grid of plots
##    – returns a ggplot object
## ================================================================
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
    multi_plot_scenarios <- function(scenarios, folder = "figs/", name) {
      plots <- list() # Initialize an empty list to store plots
      for(i in 1:nrow(scenarios)) {
        params <- scenarios[i,]
        plots[[i]] <- plot_scenario(params, print = TRUE, xpos = 3)
      }
    # Combine all plots into a grid layout
    ggsave(filename = paste0(folder, "scenario_plots_", name, ".png"), plot = wrap_plots(plots, ncol = 2), width = 10, height = 10, dpi = 300)
}

## ================================================================
## 13. Fisher Z-transformation of Correlation Coefficient
##    – converts a Pearson correlation coefficient to Fisher's Zr
##    – returns the Zr-transformed correlation value(s)
## ================================================================
  #' @title Fisher Z-transformation of Correlation Coefficient
  #' @description Converts a Pearson correlation coefficient \eqn{r} to Fisher's Zr using the inverse hyperbolic tangent.
  #' @param r A numeric vector of Pearson correlation coefficients.
  #' @return The Zr-transformed correlation value(s).
  #' @examples
  #' r.to.zr(0.5)
  r.to.zr <- # Zr estimate
    function(r) { 
      0.5 * log((1 + r) / (1 - r))
  }

## ================================================================
## 14. Variance of Fisher Z-transformed Correlation
##    – computes the approximate sampling variance of Fisher's Zr given a sample size
##    – returns the variance of the Zr-transformed correlation
## ================================================================
  #' @title Variance of Fisher Z-transformed Correlation
  #' @description Computes the approximate sampling variance of Fisher's Zr given a sample size.
  #' @param n Sample size (must be greater than 3).
  #' @return The variance of the Zr-transformed correlation.
  #' @examples
  #' zr.variance(30)
  zr.variance <- # Zr variance 
    function(n) {
      1 / (n - 3)
  }
## ================================================================

## ================================================================
## 15. # Plot function for the results
##    – creates violin plots for bias in estimates
## ================================================================

plot_bias_violin <- 
  function(data, 
           y_var, 
           y_lab, 
           title = "", 
           ylim = c(-1, 
                    1)) {
    ggplot(data,
           aes(x = factor(n),
               y = .data[[y_var]],
               fill = factor(n))) + 
      ylim(ylim) +
      geom_violin() + 
      geom_hline(aes(yintercept = 0),
                 linetype = "dashed",
                 color = "black") +
      labs(x = "Sample size",
           y = TeX(y_lab),
           title = TeX(title)) +
      scale_fill_viridis_d() +
      theme_classic() + 
      theme(legend.position = "none",
            axis.title = element_text(size = 12),
            axis.text = element_text(size = 10),
            title = element_text(size = 6))
  }

## ================================================================
## 16. # BCa Bootstrap CI function for the results
##    – BCa bootstrap confidence intervals
#================================================================

#' @title Bootstrap Estimation of Excess Kurtosis using BCa Method
#' @description Computes bootstrap estimate and BCa CIs of excess kurtosis.
#' @param data A data frame containing the variables x1, and x2 (column names must match) which is the simulated excess kurtosis for each group.
#' @param B Number of bootstrap replicates. Default is 2000.
#' @param return.replicates Logical, whether to return replicate values. Default is FALSE.
#' @return A list with estimate (`est`), bias-corrected estimate (`est_bc`), variance, and standard error.
#' @examples
#' x1 <- rgamma(25, shape = 5)
#' x2 <- rgamma(25, shape = 5)
#' data <- data.frame(x1 = x1, x2 = x2)
#' t  <- boot_kurt_bca(data)
  # x1 <- tryCatch(rpearson(1000, 
	#  			 moments = c(    mean = 0, 
	#  						 variance = 1, 
	#  						 skewness = 0, 
	#  						 kurtosis = 6)),   
  #                 error = function(e) {
  #                   message("Error in rpearson for group 1: ", e)
  #                   return(NA)
  #                 })
  # x2 <- tryCatch(rpearson(1000, 
	#  			 moments = c(    mean = 0, 
	#  						 variance = 1, 
	#  						 skewness = 0, 
	#  						 kurtosis = 3)),   
  #                 error = function(e) {
  #                   message("Error in rpearson for group 1: ", e)
  #                   return(NA)
  #                 })
  #       data <- data.frame(x1 = x1, x2 = x2)
  #  t  <- boot_kurt_bca(data)

  boot_kurt_bca <- function(data, B = 2000, return.replicates = FALSE) {

    # Excess kurtosis statistic function
    g2 <- function(x, i) {
      n1 = length(x$x1[i])
      n2 = length(x$x2[i])
      (((((n1 + 1) * n1 * (n1 - 1)) / ((n1 - 2) * (n1 - 3))) *
       (sum((x$x1[i] - mean(x$x1[i])) ^ 4) / (sum((x$x1[i] - mean(x$x1[i])) ^ 2) ^ 2))) -(3 * ((n1 - 1) ^ 2) / ((n1 - 2) * (n1 - 3)))) -  
       (((((n2 + 1) * n2 * (n2 - 1)) / ((n2 - 2) * (n2 - 3))) *
       (sum((x$x2[i] - mean(x$x2[i])) ^ 4) / (sum((x$x2[i] - mean(x$x2[i])) ^ 2) ^ 2))) -(3 * ((n2 - 1) ^ 2) / ((n2 - 2) * (n2 - 3))))
      }

    # Bootstrap resampling, non-parametric
      b.reps <- boot::boot(data, statistic = g2, R = B)

    # Point estimate
    est    <- b.reps$t0
    
    # Bias corrected estimate
    est.bc <-  2 * est - mean(b.reps$t)

    # Bias corrected accelerated (BCa) bootstrap CIs
    boot_bca <- boot::boot.ci(b.reps, type = "bca")

    out    <- list(est     = est,
                  est_bc   = est.bc,
                  est_ci_l = boot_bca$bca[4],
                  est_ci_u = boot_bca$bca[5],
                  var      = sd(b.reps$t) ^ 2,
                  se       = sd(b.reps$t))

    if (return.replicates) {out$boot <- b.reps}
    
    return(out)
}
