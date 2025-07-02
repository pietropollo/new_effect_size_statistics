###------------------------------------------------------------------------###
# Simulate Pearson distribution with specified moments
# Description: Equal sample size between groups. 
###------------------------------------------------------------------------###

# Load required packages
	pacman::p_load(moments, PearsonDS)

# Parameters
	          n <- c(10, 25, 50, 75, 100, 150) 
	    scenarios  <- data.frame(mean_g1 = c(0, 5, 10, 15, 20),
	    mean_g2 = c(0, 5, 10, 15, 20),
	variance_g1 = c(1, 2, 3, 4, 5),
	variance_g2 = c(1, 2, 3, 4, 5),
	skewness_g1 = c(1, 1.5, 2, 3, 4),
	skewness_g2 = c(1, 1.5, 2, 3, 4),
	kurtosis_g1 = c(3, 4, 5, 6, 7),
	kurtosis_g2 = c(3, 4, 5, 6, 7))
          nsims <- 5000  # Number of simulations

# Create combinations of parameters expanded by sample size vector. Each row is a scenario with a sample size which is used to set up the simulation
	params_all <- data.frame(tidyr::crossing(scenarios, n = n))

# functions for calculating effect sizes ----
## skewness ----
calc.skewness <- function(x, output = "est") {
  n <- length(x)
  
  if (output == "est") { # skewness estimate
    (sqrt(n * (n - 1)) / (n - 2)) *
      (((1 / n) * sum((x - mean(x)) ^ 3)) /
         (((1 / n) * sum((x - mean(x)) ^ 2)) ^ (3/2)))
    
  } else if (output == "var") { # skewness sampling variance
    (6 * n * (n - 1)) /
      ((n - 2) * (n + 1) * (n + 3))
  }
}

## kurtosis ----
calc.kurtosis <- function(x, output = "est") {
  n <- length(x)
  
  if (output == "est") { # kurtosis estimate
    ((((n + 1) * n * (n - 1)) / ((n - 2) * (n - 3))) *
       (sum((x - mean(x)) ^ 4) / (sum((x - mean(x)) ^ 2) ^ 2))) -
      (3 * ((n - 1) ^ 2) / ((n - 2) * (n - 3)))
  } else if (output == "var") { # kurtosis sampling variance
    (24 * n * ((n - 1) ^ 2)) /
      ((n - 3) * (n - 2) * (n + 3) * (n + 5))
  }
}


#| @title Simulate Pearson distribution with equal sample sizes 
#| @description Function to simulate Pearson distribution with equal sample sizes
#| and calculate bias and coverage
#| @param n Sample size for each group
#| @param params Data frame containing parameters for the simulation. Note that this is from a params data frame created by the expand.grid function above and so included data from a single row which is a single scenario.
#| @param nsims Number of simulations to run (default is 1000)
#| @return A data frame with bias and coverage for each simulation
#| @examples
#| sim_func(n = 50, params = params[i,], nsims = 1000)

sim_func <- function(n, params, nsims = nsims) {
  # Vectors to store results
      bias_sk <- numeric(nsims)
	  bias_ku <- numeric(nsims)
     coverage <- numeric(nsims)
  
for(i in 1:nsims) {
##---------------------------##
  # Simulate data for group 1
##---------------------------##
  x1 <- tryCatch(rpearson(n, 
				 moments = c(    mean = params$mean_g1, 
							 variance = params$variance_g1, 
							 skewness = params$skewness_g1, 
							 kurtosis = params$kurtosis_g1)),
                 error = function(e) {
                   message("Error in rpearson for group 1: ", e)
                   return(NA)
                 })
				 
##---------------------------##
  # Simulate data for group 2
##---------------------------##
  x2 <- tryCatch(rpearson(n, 
				 moments = c(    mean = params$mean_g2, 
							 variance = params$variance_g2, 
							 skewness = params$skewness_g2, 
							 kurtosis = params$kurtosis_g2)),
                 error = function(e) {
                   message("Error in rpearson for group 2: ", e)
                   return(NA)
                 })

##-------------------------------------------------##
  # Calculate bias and coverage for effect statistics
##-------------------------------------------------##
## BELOW IS JUST AN EXAMPLE WE NEED TO MODIFY FOR VARIOUS EFFECTS
# Calculate bias between groups, which is how much the difference between estimated sk and kur deviates from the true differences
  bias_sk[i] = (calc.skewness(x1) - calc.skewness(x2)) - (params$skewness_g1 - params$skewness_g2)
  bias_ku[i] = (calc.kurtosis(x1) - calc.kurtosis(x2)) - (params$kurtosis_g1 - params$kurtosis_g2)

# Calculate coverage, which is the proportion of times the confidence interval contains the true difference
           ci <- tryCatch(t.test(x1, x2)$conf.int, error = function(e) {
							message("Error in t.test: ", e)
							return(NA)
							})
  coverage[i] <- ifelse(ci[1] <= (params$mean_g1 - params$mean_g2) && 
					    ci[2] >= (params$mean_g1 - params$mean_g2), 1, 0)

}
  
##-------------------------------------------------##
  # Return data with all the simulation results
##-------------------------------------------------##
  return(data.frame(bias = mean(bias1, na.rm = TRUE), 
  				coverage = sum(coverage, na.rm = TRUE) / nsims,
					   n = length(bias1)))
}

###------------------------------------------------------------------------###
# Run simulations for all scenarios assuming equal sample size 
###------------------------------------------------------------------------###
# Initialize an empty data frame to store results
result <- data.frame()

# Loop through each scenario and run the simulation function
for(i in 1:nrow(params_all)) {
  params <- params_all[i,]
  result <- rbind(result, sim_func(n = params$n, params = params, nsims = nsims))
  print(paste("Simulation for scenario", i, "completed. Bias:", round(result$bias, 2), "Coverage:", round(result$coverage, 2)))
}
