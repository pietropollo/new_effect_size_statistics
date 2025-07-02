###------------------------------------------------------------------------###
# Simulate Pearson distribution with specified moments
# Description: Equal sample size between groups. 
###------------------------------------------------------------------------###

# Load required packages
	pacman::p_load(moments, PearsonDS)

# Parameters
	                                      n = c(10, 25, 50, 75, 100, 150) 
	    							mean_g1 = c(0, 0, 0)
									mean_g2 = c(0, 10, 20)
								variance_g1 = c(2)
								variance_g2 = c(2)
								skewness_g1 = c(1, 1.5, 2, 3, 4)
								skewness_g2 = c(1, 1.5, 2, 3, 4)
								kurtosis_g1 = c(3, 6)
								kurtosis_g2 = c(3, 6)
		scenarios <- expand.grid(mean_g1 = mean_g1, mean_g2 = mean_g2, variance_g1 = variance_g1, variance_g2 = variance_g2, skewness_g1 = skewness_g1, skewness_g2 = skewness_g2, kurtosis_g1 = kurtosis_g1, kurtosis_g2 = kurtosis_g2) # Create all combinations of scenarios
          nsims <- 1000  # Number of simulations

# Create combinations of parameters expanded by sample size vector. Each row is a scenario with a sample size which is used to set up the simulation
	params_all <- data.frame(tidyr::crossing(scenarios, n = n))

# functions for calculating effect sizes ----
## skewness ---- !!!!# Need to check because this function does not match `moments::skewness`
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

## kurtosis ---- !!!!# Need to check because this function does not match `moments::kurtosis`
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


#| @title Simulate Pearson distribution with specified moments
#| @description Function to simulate Pearson distribution with defined moments: mean, variance, skewness and kurtosis
#| and calculate bias and coverage
#| @param params Data frame containing parameters for the simulation. Note that this is from a params data frame created by the expand.grid function above and so included data from a single row which is a single scenario. Params must define the sample size and the moments for both groups (mean_g1, mean_g2, variance_g1, variance_g2, skewness_g1, skewness_g2, kurtosis_g1, kurtosis_g2).
#| @param nsims Number of simulations to run (default is 1000)
#| @return A data frame with bias with respect to point estimates of effect size and its sampling variance for each scenario.
#| @examples
#| sim_func(n = 50, params = params[i,], nsims = 1000)

sim_func <- function(params, nsims = nsims) {
  # Vectors to store results
      sk <- numeric(nsims)
	  ku <- numeric(nsims)

for(i in 1:nsims) {
##---------------------------##
  # Simulate data for group 1
##---------------------------##
  x1 <- tryCatch(rpearson(params$n, 
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
  x2 <- tryCatch(rpearson(params$n, 
				 moments = c(    mean = params$mean_g2, 
							 variance = params$variance_g2, 
							 skewness = params$skewness_g2, 
							 kurtosis = params$kurtosis_g2)),
                 error = function(e) {
                   message("Error in rpearson for group 2: ", e)
                   return(NA)
                 })

##-------------------------------------------------##
  # Calculate effect statistics
##-------------------------------------------------##

#Calculate and store the effect statistics for skewness and kurtosis
  sk[i] = tryCatch((calc.skewness(x1) - calc.skewness(x2)), error = function(e) {return(NA)}) #!!! CHECK THIS FUNCTION does not match `moments::skewness`
  ku[i] = tryCatch((calc.kurtosis(x1) - calc.kurtosis(x2)), error = function(e) {return(NA)}) #!!! CHECK THIS FUNCTION does not match `moments::kurtosis`
}
  
##-------------------------------------------------##
  # Return data with all the simulation results
##-------------------------------------------------##
  return(data.frame(bias_sk = mean(sk) - (params$skewness_g1 - params$skewness_g2), # Bias for skewness from true value
  				    bias_ku = mean(ku) - (params$kurtosis_g1 - params$kurtosis_g2), # Bias for kurtosis from true value
				 bias_sk_sv = sd(sk)^2 - (calc.skewness(x1, output = "var") + calc.skewness(x2, output = "var")), # Bias for skewness sampling variance from analytical approximation
				 bias_ku_sv = sd(ku)^2 - (calc.kurtosis(x1, output = "var") + calc.kurtosis(x2, output = "var")), # Bias for kurtosis sampling variance from analytical approximation
					 n_sims = nsims)) 
}

###------------------------------------------------------------------------###
# Run simulations for all scenarios assuming equal sample size 
###------------------------------------------------------------------------###
# Initialize an empty data frame to store results
result <- data.frame()

# Loop through each scenario and run the simulation function
for(i in 1:nrow(params_all)) {
  params <- params_all[i,]
  result <- rbind(result, sim_func(params = params, nsims = nsims))
  print(paste("Simulation for scenario", i, "completed. Bias_sk:", round(result$bias_sk[i], 2), "Bias_ku:", round(result$bias_ku[i], 2)))
}
