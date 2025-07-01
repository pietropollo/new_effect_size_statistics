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

#| @title Simulate Pearson distribution with equal sample sizes 
#| @description Function to simulate Pearson distribution with equal sample sizes
#| and calculate bias and coverage
#| @param n Sample size for each group
#| @param params Data frame containing parameters for the simulation. Note that this is from a params data frame created by the expand.grid function above and so included data from a single row which is a single scenerio.
#| @param nsims Number of simulations to run (default is 1000)
#| @return A data frame with bias and coverage for each simulation
#| @examples
#| sim_equal(n = 50, params = params[i,], nsims = 1000)

sim_equal <- function(n, params, nsims = nsims) {
  # Vectos to store results
      bias1 <- numeric(nsims)
  coverage <- numeric(nsims)
  
  for(i in 1:nsims) {
  # Simulate data for group 1
  x1 <- tryCatch(rpearson(n, 
				 moments = c(    mean = params$mean_g1, 
							 variance = params$variance_g1, 
							 skewness = params$skewness_g1, 
							 kurtosis = params$kurtosis_g1)),
                 error = function(e) {
                   message("Error in rpearson for group 1: ", e)
                   return(NA)
                 })

  # Simulate data for group 2
  x2 <- tryCatch(rpearson(n, 
				 moments = c(    mean = params$mean_g2, 
							 variance = params$variance_g2, 
							 skewness = params$skewness_g2, 
							 kurtosis = params$kurtosis_g2)),
                 error = function(e) {
                   message("Error in rpearson for group 2: ", e)
                   return(NA)
                 })

# Calculate bias between groups, which is how much the difference between means deviates from the true difference
  bias1[i] = (mean(x1) - mean(x2)) - (params$mean_g1 + params$mean_g2)  

# Calculate coverage, which is the proportion of times the confidence interval contains the true difference
           ci <- tryCatch(t.test(x1, x2)$conf.int, error = function(e) {
							message("Error in t.test: ", e)
							return(NA)
							})
  coverage[i] <- ifelse(ci[1] <= (params$mean_g1 - params$mean_g2) && 
					    ci[2] >= (params$mean_g1 - params$mean_g2), 1, 0)

  }
  
  # Return a data frame with bias and coverage
  return(data.frame(bias = mean(bias1, na.rm = TRUE), 
  				coverage = sum(coverage, na.rm = TRUE) / nsims,
					   n = length(bias1)))
}

params <- params_all[1,]  # Select the first row of parameters for testing
system.time(
t <- sim_equal(n = params[5,"n"], params = params[5,], nsims = 5000)
)

# Run simulations for all scenarios assuming equal sample size 
result <- data.frame()
for(i in 1:nrow(params_all)) {
  params <- params_all[i,]
  result <- rbind(result, sim_equal(n = params$n, params = params, nsims = nsims))
  print(paste("Simulation for scenario", i, "completed. Bias:", round(result$bias, 2), "Coverage:", round(result$coverage, 2)))
}
