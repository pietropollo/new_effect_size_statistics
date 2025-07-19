###------------------------------------------------------------------------###
# Simulate correlation-based effect sizes differences between two groups
# Description: Equal sample size between groups. 
###------------------------------------------------------------------------###

# Load required packages
	pacman::p_load(moments, PearsonDS, MASS)

# Parameters
	                n <- c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 150, 500) 
	    	     cor_g1 <- c(-0.8, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8)
             cor_g2 <- c(-0.8, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8)
             nsims <- 5000  # Number of simulations

# Create combinations of parameters expanded by sample size vector. Each row is a scenario with a sample size which is used to set up the simulation
	params_all <- expand.grid(n = n, cor_g1 = cor_g1, cor_g2 = cor_g2)

  ## Zr ----
r.to.zr <- # Zr estimate
  function(r) { 
    0.5 * log((1 + r) / (1 - r))
  }

zr.variance <- # Zr variance 
  function(n) {
    1 / (n - 3)
  }



sim_cor <- function(n, params, nsims = nsims) {
  # Vectors to store results
      bias1 <- numeric(nsims)
   coverage <- numeric(nsims)
  
for(i in 1:nsims) {
##---------------------------##
  # Simulate data for group 1
##---------------------------##
  g1 <- tryCatch(MASS::mvrnorm(n, 
				               mu = c(0, 0), 
				            Sigma = matrix(c(1, params$cor_g1, params$cor_g1, 1), nrow = 2)),
                            error = function(e) {
                        		message("Error in mvrnorm for group 1: ", e)
                        		return(NA)
                      })

##---------------------------##
  # Simulate data for group 2
##---------------------------##
  g2 <- tryCatch(MASS::mvrnorm(n, 
				               mu = c(0, 0), 
				            Sigma = matrix(c(1, params$cor_g2, params$cor_g2, 1), nrow = 2)),
                            error = function(e) {
						         message("Error in mvrnorm for group 2: ", e)
						         return(NA)
					  })

##-------------------------------------------------##
  # Calculate bias and coverage for effect statistics
##-------------------------------------------------##
## BELOW IS JUST AN EXAMPLE WE NEED TO MODIFY FOR VARIOUS EFFECTS
# Calculate bias between groups, which is how much the difference between correlations deviates from the true difference in values
  bias1[i] = (cor(g1)[1,2] - cor(g2)[1,2]) - (params$cor_g1 - params$cor_g2)

}
  
##-------------------------------------------------##
  # Return data with all the simulation results
##-------------------------------------------------##
  return(data.frame(bias = mean(bias1, na.rm = TRUE), 
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
  result <- rbind(result, sim_cor(n = params$n, params = params, nsims = nsims))
  print(paste("Simulation for scenario", i, "completed. Bias:", round(result$bias, 2)))
}

# cbind results with parameters
result <- cbind(params_all, result)