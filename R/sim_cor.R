###------------------------------------------------------------------------###
# Simulate correlation-based effect sizes differences between two groups
# Description: Equal sample size between groups. 
###------------------------------------------------------------------------###

# Clean working environment
rm(list = ls()) # Remove all objects from the environment

# Load required packages
	pacman::p_load(moments, PearsonDS, MASS)
   source("./R/func.R") # Load the functions from func.R


  ## Zr ----
r.to.zr <- # Zr estimate
  function(r) { 
    0.5 * log((1 + r) / (1 - r))
  }

zr.variance <- # Zr variance 
  function(n) {
    1 / (n - 3)
  }



sim_cor <- function(params, nsims = nsims) {
  # Vectors to store results
      bias_d_cor <- numeric(nsims)
   bias_d_cor_sv <- numeric(nsims)
  
for(i in 1:nsims) {
##---------------------------##
  # Simulate data for group 1
##---------------------------##
  g1 <- tryCatch(MASS::mvrnorm(params$n, 
				               mu = c(0, 0), 
				            Sigma = matrix(c(1, params$cor_g1, params$cor_g1, 1), nrow = 2)),
                    error = function(e) {
                        		message("Error in mvrnorm for group 1: ", e)
                        		return(NA)
                      })

##---------------------------##
  # Simulate data for group 2
##---------------------------##
  g2 <- tryCatch(MASS::mvrnorm(params$n, 
				               mu = c(0, 0), 
				            Sigma = matrix(c(1, params$cor_g2, params$cor_g2, 1), nrow = 2)),
                    error = function(e) {
						         message("Error in mvrnorm for group 2: ", e)
						         return(NA)
					  })

##-------------------------------------------------##
  # Calculate bias 
##-------------------------------------------------##
# Calculate bias between groups, which is how much the difference between correlations deviates from the true difference in values
     bias_d_cor[i] = (r.to.zr(cor(g1)[1,2]) - r.to.zr(cor(g2)[1,2])) - (r.to.zr(params$cor_g1) - r.to.zr(params$cor_g2))
  bias_d_cor_sv[i] = (1 / (params$n - 3)) + (1 / (params$n - 3)) # sampling variance of the difference in correlations
}
  
##-------------------------------------------------##
  # Return data with all the simulation results
##-------------------------------------------------##
  return(data.frame(   bias_d_cor = mean(bias_d_cor, na.rm = TRUE),
                    bias_d_cor_sv = ((mean(bias_d_cor_sv) - sd(bias_d_cor)^2) / sd(bias_d_cor)^2)*100,
                     mcse_bias_sv = sqrt(var(bias_d_cor_sv, na.rm = TRUE) / nsims),
                        mcse_bias = sqrt(var(bias_d_cor, na.rm = TRUE) / nsims),
					                  nsims = length(bias_d_cor)))
}

###------------------------------------------------------------------------###
# Run simulations for all scenarios assuming equal sample size 
###------------------------------------------------------------------------###

# Set seed for reproducibility 
  set.seed(860)

# Parameters
	                n <- c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 150, 500) 
	    	     cor_g1 <- c(-0.8, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8)
             cor_g2 <- c(-0.8, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8)
              nsims <- 5000  # Number of simulations

# Create combinations of parameters expanded by sample size vector. Each row is a scenario with a sample size which is used to set up the simulation
	       scenarios <- expand.grid(cor_g1 = cor_g1, cor_g2 = cor_g2)
         scenarios <- scenarios[!duplicated(scenarios), ] # Remove duplicate rows if any
scenarios$scenario <- 1:nrow(scenarios) # Add a scenario number column

# Params for the simulation function
	params_all <- data.frame(tidyr::crossing(scenarios, n = n))

# Initialize an empty data frame to store results
results_d_cor <- data.frame()

# Loop through each scenario and run the simulation function
system.time(
for(i in 1:nrow(params_all)) {
         params <- params_all[i,]
  results_d_cor <- rbind(results_d_cor, sim_cor(params = params, nsims = nsims))
  
  print(paste("Simulation for scenario", i, "completed. Bias:", round(results_d_cor$bias_d_cor, 2), "Bias SV:", round(results_d_cor$bias_d_cor_sv, 2), "MCSE Bias SV:", round(results_d_cor$mcse_bias_sv, 2), "MCSE Bias:", round(results_d_cor$mcse_bias, 2)))
}
) # 335.483 elapsed time which is 6 minutes 

# cbind results with parameters
result <- cbind(params_all, results_d_cor)

#Save the results to a file
saveRDS(result, file = "./output/result_cor.rds") # Save the results to a file