###------------------------------------------------------------------------###
# Simulate correlation-based effect sizes differences between two groups
# Description: Equal sample size between groups. 
###------------------------------------------------------------------------###

# Clean working environment
rm(list = ls()) # Remove all objects from the environment

# Load required packages
	pacman::p_load(moments, PearsonDS, MASS)
   source("./R/func.R") # Load the functions from func.R


# Core function for simulating correlation-based effect sizes differences
sim_cor <- function(params, nsims = nsims) {
  # Vectors to store results
      d_cor <- numeric(nsims)
   d_cor_sv <- numeric(nsims)

   boot_d_cor_bc <- numeric(nsims)
   boot_d_cor_sv <- numeric(nsims)
   jack_d_cor_bc <- numeric(nsims)
   jack_d_cor_sv <- numeric(nsims)

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
  # Calculate effects 
##-------------------------------------------------##
# Calculate the difference between correlations and combined sampling variance
     d_cor[i] = (r.to.zr(cor(g1)[1,2]) - r.to.zr(cor(g2)[1,2])) 
  d_cor_sv[i] = (1 / (params$n - 3)) + (1 / (params$n - 3)) # sampling variance of the difference in correlations

  ## Add jacknife and bootstrap methods
  g1_cor_boot <- tryCatch(boot_cor(g1), error = function(e) {return(NA)})
  g2_cor_boot <- tryCatch(boot_cor(g2), error = function(e) {return(NA)})
  g1_cor_jack <- tryCatch(jack_cor(g1), error = function(e) {return(NA)})
  g2_cor_jack <- tryCatch(jack_cor(g2), error = function(e) {return(NA)})

  ## Calculate the bias-corrected estimates and sampling variances from bootstraps and jackknife
       boot_d_cor_bc[i] <- (g1_cor_boot$est_bc - g2_cor_boot$est_bc)
       boot_d_cor_sv[i] <- (g1_cor_boot$var + g2_cor_boot$var)
  
       jack_d_cor_bc[i] <- (g1_cor_jack$est_bc - g2_cor_jack$est_bc) 
       jack_d_cor_sv[i] <- (g1_cor_jack$var + g2_cor_jack$var)
}
  
##-------------------------------------------------##
# Return data with all the simulation results
# Calculate bias, relative bias and MCSE for all methods
##-------------------------------------------------##
  return(data.frame(          bias_d_cor = (mean(d_cor, na.rm = TRUE) - (r.to.zr(params$cor_g1) - r.to.zr(params$cor_g2))),
                           bias_d_cor_sv = ((mean(d_cor_sv) - sd(d_cor)^2) / sd(d_cor)^2)*100,
                            mcse_bias_sv = sqrt(var(d_cor_sv, na.rm = TRUE) / nsims),
                               mcse_bias = sqrt(var(d_cor, na.rm = TRUE) / nsims),

                           bias_boot_d_cor = (mean(boot_d_cor_bc, na.rm = TRUE) - (r.to.zr(params$cor_g1) - r.to.zr(params$cor_g2))),
                        bias_boot_d_cor_sv = ((mean(boot_d_cor_sv) - sd(d_cor)^2) / sd(d_cor)^2)*100,
                        mcse_boot_d_cor_sv = sqrt(var(boot_d_cor_sv, na.rm = TRUE) / nsims),
                           mcse_boot_d_cor = sqrt(var(boot_d_cor_bc, na.rm = TRUE) / nsims),

                           bias_jack_d_cor = (mean(jack_d_cor_bc, na.rm = TRUE) - (r.to.zr(params$cor_g1) - r.to.zr(params$cor_g2))),
                        bias_jack_d_cor_sv = ((mean(jack_d_cor_sv) - sd(d_cor)^2) / sd(d_cor)^2)*100,
                        mcse_jack_d_cor_sv = sqrt(var(jack_d_cor_sv, na.rm = TRUE) / nsims),
                           mcse_jack_d_cor = sqrt(var(jack_d_cor_bc, na.rm = TRUE) / nsims),
                                         n = length(d_cor)))
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
              nsims <- 2500  # Number of simulations

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

  print(paste("Simulation for scenario", i, "completed. Bias:", round(results_d_cor$bias_d_cor[i], 2), "Bias SV:", round(results_d_cor$bias_d_cor_sv[i], 2), "MCSE Bias SV:", round(results_d_cor$mcse_bias_sv[i], 2), "MCSE Bias:", round(results_d_cor$mcse_bias[i], 2)))
}
) # 335.483 elapsed time which is 6 minutes 

# cbind results with parameters
result <- cbind(params_all, results_d_cor)

#Save the results to a file
saveRDS(result, file = "./output/result_cor.rds") # Save the results to a file