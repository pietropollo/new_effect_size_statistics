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

  # Point estimates
   jack_d_cor_bc <- numeric(nsims)
   jack_d_cor_sv <- numeric(nsims)

  # Coverage indicators
  coverage_d_cor <- numeric(nsims)
  coverage_d_cor_jack_bc <- numeric(nsims)  
  coverage_jack_bc_sv <- numeric(nsims)
  coverage_bc_jack_sv <- numeric(nsims)

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
  g1_cor_jack <- tryCatch(jack_cor(g1), error = function(e) {return(NA)})
  g2_cor_jack <- tryCatch(jack_cor(g2), error = function(e) {return(NA)})

  ## Calculate the bias-corrected estimates and sampling variances from bootstraps and jackknife
       jack_d_cor_bc[i] <- (g1_cor_jack$est_bc - g2_cor_jack$est_bc) 
       jack_d_cor_sv[i] <- (g1_cor_jack$var + g2_cor_jack$var)

  # Coverage indicators
  coverage_d_cor[i] <- ((r.to.zr(params$cor_g1) - r.to.zr(params$cor_g2)) >= (d_cor[i] - qt(0.975, df = params$n - 2) * sqrt(d_cor_sv[i])) & 
                        (r.to.zr(params$cor_g1) - r.to.zr(params$cor_g2)) <= (d_cor[i] + qt(0.975, df = params$n - 2) * sqrt(d_cor_sv[i])))
  coverage_d_cor_jack_bc[i] <- ((r.to.zr(params$cor_g1) - r.to.zr(params$cor_g2)) >= (jack_d_cor_bc[i] - qt(0.975, df = params$n - 2) * sqrt(jack_d_cor_sv[i])) & 
                                (r.to.zr(params$cor_g1) - r.to.zr(params$cor_g2)) <= (jack_d_cor_bc[i] + qt(0.975, df = params$n - 2) * sqrt(jack_d_cor_sv[i])))
  coverage_jack_bc_sv[i] <- ((r.to.zr(params$cor_g1) - r.to.zr(params$cor_g2)) >= (jack_d_cor_bc[i] - qt(0.975, df = params$n - 2) * sqrt(d_cor_sv[i])) & 
                                (r.to.zr(params$cor_g1) - r.to.zr(params$cor_g2)) <= (jack_d_cor_bc[i] + qt(0.975, df = params$n - 2) * sqrt(d_cor_sv[i])))
  coverage_bc_jack_sv[i] <- ((r.to.zr(params$cor_g1) - r.to.zr(params$cor_g2)) >= (d_cor[i] - qt(0.975, df = params$n - 2) * sqrt(jack_d_cor_sv[i])) & 
                                (r.to.zr(params$cor_g1) - r.to.zr(params$cor_g2)) <= (d_cor[i] + qt(0.975, df = params$n - 2) * sqrt(jack_d_cor_sv[i])))                                                            
}

##-------------------------------------------------##
# Return data with all the simulation results
# Calculate bias, relative bias and MCSE for all methods
##-------------------------------------------------##
  return(data.frame(          
         # Bias
         bias_d_cor = (mean(d_cor, na.rm = TRUE) - (r.to.zr(params$cor_g1) - r.to.zr(params$cor_g2))),
    bias_jack_d_cor = (mean(jack_d_cor_bc, na.rm = TRUE) - (r.to.zr(params$cor_g1) - r.to.zr(params$cor_g2))),

         # Relative bias
      bias_d_cor_sv = ((mean(d_cor_sv) - sd(d_cor)^2) / sd(d_cor)^2)*100,
 bias_jack_d_cor_sv = ((mean(jack_d_cor_sv, na.rm = TRUE) - sd(jack_d_cor_bc, na.rm = TRUE)^2) / sd(jack_d_cor_bc, na.rm = TRUE)^2)*100,
 bias_d_cor_jack_sv = ((mean(jack_d_cor_sv, na.rm = TRUE) - sd(d_cor, na.rm = TRUE)^2) / sd(d_cor, na.rm = TRUE)^2)*100,
 bias_jack_d_cor_sv = ((mean(d_cor_sv, na.rm = TRUE) - sd(jack_d_cor_bc, na.rm = TRUE)^2) / sd(jack_d_cor_bc, na.rm = TRUE)^2)*100,

      # Coverage
              coverage_d_cor = sum(coverage_d_cor) / nsims,
      coverage_d_cor_jack_bc = sum(coverage_d_cor_jack_bc) / nsims,
      coverage_jack_bc_sv = sum(coverage_jack_bc_sv) / nsims,
      coverage_bc_jack_sv = sum(coverage_bc_jack_sv) / nsims,

         # Monte Carlo Error
 mcse_jack_d_cor_sv = sqrt(var(jack_d_cor_sv, na.rm = TRUE) / nsims),
    mcse_jack_d_cor = sqrt(var(jack_d_cor_bc, na.rm = TRUE) / nsims),
          mcse_bias = sqrt(var(d_cor, na.rm = TRUE) / nsims),
       mcse_bias_sv = sqrt(var(d_cor_sv, na.rm = TRUE) / nsims),
       mcse_coverage_d_cor = sqrt((coverage_d_cor * (1 - coverage_d_cor)) / nsims),
mcse_coverage_d_cor_jack_bc = sqrt((mean(coverage_d_cor_jack_bc) * (1 - mean(coverage_d_cor_jack_bc))) / nsims),
   mcse_coverage_jack_bc_sv = sqrt((mean(coverage_jack_bc_sv) * (1 - mean(coverage_jack_bc_sv))) / nsims),
   mcse_coverage_bc_jack_sv = sqrt((mean(coverage_bc_jack_sv) * (1 - mean(coverage_bc_jack_sv))) / nsims),
              n_sim = length(d_cor)))
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
) # 1.36 hours

# cbind results with parameters
result <- cbind(params_all, results_d_cor)

#Save the results to a file
saveRDS(result, file = "./output/result_cor.rds") # Save the results to a file
