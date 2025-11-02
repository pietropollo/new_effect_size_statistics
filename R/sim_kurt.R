###------------------------------------------------------------------------###
# Simulate Pearson distribution with specified moments
# Description: Equal sample size between groups. 
###------------------------------------------------------------------------###

# Clear the environment
rm(list = ls())

# Load required packages
	pacman::p_load(moments, PearsonDS, tidyverse, patchwork, boot)
  source("./R/func.R") # Load the functions from func.R
  

###------------------------------------------------------------------------###
## Key Functions ----
###------------------------------------------------------------------------###

calc.kurtosis <- function(x, output = "est") {
  n <- length(x)
  
  if (output == "est") { # kurtosis estimate
    ((((n + 1) * n * (n - 1)) / ((n - 2) * (n - 3))) *
       (sum((x - mean(x)) ^ 4) / (sum((x - mean(x)) ^ 2) ^ 2))) -(3 * ((n - 1) ^ 2) / ((n - 2) * (n - 3))) # 
  } else if (output == "var") { # kurtosis sampling variance
    (24 * n * ((n - 1) ^ 2)) /
      ((n - 3) * (n - 2) * (n + 3) * (n + 5))
  }
}

kurtosis_uncorrected <- function(x) {
  n <- length(x)
  x_bar <- mean(x)
  s <- sd(x)  # sample standard deviation uses n - 1
  m4 <- mean((x - x_bar)^4)
  kurt <- m4 / s^4
  return(kurt)
}

#| @title Simulate Pearson distribution with specified moments 
#| @description Function to simulate Pearson distribution with defined moments: mean, variance, skewness and kurtosis
#| and calculate bias and coverage
#| @param params Data frame containing parameters for the simulation. Note that this is from a params data frame created by the expand.grid function above and so included data from a single row which is a single scenario. Params must define the sample size and the moments for both groups (mean_g1, mean_g2, variance_g1, variance_g2, skewness_g1, skewness_g2, kurtosis_g1, kurtosis_g2).
#| @param nsims Number of simulations to run (default is 1000)
#| @return A data frame with bias with respect to point estimates of effect size and its sampling variance for each scenario.
#| @examples
#| sim_kurt(n = 50, params = params[i,], nsims = 1000)

sim_kurt <- function(params, nsims = nsims) {
  
  # Vectors to store results. Note that we have a bunch of new functions to calculate the effect sizes and sampling variances which don't assume normality
 # Point estimates for kurtosis
  jack_kurt_bc <- numeric(nsims) # Jackknife bias-corrected method for kurtosis
         ku <- numeric(nsims)
   
  # Sampling variance for kurtosis
              ku_sv <- numeric(nsims) # Sampling variance for kurtosis
      jack_kurt_sv <- numeric(nsims) # Jackknife method sampling variance for kurtosis
   
  # Coverage indicators
              coverage_ku <- numeric(nsims) # Coverage for kurtosis effect size
			  coverage_ku_jack_bc <- numeric(nsims) # Coverage for kurtosis with jackknife bias-corrected method
			  coverage_ku_jack_sv <- numeric(nsims) # Coverage for kurtosis with jackknife bias-corrected method and sampling variance for kurtosis
			  coverage_ku_jack_adj_sv <- numeric(nsims) # Coverage for kurtosis with jackknife bias-corrected method and adjusted jackknife sampling variance for kurtosis
       coverage_jack_ku_sv <- numeric(nsims) # Coverage for kurtosis with jackknife bias-corrected method and sampling variance for kurtosis
    
    true_ku_diff <- (params$kurtosis_g1 - 3) - (params$kurtosis_g2 - 3) # True difference in EXCESS kurtosis between groups

    coverage <- function(est, sv, n, true_diff) {
      true_diff >= (est - qt(0.975, n-2) * sqrt(sv)) && true_diff <= (est + qt(0.975, n-2) * sqrt(sv))
    }
  ##-------------------------------------------------##

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

      # Do the computations and save to objects

    x1_kurt_jack <- tryCatch(jack_kurt(x1), error = function(e) {return(NA)})
    x2_kurt_jack <- tryCatch(jack_kurt(x2), error = function(e) {return(NA)})

  # Point estimates and sampling error (variances)
                         ku[i] = tryCatch((calc.kurtosis(x1) - calc.kurtosis(x2)), error = function(e) {return(NA)}) # this is EXCESS kurtosis DIFFERENCE
                      ku_sv[i] = tryCatch((calc.kurtosis(x1, output = "var") + calc.kurtosis(x2, output = "var")), error = function(e) {return(NA)}) 
               jack_kurt_bc[i] = (x1_kurt_jack$est_bc - x2_kurt_jack$est_bc) # THIS IS EXCESS KURTOSIS DIFFERENCE
               jack_kurt_sv[i] = (x1_kurt_jack$var + x2_kurt_jack$var) 

  # Coverage indicators
                      coverage_ku[i] = coverage(ku[i], ku_sv[i], params$n, true_ku_diff)
              coverage_ku_jack_bc[i] = coverage(jack_kurt_bc[i], jack_kurt_sv[i], params$n, true_ku_diff)
              coverage_ku_jack_sv[i] = coverage(ku[i], jack_kurt_sv[i], params$n, true_ku_diff) # Recommended
              coverage_jack_ku_sv[i] = coverage(jack_kurt_bc[i], ku_sv[i], params$n, true_ku_diff)

  }

  ##-------------------------------------------------##
  # Return data with all the simulation results
  # Calculate bias, relative bias and MCSE for all methods
  ##-------------------------------------------------##

  return(data.frame(
				    	 ku_est = mean(ku),
		       jack_kurt_bc_est = mean(jack_kurt_bc),
		              ku_sv_est = mean(ku_sv),
			   jack_kurt_sv_est = mean(jack_kurt_sv),

    # Bias in point estimates
                 bias_ku = mean(ku)           - true_ku_diff,          
         bias_ku_jack_bc = mean(jack_kurt_bc) - true_ku_diff,           
    
    # Store sampling variances

    # Relative bias for sampling variance, add in adjusted sampling variances
                 bias_ku_jack_sv = ((mean(jack_kurt_sv) - sd(jack_kurt_bc)^2) / sd(jack_kurt_bc)^2)*100, 
                      bias_ku_sv = ((mean(ku_sv) - sd(ku)^2) / sd(ku)^2)*100, 
              bias_ku_ku_jack_sv = ((mean(ku_sv) - sd(jack_kurt_bc)^2) / sd(jack_kurt_bc)^2)*100,
              bias_ku_jack_ku_sv = ((mean(jack_kurt_sv) - sd(ku)^2) / sd(ku)^2)*100,

     # Coverage
                coverage_ku_all = sum(coverage_ku) / nsims,
        coverage_ku_jack_bc_all = sum(coverage_ku_jack_bc) / nsims,
        coverage_ku_jack_sv_all = sum(coverage_ku_jack_sv) / nsims,
        coverage_jack_ku_sv_all = sum(coverage_jack_ku_sv) / nsims,  
    
    # Monte Carlo error
    mcse_bias_jack_ku_bc = sqrt(var(jack_kurt_bc) / length(jack_kurt_bc)), 
    mcse_bias_sv_jack_ku = sqrt(var(jack_kurt_sv) / length(jack_kurt_sv)), 
            mcse_bias_ku = sqrt(var(ku) / length(ku)),      
         mcse_bias_sv_ku = sqrt(var(ku_sv) / length(ku_sv)), 
         mcse_coverage_ku = sqrt((mean(coverage_ku) * (1 - mean(coverage_ku))) / nsims),
     mcse_coverage_ku_jack_bc = sqrt((mean(coverage_ku_jack_bc) * (1 - mean(coverage_ku_jack_bc))) / nsims),
        mcse_coverage_ku_jack_sv = sqrt((mean(coverage_ku_jack_sv) * (1 - mean(coverage_ku_jack_sv))) / nsims),
        mcse_coverage_jack_ku_sv = sqrt((mean(coverage_jack_ku_sv) * (1 - mean(coverage_jack_ku_sv))) / nsims),
         n_sims = length(ku_sv)))
    
}

###------------------------------------------------------------------------###
# # Kurtosis simulation
###------------------------------------------------------------------------###

# Set seed for reproducibility 
  set.seed(860)

	# Parameters
                      nsims = 2500  # simulations, stick with 5000 min but maybe increase to 10,000
								          n = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 150, 500)  
								    mean_g1 = c(0, 0)
								    mean_g2 = c(0, 5)
								variance_g1 = c(1, 1)
								variance_g2 = c(1, 2)
								skewness_g1 = c(0)
								skewness_g2 = c(0)
								kurtosis_g1 = c(2.5, 3, 4, 5, 6)
								kurtosis_g2 = c(2.5, 3, 4, 5, 6)
								
# Create all combinations of scenarios
		scenarios <- expand.grid(mean_g1 = mean_g1, 
                             mean_g2 = mean_g2, 
                             variance_g1 = variance_g1, 
                             variance_g2 = variance_g2, 
                             skewness_g1 = skewness_g1, 
                             skewness_g2 = skewness_g2, 
                             kurtosis_g1 = kurtosis_g1, 
                             kurtosis_g2 = kurtosis_g2) # Create all combinations of scenarios
   scenarios <- scenarios[!duplicated(scenarios), ] # Remove duplicate rows if any
   scenarios$scenario <- 1:nrow(scenarios) # Add a scenario number column

# Create combinations of parameters expanded by sample size vector. Each row is a scenario with a sample size which is used to set up the simulation
	params_all_kur <- data.frame(tidyr::crossing(scenarios, n = n))

# Initialize an empty data frame to store results
     result_kurt <- data.frame()

# Loop through each scenario and run the simulation function
system.time(
for(i in 1:nrow(params_all_kur)) {
       params <- params_all_kur[i,]
  result_kurt <- rbind(result_kurt, sim_kurt(params = params, nsims = nsims))

  print(paste("Simulation for scenario", i, "completed.", "Bias_ku:", round(result_kurt$bias_ku[i], 2), "mcse_bias_sv_ku:", round(result_kurt$mcse_bias_sv_ku[i], 2), "mcse_bias_ku:", round(result_kurt$mcse_bias_ku[i], 2)))
}
)
# 2.6 hours to run

# Merge the results with the scenario parameters
result_kurt <- cbind(params_all_kur, result_kurt)

# Save the results to a file
saveRDS(result_kurt, file = "./output/result_kurt.rds") 

###------------------------------------------------------------------------###
# Visualise the scenarios
###------------------------------------------------------------------------###

# Split into roughly equal chunks of 10 rows
#n_per_group <- 10
#   n_chunks <- ceiling(nrow(scenarios) / n_per_group)

# Split the scenarios into chunks
# df_chunks <- split(scenarios, cut(seq_len(nrow(scenarios)), breaks = n_chunks, labels = FALSE))

# Plot each chunk of scenarios. This will create a grid of plots for each scenario
# mapply(function(x,y) multi_plot_scenarios(x, folder = "output/figs/kurtosis/", name = y), x = df_chunks, y = name) 
