###------------------------------------------------------------------------###
# Simulate Pearson distribution with specified moments
# Description: Equal sample size between groups. 
###------------------------------------------------------------------------###

# Clear the environment
rm(list = ls())

# Load required packages
	pacman::p_load(moments, PearsonDS, tidyverse, patchwork)
  source("./R/func.R") # Load the functions from func.R
  
# functions for calculating effect sizes ----
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

# Testing params object # out when done. This would give normal distribution with slight skewness 
#params <- data.frame(mean_g1 = 0, mean_g2 = 0, variance_g1 = 1, variance_g2 = 1, skewness_g1 = 1, skewness_g2 = 1.5, kurtosis_g1 = 6, kurtosis_g2 = 6, n = 100)

#| @title Simulate Pearson distribution with specified moments
#| @description Function to simulate Pearson distribution with defined moments: mean, variance, skewness and kurtosis
#| and calculate bias and coverage
#| @param params Data frame containing parameters for the simulation. Note that this is from a params data frame created by the expand.grid function above and so included data from a single row which is a single scenario. Params must define the sample size and the moments for both groups (mean_g1, mean_g2, variance_g1, variance_g2, skewness_g1, skewness_g2, kurtosis_g1, kurtosis_g2).
#| @param nsims Number of simulations to run (default is 1000)
#| @return A data frame with bias with respect to point estimates of effect size and its sampling variance for each scenario.
#| @examples
#| sim_func(n = 50, params = params[i,], nsims = 1000)

sim_func <- function(params, nsims = nsims, type = c("skewness", "kurtosis")) {
       type <- match.arg(type)
  
  # Vectors to store results. Note that we have a bunch of new functions to calculate the effect sizes and sampling variances which don't assume normality
         sk <- numeric(nsims)
         ku <- numeric(nsims)
      sk_sv <- numeric(nsims)
      ku_sv <- numeric(nsims)
 
 # Point estimates for skewness and kurtosis

  jack_skew_bc <- numeric(nsims) # Jackknife bias-corrected method for skewness
  jack_kurt_bc <- numeric(nsims) # Jackknife bias-corrected method for kurtosis

# Sampling variance for skewness and kurtosis

  jack_skew_sv <- numeric(nsims) # Jackknife method sampling variance for skewness
  jack_kurt_sv <- numeric(nsims) # Jackknife method sampling variance for kurtosis

# Coverage indicators
  coverage_sk <- numeric(nsims)         # Coverage for skewness with sampling variance for skewness
  coverage_sk_jack_bc <- numeric(nsims) # Coverage for skewness with jackknife bias-corrected method and jackknife sampling variance for skewness
  coverage_sk_jack_sv <- numeric(nsims) # Coverage for skewness with jackknife sampling variance for skewness
  coverage_ku <- numeric(nsims)         # Coverage for kurtosis with sampling variance for kurtosis
  coverage_ku_jack_bc <- numeric(nsims) # Coverage for kurtosis with jackknife bias-corrected method and jackknife sampling variance for kurtosis
  coverage_ku_jack_sv <- numeric(nsims) # Coverage for kurtosis with jackknife sampling variance for kurtosis

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

  if(type == "skewness") {

    # Do the computations and save to objects
     x1_skew_jack <- tryCatch(jack_skew(x1),  error = function(e) {return(NA)})
     x2_skew_jack <- tryCatch(jack_skew(x2),  error = function(e) {return(NA)})

  # Point estimates
       sk[i] = tryCatch((calc.skewness(x1) - calc.skewness(x2)), error = function(e) {return(NA)}) # small sample size corrected
    sk_sv[i] = tryCatch((calc.skewness(x1, output = "var") + calc.skewness(x2, output = "var")), error = function(e) {return(NA)})
     jack_skew_bc[i] = (x1_skew_jack$est_bc - x2_skew_jack$est_bc)
     jack_skew_sv[i] =  (x1_skew_jack$var + x2_skew_jack$var)

  # Coverage indicators
    coverage_sk[i] = ((params$skewness_g1 - params$skewness_g2) >= (sk[i] - 1.96 * sqrt(sk_sv[i])) && (params$skewness_g1 - params$skewness_g2) <= (sk[i] + 1.96 * sqrt(sk_sv[i])))
    coverage_sk_jack_bc[i] = ((params$skewness_g1 - params$skewness_g2) >= (jack_skew_bc[i] - 1.96 * sqrt(jack_skew_sv[i])) && (params$skewness_g1 - params$skewness_g2) <= (jack_skew_bc[i] + 1.96 * sqrt(jack_skew_sv[i])))
    coverage_sk_jack_sv[i] = ((params$skewness_g1 - params$skewness_g2) >= (sk[i] - 1.96 * sqrt(jack_skew_sv[i])) && (params$skewness_g1 - params$skewness_g2) <= (sk[i] + 1.96 * sqrt(jack_skew_sv[i])))

  }

  if(type == "kurtosis") {

      # Do the computations and save to objects

    x1_kurt_jack <- tryCatch(jack_kurt(x1), error = function(e) {return(NA)})
    x2_kurt_jack <- tryCatch(jack_kurt(x2), error = function(e) {return(NA)})

  # Point estimates
       ku[i] = tryCatch((calc.kurtosis(x1) - calc.kurtosis(x2)), error = function(e) {return(NA)})
    ku_sv[i] = tryCatch((calc.kurtosis(x1, output = "var") + calc.kurtosis(x2, output = "var")), error = function(e) {return(NA)}) 
         jack_kurt_bc[i] = (x1_kurt_jack$est_bc - x2_kurt_jack$est_bc)
         jack_kurt_sv[i] =  (x1_kurt_jack$var + x2_kurt_jack$var)

  # Coverage indicators
    coverage_ku[i] = ((params$kurtosis_g1 - params$kurtosis_g2) >= (ku[i] - 1.96 * sqrt(ku_sv[i])) && (params$kurtosis_g1 - params$kurtosis_g2) <= (ku[i] + 1.96 * sqrt(ku_sv[i])))
    coverage_ku_jack_bc[i] = ((params$kurtosis_g1 - params$kurtosis_g2) >= (jack_kurt_bc[i] - 1.96 * sqrt(jack_kurt_sv[i])) && (params$kurtosis_g1 - params$kurtosis_g2) <= (jack_kurt_bc[i] + 1.96 * sqrt(jack_kurt_sv[i])))
    coverage_ku_jack_sv[i] = ((params$kurtosis_g1 - params$kurtosis_g2) >= (ku[i] - 1.96 * sqrt(jack_kurt_sv[i])) && (params$kurtosis_g1 - params$kurtosis_g2) <= (ku[i] + 1.96 * sqrt(jack_kurt_sv[i])))

  }
}

##-------------------------------------------------##
# Return data with all the simulation results
# Calculate bias, relative bias and MCSE for all methods
##-------------------------------------------------##

if(type == "skewness") {

  return(data.frame(
    # Original analytical formulas for skewness     
               bias_sk = mean(sk)           - (params$skewness_g1 - params$skewness_g2),   
       bias_sk_jack_bc = mean(jack_skew_bc) - (params$skewness_g1 - params$skewness_g2),      
    
    # Sampling variances. Relative bias
            bias_sk_sv = ((mean(sk_sv)        - sd(sk)^2) / sd(sk)^2)*100, 
       bias_sk_jack_sv = ((mean(jack_skew_sv) - sd(jack_skew_bc)^2) / sd(jack_skew_bc)^2)*100, 
    bias_sk_sk_jack_sv = ((mean(sk_sv)        - sd(jack_skew_bc)^2) / sd(jack_skew_bc)^2)*100, 
    bias_sk_jack_sk_sv = ((mean(jack_skew_sv) - sd(sk)^2) / sd(sk)^2)*100, 

    # Coverage
            coverage_sk = sum(coverage_sk) / nsims,
        coverage_sk_jack_bc = sum(coverage_sk_jack_bc) / nsims,
        coverage_sk_jack_sv = sum(coverage_sk_jack_sv) / nsims,
    
    # Monte Carlo Error
    mcse_bias_sv_jack_sk = sqrt(var(jack_skew_sv) / length(jack_skew_sv)), 
            mcse_bias_sk = sqrt(var(sk) / length(sk)), 
         mcse_bias_sv_sk = sqrt(var(sk_sv) / length(sk_sv)), 
       mcse_bias_jack_sk = sqrt(var(jack_skew_bc) / length(jack_skew_bc)),
       n_sims = length(sk_sv))) # Number of simulations
    }

if(type == "kurtosis") {

  return(data.frame(     
    # Bias in point estimates
                 bias_ku = mean(ku)           - (params$kurtosis_g1 - params$kurtosis_g2),          
         bias_ku_jack_bc = mean(jack_kurt_bc) - (params$kurtosis_g1 - params$kurtosis_g2),           
         
    # Relative bias for sampling variance
            bias_ku_jack_sv = ((mean(jack_kurt_sv) - sd(jack_kurt_bc)^2) / sd(jack_kurt_bc)^2)*100, 
                 bias_ku_sv = ((mean(ku_sv) - sd(ku)^2) / sd(ku)^2)*100, 
         bias_ku_ku_jack_sv = ((mean(ku_sv) - sd(jack_kurt_bc)^2) / sd(jack_kurt_bc)^2)*100,
         bias_ku_jack_ku_sv = ((mean(jack_kurt_sv) - sd(ku)^2) / sd(ku)^2)*100,

     # Coverage
            coverage_ku = sum(coverage_ku) / nsims,
        coverage_ku_jack_bc = sum(coverage_ku_jack_bc) / nsims,

    # Monte Carlo error
    mcse_bias_jack_ku_bc = sqrt(var(jack_kurt_bc) / length(jack_kurt_bc)), 
    mcse_bias_sv_jack_ku = sqrt(var(jack_kurt_sv) / length(jack_kurt_sv)), 
            mcse_bias_ku = sqrt(var(ku) / length(ku)),      
         mcse_bias_sv_ku = sqrt(var(ku_sv) / length(ku_sv)), 
         n_sims = length(ku_sv)))
    }
}


###------------------------------------------------------------------------###
# Run simulations for all scenarios assuming equal sample size 
###------------------------------------------------------------------------###

###------------------------------------------------------------------------###
# # Skewness simulation 
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
								skewness_g1 = c(-1, -0.5, 0, 0.5, 1)
								skewness_g2 = c(-1, -0.5, 0, 0.5, 1)
								kurtosis_g1 = c(3)
								kurtosis_g2 = c(3)

# Create all combinations of scenarios 	
scenarios <- expand.grid(    mean_g1 = mean_g1, 
                             mean_g2 = mean_g2, 
                         variance_g1 = variance_g1, 
                         variance_g2 = variance_g2, 
                         skewness_g1 = skewness_g1, 
                         skewness_g2 = skewness_g2, 
                         kurtosis_g1 = kurtosis_g1, 
                         kurtosis_g2 = kurtosis_g2) # Create all combinations of scenarios
scenarios <- scenarios[!duplicated(scenarios), ] # Remove duplicate rows if any
scenarios$scenario <- 1:nrow(scenarios) # Add a scenario number column

# Visualise the scenarios
# Split into roughly equal chunks of 10 rows
n_per_group <- 10
   n_chunks <- ceiling(nrow(scenarios) / n_per_group)

# Split the scenarios into chunks
df_chunks <- split(scenarios, cut(seq_len(nrow(scenarios)), breaks = n_chunks, labels = FALSE))
name = 1:length(df_chunks)

# Plot each chunk of scenarios. This will create a grid of plots for each scenario
mapply(function(x,y) multi_plot_scenarios(x, folder = "output/figs/skewness/", name = y), x = df_chunks, y = name) 

# Create combinations of parameters expanded by sample size vector. Each row is a scenario with a sample size which is used to set up the simulation
	params_all <- data.frame(tidyr::crossing(scenarios, n = n))               

# Initialize an empty data frame to store results
result_skewness <- data.frame()

# Loop through each scenario and run the simulation function
system.time(
for(i in 1:nrow(params_all)) {
           params <- params_all[i,]
  result_skewness <- rbind(result_skewness, sim_func(params = params, nsims = nsims, type = "skewness"))

  print(paste("Simulation for scenario", i, "completed. Bias_sk:", round(result_skewness$bias_sk[i], 2), "mcse_bias_sv_sk:", round(result_skewness$mcse_bias_sv_sk[i], 2), "mcse_bias_sk:", round(result_skewness$mcse_bias_sk[i], 2)))
}
)
# # 3 hours to run
# Merge the results with the scenario parameters
result_skewness <- cbind(params_all, result_skewness)
saveRDS(result_skewness, file = "./output/result_skewness.rds") # Save the results to a file

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

# Visualise the scenarios
# Split into roughly equal chunks of 10 rows
n_per_group <- 10
   n_chunks <- ceiling(nrow(scenarios) / n_per_group)

# Split the scenarios into chunks
df_chunks <- split(scenarios, cut(seq_len(nrow(scenarios)), breaks = n_chunks, labels = FALSE))

# Plot each chunk of scenarios. This will create a grid of plots for each scenario
mapply(function(x,y) multi_plot_scenarios(x, folder = "output/figs/kurtosis/", name = y), x = df_chunks, y = name) 

# Create combinations of parameters expanded by sample size vector. Each row is a scenario with a sample size which is used to set up the simulation
	params_all_kur <- data.frame(tidyr::crossing(scenarios, n = n))

# Initialize an empty data frame to store results
     result_kurt <- data.frame()

# Loop through each scenario and run the simulation function
system.time(
for(i in 1:nrow(params_all_kur)) {
       params <- params_all_kur[i,]
  result_kurt <- rbind(result_kurt, sim_func(params = params, nsims = nsims, type = "kurtosis"))

  print(paste("Simulation for scenario", i, "completed.", "Bias_ku:", round(result_kurt$bias_ku[i], 2), "mcse_bias_sv_ku:", round(result_kurt$mcse_bias_sv_ku[i], 2), "mcse_bias_ku:", round(result_kurt$mcse_bias_ku[i], 2)))
}
)
# 2.6 hours to run

# Merge the results with the scenario parameters
result_kurt <- cbind(params_all_kur, result_kurt)

# Save the results to a file
saveRDS(result_kurt, file = "./output/result_kurt.rds") 


# Some visuals of the scenarios
# Testing params object # out when done. This would give normal distribution with slight skewness 
params <- data.frame(mean_g1 = 0, mean_g2 = 0, variance_g1 = 1, variance_g2 = 1, skewness_g1 = 0, skewness_g2 = 0.5, kurtosis_g1 = 3, kurtosis_g2 = 3)
plot_scenario(params, print = TRUE, xpos = 4)
t <- create_dat(params)
