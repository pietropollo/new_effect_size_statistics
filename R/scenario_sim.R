
##---------------------------------------------##
## Simulation for KU to look at distribution of 
## estimates and coverage of CIs from a given scenario
##---------------------------------------------##
# Load required packages
pacman::p_load(moments, 
               PearsonDS, 
               tidyverse, 
               patchwork, 
               latex2exp)  

# Functions
  source("R/func.R") # Load effect size functions

# Set seed for reproducibility 
  set.seed(860)

##-----------------------------##
# Define simulation scenarios
##-----------------------------##
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

# Lets keep sample size at 10 for true effect size calculations. 
  params_true <- data.frame(tidyr::crossing(scenarios, n = 500))  %>% mutate(diff_ku = kurtosis_g1 - kurtosis_g2)  %>% filter(diff_ku >= 0 & mean_g1 == 0 & mean_g2 == 0 & variance_g1 == 1 & variance_g2 == 1 & skewness_g1 == 0 & skewness_g2 == 0 & skewness_g1 == 0 & skewness_g2 == 0 & kurtosis_g1 == 5) # Only keep scenarios where group 1 kurtosis is greater than or equal to group 2 kurtosis

## Sample parameters for a single simulation run
  #params <- params_true[1, ] # Change the index to run different scenarios

for(i in 1:nrow(params_all_kur)) {
  params  <- params_true[i, ]
#-----------------------------##
# Setup objects
#-----------------------------##

# Vectors to store results. Note that we have a bunch of new functions to calculate the effect sizes and sampling variances which don't assume normality
 # Point estimates for kurtosis
        jack_kurt_bc <- numeric(nsims) # Jackknife bias-corrected method for kurtosis
                  ku <- numeric(nsims)
         boot_bc_est <- numeric(nsims) # Bootstrap bias-corrected method for kurtosis
         
   
  # Sampling variance for kurtosis
             ku_sv <- numeric(nsims) # Sampling variance for kurtosis
      jack_kurt_sv <- numeric(nsims) # Jackknife method sampling variance for kurtosis
           boot_sv <- numeric(nsims) # Bootstrap sampling variance for kurtosis
   
  # Coverage indicators
                coverage_ku <- numeric(nsims) # Coverage for kurtosis effect size
			  coverage_ku_jack_bc <- numeric(nsims) # Coverage for kurtosis with jackknife bias-corrected method
			  coverage_ku_jack_sv <- numeric(nsims) # Coverage for kurtosis with jackknife bias-corrected method and sampling variance for kurtosis
			  coverage_ku_jack_adj_sv <- numeric(nsims) # Coverage for kurtosis with jackknife bias-corrected method and adjusted jackknife sampling variance for kurtosis
       coverage_jack_ku_sv <- numeric(nsims) # Coverage for kurtosis with jackknife bias-corrected method and sampling variance for kurtosis
       boot_coverage_ku <- numeric(nsims) # Coverage for kurtosis with bootstrap method

  # Boot BCA 95% CIs
   boot_sv_l95 <- numeric(nsims)
   boot_sv_u95 <- numeric(nsims)

#-----------------------------##
# Run Simulations
#-----------------------------##

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
	# Using simulated data, calculate jackknife and bootstrap estimates
	# Jackknife
	x1_kurt_jack <- tryCatch(jack_kurt(x1), error = function(e) {return(NA)})
    x2_kurt_jack <- tryCatch(jack_kurt(x2), error = function(e) {return(NA)})

	# Bootstrap. Calculate bootstrap kurtosis difference with BCa CIs within function
        data  <- data.frame(x1 = x1, x2 = x2)
    kurt_boot <- tryCatch(boot_kurt_bca(data), error = function(e) {return(NA)})
    
  # Traditional point estimators and sampling variances
      ku[i] = tryCatch((calc.kurtosis(x1) - calc.kurtosis(x2)), error = function(e) {return(NA)}) 
   ku_sv[i] = tryCatch((calc.kurtosis(x1, output = "var") + calc.kurtosis(x2, output = "var")), error = function(e) {return(NA)}) 

  # Jacknife bias-corrected estimators and sampling variances
     jack_kurt_bc[i] = (x1_kurt_jack$est_bc - x2_kurt_jack$est_bc) 
     jack_kurt_sv[i] = (x1_kurt_jack$var    + x2_kurt_jack$var) 

  # Bootstrap bias-corrected estimators and sampling variances
      boot_bc_est[i] = kurt_boot$est_bc
      boot_sv_l95[i] = kurt_boot$est_ci_l
      boot_sv_u95[i] = kurt_boot$est_ci_u

  }


# plot distributions of the point estimates AND sampling variances for jacknife, bootstrap, and the regular estimate with regular sampling variance
	main <- paste0("Scenario ", params$scenario, ": n = ", params$n, 
				   ", Kurtosis G1 = ", params$kurtosis_g1, ", Kurtosis G2 = ", params$kurtosis_g2)


s1_p1 <- ggplot(data.frame(jack_kurt_bc = jack_kurt_bc),
       aes(x = jack_kurt_bc)) +
  geom_histogram(bins = 30, colour = "black", fill = "grey80") +
  labs(
    title = main,
    x = "Jackknife Bias-Corrected Kurtosis Difference",
    y = "Count"
  ) +
  theme_classic()

s1_p2  <- ggplot(data.frame(boot_bc_est = boot_bc_est),
       aes(x = boot_bc_est)) +
  geom_histogram(bins = 30, colour = "black", fill = "grey80") +
  labs(
    title = main,
    x = "Bootstrap Bias-Corrected Kurtosis Difference",
    y = "Count"
  ) +
  theme_classic()

ggsave(paste0("output/figs/scenario_", gsub(" ", "_", main), "_kurtosis_distributions.png"), plot = s1_p1 + s1_p2 + plot_layout(ncol = 2), width = 12.269938, height = 5.341615 )
}