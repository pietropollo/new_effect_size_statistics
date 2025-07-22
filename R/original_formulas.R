###------------------------------------------------------------------------###
# Monte Carlo Simulation to Assess Estimator Bias
#
# Description:
# This script runs a Monte Carlo simulation to evaluate the accuracy and bias
# of point and variance estimators for skewness and kurtosis. It focuses on
# near-normal distributions by varying either skewness or kurtosis while
# holding the other constant at the value for a normal distribution.
#
###------------------------------------------------------------------------###

#-- 1. Load Required Packages --#
# Ensure you have these packages installed: install.packages(c("pacman"))
# pacman will handle the installation of other packages.
if (!require("pacman")) install.packages("pacman")
pacman::p_load(PearsonDS, moments, dplyr, tidyr, ggplot2)

#-- 2. Estimator Functions (As Provided) --#

#' Calculate Sample Skewness or its Sampling Variance
#'
#' @param x A numeric vector.
#' @param output A character string, either "est" for the estimate or "var" for the variance.
#' @return The calculated skewness or its sampling variance.
calc.skewness <- function(x, output = "est") {
  n <- length(x)
  if (n < 3) return(NA) # Not defined for small n
  
  if (output == "est") { # Skewness estimate (unbiased for normal population)
    (sqrt(n * (n - 1)) / (n - 2)) *
      (mean((x - mean(x))^3) / (sd(x)^3)) # Using sd() which is sqrt(var)
    
  } else if (output == "var") { # Skewness sampling variance (under normality)
    (6 * n * (n - 1)) /
      ((n - 2) * (n + 1) * (n + 3))
  }
}


#' Calculate Sample Kurtosis or its Sampling Variance
#'
#' @param x A numeric vector.
#' @param output A character string, either "est" for the estimate or "var" for the variance.
#' @return The calculated kurtosis or its sampling variance.
calc.kurtosis <- function(x, output = "est") {
  n <- length(x)
  if (n < 4) return(NA) # Not defined for small n
  
  if (output == "est") { # Kurtosis estimate (unbiased for normal population)
    ((((n + 1) * n * (n - 1)) / ((n - 2) * (n - 3))) *
       (sum((x - mean(x))^4) / (sum((x - mean(x))^2)^2))) -
      (3 * ((n - 1)^2) / ((n - 2) * (n - 3)))
    
  } else if (output == "var") { # Kurtosis sampling variance (under normality)
    (24 * n * ((n - 1)^2)) /
      ((n - 3) * (n - 2) * (n + 3) * (n + 5))
  }
}


#-- 3. Main Simulation Function --#

#' Run a single simulation scenario for moment estimators.
#'
#' @param params A named list or data frame row containing simulation parameters:
#'        `n`, `pop_mean`, `pop_var`, `pop_skew`, `pop_kurt`.
#' @param nsims The number of Monte Carlo replications to run.
#' @return A data frame containing the simulation results, including bias and MCSE.
run_moment_simulation <- function(params, nsims) {
  
  # Initialize vectors to store results from each replication
  skew_estimates <- numeric(nsims)
  kurt_estimates <- numeric(nsims)
  
  # Population moments from parameters
  pop_moments <- c(
    mean = params$pop_mean,
    variance = params$pop_var,
    skewness = params$pop_skew,
    kurtosis = params$pop_kurt
  )
  
  # Run the simulation loop
  for (i in 1:nsims) {
    # Generate a sample from the Pearson distribution with specified moments
    sample_data <- tryCatch({
      rpearson(params$n, moments = pop_moments)
    }, error = function(e) {
      # Return NA if parameters are invalid for Pearson distribution
      NA 
    })
    
    if (any(is.na(sample_data))) {
      skew_estimates[i] <- NA
      kurt_estimates[i] <- NA
    } else {
      # Calculate and store the point estimates
      skew_estimates[i] <- calc.skewness(sample_data, output = "est")
      kurt_estimates[i] <- calc.kurtosis(sample_data, output = "est")
    }
  }
  
  # Remove any failed simulation runs
  valid_skew_estimates <- na.omit(skew_estimates)
  valid_kurt_estimates <- na.omit(kurt_estimates)
  
  # Calculate the empirical (observed) variance of the estimators
  empirical_var_skew <- var(valid_skew_estimates)
  empirical_var_kurt <- var(valid_kurt_estimates)
  
  # Calculate the analytical sampling variance using the formulas
  # This is based on the assumption of normality, so we test its accuracy here.
  analytical_var_skew <- calc.skewness(rnorm(params$n), output = "var")
  analytical_var_kurt <- calc.kurtosis(rnorm(params$n), output = "var")
  
  # Calculate bias
  # Point estimate bias: (Average Estimate - True Value)
  bias_skew_point <- mean(valid_skew_estimates) - params$pop_skew
  bias_kurt_point <- mean(valid_kurt_estimates) - params$pop_kurt
  
  # Sampling variance relative bias: (Analytical Variance - Empirical Variance) / Empirical Variance
  # This shows how well the formula for sampling variance holds up.
  relative_bias_skew_var <- (analytical_var_skew - empirical_var_skew) / empirical_var_skew * 100
  relative_bias_kurt_var <- (analytical_var_kurt - empirical_var_kurt) / empirical_var_kurt * 100
  
  # Calculate Monte Carlo Standard Error for the point estimate bias
  mcse_bias_skew_point <- sd(valid_skew_estimates) / sqrt(length(valid_skew_estimates))
  mcse_bias_kurt_point <- sd(valid_kurt_estimates) / sqrt(length(valid_kurt_estimates))
  
  # Return results as a data frame
  data.frame(
    n = params$n,
    pop_skew = params$pop_skew,
    pop_kurt = params$pop_kurt,
    bias_skew_point = bias_skew_point,
    mcse_bias_skew_point = mcse_bias_skew_point,
    bias_kurt_point = bias_kurt_point,
    mcse_bias_kurt_point = mcse_bias_kurt_point,
    relative_bias_skew_var_pct = relative_bias_skew_var,
    relative_bias_kurt_var_pct = relative_bias_kurt_var,
    nsims_complete = length(valid_skew_estimates)
  )
}


#-- 4. Set Up and Run Simulation Scenarios --#

# Set a seed for reproducibility
set.seed(123)

# --- Define the grids of parameters to test ---

# Grid 1: Vary skewness, keep kurtosis normal (kurtosis = 3)
skewness_scenarios <- expand.grid(
  n        = c(20, 50, 100, 200, 500),
  pop_skew = c(0, 0.8, 1.5),
  pop_kurt = 3
)

# Grid 2: Keep skewness normal (skewness = 0), vary kurtosis
kurtosis_scenarios <- expand.grid(
  n        = c(20, 50, 100, 200, 500),
  pop_skew = 0,
  pop_kurt = c(3, 5, 8)
)

# Combine the two grids and remove duplicate rows (the n=..., skew=0, kurt=3 case)
simulation_grid <- distinct(rbind(skewness_scenarios, kurtosis_scenarios))

# Add constant mean and variance
simulation_grid$pop_mean <- 0
simulation_grid$pop_var  <- 1

# Number of replications for each scenario
num_simulations <- 10000

# Initialize a list to store results from each scenario
all_results <- list()

# Loop through each row (scenario) in the grid
cat("Starting Monte Carlo Simulation...\n")
for (i in 1:nrow(simulation_grid)) {
  params <- simulation_grid[i, ]
  
  # Run the simulation for the current scenario
  all_results[[i]] <- run_moment_simulation(params, nsims = num_simulations)
  
  # Print progress
  cat(sprintf("Finished scenario %d of %d: n=%d, skew=%.1f, kurt=%.1f\n",
              i, nrow(simulation_grid), params$n, params$pop_skew, params$pop_kurt))
}
cat("Simulation complete.\n\n")

# Combine all results into a single data frame
results_df <- bind_rows(all_results)


#-- 5. Display and Visualize Results --#

# Print the full results table
print(results_df, digits = 4)

# --- Visualize Skewness Point Estimator Bias (when Kurtosis is Normal) ---
skew_bias_plot <- results_df %>%
  filter(pop_kurt == 3) %>%
  ggplot(aes(x = factor(n), y = bias_skew_point, fill = factor(pop_skew))) +
  geom_col(position = "dodge") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  labs(
    title = "Bias of Skewness Point Estimator (Population Kurtosis = 3)",
    x = "Sample Size (n)",
    y = "Bias (Average Estimate - True Value)",
    fill = "Population Skewness"
  ) +
  theme_minimal() +
  scale_fill_viridis_d()

print(skew_bias_plot)


# --- Visualize Kurtosis Point Estimator Bias (when Skewness is Normal) ---
kurt_bias_plot <- results_df %>%
  filter(pop_skew == 0) %>%
  ggplot(aes(x = factor(n), y = bias_kurt_point, fill = factor(pop_kurt))) +
  geom_col(position = "dodge") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  labs(
    title = "Bias of Kurtosis Point Estimator (Population Skewness = 0)",
    x = "Sample Size (n)",
    y = "Bias (Average Estimate - True Value)",
    fill = "Population Kurtosis"
  ) +
  theme_minimal() +
  scale_fill_viridis_d(option = "plasma")

print(kurt_bias_plot)

# --- Visualize Relative Bias of Skewness Variance Estimator (when Kurtosis is Normal) ---
skew_var_bias_plot <- results_df %>%
  filter(pop_kurt == 3) %>%
  ggplot(aes(x = factor(n), y = relative_bias_skew_var_pct, fill = factor(pop_skew))) +
  geom_col(position = "dodge") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  labs(
    title = "Relative Bias of Skewness Sampling Variance Estimator",
    subtitle = "Population Kurtosis held at 3 (Normal)",
    x = "Sample Size (n)",
    y = "Relative Bias (%)",
    fill = "Population Skewness"
  ) +
  theme_minimal() +
  scale_fill_viridis_d()

print(skew_var_bias_plot)

# --- Visualize Relative Bias of Kurtosis Variance Estimator (when Skewness is Normal) ---
kurt_var_bias_plot <- results_df %>%
  filter(pop_skew == 0) %>%
  ggplot(aes(x = factor(n), y = relative_bias_kurt_var_pct, fill = factor(pop_kurt))) +
  geom_col(position = "dodge") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  labs(
    title = "Relative Bias of Kurtosis Sampling Variance Estimator",
    subtitle = "Population Skewness held at 0 (Normal)",
    x = "Sample Size (n)",
    y = "Relative Bias (%)",
    fill = "Population Kurtosis"
  ) +
  theme_minimal() +
  scale_fill_viridis_d(option = "plasma")

print(kurt_var_bias_plot)

