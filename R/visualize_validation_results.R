### ------------------------------------------------------------------------###
# Visualize Meta-Analysis Validation Results
#
# Description: This script loads the results from the validation simulations
# run by `ma_validation_ipd_context.R` and creates visualizations for
# power, bias, and heterogeneity assessment.
### ------------------------------------------------------------------------###

# Load required packages ----
pacman::p_load(
  tidyverse, # Data manipulation and visualization
  here,      # Path management
  patchwork, # Combining plots
  ggthemes,  # Extra themes for ggplot
  cowplot    # For theme_minimal_hgrid()
)

# --- Helper Functions for Plotting ---

#' @title Plot Power Analysis Results
#' @description Visualizes power as a function of n_studies, sample_size, and effect_size.
#' @param power_results_df A data frame with power analysis results.
#' @param show_uncertainty Logical, whether to show uncertainty bands (default TRUE)
#' @return A ggplot object.
plot_power_analysis <- function(power_results_df, show_uncertainty = TRUE) {
  
  if (!"n_studies" %in% names(power_results_df)) {
    stop("Error: `n_studies` column not found in the data frame.")
  }
  
  # Create effect_size_cat from true_effect if it doesn't exist
  if (!"effect_size_cat" %in% names(power_results_df)) {
    power_results_df <- power_results_df %>%
      mutate(
        effect_size_cat = case_when(
          abs(true_effect - 0.3) < 0.05 ~ "small",
          abs(true_effect - 0.5) < 0.05 ~ "medium",
          abs(true_effect - 1.0) < 0.05 ~ "large",
          TRUE ~ as.character(true_effect)
        )
      )
  }
  
  # Calculate approximate SE for power (binomial proportion)
  if (show_uncertainty && "nsims" %in% names(power_results_df)) {
    power_results_df <- power_results_df %>%
      mutate(
        power_se = sqrt(power * (1 - power) / nsims),
        power_lower = pmax(0, power - 1.96 * power_se),
        power_upper = pmin(1, power + 1.96 * power_se)
      )
  }
  
  power_plot <- power_results_df %>%
    mutate(effect_size_cat = factor(effect_size_cat, 
                                    levels = c("small", "medium", "large"),
                                    labels = c("Small Effect", "Medium Effect", "Large Effect"))) %>%
    ggplot(aes(x = n_studies, y = power, color = factor(sample_size), 
               fill = factor(sample_size), group = factor(sample_size))) +
    {if (show_uncertainty && "power_lower" %in% names(power_results_df)) 
      geom_ribbon(aes(ymin = power_lower, ymax = power_upper), alpha = 0.2, color = NA)} +
    geom_line(linewidth = 1.1) +
    geom_point(size = 3) +
    geom_hline(yintercept = 0.8, linetype = "dashed", color = "darkred") +
    facet_grid(metric ~ effect_size_cat) +
    scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
    scale_color_viridis_d(name = "Sample Size\n(per group)") +
    scale_fill_viridis_d(name = "Sample Size\n(per group)", guide = "none") +
    labs(
      title = "Statistical Power by Number of Studies and Sample Size",
      subtitle = "Power to detect a true effect, based on subsampling from IMPC data",
      x = "Number of Studies in Meta-Analysis",
      y = "Power (Rejection Rate)"
    ) +
    theme_minimal_hgrid() +
    theme(
      strip.text = element_text(face = "bold", size = 12),
      legend.position = "bottom",
      plot.title = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 12)
    )
  
  return(power_plot)
}

#' @title Plot Bias Analysis Results
#' @description Visualizes bias as a function of n_studies, sample_size, and effect_size.
#' @param power_results_df A data frame with power and bias analysis results.
#' @param show_uncertainty Logical, whether to show uncertainty bands (default TRUE)
#' @return A ggplot object.
plot_bias_analysis <- function(power_results_df, show_uncertainty = TRUE) {
  
  if (!"n_studies" %in% names(power_results_df)) {
    stop("Error: `n_studies` column not found in the data frame.")
  }
  
  # Create effect_size_cat from true_effect if it doesn't exist
  if (!"effect_size_cat" %in% names(power_results_df)) {
    power_results_df <- power_results_df %>%
      mutate(
        effect_size_cat = case_when(
          abs(true_effect - 0.3) < 0.05 ~ "small",
          abs(true_effect - 0.5) < 0.05 ~ "medium",
          abs(true_effect - 1.0) < 0.05 ~ "large",
          TRUE ~ as.character(true_effect)
        )
      )
  }
  
  # Use RMSE to create uncertainty bands around bias
  if (show_uncertainty && "rmse" %in% names(power_results_df)) {
    power_results_df <- power_results_df %>%
      mutate(
        bias_lower = bias - rmse,
        bias_upper = bias + rmse
      )
  }
  
  bias_plot <- power_results_df %>%
    mutate(effect_size_cat = factor(effect_size_cat, 
                                    levels = c("small", "medium", "large"),
                                    labels = c("Small Effect", "Medium Effect", "Large Effect"))) %>%
    ggplot(aes(x = n_studies, y = bias, color = factor(sample_size), 
               fill = factor(sample_size), group = factor(sample_size))) +
    {if (show_uncertainty && "bias_lower" %in% names(power_results_df)) 
      geom_ribbon(aes(ymin = bias_lower, ymax = bias_upper), alpha = 0.2, color = NA)} +
    geom_line(linewidth = 1.1) +
    geom_point(size = 3) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    facet_grid(metric ~ effect_size_cat) +
    scale_y_continuous(labels = scales::percent) +
    scale_color_viridis_d(name = "Sample Size\n(per group)") +
    scale_fill_viridis_d(name = "Sample Size\n(per group)", guide = "none") +
    labs(
      title = "Estimation Bias by Number of Studies and Sample Size",
      subtitle = "Relative bias of the meta-analytic estimate compared to the true effect",
      x = "Number of Studies in Meta-Analysis",
      y = "Relative Bias"
    ) +
    theme_minimal_hgrid() +
    theme(
      strip.text = element_text(face = "bold", size = 12),
      legend.position = "bottom",
      plot.title = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 12)
    )
  
  return(bias_plot)
}

#' @title Plot Heterogeneity Analysis Results
#' @description Visualizes MA performance under different levels of heterogeneity (tau2)
#' @param het_results_df A data frame with heterogeneity test results
#' @return A ggplot object
plot_heterogeneity_analysis <- function(het_results_df) {
  
  if (!"true_tau2" %in% names(het_results_df)) {
    stop("Error: `true_tau2` column not found in the data frame.")
  }
  
  # Create combined plot showing bias and tau2 estimation
  p1 <- het_results_df %>%
    ggplot(aes(x = true_tau2, y = bias, color = metric, group = metric)) +
    geom_line(linewidth = 1.1) +
    geom_point(size = 3) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    facet_wrap(~ trait, scales = "free_y") +
    scale_color_viridis_d(name = "Metric") +
    labs(
      title = "Estimation Bias Under Heterogeneity",
      subtitle = "How bias changes with between-study variance (tau²)",
      x = "True Between-Study Variance (tau²)",
      y = "Bias in Effect Estimate"
    ) +
    theme_minimal_hgrid() +
    theme(
      strip.text = element_text(face = "bold", size = 11),
      legend.position = "bottom",
      plot.title = element_text(size = 14, face = "bold")
    )
  
  p2 <- het_results_df %>%
    ggplot(aes(x = true_tau2, y = mean_tau2, color = metric, group = metric)) +
    geom_line(linewidth = 1.1) +
    geom_point(size = 3) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "darkred") +
    facet_wrap(~ trait, scales = "free_y") +
    scale_y_log10() +
    scale_color_viridis_d(name = "Metric") +
    labs(
      title = "Heterogeneity Estimation Accuracy",
      subtitle = "Estimated vs. true tau² (dashed line = perfect estimation, log scale)",
      x = "True Between-Study Variance (tau²)",
      y = "Mean Estimated tau² (log scale)"
    ) +
    theme_minimal_hgrid() +
    theme(
      strip.text = element_text(face = "bold", size = 11),
      legend.position = "bottom",
      plot.title = element_text(size = 14, face = "bold")
    )
  
  combined_plot <- p1 / p2
  
  return(combined_plot)
}

#' @title Plot IPD vs Summary MA Comparison
#' @description Visualizes agreement between summary-level and IPD meta-analysis
#' @param ipd_comparison_df A data frame with IPD comparison results
#' @return A ggplot object
plot_ipd_comparison <- function(ipd_comparison_df) {
  
  if (!"summary_estimate" %in% names(ipd_comparison_df) || 
      !"ipd_estimate" %in% names(ipd_comparison_df)) {
    stop("Error: Required columns not found in the data frame.")
  }
  
  # Calculate overall limits for consistent scaling
  all_estimates <- c(ipd_comparison_df$summary_estimate, ipd_comparison_df$ipd_estimate)
  lim_range <- range(all_estimates, na.rm = TRUE)
  lim_buffer <- diff(lim_range) * 0.1
  limits <- c(lim_range[1] - lim_buffer, lim_range[2] + lim_buffer)
  
  # Scatterplot comparing estimates
  p1 <- ipd_comparison_df %>%
    ggplot(aes(x = ipd_estimate, y = summary_estimate, color = metric, shape = trait)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "darkred", linewidth = 1) +
    geom_point(size = 4, alpha = 0.7) +
    geom_errorbar(aes(ymin = summary_ci_lb, ymax = summary_ci_ub), alpha = 0.3, width = 0) +
    geom_errorbarh(aes(xmin = ipd_ci_lb, xmax = ipd_ci_ub), alpha = 0.3, height = 0) +
    coord_fixed(xlim = limits, ylim = limits) +
    scale_color_viridis_d(name = "Metric") +
    scale_shape_manual(name = "Trait", values = c(16, 17, 15, 18)) +
    labs(
      title = "Summary-Level vs. IPD Meta-Analysis Agreement",
      subtitle = "Dashed line indicates perfect agreement",
      x = "IPD Meta-Analysis Estimate",
      y = "Summary-Level Meta-Analysis Estimate"
    ) +
    theme_minimal_grid() +
    theme(
      legend.position = "right",
      plot.title = element_text(size = 14, face = "bold"),
      aspect.ratio = 1
    )
  
  # Bland-Altman plot showing differences
  p2 <- ipd_comparison_df %>%
    mutate(
      mean_estimate = (summary_estimate + ipd_estimate) / 2,
      diff_estimate = summary_estimate - ipd_estimate
    ) %>%
    ggplot(aes(x = mean_estimate, y = diff_estimate, color = metric, shape = trait)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "darkred", linewidth = 1) +
    geom_point(size = 4, alpha = 0.7) +
    scale_color_viridis_d(name = "Metric") +
    scale_shape_manual(name = "Trait", values = c(16, 17, 15, 18)) +
    labs(
      title = "Bland-Altman Plot: Agreement Between Methods",
      subtitle = "Difference vs. average of both estimates",
      x = "Mean of Both Estimates",
      y = "Difference (Summary - IPD)"
    ) +
    theme_minimal_hgrid() +
    theme(
      legend.position = "right",
      plot.title = element_text(size = 14, face = "bold")
    )
  
  combined_plot <- p1 / p2
  
  return(combined_plot)
}


### ------------------------------------------------------------------------###
# Main Execution Block
#
# This section shows how to run new simulations varying `n_studies`
# and then generate the plots.
### ------------------------------------------------------------------------###

# --- STEP 1: Run Simulations (if results don't exist) ---

# To visualize the effect of `n_studies`, we need to run simulations across a range of values.
# The following code block demonstrates how to do this.
# NOTE: This can be time-consuming.

run_power_sims_for_n_studies <- function() {
  
  # Source the main validation script to access its functions
  source(here("R", "ma_validation_ipd_context.R"))
  
  # --- Simulation Parameters ---
  nsims_per_run <- 200 # Number of simulations for each setting (200 is good for stability)
  n_studies_vec <- c(5, 10, 20, 40, 60, 80) # Range of 'number of studies' to test
  sample_sizes <- c(50, 100) # Sample sizes per study group
  effect_sizes <- c("small", "medium", "large")
  metrics <- c("skewness", "kurtosis")
  
  # We'll use one trait to speed things up
  trait_to_test <- "glucose"
  
  # --- Run Simulation Loop ---
  message("Starting power and bias simulations across different n_studies...")
  
  all_power_results <- list()
  
  for (n_st in n_studies_vec) {
    for (ss in sample_sizes) {
      for (es in effect_sizes) {
        for (m in metrics) {
          
          message(paste0("\nRunning: n_studies=", n_st, ", sample_size=", ss, 
                         ", effect=", es, ", metric=", m, "\n"))
          
          # This function is from `ma_validation_ipd_context.R`
          result <- test_ma_known_effects(
            nsims = nsims_per_run,
            n_studies = n_st,
            sample_size = ss,
            trait = trait_to_test,
            metric = m,
            effect_size = es
          )
          
          if (!is.null(result)) {
            all_power_results[[length(all_power_results) + 1]] <- result
          }
        }
      }
    }
  }
  
  # Combine and save results
  power_results_df <- bind_rows(all_power_results)
  
  # Ensure output directory exists
  output_dir <- here("output", "ma_validation")
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  saveRDS(
    power_results_df,
    file.path(output_dir, "power_bias_by_n_studies.rds")
  )
  
  write_csv(
    power_results_df,
    file.path(output_dir, "power_bias_by_n_studies.csv")
  )
  
  message("\nSimulations complete. Results saved to 'output/ma_validation/'.")
  
  return(power_results_df)
}


# --- STEP 2: Load Data and Generate Plots ---

# To run this part, you can either:
# 1. Use existing results from run_complete_validation() in output/ma_validation/
# 2. Run new simulations with run_power_sims_for_n_studies()

visualize_results <- function(run_new_sims = FALSE, show_uncertainty = TRUE) {
  
  # --- Check for existing data files ---
  output_dir <- here("output", "ma_validation")
  
  # Files created by run_power_sims_for_n_studies()
  new_sims_file <- file.path(output_dir, "power_bias_by_n_studies.rds")
  
  # Files created by run_complete_validation() from ma_validation_ipd_context.R
  existing_power_file <- file.path(output_dir, "power_bias_results.rds")
  existing_type1_file <- file.path(output_dir, "type1_error_results.rds")
  het_file <- file.path(output_dir, "heterogeneity_results.rds")
  ipd_file <- file.path(output_dir, "ipd_comparison.rds")
  
  # --- Decide which data to use for power/bias plots ---
  if (run_new_sims) {
    message("Running new simulations as requested...")
    power_data <- run_power_sims_for_n_studies()
    data_source <- "new simulations"
    
  } else if (file.exists(new_sims_file)) {
    message("Loading results from previous run_power_sims_for_n_studies()...")
    power_data <- readRDS(new_sims_file)
    data_source <- "power_bias_by_n_studies.rds"
    
  } else if (file.exists(existing_power_file)) {
    message("Loading results from run_complete_validation()...")
    power_data <- readRDS(existing_power_file)
    data_source <- "power_bias_results.rds"
    
  } else {
    message("No existing results found. Running new simulations...")
    power_data <- run_power_sims_for_n_studies()
    data_source <- "new simulations"
  }
  
  # --- Check data validity ---
  if (nrow(power_data) == 0) {
    stop("The results data frame is empty. Cannot generate plots.")
  }
  
  if (!"n_studies" %in% names(power_data)) {
    stop("Error: 'n_studies' column not found. This data may not be compatible with these plots.")
  }
  
  # --- Generate Power and Bias Plots ---
  message("Generating power and bias plots from: ", data_source)
  if (show_uncertainty) {
    message("  Including uncertainty bands based on simulation variability")
  }
  
  power_plot <- plot_power_analysis(power_data, show_uncertainty = show_uncertainty)
  bias_plot <- plot_bias_analysis(power_data, show_uncertainty = show_uncertainty)
  
  # --- Load and plot heterogeneity results if available ---
  if (file.exists(het_file)) {
    message("Loading and plotting heterogeneity results...")
    het_data <- readRDS(het_file)
    
    if (nrow(het_data) > 0) {
      het_plot <- plot_heterogeneity_analysis(het_data)
    } else {
      message("  Heterogeneity data is empty, skipping plot.")
      het_plot <- NULL
    }
  } else {
    message("No heterogeneity results file found, skipping heterogeneity plots.")
    het_plot <- NULL
  }
  
  # --- Load and plot IPD comparison if available ---
  if (file.exists(ipd_file)) {
    message("Loading and plotting IPD comparison results...")
    ipd_data <- readRDS(ipd_file)
    
    if (nrow(ipd_data) > 0) {
      ipd_plot <- plot_ipd_comparison(ipd_data)
    } else {
      message("  IPD comparison data is empty, skipping plot.")
      ipd_plot <- NULL
    }
  } else {
    message("No IPD comparison file found, skipping IPD comparison plots.")
    ipd_plot <- NULL
  }
  
  # --- Save Plots ---
  plots_dir <- here("output", "ma_validation", "plots")
  if (!dir.exists(plots_dir)) {
    dir.create(plots_dir, recursive = TRUE)
  }
  
  # Save power and bias plots
  ggsave(
    filename = file.path(plots_dir, "power_vs_n_studies.png"),
    plot = power_plot,
    width = 10,
    height = 8,
    dpi = 300,
    bg = "white"
  )
  
  ggsave(
    filename = file.path(plots_dir, "bias_vs_n_studies.png"),
    plot = bias_plot,
    width = 10,
    height = 8,
    dpi = 300,
    bg = "white"
  )
  
  # Save heterogeneity plot if generated
  if (!is.null(het_plot)) {
    ggsave(
      filename = file.path(plots_dir, "heterogeneity_performance.png"),
      plot = het_plot,
      width = 10,
      height = 10,
      dpi = 300,
      bg = "white"
    )
  }
  
  # Save IPD comparison plot if generated
  if (!is.null(ipd_plot)) {
    ggsave(
      filename = file.path(plots_dir, "ipd_vs_summary_comparison.png"),
      plot = ipd_plot,
      width = 10,
      height = 10,
      dpi = 300,
      bg = "white"
    )
  }
  
  message(paste("\nAll plots saved to:", plots_dir))
  
  # Display plots
  print(power_plot)
  print(bias_plot)
  
  if (!is.null(het_plot)) {
    print(het_plot)
  }
  
  if (!is.null(ipd_plot)) {
    print(ipd_plot)
  }
  
  # Return all data invisibly for further analysis if needed
  invisible(list(
    power_data = power_data,
    het_data = if (exists("het_data")) het_data else NULL,
    ipd_data = if (exists("ipd_data")) ipd_data else NULL
  ))
}

# --- HOW TO USE ---
#
# 1. Use existing results (from previous runs) or run new simulations if none exist:
visualize_results()
visualize_results(show_uncertainty = TRUE)
#
# 2. Force new simulations even if results exist:
#    visualize_results(run_new_sims = TRUE)
#
# 3. Generate plots without uncertainty bands:
#    visualize_results(show_uncertainty = FALSE)
#
# 4. If you only want to run the simulations without plotting:
#    # run_power_sims_for_n_studies()
#
### ------------------------------------------------------------------------###
