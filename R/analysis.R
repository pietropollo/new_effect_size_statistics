###------------------------------------------------------------------------###
# Simulation Results for New Effect Size Statistics
# Description: Plotting the results of the simulations for new effect size statistics.
###------------------------------------------------------------------------###

# Clean working environment
rm(list = ls()) # Remove all objects from the environment

# Load required packages
	pacman::p_load(moments, PearsonDS, tidyverse, patchwork, latex2exp)  
  	source("./R/func.R") # Load the functions from func.R

	# Load the simulation results
	result_kurt <- readRDS("./output/result_kurt.rds") # Load the kurtosis results
	 result_cor <- readRDS("./output/result_cor.rds") # Load the correlation results	
	result_skew <- readRDS("./output/result_skewness.rds") # Load the skewness results

# Rename second n column to n_sim 
     colnames(result_cor)[17] <- "n_sim"

# Plot function for the results
plot_bias_violin <- function(data, y_var, y_lab) {
  ggplot(data, aes(x = factor(n), y = .data[[y_var]], fill = factor(n))) + 
    geom_violin() + 
    geom_hline(aes(yintercept = 0), linetype = "dashed", color = "black") +
    labs(x = "Sample Size", y = TeX(y_lab)) +
    scale_fill_viridis_d() +
    theme_classic() + 
    theme(
      legend.position = "none",
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12),
      plot.tag = element_text(size = 16, face = "bold")
    )
}

##------------------------------------------------------------------------##
## Bias in Estimates
##------------------------------------------------------------------------##
## Skewness

# Skewness plots
bias_sk_plot        <- plot_bias_violin(result_skew, "bias_sk",         "Bias $\\Delta sk$")
bias_sk_plot_delta  <- plot_bias_violin(result_skew, "bias_sk_delta",   "Bias $\\Delta sk$ (Delta Method)")
bias_sk_plot_boot   <- plot_bias_violin(result_skew, "bias_sk_boot_bc", "Bias $\\Delta sk$ (Bootstrap)")
bias_sk_plot_jack   <- plot_bias_violin(result_skew, "bias_sk_jack_bc", "Bias $\\Delta sk$ (Jackknife)")

# Kurtosis
# Kurtosis plots
bias_ku_plot        <- plot_bias_violin(result_kurt, "bias_ku",         "Bias $\\Delta ku$")
bias_ku_plot_delta  <- plot_bias_violin(result_kurt, "bias_ku_delta",   "Bias $\\Delta ku$ (Delta Method)")
bias_ku_plot_boot   <- plot_bias_violin(result_kurt, "bias_ku_boot_bc", "Bias $\\Delta ku$ (Bootstrap)")
bias_ku_plot_jack   <- plot_bias_violin(result_kurt, "bias_ku_jack_bc", "Bias $\\Delta ku$ (Jackknife)")

# Correlation plots
bias_z_plot         <- plot_bias_violin(result_cor, "bias_d_cor",        "Bias $\\Delta Z_{r}$")
bias_z_plot_boot    <- plot_bias_violin(result_cor, "bias_boot_d_cor",   "Bias $\\Delta Z_{r}$ (Bootstrapped)")
bias_z_plot_jack    <- plot_bias_violin(result_cor, "bias_jack_d_cor",   "Bias $\\Delta Z_{r}$ (Jackknife)")

# Combine all bias plots
est_plot <- (bias_sk_plot| bias_sk_plot_delta | bias_sk_plot_boot | bias_sk_plot_jack) / (bias_ku_plot | bias_ku_plot_delta | bias_ku_plot_boot | bias_ku_plot_jack)  +
plot_annotation(tag_levels = 'A', tag_suffix = ")") &
theme(plot.tag = element_text(size = 16, face = "bold"))



##------------------------------------------------------------------------##
## Relative Bias in Sampling Error of Estimates
##------------------------------------------------------------------------##

# Skewness
bias_sv_sk       <- plot_bias_violin(result_skew, "bias_sk_sv",         "Relative Bias $SV_{\\Delta sk}$ (%)")
bias_sv_sk_delta <- plot_bias_violin(result_skew, "bias_sk_delta_sv",   "Relative Bias $SV_{\\Delta sk}$ (%) (Delta method)")
bias_sv_sk_boot  <- plot_bias_violin(result_skew, "bias_sk_boot_sv",    "Relative Bias $SV_{\\Delta sk}$ (%) (Bootstrap)")
bias_sv_sk_jack  <- plot_bias_violin(result_skew, "bias_sk_jack_sv",    "Relative Bias $SV_{\\Delta sk}$ (%) (Jackknife)")

# Kurtosis
bias_sv_ku       <- plot_bias_violin(result_kurt, "bias_ku_sv",         "Relative Bias $SV_{\\Delta ku}$ (%)")
bias_sv_ku_delta <- plot_bias_violin(result_kurt, "bias_ku_delta_sv",   "Relative Bias $SV_{\\Delta ku}$ (%) (Delta method)")
bias_sv_ku_boot  <- plot_bias_violin(result_kurt, "bias_ku_boot_sv",    "Relative Bias $SV_{\\Delta ku}$ (%) (Bootstrap)")
bias_sv_ku_jack  <- plot_bias_violin(result_kurt, "bias_ku_jack_sv",    "Relative Bias $SV_{\\Delta ku}$ (%) (Jackknife)")

# Correlation
bias_sv_z        <- plot_bias_violin(result_cor,  "bias_d_cor_sv",      "Relative Bias $SV_{\\Delta Z_{r}}$ (%)")
bias_sv_z_boot   <- plot_bias_violin(result_cor,  "bias_boot_d_cor_sv", "Relative Bias $SV_{\\Delta Z_{r}}$ (%) (Bootstrapped)")
bias_sv_z_jack   <- plot_bias_violin(result_cor,  "bias_jack_d_cor_sv", "Relative Bias $SV_{\\Delta Z_{r}}$ (%) (Jackknife)")

# Combine all plots
final_rel_bias_plot <- (bias_sv_sk | bias_sv_sk_delta | bias_sv_sk_boot | bias_sv_sk_jack) / (bias_sv_ku | bias_sv_ku_delta | bias_sv_ku_boot | bias_sv_ku_jack) +
  plot_annotation(tag_levels = 'A', tag_suffix = ")") &
  theme(plot.tag = element_text(size = 16, face = "bold"))
  
# Show the final plot
print(final_rel_bias_plot)

# Make the correlation plots
cor_plot <- (bias_z_plot | bias_z_plot_boot   | bias_z_plot_jack) / (bias_sv_z| bias_sv_z_boot   | bias_sv_z_jack)+
  plot_annotation(tag_levels = 'A', tag_suffix = ")") &
  theme(plot.tag = element_text(size = 16, face = "bold"))
