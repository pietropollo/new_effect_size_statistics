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


  # Lets find out what is the best estimator
      result_kurt_least_bias <- result_kurt  %>% summarise(across(starts_with("bias_"), ~ sum(.x^2, na.rm = TRUE), .names = "least_bias_{col}"))  %>% t()  %>%  data.frame() %>% arrange() 
      colnames(result_kurt_least_bias) <- "Least Bias"
      result_skew_least_bias <- result_skew  %>% summarise(across(starts_with("bias_"), ~ sum(.x^2, na.rm = TRUE), .names = "least_bias_{col}"))  %>% t()  %>%  data.frame() %>% arrange() 
      colnames(result_skew_least_bias) <- "Least Bias"
      result_cor_least_bias <- result_cor  %>% summarise(across(starts_with("bias_"), ~ sum(.x^2, na.rm = TRUE), .names = "least_bias_{col}"))  %>% t()  %>%  data.frame() %>% arrange() 
      colnames(result_cor_least_bias) <- "Least Bias"

# Plot function for the results
plot_bias_violin <- function(data, y_var, y_lab, title = "", ylim = c(-1, 1)) {
  ggplot(data, aes(x = factor(n), y = .data[[y_var]], fill = factor(n))) + 
    ylim(ylim) +
    geom_violin() + 
    geom_hline(aes(yintercept = 0), linetype = "dashed", color = "black") +
    labs(x = "Sample Size", y = TeX(y_lab), title = TeX(title)) +
    scale_fill_viridis_d() +
    theme_classic() + 
    theme(
      legend.position = "none",
           axis.title = element_text(size = 12, face = "bold"),
            axis.text = element_text(size = 12),
             plot.tag = element_text(size = 16, face = "bold")
    )
}

##------------------------------------------------------------------------##
## Bias in Estimates
##------------------------------------------------------------------------##

# Skewness plots
bias_sk_plot        <- plot_bias_violin(result_skew, "bias_sk",         "Bias $\\Delta sk$", ylim = c(-0.5, 0.5))
bias_sk_plot_jack   <- plot_bias_violin(result_skew, "bias_sk_jack_bc", "Bias $\\Delta sk$ (Jackknife)", ylim = c(-0.5, 0.5))

# Kurtosis plots
bias_ku_plot        <- plot_bias_violin(result_kurt, "bias_ku",         "Bias $\\Delta ku$", ylim = c(-3, 3))
bias_ku_plot_jack   <- plot_bias_violin(result_kurt, "bias_ku_jack_bc", "Bias $\\Delta ku$ (Jackknife)", ylim = c(-3, 3))

# Correlation plots
bias_z_plot         <- plot_bias_violin(result_cor, "bias_d_cor",        "Bias $\\Delta Z_{r}$", ylim = c(-0.1, 0.1))
bias_z_plot_jack    <- plot_bias_violin(result_cor, "bias_jack_d_cor",   "Bias $\\Delta Z_{r}$ (Jackknife)", ylim = c(-0.1, 0.1))

# Combine all bias plots
est_plot <- (bias_sk_plot| bias_sk_plot_jack) / (bias_ku_plot | bias_ku_plot_jack) / (bias_z_plot | bias_z_plot_jack) +
  plot_annotation(tag_levels = 'A', tag_suffix = ")") &
  theme(plot.tag = element_text(size = 16, face = "bold"))
ggsave("./output/figs/bias.png", plot = est_plot, width = 13, height = 12)

##------------------------------------------------------------------------##
## Relative Bias in Sampling Error of Estimates
##------------------------------------------------------------------------##

# Skewness
bias_sv_sk       <- plot_bias_violin(result_skew, "bias_sk_sv",         "Relative Bias $SV_{\\Delta sk}$ (%)", title = "((mean(sk_sv)- sd(sk)^2) / sd(sk)^2)*100", ylim = c(-80, 80))
bias_sv_sk_jack  <- plot_bias_violin(result_skew, "bias_sk_jack_sv",    "Relative Bias $SV_{\\Delta sk}$ (%)", title = "((mean(jack_skew_sv) - sd(jack_skew_bc)^2) / sd(jack_skew_bc)^2)*100", ylim = c(-80, 80))
bias_sk_sk_jack_sv  <- plot_bias_violin(result_skew, "bias_sk_sk_jack_sv",    "Relative Bias $SV_{\\Delta sk}$ (%)", title = "((mean(sk_sv)- sd(jack_skew_bc)^2) / sd(jack_skew_bc)^2)*100", ylim = c(-80, 80))
bias_sk_jack_sk_sv  <- plot_bias_violin(result_skew, "bias_sk_jack_sk_sv",    "Relative Bias $SV_{\\Delta sk}$ (%)", title = "((mean(jack_skew_sv) - sd(sk)^2) / sd(sk)^2)*100", ylim = c(-80, 80))

# Kurtosis
bias_ku_sv       <- plot_bias_violin(result_kurt, "bias_ku_sv",         "Relative Bias $SV_{\\Delta ku}$ (%)", title = "((mean(ku_sv) - sd(ku)^2) / sd(ku)^2)*100", ylim = c(-100, 100)) ## Problem as excludes 40 rows
bias_sv_ku_jack  <- plot_bias_violin(result_kurt, "bias_ku_jack_sv",    "Relative Bias $SV_{\\Delta ku}$ (%)", title = "((mean(jack_ku_sv) - sd(jack_ku_bc)^2) / sd(jack_ku_bc)^2)*100", ylim = c(-100, 100))
bias_ku_ku_jack_sv  <- plot_bias_violin(result_kurt, "bias_ku_ku_jack_sv",    "Relative Bias $SV_{\\Delta ku}$ (%)", title = "((mean(ku_sv) - sd(jack_ku_bc)^2) / sd(jack_ku_bc)^2)*100", ylim = c(-100, 100))
bias_ku_jack_ku_sv  <- plot_bias_violin(result_kurt, "bias_ku_jack_ku_sv",    "Relative Bias $SV_{\\Delta ku}$ (%)", title = "((mean(jack_ku_sv) - sd(ku)^2) / sd(ku)^2)*100", ylim = c(-100, 100))

# Correlation
bias_sv_z        <- plot_bias_violin(result_cor,  "bias_d_cor_sv",      "Relative Bias $SV_{\\Delta Z_{r}}$ (%)", title = "((mean(d_cor_sv) - sd(d_cor)^2) / sd(d_cor)^2)*100", ylim = c(-20, 80))
bias_sv_z_jack   <- plot_bias_violin(result_cor,  "bias_jack_d_cor_sv", "Relative Bias $SV_{\\Delta Z_{r}}$ (%)", title = "(mean(jack_d_cor_sv) - sd(jack_d_cor_bc)^2) / sd(jack_d_cor_bc)^2)*100", ylim = c(-20, 80))
bias_d_cor_jack_sv   <- plot_bias_violin(result_cor,  "bias_d_cor_jack_sv", "Relative Bias $SV_{\\Delta Z_{r}}$ (%)", title = "(mean(jack_d_cor_sv) - sd(d_cor)^2) / sd(d_cor)^2)*100", ylim = c(-20, 80))
bias_jack_d_cor_sv.1   <- plot_bias_violin(result_cor,  "bias_jack_d_cor_sv.1", "Relative Bias $SV_{\\Delta Z_{r}}$ (%)", title = "(mean(d_cor_sv) - sd(jack_d_cor_bc)^2) / sd(jack_d_cor_bc)^2)*100", ylim = c(-20, 80))

# Combine all plots
final_rel_bias_plot <- (bias_sv_sk | bias_sv_sk_jack | bias_sk_sk_jack_sv | bias_sk_jack_sk_sv) / (bias_ku_sv | bias_sv_ku_jack | bias_ku_ku_jack_sv |bias_ku_jack_ku_sv) / (bias_sv_z | bias_sv_z_jack | bias_d_cor_jack_sv | bias_jack_d_cor_sv.1) +
  plot_annotation(tag_levels = 'A', tag_suffix = ")") &
  theme(plot.tag = element_text(size = 16, face = "bold"))
  
# Show the final plot
final_rel_bias_plot
ggsave("./output/figs/relativebias.png", plot = final_rel_bias_plot, width = 23, height = 12)


# Main MS plot

bias_ku_jack_ku_sv2  <- plot_bias_violin(result_kurt, "bias_ku_jack_ku_sv",    "Relative Bias $SV_{\\Delta ku}$ (%)", title = "", ylim = c(-50, 50))
bias_sk_jack_sk_sv2  <- plot_bias_violin(result_skew, "bias_sk_jack_sk_sv",    "Relative Bias $SV_{\\Delta sk}$ (%)", title = "", ylim = c(-50, 50))
bias_sv_z2        <- plot_bias_violin(result_cor,  "bias_d_cor_sv",      "Relative Bias $SV_{\\Delta Z_{r}}$ (%)", title = "", ylim = c(-30, 30))

ms_plot <- (bias_sk_plot| bias_sk_jack_sk_sv2) / (bias_ku_plot | bias_ku_jack_ku_sv2) / (bias_z_plot | bias_sv_z2) + plot_annotation(tag_levels = 'A', tag_suffix = ")") & theme(plot.tag = element_text(size = 16, face = "bold"))
