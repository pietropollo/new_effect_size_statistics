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

##------------------------------------------------------------------------##
## Bias in Estimates
##------------------------------------------------------------------------##
## Skewness
bias_sk_plot <- ggplot(result_skew, aes(x = factor(n), y = bias_sk, fill = factor(n))) + 
  geom_violin() + geom_hline(aes(yintercept = 0), linetype = "dashed", color = "black") +
  labs(x = "Sample Size",
       y = TeX("Bias $\\Delta sk$")) +
  theme_classic() + theme(legend.position = "none", axis.title=element_text(size = 14), axis.text=element_text(size = 12), plot.tag = element_text(size = 16, face = "bold")) + scale_fill_viridis_d() 

bias_sk_plot_delta <- ggplot(result_skew, aes(x = factor(n), y = bias_sk_delta, fill = factor(n))) + 
  geom_violin() + geom_hline(aes(yintercept = 0), linetype = "dashed", color = "black") +
  labs(x = "Sample Size",
       y = TeX("Bias $\\Delta sk$ (Delta Method)")) +
  theme_classic() + theme(legend.position = "none", axis.title=element_text(size = 14), axis.text=element_text(size = 12), plot.tag = element_text(size = 16, face = "bold")) + scale_fill_viridis_d() 

bias_sk_plot_boot <- ggplot(result_skew, aes(x = factor(n), y = bias_sk_boot_bc, fill = factor(n))) + 
  geom_violin() + geom_hline(aes(yintercept = 0), linetype = "dashed", color = "black") +
  labs(x = "Sample Size",
       y = TeX("Bias $\\Delta sk$ (Bootstrap)")) +
  theme_classic() + theme(legend.position = "none", axis.title=element_text(size = 14), axis.text=element_text(size = 12), plot.tag = element_text(size = 16, face = "bold")) + scale_fill_viridis_d() 

bias_sk_plot_jack <- ggplot(result_skew, aes(x = factor(n), y = bias_sk_jack_bc, fill = factor(n))) + 
  geom_violin() + geom_hline(aes(yintercept = 0), linetype = "dashed", color = "black") +
  labs(x = "Sample Size",
       y = TeX("Bias $\\Delta sk$ (Jackknife)")) +
  theme_classic() + theme(legend.position = "none", axis.title=element_text(size = 14), axis.text=element_text(size = 12), plot.tag = element_text(size = 16, face = "bold")) + scale_fill_viridis_d() 

# Kurtosis
bias_kurt_plot <- ggplot(result_kurt, aes(x = factor(n), y = bias_ku, fill = factor(n))) + 
  geom_violin() + geom_hline(aes(yintercept = 0), linetype = "dashed", color = "black") +
  labs(x = "Sample Size",
       y = TeX("Bias $\\Delta ku$")) +
  theme_classic() + theme(legend.position = "none", axis.title=element_text(size = 14), axis.text=element_text(size = 12), plot.tag = element_text(size = 16, face = "bold")) + scale_fill_viridis_d() 

bias_kurt_plot_delta <- ggplot(result_kurt, aes(x = factor(n), y = bias_ku_delta, fill = factor(n))) + 
  geom_violin() + geom_hline(aes(yintercept = 0), linetype = "dashed", color = "black") +
  labs(x = "Sample Size",
       y = TeX("Bias $\\Delta ku$ (Delta Method)")) +
  theme_classic() + theme(legend.position = "none", axis.title=element_text(size = 14), axis.text=element_text(size = 12), plot.tag = element_text(size = 16, face = "bold")) + scale_fill_viridis_d()

bias_kurt_plot_boot <- ggplot(result_kurt, aes(x = factor(n), y = bias_ku_boot_bc, fill = factor(n))) + 
  geom_violin() + geom_hline(aes(yintercept = 0), linetype = "dashed", color = "black") +
  labs(x = "Sample Size",
       y = TeX("Bias $\\Delta ku$ (Bootstrap)")) +
  theme_classic() + theme(legend.position = "none", axis.title=element_text(size = 14), axis.text=element_text(size = 12), plot.tag = element_text(size = 16, face = "bold")) + scale_fill_viridis_d()

bias_kurt_plot_jack <- ggplot(result_kurt, aes(x = factor(n), y = bias_ku_jack_bc, fill = factor(n))) + 
  geom_violin() + geom_hline(aes(yintercept = 0), linetype = "dashed", color = "black") +
  labs(x = "Sample Size",
       y = TeX("Bias $\\Delta ku$ (Jackknife)")) +
  theme_classic() + theme(legend.position = "none", axis.title=element_text(size = 14), axis.text=element_text(size = 12), plot.tag = element_text(size = 16, face = "bold")) + scale_fill_viridis_d()                           


# Correlation
bias_cor_plot <- ggplot(result_cor, aes(x = factor(n), y = bias_d_cor, fill = factor(n))) + 
 geom_violin() + geom_hline(aes(yintercept = 0), linetype = "dashed", color = "black") +
  labs(x = "Sample Size",
       y = TeX("Bias $\\Delta Z_{r}$")) +
  theme_classic() + theme(legend.position = "none", axis.title=element_text(size = 14), axis.text=element_text(size = 12), plot.tag = element_text(size = 16, face = "bold")) + scale_fill_viridis_d() 

bias_cor_plot_boot <- ggplot(result_cor, aes(x = factor(n), y = bias_boot_d_cor, fill = factor(n))) + 
 geom_violin() + geom_hline(aes(yintercept = 0), linetype = "dashed", color = "black") +
  labs(x = "Sample Size",
       y = TeX("Bias $\\Delta Z_{r}$ (Bootstrapped)")) +
  theme_classic() + theme(legend.position = "none", axis.title=element_text(size = 14), axis.text=element_text(size = 12), plot.tag = element_text(size = 16, face = "bold")) + scale_fill_viridis_d()

bias_cor_plot_jack <- ggplot(result_cor, aes(x = factor(n), y = bias_jack_d_cor, fill = factor(n))) + 
 geom_violin() + geom_hline(aes(yintercept = 0), linetype = "dashed", color = "black") +
  labs(x = "Sample Size",
       y = TeX("Bias $\\Delta Z_{r}$ (Jackknife)")) +
  theme_classic() + theme(legend.position = "none", axis.title=element_text(size = 14), axis.text=element_text(size = 12), plot.tag = element_text(size = 16, face = "bold")) + scale_fill_viridis_d()  

est_plot <- (bias_sk_plot | bias_sk_plot_delta | bias_sk_plot_boot | bias_sk_plot_jack) / (bias_kurt_plot | bias_kurt_plot_delta | bias_kurt_plot_boot | bias_kurt_plot_jack) / (bias_cor_plot |  bias_cor_plot_boot | bias_cor_plot_jack) + 
  plot_annotation(tag_levels = 'A', tag_suffix = ")") & 
  theme(plot.tag = element_text(size = 16, face = "bold"))

##------------------------------------------------------------------------##
## Relative Bias in Sampling Error of Estimates
##------------------------------------------------------------------------##

# Skewness

bias_sk_sv_plot <- ggplot(result_skew, aes(x = factor(n), y = bias_sk_sv, fill = factor(n))) + 
 geom_violin() + geom_hline(aes(yintercept = 0), linetype = "dashed", color = "black") +
  labs(x = "Sample Size",
       y = TeX("Relative Bias $SV_{\\Delta sk}$ (%)")) +
  theme_classic() + theme(legend.position = "none", axis.title=element_text(size = 14), axis.text=element_text(size = 12), plot.tag = element_text(size = 16, face = "bold")) + scale_fill_viridis_d()

bias_sk_sv_plot_delta <- ggplot(result_skew, aes(x = factor(n), y = bias_sk_delta_sv, fill = factor(n))) + 
 geom_violin() + geom_hline(aes(yintercept = 0), linetype = "dashed", color = "black") +
  labs(x = "Sample Size",
       y = TeX("Relative Bias $SV_{\\Delta sk}$ (%) (Delta method)")) +
  theme_classic() + theme(legend.position = "none", axis.title=element_text(size = 14), axis.text=element_text(size = 12), plot.tag = element_text(size = 16, face = "bold")) + scale_fill_viridis_d()

bias_sk_sv_plot_boot <- ggplot(result_skew, aes(x = factor(n), y = bias_sk_boot_sv, fill = factor(n))) + 
 geom_violin() + geom_hline(aes(yintercept = 0), linetype = "dashed", color = "black") +
  labs(x = "Sample Size",
       y = TeX("Relative Bias $SV_{\\Delta sk}$ (%) (Bootstrap)")) +
  theme_classic() + theme(legend.position = "none", axis.title=element_text(size = 14), axis.text=element_text(size = 12), plot.tag = element_text(size = 16, face = "bold")) + scale_fill_viridis_d()

bias_sk_sv_plot_jack <- ggplot(result_skew, aes(x = factor(n), y = bias_sk_jack_sv, fill = factor(n))) + 
 geom_violin() + geom_hline(aes(yintercept = 0), linetype = "dashed", color = "black") +
  labs(x = "Sample Size",
       y = TeX("Relative Bias $SV_{\\Delta sk}$ (%) (Jackknife)")) +
  theme_classic() + theme(legend.position = "none", axis.title=element_text(size = 14), axis.text=element_text(size = 12), plot.tag = element_text(size = 16, face = "bold")) + scale_fill_viridis_d()

# Kurtosis
bias_kurt_sv_plot <- ggplot(result_kurt, aes(x = factor(n), y = bias_ku_sv, fill = factor(n))) + 
 geom_violin() + geom_hline(aes(yintercept = 0), linetype = "dashed", color = "black") +
  labs(x = "Sample Size",
       y = TeX("Relative Bias $SV_{\\Delta ku}$ (%)")) +
  theme_classic() + theme(legend.position = "none", axis.title=element_text(size = 14), axis.text=element_text(size = 12), plot.tag = element_text(size = 16, face = "bold")) + scale_fill_viridis_d()

bias_kurt_sv_plot_delta <- ggplot(result_kurt, aes(x = factor(n), y = bias_ku_delta_sv, fill = factor(n))) + 
 geom_violin() + geom_hline(aes(yintercept = 0), linetype = "dashed", color = "black") +
  labs(x = "Sample Size",
       y = TeX("Relative Bias $SV_{\\Delta ku}$ (%) (Delta method)")) +
  theme_classic() + theme(legend.position = "none", axis.title=element_text(size = 14), axis.text=element_text(size = 12), plot.tag = element_text(size = 16, face = "bold")) + scale_fill_viridis_d()

bias_kurt_sv_plot_boot <- ggplot(result_kurt, aes(x = factor(n), y = bias_ku_boot_sv, fill = factor(n))) + 
 geom_violin() + geom_hline(aes(yintercept = 0), linetype = "dashed", color = "black") +
  labs(x = "Sample Size",
       y = TeX("Relative Bias $SV_{\\Delta ku}$ (%) (Bootstrap)")) +
  theme_classic() + theme(legend.position = "none", axis.title=element_text(size = 14), axis.text=element_text(size = 12), plot.tag = element_text(size = 16, face = "bold")) + scale_fill_viridis_d()

bias_kurt_sv_plot_jack <- ggplot(result_kurt, aes(x = factor(n), y = bias_ku_jack_sv, fill = factor(n))) + 
 geom_violin() + geom_hline(aes(yintercept = 0), linetype = "dashed", color = "black") +
  labs(x = "Sample Size",
       y = TeX("Relative Bias $SV_{\\Delta ku}$ (%) (Jackknife)")) +
  theme_classic() + theme(legend.position = "none", axis.title=element_text(size = 14), axis.text=element_text(size = 12), plot.tag = element_text(size = 16, face = "bold")) + scale_fill_viridis_d()

# Correlation
bias_cor_sv_plot <- ggplot(result_cor, aes(x = factor(n), y = bias_d_cor_sv, fill = factor(n))) + 
 geom_violin() + geom_hline(aes(yintercept = 0), linetype = "dashed", color = "black") +
  labs(x = "Sample Size",
       y = TeX("Relative Bias $SV_{\\Delta Z_{r}}$ (%)")) +
  theme_classic() + theme(legend.position = "none", axis.title=element_text(size = 14), axis.text=element_text(size = 12), plot.tag = element_text(size = 16, face = "bold")) + scale_fill_viridis_d()

bias_cor_sv_boot_plot <- ggplot(result_cor, aes(x = factor(n), y = bias_boot_d_cor_sv, fill = factor(n))) + 
 geom_violin() + geom_hline(aes(yintercept = 0), linetype = "dashed", color = "black") +
  labs(x = "Sample Size",
       y = TeX("Relative Bias $SV_{\\Delta Z_{r}}$ (%) (Bootstrapped)")) +
  theme_classic() + theme(legend.position = "none", axis.title=element_text(size = 14), axis.text=element_text(size = 12), plot.tag = element_text(size = 16, face = "bold")) + scale_fill_viridis_d()  

bias_cor_sv_jacknife_plot <- ggplot(result_cor, aes(x = factor(n), y = bias_jack_d_cor_sv, fill = factor(n))) + 
 geom_violin() + geom_hline(aes(yintercept = 0), linetype = "dashed", color = "black") +
  labs(x = "Sample Size",
       y = TeX("Relative Bias $SV_{\\Delta Z_{r}}$ (%) (Jackknife)")) +
  theme_classic() + theme(legend.position = "none", axis.title=element_text(size = 14), axis.text=element_text(size = 12), plot.tag = element_text(size = 16, face = "bold")) + scale_fill_viridis_d() 

# Combine all plots
final_rel_bias_plot <- (bias_sk_sv_plot | bias_sk_sv_plot_delta | bias_sk_sv_plot_boot | bias_sk_sv_plot_jack) / (bias_kurt_sv_plot | bias_kurt_sv_plot_delta | bias_kurt_sv_plot_boot | bias_kurt_sv_plot_jack) / (bias_cor_sv_plot | bias_cor_sv_boot_plot | bias_cor_sv_jacknife_plot) + plot_annotation(tag_levels = 'A', tag_suffix = ")")

# Show the final plot
print(final_rel_bias_plot)
