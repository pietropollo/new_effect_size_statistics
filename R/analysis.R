###------------------------------------------------------------------------###
# Simulation Results for New Effect Size Statistics
# Description: Plotting the results of the simulations for new effect size statistics.
###------------------------------------------------------------------------###

# Clean working environment ----
rm(list = ls()) # Remove all objects from the environment

# Load required packages ----
pacman::p_load(moments, 
               PearsonDS, 
               tidyverse, 
               patchwork, 
               latex2exp)  
source("./R/func.R") # Load the functions from func.R

# Load the simulation results ----
result_kurt <- readRDS("./output/result_kurt.rds") # Load the kurtosis results
result_cor <- readRDS("./output/result_cor.rds") # Load the correlation results	
result_skew <- readRDS("./output/result_skewness.rds") # Load the skewness results

# Best estimators ----
## Skewness ----
result_skew_least_bias <- 
  result_skew %>% 
  summarise(across(starts_with("bias_"), 
                   ~ sum(.x^2, 
                         na.rm = TRUE),
                   .names = "least_bias_{col}")) %>%
  t() %>%
  data.frame() %>%
  arrange() 

colnames(result_skew_least_bias) <- "Least bias"

## Kurtosis ----
result_kurt_least_bias <- 
  result_kurt %>% 
  summarise(across(starts_with("bias_"), 
                   ~ sum(.x^2, 
                         na.rm = TRUE), 
                   .names = "least_bias_{col}")) %>%
  t() %>%
  data.frame() %>%
  arrange()

colnames(result_kurt_least_bias) <- "Least bias"

## Correlation ----
result_cor_least_bias <- 
  result_cor %>%
  summarise(across(starts_with("bias_"),
                   ~ sum(.x^2, 
                         na.rm = TRUE),
                   .names = "least_bias_{col}")) %>%
  t() %>%
  data.frame() %>%
  arrange()

colnames(result_cor_least_bias) <- "Least bias"

# Plot function for the results ----
plot_bias_violin <- 
  function(data, 
           y_var, 
           y_lab, 
           title = "", 
           ylim = c(-1, 
                    1)) {
    ggplot(data,
           aes(x = factor(n),
               y = .data[[y_var]],
               fill = factor(n))) + 
      ylim(ylim) +
      geom_violin() + 
      geom_hline(aes(yintercept = 0),
                 linetype = "dashed",
                 color = "black") +
      labs(x = "Sample size",
           y = TeX(y_lab),
           title = TeX(title)) +
      scale_fill_viridis_d() +
      theme_classic() + 
      theme(legend.position = "none",
            axis.title = element_text(size = 12),
            axis.text = element_text(size = 10),
            title = element_text(size = 6))
  }

##------------------------------------------------------------------------##
## Bias in Estimates
##------------------------------------------------------------------------##

# Plots bias ----
## Skewness ----
bias_sk_plot <- 
  plot_bias_violin(result_skew,
                   "bias_sk",
                   "Bias $\\Delta \\textit{sk}$",
                   ylim = c(-0.5,
                            0.5))
bias_sk_plot_jack <- 
  plot_bias_violin(result_skew, "bias_sk_jack_bc", 
                   "Bias $\\Delta \\textit{sk}$ (jackknife)", 
                   ylim = c(-0.5, 
                            0.5))

## Kurtosis ----
bias_ku_plot <-
  plot_bias_violin(result_kurt,
                   "bias_ku",
                   "Bias $\\Delta \\textit{ku}$",
                   ylim = c(-3, 
                            3))
bias_ku_plot_jack <- plot_bias_violin(result_kurt,
                                      "bias_ku_jack_bc",
                                      "Bias $\\Delta \\textit{ku}$ (jackknife)",
                                      ylim = c(-3, 
                                               3))

## Correlation ----
bias_z_plot <- plot_bias_violin(result_cor,
                                "bias_d_cor",
                                "Bias $\\Delta \\textit{Zr}$",
                                ylim = c(-0.1,
                                         0.1))
bias_z_plot_jack <- plot_bias_violin(result_cor, 
                                     "bias_jack_d_cor",
                                     "Bias $\\Delta \\textit{Zr}$ (jackknife)", 
                                     ylim = c(-0.1, 
                                              0.1))

## Combine ----
est_plot <- 
  bias_sk_plot +
  bias_sk_plot_jack +
  bias_ku_plot +
  bias_ku_plot_jack +
  bias_z_plot +
  bias_z_plot_jack +
  plot_layout(nrow = 3) +
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(size = 12))

ggsave("./output/figs/bias.png", 
       plot = est_plot, 
       bg = "white",
       dpi = 600,
       width = 9,
       height = 7,
       units = "in")

##------------------------------------------------------------------------##
## Relative bias in Sampling Error of Estimates
##------------------------------------------------------------------------##

# Plots sampling error ----
## Skewness ----
bias_sv_sk <- 
  plot_bias_violin(result_skew,
                   "bias_sk_sv",
                   "Relative bias $\\textit{s^2}_{\\Delta \\textit{sk}}$ (%)",
                   title = "((mean(sk_sv)- sd(sk)^2) / sd(sk)^2)*100",
                   ylim = c(-80,
                            80))

bias_sv_sk_jack <- 
  plot_bias_violin(result_skew,
                   "bias_sk_jack_sv",
                   "Relative bias $\\textit{s^2}_{\\Delta \\textit{sk}}$ (%)",
                   title = "((mean(jack_skew_sv) - sd(jack_skew_bc)^2) / sd(jack_skew_bc)^2)*100",
                   ylim = c(-80,
                            80))

bias_sk_sk_jack_sv <- 
  plot_bias_violin(result_skew,
                   "bias_sk_sk_jack_sv",
                   "Relative bias $\\textit{s^2}_{\\Delta \\textit{sk}}$ (%)",
                   title = "((mean(sk_sv)- sd(jack_skew_bc)^2) / sd(jack_skew_bc)^2)*100",
                   ylim = c(-80, 
                            80))

bias_sk_jack_sk_sv <- 
  plot_bias_violin(result_skew,
                   "bias_sk_jack_sk_sv",
                   "Relative bias $\\textit{s^2}_{\\Delta \\textit{sk}}$ (%)",
                   title = "((mean(jack_skew_sv) - sd(sk)^2) / sd(sk)^2)*100",
                   ylim = c(-80,
                            80))

## Kurtosis ----
bias_ku_sv <- 
  plot_bias_violin(result_kurt,
                   "bias_ku_sv",
                   "Relative bias $\\textit{s^2}_{\\Delta \\textit{ku}}$ (%)",
                   title = "((mean(ku_sv) - sd(ku)^2) / sd(ku)^2)*100",
                   ylim = c(-100,
                            100)) ## Problem as excludes 40 rows

bias_sv_ku_jack <- 
  plot_bias_violin(result_kurt,
                   "bias_ku_jack_sv",
                   "Relative bias $\\textit{s^2}_{\\Delta \\textit{ku}}$ (%)",
                   title = "((mean(jack_ku_sv) - sd(jack_ku_bc)^2) / sd(jack_ku_bc)^2)*100",
                   ylim = c(-100,
                            100))

bias_ku_ku_jack_sv <- 
  plot_bias_violin(result_kurt,
                   "bias_ku_ku_jack_sv",
                   "Relative bias $\\textit{s^2}_{\\Delta \\textit{ku}}$ (%)",
                   title = "((mean(ku_sv) - sd(jack_ku_bc)^2) / sd(jack_ku_bc)^2)*100",
                   ylim = c(-100,
                            100))

bias_ku_jack_ku_sv <- 
  plot_bias_violin(result_kurt, 
                   "bias_ku_jack_ku_sv",
                   "Relative bias $\\textit{s^2}_{\\Delta \\textit{ku}}$ (%)",
                   title = "((mean(jack_ku_sv) - sd(ku)^2) / sd(ku)^2)*100",
                   ylim = c(-100,
                            100))

## Correlation ----
bias_sv_z <-
  plot_bias_violin(result_cor,
                   "bias_d_cor_sv",
                   "Relative bias $\\textit{s^2}_{\\Delta \\textit{Zr}}$ (%)",
                   title = "((mean(d_cor_sv) - sd(d_cor)^2) / sd(d_cor)^2)*100",
                   ylim = c(-20,
                            80))
bias_sv_z_jack <-
  plot_bias_violin(result_cor,
                   "bias_jack_d_cor_sv",
                   "Relative bias $\\textit{s^2}_{\\Delta \\textit{Zr}}$ (%)",
                   title = "(mean(jack_d_cor_sv) - sd(jack_d_cor_bc)^2) / sd(jack_d_cor_bc)^2)*100",
                   ylim = c(-20, 
                            80))

bias_d_cor_jack_sv <-
  plot_bias_violin(result_cor,
                   "bias_d_cor_jack_sv",
                   "Relative bias $\\textit{s^2}_{\\Delta \\textit{Zr}}$ (%)",
                   title = "(mean(jack_d_cor_sv) - sd(d_cor)^2) / sd(d_cor)^2)*100", 
                   ylim = c(-20,
                            80))

bias_jack_d_cor_sv.1 <-
  plot_bias_violin(result_cor,
                   "bias_jack_d_cor_sv.1", 
                   "Relative bias $\\textit{s^2}_{\\Delta \\textit{Zr}}$ (%)",
                   title = "(mean(d_cor_sv) - sd(jack_d_cor_bc)^2) / sd(jack_d_cor_bc)^2)*100",
                   ylim = c(-20, 
                            80))

## Combine ----
final_rel_bias_plot <- 
  bias_sv_sk + 
  bias_sv_sk_jack +
  bias_sk_sk_jack_sv + 
  bias_sk_jack_sk_sv +
  bias_ku_sv +
  bias_sv_ku_jack +
  bias_ku_ku_jack_sv +
  bias_ku_jack_ku_sv +
  bias_sv_z +
  bias_sv_z_jack +
  bias_d_cor_jack_sv +
  bias_jack_d_cor_sv.1 +
  plot_layout(nrow = 3) +
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(size = 12))

# Final plot ----
ggsave("./output/figs/relativebias.png",
       plot = final_rel_bias_plot,
       bg = "white",
       dpi = 600,
       width = 16,
       height = 7,
       units = "in")

# Main MS plot ----
bias_ku_jack_ku_sv2 <-
  plot_bias_violin(result_kurt,
                   "bias_ku_jack_ku_sv",
                   "Relative bias $\\textit{s^2}_{\\Delta \\textit{ku}}$ (%)",
                   ylim = c(-50,
                            50))

bias_sk_jack_sk_sv2 <- 
  plot_bias_violin(result_skew,
                   "bias_sk_jack_sk_sv",
                   "Relative bias $\\textit{s^2}_{\\Delta \\textit{sk}}$ (%)",
                   ylim = c(-50,
                            50))

bias_sv_z2 <-
  plot_bias_violin(result_cor,
                   "bias_d_cor_sv",
                   "Relative bias $\\textit{s^2}_{\\Delta \\textit{Zr}}$ (%)",
                   ylim = c(-30,
                            30))

ms_plot <- 
  bias_sk_plot + 
  bias_sk_jack_sk_sv2 +
  bias_ku_plot + 
  bias_ku_jack_ku_sv2 +
  bias_z_plot + 
  bias_sv_z2 +
  plot_layout(nrow = 3) +
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(size = 12))

ggsave("./output/figs/ms_plot.png",
       plot = ms_plot,
       bg = "white",
       dpi = 600,
       width = 10,
       height = 7,
       units = "in")