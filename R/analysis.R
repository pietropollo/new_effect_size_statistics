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
result_kurt <- readRDS("./output/result_kurt.rds") # Load the kurtosis results, has adjustements (x1_kurt_jack$var + x2_kurt_jack$var) + ku[i]^2 / (2*(params$n + params$n -2))
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
                                ylim = c(-0.25,
                                         0.25))
bias_z_plot_jack <- plot_bias_violin(result_cor, 
                                     "bias_jack_d_cor",
                                     "Bias $\\Delta \\textit{Zr}$ (jackknife)", 
                                     ylim = c(-0.25, 
                                              0.25))

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

list_sk <- c("bias_sk_sv",
               "bias_sk_jack_sv",
               "bias_sk_sk_jack_sv",
               "bias_sk_jack_sk_sv")

plot_labels <- c("((mean(sk_sv) - sd(sk)^2) / sd(sk)^2)*100",
                 "((mean(jack_skew_sv) - sd(jack_skew_bc)^2) / sd(jack_skew_bc)^2)*100",
                 "((mean(sk_sv) - sd(jack_skew_bc)^2) / sd(jack_skew_bc)^2)*100",
                 "((mean(jack_skew_sv) - sd(sk)^2) / sd(sk)^2)*100")

list_sk_plots <- list()
for (i in seq_along(list_sk)) {
  list_sk_plots[[i]] <- plot_bias_violin(result_skew,
                                          list_sk[i],
   "Relative bias $\\textit{s^2}_{\\Delta \\textit{sk}}$ (%)",
                                        title = plot_labels[i],
                                        ylim = c(-80,
                                                 80))
  names(list_sk_plots)[i] <- list_sk[i]
}

## Kurtosis ----

list_ku  <- c("bias_ku_sv",
               "bias_ku_jack_sv",
               "bias_ku_ku_jack_sv",
               "bias_ku_jack_ku_sv") # replace with adjusted jackknife relative bias "adj_jack_sv_w_ku_rel_bias" if needed

plot_labels <- c("((mean(ku_sv) - sd(ku)^2) / sd(ku)^2)*100",
                 "((mean(jack_ku_sv) - sd(jack_ku_bc)^2) / sd(jack_ku_bc)^2)*100",
                 "((mean(ku_sv) - sd(jack_ku_bc)^2) / sd(jack_ku_bc)^2)*100",
                 "((mean(jack_ku_sv) - sd(ku)^2) / sd(ku)^2)*100")

list_ku_plots <- list()
for (i in seq_along(list_ku)) {
  list_ku_plots[[i]] <- plot_bias_violin(result_kurt,
                                        list_ku[i],
                                        "Relative bias $\\textit{s^2}_{\\Delta \\textit{ku}}$ (%)",
                                        title = plot_labels[i],
                                        ylim = c(-100,
                                                 100))
  names(list_ku_plots)[i] <- list_ku[i]
}

## Correlation ----

list_cor  <- c("bias_d_cor_sv",
                  "bias_jack_d_cor_sv",
                  "bias_d_cor_jack_sv",
                  "bias_jack_d_cor_sv.1")

plot_labels <- c("((mean(d_cor_sv) - sd(d_cor)^2) / sd(d_cor)^2)*100",
                 "(mean(jack_d_cor_sv) - sd(jack_d_cor_bc)^2) / sd(jack_d_cor_bc)^2)*100",
                 "(mean(jack_d_cor_sv) - sd(d_cor)^2) / sd(d_cor)^2)*100",
                 "(mean(d_cor_sv) - sd(jack_d_cor_bc)^2) / sd(jack_d_cor_bc)^2)*100")

list_cor_plots <- list()
for (i in seq_along(list_cor)) {
  list_cor_plots[[i]] <- plot_bias_violin(result_cor,
                                         list_cor[i],
                                         "Relative bias $\\textit{s^2}_{\\Delta \\textit{Zr}}$ (%)",
                                         title = plot_labels[i],
                                         ylim = c(-80,
                                                  80))
  names(list_cor_plots)[i] <- list_cor[i]
}


## Combine ----
final_rel_bias_plot <- 
  list_sk_plots[["bias_sk_sv"]] + 
  list_sk_plots[["bias_sk_jack_sv"]] +
  list_sk_plots[["bias_sk_sk_jack_sv"]] + 
  list_sk_plots[["bias_sk_jack_sk_sv"]] +
  list_ku_plots[["bias_ku_sv"]] +
  list_ku_plots[["bias_ku_jack_sv"]] + 
  list_ku_plots[["bias_ku_ku_jack_sv"]] +
  list_ku_plots[["bias_ku_jack_ku_sv"]] +
  list_cor_plots[["bias_d_cor_sv"]] + 
  list_cor_plots[["bias_jack_d_cor_sv"]] + 
  list_cor_plots[["bias_d_cor_jack_sv"]] +
  list_cor_plots[["bias_jack_d_cor_sv.1"]] +
  plot_layout(nrow = 3, ncol = 4) +
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

# Coverage plots ----

## Skewness ----
lims_cov <- c(0, 1)

coverage_sk  <- c("coverage_sk_all",
                  "coverage_sk_jack_bc_all",
                  "coverage_sk_jack_sv_all",
                  "coverage_jack_sk_sv_all")

plot_labels_cov <- c("Traditional (estimator & Sampling Variance)",
                     "Jackknife (estimator & Sampling Variance)",
                     "Traditional estimator with Jackknife Sampling Variance",
                     "Jackknife estimator with Traditional Sampling Variance")

list_cov_skew_plots <- list()
for(i in seq_along(coverage_sk)) {
  list_cov_skew_plots[[i]] <- plot_bias_violin(result_skew,
                          coverage_sk[i],
                          "Coverage (%)",
                          title = plot_labels_cov[i],
                          ylim = lims_cov) + geom_hline(aes(yintercept = 0.95),
                                                linetype = "dashed",
                                                color = "red")
        names(list_cov_skew_plots)[i]  <- coverage_sk[i]
}

## Kurtosis ----
coverage_ku  <- c("coverage_ku_all",
                  "coverage_ku_jack_bc_all",
                  "coverage_ku_jack_sv_all",
                  "coverage_jack_ku_sv_all",
                  "coverage_adj_ku_sv_all",
                  "coverage_adj_jack_ku_sv_all", 
                  "boot_coverage_ku_all")

plot_labels_cov_ku <- c("Traditional (estimator & Sampling Variance)",
                         "Jackknife (estimator & Sampling Variance)",
                         "Traditional estimator with Jackknife Sampling Variance",
                         "Jackknife estimator with Traditional Sampling Variance",
                         "Traditional estimator & Adjusted Traditional Sampling Variance",
                         "Jackknife estimator & Adjusted Jackknife Sampling Variance",
                         "Bootstrap BCA 95% CI")

list_cov_kurt_plots <- list()

for(i in seq_along(coverage_ku)) {
  list_cov_kurt_plots[[i]] <- plot_bias_violin(result_kurt,
                                              coverage_ku[i],
                                              "Coverage (%)",
                                              title = plot_labels_cov_ku[i],
                                              ylim = lims_cov) + geom_hline(aes(yintercept = 0.95),
                                                                    linetype = "dashed",
                                                                    color = "red")
  names(list_cov_kurt_plots)[i]  <- coverage_ku[i]
}

## Correlation ----

coverage_cor  <- c("coverage_d_cor_all",
                   "coverage_d_cor_jack_bc_all",
                   "coverage_jack_bc_sv_all",
                   "coverage_bc_jack_sv_all")

plot_labels_cov_cor <- c("Traditional (estimator & Sampling Variance)",
                          "Jackknife (estimator & Sampling Variance)",
                          "Traditional estimator with Jackknife Sampling Variance",
                          "Jackknife estimator with Traditional Sampling Variance")

list_cov_cor_plots <- list()
for(i in seq_along(coverage_cor)) {
  list_cov_cor_plots[[i]] <- plot_bias_violin(result_cor,
                                              coverage_cor[i],
                                              "Coverage (%)",
                                              title = plot_labels_cov_cor[i],
                                              ylim = lims_cov) + geom_hline(aes(yintercept = 0.95),
                                                                    linetype = "dashed",
                                                                    color = "red")
  names(list_cov_cor_plots)[i]  <- coverage_cor[i]
}


                  
## Combine ----
final_coverage_plot <- 
  list_cov_skew_plots[["coverage_sk_all"]] +
  list_cov_skew_plots[["coverage_sk_jack_bc_all"]] +    
  list_cov_skew_plots[["coverage_sk_jack_sv_all"]] + 
  list_cov_skew_plots[["coverage_jack_sk_sv_all"]] +
  list_cov_kurt_plots[["coverage_ku_all"]] +
  list_cov_kurt_plots[["coverage_ku_jack_bc_all"]] +
  list_cov_kurt_plots[["coverage_ku_jack_sv_all"]] + # "coverage_adj_jack_ku_sv_all" & "coverage_ku_jack_sv_all"
  list_cov_kurt_plots[["boot_coverage_ku_all"]] + # "coverage_adj_ku_sv_all" & "coverage_jack_ku_sv_all" & coverage_adj_ku_sv_all
  list_cov_cor_plots[["coverage_d_cor_all"]] +
  list_cov_cor_plots[["coverage_d_cor_jack_bc_all"]] +
  list_cov_cor_plots[["coverage_jack_bc_sv_all"]] +
  list_cov_cor_plots[["coverage_bc_jack_sv_all"]] +
  plot_layout(nrow = 3, ncol = 4) +
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(size = 12))

# Final plot ----
ggsave("./output/figs/coverage.png",
       plot = final_coverage_plot,
       bg = "white",
       dpi = 600,
       width = 16,
       height = 7,
       units = "in")

#### COVERAGE ACROSS SCENARIOS ######

results_abs_ku2_N <- readRDS("./output/result_kurt_abs(ku)^2_N.rds") # Load the kurtosis results, has adjustements (x1_kurt_jack$var + x2_kurt_jack$var) + ku[i]^2 / (2*(params$n + params$n -2))
    results_ku2_N <- readRDS("./output/result_kurt_ku^2_N.rds") # Load the kurtosis results, has adjustements (x1_kurt_jack$var + x2_kurt_jack$var) + ku[i]^2 / (2*(params$n + params$n -2))
     results_ku2 <- readRDS("./output/result_kurt_ku^2.rds") # Load the kurtosis results, has adjustements (x1_kurt_jack$var + x2_kurt_jack$var) + ku[i]^2 / (2*(params$n + params$n -2))

ku_data_list <- list(results_ku2,
                  results_abs_ku2_N,
                     results_ku2_N,
                     results_ku2)

coverage_ku_adj  <- c("coverage_ku_jack_sv_all",
                  "coverage_adj_jack_ku_sv_all",
                  "coverage_adj_jack_ku_sv_all",
                  "coverage_adj_jack_ku_sv_all")

plot_labels_cov_ku_adj <- c("(x1_kurt_jack$var + x2_kurt_jack$var)",
                         "(x1_kurt_jack$var + x2_kurt_jack$var) + abs(ku[i])^2 / (2*(params$n + params$n -2))",
                         "(x1_kurt_jack$var + x2_kurt_jack$var) + ku[i]^2 / (2*(params$n + params$n -2))",
                         "(x1_kurt_jack$var + x2_kurt_jack$var) + ku[i]^2")

list_cov_kurt_adj_plots <- list()

for(i in seq_along(coverage_ku_adj)) {
  list_cov_kurt_adj_plots[[i]] <- plot_bias_violin(ku_data_list[[i]],
                                              coverage_ku_adj[i],
                                              "Coverage (%)",
                                              title = plot_labels_cov_ku_adj[i],
                                              ylim = lims_cov) + geom_hline(aes(yintercept = 0.95),
                                                                    linetype = "dashed",
                                                                    color = "red")
}

ku_coverage <- list_cov_kurt_adj_plots[[1]] +
  list_cov_kurt_adj_plots[[2]] +
  list_cov_kurt_adj_plots[[3]] + 
  list_cov_kurt_adj_plots[[4]] + 
  plot_layout(nrow = 2, ncol = 2) +
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(size = 12))                  



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

coverage_sk_jack_sv_all <-
  plot_bias_violin(result_skew,
                   "coverage_sk_jack_sv_all",
                   "Coverage in $\\Delta \\textit{sk}$",
                   ylim = c(0,
                            1)) + geom_hline(aes(yintercept = 0.95),
                                              linetype = "dashed",
                                              color = "red")

coverage_ku_jack_sv_all <-
  plot_bias_violin(result_kurt,
                   "coverage_ku_jack_sv_all",
                   "Coverage in $\\Delta \\textit{ku}$",
                   ylim = c(0,
                            1)) + geom_hline(aes(yintercept = 0.95),
                                              linetype = "dashed",
                                              color = "red")  

coverage_bc_jack_sv_all <-
  plot_bias_violin(result_cor,
                   "coverage_bc_jack_sv_all",
                   "Coverage in $\\Delta \\textit{Zr}$",
                   ylim = c(0,
                            1)) + geom_hline(aes(yintercept = 0.95),
                                              linetype = "dashed",
                                              color = "red")                                            
ms_plot <- 
  bias_sk_plot + 
  bias_sk_jack_sk_sv2 +
  coverage_sk_jack_sv_all +
  bias_ku_plot + 
  bias_ku_jack_ku_sv2 +
  coverage_ku_jack_sv_all +
  bias_z_plot + 
  bias_sv_z2 +
  coverage_bc_jack_sv_all +
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

# Exploring kurtosis results ----

# First, is bias, coverage and relative bias related to mean differences?

result_kurt <- result_kurt  %>% mutate(diff = ifelse(mean_g1 == mean_g2, 0, abs(mean_g1 - mean_g2)),
                                      var_diff = ifelse(variance_g1 == variance_g2, 0, abs(variance_g1 - variance_g2)),
                                      kurt_diff = ifelse(kurtosis_g1 == kurtosis_g2, 0, abs(kurtosis_g1 - kurtosis_g2)),
                                      abs_ku_est = abs(ku_est))

result_skew <- result_skew  %>% mutate(diff = ifelse(mean_g1 == mean_g2, 0, abs(mean_g1 - mean_g2)),
                                      var_diff = ifelse(variance_g1 == variance_g2, 0, abs(variance_g1 - variance_g2)),
                                      sk_diff = ifelse(skewness_g1 == skewness_g2, 0, abs(skewness_g1 - skewness_g2)),
                                      abs_sk_est = abs(sk_diff))                                      

model <- lm(coverage_ku_jack_sv_all ~ abs_ku_est + log(1/n) + abs_ku_est*log(1/n), data = result_kurt) 
summary(model)
hist(residuals(model)) # Check normality of residuals

cor_factors <- coefficients(model)

# When the groups differ in their means does bias vary?
mean_diff_ku <- result_kurt %>%
  ggplot() +
  geom_violin(aes(x = diff, y = bias_ku, group = diff, fill = diff)) + labs(x = "Mean difference between groups", y = TeX("Bias $\\Delta \\textit{ku}$")) + theme_classic() + theme(legend.position = "none")

# When groups differ in their variances does bias vary?
var_diff_ku <- result_kurt %>%
  ggplot() +
  geom_violin(aes(x = var_diff, y = bias_ku, group = var_diff, fill = var_diff)) + labs(x = "Variance difference between groups", y = TeX("Bias $\\Delta \\textit{ku}$")) + theme_classic() + theme(legend.position = "none")

mean_diff_sk  <- result_skew %>%
  ggplot() +
  geom_violin(aes(x = diff, y = bias_sk, group = diff, fill = diff)) + labs(x = "Mean difference between groups", y = TeX("Bias $\\Delta \\textit{sk}$")) + theme_classic() + theme(legend.position = "none") 

var_diff_sk  <- result_skew %>%
  ggplot() +
  geom_violin(aes(x = var_diff, y = bias_sk, group = var_diff, fill = var_diff)) + labs(x = "Variance difference between groups", y = TeX("Bias $\\Delta \\textit{sk}$")) + theme_classic() + theme(legend.position = "none")

combined_mean_var <- (mean_diff_sk +
  var_diff_sk)/ (mean_diff_ku +
  var_diff_ku) +
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(size = 12))
ggsave("./output/figs/combined_mean_var.png", plot = combined_mean_var, width = 10, height = 5)

# Plot difference in kurtosis when groups are either DIFERENT or NOT and when the sample size is SMALL or LARGE
kurt_explore<- result_kurt %>% mutate(diff = ifelse(kurtosis_g1 == kurtosis_g2, "No", "Yes")) %>% filter(n %in% c(10, 500)) %>% mutate(n_diff = interaction(as.factor(n), as.factor(diff)))  %>% 
  ggplot() +
  geom_violin(aes(y = kurt_diff, x = n_diff, group = n_diff, fill = n_diff)) + theme_classic() + labs(y = TeX("$\\Delta \\textit{ku}$"), x = TeX("Scenario (no true diff vs true diff) and Sample Size (n)")) + theme(legend.position = "none") 

# When groups differ in their kurtosis does bias vary?
# Use the estimated Kurtosis difference **IMPORTANT PLOT**
ku_est_coverage <- result_kurt %>%
  ggplot() +
  geom_point(aes(x = abs(ku_est), y = coverage_ku_jack_sv_all, col = kurt_diff, size = n)) + labs(x = TeX("Absolute Estimated $\\Delta \\textit{ku}$"), y = TeX("Coverage in $\\Delta \\textit{ku}$"), col = TeX("True Absolute $\\Delta \\textit{ku}$")) + theme_classic() + theme(legend.position = "right", axis.title = element_text(size = 14)) + geom_hline(aes(yintercept = 0.95),
             linetype = "dashed",
             color = "red")

ku_est_coverage_var <- result_kurt %>%
  ggplot() +
  geom_point(aes(x = abs(ku_est), y = coverage_adj_ku_sv_all, col = var_diff, size = n)) + labs(x = TeX("Absolute Estimated $\\Delta \\textit{ku}$"), y = TeX("Coverage in $\\Delta \\textit{ku}$"), col = TeX("True Absolute $\\Delta \\textit{v}$")) + theme_classic() + theme(legend.position = "right", axis.title = element_text(size = 14)) + geom_hline(aes(yintercept = 0.95),
             linetype = "dashed",
             color = "red")             

explore_ku <- ku_est_coverage + ku_est_coverage_var +
  plot_layout(nrow = 1) +
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(size = 12))
# -----------------------

kurt_diff_bias<- result_kurt %>%
  ggplot() +
  geom_violin(aes(x = kurt_diff, y = bias_ku_jack_bc, group = kurt_diff, fill = kurt_diff)) + labs(x = "Absolute Kurtosis difference between groups", y = "Bias in kurtosis estimate") + theme_classic() + theme(legend.position = "none") + geom_hline(aes(yintercept = 0),
             linetype = "dashed",
             color = "black")

kurt_diff_relbias<- result_kurt %>%
  ggplot() +
  geom_violin(aes(x = kurt_diff, y = bias_ku_jack_ku_sv, group = kurt_diff, fill = kurt_diff)) + labs(x = "Absolute Kurtosis difference between groups", y = "Relative Bias") + theme_classic() + theme(legend.position = "none") + geom_hline(aes(yintercept = 0),
             linetype = "dashed",
             color = "black")

final_kurtosis_explore <- kurt_diff_bias +
  kurt_diff_relbias +
  kurt_diff_coverage +
  plot_layout(nrow = 1) +
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(size = 12))
