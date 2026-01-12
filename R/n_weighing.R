# test

# Clean working environment ----
rm(list = ls()) # Remove all objects from the environment


# latest ochaRd - 

#install.packages("pacman")
#pacman::p_load(devtools, tidyverse, metafor, patchwork, R.rsp, emmeans)

#devtools::install_github("daniel1noble/orchaRd", force = TRUE)
#library(orchaRd)

# Load required packages ----
pacman::p_load(moments, 
               PearsonDS,
               tidyverse,
               patchwork,
               latex2exp,
               metafor,
               here,
               orchaRd)
source("./R/func.R") # Load the functions from func.R


# extra function

# calculate established effect sizes ----
calc.effect <- function(data = raw_data, 
                        m) { # calculates other already established effect size statistics 
  escalc(measure = m,
         m1i = data$mean_male,
         m2i = data$mean_female,
         sd1i = data$sd_male,
         sd2i = data$sd_female,
         n1i = data$n_male,
         n2i = data$n_female,
         var.names = c(paste0(m,
                              "_est"),
                       paste0(m,
                              "_var")))
}

# read csv mice_data_sample.csv

# raw data ----
df_raw <- 
  read_csv("mice_data_sample.csv") %>% 
  # small adjustments to make plots more readable:
  mutate(phenotyping_center = 
           ifelse(phenotyping_center == "MRC Harwell",
                  "MRC H",
                  phenotyping_center),
         strain_fig = case_when(strain_accession_id == "MGI:2159965" ~ 
                                  "N",
                                strain_accession_id == "MGI:2683688" ~ 
                                  "NCrl",
                                strain_accession_id == "MGI:2164831" ~ 
                                  "NTac",
                                strain_accession_id == "MGI:3056279" ~ 
                                  "NJ",
                                strain_accession_id == "MGI:2160139" ~ 
                                  "NJcl"))

df_raw_wide <-
  df_raw %>% 
  select(- strain_accession_id) %>% 
  pivot_wider(id_cols = c(specimen_id,
                          trait_name,
                          phenotyping_center,
                          strain_fig),
              names_from = sex,
              values_from = value)

df_meta_analysed <- # takes 30s to a minute to run
  df_raw %>% 
  group_by(sex,
           trait_name,
           phenotyping_center,
           strain_fig) %>% 
  summarize(mean = mean(value,
                        na.rm = T),
            sd = sd(value,
                    na.rm = T),
            n = n()) %>% 
  pivot_wider(id_cols = c(trait_name,
                          phenotyping_center,
                          strain_fig),
              names_from = sex,
              values_from = c(mean:n)) %>% 
  bind_cols(calc.effect(., m = "ROM")) %>% # lnRR
  bind_cols(calc.effect(., m = "CVR")) %>% # lnCVR
  bind_cols(calc.effect(., m = "VR")) %>%  # lnVR
  left_join(df_raw_wide %>% 
              group_by(trait_name,
                       phenotyping_center,
                       strain_fig) %>% 
              reframe(orchaRd::moment_effects(x1 = na.omit(male),
                                              x2 = na.omit(female),
                                              type = "skew"))) %>% 
  left_join(df_raw_wide %>% 
              group_by(trait_name,
                       phenotyping_center,
                       strain_fig) %>% 
              reframe(orchaRd::moment_effects(x1 = na.omit(male),
                                              x2 = na.omit(female),
                                              type = "kurt"))) %>% 
  rename(SK_delta_est = d_skew,
         SK_delta_var = d_skew_v,
         KU_delta_est = d_kurt,
         KU_delta_var = d_kurt_v) %>% 
  mutate(n_total = n_female + n_male,
         prop_females = n_female / (n_female + n_male)) %>% 
  select(trait_name,
         phenotyping_center,
         strain_fig,
         n_total,
         n_male,
         n_female,
         prop_females,
         ROM_est,
         ROM_var,
         CVR_est,
         CVR_var,
         VR_est,
         VR_var,
         SK_delta_est,
         SK_delta_var,
         KU_delta_est,
         KU_delta_var) %>% 
  mutate(ROM_upper = ROM_est + qt(0.975, 
                                  n_total - 1) * sqrt(ROM_var),
         ROM_lower = ROM_est - qt(0.975, 
                                  n_total - 1) * sqrt(ROM_var),
         CVR_upper = CVR_est + qt(0.975, 
                                  n_total - 1) * sqrt(CVR_var),
         CVR_lower = CVR_est - qt(0.975, 
                                  n_total - 1) * sqrt(CVR_var),
         VR_upper = VR_est + qt(0.975, 
                                n_total - 1) * sqrt(VR_var),
         VR_lower = VR_est - qt(0.975, 
                                n_total - 1) * sqrt(VR_var),
         SK_delta_upper = SK_delta_est + qt(0.975, 
                                            n_total - 1) * sqrt(SK_delta_var),
         SK_delta_lower = SK_delta_est - qt(0.975, 
                                            n_total - 1) * sqrt(SK_delta_var),
         KU_delta_upper = KU_delta_est + qt(0.975, 
                                            n_total - 1) * sqrt(KU_delta_var),
         KU_delta_lower = KU_delta_est - qt(0.975, 
                                            n_total - 1) * sqrt(KU_delta_var))

# getting n0
dat <- df_meta_analysed

dat$n0 <- dat$n_male*dat$n_female / (dat$n_male + dat$n_female)

vtilde <- 1/dat$n0
dat$effect_size_id <- 1:dim(dat)[1]
dat$sampling_error_id <- dat$effect_size_id

Vf <- diag(as.numeric(vtilde))
levs <- levels(factor(dat$sampling_error_id))
rownames(Vf) <- levs
colnames(Vf) <- levs

model <- 
  rma.mv(yi = KU_delta_est,
         V = 0,
         random = list(~ 1|effect_size_id,
                       # this is like sampling error bit which usually goes to V
                       ~ 1|sampling_error_id, 
                       ~ 1|phenotyping_center, 
                       ~ 1|strain_fig),
         data = dat,
         method = "REML",
         test = "t",
         R = list(sampling_error_id = Vf),
         Rscale = FALSE,
         # new optimizer
         control = list(rel.tol = 1e-8))

summary(model)

# how do I get I2 from this - 
# see this paper - https://onlinelibrary.wiley.com/doi/10.1111/ele.14144
sampling_variance <- model$sigma2[2]*vtilde

# from Higgins - https://doi.org/10.1002/sim.1186
av_sv <- sum(1 / sampling_variance) / 
  sum((1 / sampling_variance) ^ 2) * sum(1 / sampling_variance) ^ 2 / (sum(1 / sampling_variance) * (length(sampling_variance) - 1))

# I2
(model$sigma2[1] + model$sigma2[3] + model$sigma2[4])/(av_sv + model$sigma2[1] + model$sigma2[3] + model$sigma2[4])
