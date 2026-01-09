### ------------------------------------------------------------------------###
# Meta-Analysis Validation Study - IPD Context
# Description: Validates meta-analytic performance of Δsk, Δku, and ΔZr
# in the context of the analysis presented in supp_material.qmd
#
# Reviewer Request:
# "Comparing the meta-analyzed effect estimate against an individual-rodent-data
# meta-analysis would inform how well the meta-analyzed metrics match the real data.
# Additional simulations could inform the performance: How does a meta-analysis
# perform under the null? Under different known effect sizes, variances, and sample
# sizes? Under heterogeneity of effects? This could be assessed through simulations
# based on the raw, individual-rodent data by creating subsamples, changing the
# effect size, variance, and sample size, and introducing heterogeneity."
### ------------------------------------------------------------------------###

# Clean working environment ----
rm(list = ls())

# Load required packages ----
pacman::p_load(
  tidyverse, # Data manipulation and visualization
  metafor, # Meta-analysis (rma.mv)
  orchaRd, # For moment_effects function
  lme4, # For IPD meta-analysis
  boot, # Bootstrap methods
  moments, # Moment calculations
  PearsonDS, # Pearson distribution simulation
  patchwork, # Combining plots
  latex2exp, # LaTeX in plots
  here, # Path management
  janitor # Data cleaning
)

# Source functions from existing codebase
source(here("R", "func.R"))

### ------------------------------------------------------------------------###
## Helper Functions - Using Same Formulas as Simulation Scripts
### ------------------------------------------------------------------------###

#' @title Calculate Sample Skewness or its Sampling Variance
#' @description Same function as in sim_skew.R for consistency
#' @param x A numeric vector
#' @param output "est" for estimate or "var" for sampling variance
#' @return Skewness estimate or its sampling variance
calc.skewness <- function(x, output = "est") {
  n <- length(x)
  if (n < 3) {
    return(NA)
  }

  if (output == "est") { # Skewness estimate
    (sqrt(n * (n - 1)) / (n - 2)) *
      (((1 / n) * sum((x - mean(x))^3)) /
        (((1 / n) * sum((x - mean(x))^2))^(3 / 2)))
  } else if (output == "var") { # Skewness sampling variance
    (6 * n * (n - 1)) /
      ((n - 2) * (n + 1) * (n + 3))
  }
}

#' @title Calculate Sample Kurtosis or its Sampling Variance
#' @description Same function as in sim_kurt.R for consistency
#' @param x A numeric vector
#' @param output "est" for estimate or "var" for sampling variance
#' @return Kurtosis estimate or its sampling variance
calc.kurtosis <- function(x, output = "est") {
  n <- length(x)
  if (n < 4) {
    return(NA)
  }

  if (output == "est") { # Kurtosis estimate (returns EXCESS kurtosis)
    ((((n + 1) * n * (n - 1)) / ((n - 2) * (n - 3))) *
      (sum((x - mean(x))^4) / (sum((x - mean(x))^2)^2))) -
      (3 * ((n - 1)^2) / ((n - 2) * (n - 3)))
  } else if (output == "var") { # Kurtosis sampling variance
    (24 * n * ((n - 1)^2)) /
      ((n - 3) * (n - 2) * (n + 3) * (n + 5))
  }
}

### ------------------------------------------------------------------------###
## PART 1: Individual Participant Data (IPD) Meta-Analysis
## Compare summary-level meta-analysis (as in supp_material.qmd) to IPD MA
### ------------------------------------------------------------------------###

#' @title Load and Prepare IMPC Data
#' @description Loads the mice data and prepares it following supp_material.qmd
#' @return List with raw data and summary data
load_and_prepare_impc_data <- function() {
  # Load raw data - same as in supp_material.qmd
  df_raw <-
    read_csv(here("mice_data_sample.csv")) %>%
    mutate(
      phenotyping_center =
        ifelse(phenotyping_center == "MRC Harwell",
          "MRC H",
          phenotyping_center
        ),
      strain_fig = case_when(
        strain_accession_id == "MGI:2159965" ~ "N",
        strain_accession_id == "MGI:2683688" ~ "NCrl",
        strain_accession_id == "MGI:2164831" ~ "NTac",
        strain_accession_id == "MGI:3056279" ~ "NJ",
        strain_accession_id == "MGI:2160139" ~ "NJcl"
      )
    )

  # Wide format for orchaRd::moment_effects
  df_raw_wide <-
    df_raw %>%
    select(-strain_accession_id) %>%
    pivot_wider(
      id_cols = c(specimen_id, trait_name, phenotyping_center, strain_fig),
      names_from = sex,
      values_from = value
    )

  # Summary data with moment effects - following supp_material.qmd approach
  df_meta_analysed <-
    df_raw %>%
    group_by(sex, trait_name, phenotyping_center, strain_fig) %>%
    summarize(
      mean = mean(value, na.rm = TRUE),
      sd = sd(value, na.rm = TRUE),
      n = n(),
      .groups = "drop"
    ) %>%
    pivot_wider(
      id_cols = c(trait_name, phenotyping_center, strain_fig),
      names_from = sex,
      values_from = c(mean:n)
    ) %>%
    # Calculate moment effects using orchaRd
    left_join(
      df_raw_wide %>%
        group_by(trait_name, phenotyping_center, strain_fig) %>%
        reframe(orchaRd::moment_effects(
          x1 = na.omit(male),
          x2 = na.omit(female),
          type = "skew"
        ))
    ) %>%
    left_join(
      df_raw_wide %>%
        group_by(trait_name, phenotyping_center, strain_fig) %>%
        reframe(orchaRd::moment_effects(
          x1 = na.omit(male),
          x2 = na.omit(female),
          type = "kurt"
        ))
    ) %>%
    rename(
      SK_delta_est = d_skew,
      SK_delta_var = d_skew_v,
      KU_delta_est = d_kurt,
      KU_delta_var = d_kurt_v
    ) %>%
    mutate(
      n_total = n_female + n_male,
      SK_delta_upper = SK_delta_est + qt(0.975, n_total - 1) * sqrt(SK_delta_var),
      SK_delta_lower = SK_delta_est - qt(0.975, n_total - 1) * sqrt(SK_delta_var),
      KU_delta_upper = KU_delta_est + qt(0.975, n_total - 1) * sqrt(KU_delta_var),
      KU_delta_lower = KU_delta_est - qt(0.975, n_total - 1) * sqrt(KU_delta_var)
    )

  return(list(
    raw_data = df_raw,
    raw_wide = df_raw_wide,
    summary_data = df_meta_analysed
  ))
}

#' @title Conduct Summary-Level Meta-Analysis (as in supp_material.qmd)
#' @description Performs the multi-level meta-analysis exactly as in the supplement
#' @param summary_data The summary statistics data frame
#' @param trait The trait to analyze
#' @param metric The metric ("SK_delta" or "KU_delta")
#' @param centers Centers to include (default: main 5)
#' @return Meta-analysis results using rma.mv
conduct_summary_ma_impc <- function(summary_data,
                                    trait,
                                    metric = "SK_delta",
                                    centers = c(
                                      "CCP-IMG", "HMGU", "JAX",
                                      "MRC H", "TCP"
                                    )) {
  # Filter and prepare data - same as process.ind_effects in supp_material.qmd
  ind_effects <- summary_data %>%
    filter(
      trait_name == trait,
      phenotyping_center %in% centers
    ) %>%
    mutate(type = "individual") %>%
    select(
      phenotyping_center,
      strain_fig,
      n = n_total,
      est = all_of(paste0(metric, "_est")),
      var = all_of(paste0(metric, "_var"))
    ) %>%
    rowid_to_column("effect_size_id")

  if (nrow(ind_effects) == 0) {
    return(NULL)
  }

  # Multi-level meta-analysis as in supp_material.qmd
  model <- tryCatch(
    {
      rma.mv(
        data = ind_effects,
        yi = est,
        V = var,
        test = "t",
        random = list(
          ~ 1 | effect_size_id,
          ~ 1 | phenotyping_center,
          ~ 1 | strain_fig
        )
      )
    },
    error = function(e) {
      message("Summary MA failed: ", e$message)
      return(NULL)
    }
  )

  if (is.null(model)) {
    return(NULL)
  }

  return(list(
    estimate = as.numeric(model$beta),
    se = as.numeric(model$se),
    ci_lb = as.numeric(model$ci.lb),
    ci_ub = as.numeric(model$ci.ub),
    tau2 = sum(model$sigma2), # Total heterogeneity
    I2 = orchaRd::i2_ml(model),
    n_studies = nrow(ind_effects),
    method = "Summary MA (rma.mv)",
    metric = metric,
    trait = trait
  ))
}

#' @title Conduct IPD Meta-Analysis
#' @description TRUE IPD approach: analyzes individual rodent values directly
#' @details In true IPD meta-analysis, each row (y_ij) is an individual rodent's value.
#' For skewness/kurtosis, we analyze distributional moments using individual-level data:
#' - For skewness: analyze cubed standardized residuals (captures 3rd moment)
#' - For kurtosis: analyze 4th power of standardized residuals (captures 4th moment)
#' The lmer model accounts for clustering within centers and strains while
#' estimating sex differences in distributional shape from ALL individual rodents.
#' @param raw_data The individual rodent data (each row = one rodent)
#' @param trait The trait to analyze
#' @param metric Type of analysis ("skewness" or "kurtosis")
#' @param centers Centers to include
#' @return IPD meta-analysis results
conduct_ipd_ma_impc <- function(raw_data,
                                trait,
                                metric = "skewness",
                                centers = c(
                                  "CCP-IMG", "HMGU", "JAX",
                                  "MRC H", "TCP"
                                )) {
  # Filter data - each row is an individual rodent (y_ij = individual value)
  trait_data <- raw_data %>%
    filter(
      trait_name == trait,
      phenotyping_center %in% centers
    ) %>%
    drop_na(value, sex) %>%
    mutate(
      # Create study identifier (center/strain combination)
      study_id = paste0(phenotyping_center, "_", strain_fig)
    )

  if (nrow(trait_data) < 20) {
    return(NULL)
  }

  # TRUE IPD meta-analysis: analyze individual rodent values directly
  ipd_result <- tryCatch(
    {
      if (metric == "skewness") {
        # IPD approach for skewness: analyze cubed standardized residuals
        # Step 1: Standardize within sex-study groups
        trait_data <- trait_data %>%
          group_by(study_id, sex) %>%
          mutate(
            value_centered = value - mean(value),
            value_sd = sd(value),
            value_standardized = value_centered / value_sd,
            # Cubed standardized values capture skewness (3rd moment)
            value_cubed = value_standardized^3
          ) %>%
          ungroup()

        # Step 2: Fit mixed-effects model on cubed values
        # This analyzes sex differences in distributional skewness
        # using ALL individual rodents (y_ij = cubed standardized value)
        ipd_model <- lmer(
          value_cubed ~ sex + (1 | phenotyping_center / strain_fig),
          data = trait_data,
          REML = TRUE
        )

        # Extract sex difference in skewness
        sex_effect <- fixef(ipd_model)["sexmale"]
        sex_se <- sqrt(vcov(ipd_model)["sexmale", "sexmale"])

        # The coefficient represents difference in E[Z^3] between sexes
        # For standardized values, E[Z^3] IS the skewness
        # So sex_effect is already a difference in skewness (Δsk)
        # No further scaling needed - it's already in the right units
        
        list(
          estimate = sex_effect,
          se = sex_se,
          ci_lb = sex_effect - 1.96 * sex_se,
          ci_ub = sex_effect + 1.96 * sex_se,
          method = "IPD (individual rodent-level analysis)",
          metric = "skewness",
          model_type = "lmer on cubed standardized values",
          n_individuals = nrow(trait_data),
          n_studies = n_distinct(trait_data$study_id),
          y_ij_description = "cubed standardized value for each rodent"
        )
      } else if (metric == "kurtosis") {
        # IPD approach for kurtosis: analyze 4th power of standardized residuals
        # Step 1: Standardize within sex-study groups
        trait_data <- trait_data %>%
          group_by(study_id, sex) %>%
          mutate(
            value_centered = value - mean(value),
            value_sd = sd(value),
            value_standardized = value_centered / value_sd,
            # 4th power captures kurtosis (4th moment)
            value_fourth = value_standardized^4
          ) %>%
          ungroup()

        # Step 2: Fit mixed-effects model on 4th power values
        # This analyzes sex differences in distributional kurtosis
        # using ALL individual rodents (y_ij = 4th power of standardized value)
        ipd_model <- lmer(
          value_fourth ~ sex + (1 | phenotyping_center / strain_fig),
          data = trait_data,
          REML = TRUE
        )

        # Extract sex difference
        sex_effect <- fixef(ipd_model)["sexmale"]
        sex_se <- sqrt(vcov(ipd_model)["sexmale", "sexmale"])

        # Under normality, E[Z^4] = 3 (raw 4th moment)
        # The coefficient represents difference in E[Z^4] between sexes
        # For excess kurtosis, we need the difference in (E[Z^4] - 3)
        # Since we're looking at differences, the -3 cancels out
        # So sex_effect is already Δ(excess kurtosis)
        
        list(
          estimate = sex_effect,
          se = sex_se,
          ci_lb = sex_effect - 1.96 * sex_se,
          ci_ub = sex_effect + 1.96 * sex_se,
          method = "IPD (individual rodent-level analysis)",
          metric = "kurtosis",
          model_type = "lmer on 4th power standardized values",
          n_individuals = nrow(trait_data),
          n_studies = n_distinct(trait_data$study_id),
          y_ij_description = "4th power of standardized value for each rodent"
        )
      }
    },
    error = function(e) {
      message("IPD MA failed: ", e$message)
      return(NULL)
    }
  )

  return(ipd_result)
}

### ------------------------------------------------------------------------###
## KEY DISTINCTION: Summary-level MA vs TRUE IPD MA
### ------------------------------------------------------------------------###
##
## SUMMARY-LEVEL MA (as in supp_material.qmd):
## - Input: Summary statistics per study (mean, SD, skewness, kurtosis)
## - Process: Calculate Δsk or Δku for each study, then meta-analyze
## - Unit of analysis: Study-level effect sizes (one per center/strain)
## - Meta-analysis: rma.mv on study-level Δsk or Δku values
##
## TRUE IPD MA (conduct_ipd_ma_impc function above):
## - Input: Individual rodent values (y_ij for each rodent i in study j)
## - Process: Transform each individual value (cube or 4th power),
##            then fit mixed model on ALL individual rodents
## - Unit of analysis: Individual rodent (thousands of data points)
## - Model: lmer on transformed individual values
##          y*_ij = f(y_ij) where f = cubing or 4th power
## - Benefits:
##   1. Uses all individual-level information
##   2. Accounts for within-study correlation structure
##   3. More precise estimates with correct uncertainty
##   4. Gold standard for assessing summary-level MA validity
##
### ------------------------------------------------------------------------###

#' @title Compare Summary MA vs IPD MA
#' @description Compares summary-level MA (as done in supplement) vs true IPD MA
#' @details This comparison assesses how well the summary-level approach
#' (calculating Δsk/Δku per study then meta-analyzing) matches the gold-standard
#' IPD approach (analyzing all individual rodent values in one model).
#' Differences indicate potential information loss or bias in summary approach.
compare_summary_vs_ipd <- function() {
  message("\n=== Loading IMPC Data ===\n")
  data_list <- load_and_prepare_impc_data()

  # Traits to analyze (from supp_material.qmd)
  traits <- c("fat_mass", "heart_weight", "glucose", "total_cholesterol")
  metrics_summary <- c("SK_delta", "KU_delta")
  metrics_ipd <- c("skewness", "kurtosis")

  comparison_results <- list()
  idx <- 1

  for (i in seq_along(traits)) {
    for (j in seq_along(metrics_summary)) {
      trait <- traits[i]
      metric_sum <- metrics_summary[j]
      metric_ipd <- metrics_ipd[j]

      message("Analyzing: ", trait, " - ", metric_sum)

      # Summary-level MA (as in supp_material.qmd)
      summary_ma <- conduct_summary_ma_impc(
        data_list$summary_data,
        trait = trait,
        metric = metric_sum
      )

      # IPD MA (gold standard)
      ipd_ma <- conduct_ipd_ma_impc(
        data_list$raw_data,
        trait = trait,
        metric = metric_ipd
      )

      if (!is.null(summary_ma) && !is.null(ipd_ma)) {
        comparison_results[[idx]] <- tibble(
          trait = trait,
          metric = metric_sum,
          summary_estimate = summary_ma$estimate,
          summary_se = summary_ma$se,
          summary_ci_lb = summary_ma$ci_lb,
          summary_ci_ub = summary_ma$ci_ub,
          ipd_estimate = ipd_ma$estimate,
          ipd_se = ipd_ma$se,
          ipd_ci_lb = ipd_ma$ci_lb,
          ipd_ci_ub = ipd_ma$ci_ub,
          difference = summary_ma$estimate - ipd_ma$estimate,
          # Avoid division by zero - only calculate if IPD estimate is non-zero
          relative_difference = ifelse(abs(ipd_ma$estimate) < 1e-10,
            NA_real_,
            (summary_ma$estimate - ipd_ma$estimate) / abs(ipd_ma$estimate) * 100
          ),
          overlap = (summary_ma$ci_ub >= ipd_ma$ci_lb) &
            (summary_ma$ci_lb <= ipd_ma$ci_ub)
        )
        idx <- idx + 1
      }
    }
  }

  return(bind_rows(comparison_results))
}

### ------------------------------------------------------------------------###
## PART 2: Simulation-Based Validation Using Real IMPC Data
## Create subsamples, manipulate effects, test MA performance
### ------------------------------------------------------------------------###

#' @title Subsample IMPC Data for Simulations
#' @description Creates subsamples from real data with controlled properties
#' @details Uses sampling WITH REPLACEMENT to enable any requested sample size.
#' As long as there's some data for each sex, we can create studies of any size
#' by resampling. This allows testing different sample size scenarios even when
#' the original dataset is limited.
#' @param raw_data The individual rodent data
#' @param trait Trait to subsample
#' @param n_studies Number of "studies" to create
#' @param sample_size_per_study Sample size per study (per sex)
#' @param effect_manipulation How to manipulate effect ("null", "small", "medium", "large")
#' @param heterogeneity Level of heterogeneity to introduce (0 = none, >0 = tau^2)
#' @return Subsampled data ready for MA
subsample_impc_data <- function(raw_data,
                                trait,
                                n_studies = 10,
                                sample_size_per_study = 50,
                                effect_manipulation = "null",
                                heterogeneity = 0,
                                metric = "skewness") {
  # Get all data for trait
  trait_data <- raw_data %>%
    filter(trait_name == trait) %>%
    drop_na(value, sex)

  # Check we have at least SOME data for each sex
  # (with replacement, we can create any size sample from even 1 observation)
  n_male <- sum(trait_data$sex == "male")
  n_female <- sum(trait_data$sex == "female")
  
  min_required <- 10  # Minimum for reasonable distribution estimation
  
  if (n_male < min_required || n_female < min_required) {
    message("Insufficient data for trait: ", trait, 
            " (male: ", n_male, ", female: ", n_female, ")")
    return(NULL)
  }
  
  # Information message if heavily resampling
  total_requested <- n_studies * sample_size_per_study
  if (total_requested > n_male || total_requested > n_female) {
    message("Note: Resampling with replacement for trait '", trait, 
            "' (requested: ", total_requested, " per sex, ",
            "available: male=", n_male, ", female=", n_female, ")")
  }

  # Generate heterogeneous true effects if needed
  if (heterogeneity > 0) {
    if (effect_manipulation == "null") {
      true_effects <- rnorm(n_studies, mean = 0, sd = sqrt(heterogeneity))
    } else if (effect_manipulation == "small") {
      true_effects <- rnorm(n_studies, mean = 0.3, sd = sqrt(heterogeneity))
    } else if (effect_manipulation == "medium") {
      true_effects <- rnorm(n_studies, mean = 0.5, sd = sqrt(heterogeneity))
    } else if (effect_manipulation == "large") {
      true_effects <- rnorm(n_studies, mean = 1.0, sd = sqrt(heterogeneity))
    }
  } else {
    # No heterogeneity - same effect for all
    if (effect_manipulation == "null") {
      true_effects <- rep(0, n_studies)
    } else if (effect_manipulation == "small") {
      true_effects <- rep(0.3, n_studies)
    } else if (effect_manipulation == "medium") {
      true_effects <- rep(0.5, n_studies)
    } else if (effect_manipulation == "large") {
      true_effects <- rep(1.0, n_studies)
    }
  }

  # Create subsamples
  study_data <- list()

  for (study in 1:n_studies) {
    
    # Safely sample with replacement from real data
    # This works even if requested sample size > available data
    tryCatch({
      male_sample <- trait_data %>%
        filter(sex == "male") %>%
        slice_sample(n = sample_size_per_study, replace = TRUE)

      female_sample <- trait_data %>%
        filter(sex == "female") %>%
        slice_sample(n = sample_size_per_study, replace = TRUE)
      
      # Verify we got the requested sample size
      if (nrow(male_sample) != sample_size_per_study || 
          nrow(female_sample) != sample_size_per_study) {
        warning("Study ", study, ": Got unexpected sample size")
        next
      }

      # Apply effect manipulation
      if (metric == "skewness" && effect_manipulation != "null") {
        # Manipulate distribution to create skewness difference
        # Apply transformation to male data
        male_values <- male_sample$value
        # Box-Cox-like transformation to induce skewness
        lambda <- 1 + true_effects[study] * 0.5
        if (lambda != 0 && all(is.finite(male_values))) {
          # Ensure transformation is valid (positive values for power transform)
          if (any(male_values <= 0)) {
            # Shift to positive range
            male_values_shifted <- male_values - min(male_values) + 1
            male_sample$value <- sign(male_values_shifted) *
              abs(male_values_shifted)^lambda
          } else {
            male_sample$value <- male_values^lambda
          }
        }
      } else if (metric == "kurtosis" && effect_manipulation != "null") {
        # Manipulate to create kurtosis difference
        # Add outliers proportional to effect size
        n_outliers <- round(sample_size_per_study * 0.1 * abs(true_effects[study]))
        if (n_outliers > 0 && n_outliers < sample_size_per_study) {
          outlier_indices <- sample(1:sample_size_per_study, n_outliers)
          male_sample$value[outlier_indices] <-
            male_sample$value[outlier_indices] *
              (1 + 3 * sign(true_effects[study]))
        }
      }

      # Calculate moment effect using orchaRd (as in supp_material.qmd)
      moment_result <- orchaRd::moment_effects(
        x1 = male_sample$value,
        x2 = female_sample$value,
        type = ifelse(metric == "skewness", "skew", "kurt")
      )

      study_data[[study]] <- tibble(
        study_id = study,
        true_effect = true_effects[study],
        est = ifelse(metric == "skewness",
          moment_result$d_skew,
          moment_result$d_kurt
        ),
        var = ifelse(metric == "skewness",
          moment_result$d_skew_v,
          moment_result$d_kurt_v
        ),
        n_male = sample_size_per_study,
        n_female = sample_size_per_study
      )
      
    }, error = function(e) {
      warning("Study ", study, " failed: ", e$message)
      # Return NULL for this study - will be filtered out
      return(NULL)
    })
  }

  # Filter out any NULL entries (failed studies) and bind
  study_data <- study_data[!sapply(study_data, is.null)]
  
  if (length(study_data) == 0) {
    message("All studies failed for trait: ", trait)
    return(NULL)
  }
  
  if (length(study_data) < n_studies) {
    message("Warning: Only ", length(study_data), "/", n_studies, 
            " studies succeeded")
  }

  return(bind_rows(study_data))
}

#' @title Test MA Performance Under Null
#' @description Tests Type I error rate
#' @param nsims Number of simulation iterations
#' @param n_studies Number of studies per MA
#' @param sample_size Sample size per study
#' @param trait Trait to use for subsampling
#' @param metric Metric to test
#' @return Performance metrics
test_ma_under_null <- function(nsims = 100,
                               n_studies = 10,
                               sample_size = 50,
                               trait = "glucose",
                               metric = "skewness") {
  message("\nTesting MA under null hypothesis...")
  message("Trait: ", trait, ", Metric: ", metric)

  # Load data once
  data_list <- load_and_prepare_impc_data()

  # Storage
  reject_h0 <- numeric(nsims)
  estimates <- numeric(nsims)
  ses <- numeric(nsims)

  for (sim in 1:nsims) {
    if (sim %% 20 == 0) message("  Simulation ", sim, "/", nsims)

    # Subsample with null effect
    study_data <- subsample_impc_data(
      data_list$raw_data,
      trait = trait,
      n_studies = n_studies,
      sample_size_per_study = sample_size,
      effect_manipulation = "null",
      heterogeneity = 0,
      metric = metric
    )

    if (is.null(study_data)) next

    # Meta-analysis using rma.mv (as in supp_material.qmd)
    ma <- tryCatch(
      {
        rma.mv(
          data = study_data,
          yi = est,
          V = var,
          test = "t",
          random = ~ 1 | study_id
        )
      },
      error = function(e) NULL
    )

    if (!is.null(ma)) {
      estimates[sim] <- as.numeric(ma$beta)
      ses[sim] <- as.numeric(ma$se)
      # Test H0: estimate = 0
      # For rma.mv with test="t", use the p-value directly
      reject_h0[sim] <- ma$pval < 0.05
    }
  }

  # Calculate performance
  n_successful <- sum(!is.na(estimates))
  
  if (n_successful < nsims * 0.5) {
    warning("Only ", n_successful, "/", nsims, " simulations succeeded. ",
            "Results may be unreliable.")
  }
  
  return(tibble(
    metric = metric,
    trait = trait,
    n_studies = n_studies,
    sample_size = sample_size,
    nsims = nsims,
    n_successful = n_successful,
    type1_error = mean(reject_h0, na.rm = TRUE),
    mean_estimate = mean(estimates, na.rm = TRUE),
    sd_estimates = sd(estimates, na.rm = TRUE),
    mean_se = mean(ses, na.rm = TRUE)
  ))
}

#' @title Test MA Performance Under Known Effects
#' @description Tests power and bias under various effect sizes
test_ma_known_effects <- function(nsims = 100,
                                  n_studies = 10,
                                  sample_size = 50,
                                  trait = "glucose",
                                  metric = "skewness",
                                  effect_size = "medium") {
  message("\nTesting MA under known effect: ", effect_size)
  message("Trait: ", trait, ", Metric: ", metric)

  # Load data
  data_list <- load_and_prepare_impc_data()

  # Storage
  reject_h0 <- numeric(nsims)
  estimates <- numeric(nsims)
  ses <- numeric(nsims)
  coverage <- numeric(nsims)

  # True effect value
  true_effect <- switch(effect_size,
    "small" = 0.3,
    "medium" = 0.5,
    "large" = 1.0
  )

  for (sim in 1:nsims) {
    if (sim %% 20 == 0) message("  Simulation ", sim, "/", nsims)

    # Subsample with known effect
    study_data <- subsample_impc_data(
      data_list$raw_data,
      trait = trait,
      n_studies = n_studies,
      sample_size_per_study = sample_size,
      effect_manipulation = effect_size,
      heterogeneity = 0,
      metric = metric
    )

    if (is.null(study_data)) next

    # Meta-analysis
    ma <- tryCatch(
      {
        rma.mv(
          data = study_data,
          yi = est,
          V = var,
          test = "t",
          random = ~ 1 | study_id
        )
      },
      error = function(e) NULL
    )

    if (!is.null(ma)) {
      estimates[sim] <- as.numeric(ma$beta)
      ses[sim] <- as.numeric(ma$se)

      # Use p-value directly from model
      reject_h0[sim] <- ma$pval < 0.05

      # Coverage
      ci_lb <- as.numeric(ma$ci.lb)
      ci_ub <- as.numeric(ma$ci.ub)
      coverage[sim] <- (true_effect >= ci_lb) & (true_effect <= ci_ub)
    }
  }

  # Performance metrics
  n_successful <- sum(!is.na(estimates))
  
  if (n_successful < nsims * 0.5) {
    warning("Only ", n_successful, "/", nsims, " simulations succeeded. ",
            "Results may be unreliable.")
  }
  
  return(tibble(
    metric = metric,
    trait = trait,
    n_studies = n_studies,
    sample_size = sample_size,
    true_effect = true_effect,
    nsims = nsims,
    n_successful = n_successful,
    mean_estimate = mean(estimates, na.rm = TRUE),
    bias = mean(estimates, na.rm = TRUE) - true_effect,
    rmse = sqrt(mean((estimates - true_effect)^2, na.rm = TRUE)),
    power = mean(reject_h0, na.rm = TRUE),
    coverage = mean(coverage, na.rm = TRUE)
  ))
}

#' @title Test MA Performance Under Heterogeneity
#' @description Tests how MA handles between-study variability
test_ma_heterogeneity <- function(nsims = 100,
                                  n_studies = 10,
                                  sample_size = 50,
                                  trait = "glucose",
                                  metric = "skewness",
                                  tau2 = 0.1,
                                  mean_effect = 0.5) {
  message("\nTesting MA under heterogeneity (tau2 = ", tau2, ")")
  message("Trait: ", trait, ", Metric: ", metric)

  # Load data
  data_list <- load_and_prepare_impc_data()

  # Storage
  estimates <- numeric(nsims)
  tau2_estimates <- numeric(nsims)
  coverage <- numeric(nsims)

  for (sim in 1:nsims) {
    if (sim %% 20 == 0) message("  Simulation ", sim, "/", nsims)

    # Subsample with heterogeneity
    study_data <- subsample_impc_data(
      data_list$raw_data,
      trait = trait,
      n_studies = n_studies,
      sample_size_per_study = sample_size,
      effect_manipulation = ifelse(mean_effect > 0, "medium", "null"),
      heterogeneity = tau2,
      metric = metric
    )

    if (is.null(study_data)) next

    # Meta-analysis
    ma <- tryCatch(
      {
        rma.mv(
          data = study_data,
          yi = est,
          V = var,
          test = "t",
          random = ~ 1 | study_id
        )
      },
      error = function(e) NULL
    )

    if (!is.null(ma)) {
      estimates[sim] <- as.numeric(ma$beta)
      tau2_estimates[sim] <- sum(ma$sigma2)

      # Coverage for mean effect
      ci_lb <- ma$ci.lb
      ci_ub <- ma$ci.ub
      coverage[sim] <- (mean_effect >= ci_lb) & (mean_effect <= ci_ub)
    }
  }

  # Performance
  n_successful <- sum(!is.na(estimates))
  
  if (n_successful < nsims * 0.5) {
    warning("Only ", n_successful, "/", nsims, " simulations succeeded. ",
            "Results may be unreliable.")
  }
  
  return(tibble(
    metric = metric,
    trait = trait,
    n_studies = n_studies,
    sample_size = sample_size,
    true_tau2 = tau2,
    mean_effect = mean_effect,
    nsims = nsims,
    n_successful = n_successful,
    mean_estimate = mean(estimates, na.rm = TRUE),
    bias = mean(estimates, na.rm = TRUE) - mean_effect,
    mean_tau2 = mean(tau2_estimates, na.rm = TRUE),
    tau2_bias = mean(tau2_estimates, na.rm = TRUE) - tau2,
    coverage = mean(coverage, na.rm = TRUE)
  ))
}

### ------------------------------------------------------------------------###
## PART 3: Comprehensive Simulation Grid
### ------------------------------------------------------------------------###

#' @title Run Complete Validation Study
#' @description Runs all validation analyses
run_complete_validation <- function(nsims = 100, save_results = TRUE, n_studies = 100) {
  message("\n", rep("=", 70))
  message("META-ANALYSIS VALIDATION STUDY - IMPC CONTEXT")
  message(rep("=", 70), "\n")

  # Part 1: IPD vs Summary comparison
  message("\n### PART 1: Summary MA vs IPD MA Comparison ###\n")
  ipd_comparison <- tryCatch(
    {
      compare_summary_vs_ipd()
    },
    error = function(e) {
      message("IPD comparison failed: ", e$message)
      NULL
    }
  )

  if (!is.null(ipd_comparison)) {
    print(ipd_comparison)
  }

  # Part 2: Simulations
  message("\n### PART 2: Simulation-Based Validation ###\n")

  traits <- c("glucose", "fat_mass")
  metrics <- c("skewness", "kurtosis")
  sample_sizes <- c(50, 80, 200)
  effect_sizes <- c("null", "small", "medium", "large")
  tau2_values <- c(0, 0.05, 0.1, 0.25)

  # 2a. Type I error
  message("\n2a. Testing Type I Error Rates...\n")
  null_results <- list()
  idx <- 1
  for (trait in traits) {
    for (metric in metrics) {
      for (ss in sample_sizes) {
        null_results[[idx]] <- test_ma_under_null(
          nsims = nsims,
          n_studies = n_studies,
          sample_size = ss,
          trait = trait,
          metric = metric
        )
        idx <- idx + 1
      }
    }
  }
  null_results_df <- bind_rows(null_results)
  print(null_results_df)

  # 2b. Power and bias
  message("\n2b. Testing Power and Bias Under Known Effects...\n")
  power_results <- list()
  idx <- 1
  for (trait in traits[1]) { # Just first trait for time
    for (metric in metrics) {
      for (effect in effect_sizes[-1]) { # Skip null
        for (ss in sample_sizes) {
          power_results[[idx]] <- test_ma_known_effects(
            nsims = nsims,
            n_studies = n_studies,
            sample_size = ss,
            trait = trait,
            metric = metric,
            effect_size = effect
          )
          idx <- idx + 1
        }
      }
    }
  }
  power_results_df <- bind_rows(power_results)
  print(power_results_df)

  # 2c. Heterogeneity
  message("\n2c. Testing Performance Under Heterogeneity...\n")
  het_results <- list()
  idx <- 1
  for (trait in traits[1]) {
    for (metric in metrics) {
      for (tau2 in tau2_values) {
        het_results[[idx]] <- test_ma_heterogeneity(
          nsims = nsims,
          n_studies = n_studies,
          sample_size = 50,
          trait = trait,
          metric = metric,
          tau2 = tau2,
          mean_effect = 0.5
        )
        idx <- idx + 1
      }
    }
  }
  het_results_df <- bind_rows(het_results)
  print(het_results_df)

  # Save results
  if (save_results) {
    message("\nSaving results...")

    if (!dir.exists(here("output", "ma_validation"))) {
      dir.create(here("output", "ma_validation"), recursive = TRUE)
    }

    saveRDS(
      ipd_comparison,
      here("output", "ma_validation", "ipd_comparison.rds")
    )
    saveRDS(
      null_results_df,
      here("output", "ma_validation", "type1_error_results.rds")
    )
    saveRDS(
      power_results_df,
      here("output", "ma_validation", "power_bias_results.rds")
    )
    saveRDS(
      het_results_df,
      here("output", "ma_validation", "heterogeneity_results.rds")
    )

    # Summary CSV files
    write_csv(
      ipd_comparison,
      here("output", "ma_validation", "ipd_comparison.csv")
    )
    write_csv(
      null_results_df,
      here("output", "ma_validation", "type1_error_results.csv")
    )
    write_csv(
      power_results_df,
      here("output", "ma_validation", "power_bias_results.csv")
    )
    write_csv(
      het_results_df,
      here("output", "ma_validation", "heterogeneity_results.csv")
    )

    message("Results saved to ./output/ma_validation/")
  }

  message("\n", rep("=", 70))
  message("VALIDATION COMPLETE!")
  message(rep("=", 70), "\n")

  return(list(
    ipd_comparison = ipd_comparison,
    type1_error = null_results_df,
    power_bias = power_results_df,
    heterogeneity = het_results_df
  ))
}

### ------------------------------------------------------------------------###
##
## # Run complete validation (may take 15-30 minutes with nsims=200)
results <- run_complete_validation(
  nsims = 100,
  save_results = FALSE, n_studies = 20
)
results2 <- run_complete_validation(
  nsims = 200,
  save_results = FALSE,
  n_studies = 100
)
##
## # Run specific tests:
## # 1. IPD vs Summary comparison
ipd_comp <- compare_summary_vs_ipd()
##
## # 2. Type I error test
null_test <- test_ma_under_null(
  nsims = 200, trait = "glucose",
  metric = "skewness"
)
##
## # 3. Power test
power_test <- test_ma_known_effects(
  nsims = 100, trait = "glucose",
  metric = "skewness",
  effect_size = "medium"
)
##
## # 4. Heterogeneity test
het_test <- test_ma_heterogeneity(
  nsims = 100, trait = "glucose",
  metric = "skewness", tau2 = 0.1
)
### ------------------------------------------------------------------------###
