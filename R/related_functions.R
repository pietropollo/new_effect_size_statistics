## ================================================================
## 1. delta-method plug-in variance for sample skewness  (g1)
##    – returns a list:  estimate, var, se
## ================================================================
#' @title Delta-method Variance for Sample Skewness
#' @description Computes the sample skewness and estimates its variance using the delta method.
#' @param x A numeric vector.
#' @return A list with the skewness estimate (`est`), its variance (`var`), and standard error (`se`).
#' @examples
#' x <- rgamma(25, shape = 5)
#' skew_delta(x)
skew_delta <- function(x) {
  n  <- length(x)
  xc <- x - mean(x)
  
  ## central moments up to order 6
  m2 <- mean(xc^2); m3 <- mean(xc^3); m4 <- mean(xc^4)
  m5 <- mean(xc^5); m6 <- mean(xc^6)
  
  ## standardised moments (λr = μr / σ^r)
  l3 <- m3 / m2^(3/2)
  l4 <- m4 / m2^2
  l5 <- m5 / m2^(5/2)
  l6 <- m6 / m2^3
  
  ## Δ-method variance  (Eq. 1 in the previous reply)
  var_g1 <- ( l6
              - 3 * l3 * l5
              + 2.25 * l3^2 * (l4 - 1)
              - 0.25 * l3^2 ) / n
  
  list(est = l3,
       var = var_g1,
       se  = sqrt(var_g1))
}


## ================================================================
## 2. delta-method plug-in variance for sample *excess* kurtosis  (g2)
##    – returns a list:  estimate, var, se
## ================================================================
#' @title Delta-method Variance for Sample Excess Kurtosis
#' @description Computes the sample excess kurtosis and estimates its variance using the delta method.
#' @param x A numeric vector.
#' @return A list with the excess kurtosis estimate (`est`), its variance (`var`), and standard error (`se`).
#' @examples
#' x <- rgamma(25, shape = 5)
#' kurt_delta(x)
kurt_delta <- function(x) {
  n  <- length(x)
  xc <- x - mean(x)
  
  ## central moments up to order 8
  m2 <- mean(xc^2); m4 <- mean(xc^4)
  m6 <- mean(xc^6); m8 <- mean(xc^8)
  
  l4 <- m4 / m2^2        # raw kurtosis
  l6 <- m6 / m2^3
  l8 <- m8 / m2^4
  
  ## Δ-method variance  (Eq. 2 in the previous reply)
  var_g2 <- ( l8
              - 4 * l4 * l6
              + 4 * l4^3
              - l4^2 ) / n
  
  g2 <- l4 - 3            # convert to *excess* kurtosis
  
  list(est = g2,
       var = var_g2,
       se  = sqrt(var_g2))
}


## ================================================================
## 3. Non-parametric bootstrap for skewness
##    – returns a list:  point est, bias-corrected est, se, replicates*
## ================================================================
#' @title Bootstrap Estimation of Skewness
#' @description Computes bootstrap estimate and standard error of skewness.
#' @param x A numeric vector.
#' @param B Number of bootstrap replicates. Default is 2000.
#' @param bias.correct Logical, whether to apply bias correction. Default is TRUE.
#' @param return.replicates Logical, whether to return replicate values. Default is FALSE.
#' @return A list with estimate (`est`), bias-corrected estimate (`est_bc`), variance, and standard error.
#' @examples
#' x <- rgamma(25, shape = 5)
#' boot_skew(x)
boot_skew <- function(x, B = 2000, bias.correct = TRUE,
                      return.replicates = FALSE) {
  g1 <- function(z) mean((z - mean(z))^3) /
    mean((z - mean(z))^2)^(3/2)
  
  b.reps <- replicate(B, g1(sample(x, replace = TRUE)))
  
  est    <- g1(x)
  est.bc <- if (bias.correct) 2 * est - mean(b.reps) else est
  out    <- list(est     = est,
                 est_bc  = est.bc,
                 var     = sd(b.reps)^2,
                 se      = sd(b.reps))
  if (return.replicates) out$boot <- b.reps
  out
}


## ================================================================
## 4. Non-parametric bootstrap for excess kurtosis
##    – returns a list:  point est, bias-corrected est, se, replicates*
## ================================================================
#' @title Bootstrap Estimation of Excess Kurtosis
#' @description Computes bootstrap estimate and standard error of excess kurtosis.
#' @param x A numeric vector.
#' @param B Number of bootstrap replicates. Default is 2000.
#' @param bias.correct Logical, whether to apply bias correction. Default is TRUE.
#' @param return.replicates Logical, whether to return replicate values. Default is FALSE.
#' @return A list with estimate (`est`), bias-corrected estimate (`est_bc`), variance, and standard error.
#' @examples
#' x <- rgamma(25, shape = 5)
#' boot_kurt(x)
boot_kurt <- function(x, B = 2000, bias.correct = TRUE,
                      return.replicates = FALSE) {
  g2 <- function(z) mean((z - mean(z))^4) /
    mean((z - mean(z))^2)^2 - 3
  
  b.reps <- replicate(B, g2(sample(x, replace = TRUE)))
  
  est    <- g2(x)
  est.bc <- if (bias.correct) 2 * est - mean(b.reps) else est
  out    <- list(est     = est,
                 est_bc  = est.bc,
                 var     = sd(b.reps)^2,
                 se      = sd(b.reps))
  if (return.replicates) out$boot <- b.reps
  out
}


## ================================================================
## 5. Jack-knife SE (and bias-correction) for skewness
## ================================================================
#' @title Jackknife Estimation of Skewness
#' @description Computes jackknife estimate and standard error for skewness.
#' @param x A numeric vector.
#' @param bias.correct Logical, whether to apply bias correction. Default is TRUE.
#' @param return.replicates Logical, whether to return jackknife replicates. Default is FALSE.
#' @return A list with estimate (`est`), bias-corrected estimate (`est_bc`), variance, and standard error.
#' @examples
#' x <- rgamma(25, shape = 5)
#' jack_skew(x)
jack_skew <- function(x, bias.correct = TRUE,
                      return.replicates = FALSE) {
  
  n   <- length(x)
  g1  <- function(z) mean((z - mean(z))^3) /
    mean((z - mean(z))^2)^(3/2)
  
  ## leave-one-out replicates
  theta_i <- vapply(seq_len(n),
                    function(i) g1(x[-i]),
                    numeric(1))
  
  theta   <- g1(x)                 # full-sample estimate
  theta_bar <- mean(theta_i)
  
  ## jack-knife standard error
  se_jack <- sqrt((n - 1) * mean((theta_i - theta_bar)^2))
  
  ## bias correction
  theta_bc <- if (bias.correct) n * theta - (n - 1) * theta_bar else theta
  
  out <- list(est    = theta,
              est_bc = theta_bc,
              var    = se_jack^2,
              se     = se_jack)
  
  if (return.replicates) out$jack <- theta_i
  out
}

## ================================================================
## 6. Jack-knife SE (and bias-correction) for *excess* kurtosis
## ================================================================
#' @title Jackknife Estimation of Excess Kurtosis
#' @description Computes jackknife estimate and standard error for excess kurtosis.
#' @param x A numeric vector.
#' @param bias.correct Logical, whether to apply bias correction. Default is TRUE.
#' @param return.replicates Logical, whether to return jackknife replicates. Default is FALSE.
#' @return A list with estimate (`est`), bias-corrected estimate (`est_bc`), variance, and standard error.
#' @examples
#' x <- rgamma(25, shape = 5)
#' jack_kurt(x)
jack_kurt <- function(x, bias.correct = TRUE,
                      return.replicates = FALSE) {
  
  n   <- length(x)
  g2  <- function(z) mean((z - mean(z))^4) /
    mean((z - mean(z))^2)^2 - 3
  
  theta_i <- vapply(seq_len(n),
                    function(i) g2(x[-i]),
                    numeric(1))
  
  theta     <- g2(x)
  theta_bar <- mean(theta_i)
  
  se_jack   <- sqrt((n - 1) * mean((theta_i - theta_bar)^2))
  
  theta_bc  <- if (bias.correct) n * theta - (n - 1) * theta_bar else theta
  
  out <- list(est    = theta,
              est_bc = theta_bc,
              var    = se_jack^2,
              se     = se_jack)
  
  if (return.replicates) out$jack <- theta_i
  out
}

## ================================================================
## 7. Jack-knife SE (and bias-correction) for correlation
## ================================================================
#' @title Jackknife Estimation of Correlation
#' @description Computes jackknife estimate and standard error for the Pearson correlation between two variables.
#' @param dat A two-column data frame or matrix.
#' @param bias.correct Logical, whether to apply bias correction. Default is TRUE.
#' @param return.replicates Logical, whether to return jackknife replicates. Default is FALSE.
#' @return A list with estimate (`est`), bias-corrected estimate (`est_bc`), variance, and standard error.
#' @examples
#' g2 <- MASS::mvrnorm(1000, mu = c(0, 0), Sigma = matrix(c(1, -0.4, -0.4, 1), 2))
#' jack_cor(g2)
jack_cor <- function(dat, bias.correct = TRUE,
                      return.replicates = FALSE) {
  
  n   <- nrow(dat)
  g2  <- function(z) cor(z)[1,2]
  
  theta_i <- vapply(seq_len(n),
                    function(i) g2(dat[-i, ]),
                    numeric(1))
  
  theta     <- g2(dat)
  theta_bar <- mean(theta_i)
  
  se_jack   <- sqrt((n - 1) * mean((theta_i - theta_bar)^2))
  
  theta_bc  <- if (bias.correct) n * theta - (n - 1) * theta_bar else theta
  
  out <- list(est    = theta,
              est_bc = theta_bc,
              var    = se_jack^2,
              se     = se_jack)
  
  if (return.replicates) out$jack <- theta_i
  out
}

## ================================================================
## 8. Non-parametric bootstrap for correlation
##    – returns a list:  point est, bias-corrected est, se, replicates*
## ================================================================
#' @title Bootstrap Estimation of Correlation
#' @description Computes bootstrap estimate and standard error for the Pearson correlation between two variables.
#' @param dat A two-column data frame or matrix.
#' @param B Number of bootstrap replicates. Default is 2000.
#' @param bias.correct Logical, whether to apply bias correction. Default is TRUE.
#' @param return.replicates Logical, whether to return replicate values. Default is FALSE.
#' @return A list with estimate (`est`), bias-corrected estimate (`est_bc`), variance, and standard error.
#' @examples
#' g2 <- MASS::mvrnorm(1000, mu = c(0, 0), Sigma = matrix(c(1, -0.4, -0.4, 1), 2))
#' boot_cor(g2)
boot_cor <- function(dat, B = 2000, bias.correct = TRUE,
                     return.replicates = FALSE) {
  g2 <- function(z) cor(z)[1,2]

  b.reps <- replicate(B, {
    sample_rows <- dat[sample(nrow(dat), replace = TRUE), ]
    g2(sample_rows)
  })

  est    <- g2(dat)
  est.bc <- if (bias.correct) 2 * est - mean(b.reps) else est
  out    <- list(est     = est,
                 est_bc  = est.bc,
                 var     = sd(b.reps)^2,
                 se      = sd(b.reps))
  if (return.replicates) out$boot <- b.reps
  out
}

# How to use it

set.seed(1)
x <- rgamma(25, shape = 5)  # moderately skewed, kurtotic data

sk_res  <- skew_delta(x)
ku_res  <- kurt_delta(x)

bs_sk   <- boot_skew(x, 2000)      # default bias correction on
bs_ku   <- boot_kurt(x, 2000)

jk_sk <- jack_skew(x)
jk_ku <- jack_kurt(x)


g2 <- tryCatch(MASS::mvrnorm(1000, 
				               mu = c(0, 0), 
				            Sigma = matrix(c(1, -0.4, -0.4, 1), nrow = 2)),
                    error = function(e) {
						         message("Error in mvrnorm for group 2: ", e)
						         return(NA)
					  })
jack_cor(g2, bias.correct = TRUE, return.replicates = FALSE)
boot_cor(g2, B = 10000, bias.correct = TRUE, return.replicates = FALSE)

sk_res
bs_sk
jk_sk

ku_res
bs_ku
jk_ku