

# Install packages
pacman::p_load(moments, PearsonDS)

## target central moments
m      <- c(    mean = 10,           # μ
            variance = 1,            # σ²
            skewness = 1,            # γ1  (0 ⇒ symmetry)
            kurtosis = 6)            # γ2  (6 ⇒ excess kurtosis = 3)

x <- rpearson(1e5, moments = m)      # draw 100,000 values
moments::kurtosis(x, na.rm = TRUE)   # sample estimate (should be ≈ 6)
moments::skewness(x, na.rm = TRUE)   # sample estimate (should be ≈ 0)
mean(x, na.rm = TRUE)                # sample estimate (should be ≈ 0)
var(x, na.rm = TRUE)                 # sample estimate (should be ≈ 1)
hist(x)
