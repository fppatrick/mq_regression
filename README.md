# M-Quantile Regression (`mquantile.R`)

An R function for robust M-quantile regression, combining quantile regression with M-estimation to handle outliers, heteroscedasticity and non-Gaussian conditional distribution.

---

## Description
The `mqlm` function fits M-quantile regression models using iterative reweighted least squares (IRLS). It extends traditional quantile regression by incorporating robustness through M-estimation, making it suitable for data with outliers or non-constant variance.

Key features:
- Estimates conditional quantiles (e.g., median, 90th percentile).
- Supports robustness tuning via influence functions (`psi.huber`, `psi.bisquare`).
- Allows case weights and variance weights for flexible modeling. 
- Provides convergence control via tolerance (`acc`) and maximum iterations (`maxit`).

---

## Installation
Ensure the `MASS` package is installed:
```R
install.packages("MASS")

source("mquantile.R")

library(MASS)
set.seed(123)
x <- rnorm(100)
y <- 2 * x + rnorm(100, sd = 0.5)  # Linear relationship with noise

# Fit M-quantile regression for the median (q = 0.5)
fit <- mqlm(x = cbind(x), y = y, q = 0.5)

# Print coefficients
print(fit$coefficients)

# Plot results
plot(x, y, main = "M-Quantile Regression (Median)")
lines(x, fit$fitted.values, col = "red", lwd = 2)
