# Created by Yixin Sun and Eric Karsten in October 2018
# code to calculate bias-corrected confidence intervals
library(tidyverse)
library(furrr)

# bca steps
# 1. calculate theta_0
# 2. sample B times from df - make this into a list and a=pass into futures_map
# 3. Calculate B empirical thetas by passing samples into est(formula, data = df)
# 4. compute z_0 from {theta_B}
# 5. compute sigma_boot, empirical standard deviation from empirical thetas
# 6. compute a from jackknife
# 7. compute sigma_jackknife
# 8. compute internal standard error

# split up the bias-correction process from the acceleration calculation
uggs <- function(df, B, formula, est, ..., m = nrow(df), alpha = c(0.025, 0.05, 0.1)){
	t0 <- est(formula, data = df)

}
