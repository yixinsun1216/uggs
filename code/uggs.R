# Created by Yixin Sun and Eric Karsten in October 2018
# code to calculate bias-corrected confidence intervals
library(tidyverse)
library(furrr)

# bca steps
# 1. calculate t_0
# 2. sample B times from df - make this into a list and a=pass into futures_map
# 3. Calculate B empirical ts by passing samples into est(formula, data = df)
# 4. compute z_0 from {t_B}
# 5. compute sigma_boot, empirical standard deviation from empirical ts
# 6. compute a from jackknife
# 7. compute sigma_jackknife
# 8. compute internal standard error

num.workers <- 4
plan(multiprocess(workers = eval(num.workers)))

bias_ugg <- function(sample, formula, est){
	t_star <- future_map_dbl(sample, function(x) est(formula, data = x))
}

# function for calculating bca level-alpha confidence interval endpoint
ci_ugg <- function(alpha, t_star, t0, a){
	alpha <- alpha[alpha < 0.5]
    alpha <- c(alpha, 0.5, rev(1 - alpha))
    zalpha <- qnorm(alpha)
    sdboot0 <- sd(t_star)

    # compute bca confidence limits
	phi <- stats::pnorm(z0 + (z0 + zalpha)/(1 - a * (z0 + zalpha)))
	perc <- round(phi * B)
	t_bca <- sort(t_star)[perc]

	# also compute standard errors
	t_std <- t0 + sdboot0 * stats::qnorm(alpha)
    ugg_ci <- cbind(t_alpha, t_std)
    dimnames(ugg_ci) <- list(alpha, c("bca", "std"))

    return(ugg_ci)
}

# split up the bias-correction process from the acceleration calculation
uggs <- function(df, B, formula, est, ..., m = nrow(df), alpha = c(0.025, 0.05, 0.1), progress = TRUE){
	n <- nrow(df)
	t0 <- est(formula, data = df)

	# create B samples from df
	dfx <- rerun(B, sample_n(df, n, replace = TRUE))

	# pass each sample into estimator
	t_star <- future_map_dbl(dfx, function(x) est(formula, data = x), .progress = progress)

	# figure out how to calculate sampling error


}
