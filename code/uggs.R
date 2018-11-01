# Created by Yixin Sun and Eric Karsten in October 2018
# code to calculate bias-corrected confidence intervals
library(tidyverse)
library(furrr)

if(Sys.getenv("USERNAME") == "Yixin Sun"){
	root <- "C:/Users/Yixin Sun/Documents/Github/uggs"
	ddir <- "C:/Users/Yixin Sun/Documents/Dropbox/texas"
}

source(file.path(root, 'code/jack.R'))
# bca steps
# 1. calculate theta_0
# 2. sample B times from df - make this into a list and pass into futures_map
# 3. Calculate B empirical ts by passing samples into est(formula, data = df)
# 4. compute z_0 from {theta_B}
# 5. compute sigma_boot, empirical standard deviation from empirical ts
# 6. compute a from jackknife
# 7. compute sigma_jackknife
# 8. compute internal standard error

num.workers <- 4
plan(multiprocess(workers = eval(num.workers)))


# function for calculating bca level-alpha confidence interval endpoint
ci_ugg <- function(alpha, theta_star, t0, a){
	alpha <- alpha[alpha < 0.5]
    alpha <- c(alpha, 0.5, rev(1 - alpha))
    zalpha <- qnorm(alpha)
    sdboot0 <- sd(theta_star)

    # compute bca confidence limits
	phi <- stats::pnorm(z0 + (z0 + zalpha)/(1 - a * (z0 + zalpha)))
	perc <- round(phi * B)
	theta_bca <- sort(theta_star)[perc]

	# also compute usual standard errors
	theta_std <- t0 + sdboot0 * stats::qnorm(alpha)
    ugg_ci <- cbind(theta_bca, theta_std)
    dimnames(ugg_ci) <- list(alpha, c("bca", "std"))

    return(ugg_ci)
}


# split up the bias-correction process from the acceleration calculation
uggs <- function(df, B, formula, est, ..., jcount = nrow(df), jreps = 5,
	alpha = c(0.025, 0.05, 0.1), progress = TRUE){
	## Save rng state
    if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        stats::runif(1)
    seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)

	n <- nrow(df)
	t0 <- est(formula, df)

	# create B samples from df
	dfx <- rerun(B, sample_n(df, n, replace = TRUE))

	# pass each sample into estimator
	theta_boot <- future_map_dbl(dfx, function(x) est(formula, x), .progress = progress)
	se_boot <- sd(theta_boot)

	# calculate a
	theta_j <- furrrjack(df, formula, est, jcount = jcount)
	ajack <- a(theta_j)

	# use a to calculate bca level-alpha CI
	limits <- ci_ugg(alpha, theta_star, t0, ajack)


	# calculate bias-corrected estimator
	ustat <- 2 * t0 - mean(theta_boot)
	# figure out how to calculate sampling error sdu????
	stats <- ustat

	# calculate B.mean
	B.mean <- c(B, mean(theta_boot))

	bca_output <- list(limits = limits, stats = stats, B.mean = B.mean)
	return(bca_output)
}
