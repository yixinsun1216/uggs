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

# passing bootstrapped sample into estimating function
furrrboot <- function(df, i, formula, est){
	dfx <- df[i,]
	theta_star <- est(formula, dfx)
	n <- nrow(dfx)

	# indicidence matrix - count of data used
	# want to make sure the tabulation contains every index, so we first add in 
	# vector 1:n and then subtract 1 from count of each index 
	Yj <- table(c(i, 1:n)) - 1

	return(list(theta = theta_star, Yj = Yj))
}

# function for calculating bca level-alpha confidence interval endpoint
furrr_bca <- function(alpha, theta_star, t0, a){
	alpha <- alpha[alpha < 0.5]
    alpha <- c(alpha, 0.5, rev(1 - alpha))
    zalpha <- qnorm(alpha)
    sdboot0 <- sd(theta_star)

    # compute bca confidence limits
    z0 <- stats::qnorm(sum(theta_star < t0)/B)
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
	dfx_indicies <- rerun(B, sample(n, n, replace = TRUE))

	# pass each sample into estimator
	boot <- 
	  dfx_indicies %>%
	  future_map(function(x) furrrboot(reg_data, x, formula, est), 
	  								   .progress = progress)

	# furrrboot returns theta and count vector Y - extract these theta
	theta_boot <- map_dbl(boot, function(x) x[[1]])
	se_boot <- sd(theta_boot)

	# use total count vector to calculate sample error
	Y <-
	  map(boot, function(x) x[[2]]) %>% 
	  reduce(function(x, y) rowSums(cbind(x, y)))

	# weight theta by count vector
	tY <- 
	  map(boot, function(x) x[[1]] * x[[2]]) %>%
	  reduce(function(x, y) rowSums(cbind(x, y))) 

	theta_boot_mean <- mean(theta_boot)
	tY_avg <- tY / B
	Y_avg <- Y / B
	s <- n * (tY_avg - theta_boot_mean * Y_avg)
	u <- 2 * theta_boot_mean - s
	sdu <- sqrt(sum(u ^ 2)) / n

	# calculate bias corrected estimator and bind with sdu
	ustat <- 2 * t0 - mean(theta_boot) 
	bootstat <- c(ustat, sdu)
	names(bootstat) <- c("ustat", "sdu")

	# calculate a
	theta_j <- furrrjack(df, formula, est, jcount = jcount)
	ajack <- a(theta_j)

	# use a to calculate bca level-alpha CI
	limits <- furrr_bca(alpha, theta_boot, t0, ajack)


	# calculate B.mean
	B.mean <- c(as.integer(B), mean(theta_boot))

	stats <- bootstat
	bca_output <- list(limits = limits, stats = stats, B.mean = B.mean)
	return(bca_output)
}
