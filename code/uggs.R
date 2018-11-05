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


# split up the bias-correction process from the acceleration calculation
uggs <- function(df, B, formula, est, ..., jcount = nrow(df), jreps = 5, K =2, J = 10,
	alpha = c(0.025, 0.05, 0.1), progress = TRUE){
	## Save rng state
    if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        runif(1)
    seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)

	n <- nrow(df)
	t0 <- est(formula, df)

	# create B samples from df and passes each sample into estimator
	dfx_indicies <- rerun(B, sample(n, n, replace = TRUE))
	boot <- 
	  dfx_indicies %>%
	  future_map(function(x) furrrboot(reg_data, x, formula, est), 
	  								   .progress = progress)

	# furrrboot returns list with theta and count vector Y - extract theta
	theta_boot <- map_dbl(boot, function(x) x[[1]])
	se_boot <- sd(theta_boot)

	# calculate a
	jackoutput <- furrrjack(df, formula, est, jcount = jcount, jreps = jreps, 
						 	progress = progress)
	ajack <- jackoutput$`a`


	# use a to calculate bca level-alpha CI
	limits <- furrrgi(alpha, theta_boot, t0, ajack)

	# still have to fill in sdboot
    z0 <- qnorm(sum(theta_boot < t0)/B)
    sdboot <- sd(theta_boot)

	# calculate internal errors and average across K calculations
	ie <- rerun(K, internal_error(theta_boot, B, J, ajack, alpha, t0))
	jacksd <- 
	  ie %>%
	  map(function(x) x$`jacksd`) %>%
	  reduce(cbind) %>%
	  split(seq(nrow(.))) %>%
	  map_dbl(mean)
	limits <- cbind(limits, jacksd)

	jsd <- 
	  ie %>%
	  map(function(x) x$jsd) %>%
	  reduce(cbind) %>%
	  split(seq(nrow(.))) %>%
	  map_dbl(mean)
	jsd <- c(0, jsd, 0, 0)
	stats <- t(cbind(c(t0, sdboot, z0, ajack, jackoutput$sdjack), jsd))
	rownames(stats) <- c("est", "jsd")
	colnames(stats) <- c("theta", "sdboot", "z0", "a", "sdjack")

	
    # use total count vector to calculate sample error aka black magic sdu
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
	ustats <- c(ustat, sdu)
	names(ustats) <- c("ustat", "sdu")

	options(scipen=999)
	B.mean <- c(B, mean(theta_boot))


	bca_output <- list(limits = limits, stats = stats, B.mean = B.mean, ustats = ustats)
	return(bca_output)
}
