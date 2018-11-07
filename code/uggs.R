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
uggs <- function(df, B, formula, est, ..., jcount = nrow(df), jreps = 5,
                 iereps = 2, J = 10, alpha = c(0.025, 0.05, 0.1),
                 progress = TRUE){
	## Save rng state
  if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      runif(1)
  seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)

	n <- nrow(df)
	t0 <- est(formula, df)

	# create B samples from df and pass each sample into estimator
	bootstrap_indicies <- rerun(B, sample(n, n, replace = TRUE))
	theta_boot <- 
	  bootstrap_indicies %>%
	  future_map_dbl(function(i) {est(formula, df[i,])}, .progress = progress)
	z0 <- qnorm(sum(theta_boot < t0)/B)
	sdboot <- sd(theta_boot)
	theta_boot_mean <- mean(theta_boot)

	# calculate a and jacknife standard error
	jackoutput <-
	  rerun(jreps, furrrjack(df, formula, est, jcount, n)) %>%
	  bind_rows() %>%
	  colMeans() %>%
	  as.list()
	ajack <- jackoutput$a
	
	# use a to calculate bca level-alpha CI
	limits <- furrrgi(alpha, theta_boot, t0, ajack)

	# calculate internal errors and average across iereps calculations
	ie <-
	  rerun(iereps, furrrie(theta_boot, B, J, ajack, alpha, t0)) %>%
	  bind_rows() %>%
	  map_df(mean)
	
	jacksd <- unlist(ie)[1:length(alpha)]
	limits <- cbind(limits, jacksd)

	jsd <- c(0, ie$sdboot, ie$z0, 0, 0)
	stats <- t(cbind(c(t0, sdboot, z0, ajack, jackoutput$se), jsd))
	rownames(stats) <- c("est", "jsd")
	colnames(stats) <- c("theta", "sdboot", "z0", "a", "sdjack")
	
  # use total count vector to calculate sample error
	Y <-
	  bootstrap_indicies %>% 
	  unlist() %>%
	  c(1:n) %>%
	  table() - 1

	# weight theta by count vector
	tY <-
	  map2_df(theta_boot, bootstrap_indicies,
	          function(x,y) tibble(theta = x, index = y)) %>%
	  group_by(index) %>%
	  summarise(t_sum = sum(theta)) %>%
	  right_join(tibble(index = 1:n), by = "index") %>%
	  replace_na(list(t_sum = 0))

	tY_avg <- tY$t_sum / B
	Y_avg <- Y / B
	s <- n * (tY_avg - theta_boot_mean * Y_avg)
	u <- 2 * theta_boot_mean - s
	sdu <- sqrt(sum(u ^ 2)) / n

	# calculate bias corrected estimator and bind with sdu
	ustat <- 2 * t0 - theta_boot_mean 
	ustats <- c(ustat, sdu)
	names(ustats) <- c("ustat", "sdu")

	options(scipen=999)
	B.mean <- c(B, theta_boot_mean)

	bca_output <- list(limits = limits, stats = stats,
	                   B.mean = B.mean, ustats = ustats)
	return(bca_output)
}
