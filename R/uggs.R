## Version of November 8, 2018


#' @title Nonparametric bias-corrected and accelerated bootstrap confidence
#'   limits
#'
#' @description This routine computes nonparametric confidence limits for
#'   bootstrap estimates.
#'
#' @details Bootstrap confidence limits correct the standard method of
#' confidence intervals in three ways
#'
#' - the bootstrap cdf, \eqn{G} - the bias-correction number \eqn{z_{0} - the
#' acceleration number \eqn{a} that measures the rate of change in
#' \eqn{\sigma_{t_0}} as the data changes.
#'


#' @importFrom magrittr %>%
#' @importFrom tibble tibble
#' @importFrom purrr map map2_df map_df
#' @importFrom stats cov dnorm lm pnorm qnorm runif sd var

# Main function
#' @export
uggs <- function(df, B, est, ..., jcount = nrow(df), jreps = 5,
                 iereps = 2, J = 10, alpha = c(0.025, 0.05, 0.1),
                 progress = TRUE, num.workers = 4){

	plan(multiprocess(workers = eval(num.workers)))

	n <- nrow(df)
	t0 <- est(df, ...)

	# create B samples from df and pass each sample into estimator
	bootstrap_indicies <- rerun(B, sample(n, n, replace = TRUE))
	theta_boot <-
	  bootstrap_indicies %>%
	  future_map_dbl(function(i) {est(df[i,])}, .progress = progress)
	z0 <- qnorm(sum(theta_boot < t0)/B)
	sdboot <- sd(theta_boot)
	theta_boot_mean <- mean(theta_boot)

	# calculate a and jacknife standard error
	jackoutput <-
	  rerun(jreps, furrrjack(df, est, jcount, n)) %>%
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

	jacksd <- unlist(ie)[1:(length(alpha) * 2 + 1)]
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

# ===========================================================================
# Helper functions
# ===========================================================================
# function for calculating bca level-alpha confidence interval endpoint: G^{-1}
furrrgi <- function(alpha, theta_star, t0, a){
	B <- length(theta_star)
	alpha <- alpha[alpha < 0.5]
	alpha <- c(alpha, 0.5, rev(1 - alpha))
	zalpha <- qnorm(alpha)
	sdboot0 <- sd(theta_star)

	# compute bca confidence limits
	z0 <- qnorm(sum(theta_star < t0)/B)
	phi <- pnorm(z0 + (z0 + zalpha)/(1 - a * (z0 + zalpha)))
	theta_bca <- sort(theta_star)[round(phi * B)]

	# also compute usual standard errors
	theta_std <- t0 + sdboot0 * qnorm(alpha)

	# compute proportion of bootstrapped estimates which are less than conf limit
	pct <-
	  theta_bca %>%
	  map(function(x) sum(theta_star <= x) / B)

	ugg_ci <- cbind(theta_bca, theta_std, pct)
	dimnames(ugg_ci) <- list(alpha, c("bca", "std", "pct"))

	return(ugg_ci)
}


# Function for calculating the internal error of the confidence limits
# These use a different jackknife calculation, which are based on the
# original B bootstrap replications.
# The B-vector theta_star is divided into J groups, and each group is deleted
# in turn to recompute the limits.
furrrie <- function(theta_star, B, J, ajack, alpha, t0){
	indicies <-
	  sample(B, B - B %% J) %>%
	  matrix(J) %>%
	  t() %>%
	  as_tibble() %>%
	  as.list()

	internal_jack <- function(drop_index){
		ttj <- theta_star[-drop_index]
		Bj <- length(ttj)
		sdboot <- sd(ttj)
		z0 <- qnorm(sum(ttj < t0)/B)
		limit_j <- unlist(furrrgi(alpha, ttj, t0, ajack)[,"bca"]) %>% as.numeric()
		out <- c(limit_j, sdboot, z0)

		nalpha <- alpha[alpha < 0.5]
		nalpha <- c(nalpha, 0.5, rev(1 - nalpha))
		names(out) <- c(as.character(nalpha),"sdboot", "z0")

		return(as.list(out))
	}

	int_sd <-
	  future_map(indicies, function(x) internal_jack(x)) %>%
	  bind_rows %>%
	  map_df(function(x) {sd(x) * (J-1)/(sqrt(J))})

	return(int_sd)
}


# function that takes in data, function to estimate, jcount, jreps,
#   and spits out the value accelaration value a and jacknife standard deviation
# jcount is how many jacknifes we want to do, defaulting to number of rows
# jreps is the number of times we want to do the random jacknife
furrrjack <- function(df, est, jcount,...) {
	n <- nrow(df)
	r <- n %% jcount

	# create list of index sets to be jacknifed
	indicies <-
	  sample(1:n, n-r) %>%
	  matrix(jcount) %>%
	  t() %>%
	  as_tibble() %>%
	  as.list()

	# estimates theta without sampled subgroup
	theta_j <-
	  indicies %>%
	  future_map_dbl(function(x)
	    est(df[-unlist(x),])) %>%
	  unname()

	# calculate acceleration a and jacknife standard error
	theta_dot <- (mean(theta_j) - theta_j) * (jcount - 1)
	a_j <- 1/6 * sum((theta_j - theta_dot)^3) /
	((sum((theta_j - theta_dot)^2))^1.5)
	se_j <- sqrt(sum(theta_dot^2)) / sqrt(jcount * (jcount - 1))

	return(list(a = a_j, se = se_j))
}
