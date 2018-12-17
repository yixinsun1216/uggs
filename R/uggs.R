## Version of November 8, 2018


#' @title Nonparametric bias-corrected and accelerated bootstrap confidence
#'   limits
#'
#' @description This routine computes nonparametric confidence limits for
#'   bootstrap estimates. 
#' 
#' @details 
#' Bootstrap confidence limits correct the standard method of confidence 
#' intervals in three ways:
#' 
#' \enumerate{
	#' \item param the bootstrap cdf, \eqn{G}
	#' \item the bias-correction number \eqn{z_{0}} 
	#' \item the acceleration number \eqn{a} that measures the rate of change in 
	#'   \eqn{\sigma_{t_0}} as the data changes.
#' }
#' 
#' @param df a dataframe with \eqn{n} rows, assumed to be independently sampled 
#'  	from the target population.
#' @param B number of bootstrap replications
#' @param est function of the estimating equation, \eqn{\hat{\theta} = est(x)}, 
#' 		which returns a real value for the parameter of interest
#' @param ... additional arguments for est
#' @param jcount value used in calculating a. Because n can get very large, 
#'  	calculating \eqn{n} jackknife values can be slow. A way to speed up the 
#'  	calculation is to collect the \eqn{n} observations into \eqn{jcount} 
#'  	groups and deleting each group in turn. Thus we only 
#'  	evaluate \eqn{jcount} calculations instead of \eqn{n}.
#' @param jreps number of repetitions of grouped jackknives. These \eqn{jreps}
#'  	calculations are averaged to obtain our \eqn{\hat{a}}.
#' @param iereps a separate jackknife calculation estimates the 
#'  	internal standard error. The \eqn{B}-length vector \eqn{theta^*} is 
#' 		randomly grouped into \eqn{J} groups, and each group is deleted in turn
#' 		to recompute our estimates. This is done \eqn{iereps} times and 
#' 		averaged to compute the final jackknife estimates. 
#' @param J the number of groups \eqn{B}-length vector \eqn{theta^*} is 
#' 		partitioned into to calculate internal standard error
#' @param alpha percentiles to be computed for the confidence limits. Providing
#' 		alpha values below 0.5 is sufficient; upper limits are automatically 
#' 		computed
#' @param progress logical for a progress bar in bootstrap calculations
#' @param num_workers the number of workers used for parallel processing
#' @param ie_calc logical for whether or not we should calculate internal 
#' 		standard errors
#' 
#' 
#' @return
#' 
#' \itemize{
	#' \item limits: four columns housing information on confidence limits
	#' \itemize{	
	#'    \item `bca` shows the empirical bca confidence limits
	#' 		at the alpha percentiles
	#'    \item `std` shows the the standard confidence limits, 
	#'      \eqn{\hat{\theta} + \hat{\sigma}z_{\alpha}}
	#'    \item `pct` gives the percentiles of the sorted B bootstrap replications 
	#'      that correspond to `bca`
	#'    \item `pct`, gives the percentiles of the ordered B bootstrap replications
	#'      corresponding to the bca limits
	#'    \item `jacksd` is internal standard errors for the bca limits
	#' }
	#' \item stats: top line of stats shows 5 estimates, and bottom line gives the 
	#' 		internal standard errors for the five quantities below:
	#' \itemize{
	#'    \item theta: \eqn{f(x)}, original point estimate of the parameter of
	#'     interest
	#' 	  \item `sdboot` is the bootstrap estimate of standard error;
	#' 	  \item `z0` is the bca bias correction value, in this case quite
	#'     negative
	#' 	  \item `a` is the _acceleration_, a component of the bca
	#'     limits (nearly zero here)
	#' 	  \item `sdjack` is the jackknife estimate
	#'     of standard error for theta. 
	#'}
	#' \item B.mean: bootstrap sample size B, and the mean of the B
	#'     bootstrap replications \eqn{\hat{\theta^*}}
	#'
	#' \item ustats: The bias-corrected estimator `2 * t0 - mean(tt)`,
	#'     and an estimate, `sdu`, of its sampling error
#' }
#' 
#' @examples
#' library(lfe)
#' library(uggs)
#' 
#' ## create covariates
#' x1 <- rnorm(1000)
#' x2 <- rnorm(length(x1))
#' 
#' ## fixed effects
#' fe <- factor(sample(20, length(x1), replace=TRUE))
#' 
#' ## effects for fe
#' fe_effs <- rnorm(nlevels(fe))
#' 
#' ## creating left hand side y
#' u <- rnorm(length(x1))
#' y <- 2 * x1 + x2 + fe_effs[fe] + u
#' 
#' # create dataframe to pass into uggs
#' df_test <- as.data.frame(cbind(y, x1, x2, fe))
#' 
#' # function that returns parameter of interest, x1
#' est_test <- function(df){
#' 	m <- felm(y ~ x1 + x2 | fe, df)
#' 	as.numeric(coef(m)["x1"])
#' }
#' 
#' x1_boot <- uggs(df_test, 1000, est_test, jcount = 40, jreps = 5)
#' 
#' @references Efron, Bradley, and Trevor J. Hastie. Computer Age Statistical 
#' 		Inference: Algorithms, Evidence, and Data Science. Cambridge University 
#' 		Press, 2017.
#' @references Efron, Brad, and Balasubramanian Narasimhan. The Automatic 
#' 		Construction of Bootstrap Confidence Intervals. 2018.
#' @importFrom magrittr %>%
#' @importFrom tibble tibble as_tibble
#' @importFrom purrr map map2_df map_df rerun
#' @importFrom stats cov dnorm lm pnorm qnorm runif sd var
#' @importFrom furrr future_map_dbl future_map 
#' @importFrom future plan multiprocess
#' @importFrom dplyr bind_rows group_by summarise right_join left_join
#' @importFrom tidyr replace_na
#' 



# Main function
#' @export
uggs <- function(df, B, est, ..., jcount = nrow(df), jreps = 5,
                 iereps = 2, J = 10, alpha = c(0.025, 0.05, 0.1),
                 progress = TRUE, num_workers = 4, ie_calc = FALSE){

    # this is a silly, but necessary workaround
    num_workers <<- num_workers
	plan(future::multiprocess(workers = num_workers))

	n <- nrow(df)
	t0 <- est(df, ...)

	# create B samples from df and pass each sample into estimator
	print("Bootstrapping B samples and re-estimating theta")
	bootstrap_indicies <- rerun(B, sample(n, n, replace = TRUE))
	theta_boot <-
	  bootstrap_indicies %>%
	  future_map_dbl(function(i) {est(df[i,], ...)}, .progress = progress)
	z0 <- qnorm(sum(theta_boot < t0)/B)
	sdboot <- sd(theta_boot)
	theta_boot_mean <- mean(theta_boot)

	# calculate a and jacknife standard error
	print(paste("Calculating acceleration value using", jreps, 
				"rounds of jackknifing"))
	jackoutput <-
	  rerun(jreps, furrrjack(df, est, jcount, progress, ...)) %>%
	  bind_rows() %>%
	  colMeans() %>%
	  as.list()
	ajack <- jackoutput$a

	# use a to calculate bca level-alpha CI
	limits <- furrrgi(alpha, theta_boot, t0, ajack)

	# calculate internal errors and average across iereps calculations
	# bind together stats and limits for outputting
	if(ie_calc){
		print(paste("Estimating internal error of confidence limits using", 
					iereps, "rounds of jackknifing"))
		ie <-
		  rerun(iereps, furrrie(theta_boot, B, J, ajack, alpha, t0, progress)) %>%
		  bind_rows() %>%
		  map_df(mean)
		jacksd <- unlist(ie)[1:(length(alpha) * 2 + 1)]
		limits <- cbind(limits, jacksd)
		jsd <- c(0, ie$sdboot, ie$z0, 0, 0)
		stats <- t(cbind(c(t0, sdboot, z0, ajack, jackoutput$se), jsd))
		rownames(stats) <- c("est", "jsd")	
	}else{
		stats <- cbind(t0, sdboot, z0, ajack, jackoutput$se)
		rownames(stats) <- "est"
	}
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

	# compute proportion of bootstrapped estimates that are less than conf limit
	pct <-
	  theta_bca %>%
	  map(function(x) sum(theta_star <= x) / B)

	ugg_ci <- cbind(theta_bca, theta_std, pct)
	dimnames(ugg_ci) <- list(alpha, c("bca", "std", "pct"))

	return(ugg_ci)
}


# function that takes in data, function to estimate, jcount, jreps,
# and spits out the value accelaration value a and jacknife standard deviation
# jcount is how many jacknifes we want to do, defaulting to number of rows
# jreps is the number of times we want to do the random jacknife
furrrjack <- function(df, est, jcount, progress, ...) {
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
	    est(df[-unlist(x),], ...), .progress = progress) %>%
	  unname()

	# calculate acceleration a and jacknife standard error
	theta_dot <- (mean(theta_j) - theta_j) * (jcount - 1)
	a_j <- 1/6 * sum((theta_j - theta_dot)^3) /
	((sum((theta_j - theta_dot)^2))^1.5)
	se_j <- sqrt(sum(theta_dot^2)) / sqrt(jcount * (jcount - 1))

	return(list(a = a_j, se = se_j))
}


# Function for calculating the internal error of the confidence limits
# These use a different jackknife calculation, which are based on the
# original B bootstrap replications.
# The B-vector theta_star is divided into J groups, and each group is deleted
# in turn to recompute the limits.
furrrie <- function(theta_star, B, J, ajack, alpha, t0, progress){
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
		limit_j <- as.numeric(unlist(furrrgi(alpha, ttj, t0, ajack)[,"bca"]))
		out <- c(limit_j, sdboot, z0)

		nalpha <- alpha[alpha < 0.5]
		nalpha <- c(nalpha, 0.5, rev(1 - nalpha))
		names(out) <- c(as.character(nalpha),"sdboot", "z0")

		return(as.list(out))
	}

	int_sd <-
	  future_map(indicies, function(x) internal_jack(x), .progress = progress) %>%
	  bind_rows() %>%
	  map_df(function(x) {sd(x) * (J-1)/(sqrt(J))})

	return(int_sd)
}
