# Created by Yixin Sun and Eric Karsten in October 2018
# code to calculate bias-corrected confidence intervals
library(tidyverse)
library(furrr)

plan(multiprocess)


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

# function for calculating bca level-alpha confidence interval endpoint - aka G inverse
furrrgi <- function(alpha, theta_star, t0, a){
  B <- length(theta_star)
  alpha <- alpha[alpha < 0.5]
  alpha <- c(alpha, 0.5, rev(1 - alpha))
  zalpha <- qnorm(alpha)
  sdboot0 <- sd(theta_star)

  # compute bca confidence limits
  z0 <- qnorm(sum(theta_star < t0)/B)
  phi <- pnorm(z0 + (z0 + zalpha)/(1 - a * (z0 + zalpha)))
  perc <- round(phi * B)
  theta_bca <- sort(theta_star)[perc]

  # also compute usual standard errors
  theta_std <- t0 + sdboot0 * qnorm(alpha)
    
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
# This is done K times, and averaged to get jackknife estimates 
furrrie <- function(theta_star, B, J, ajack, alpha, t0){
  indicies <-
        sample(B, size = B) %>%
        matrix(ncol = J) %>%
        as_tibble() %>%
        as.list()

  internal_jack <- function(alpha, t0, theta_star, drop_index, ajack){
    ttj <- theta_star[-drop_index]
    Bj <- length(ttj)
    sdboot <- sd(ttj)
    
    limit_j <- unlist(furrrgi(alpha, ttj, t0, ajack)[,"bca"])
    z0 <- qnorm(sum(ttj < t0)/B)

    stats_j <- c(sdboot, z0)

    return(list(limit = limit_j, stats = stats_j))
  }

  int_sd <- future_map(indicies, function(x) internal_jack(alpha, t0, theta_star, x, ajack))

  limits <- 
    int_sd %>%
    map(function(x) x$`limit`) %>%
    reduce(cbind) %>%
    split(seq(nrow(.))) %>%
    map_dbl(sd)
  limits <- limits * (J - 1) / sqrt(J)

  stats <- 
    int_sd %>%
    map(function(x) x$stats) %>%
    reduce(cbind) %>%
    split(seq(nrow(.))) %>%
    map_dbl(sd)
    stats <- stats * (J - 1) / sqrt(J)
    names(stats) <- c("sdboot", "z0")

    return(list(jacksd = limits, jsd = stats))

}


# function that takes in data, formula, function to estimate, jcount, jreps,
# and spits out a jacknife theta estimates

# jcount is how many jacknifes we want to do, defaulting to number of rows,
# in the old code this was `m`

# jreps is the number of times we want to do the random jacknife 

furrrjack <- function(df, formula, est, ..., jcount, jreps, progress = progress) {
  jack <- function(df, jcount) {
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
                     est(formula, df[-unlist(x),]), .progress = progress) %>% 
      unname()

    # calculate a
    aj <- a(theta_j, jcount)
    ssj <- se_jack(theta_j, jcount)
  
    return(list(a = aj, ssj = ssj))
  }

  # run calculation of jacknife theta estimates jrep times
  aij <- rerun(jreps, jack(df, jcount))

  # jack() returns the a and jackknife estimate of standard error
  # extract these 2 elements and average the jreps 
  a <- 
    aij %>%
    map_dbl(function(x) x[[1]]) %>%
    mean

  sdjack <- 
    aij %>%
    map_dbl(function(x) x[[2]]) %>%
    mean
  
  return(list(a = a, sdjack = sdjack))
}


# function that turns jacknife thetas into `a`
a <- function(theta_j, jcount) {
  theta_dot <- (mean(theta_j) - theta_j) * (jcount - 1)
  aj <- 1/6 * sum((theta_j - theta_dot)^3) / ((sum((theta_j - theta_dot)^2))^1.5)
  return(aj)
}

# function that gives the standard error
se_jack <- function(theta_j, jcount) {
  n <- length(theta_j)
  theta_dot <- (mean(theta_j) - theta_j) * (jcount - 1)
  se <- sqrt(sum(theta_dot^2)) / sqrt(jcount * (jcount - 1))
  return(se)
}
