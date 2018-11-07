# Created by Yixin Sun and Eric Karsten in October 2018
# code to calculate bias-corrected confidence intervals
library(tidyverse)
library(furrr)

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
    names(out) <- c(as.character(alpha),"sdboot", "z0")
    return(as.list(out))
  }

  int_sd <-
    future_map(indicies, function(x) internal_jack(x)) %>%
    bind_rows %>%
    map_df(function(x) {sd(x) * (J-1)/(sqrt(J))})

  return(int_sd)
}

# function that takes in data, formula, function to estimate, jcount, jreps,
#   and spits out the value accelaration value a and jacknife standard deviation
# jcount is how many jacknifes we want to do, defaulting to number of rows
# jreps is the number of times we want to do the random jacknife 
furrrjack <- function(df, formula, est, jcount,...) {
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
  
  # calculate acceleration a and jacknife standard error
  theta_dot <- (mean(theta_j) - theta_j) * (jcount - 1)
  a_j <- 1/6 * sum((theta_j - theta_dot)^3) /
    ((sum((theta_j - theta_dot)^2))^1.5)
  se_j <- sqrt(sum(theta_dot^2)) / sqrt(jcount * (jcount - 1))
  
  return(list(a = a_j, se = se_j))
}
