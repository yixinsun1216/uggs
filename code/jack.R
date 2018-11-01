# Created by Yixin Sun and Eric Karsten in October 2018
# code to calculate bias-corrected confidence intervals
library(tidyverse)
library(furrr)

plan(multiprocess)

# function that takes in data, formula, function to estimate, jcount, jreps,
# and spits out a jacknife theta estimates

# jcount is how many jacknifes we want to do, defaulting to number of rows,
# in the old code this was `m`

# jreps is the number of times we want to do the random jacknife 

furrrjack <- function(df, formula, est, ..., jcount) {
  jack <- function(drop_index) {
    return(est(formula, df[-unlist(drop_index),]))
  }
  
  # create list of index sets to be jacknifed
  n <- nrow(df)
  r <- n %% jcount
  indicies <-
    sample(1:n, n-r) %>%
    matrix(jcount) %>%
    t() %>%
    as_tibble() %>%
    as.list()

  # compute theta for each jacknife sample
  theta_j <- future_map_dbl(indicies , jack) %>% unname()
  return(theta_j)
}

# function that turns jacknife thetas into `a`
a <- function(theta_j) {
  theta_dot <- mean(theta_j)
  return(1/6 * sum((theta_j - theta_dot)^3)/
    ((sum((theta_j - theta_dot)^2))^1.5))
}

# function that gives the standard error
se_jack <- function(theta_j) {
  n <- length(theta_j)
  theta_dot <- mean(theta_j)
  se <- sqrt( (n-1)/n * sum((theta_j - theta_dot)^2))
  return(se)
}


