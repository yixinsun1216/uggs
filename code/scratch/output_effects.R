#===========================================================================
# code to test the hypothesese that
# (a) output is (conditionally) higher in auctioned leases
# (b) the difference can be explained by differences in firm assignment
#===========================================================================
library(tidyverse)
library(lubridate)
library(knitr)
library(kableExtra)
library(bcaboot)
library(lfe)
library(alpaca)
library(furrr)

#===========================================================================
# BASIC TEXAS SETUP
#===========================================================================
root <- getwd()
while(basename(root) != "texas") {
  root <- dirname(root)
}

source(file.path(root, "code", "paths.R"))

if(length(commandArgs(TRUE)) < 4) {
  num.workers <- 4
  estimator <- NA_real_
  control_spec <- NA_character_
  n_boots <- 500
} else {
  args <- commandArgs(TRUE)
  num.workers <- as.integer(args[1])
  estimator <- args[2]
  control_spec <- as.integer(args[3])
  n_boots <- as.integer(args[4])  
}

plan(multiprocess(workers = eval(num.workers)))

#===========================================================================
# LOAD DATA CLEANED IN SAMPLE_SELECTION.R
#===========================================================================
load(file.path(ddir, "generated_data/final_sample.Rda"))

reg_data <-
  final_sample %>%
  mutate(acres = Polygon_Acres,
         shape_quality = Polygon_Acres / Hull_Acres,
         Auction = if_else(Auction, 1, 0),
         Firm = NewFirm,
         EffDate = as.numeric(Effective_Date),
         Rental_Total = if_else(is.na(Rental_Total), 0, Rental_Total),
         Output_per_acre = total_deur_boe / acres) %>%
  mutate_at(vars(DistRoad:DistWaterbodies),
            funs(if_else(is.na(.), 100, ./1000))) %>%
  mutate(Developed_HI = Developed_HI + Developed_MI,
         Developed_LO = Developed_LI + Developed_OS) %>%
  filter(in_sample) %>%
  filter(!censored) %>%
  filter(Type != "UT") %>%
  mutate(NewFirm = as.integer(factor(NewFirm)),
         Grid10Yr = as.integer(Grid10Yr),
         Grid20Yr = as.integer(Grid20Yr),
         CountyYr = as.integer(CountyYr),
         YearQtr = as.integer(YearQtr))

#===========================================================================
# define specifications
#===========================================================================
base_controls <- "log(acres)"

shape_controls <- "shape_quality"

distance_controls <-
  c("DistRoad", "DistRiverStreams", "DistWaterbodies") %>%
  paste(collapse = " + ")

# here I am going to drop the first landcover category, Shrub_Scrub, because
# we need to drop one category and this guy is the bulk of the coverage
landcover_controls <-
  reg_data %>%
  select(Developed_HI, Developed_LO, CoverWater,
         CoverCultivated, CoverForest, CoverWetlands) %>%
  names %>%
  paste(collapse = " + ")

#gridsize <- "County"
gridsize <- "Grid10"
loctime_felm <- paste(gridsize, "Yr + YearQtr", sep = "")

loctime_feglm <- "CountyYr + YearQtr"


#===========================================================================
# helper functions
#===========================================================================
# this is a functor which generates a function g(X), corresponding to
# estimator(f, as.data.frame(X))
matrix_call <- function(estimator, f, fvars, ...) {
  return(function(X) estimator(f, as.data.frame(X[, fvars]), ...))
}

# a version which does it for two-formula models
matrix_call2 <- function(estimator, f1, f2, fvars, ...) {
  return(function(X) estimator(f1, f2, as.data.frame(X[, fvars]), ...))
}

# this routine generates a bootstrap sample of d, taking care to re-define the
# variable groupid to have an omitted small groups category that is computed
# within the bootstrap sample for any group with a singleton number of
# observations, or which has no variation in the varvar variable.  also, since
# we know we'll eventually have to have all columns be numeric, and factorized
# to be useful, the groupid column will be returned as a factor.  if
# do_bootstrap is FALSE, this routine still "cleans" the groupid of singleton
# values and zero varvar variation.  if incidence = TRUE, the return value is a
# list: the (possibly) bootstrapped data, and an incidence vector y which
# indicates how many times each observation in the original dataset appears in
# the bootstrap.  this second piece is useful if you want to use the bcajack2
# routine.  Note: at present, all the code that uses BCa bootstraps below does
# *NOT* use bcajack2, so this is probably not useful right now
bootstrap_fe_data <- function(d,
                              groupid,
                              varvar,
                              do_bootstrap = TRUE,
                              incidence = FALSE) {
  quo_groupid <- enquo(groupid)
  string_groupid <- quo_name(quo_groupid)
  quo_varvar <- enquo(varvar)
  
  if(do_bootstrap) {
    # bootstrap the rows and generate the incidence vector
    inds <- sample.int(nrow(d), replace = TRUE)
    y <- tabulate(inds, nrow(d))
    d <- d[inds, ]
  } else {
    # generate dummy incidence vector
    y <- rep(1, nrow(d))
  }

  # do several simulation friendly steps on this (possibly) bootstrapped data
  # 1. replace groupid with an equivalent integer
  # 2. define an outside group (i.e., 0) in groupid which collects all rows
  # that are singletons at the groupid level or have no variation in varvar
  d <-
    d %>%
    mutate(!! string_groupid := as.numeric(factor(!! quo_groupid))) %>%
    group_by(!! quo_groupid) %>%
    mutate(temp_N = n()) %>%
    mutate(temp_varvar = sd(!! quo_varvar)) %>%
    ungroup %>%
    replace_na(list(temp_varvar = 0)) %>%
    mutate(!! string_groupid := if_else(temp_N == 1 |
                                        temp_varvar == 0,
                                        0,
                                        !! quo_groupid)) %>%
    mutate(!! string_groupid := factor(!! quo_groupid)) %>%
    select(-temp_N, -temp_varvar)

  if(incidence) {
    return(list(d, y))
  } else {
    return(d)
  }
}

# this function replaces any NaN, NA, or Inf value with the finite, nomissing
# average in the vector
replace_na2 <- function(v) {
  n <- sum(is.na(v) | is.nan(v) | is.infinite(v))
  if(n > 0) {
    avg <- mean(v[!(is.na(v) | is.nan(v) | is.infinite(v))])
    newv <- v
    newv[(is.na(v) | is.nan(v) | is.infinite(v))] <- avg
    print(paste("replaced", n, "values with average", sep = " "))
    return(newv)
  } else {
    return(v)
  }
}

twopart <- function(f_drilled, f_log_output, d, B = 800) {
  bootstraps <- twopart_bootstraps(f_drilled, f_log_output, d, B = B)
  return(tibble(estimate = bootstraps[[1]],
                ci_lo = bootstraps[[2]]$lims["0.025", "bca"],
                ci_hi = bootstraps[[2]]$lims["0.975", "bca"]))
}

logit <- function(x) exp(x) / (1 + exp(x))
twopart_work <- function(f_drilled, f_log_output, d, run = 0) {

  # clean groups and bootstrap the data if need be
  do_bootstrap <- run > 0
  d <-
    bootstrap_fe_data(d,
                      NewFirm,
                      drilled, 
                      do_bootstrap = do_bootstrap)
  
  m_drilled <-
    feglm(f_drilled,
          d,
          family = binomial())

  logodds <- predict(m_drilled, type = "link")
  newprobs <-
    m_drilled$data %>%
    mutate(logodds = logodds,
           logodds1 = if_else(Auction == 1,
                              logodds,
                              logodds + coef(m_drilled)["Auction"]),
           logodds0 = if_else(Auction == 0,
                              logodds,
                              logodds - coef(m_drilled)["Auction"]))

  p1 <- with(newprobs, logit(logodds1))
  p0 <- with(newprobs, logit(logodds0))

  ape_drilled <- mean(log(p1) - log(p0))
  
  # estimate log output APE
  m_log_output <- felm(f_log_output, filter(d, drilled == 1))
  ape_log_output <- as.numeric(coef(m_log_output)["Auction"])
  
  theta <- ape_drilled + ape_log_output
  return(theta)
}

twopart_bootstraps <- function(f_drilled, f_log_output, d, B = 800) {
  bootstraps <-
    seq(1, B) %>%
    future_map_dbl(possibly(function(x) twopart_work(f_drilled,
                                                 f_log_output,
                                                 d,
                                                 x),
                            NA_real_),
               .progress = T,
               .options = future_options(seed = T))

  tt <- replace_na2(bootstraps)

  # compute bca bootstrap CIs
  t0 <- twopart_work(f_drilled, f_log_output, d, 0)
  fvars <- union(all.vars(f_drilled), all.vars(f_log_output))
  fvars <- union(fvars, "NewFirm")
  func <-
    function(x) matrix_call2(twopart_work,
                             f_drilled,
                             f_log_output,
                             fvars,
                             run = 0)(x)
  
  bca <- bcajack(as.matrix(d[, fvars]), tt, func, m = 80)
  return(list(t0, bca))
}

poisson2 <- function(f, d, B = 800) {
  bootstraps <- poisson2_bootstraps(f, d, B = B)
  return(tibble(estimate = bootstraps[[1]],
                ci_lo = bootstraps[[2]]$lims["0.025", "bca"],
                ci_hi = bootstraps[[2]]$lims["0.975", "bca"]))
}

poisson2_work <- function(f, d, run = 0) {
  # clean groups and bootstrap the data if need be
  do_bootstrap <- run > 0
  d <-
    bootstrap_fe_data(d,
                      NewFirm,
                      total_deur_boe,
                      do_bootstrap = do_bootstrap)
  
  m <-
    feglm(f,
          d,
          family = poisson())
  theta <- as.numeric(coef(m)["Auction"])
  return(theta)
}

poisson2_bootstraps <- function(f, d, B = 800) {
  # note the hack below to replace NAs with average values
  # deals with obnoxious warning messages in bcajack later on
  bootstraps <-
    seq(1, B) %>%
    future_map_dbl(possibly(function(x) poisson2_work(f, d, x), NA_real_),
                   .progress = T,
                   .options = future_options(seed = T))

  tt <- replace_na2(bootstraps)

  # compute bca bootstrap CIs
  t0 <- poisson2_work(f, d, 0)
  fvars <- union(all.vars(f), "NewFirm")
  func <- function(x) matrix_call(poisson2_work, f, fvars, run = 0)(x)
  bca <- bcajack(as.matrix(d[, fvars]), tt, func, m = 80)
  
  return(list(t0, bca))
}

make_table <- function(results, filename) {
  results %>%
    select(estimator, outcome, estimate, ci_lo, ci_hi, 
           firms, controls, spacetime) %>%
    rename(Estimator = estimator,
           Outcome = outcome,
           Estimate = estimate,
           `CI Lower` = ci_lo,
           `CI Upper` = ci_hi,           
           `Firm FE` = firms,
           Controls = controls,
           `Loc/Time` = spacetime) %>%
    kable(format = "latex",
          booktabs = TRUE,
          linesep = "",
          escape = F,
          align = c('l', 'c', 'c', 'c', 'c', 'c', 'l', 'l'),
          digits = 3) %>%
    writeLines(file.path(tdir, filename))
}

coefsummary_felm <- function(m) {
  coef(summary(m))
}

# come back here and refactor this to do:
# 1) basic linear FE specs, which we already have
# 2) poisson FE specs, which there is code for above but not a "models" method
# 3) revise twopart work to do a (consistent) set of FE models instead of
# polynomials.
# use bootstraps for versions 2 + 3
#===========================================================================
# estimate basic linear FE models
#===========================================================================
fe_models <- function(controls, cname = NULL, stname = NULL) {
  print("running fe models")
  controls <- paste(controls, collapse = " + ")
  base_rhs <-
    paste("Auction", controls, sep = " + ") %>%
    paste(loctime_felm, sep = " | ")
  firms_rhs <- paste(base_rhs, "NewFirm", sep = " + ")

  base_rhs <- paste(base_rhs, "0", gridsize, sep = " | ")
  firms_rhs <- paste(firms_rhs, "0", gridsize, sep = " | ")  

  fe_drilled <-
    paste("drilled", base_rhs, sep = " ~ ") %>%
    as.formula %>%
    felm(reg_data)

  fe_drilled_firms <-
    paste("drilled", firms_rhs, sep = " ~ ") %>%
    as.formula %>%
    felm(reg_data, exactDOF = T)

  fe_output <-
    paste("log(total_deur_boe)", base_rhs, sep = " ~ ") %>%
    as.formula %>%
    felm(filter(reg_data, drilled == 1))

  fe_output_firms <-
    paste("log(total_deur_boe)", firms_rhs, sep = " ~ ") %>%
    as.formula %>%
    felm(filter(reg_data, drilled == 1), exactDOF = T)

  fe_output_1 <-
    paste("log(total_deur_boe+1)", base_rhs, sep = " ~ ") %>%
    as.formula %>%
    felm(reg_data)

  fe_output_1_firms <-
    paste("log(total_deur_boe+1)", firms_rhs, sep = " ~ ") %>%
    as.formula %>%
    felm(reg_data, exactDOF = T)

  if(is.null(cname)) {
    cname <- as.character(str_length(controls))
  }

  if(is.null(stname)) {
    stname <- str_replace_all(location_time_FE, fixed(" "), "")
  }

  results <-
    list(fe_drilled, fe_drilled_firms,
         fe_output, fe_output_firms,
         fe_output_1, fe_output_1_firms) %>%
    map_dfr(~ as_tibble(coefsummary_felm(.), rownames = "var")) %>%
    filter(var == "Auction") %>%
    select(estimate = 2, stderr = 3) %>%
    mutate(ci_lo = estimate - 1.96 * stderr,
           ci_hi = estimate + 1.96 * stderr) %>%
    select(-stderr) %>%
    bind_cols(tibble(outcome = c("drilled", "drilled",
                                 "log(Output)|D", "log(Output)|D",
                                 "log(Output+1)", "log(Output+1)"),
                     estimator = rep("FE", 6),
                     firms = c("No", "Yes", "No", "Yes", "No", "Yes"),
                     controls = cname,
                     spacetime = stname))
  
  return(results)
}

#===========================================================================
# estimate a parametric log Pr + log E[Y|x] model, with bootstrap SEs
#===========================================================================
twopart_models <- function(controls,
                           cname = NULL,
                           stname = NULL,
                           B = n_boots) {
  print("running twopart models")
  controls <- paste(controls, collapse = " + ")
  
  base_drilled_rhs <-
    paste("Auction", controls, sep = " + ") %>%
    paste(loctime_feglm, sep = " | ")
  firms_drilled_rhs <- paste(base_drilled_rhs, "NewFirm", sep = " + ")

  base_log_output_rhs <-
    paste("Auction", controls, sep = " + ") %>%
    paste(loctime_feglm, sep = " | ")
  firms_log_output_rhs <- paste(base_log_output_rhs, "NewFirm", sep = " + ")

  nofirms <-
    twopart(as.formula(paste("drilled",
                             base_drilled_rhs,
                             sep = " ~ ")),
            as.formula(paste("log(total_deur_boe)",
                             base_log_output_rhs,
                             sep = " ~ ")),
            reg_data,
            B = B)

  firms <-
    twopart(as.formula(paste("drilled",
                             firms_drilled_rhs,
                             sep = " ~ ")),
            as.formula(paste("log(total_deur_boe)",
                             firms_log_output_rhs,
                             sep = " ~ ")),
            reg_data,
            B = B)

  if(is.null(cname)) {
    cname <- as.character(str_length(controls))
  }

  if(is.null(stname)) {
    stname <- str_replace_all(location_time_FE, fixed(" "), "")
  }

  results <-
    bind_rows(mutate(nofirms,
                     outcome = "log E[Output]",
                     estimator = "twopart",
                     firms = "No",
                     controls = cname,
                     spacetime = stname),
              mutate(firms,
                     outcome = "log E[Output]",
                     estimator = "twopart",
                     firms = "Yes",
                     controls = cname,
                     spacetime = stname))
  return(results)
}

#===========================================================================
# estimate bootstrapped poisson models
#===========================================================================
poisson2_models <- function(controls,
                            cname = NULL,
                            stname = NULL,
                            B = n_boots) {
  print("running poisson models")
  controls <- paste(controls, collapse = " + ")
  base_rhs <-
    paste("Auction", controls, sep = " + ") %>%
    paste(loctime_feglm, sep = " | ")
  firms_rhs <- paste(base_rhs, "factor(NewFirm)", sep = " + ")

  nofirms <-
    base_rhs %>%
    paste("total_deur_boe", ., sep = " ~ ") %>%
    as.formula %>%
    poisson2(reg_data, B = B)

  firms <-
    firms_rhs %>%
    paste("total_deur_boe", ., sep = " ~ ") %>%
    as.formula %>%
    poisson2(reg_data, B = B)

  if(is.null(cname)) {
    cname <- as.character(str_length(controls))
  }

  if(is.null(stname)) {
    stname <- str_replace_all(loctime_feglm, fixed(" "), "")
  }

  results <-
    bind_rows(mutate(nofirms,
                     outcome = "log E[Output]",
                     estimator = "poisson",
                     firms = "No",
                     controls = cname,
                     spacetime = stname),
              mutate(firms,
                     outcome = "log E[Output]",
                     estimator = "poisson",
                     firms = "Yes",
                     controls = cname,
                     spacetime = stname))
  return(results)
}


# omnibus routine that we can pass a high level model descriptor to for use in
# a larger SLURM batch file
run_models <- function(estimator, control_spec) {
  print(paste("estimator is:", estimator))
  print(paste("control spec is:", control_spec))

  controls <-
    case_when(
      control_spec == 1 ~ list(base_controls, "Base"),
      control_spec == 2 ~ list(list(base_controls, shape_controls),
                               "Shape"),
      control_spec == 3 ~ list(list(base_controls,
                                    shape_controls,
                                    distance_controls),
                               "Distance"),
      control_spec == 4 ~ list(list(base_controls,
                                    shape_controls,
                                    distance_controls,
                                    landcover_controls),
                               "Landcover"))

  if(estimator == "fe") {
    results <- fe_models(controls[[1]],
                         cname = controls[[2]],
                         stname = paste(gridsize, "Yr+YQ", sep = ""))
  } else if(estimator == "twopart") {
    results <- twopart_models(controls[[1]],
                              cname = controls[[2]],
                              stname = "CountyYr+YQ")
  } else if(estimator == "poisson") {
    results <- poisson2_models(controls[[1]],
                               cname = controls[[2]],
                               stname = "CountyYr+YQ")
  } 
  return(results)
}

#===========================================================================
# run models and make some tables
#===========================================================================

# if running in batch mode, estimate the model/spec requested
if(!is.na(estimator) & !is.na(control_spec)) {
  results <- run_models(estimator, control_spec)
  print(results)
  
  filename <-
    paste("output_effects", estimator, control_spec, sep = "_") %>%
    paste(".Rda", sep = "")
  
  save(results, file = file.path(edir, filename))
}




