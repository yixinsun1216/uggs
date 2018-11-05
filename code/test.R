library(tidyverse)
library(furrr)
library(alpaca)
library(microbenchmark)
library(bcaboot)

if(Sys.getenv("USERNAME") == "Yixin Sun"){
	root <- "C:/Users/Yixin Sun/Documents/Github/uggs"
	ddir <- "C:/Users/Yixin Sun/Documents/Dropbox/texas"
}
setwd(root)

source(file.path(root, "code/uggs.R"))

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


f <- as.formula("drilled ~ Auction + log(acres) | CountyYr + YearQtr")
poisson_test <- function(formula, df){
  m <- feglm(formula, df, family = poisson())
  return(as.numeric(coef(m)["Auction"]))
}

check <- uggs(reg_data, 2000, f, poisson_test, jcount = 80, jreps = 5,
	alpha = c(0.025, 0.05, 0.1, .16), progress = TRUE)


# make data into matrix and use in bcajack function
reg_matrix <-
  reg_data %>%
  select(drilled, Auction, acres, CountyYr, YearQtr) %>%
  as.matrix

poisson_test_og <- function(x){
	f <- "drilled ~ Auction + log(acres) | CountyYr + YearQtr"
	df <- 
	  as_tibble(x) %>%
	  setNames(c("drilled", "Auction", "acres", "CountyYr", "YearQtr"))
	m <- feglm(as.formula(f), df, family = poisson())
  	return(as.numeric(coef(m)["Auction"]))
}


print(system.time(bca_check <- bcajack(reg_matrix, 2000, poisson_test_og, m = 80)))
print(system.time(check <- uggs(reg_data, 2000, f, poisson_test, jcount = 80)))
