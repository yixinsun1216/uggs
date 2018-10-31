library(tidyverse)
library(furrr)
library(alpaca)

if(Sys.getenv("USERNAME") == "Yixin Sun"){
	root <- "C:/Users/Yixin Sun/Documents/Github/uggs"
	ddir <- "C:/Users/Yixin Sun/Documents/Dropbox/texas"
}

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

controls <- c("log(acres)")
loctime_feglm <- "CountyYr + YearQtr"

rhs <-
  paste("Auction", controls, sep = " + ") %>%
  paste(loctime_feglm, sep = " | ")


f <- 
  paste("drilled", rhs, sep = " ~ ") %>%
  as.formula

poisson_test <- function(formula, df){
  m <- feglm(formula, df, family = poisson())
  return(as.numeric(coef(m)["Auction"]))
}
