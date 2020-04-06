setwd("/Users/jianchen/Documents/GitHub/735_g5/Bayesian method")
library(rstan)
library(JM)
library(tidyverse)
#options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

md <- stan_model('Bayesian.stan')

aids$drug = ifelse(aids$drug == "ddC", 1, 0)
aids$gender = ifelse(aids$gender == "male", 1, 0)
aids$prevOI = ifelse(aids$prevOI == "AIDS", 1, 0)
aids$AZT = ifelse(aids$AZT == "intolerance", 1, 0)

aids_short <- aids[1:9] %>%
  spread(key = obstime, value = CD4)


set.seed(1)
n <- length(unique(aids$patient))
nobs <- dim(aids)[1]
ID <- as.numeric(aids$patient)
ncov <- 4
J <- c(match(1:n, aids$patient), nobs+1)
time <- aids$obstime
surv_time <- aids_short$Time
status <- aids_short$death
w <- aids_short$drug
x <- as.matrix(aids[,6:9])
Y <- aids$CD4

## Specify stan data, parameters, and initial values for parameters
stanDat <- list(n = n, nobs = nobs, ID = ID, ncov = ncov, J = J, time = time,
                w = w, surv_time = surv_time, status = status, x = x, Y = Y)

## sampling from stan
fitStan <- sampling(md, data = stanDat)
fitStan

