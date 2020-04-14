library(tidyverse)
library(knitr)
library(devtools)
load_all('jm5')

## MCEM and Bayesian pre-processing 

l.formula <- CD4 ~ obstime + drug + prevOI + AZT
s.formula <- Surv(Time, death) ~ drug
dat <- JM::aids
l <- jm_filter(l.formula, s.formula, dat, dat$patient)
long.dat <- l[[1]]
surv.dat <- l[[2]]

## RSF pre-processing
nearestCD4 = long.dat %>%
  select(ID, obstime, CD4, drugddI) %>%
  group_by(ID) %>%
  filter(obstime==max(obstime)) %>%
  select(ID, CD4)

surv.dat2 = merge(surv.dat, nearestCD4, by="ID")

## summary statistics (need to edit)


## Model fit: MCEM; Bayesian;

l.formula2 <- CD4 ~ obstime + drugddI + prevOIAIDS + AZTfailure
s.formula2 <- Surv(Time, death) ~ drugddI

mcem.result <- MCEM(l.formula2, s.formula2, long.dat, surv.dat, max.iter = 10)
mcem.result

fit.result <- fit_stan(l.formula2, s.formula2, long.dat, surv.dat)
summary.mean <- summary(fit.result)$summary[, 1]
summary.mean <- summary.mean[1:(length(summary.mean) - 1)]
summary.mean

## model fit: RSF
s.formula2 <- Surv(Time, death) ~ drugddI + CD4
fit = RandomSurvivalForest(surv.dat2, s.formula2)

head(fit$predicted) ## predicted value for each patient
head(fit$survival[, 1:10]) ## In-bag survival function
plot.survival(fit) ## generates Survival plots
##

## Cross-validation to estimate accuracy of MCEM and Bayesian method
cv.dat <- create_cv(long.dat, surv.dat)
st <- c(2, 4, 6) ## starting time
dt <- c(0.5, 1) ## prediction time window
mcem.auc.bs.mean <- matrix(0, length(st)*length(dt), 2)
stan.auc.bs.mean <- matrix(0, length(st)*length(dt), 2)
for (i in 1:length(cv.dat)){
  ## training data and testing data
  long.train <- cv.dat[[i]][[1]]
  long.test <- cv.dat[[i]][[2]]
  surv.train <- cv.dat[[i]][[3]]
  surv.test <- cv.dat[[i]][[4]]
  
  ## run MCEM and Stan on training dataset
  mcem.train.mean <- MCEM(l.formula2, s.formula2, long.train, surv.train, max.iter = 10)
  stan.train <- fit_stan(l.formula2, s.formula2, long.train, surv.train)
  stan.train.mean <- summary(stan.train)$summary[, 1]
  stan.train.mean <- stan.train.mean[-length(stan.train.mean)]
  
  ## Calculate AUC and BS ##
  mcem.auc.bs <- matrix(NA, length(st)*length(dt), 2)
  stan.auc.bs <- matrix(NA, length(st)*length(dt), 2)
  index <- 0
  for (Tstart in st){
    for (Tdelta in dt){
      index <- index + 1
      mcem.auc.bs[index, ] <- calc_AUC_BS(l.formula2, s.formula2, 
                                          long.test, surv.test,
                                          Tstart, Tstart+Tdelta, 
                                          mcem.train.mean)
      stan.auc.bs[index, ] <- calc_AUC_BS(l.formula2, s.formula2, 
                                          long.test, surv.test,
                                          Tstart, Tstart+Tdelta, 
                                          stan.train.mean)
    }
  }
  mcem.auc.bs.mean <- mcem.auc.bs.mean + mcem.auc.bs
  stan.auc.bs.mean <- stan.auc.bs.mean + stan.auc.bs
}

mcem.auc.bs.mean/length(cv.dat)
stan.auc.bs.mean/length(cv.dat)

## Cross-validation: RSF
cv.dat <- create_cv(long.dat, surv.dat)
surv.train2 = left_join(cv.dat[[2]]$surv.train, nearestCD4, by="ID")
surv.test2 = left_join(cv.dat[[2]]$surv.test, nearestCD4, by="ID")

fit = RandomSurvivalForest(surv.train2, s.formula2)
res = predict(fit, newdata = surv.test2)
res$err.rate[!is.na(res$err.rate)] ## prediction error rate

## Simulation
sim.dat <- SimulateDataset(seed = 1)
