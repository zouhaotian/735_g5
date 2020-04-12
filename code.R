library(knitr)
library(devtools)
load_all('jm5')

l.formula <- CD4 ~ obstime + drug + prevOI + AZT
s.formula <- Surv(Time, death) ~ drug
dat <- JM::aids
l <- jm_filter(l.formula, s.formula, dat, dat$patient)
long.dat <- l[[1]]
surv.dat <- l[[2]]

## summary statistics (need to edit)

l.formula2 <- CD4 ~ obstime + drugddI + prevOIAIDS + AZTfailure
s.formula2 <- Surv(Time, death) ~ drugddI

mcem.result <- MCEM(l.formula2, s.formula2, long.dat, surv.dat, max.iter = 10)

fit.result <- fit_stan(l.formula2, s.formula2, long.dat, surv.dat)
summary.mean <- summary(fit.result)$summary[, 1]
summary.mean <- summary.mean[1:(length(summary.mean) - 1)]

#save.image('1.RData')
cv.dat <- create_cv(long.dat, surv.dat)
st <- c(6, 12) ## starting time
dt <- c(2, 4, 6) ## prediction time window
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
}