library(knitr)
library(devtools)
load_all('jm5')

l.formula <- CD4 ~ obstime + drug + gender + prevOI + AZT
s.formula <- Surv(Time, death) ~ drug
dat <- JM::aids
l <- jm_filter(l.formula, s.formula, dat, dat$patient)
long.dat <- l[[1]]
surv.dat <- l[[2]]

kable(summary(l[[1]]))
kable(summary(l[[2]]))

l.formula2 <- CD4 ~ obstime + drugddI + gendermale + prevOIAIDS + AZTfailure
s.formula2 <- Surv(Time, death) ~ drugddI

mcem.result <- MCEM(l.formula2, s.formula2, long.dat, surv.dat)

fit.result <- fit_stan(l.formula2, s.formula2, long.dat, surv.dat)
summary.mean <- summary(fit.result)$summary[, 1]
summary.mean <- summary.mean[1:(length(summary.mean) - 1)]

#save.image('1.RData')
