#Data Simulation based on Propsal_v3
# We compare the model fitting results using the simulation set-up as follows:
  
# We simulate N=500 subjects, and time points as (0, 2, 6, 12, 18)
# as from the `aids' data in the longitudinal model, 
# and another 200 subjects as the testing dataset.

#Set Seed
set.seed(2020)

#number of individuals
n = 500

#number of observations
time = c(0, 2, 6, 12, 18)
p = length(time)

# First, we simulate the longitudinal responses using the following formula: 
# Yij = 20 + (-1) x tij + (-4) x xi + ui + eij , 
# where ui ~ N(0; 2); eij ~ N(0; 1),
# xi ~ binom(1; 0.4).

### set up data frame
dat = data.frame(subject=rep(1:n, times = p), time=rep(time,each=n))
#dat3 <- dat[order(dat$subject),]
ui = as.vector(rnorm(n,0,2))
eij = as.vector(rnorm(n*p,0,1))
xi = as.vector(rbinom(n, 1, 0.4))
dat1 = cbind(dat,rep(xi,p),rep(ui,p),eij)
dat2 = dat1[order(dat1$subject,dat1$time),]
colnames(dat2)[3] = "xi"
colnames(dat2)[4] = "ui"
dat2$Yij = 20-1*dat2$time -4*dat2$xi + dat2$ui + dat2$eij
head(dat2)

# Second, we simulate N random uniform distribution, which corresponds to the survival probability. 
# We calculate the failure time as: hi(Ti) = exp(-2) exp(-1*wi + 0:2 *mi(Ti)), 
# where mi(Ti) = 20 + (-1) * Ti + (-4) * xi + ui.
# And we can calculate the survival time using an explicit formula.

dat_tmp = dat2[ which(dat2$time==0), ]
dat_tmp$time <- NULL
dat_tmp$eij <- NULL
dat_tmp$Yij <- NULL
dat_tmp$N = runif(n,0,1) #survival probability
#define Ti and wi
dat_tmp$mi_ti = 20 + (-1)*dat_tmp$Ti + (-4)*dat_tmp$xi + dat_tmp$ui
dat_tmp$hi_ti = exp(-2)*exp(-1*dat_tmp$wi+0.2*dat_tmp$mi_ti) #actual failure time
head(dat_tmp,10)
