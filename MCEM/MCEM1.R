library(data.table)
library(JM)
data(aids)
aids = data.table(aids)

# xi: drug, gender, prevOI, AZT
# wi: drug
# ti: obstime
# Ti: Time
# yi: CD4
# deltai: death

# latent mean
m = function(ti, xi, ui, beta0, beta1, beta){
  # beta0: intercept
  # beta1: coefficient of t_i
  # beta: coeffieicent of x_i
  return(beta0 + ti * beta1 + xi %*% beta+ ui)
}

# log hazard function
logh = function(wi, Ti, xi, ui, h0, beta0, beta1, beta, gamma, alpha) {
  # h0: baseline hazard function
  # gamma: coefficient of wi
  # alpha: coefficient of latent mean
  return(log(h0) + wi * gamma + alpha * m(Ti, xi[1, ], ui, beta0, beta1, beta))
}

# log survival function
logS = function(wi, Ti, xi, ui, h0, beta0, beta1, beta, gamma, alpha) {
  return(-h0 * exp(wi * gamma + alpha * (beta0 + xi[1, ] %*% beta + ui)) * (exp(alpha * beta1 * Ti) - 1) / (alpha * beta1))
}

# complete data log likelihood for subject i or the log posterior distribution of ui
lli = function(ui, yi, wi, Ti, xi, ti, deltai, h0, beta0, beta1, beta, gamma, alpha, s2e, s2u) {
  ll = sum(dnorm(yi, mean = m(ti, xi, ui, beta0, beta1, beta), sd = sqrt(s2e), log = TRUE)) + 
    deltai * logh(wi, Ti, xi, ui, h0, beta0, beta1, beta, gamma, alpha) +
    logS(wi, Ti, xi, ui, h0, beta0, beta1, beta, gamma, alpha) + 
    dnorm(ui, mean = 0, sd = sqrt(s2u), log = TRUE)
  return(ll)
}

# MH ratio
R = function(u, ut, yi, wi, Ti, xi, ti, deltai, h0, beta0, beta1, beta, gamma, alpha, s2e, s2u){
  # u: proposal
  # ut: current state
  logR = lli(u, yi, wi, Ti, xi, ti, deltai, h0, beta0, beta1, beta, gamma, alpha, s2e, s2u) - 
    lli(ut, yi, wi, Ti, xi, ti, deltai, h0, beta0, beta1, beta, gamma, alpha, s2e, s2u)
  
  return(exp(logR))
}

# mcmc sampler for ui
mcmc.sampler.i = function(yi, wi, Ti, xi, ti, deltai, h0, beta0, beta1, beta, gamma, alpha, s2e, s2u, M = 10000, burn_in = 2000, sigma = 1, trace = 0){
  # sigma is the standard deviation of normal distribution proposal
  u.rw.chain = rep(0, M)
  accept = 0
  for (i in 1:(M-1)) {
    ut = u.rw.chain[i]
    u = ut + rnorm(1, sd = sigma)
    r = R(u, ut, yi, wi, Ti, xi, ti, deltai, h0, beta0, beta1, beta, gamma, alpha, s2e, s2u)
    r = min(r, 1)
    keep = rbinom(1, 1, r)
    if (keep == 1){
      u.rw.chain[i + 1] = u
      accept = accept + 1
    }else{
      u.rw.chain[i + 1] = ut
    }
  }
  
  if (trace != 0) {
    cat(sprintf("Acceptance rate: %.3f", accept / M))
  }
  return(u.rw.chain[-(1:burn_in)])
}



# data augmentation
augment = function(data, M = 10000){
  if(!is.data.table(data)) data = data.table(data)
  data_aug = data[rep(1:nrow(data),  M),]
  data_aug$M = rep(1:M, each = nrow(data))
  data_aug = data_aug[order(patient, M),]
  data_aug$samples = rep(0, nrow(data_aug))
  return(data_aug)
}

mcmc.sampler.all = function(data, data_aug, h0, beta0, beta1, beta, gamma, alpha, s2e, s2u, M = 10000, burn_in = 2000, sigma = 1) {
  n = length(unique(data$patient))
  for (i in 1:n) {
    yi = data[patient == i, CD4]
    xi = as.matrix(data[patient == i, .(drug, gender, prevOI, AZT)])
    ti = data[patient == i, obstime]
    deltai = data[patient == i, death][1]
    Ti = data[patient == i, Time][1]
    wi = data[patient == i, drug][1]
    sample_i = mcmc.sampler.i(yi, wi, Ti, xi, ti, deltai, h0, beta0, beta1, beta, gamma, alpha, s2e, s2u, M + burn_in, burn_in, sigma)
    sample_i = sample_i[rep(1:M, each = length(yi))]
    data_aug[patient == i, samples := sample_i]
  }
  return(data_aug)
}


aids$drug = as.numeric(aids$drug) - 1
aids$gender = as.numeric(aids$gender) - 1
aids$prevOI = as.numeric(aids$prevOI) - 1
aids$AZT = as.numeric(aids$AZT) - 1
aids$death = as.numeric(aids$death)

data_aug = augment(aids, M=10000)

h0=1
beta0 = 0.5
beta1 = 1.5
beta = c(0.6,0.7,0.8,0.9)
gamma = 0.5
alpha = 0.8
s2e = 1
s2u = 1

aids_M_aug = mcmc.sampler.all(aids, data_aug, h0, beta0, beta1, beta, gamma, alpha, s2e, s2u, M = 10000, burn_in = 2000, sigma = 1)
X_aug = model.matrix(~ 1 + obstime + drug + gender + prevOI + AZT, data = aids_M_aug)
X_aug_surv = model.matrix(~ 1 + Time + drug + gender + prevOI + AZT, data = aids_M_aug)

m = X_aug %*% c(beta0, beta1, beta) + aids_M_aug[, samples]
m_surv = X_aug_surv %*% c(beta0, beta1, beta) + aids_M_aug[, samples]
log_h =  log(h0) + aids_M_aug$drug * gamma + alpha * m_surv
log_s = -exp(log_h) * (exp(alpha * beta1 * aids_M_aug$Time) - 1) / (alpha * beta1)
qfunction = sum(dnorm(aids_M_aug$CD4, mean = m, sd = sqrt(s2e), log = TRUE))/10000 + 
  sum(aids_M_aug$death * log_h)/30000 + sum(log_s)/30000 + 
  sum(dnorm(aids_M_aug$samples, mean = 0, sd = sqrt(s2u), log = TRUE))/30000

  
  

