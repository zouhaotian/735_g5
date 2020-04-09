library(data.table)

# latent mean
m = function(ti, xi, ui, beta0, beta1, beta){
  # beta0: intercept
  # beta1: coefficient of t_i
  # beta: coeffieicent of x_i
  return(beta0 + ti * beta1 + xi %*% beta)+ ui)
}

# log hazard function
logh = function(wi, ti, xi, ui, h0, beta0, beta1, beta, gamma, alpha) {
  # h0: baseline hazard function
  # gamma: coefficient of wi
  # alpha: coefficient of latent mean
  return(log(h0) + wi * gamma + alpha * m(ti, xi[1, ], ui, beta0, beta1, beta))
}

# log survival function
logS = function(wi, Ti, xi, ui, h0, beta0, beta1, beta, gamma, alpha) {
  # xi: covariates of subject i at Ti
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
  # u is the proposal, ut is the current state
  logR = lli(u, yi, wi, Ti, xi, ti, deltai, h0, beta0, beta1, beta, gamma, alpha, s2e, s2u) - 
    lli(ut, yi, wi, Ti, xi, ti, deltai, h0, beta0, beta1, beta, gamma, alpha, s2e, s2u)
  
  return(exp(logR))
}

# mcmc sampler for ui
mcmc.sampler.i = function(yi, wi, Ti, xi, ti, deltai, h0, beta0, beat1, beta, gamma, alpha, s2e, s2u, M = 10000, burn_in = 2000, sigma = 1){
  # sigma is the standard deviation of normal distribution proposal
  u.rw.chain = rep(0, M)
  for (i in 1:(M-1)) {
    ut = u.rw.chain[i]
    u = ut + rnorm(1, sd = sigma)
    r = R(u, ut, yi, wi, Ti, xi, ti, deltai, h0, beta0, beta1, beta, gamma, alpha, s2e, s2u)
    r = min(r, 1)
    keep = rbinom(1, 1, r)
    if (keep == 1){
      u.rw.chain[i + 1] = u
    }else{
      u.rw.chain[i + 1] = ut
    }
  }
  return(u.rw.chain[-(1:burn_in)])
}

augment = function(data, M = 10000){
  if(!is.data.table(data)) data = data.table(data)
  data_aug = data[rep(1:nrow(data),  M),]
  data_aug$M = rep(1:M, each = nrow(data))
  data_aug = data_aug[order(patient, M),]
  data_aug$samples = rep(0, nrow(data_aug))
  return(data_aug)
}

mcmc.sampler.all = function(data, h0, beta0, beta1, beta, gamma, alpha, s2e, s2u, M = 10000, burn_in = 2000, sigma = 1) {
  n = length(unique(data$patient))
  data_aug = augment(data, M)
  for (i in 1:n) {
    yi = data[patient == i, CD4]
    xi = as.matrix(data[patient == i, .(drug, gender, prevOI, AZT)])
    ti = data[patient == i, obstime]
    deltai = data[patient == i, death][1]
    Ti = data[patient == i, Time][1]
    wi = data[patient == i, drug][1]
    sample_i = mcmc.sampler.i(yi, wi, Ti, xi, ti, deltai, h0, beta0, beat1, beta, gamma, alpha, s2e, s2u, M + burn_in, burn_in, sigma)
    sample_i = sample_i[rep(1:M, each = length(yi))]
    data_aug[patient == i, samples := sample_i]
  }
  return(data_aug)
}

aids$drug = as.numeric(aids$drug) - 1
aids$gender = as.numeric(aids$gender) - 1
aids$prevOI = as.numeric(aids$prevOI) - 1
aids$AZT = as.numeric(aids$AZT) - 1
aids_aug = augment(aids, M = 1000)
h0 = exp(-2)
beta0 = 20
alpha = 0.2
beta1 = -1
beta = rep(0.1, 4)
gamma = -1
s2e = 1
s2u = 1
yi = aids[patient == 1, CD4]
xi = as.matrix(aids[patient == 1, .(drug, gender, prevOI, AZT)])
ti = aids[patient == 1, obstime]
deltai = aids[patient == 1, death][1]
Ti = aids[patient == 1, Time][1]
wi = aids[patient == 1, drug][1]
ui = 0.5
u = 0.9

chain = mcmc.sampler.i(yi, wi, Ti, xi, ti, deltai, h0, beta0, beat1, beta, gamma, alpha, s2e, s2u, M = 10000, burn_in = 2000, sigma = 1)

