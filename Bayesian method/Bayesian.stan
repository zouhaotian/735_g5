data{
  int<lower=0> n;
  int<lower=0> nobs;
  int<lower=0> ncov;
  int<lower=0> ID[nobs];
  int<lower=0> J[n+1];
  real<lower=0> time[nobs];
  real surv_time[n];
  real status[n];
  real w[n];
  matrix[nobs, ncov] x;
  real Y[nobs];
}

parameters{
  real beta[2];
  vector[ncov] boldbeta;
  real logh0;
  real gamma;
  real alpha;
  real <lower=0> sigma_u;
  real <lower=0> sigma_e;
  
  real u[n]; //random effects
}

transformed parameters{
  vector[nobs] m;
  vector[nobs] mT;
  real h[n];
  real S[n];
  real LL[n];
  
  for (i in 1:nobs){
    m[i] = beta[1] + beta[2]*time[i] + u[ID[i]];
    mT[i] = beta[1] + beta[2]*surv_time[ID[i]] + u[ID[i]];
  }
  m = m + x*boldbeta;
  mT = mT + x*boldbeta;
  
  for (i in 1:n){
    h[i] = exp(logh0 + w[i]*gamma + alpha*mT[ID[i]]);
    S[i] = exp(-h[i]*(1-exp(-alpha*beta[2]*surv_time[i]))/(alpha*beta[2]));
    LL[i] = log(pow(h[i], status[i])*S[i]) + normal_lpdf(u[i] | 0, sigma_u);
    
    for(j in J[i]:(J[i+1]-1)){
      LL[i] = LL[i] + normal_lpdf(Y[j] | m[j], sigma_e);
    }
  }
}

model{
  beta ~ normal(0, 10);
  boldbeta ~ normal(0, 10);
  gamma ~ normal(0, 10);
  alpha ~ normal(0, 10);
  logh0 ~ normal(0, 10);
  sigma_u ~ inv_gamma(0.1, 0.1);
  sigma_e ~ inv_gamma(0.1, 0.1);
  u ~ normal(0, sigma_u);
  
  target += LL;
}
