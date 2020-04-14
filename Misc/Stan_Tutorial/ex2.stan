data{
  int n;
  int nobs;
  int ID[nobs];
  real time[nobs];
  real Y[nobs];
}
parameters{
  real beta[2];
  real <lower=0> sigma_u;
  real <lower=0> sigma_e;
  
  real u[n]; //random effects
}
transformed parameters{
  real mu[nobs];
  for (i in 1:nobs){
    mu[i] = beta[1] + beta[2]*time[i] + u[ID[i]];
  }
}
model{
  beta ~ normal(0, 10);
  sigma_u ~ inv_gamma(0.1, 0.1);
  sigma_e ~ inv_gamma(0.1, 0.1);
  
  u ~ normal(0, sigma_u);
  Y ~ normal(mu, sigma_e);
}
