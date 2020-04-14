data{
  int n;
  real x[n];
  real y[n];
}
parameters{
  real beta[2];
  real <lower=0> sigma_e;
}
transformed parameters{
  real mu[n];
  for (i in 1:n){
    mu[i] = beta[1] + beta[2]*x[i];
  }
}
model{
  beta ~ normal(0, 10);
  sigma_e ~ inv_gamma(0.1, 0.1);
  
  y ~ normal(mu, sigma_e);
}
