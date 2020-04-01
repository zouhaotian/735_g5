data{
  int n;
  real w[n];
  real surv_time[n];
  real status[n];
}
parameters{
  real logh0;
  real gamma;
}
transformed parameters{
  real h[n];
  real S[n];
  real LL[n];
  for (i in 1:n){
    h[i] = exp(logh0 + w[i]*gamma);
    S[i] = exp(-exp(logh0 + w[i]*gamma)*surv_time[i]);
    LL[i] = log(pow(h[i], status[i])*S[i]);
  }
}
model{
  logh0 ~ normal(0, 10);
  gamma ~ normal(0, 10);
 
  target+=LL;
}
