check_exist <- function(variable){
  for (i in 1:length(variable)){
    if (!exists(variable[i])){
      err_message <- paste0('Variable ', l.variable, ' does not exist!')
      stop(err_message)
    }
  }
}

#' Joint modeling using Stan
#'
#' This function allows you to pass in longitudinal and survival
#' formula, with the longitudinal and survival dataset, 
#' to perform NUTS sampling using Stan
#' 
#' @param l.formula the formula for longitudinal model
#' @param s.formula the formula for survival model
#' @param long.dat the longitudinal dataset
#' @param surv.dat the survival dataset
#' @param n.iter total number of iterations in NUTS (default = 2000)
#' @param n.burnin number of burnins in NUTS (default = n.iter/2)
#' @param seed random seed (default = 123)
#' @param progress whether to show progress in sampling (default = TRUE)
#' 
#' @return a Stan fit object
#'  
#' @examples
#' 
#' fit_stan(CD4 ~ obstime, Surv(Time, death) ~ drugddI, 
#' jm_filter(CD4 ~ obstime, Surv(Time, death) ~ drug, JM::aids, JM::aids$patient)[[1]],
#' jm_filter(CD4 ~ obstime, Surv(Time, death) ~ drug, JM::aids, JM::aids$patient)[[2]])
#' 
#' @importFrom nlme lme VarCorr
#' @importFrom survival Surv survreg
#' @importFrom rstan rstan_options stan_model sampling 
#' 
#' @export
fit_stan <- function(l.formula, s.formula, long.dat, surv.dat,
                     n.iter = 2000, 
                     n.burnin = floor(n.iter/2),
                     seed = 123,
                     progress = TRUE){
  options(mc.cores = parallel::detectCores())
  rstan_options(auto_write = TRUE)
  
  ## Check existence of variable
  attach(long.dat)
  l.variable <- c(all.vars(l.formula), 'ID')
  check_exist(l.variable)
  detach(long.dat)
  
  attach(surv.dat)
  s.variable <- c(all.vars(s.formula), 'ID')
  check_exist(s.variable)
  detach(surv.dat)
  
  ## first, we determine number of subjects, number of covariates
  n <- nrow(surv.dat) ## number of subjects
  nobs <- nrow(long.dat) ## number of observations
  count <- as.vector(table(long.dat$ID)) ## number of observation each subject
  index <- 1:nobs
  
  ID.start <- rep(NA, n) ## starting index of each subject
  ID.end <- rep(NA, n) ## ending index of each subject
  attach(long.dat)
  for (i in unique(ID)){
    tmp.index <- index[which(ID==i)]
    ID.start[i] <- min(tmp.index)
    ID.end[i] <- max(tmp.index)
  }
  detach(long.dat)
  
  p.long <- ncol(long.dat) - 1 ## number of covariates for longitudinal
  p.surv <- ncol(surv.dat) - 3 ## number of covariates for survival
  
  ## First, we initialize our parameters:
  ## longitudinal part
  fit.lme <- lme(l.formula, random = ~1|ID, data = long.dat)
  coef.lme <- unname(fit.lme$coefficients$fixed)
  random.effects <- as.vector(fit.lme$coefficients$random[[1]])
  Sigma <- as.numeric(VarCorr(fit.lme))
  sigma_u_hat <- Sigma[2]
  sigma_e_hat <- Sigma[4]
  
  ## survival part
  latent.dat <- as.matrix(long.dat[ID.end, 3:(ncol(long.dat)-1)])
  latent.dat <- cbind(rep(1, n), latent.dat)
  mu_hat <- as.vector((latent.dat %*% coef.lme[-2])[, 1] + random.effects)
  surv.dat2 <- cbind(surv.dat, mu_hat)
  s.formula2 <- update(s.formula, ~. + mu_hat)
  
  fit.surv <- survreg(s.formula2, data = surv.dat2, dist = 'exponential')
  coef.surv <- unname(fit.surv$coefficients)
  logh0_hat <- -coef.surv[1]
  gamma_hat <- -coef.surv[2:(p.surv+1)]
  alpha_hat <- -coef.surv[length(coef.surv)]
  
  ## Second, we put all required data into a list, and define parameters
  stanDat <- list(n = n, nobs = nobs, P1 = p.long, P2 = p.surv,
                  ID = as.numeric(long.dat$ID), 
                  ID_start = ID.start, ID_end = ID.end,
                  X = as.matrix(long.dat[, 2:(ncol(long.dat) - 1)]),
                  Y = long.dat[, 1],
                  W = as.matrix(surv.dat[, 3:(ncol(surv.dat) - 1)]),
                  TIME = surv.dat[, 1],
                  STATUS = surv.dat[, 2])
  pars <- c('beta', 'sigma_u', 'sigma_e','logh0', 'gamma', 'alpha')
  inits <- list(beta = coef.lme,
                sigma_u = sigma_u_hat,
                sigma_e = sigma_e_hat,
                logh0 = logh0_hat,
                gamma = as.array(gamma_hat), 
                alpha = alpha_hat)
  inits <- list(c1 = inits, c2 = inits, c3 = inits, c4 = inits)
  
  ## Third, we perform sampling
  model_code <- "
    data{
      int n;
      int nobs;
      int P1;
      int P2;
      int ID[nobs];
      int ID_start[n];
      int ID_end[n];
      
      matrix[nobs, P1-1] X;
      real Y[nobs];
      
      matrix[n, P2] W;
      real TIME[n];
      real STATUS[n];
    }
    parameters{
      vector[P1] beta;
      real <lower=0> sigma_u;
      real <lower=0> sigma_e;
      
      real logh0;
      vector[P2] gamma;
      real alpha;
      
      real u[n];
    }
    transformed parameters{
      real mu[nobs];
      real lp_long[nobs];
      real lp_surv[n];
      real h[n];
      real S[n];
      real LL[n];
      
      for (i in 1:nobs){
        lp_long[i] = beta[1] + X[i, 2:(P1-1)]*beta[3:P1] + u[ID[i]];
        mu[i] = lp_long[i] + beta[2]*X[i, 1];
      }
      for (i in 1:n){
        lp_surv[i] = W[i]*gamma;
        h[i] = exp(logh0 + lp_surv[i] + alpha*(lp_long[ID_end[i]] + beta[2]*TIME[i]));
        S[i] = exp(-exp(logh0 + lp_surv[i] + alpha*lp_long[ID_end[i]])*
                   (exp(alpha*beta[2]*TIME[i]) - 1)/(alpha*beta[2]));
        LL[i] = log(pow(h[i], STATUS[i])*S[i]);
      }
    }
    model{
      beta ~ normal(0, 10);
      logh0 ~ normal(0, 10);
      gamma ~ normal(0, 10);
      alpha ~ normal(0, 10);
      
      sigma_u ~ inv_gamma(0.1, 0.1);
      sigma_e ~ inv_gamma(0.1, 0.1);
      u ~ normal(0, sigma_u);
      Y ~ normal(mu, sigma_e);
      target+=LL;
    }
  "
  md = stan_model(model_code = model_code)
  fitStan <- sampling(md, data = stanDat, pars = pars, init = inits,
                      chains = 4, iter = n.iter, warmup = n.burnin,
                      seed = seed, open_progress = progress)
  return(fitStan)
}
