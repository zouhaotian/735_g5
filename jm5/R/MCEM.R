check_exist <- function(variable){
  for (i in 1:length(variable)){
    if (!exists(variable[i])){
      err_message <- paste0('Variable ', l.variable, ' does not exist!')
      stop(err_message)
    }
  }
}

log.lik <- function(yi, ti, xi, Ti, di, wi, theta, ut, fixed.eff){
  ll <- sum(dnorm(yi, fixed.eff$long.mui.nore + ut, theta$sigma_e, log = T)) + 
    di*(fixed.eff$loghi.nore + theta$alpha*ut) + fixed.eff$logSi.nore*exp(theta$alpha*ut) + 
    dnorm(ut, 0, theta$sigma_u, log = T)
  return(ll)
}

MH.rw <- function(yi, ti, xi, Ti, di, wi, theta, fixed.eff, M, burnin){
  u <- rep(NA, M + burnin)
  u[1] <- rnorm(1, sd = theta$sigma_u)
  for (m in 1:(M + burnin - 1)){
    u[m+1] <- u[m] + runif(1, -1, 1)
    logR.numerator <- log.lik(yi, ti, xi, Ti, di, wi, theta, u[m+1], fixed.eff)
    logR.denominator <- log.lik(yi, ti, xi, Ti, di, wi, theta, u[m], fixed.eff)
    R <- exp(logR.numerator - logR.denominator)
    keep <- rbinom(1, 1, min(R, 1))
    if (keep!=1){
      u[m+1] <- u[m]
    }
  }
  return(u[(burnin+1):(M+burnin)])
}

tolerance <- function(theta.1, theta.2){
  p <- length(theta.1)
  max.tol <- 0
  for (i in 1:p){
    rel.tol <- abs(max(theta.1[[i]] - theta.2[[i]]))
    max.tol <- max(max.tol, rel.tol)
  }
  return(max.tol)
}


#' Joint modeling using MCEM
#'
#' This function allows you to pass in longitudinal and survival
#' formula, with the longitudinal and survival dataset, 
#' to perform EM algorithm using Metropolis Hastings random walk sampler 
#' 
#' @param l.formula the formula for longitudinal model
#' @param s.formula the formula for survival model
#' @param long.dat the longitudinal dataset
#' @param surv.dat the survival dataset
#' @param max.iter maximum iterations allowed (default = 30)
#' @param tol tolerance (default = 1e-4)
#' @param seed random seed (default = 123)
#' 
#' @return a list of estimated parameters
#'  
#' @examples
#' 
#' MCEM(CD4 ~ obstime, Surv(Time, death) ~ drugddI, 
#' jm_filter(CD4 ~ obstime, Surv(Time, death) ~ drug, JM::aids, JM::aids$patient)[[1]],
#' jm_filter(CD4 ~ obstime, Surv(Time, death) ~ drug, JM::aids, JM::aids$patient)[[2]])
#' 
#' @importFrom nlme lme VarCorr
#' @importFrom survival Surv coxph survreg
#' 
#' @export

MCEM <- function(l.formula, s.formula, long.dat, surv.dat, 
                 max.iter = 30, tol = 1e-4, seed = 123){
  
  set.seed(seed)
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
  
  ## We obtain covariates and response from dataset ##
  Y <- long.dat[, 1]
  obstime <- long.dat[, 2]
  X <- as.matrix(long.dat[, 3:(ncol(long.dat)-1)])
  
  surv_time <- surv.dat[, 1]
  delta <- surv.dat[, 2]
  W <- as.matrix(surv.dat[, 3:(ncol(surv.dat)-1)])
  
  ## First, we initialize our parameters:
  ## longitudinal part
  fit.lme <- lme(l.formula, random = ~1|ID, data = long.dat)
  coef.lme <- unname(fit.lme$coefficients$fixed)
  beta0 <- coef.lme[1]
  beta1 <- coef.lme[2]
  beta <- coef.lme[3:p.long]
  
  random.effects <- as.vector(fit.lme$coefficients$random[[1]])
  Sigma <- as.numeric(VarCorr(fit.lme))
  sigma_u <- Sigma[2]
  sigma_e <- Sigma[4]
  
  ## survival part
  latent.dat <- as.matrix(long.dat[ID.end, 3:(ncol(long.dat)-1)])
  latent.dat <- cbind(rep(1, n), latent.dat)
  mu_hat <- as.vector((latent.dat %*% coef.lme[-2])[, 1] + random.effects)
  surv.dat2 <- cbind(surv.dat, mu_hat)
  s.formula2 <- update(s.formula, ~. + mu_hat)
  
  fit.surv <- coxph(s.formula2, method = 'breslow', data = surv.dat2)
  coef.surv <- unname(fit.surv$coefficients)
  gamma <- coef.surv[1:p.surv]
  alpha <- coef.surv[p.surv+1]
  
  inits <- list(beta0 = beta0, beta1 = beta1, beta = beta, 
                sigma_u = sigma_u, sigma_e = sigma_e,
                logh0 = 0, gamma = gamma, alpha = alpha)
  theta <- inits
  
  ## specify tuning parameters for EM algorithm
  iter <- 0 ## current iteration
  curr.tol <- Inf ## current tolerance
  M <- 5000 ## total 5000 samples
  burnin <- 2000 ## with 2000 burnin period
  
  ## augment data ##
  long.dat.aug <- matrix(NA, nobs*M, ncol(long.dat))
  surv.dat.aug <- matrix(NA, n*M, ncol(surv.dat))
  latent.dat.aug <- matrix(NA, n*M, ncol(latent.dat))
  for (i in 1:nobs){
    long.dat.aug[((i-1)*M+1):(i*M), ] <- matrix(rep(as.numeric(long.dat[i, ]), M), M, ncol(long.dat), byrow = T)
  }
  for (i in 1:n){
    surv.dat.aug[((i-1)*M+1):(i*M), ] <- matrix(rep(as.numeric(surv.dat[i, ]), M), M, ncol(surv.dat), byrow = T)
    latent.dat.aug[((i-1)*M+1):(i*M), ] <- matrix(rep(as.numeric(latent.dat[i, ]), M), M, ncol(latent.dat), byrow = T)
  }
  colnames(long.dat.aug) <- colnames(long.dat)
  colnames(surv.dat.aug) <- colnames(surv.dat)
  long.dat.aug <- data.frame(long.dat.aug)
  surv.dat.aug <- data.frame(surv.dat.aug)
  
  while (iter<max.iter & curr.tol>tol){
    iter <- iter + 1
    theta.old <- theta
    
    ID <- surv.dat$ID
    re.sample <- rep(NA, M*n)
    re.sample.aug <- rep(NA, M*nobs)
    start <- 0
    end <- 0
    for (i in ID){
      ## obtain the response and covariates for subject i
      yi <- Y[ID.start[i]:ID.end[i]]
      ti <- obstime[ID.start[i]:ID.end[i]]
      xi <- X[ID.start[i]:ID.end[i], ]
      long.const.mui <- theta.old$beta0 + (xi %*% theta.old$beta)[, 1]
      long.mui.nore <- long.const.mui + theta.old$beta1*ti
      
      Ti <- surv_time[i]
      di <- delta[i]
      wi <- W[i, ]
      loghi.nore <- theta$logh0 + (wi %*% theta.old$gamma)[, 1] + 
        theta.old$alpha*(long.const.mui[1] + theta.old$beta1*Ti)
      logSi.nore <- -exp(theta$logh0 + (wi %*% theta.old$gamma)[, 1] + theta.old$alpha*long.const.mui[1])*
        exp(theta$alpha*theta$beta1*Ti - 1)/(theta$alpha*theta$beta1*Ti)
      fixed.eff <- list(long.mui.nore = long.mui.nore,
                        loghi.nore = loghi.nore,
                        logSi.nore = logSi.nore)
      u <- MH.rw(yi, ti, xi, Ti, di, wi, theta.old, fixed.eff, M, burnin)
      re.sample[((i-1)*M+1):(i*M)] <- u
      start <- end + 1
      end <- end + count[i]*M
      re.sample.aug[start:end] <- rep(u, count[i]) 
    }
    
    long.dat.aug2 <- cbind(long.dat.aug, u = re.sample.aug)
    
    ## M-step ##
    lm.fit <- as.numeric(lm(l.formula, data = long.dat.aug2, offset = u)$coefficients)
    sigma_e_new <- sigma(lm(l.formula, data = long.dat.aug2, offset = u))
    sigma_u_new <- sd(re.sample)
    
    mu_hat_aug <- as.vector((latent.dat %*% lm.fit[-2])[, 1] + re.sample)
    surv.dat.aug2 <- cbind(surv.dat.aug, mu_hat_aug)
    s.formula3 <- update(s.formula, ~. + mu_hat_aug)
    surv.fit <- survreg(s.formula3, data = surv.dat.aug2, dist = 'exponential')
    surv.coef <- as.numeric(surv.fit$coefficients)
    
    theta$beta0 <- lm.fit[1]
    theta$beta1 <- lm.fit[2]
    theta$beta <- lm.fit[3:length(lm.fit)]
    theta$sigma_u <- sigma_u_new
    theta$sigma_e <- sigma_e_new
    theta$logh0 <- -surv.coef[1]
    theta$gamma <- -surv.coef[2:(1+p.surv)]
    theta$alpha <- -surv.coef[length(surv.coef)]
    
    ## Calculate tolerance
    curr.tol <- tolerance(theta, theta.old)
    cat(sprintf('Iter:%d curr.tol:%.6f\n', iter, curr.tol))
    if (iter==max.iter & curr.tol>tol){
      warning('Iteration reached maximum without convergence!')
    }
  }
  return(theta)
}
