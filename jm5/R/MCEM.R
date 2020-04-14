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

log.lik.aug <- function(y.aug, t.aug, x.aug, T.aug, d.aug, w.aug, tmp.paras, p.long, p.surv, re.sample.aug, re.sample, latent.dat.aug){
  beta0 <- tmp.paras[1] 
  beta1 <- tmp.paras[2]
  beta <- tmp.paras[3:p.long]
  sigma_e <- tmp.paras[p.long+2]
  logh0 <- tmp.paras[p.long+3]
  gamma <- tmp.paras[(p.long+4):(p.long+p.surv+3)] 
  alpha <- tmp.paras[length(tmp.paras)]
  
  long.mu.aug <- beta0 + beta1*t.aug + (x.aug %*% beta)[, 1] + re.sample.aug
  logh.aug <- logh0 + (w.aug %*% gamma)[, 1] + 
    alpha*((latent.dat.aug %*% c(beta0, beta))[, 1] + beta1*T.aug + re.sample)
  logS.aug <- -exp(logh0 + (w.aug %*% gamma)[, 1] + alpha*((latent.dat.aug %*% c(beta0, beta))[, 1] + re.sample))*
    (exp(alpha*beta1*T.aug) - 1)/(alpha*beta1)
  
  ll <- sum(dnorm(y.aug, long.mu.aug, sigma_e, log = T)) + 
    sum(d.aug*logh.aug + logS.aug)
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
#' @param tol tolerance (default = 1e-3)
#' @param seed random seed (default = 123)
#' @param progress whether to show progress for each iteration (default = TRUE)
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
                 max.iter = 30, tol = 1e-3, seed = 123, progress = TRUE){
  
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
  
  fit.surv <- survreg(s.formula2, data = surv.dat2, dist = 'exponential')
  coef.surv <- unname(fit.surv$coefficients)
  logh0 <- -coef.surv[1]
  gamma <- -coef.surv[2:(p.surv+1)]
  alpha <- -coef.surv[length(coef.surv)]
  
  inits <- list(beta0 = beta0, beta1 = beta1, beta = beta, 
                sigma_u = sigma_u, sigma_e = sigma_e,
                logh0 = logh0, gamma = gamma, alpha = alpha)
  theta <- inits
  
  ## specify tuning parameters for EM algorithm
  iter <- 0 ## current iteration
  curr.tol <- 1 ## current tolerance
  M <- 100 ## total random effect samples
  M.paras <- 1000 ## total parameter samples
  burnin <- 1000 ## with burnin period
  
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
  
  y.aug <- long.dat.aug[, 1]
  t.aug <- long.dat.aug[, 2]
  x.aug <- as.matrix(long.dat.aug[, 3:(ncol(long.dat.aug)-1)])
  
  T.aug <- surv.dat.aug[, 1]
  d.aug <- surv.dat.aug[, 2]
  w.aug <- as.matrix(surv.dat.aug[, 3:(ncol(surv.dat.aug)-1)])
  
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
        (exp(theta$alpha*theta$beta1*Ti) - 1)/(theta$alpha*theta$beta1)
      fixed.eff <- list(long.mui.nore = long.mui.nore,
                        loghi.nore = loghi.nore,
                        logSi.nore = logSi.nore)
      u <- MH.rw(yi, ti, xi, Ti, di, wi, theta.old, fixed.eff, M, burnin)
      re.sample[((i-1)*M+1):(i*M)] <- u
      start <- end + 1
      end <- end + count[i]*M
      re.sample.aug[start:end] <- rep(u, count[i]) 
    }
    
    ## M-step: we use MH.rw rather than directly maximizing the likelihood
    ## paras.est: beta0, beta1, beta, sigma_u, sigma_e, logh0
    surv.paras <- as.numeric(unlist(theta))[-(1:(p.long+3))]
    paras.est <- matrix(NA, (M.paras+burnin), p.long+3)
    paras.est[1, ] <- as.numeric(unlist(theta))[1:(p.long+3)]
    paras.est[, p.long+1] <- sd(re.sample)
    long.dat.aug2 <- cbind(long.dat.aug, re.sample.aug)
    paras.est[, p.long+2] <- sigma(lm(l.formula, data=long.dat.aug2, offset = re.sample.aug))
    
    for (m in 1:(M.paras+burnin-1)){
      paras.est[m+1, -((p.long+1):(p.long+2))] <- paras.est[m, -((p.long+1):(p.long+2))] + 0.1*curr.tol*runif(1, -1, 1)
      full.paras.est.new <- c(paras.est[m+1, ], surv.paras)
      full.paras.est.old <- c(paras.est[m, ], surv.paras)
      logR.numerator <- log.lik.aug(y.aug, t.aug, x.aug, T.aug, d.aug, w.aug, full.paras.est.new, p.long, p.surv, re.sample.aug, re.sample, latent.dat.aug)
      logR.denominator <- log.lik.aug(y.aug, t.aug, x.aug, T.aug, d.aug, w.aug, full.paras.est.old, p.long, p.surv, re.sample.aug, re.sample, latent.dat.aug)
      R <- exp(logR.numerator - logR.denominator)
      keep <- rbinom(1, 1, min(1, R))
      if (keep!=1){
        paras.est[m+1, ] <- paras.est[m, ]
      }
    }
    
    theta$beta0 <- mean(paras.est[(M.paras+1):(M.paras+burnin), 1])
    theta$beta1 <- mean(paras.est[(M.paras+1):(M.paras+burnin), 2])
    theta$beta <- colMeans(paras.est[(M.paras+1):(M.paras+burnin), 3:p.long])
    theta$sigma_u <- mean(paras.est[(M.paras+1):(M.paras+burnin), p.long + 1])
    theta$sigma_e <- mean(paras.est[(M.paras+1):(M.paras+burnin), p.long + 2])
    theta$logh0 <- mean(paras.est[(M.paras+1):(M.paras+burnin), p.long + 3])
    
    ## maximize survival part
    mu_hat <- (latent.dat.aug %*% c(theta$beta0, theta$beta))[, 1] + re.sample
    surv.dat.aug2 <- cbind(surv.dat.aug, mu_hat)
    s.formula2 <- update(s.formula, ~. + mu_hat)
    fit.surv <- coxph(s.formula2, data = surv.dat.aug2, method = 'breslow')
    coef.surv <- unname(fit.surv$coefficients)
    theta$gamma <- coef.surv[1:p.surv]
    theta$alpha <- coef.surv[p.surv+1]

    ## Calculate tolerance
    curr.tol <- tolerance(theta, theta.old)
    if (progress==TRUE){
      cat(sprintf('Iter:%d curr.tol:%.6f\n', iter, curr.tol)) 
    }
    if (iter==max.iter){
      warning('Iteration reached maximum without convergence!')
    }
  }
  k <- unlist(theta)
  return(k)
}
