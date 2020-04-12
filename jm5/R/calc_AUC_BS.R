check_exist <- function(variable){
  for (i in 1:length(variable)){
    if (!exists(variable[i])){
      err_message <- paste0('Variable ', l.variable, ' does not exist!')
      stop(err_message)
    }
  }
}

filter_dat <- function(long.test, surv.test, Tstart){
  surv.ID <- surv.test[surv.test$Time>=Tstart, ]$ID
  surv.filtered <- surv.test[surv.test$ID %in% surv.ID, ]
  long.filtered <- long.test[long.test$ID %in% surv.ID, ]
  l <- list(long.filtered = long.filtered, surv.filtered = surv.filtered)
  return(l)
}

log.lik <- function(yi, ti, xi, Ti, di, wi, sigma_u, sigma_e, alpha, ut, fixed.eff){
  ll <- sum(dnorm(yi, fixed.eff$long.mui.nore + ut, sigma_e, log = T)) + 
    di*(fixed.eff$loghi.nore + alpha*ut) + fixed.eff$logSi.nore*exp(alpha*ut) + 
    dnorm(ut, 0, sigma_u, log = T)
  return(ll)
}

MH.rw <- function(yi, ti, xi, Ti, di, wi, sigma_u, sigma_e, alpha, fixed.eff, M, burnin){
  u <- rep(NA, M + burnin)
  u[1] <- rnorm(1, sd = sigma_u)
  for (m in 1:(M + burnin - 1)){
    u[m+1] <- u[m] + runif(1, -1, 1)
    logR.numerator <- log.lik(yi, ti, xi, Ti, di, wi, sigma_u, sigma_e, alpha, u[m+1], fixed.eff)
    logR.denominator <- log.lik(yi, ti, xi, Ti, di, wi, sigma_u, sigma_e, alpha, u[m], fixed.eff)
    R <- exp(logR.numerator - logR.denominator)
    keep <- rbinom(1, 1, min(R, 1))
    if (keep!=1){
      u[m+1] <- u[m]
    }
  }
  return(u[(burnin+1):(M+burnin)])
}

calc_AUC_BS <- function(l.formula, s.formula, long.test, surv.test, 
                        Tstart, Tstop, summary.mean){
  
  ## Check existence of variable
  attach(long.test)
  l.variable <- c(all.vars(l.formula), 'ID')
  p.long <- length(all.vars(l.formula))
  check_exist(l.variable)
  detach(long.test)
  
  attach(surv.test)
  s.variable <- c(all.vars(s.formula), 'ID')
  p.surv <- length(all.vars(s.formula)) - 2
  check_exist(s.variable)
  detach(surv.test)
  
  ## define MH sampler ##
  M <- 2000
  burnin <- 5000
  
  ## get estimated parameters ##
  beta0 <- summary.mean[1]
  beta1 <- summary.mean[2]
  beta <- summary.mean[3:p.long]
  sigma_u <- summary.mean[p.long+1]
  sigma_e <- summary.mean[p.long+2]
  logh0 <- summary.mean[p.long+3]
  gamma <- summary.mean[(p.long+4):(p.long+p.surv+3)]
  alpha <- summary.mean[length(summary.mean)]
  
  ## filter the longitudinal data and survival data who survives up to Tstart
  l <- filter_dat(long.test, surv.test, Tstart)
  long.filtered <- l[[1]]
  surv.filtered <- l[[2]]
  ID.unique <- surv.filtered$ID
  surv.prob.est <- rep(NA, length(ID.unique))
  ID.index <- 0
  for (test.ID in ID.unique){
    ID.index <- ID.index + 1
    tmp.long <- long.filtered[long.filtered$ID==test.ID & long.filtered$obstime<=Tstart, ]
    tmp.surv <- surv.filtered[surv.filtered$ID==test.ID, ]
    
    ## MH.rw sample the random effect u
    yi <- tmp.long[, 1]
    ti <- tmp.long[, 2]
    xi <- as.matrix(tmp.long[, 3:p.long])
    long.const.mui <- beta0 + (xi %*% beta)[, 1]
    long.mui.nore <- long.const.mui + beta1*ti
    
    Ti <- Tstart
    di <- 0
    wi <- as.matrix(tmp.surv[, 3:(p.surv+2)])
    loghi.nore <- logh0 + (wi %*% gamma)[, 1] + alpha*(long.const.mui[1] + beta1*Ti)
    logSi.nore <- -exp(logh0 + (wi %*% gamma)[, 1] + alpha*long.const.mui[1])*
      (exp(alpha*beta1*Ti) - 1)/(alpha*beta1)
    fixed.eff <- list(long.mui.nore = long.mui.nore,
                      loghi.nore = loghi.nore,
                      logSi.nore = logSi.nore)
    
    u <- MH.rw(yi, ti, xi, Ti, di, wi, sigma_u, sigma_e, alpha, fixed.eff, M, burnin)
    logSi.Tstart <- -exp(logh0 + (wi %*% gamma)[, 1] + alpha*long.const.mui[1])*
      (exp(alpha*beta1*Tstart) - 1)/(alpha*beta1)*exp(alpha*u)
    logSi.Tstop <- -exp(logh0 + (wi %*% gamma)[, 1] + alpha*long.const.mui[1])*
      (exp(alpha*beta1*Tstop) - 1)/(alpha*beta1)*exp(alpha*u)
    Si <- exp(logSi.Tstop - logSi.Tstart)
    surv.prob.est[ID.index] <- mean(Si)
  }
  
  
}