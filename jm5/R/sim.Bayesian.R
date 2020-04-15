get.estimate = function(result.sim){
  parname = names(result.sim)
  n = nrow(result.sim)
  estimate = colMeans(result.sim)
  SE = rep(NA, length(parname))
  Lower = rep(NA, length(parname))
  Upper = rep(NA, length(parname))
  for (i in 1:length(parname)) {
    SE[i] = sqrt(var(result.sim[,i]-estimate[i]))
    sort_par = sort(result.sim[,i])
    Lower[i] = sort_par[floor(0.025*n)]
    Upper[i] = sort_par[floor(0.975*n)]
  }
  return(remove_rownames(data.frame(Parameter = parname,
                                    Estimate = estimate,
                                    StdErr = SE,
                                    LowerLimit = Lower,
                                    UpperLimit = Upper)))
}

#' Simulation for Bayesian method
#' 
#' This function allows you to evaluate the correctness of Bayesian method
#' using simulated data
#' 
#' @param sim.num number of simulations (default = 10)
#' @param n.iter total number of iterations in NUTS (default = 2000)
#' @param n.burnin number of burnins in NUTS (default = n.iter/2)
#' @param seed random seed (default = 123)
#' 
#' @return a data frame containing parameter name, true values of parameters, 
#' mean of estimated parameters, SD of estimated parameters, square root of 
#' mean of SE^2 and coverage probability.
#' 
#' @example 
#' 
#' sim.Bayesian(sim.num = 2)
#' 
#' @export
sim.Bayesian = function(sim.num = 10,
                        n.iter = 2000, 
                        n.burnin = floor(n.iter/2),
                        seed = 123){
  
  l.formula.sim <- Y ~ obstime + x
  s.formula.sim <- Surv(Time, death) ~ w
  true_value = c(20, -1, -4, 2, 1, -2, -1, -0.2)
  estimate = matrix(NA, length(true_value), sim.num)
  SE = matrix(NA, length(true_value), sim.num)
  cover = matrix(NA, length(true_value), sim.num)
  
  for (i in 1:sim.num) {
    dat.sim = SimulateDataset(seed = i)
    fit.Bayesian = fit_stan(l.formula.sim, s.formula.sim, dat.sim$long.train, dat.sim$surv.train,
                            n.iter = n.iter, 
                            n.burnin = n.burnin,
                            seed = seed,
                            progress = FALSE)
    result.sim = as.data.frame(fit.Bayesian)
    est = get.estimate(result.sim)[1:length(true_value),]
    estimate[,i] = est$Estimate
    SE[,i] = est$StdErr
    cover[, i] = ifelse((true_value > est$LowerLimit) & (true_value < est$UpperLimit), 1, 0)
  }
  
  Parameter = est$Parameter
  mean_estimate = rowMeans(estimate)
  StdDev = sqrt(rowMeans((estimate - true_value)^2))
  sr_MSE2 = sqrt(rowMeans(SE^2))
  coverage = rowMeans(cover)
  result = data.frame(Parameter = Parameter,
                      TrueVal = true_value,
                      Estimate = mean_estimate,
                      StdDev = StdDev,
                      StdDev_hat = sr_MSE2,
                      Coverage = coverage)
  
  return(result)
}