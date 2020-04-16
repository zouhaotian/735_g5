#' Simulation for Bayesian method
#' 
#' This function allows you to evaluate the correctness of Bayesian method
#' using simulated data
#' 
#' @param sim.num number of simulations (default = 10)
#' @param n.iter total number of iterations in NUTS (default = 2000)
#' @param n.burnin number of burnins in NUTS (default = n.iter/2)
#' @param seed random seed (default = 123)
#' @param progress whether to show progress for sampling (default = TRUE)
#' 
#' @return a data frame containing parameter name, true values of parameters, 
#' mean of estimated parameters, SD of estimated parameters, square root of 
#' mean of SE^2 and coverage probability.
#' 
#' @examples 
#' 
#' sim.Bayesian(sim.num = 2)
#' 
#' @export
sim.Bayesian = function(sim.num = 10,
                        n.iter = 2000, 
                        n.burnin = floor(n.iter/2),
                        seed = 123, progress = TRUE){
  
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
                            progress = progress)
    result.sim = as.data.frame(fit.Bayesian[1:length(true_value), ])
    estimate[,i] = result.sim$mean
    SE[,i] = result.sim$sd
    cover[, i] = ifelse((true_value > result.sim$`2.5%`) & (true_value < result.sim$`97.5%`), 1, 0)
  }
  
  Parameter = rownames(result.sim)
  mean_estimate = rowMeans(estimate)
  
  StdDev <- rep(NA, length(true_value))
  for (i in 1:length(true_value)){
    StdDev[i] <- sd(estimate[i, ])
  }

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