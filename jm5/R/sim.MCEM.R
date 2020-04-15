#' Simulation for MCEM method
#' 
#' This function allows you to evaluate the correctness of MCEM method
#' using simulated data
#' 
#' @param sim.num number of simulations (default = 10)
#' @param max.iter maximum iterations allowed (default = 30)
#' @param tol tolerance (default = 1e-3)
#' @param seed random seed (default = 123)
#' @param progress whether to show progress for each iteration (default = TRUE)
#' 
#' @return a data frame containing parameter name, true values of parameters, 
#' mean of estimated parameters, SD of estimated parameters.
#' 
#' @example 
#' 
#' sim.MCEM(sim.num = 2)
#' 
#' @export
sim.MCEM = function(sim.num = 10, max.iter = 30, tol = 1e-3, seed = 123, progress = TRUE){
  
  l.formula.sim <- Y ~ obstime + x
  s.formula.sim <- Surv(Time, death) ~ w
  true_value = c(20, -1, -4, 2, 1, -2, -1, -0.2)
  estimate = matrix(NA, length(true_value), sim.num)
  
  for (i in 1:sim.num) {
    dat.sim = SimulateDataset(seed = i)
    fit.MCEM = MCEM(l.formula.sim, s.formula.sim, dat.sim$long.train, dat.sim$surv.train,
                    max.iter = max.iter, tol = tol, seed = seed, progress = progress)
    estimate[,i] = fit.MCEM
  }
  
  Parameter = names(fit.MCEM)
  mean_estimate = rowMeans(estimate)
  StdDev = sqrt(rowMeans((estimate - true_value)^2))
  
  result = data.frame(Parameter = Parameter,
                      TrueVal = true_value,
                      Estimate = mean_estimate,
                      StdDev = StdDev)
  
  return(result)
}