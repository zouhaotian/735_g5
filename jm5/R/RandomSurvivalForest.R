#' Fitting survival model using Random Survival Forest
#'
#' This function allows you to pass in survival dataset, using
#' Random Survival Forest to fit the survival model
#' 
#' @param surv.dat2 the survival dataset
#' @param s.formula the formula for survival model
#' @param ntree number of trees (default = 100)
#' @param tol tolerance (default = 1e-3)
#' @param seed random seed (default = 12)
#' 
#' @return a RSF fit object
#' 
#' @importFrom randomForestSRC rfsrc
#' 
#' @export
RandomSurvivalForest = function(surv.dat, s.formula, ntree=100,seed=12,...)  {
  # set seed
  set.seed(seed)
  
  fit = rfsrc(s.formula, data = surv.dat, ntree = ntree,statistics=T, ...)
  return(fit)
}
