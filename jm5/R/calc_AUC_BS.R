check_exist <- function(variable){
  for (i in 1:length(variable)){
    if (!exists(variable[i])){
      err_message <- paste0('Variable ', l.variable, ' does not exist!')
      stop(err_message)
    }
  }
}

calc_AUC_BS <- function(l.formula, s.formula, long.test, surv.test, 
                        Tstart, Tstop, summary.mean,
                        n.samples, n.burnin){
  
  ## Check existence of variable
  attach(long.test)
  l.variable <- c(all.vars(l.formula), 'ID')
  check_exist(l.variable)
  detach(long.test)
  
  attach(surv.test)
  s.variable <- c(all.vars(s.formula), 'ID')
  check_exist(s.variable)
  detach(surv.test)
  
  ## filter the longitudinal data and survival data who survives up to Tstart
  
  
  
}