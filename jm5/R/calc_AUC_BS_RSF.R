filter_surv <- function(surv.test, Tstart){
  Time <- surv.test[, 1]
  surv.ID <- surv.test[Time>=Tstart, ]$ID
  surv.filtered <- surv.test[surv.test$ID %in% surv.ID, ]
  return(surv.filtered)
}

extract_prob <- function(surv.prob.test.ID, toi, Tstart){
  return(surv.prob.test.ID[min(which(Tstart<=toi))])
}

#' Calculate time-dependent AUC and BS for RSF
#'
#' This function allows you to pass in survival test dataset, 
#' start time, stop time and predicted RSF result on test dataset, 
#' to calculate time-dependent Area Under Curve for ROC, and Brier Score
#' 
#' @param surv.dat2 the survival test set
#' @param Tstart starting time
#' @param Tstop stop time
#' @param res predicted result using RSF on test dataset 
#' 
#' @return a length-2 vector: estimated AUC and BS
#'  
#' 
#' @importFrom survival Surv
#' @importFrom tdROC tdROC
#' @importFrom ipred sbrier
#' 
#' @export
calc_AUC_BS_RSF <- function(surv.test2, Tstart, Tstop, res){
  ## extract the predicted survival probability at time of interest
  toi <- res$time.interest
  surv.prob <- res$survival
  
  surv.filtered <- filter_surv(surv.test2, Tstart)
  ID.unique <- surv.filtered$ID
  surv.prob.est <- rep(NA, length(ID.unique))
  ID.index <- 0
  for (test.ID in ID.unique){
    ID.index <- ID.index + 1
    surv.prob.test.ID <- surv.prob[surv.test2$ID==test.ID, ]
    surv.prob.est.Tstart <- extract_prob(surv.prob.test.ID, toi, Tstart)
    surv.prob.est.Tstop <- extract_prob(surv.prob.test.ID, toi, Tstop)
    surv.prob.est[ID.index] <- surv.prob.est.Tstop/surv.prob.est.Tstart
  }
  
  ROC.est <- tdROC(X = 1 - surv.prob.est, Y = surv.filtered[, 1], 
                   delta = surv.filtered[, 2], tau = Tstop,
                   span = 0.1, alpha = 0.05,
                   n.grid = 1000, cut.off = 0.5)
  AUC <- ROC.est$AUC$value
  surv.obj <- Surv(surv.filtered[, 1], surv.filtered[, 2])
  BS <- sbrier(surv.obj, surv.prob.est, btime = Tstop)
  return(c(AUC, BS))
}