#' Create cross-validation data
#'
#' This function allows you to pass in longitudinal and survival
#' datasets with the cv ratio and replicates, 
#' to create cross-validation partitions (training data and testing data).
#' 
#' @param long.dat the dataset for longitudinal model
#' @param surv.dat the dataset for survival model
#' @param ratio ratios of subjects in training dataset and testing dataset (default=3)
#' @param reps how many replicates needed (default=10)
#' @param seed random seed (default = 123)
#' 
#' @return a length-reps list, each element contains 
#' longitudianl and survival training and testing dataset
#'  
#' @examples
#' 
#' create_cv(
#' jm_filter(CD4 ~ obstime, Surv(Time, death) ~ drug, JM::aids, JM::aids$patient)[[1]],
#' jm_filter(CD4 ~ obstime, Surv(Time, death) ~ drug, JM::aids, JM::aids$patient)[[2]])
#' 
#' 
#' @export
create_cv <- function(long.dat, surv.dat, ratio = 3, reps = 10, seed = 123){
  set.seed(seed)
  
  n <- nrow(surv.dat) ## subjects
  n.train <- floor(n*ratio/(ratio+1)) ## training data subjects
  n.test <- n - n.train ## testing dataset subjects
  cv.dat <- vector('list', length = 10)
  for (i in 1:reps){
    ID.train <- sort(sample(1:n, n.train))
    ID.test <- (1:n)[-ID.train]
    
    long.train <- long.dat[long.dat$ID %in% ID.train, ]
    long.train.ID.count <- as.numeric(table(long.train$ID))
    long.test <- long.dat[long.dat$ID %in% ID.test, ]
    long.test.ID.count <- as.numeric(table(long.test$ID))
    surv.train <- surv.dat[surv.dat$ID %in% ID.train, ]
    surv.test <- surv.dat[surv.dat$ID %in% ID.test, ]
    
    ## rename the ID ##
    long.train$ID <- rep(1:n.train, long.train.ID.count)
    long.test$ID <- rep((n.train+1):n, long.test.ID.count)
    surv.train$ID <- 1:n.train
    surv.test$ID <- (n.train+1):n
    
    l <- list(long.train = long.train, long.test = long.test,
              surv.train = surv.train, surv.test = surv.test)
    cv.dat[[i]] <- l
  }
  return(cv.dat)
}