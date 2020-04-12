
create_cv <- function(long.dat, surv.dat, ratio = 3, reps = 10, seed = 123){
  set.seed(seed)
  
  n <- nrow(surv.dat) ## subjects
  n.train <- floor(n*ratio/(ratio+1)) ## training data subjects
  n.test <- n - n.train ## testing dataset subjects
  
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
  }
  
}