#' Simulate the longitudinal and survival datasets
#'
#' This function allows you to pass in subjects in training dataset
#' and testing dataset to create four separate datasets: 
#' longitudinal training, longitudinal testing, 
#' survival training, survival testing.
#' 
#' @param n.train number of subjects in training dataset (default=500)
#' @param n.test number of subjects in testing dataset (default=200)
#' @param seed random seed
#' 
#' @return a length-4 list: longitudinal training, longitudinal testing, 
#' survival training, survival testing.
#'  
#' @examples
#' 
#' l = SimulateDataset(seed = 123)
#' head(l[[1]])
#' 
#' @export
SimulateDataset = function(n.train = 500, n.test = 200, seed){
  
  #Set Seed
  set.seed(seed)
  
  #number of individuals
  N = n.train+n.test #TOTAL; split into N_train = 500 and N_test = 200 later;
  
  #number of observations
  time = c(0, 2, 6, 12, 18)
  p = length(time) #each person starts with 5 observations of data (later, changes when event time is simulated)
  
  # First, we simulate the longitudinal responses using the following formula: 
  # Yij = 20 + (-1) x tij + (-4) x xi + ui + eij , 
  # where ui ~ N(0; 2); eij ~ N(0; 1),
  # xi ~ binom(1; 0.4).
  
  ### set up data frame
  dat = data.frame(patient=rep(1:N, times = p), time=rep(time,each=N))
  #dat3 <- dat[order(dat$subject),]
  ui = as.vector(rnorm(N,0,2))
  eij = as.vector(rnorm(N*p,0,1))
  xi = as.vector(rbinom(N, 1, 0.4))
  drug = as.vector(rbinom(N, 1, 0.5)) #ddI = 1 #this is wi
  dat1 = cbind(dat,rep(xi,p),rep(ui,p),eij,rep(drug,p))
  dat2 = dat1[order(dat1$patient,dat1$time),]
  colnames(dat2)[3] = "xi"
  colnames(dat2)[4] = "ui"
  colnames(dat2)[6] = "drug"
  dat2$Yij = 20-1*dat2$time -4*dat2$xi + dat2$ui + dat2$eij
  #head(dat2) #use again later for longitudinal dataset
  
  # Second, we simulate N random uniform distribution, which corresponds to the survival probability. 
  # We calculate the failure time as: hi(Ti) = exp(-2) exp(-1*wi + 0.2 *mi(Ti)), 
  # where mi(Ti) = 20 + (-1) * Ti + (-4) * xi + ui.
  # And we can calculate the survival time using an explicit formula.
  
  #take t = 0 for now to get one row per person
  dat3 = dat2[which(dat2$time==0),]
  
  # T_i^*: the person's death time. 
  # C_i: the person's censoring time. 
  # We simulate T_i^* as stated today. 
  # Determine C_i independent of T_i^*, say constant. 
  # Event time T_i = min(T_i^*, C_i)
  # Delta_i = I(T_i^*<=C_i)
  #mi_ti = 20 + (-1)*dat3$T_i + (-4)*dat3$xi + dat3$ui
  #hi_ti = exp(-2)*exp(-1*dat3$drug - 0.2*mi_ti) #hazard function for subject i 
  
  S_t = runif(N) #survival probability (random uniform distribution)
  alpha = -0.2
  beta0 = 20
  beta1 = -1
  beta2 = -4
  gamma = -1
  h0 = exp(-2)
  wi = dat3$drug
  Ti_star = log(1+(alpha*beta1*log(S_t)/(-h0*exp(gamma*wi+alpha*(beta0+beta2*xi+ui)))))/(alpha*beta1)
  #We want to death rate as 40.25% --> try different values of c and see which c makes the death rate close to 40.25%.
  c = 27 #choose c
  Ci = runif(N, 0, c)
  Ti = pmin(Ti_star, Ci)
  delta_i = ifelse((Ti_star <= Ci),1,0)
  death_rate = sum(delta_i)/N #is this close to 40.25%?
  #death_rate #0.42
  
  surv.dat.sim = cbind(dat3[,-c(2,5,7)], Ti_star, Ci, Ti, delta_i)
  surv.dat <- data.frame(Time = surv.dat.sim$Ti, 
                         death = surv.dat.sim$delta_i,
                         w = surv.dat.sim$drug,
                         ID = surv.dat.sim$patient)
  
  #head(surv.dat.sim, 10) #This is our simulated survival dataset
  #assign(paste("surv.dat.sim", num, sep=""),surv.dat.sim)
  
  #longitudinal
  Ti = as.vector(Ti)
  dat4 = cbind(dat2, rep(Ti,each = p))
  colnames(dat4)[8] = "Ti"
  long.dat.sim <- dat4[!(dat4$Ti < dat4$time),]
  long.dat <- data.frame(Y = long.dat.sim$Yij, 
                         obstime = long.dat.sim$time,
                         x = long.dat.sim$xi, 
                         ID = long.dat.sim$patient)
  #head(long.dat.sim,10)  #This is our corresponding simulated longitudinal dataset
  #View(long.dat.sim)
  #assign(paste("long.dat.sim", num, sep=""),long.dat.sim)
  
  #datasets_final = list(noquote(paste0("long.dat.sim", num)) = long.dat.sim, noquote(paste0("surv.dat.sim", num)) = surv.dat.sim)
  
  datasets_final = list(
    long.train = long.dat[long.dat$ID %in% (1:n.train), ], 
    long.test = long.dat[long.dat$ID %in% ((n.train+1):N), ],
    surv.train = surv.dat[surv.dat$ID %in% (1:n.train), ],
    surv.test = surv.dat[surv.dat$ID %in% ((n.train+1):N), ])
  return(datasets_final)
}

