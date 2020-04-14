#Data Simulation based on Proposal_v3
# We compare the model fitting results using the simulation set-up as follows:
  
# We simulate n=500 subjects, and time points as (0, 2, 6, 12, 18)
# as from the `aids' data in the longitudinal model, 
# and another 200 subjects as the testing dataset.

SimulateDataset = function(seed, num){
  
  #Set Seed
  set.seed(seed)
  
  #number of individuals
  N = 700 #TOTAL; split into N_train = 500 and N_test = 200 later;
  
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
  drug = as.vector(rbinom(N, 1, 0.49)) #ddI = 1 #this is wi
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
  #hi_ti = exp(-2)*exp(-1*dat3$drug + 0.2*mi_ti) #hazard function for subject i 
  
  S_t = runif(N) #survival probability (random uniform distribution)
  alpha = -0.2
  beta0 = 20
  beta1 = -1
  beta2 = -4
  gamma = -1
  h0 = exp(-1)
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
  #head(surv.dat.sim, 10) #This is our simulated survival dataset
  assign(paste("surv.dat.sim", num, sep=""),surv.dat.sim)
  
  #longitudinal
  Ti = as.vector(Ti)
  dat4 = cbind(dat2, rep(Ti,each = p))
  colnames(dat4)[8] = "Ti"
  long.dat.sim<-dat4[!(dat4$Ti < dat4$time),]
  #head(long.dat.sim,10)  #This is our corresponding simulated longitudinal dataset
  #View(long.dat.sim)
  assign(paste("long.dat.sim", num, sep=""),long.dat.sim)
  
  #datasets_final = list(noquote(paste0("long.dat.sim", num)) = long.dat.sim, noquote(paste0("surv.dat.sim", num)) = surv.dat.sim)
  
  datasets_final = list(Longitudinal = long.dat.sim, Survival= surv.dat.sim)
  return(datasets_final)
}

SimulateDataset(20201, 1)
SimulateDataset(20202, 2)
SimulateDataset(20203, 3)
SimulateDataset(20204, 4)
#etc

#500 training
#200 testing
#Look at create_cv.R for how to split (JM5 package)
