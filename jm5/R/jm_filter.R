check_exist <- function(variable){
  for (i in 1:length(variable)){
    if (!exists(variable[i])){
      err_message <- paste0('Variable ', l.variable, ' does not exist!')
      stop(err_message)
    }
  }
}

get_variable <- function(variable){
  ## prespecify the data matrix as n.row = length of first variable
  ## n.col = 100 (we discard those empty columns later)
  dat.mat <- as.data.frame(matrix(NA, length(get(variable[1])), 100))
  total.cols <- 0
  for (i in 1:length(variable)){
    tmp <- get(variable[i])
    vars <- model.matrix(~ tmp)
    ## if the factor has more than 2 levels, then we need to loop
    for (j in 2:ncol(vars)){
      vname <- str_remove(colnames(vars)[j], 'tmp')
      vname <- paste0(variable[i], vname)
      total.cols <- total.cols + 1
      colnames(dat.mat)[total.cols] <- vname
      dat.mat[, total.cols] <- vars[, j]
    }
  }
  return(dat.mat[, 1:total.cols])
}

#' Create the longitudinal and survival datasets
#'
#' This function allows you to pass in longitudinal and survival
#' formula, with the full dataset and ID variable, to create two
#' separate datasets: longitudinal dataset and survival dataset.
#' 
#' @param l.formula the formula for longitudinal model
#' @param s.formula the formula for survival model
#' @param dat original dataset
#' @param ID ID variable
#' 
#' @return a length-2 list: longitudinal and survival dataframe
#'  
#' @examples
#' 
#' jm_filter(CD4 ~ obstime, Surv(Time, death) ~ drug, JM::aids, JM::aids$patient)
#' 
#' @importFrom stringr str_remove
#' 
#' @export
jm_filter <- function(l.formula, s.formula, dat, ID){
  ## attach the dataset, so that we can check existence of variables
  attach(dat) 
  
  ## check existence of variable in longitudinal and survival formula
  l.variable <- all.vars(l.formula)
  s.variable <- all.vars(s.formula)
  
  check_exist(l.variable)
  check_exist(s.variable)
  
  ## if variables all exist in dataset, then we filter the variables to a matrix
  
  l.mat <- get_variable(l.variable)
  l.mat$ID <- ID
  s.mat <- get_variable(s.variable)
  s.mat$ID <- ID
  s.mat <- unique(s.mat)
  
  if (sum(!complete.cases(l.mat))!=0){
    stop('There are missing values in longitudinal data!')
  }
  if (sum(!complete.cases(s.mat))!=0){
    stop('There are missing values in survival data!')
  }
  
  detach(dat)
  l <- list(long.dat = l.mat, surv.dat = s.mat)
  return(l)
}