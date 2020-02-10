glm2 <- function(formula, dat, family, link, inits, intercept=T, tol, max.iter=50){
  iter=0
  curr.tol = Inf
  while (iter<max.iter & curr.tol>tol){
    if (family=='Poisson'){
      b <- function()
    }
    if (link=='log'){
      g <- function()
    }
  }
}