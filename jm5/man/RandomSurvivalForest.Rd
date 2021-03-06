% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RandomSurvivalForest.R
\name{RandomSurvivalForest}
\alias{RandomSurvivalForest}
\title{Fitting survival model using Random Survival Forest}
\usage{
RandomSurvivalForest(surv.dat, s.formula, ntree = 100, seed = 12, ...)
}
\arguments{
\item{s.formula}{the formula for survival model}

\item{ntree}{number of trees (default = 100)}

\item{seed}{random seed (default = 12)}

\item{surv.dat2}{the survival dataset}

\item{tol}{tolerance (default = 1e-3)}
}
\value{
a RSF fit object
}
\description{
This function allows you to pass in survival dataset, using
Random Survival Forest to fit the survival model
}
