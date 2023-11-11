\name{calipeRmatching-package}
\alias{calipeRmatching-package}
\alias{calipeRmatching}
\docType{package}
\title{
  Estimates average treatment effects using the caliper matching estimator
}
\description{
  Caliper matching estimator of the Average Treatment Effect (ATE) and the Average Treatment Effect on the Treated (ATT)
  of a binary treatment. These effects are estimated by matching on the \itemize{
    \item known or 
    \item estimated parametric
  } propensity score.
  See \href{https://arxiv.org/abs/2304.08373}{[KPV23]} for more details. This package is an R wrapper to the \href{https://github.com/kmmate/caliper_matching}{C implementation of the caliper matching estimator}.
}
\details{

  The two main functions of the package are \itemize{
    \item \code{\link{cm_cm_known_propscore}}
    \item \code{\link{cm_cm_estimated_propscore}}
  }
  that estimate the ATE and ATT via caliper matching, and, if specified, their asymptotic variances.
  The ATE and ATT estimates are constructed by comparing the outcome of matched treated and control units: a treated and 
  a control unit are matched if their propensity scores are closer than a specified distance called the caliper. 

  \code{\link{cm_cm_known_propscore}} and \code{\link{cm_cm_estimated_propscore}} use the known and the estimated parametric propensity score, respectively, for matching. 
  Correspondingly, \code{\link{cm_cm_known_propscore}} requires the known propensity score as an input, 
  while \code{\link{cm_cm_estimated_propscore}} needs the covariate data, the type of the parametric propensity score model
  and the estimated propensity score parameter. In addition, the functions produce statistics on the matches.
  The caliper can be set manually or automatically.

  The main functions heavily exploit the properties of caliper matching to compute the estimates very efficiently 
  by processing observations simultaneously in a parallel fashion on multiple threads. The package has a setup function \code{\link{cm_set_number_of_threads}} that can be used to set the number of threads for 
  parallel execution. Unless a different number is specified, the number of parallel threads used is equal to the 
  number of logical cores on the machine.

  The function \code{\link{test_cm}} can be used to run the extensive test suit of the underlying C library.

  Functions prepended by `cpp_` are for internal use only.
}
\author{
Mate Kormos

Maintainer: Mate Kormos <mate.kormos@runbox.com>
}
\references{
  [KPV23] Kormos, V. d. Pas, V. d. Vaart (2023): Asymptotics of Caliper Matching Estimators for Average Treatment Effects, https://arxiv.org/abs/2304.08373
}
\keyword{ caliper matching }
%\examples{
%  \dontrun{
%     ## Optional simple examples of the most important functions
%     ## These can be in \dontrun{} and \donttest{} blocks.   
%  }
%}
