% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/update.FluHMM.R
\name{update.FluHMM}
\alias{update.FluHMM}
\title{Generate posterior samples from a `FluHMM' object}
\usage{
\method{update}{FluHMM}(x, iter = 5000, thin = 1, Rhat = 1.1,
  enlarge = FALSE)
}
\arguments{
\item{x}{An object of class `FluHMM', from which to generate posterior samples.}

\item{iter}{Number of iterations to run}

\item{thin}{Thinning interval between consequtive iterations to keep.}

\item{Rhat}{Gelman-Rubin diagnostic cutoff value to determine convergence.}

\item{enlarge}{If \code{FALSE} (the default), a new posterior sample is generated and the old
one discarded. If \code{TRUE}, the new posterior sample is appended to the old one to
increase its size. If you do that you must make sure you use the same thinning interval.}
}
\value{
None. The function mutates its argument `x' directly.
}
\description{
This function generates new posterior samples from a `FluHMM' object, or runs more iterations
to increase the length of the existing ones.
}
\details{
After reaching initial convergence (i.e. convergence for the sigma[1] parameter, which
   is the standard deviation of the pre- and post-epidemic phases), \code{update} should be run
   to generate a fresh posterior sample, for as many iterations as required based on the desired
   level of precision.

   Afterwards, convergence should be checked by examining the \code{converged} and \code{gelman}
   elements of the FluHMM object. If the chains have not converged, \code{update} should be run
   again (or alternatively \code{\link{autoUpdate}}). However, parameters beta[2] and beta[3]
   (the slopes of the epidemic growth and epidemic decline phases) are naturally slow-mixing and
   thus convergence may take longer to achieve. If all other parameters have converged (i.e. the
   Gelman-Rubin statistic is below 1.1) but beta[2] and beta[3] have not, it is recommended to
   *extend* the sample by running \code{update} with option \code{enlarge=TRUE}.
}
