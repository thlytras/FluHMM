#' Generate posterior samples from a `FluHMM' object
#'
#' This function generates new posterior samples from a `FluHMM' object, or runs more iterations
#' to increase the length of the existing ones.
#'
#' @param x An object of class `FluHMM', from which to generate posterior samples.
#' @param iter Number of iterations to run
#' @param thin Thinning interval between consequtive iterations to keep.
#' @param Rhat Gelman-Rubin diagnostic cutoff value to determine convergence.
#' @param enlarge If \code{FALSE} (the default), a new posterior sample is generated and the old
#'    one discarded. If \code{TRUE}, the new posterior sample is appended to the old one to
#'    increase its size. If you do that you must make sure you use the same thinning interval.
#'
#' @details After reaching initial convergence (i.e. convergence for the sigma[1] parameter, which
#'    is the standard deviation of the pre- and post-epidemic phases), \code{update} should be run
#'    to generate a fresh posterior sample, for as many iterations as required based on the desired
#'    level of precision.
#'
#'    Afterwards, convergence should be checked by examining the \code{converged} and \code{gelman}
#'    elements of the FluHMM object. If the chains have not converged, \code{update} should be run
#'    again (or alternatively \code{\link{autoUpdate}}). However, parameters beta[2] and beta[3]
#'    (the slopes of the epidemic growth and epidemic decline phases) are naturally slow-mixing and
#'    thus convergence may take longer to achieve. If all other parameters have converged (i.e. the
#'    Gelman-Rubin statistic is below 1.1) but beta[2] and beta[3] have not, it is recommended to
#'    *extend* the sample by running \code{update} with option \code{enlarge=TRUE}.
#'
#' @return None. The function mutates its argument `x' directly.
#'
#' @export
update.FluHMM <- function(x, iter=5000, thin=1, Rhat=1.1, enlarge=FALSE) {
  if (is.null(x$model)) stop("Object has been \"compressed\" for storage. No further MCMC sampling possible")
  t <- unname(system.time({
    varnamesA <- c("muPre", "beta", "sigma", "P12", "P23", "P24", "P34", "P45", "state", "mu")
    varnamesB <- c("muPre", "beta[1]", "beta[2]", "beta[3]", "beta[4]", "sigma[1]", "sigma[2]",
            "P12", "P23", "P24", "P34", "P45")
    cSample <- coda.samples(x$model, var=varnamesA, n.iter=iter, thin=thin)
    if (enlarge) {
      a <<- x$cSample
      b <<- cSample
      cSample <- c(x$cSample, cSample)
    }
    states <- calculateStates(cSample)
    mu <- summary(cSample[,grep("mu\\[", varnames(cSample))])[[1]][,1]
    params <- summary(cSample[, varnamesB])[[1]][,1:2]
    iterCount <- attributes(cSample[[1]])$mcpar[2]
    gelman <- gelman(cSample)
    converged <- (sum(gelman[,1]>Rhat)==0)
    initConv <- (gelman["sigma[1]",1]<Rhat)
    eval.parent(substitute( x[["cSample"]] <- cSample ))
    eval.parent(substitute( x[["states"]] <- states ))
    eval.parent(substitute( x[["mu"]] <- mu ))
    eval.parent(substitute( x[["params"]] <- params ))
    eval.parent(substitute( x[["iterCount"]] <- attributes(cSample[[1]])$mcpar[2] ))
    eval.parent(substitute( x[["gelman"]] <- gelman ))
    eval.parent(substitute( x[["converged"]] <- converged ))
    eval.parent(substitute( x[["initConv"]] <- initConv ))
  })[3])
  eval.parent(substitute(x[["elapsedTime"]] <- x$elapsedTime + t ))
  invisible()
}

