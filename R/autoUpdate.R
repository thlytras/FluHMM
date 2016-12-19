#' Run MCMC iterations until full convergence
#'
#' This function runs further MCMC iteration to extend the current posterior sample
#' until convergence is reached for all parameters.
#'
#' @param x An object of class `FluHMM', from which to generate posterior samples.
#' @param iter Number of iterations to run.
#' @param maxit Maximum number of iterations to run before giving up.
#' @param Rhat Gelman-Rubin diagnostic cutoff value to determine convergence
#'    (only for the sigma[1] parameter).
#' @param thin Thinning interval. If the beta[2] and beta[3] parameters are very slow to mix,
#'    and your RAM memory is insufficient, it is advisable to use a larger thinning interval
#'    (e.g. 10 or 50), i.e. keep in the posterior sample every \code{thin} MCMC iteration.
#'
#' @details Most parameters in the model converge fairly quickly, but mixing is slower for parameters
#'    beta[2] and beta[3] (i.e. the slopes of the epidemic growth and epidemic decline phases). Initial
#'    convergence of the model is judged by monitoring the sigma[1] parameter and is performed by the
#'    function \code{\link{autoInitConv}}, which is usually run in the \code{\link{FluHMM}} constructor.
#'    If all other parameters have converged but beta[2] and beta[3] have not, due to slow mixing, then
#'    it is advisable to run the MCMC chain for longer, until beta[2] and beta[3] cover their posterior
#'    distributions adequately. This can be done manually using \code{\link{update}} with argument
#'    \code{enlarge=TRUE}. Alternatively \code{autoUpdate} can be used, which extends the sample by
#'    \code{iter} iterations, checks the Gelman-Rubin diagnostic for all parameters, and if convergence
#'    has not been reached it extends again the sample until convergence or until \code{maxit} iterations.
#'
#'    In some cases (usually near the end of the influenza surveillance period, when all five model phases
#'    have occured) mixing is very slow and the chains need to be run for very long. This requires lots of
#'    RAM memory and a fair amount of CPU power to handle the very large posterior sample. In such cases
#'    it is advisable to thin the chain (i.e. by 10 or 20, i.e. keep every 10th or 20th element in the MCMC
#'    chain), to generate a smaller posterior sample with good coverage of the posterior distribution.
#'
#' @return None. The function mutates its argument `x' directly.
#'
#' @export
autoUpdate <- function(x, iter=20000, maxit=200000, Rhat=1.1, thin=1) {
    if (is.null(x$model)) stop("Object has been \"compressed\" for storage. No further MCMC sampling possible")
    cat(sprintf("Auto-updating: sampling from the MCMC chain until\nconvergence or until %d iterations. Please wait...\n", maxit))
    xx <- x
    iter.remaining <- maxit
    while (!xx$converged && iter.remaining>0) {
        iter.next <- ifelse(iter<iter.remaining, iter, iter.remaining)
        iter.remaining <- iter.remaining - iter.next
        cat(sprintf("Not yet converged: sampling for another %d iterations...\n", iter.next))
        # Make sure that sigma[1] is still converged....
        update(xx, iter=iter.next, Rhat=Rhat, thin=thin, enlarge=TRUE)
        eval.parent(substitute( x <- xx ))
    }
    if (xx$converged) {
        cat(sprintf("Model has now converged, good! (After %d total iterations.\n\n", xx$iterCount))
    } else {
        cat(sprintf("Model has still not converged (after %d total iterations).\n", xx$iterCount))
        cat("Run more iterations or check MCMC diagnostics...\n\n")
    }
}


#' Run MCMC iterations until initial convergence
#'
#' This function is not normally called by the user; it is called from \code{\link{FluHMM}}} provided
#' \code{initConv=TRUE} (the default). It generates posterior samples from the model repeatedly
#' until convergence is reached for the sigma[1] parameter (this is called "initial convergence").
#'
#' @param x An object of class `FluHMM', from which to generate posterior samples.
#' @param iter Number of iterations to run.
#' @param maxit Maximum number of iterations to run before giving up.
#' @param Rhat Gelman-Rubin diagnostic cutoff value to determine convergence
#'    (only for the sigma[1] parameter).
#'
#' @details The sigma[1] parameter (standard deviation of the pre-epidemic phase) is of primary
#'    importance in the model, since the pre-epidemic phase comes first and its correct identification
#'    is the basis on which to estimate the other phases. If convergence for sigma[1] has been reached,
#'    the other parameters in the model are very likely to have converged too, with the exception of
#'    beta[2] and beta[3] (slopes of the epidemic growth and epidemic decline phase); the latter mix
#'    more slowly and may necessitate a longer, and possibly thinned, posterior sample.
#'
#'    Therefore "initial convergence" is defined as convergence for the sigma[1] parameter. Unless
#'    this is achieved, it is inadvisable to use the posterior samples for any inference at all. Only
#'    _after_ this has been achieved can a new posterior sample be generated (using \code{\link{update}}).
#'    Then convergence for all parameters is checked again and, if not achieved, a new sample is generated
#'    from scratch or the current one further is further extended.
#'
#' @return None. The function mutates its argument `x' directly.
#'
#' @export
autoInitConv <- function(x, iter=5000, maxit=95000, Rhat=1.1) {
  if (is.null(x$model)) stop("Object has been \"compressed\" for storage. No further MCMC sampling possible!")
  cat("Checking for initial convergence... ")
  if (x$gelman["sigma[1]",1]<Rhat) {
    cat("initial convergence has been reached.\n")
    cat("You can now sample the chain long enough to get reliable estimates!\n\n")
    return()
  }
  cat("initial convergence NOT reached.\n")
  cat(sprintf("Updating in increments of %d iterations until initial convergence,\n", iter))
  cat(sprintf("   or until we do %d iterations. Please wait...", maxit))
  xx <- x
  iter.remaining <- maxit
  while (!(xx$gelman["sigma[1]",]<Rhat) && iter.remaining>0) {
    iter.next <- ifelse(iter<iter.remaining, iter, iter.remaining)
    if (iter.remaining != maxit) {
      cat(sprintf("Not yet reached initial convergence: sampling for another %d iterations...\n", iter.next))
    }
    iter.remaining <- iter.remaining - iter.next
    update(xx, iter=iter.next, Rhat=Rhat, thin=1)
    eval.parent(substitute( x <- xx ))
  }
  if (xx$gelman["sigma[1]",1]<Rhat) {
    eval.parent(substitute( x[["initConv"]] <- TRUE ))
    cat(sprintf("Reached initial convergence, good! (After %d total iterations.\n\n", xx$iterCount))
    cat("You can now sample the chain long enough to get reliable estimates.\n\n")
  } else {
    eval.parent(substitute( x[["initConv"]] <- FALSE ))
    cat(sprintf("Initial convergence still not reached (after %d total iterations).\n", xx$iterCount))
    cat("You should try more iterations before you do any parameter sampling...\n\n")
  }
}


