#' Combine an mcmc object
#'
#' This method combines `mcmc' objects into a single `mcmc' object.
#' Objects must have come from the same model and have the same variable names.
#'
#' @param ... `mcmc' objects to be concatenated
#' @param start Iteration number of first observation in the chain, for the resulting object.
#' If \code{NA}, the first iteration number of the first `mcmc' object is used.
#' @param thin Thinning interval between consecutive observations, for the resulting object.
#' If \code{NA}, the thinning interval of the first `mcmc' object is used.
#'
#' @return Returns the combined `mcmc' object.
#'
#' @export
c.mcmc <- function(..., start=NA, thin=NA) {
  l <- list(...)
  if (is.na(start)) {
    start <- attr(l[[1]], "mcpar")[1]
  }
  if (is.na(thin)) {
    thin <- attr(l[[1]], "mcpar")[3]
  }
  v <- lapply(l, varnames)
  if (!all(sapply(2:length(v), function(i) identical(v[1], v[i])))) {
    stop("Not all 'mcmc' objects have the same varnames!")
  }
  mcmc(do.call(rbind, lapply(l, as.matrix)), start=start, thin=thin)
}


#' Combine an mcmc.list object
#'
#' This method combines `mcmc.list' objects into a single `mcmc.list' object.
#' Objects must have come from the same model, and have the same
#' variable names and same number of chains.
#'
#' @param ... `mcmc.list' objects to be concatenated
#'
#' @return Returns the combined `mcmc.list' object.
#'
#' @export
c.mcmc.list <- function(...) {
  l <- list(...)
  nchains <- unique(sapply(l, length))
  if (length(nchains)>1) {
    stop("Not all 'mcmc.list' objects have the same number of chains!")
  }
  do.call(mcmc.list, lapply(1:nchains, function(i) {
    do.call(c, lapply(l, function(x) x[[i]]))
  }))
}

