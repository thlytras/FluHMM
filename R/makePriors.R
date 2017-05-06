#' Create informative prior from other FluHMM objects
#'
#' This function uses past end-of-season FluHMM objects and constructs
#' informative priors for the \code{beta[1]}, \code{binc}, \code{bprop} 
#' and \code{sigma[1]} model parameters.
#'
#' @param ... `FluHMM' objects to be used
#' @param par Which parameters to create priors for? A character vector.
#'
#' @return Returns a named list with the prior distribution parameters for 
#' each model parameter.
#'
#' @export
makePriors <- function(..., par=c("beta[1]", "beta[4]", "binc", "bprop")) {
  objs <- list(...)
  if (sum(sapply(objs, class)!="FluHMM")>0) {
    stop("Not all provided objects are of class 'FluHMM'!")
  }
  objs <- c(unlist(objs[which(sapply(objs, class)=="list")], recursive=FALSE), 
        objs[which(sapply(objs, class)=="FluHMM")])
  par <- par[par %in% c("beta[1]", "beta[4]", "binc", "bprop")]
  priors <- lapply(par, function(v) {
    val <- as.data.frame.matrix(t(sapply(objs, function(x) x$params[v,])))
    list(mu = with(val, mean(Mean)),
        sigma = with(val, sqrt((sum(SD^2) + sum(Mean^2))/nrow(val) - mean(Mean)^2)))
  })
  names(priors) <- par
  for (x in c("bprop", "beta[4]")) {
    priors[[x]]$alpha <- with(priors[[x]], beta_a(mu, sigma))
    priors[[x]]$beta <- with(priors[[x]], beta_b(mu, sigma))
    priors[[x]]$mu <- NULL
    priors[[x]]$sigma <- NULL
  }
  return(priors)  
}


beta_a <- function(m,s) m*(m*(1-m)/s^2 - 1)
beta_b <- function(m,s) beta_a(m,s)*(1-m)/m

