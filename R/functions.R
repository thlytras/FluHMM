#' Compress a FluHMM object
#'
#' This methods removes the MCMC samples and model object from the FluHMM object (elements
#' \code{cSample} and \code{model} of the FluHMM object), thus substantially reducing its size.
#' After this operation results are still available and the object can be plotted, but it is
#' not possible to generate further MCMC samples or run MCMC diagnostics. Use this function
#' to compress the object for long-term storage, when you are satisfied with the output.
#'
#' @param x Object (of class `FluHMM') to be compressed.
#'
#' @return The function does NOT return anything, it mutates its argument directly.
#'
#' @export
compress <- function(x) {
  eval.parent(substitute( x[["cSample"]] <- NULL ))
  eval.parent(substitute( x[["model"]] <- NULL ))
}


calculateStates <- function(cSample) {
  states <- cSample[,grep("state", varnames(cSample))]
  states <- t(sapply(1:length(varnames(states)), function(i) table(unlist(states[,i]))[as.character(1:5)]))
  states[is.na(states)] <- 0
  states <- 100 * states / sum(states[1,])
  colnames(states) <- NULL
  states
}


gelman <- function(x, i=1:length(x)) {
    vars <- c("muPre", "sigma[1]", "sigma[2]", "beta[1]", "beta[2]", "beta[3]", "beta[4]",
                "P12", "P23", "P24", "P34", "P45",
                "muPreIsol", "sigmaIsol[1]", "sigmaIsol[2]", 
                "betaIsol[1]", "betaIsol[2]", "betaIsol[3]", "betaIsol[4]")
    vars <- vars[vars %in% varnames(x)]
    res <- gelman.diag(x[i][,vars])$psrf
    colnames(res) <- c("Point", "UpperCI")
    res
}

