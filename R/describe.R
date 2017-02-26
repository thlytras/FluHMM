#' Describe a `FluHMM' object
#'
#' This method prints various interesting information about the fitted FluHMM object: the
#' number of initial weeks fixed to the pre-epidemic phase (K), the probability that the
#' epidemic wave has begun, the three most likely first epidemic weeks, the probability that
#' peak intensity has been reached, and the three most likely peak intensity weeks.
#'
#' @aliases describe
#'
#' @param x An object of class `FluHMM', to be described.
#' @param recalc If \code{TRUE}, recalculate the most likely first epidemic weeks and most
#'    likely peak intensity weeks, and store them in the object.
#'
#' @details The current probability that the epidemic wave has begun is equal to the sum of the
#'    probabilities of phases 2 to 5 for the current (i.e. latest) week. Similarly, the current
#'    probability that peak intensity has been reached is the sum of the probabilities of phases
#'    3 to 5 for the current week.
#'
#'    The most likely first epidemic week is calculated by iterating over the MCMC chains and
#'    finding, for each iteration, the first week where Z=2 (first week of the epidemic growth
#'    phase). Similarly, the most likely peak intensity week is calculated by iterating and
#'    finding the latest week where Z=2.
#'
#' @return None. The information are printed on the screen, and if \code{recalc=TRUE} they
#'    are also stored in the \code{descr} element of the FluHMM object. Thus the method
#'    mutates its argument `x' directly.
#'
#' @export
describe.FluHMM <- function(x, recalc=FALSE) {
  wkNames <- names(x$rates)
  if (is.null(wkNames)) wkNames <- 1:length(x$rates)

  if (recalc && !is.null(x$cSample)) {
    descr <- calcDescr(x$cSample, wkNames)
    eval.parent(substitute( x[["descr"]] <- descr ))
  } else {
    descr <- x$descr
  }

  cat(sprintf("Current week: i=%s (week %s)\n",
    length(x$rates), wkNames[length(x$rates)]))
  cat(sprintf("Weeks fixed on phase 1 (pre-epidemic): K=%s\n", x$K))
  cat(sprintf("\nProbability that the epidemic has started = %s%%\n",
    round(sum(x$states[nrow(x$states), 2:5]) ,1)))
  cat("Most likely first epidemic weeks:\n")
  print(descr$firstEpiWeek)
  cat(sprintf("\nProbability that peak intensity has been reached = %s%%\n",
    round(sum(x$states[nrow(x$states), 3:5]) ,1)))
  cat("Most likely peak weeks:\n")
  print(descr$peakWeek)

  invisible()
}


#' @export
describe <- function(x) {
  UseMethod("describe",x)
}


calcDescr <- function(cSample, wkNames) {
  stateMat <- sapply(grep("state", varnames(cSample)), function(i) unlist(cSample[,i]))

  peakWeek <- apply(stateMat, 1, function(x) rev(which(x==2))[1])
  peakWeek <- summProb(peakWeek, wkNames, "peakWeek")

  firstEpiWeek <- apply(stateMat, 1, function(x) which(x==2)[1])
  firstEpiWeek <- summProb(firstEpiWeek, wkNames, "firstEpiWeek")

  list(firstEpiWeek=firstEpiWeek, peakWeek=peakWeek)
}

summProb <- function(x, wkNames, label) {
  res <- round(sort(table(x), decreasing=TRUE)[1:3]*100/length(x), 1)
  res <- as.data.frame(res, responseName="Probability")
  res[,1] <- as.integer(as.character(res[,1]))
  res <- cbind(i=res[,1], res)
  res[,2] <- wkNames[res$i]
  names(res)[2] <- label
  res
}

