#' Plot a FluJointHMM object
#'
#' This function plots the ILI/ARI rates and number of influenza-positive isolates stored
#' in a FluJointHMM object, superimposing the model results. The input and output is
#' mostly the same as \code{\link{plot.FluHMM}}, with an additional smaller sub-plot
#' of the isolates at the bottom.
#'
#' @param x An object of class `FluJointHMM' to be plotted.
#' @param xlab Label for the x-axis.
#' @param ylab Label for the y-axis of the ILI/ARI rate sub-plot.
#' @param ylabIsol Label for the y-axis of the isolates sub-plot.
#' @param main Main label of the plot.
#' @param xaxis How to annotate the x-axis of the plot. If \code{NULL}, no axis is being drawn
#'    (equivalent to \code{xaxt='n'}). If \code{xaxis=0}, the axis names are taken from the
#'    \code{names} attribute of the \code{seasonRates} element of the FluHMM object. If
#'    \code{xaxis=NA}, the \code{seasonRates} names if treated as strings of format YYYYWW, and
#'    the week number is plotted. Otherwise, if a vector of length>1 is supplied, it is treated
#'    as the week numbers to be plotted at the x-axis.
#' @param showPs If \code{TRUE} (the default), show the weekly posterior probabilities of each phase
#'    as a series of numbers too, instead of only a color representation.
#' @param cexPs Character expansion factor (cex), i.e. size for plotting the posterior probabilities
#' @param yexpand Factor (as a percentage of the original y-axis size) to expand the y-axis, so that the
#'    numbers of the posterior probabilities do not overlap with the rest of the plot.
#' @param col Color for plotting the weekly ILI/ARI rates
#' @param mucol Color for plotting the fitted weekly mean rates.
#' @param hues A numeric vector of length 5, with values between 0 and 1, containing the hue values
#'    for each of the five phases in the model.
#' @param rainbow If \code{TRUE}, the mean rates and mean number of isolates for each of the
#'    six chains are individually plotted as well.
#' @param ci If \code{TRUE}, semi-transparent 95% confidence bands are plotted around the rates,
#'    provided the `FluJointHMM' object includes a logSE element (the log Standard Errors of the rates).
#' @param alpha Alpha transparency value for the confidence bands (a number between 0 and 1).
#'
#' @details This function plots the ILI/ARI rates with the week number on the x-axis. If the FluJointHMM
#'    object has a \code{seasonRates} element, this is plotted as a thin grey line and the \code{rates}
#'    are overlaid as a thick line with points, of color \code{col}.
#'
#'    A similar subplot (a bar plot) of the influenza-positive isolates is drawn under the ILI/ARI rates plot.
#'
#'    The posterior probabilities of the five epidemic phases per week (pre-epidemic, epidemic growth,
#'    epidemic plateau, epidemic decline and post-epidemic) are displayed with colored bars on the top
#'    of the plot and optionally as a vertical stack of numbers (if \code{showPs=TRUE}). The hue of
#'    each phase in the color representation is given in the \code{hues} argument, and the saturation
#'    is dependent on the posterior probability of each phase. If the vertical space of the plot is
#'    insufficient (especially if \code{showPs=TRUE}), expand it as necessary by increasing the
#'    \code{yexpand} argument.
#'
#'    The mean weekly rates fitted by the model is plotted as a thick dotted line. If \code{rainbow=TRUE},
#'    the mean rates for each MCMC chain are plotted individually as well. These should be very close to
#'    each other; if very different, then the chains have probably not converged enough.
#'
#'    In any case, if full convergence (for all model parameters) has not been reached, a clear warning
#'    will be displayed on the plot.
#'
#' @return None
#'
#' @export
plot.FluJointHMM <- function(x, xlab="Week", ylab="ILI rate", ylabIsol="Influenza(+) samples",
            main=NA, xaxis=NA, showPs=TRUE, cexPs=0.7, yexpand=0.3,
            col="red", mucol="limegreen", hues=c(4,0,2,5,3)/6, rainbow=FALSE, ci=TRUE, alpha=0.1) {
  layout(matrix(1:2), heights=c(3,2))
  par(mar=c(2,4,4,2))
  plot.FluHMM(x, xlab="", ylab=ylab, main=main, xaxis=xaxis, showPs=showPs, cexPs=cexPs, yexpand=yexpand,
            col=col, mucol=mucol, hues=hues, rainbow=rainbow, ci=ci, alpha=alpha)
  abline(v=0.5+(0:length(x$seasonRates)), col="lightgrey", lty="dotted")
  pos <- c(barplot(x$isolates[1:length(x$seasonRates)], plot=FALSE))
  par(mar=c(5,4,2,2))
  barplot(x$isolates[1:length(x$seasonRates)], col=col, border="white", xlim=range(pos),
    ylim=c(0, ceiling(max(qchisq(0.975, 2*(x$isolates+1))/2)*1.1)), xlab=xlab, ylab=ylabIsol)
  abline(v=c(pos[1]-unique(diff(pos))[1]/2, pos+unique(diff(pos))[1]/2),
    col="lightgrey", lty="dotted")
  if (ci) {
    drawBand(x=pos[1:length(x$isolates)], col=addalpha(col, alpha),
        y.lo=qchisq(0.025, 2*x$isolates)/2,
        y.hi=qchisq(0.975, 2*(x$isolates+1))/2)
  }
  barplot(x$isolates[1:length(x$seasonRates)], col=col, border=col, add=TRUE)
  points(x=pos[1:length(x$muIsol)], y=x$muIsol, type="l", col=mucol, lwd=3, lty="dotted")
  if (!is.null(xaxis)) {
    if (length(xaxis)==1) {
      if (is.na(xaxis)) {
        wklab <- as.integer(names(x$seasonRates)) %% 100
      } else if (xaxis==0) {
        wklab <- names(x$seasonRates)
      }
    } else {
      wklab <- xaxis[1:length(x$seasonRates)]
    }
    axis(1, labels=wklab, at=pos)
  }
  if (rainbow) {
    for (i in 1:6) points(x=pos[1:length(x$muIsol)],
      y=summary(x$cSample[[i]][,grep("muIsol\\[", varnames(x$cSample))])[[1]][,1],
      type="l", col=rainbow(6)[i])
  }
  if (!is.null(x$descr)) {
    # Mark the most likely first epidemic week, if probability > 0.5
    if (sum(x$states[nrow(x$states), 2:5])>50 && x$descr$firstEpiWeek[1,3]>50) {
      text(pos[x$descr$firstEpiWeek[1,"i"]],
        y=ceiling(max(qchisq(0.975, 2*(x$isolates+1))/2)*1.1),
        "*", srt=90, adj=1, cex=2)
    }
    # Mark the most likely peak intensity week, if probability > 0.5
    if (sum(x$states[nrow(x$states), 3:5])>50 && x$descr$peakWeek[1,3]>50) {
      text(pos[x$descr$peakWeek[1,"i"]],
        y=ceiling(max(qchisq(0.975, 2*(x$isolates+1))/2)*1.1),
        "***", srt=90, adj=1, cex=2)
    }
  }
}



addalpha <- function(colors, alpha=1.0) {
  r <- col2rgb(colors, alpha=T)
  # Apply alpha
  r[4,] <- alpha*255
  r <- r/255.0
  return(rgb(r[1,], r[2,], r[3,], r[4,]))
}



drawBand <- function(x, y.lo, y.hi, col) {
  if (length(x)!=length(y.lo) || length(x)!=length(y.hi)) {
    stop("Arguments x, y.lo, y.hi must have the same length")
  }
  for (i in 1:length(x)) {
    if (sum(is.na(c(y.lo[i:(i+1)], y.hi[i:(i+1)])))==0) {
      polygon(x=x[i+c(0,1,1,0)], y=c(y.lo[i:(i+1)], y.hi[(i+1):i]), col=col, border=NA)
    } else if (sum(is.na(c(y.lo[i], y.hi[i])))==0 && (i==1 || sum(is.na(c(y.lo[i-1], y.hi[i-1])))>0)) {
      polygon(x=x[c(i,i)], y=c(y.lo[i], y.hi[i]), col=col, border=col)
    }
  }
}

