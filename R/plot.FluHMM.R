#' Plot a FluHMM object
#'
#' This function plots the ILI/ARI rates stored in a FluHMM object,
#' superimposing the model results.
#'
#' @param x An object of class `FluHMM' to be plotted.
#' @param xlab Label for the x-axis.
#' @param ylab Label for the y-axis.
#' @param main Main label of the plot.
#' @param xaxis How to annotate the x-axis of the plot. If \code{NULL}, no axis is being drawn
#'    (equivalent to \code{xaxt='n'}). If \code{xaxis=0}, the axis names are taken from the
#'    \code{names} attribute of the \code{seasonRates} element of the FluHMM object. If
#'    \code{xaxis=NA}, the \code{seasonRates} names if treated as strings of format YYYYWW, and
#'    the week number is plotted. Otherwise, if a vector of length>1 is supplied, it is treated
#'    as the week numbers to be plotted at the x-axis.
#' @param showPs If \code{TRUE} (the default), show the weekly posterior probabilities of each phase
#'    as a series of numbers too, instead of only a color representation.
#' @param yexpand Factor (as a percentage of the original y-axis size) to expand the y-axis, so that the
#'    numbers of the posterior probabilities do not overlap with the rest of the plot.
#' @param col Color for plotting the weekly ILI/ARI rates
#' @param mucol Color for plotting the fitted weekly mean rates.
#' @param hues A numeric vector of length 5, with values between 0 and 1, containing the hue values
#'    for each of the five phases in the model.
#' @param rainbow If \code{TRUE}, the mean rates for each of the six chains are
#'    individually plotted as well.
#' @param ci If \code{TRUE}, semi-transparent 95% confidence bands are plotted around the rates,
#'    provided the `FluHMM' object includes a logSE element (the log Standard Errors of the rates).
#' @param alpha Alpha transparency value for the confidence bands (a number between 0 and 1).
#'
#' @details This function plots the ILI/ARI rates with the week number on the x-axis. If the FluHMM object
#'    object has a \code{seasonRates} element, this is plotted as a thin grey line and the \code{rates}
#'    are overlaid as a thick line with points, of color \code{col}.
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
#'    will be displayed at the bottom left margin of the plot.
#'
#' @return None
#'
#' @export
plot.FluHMM <- function(x, xlab="Week", ylab="ILI rate", main=NA, xaxis=NA, showPs=TRUE, yexpand=0.3,
            col="red", mucol="limegreen", hues=c(4,0,2,5,3)/6, rainbow=FALSE, ci=TRUE, alpha=0.1) {
  plot(x$seasonRates, type="l", col="grey", bty="l", xaxt="n",
            ylim=c(0,(max(x$seasonRates, na.rm=TRUE))*(1.15+yexpand)), xlab=xlab, ylab=ylab, main=main)
  ymax <- par("usr")[4]
  y0 <- (par("usr")[4]-par("usr")[3])*0.9 + par("usr")[3]
  yrow <- (ymax-y0)/5
  abline(v=0.5+(0:length(x$seasonRates)), col="lightgrey", lty="dotted")
  points(x$rates[1:nrow(x$states)], type="o", col=col, cex=0.8, lwd=2, pch=19)
  points(x$mu, type="l", col=mucol, lwd=3, lty="dotted")
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
    axis(1, labels=wklab, at=1:length(wklab))
  }
  if (showPs) {
    text(x=1:nrow(x$states), y=y0-yrow, labels=apply(round(x$states), 1, paste, collapse="\n"), cex=0.7, pos=1)
  }
  for (st in 1:5) {
    for (i in 1:nrow(x$states)) {
      polygon(y=(5-st)*yrow + c(y0,y0+yrow,y0+yrow,y0),
                x=c(-0.5,-0.5,0.5,0.5)+i, border=NA,
                col=hsv(hues[st], x$states[i,st]/100, 1))
    }
  }
  if (!is.null(x$logSE) && ci) {
    drawBand(x=1:length(x$rates), col=addalpha(col, alpha),
    y.lo=exp(log(x$rates) - 1.96*x$logSE),
    y.hi=exp(log(x$rates) + 1.96*x$logSE))
  }
  if (rainbow) {
    for (i in 1:6) points(
      summary(x$cSample[[i]][,grep("mu\\[", varnames(x$cSample))])[[1]][,1],
      type="l", col=rainbow(6)[i])
  }
  if (!x$converged) mtext("WARNING: model has NOT converged", side=1, adj=0, cex=0.8, line=2.5, col="darkred")
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

