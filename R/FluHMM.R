#' Create a FluHMM object
#'
#' This function takes the set of rates as input, fits the model in JAGS, and constructs
#' a `FluHMM' object. The object can then be used to generate MCMC samples, summarize and
#' plot the results.
#'
#' @param rates The set of weekly influenza-like illness / acute respiratory infection (ILI/ARI)
#' rates obtained from sentinel surveillance, up to the current week, as a numeric vector.
#' @param seasonRates The set of weekly ILI/ARI rates for the whole season, if available (i.e. if
#' season has been compleated). This allows fitting the model for a partial season, but
#' plotting it overlaid on the whole season, see \code{plot.FluHMM}.
#' @param isolates A optional set of weekly numbers of influenza-positive lab isolates. Does not
#' have to be of equal length with the set of rates. If specified, an object of class
#' `FluJointHMM' is produced (inheriting also from class `FluHMM'), which jointly models both
#' series (the rates and the number of isolates) as observations from the same Hidden Markov
#' states chain.
#' @param weights A vector of of length equal to \code{length(rates)} containing observation
#' weights for the rates. If \code{NULL}, all weights are set equal to 1.
#' @param logSE An optional vector of length equal to \code{length(rates)} containing
#' log standard errors for the rates. If \code{NULL}, the rates are treated in the model as
#' "true" rates, i.e. without measurement error.
#' @param K The first K observations (weeks) of the rates are considered a priori to belong in the
#' pre-epidemic phase of the model. Set this to a higher level if you have lots of observations
#' (more than 25) to speed up fitting of the model, as long as you are confident that the weeks
#' really belong to the pre-epidemic phase.
#' @param initConv If \code{TRUE} (the default), MCMC samples are generated from the chains until
#' initial convergence, see details.
#' @param maxit Maximum number of iterations performet for initial convergence, see \code{\link{autoInitConv}}.
#'
#' @details The function constructs an object of class `FluHMM', which contains all the input,
#' model information and results, and can be processed further as required. The minimum input is
#' the set of weekly ILI/ARI rates (argument \code{rates}) up to the current week. The function fits
#' the appropriate model in JAGS (with or without measurement error, depending on the argument
#' \code{logSE}, and with or without a submodel for the isolates, depending on the argument
#' \code{isolates}), and generates posterior samples for 5000 iterations. Six MCMC chains are used.
#'
#' Then, provided the argument \code{initConv} is \code{TRUE} (the default), the sample for
#' sigma[1] (i.e. the standard deviation of the pre- and post-epidemic phases) is checked for
#' convergence using the Gelman and Rubin diagnostic. This is defined as "initial convergence".
#' If initial convergence has not been reached, the posterior sample is discarded, the chains are
#' sampled again for 5000 iterations, and a new check is made. The process is repeated again until
#' initial convergence is reached or after 95000 iterations. See \code{\link{autoInitConv}} for
#' details.
#'
#' After initial convergence is reached, a *new sample* should be generated for inference using
#' \code{\link{update.FluHMM}}, with the number of iterations dependent on the
#' desired precision. If full convergence is not reached, the object can be \code{\link{autoUpdate}}d
#' until full convergence.
#'
#' @return An object of class `FluHMM', which is a list with the following components:
#' \describe{
#'  \item{model}{The fitted model; an object of class `jags'}
#'  \item{cSample}{An `mcmc.list' object containing the posterior samples for the variables
#'  in the model.}
#'  \item{params}{Mean and standard deviation of the parameters of interest.}
#'  \item{states}{A Nx5 matrix, where N==length(rates), containing the probabilities of each phase per week}
#'  \item{mu}{A vector with the fitted mean rates per week}
#'  \item{elapsedTime}{Total processing time spent fitting the model and sampling from the chains.}
#'  \item{gelman}{The Gelman-Rubin diagnostic for the main parameters in the model}
#'  \item{converged}{\code{TRUE} if full convergence has been reached, i.e. if the Gelman-Rubin diagnostic
#'  is less than 1.1 for all parameters in the model.}
#'  \item{initConv}{\code{TRUE} if "initial convergence" has been reached, i.e. if the Gelman-Rubin
#'  diagnostic for sigma[1] (the standard deviation of the pre- and post-epidemic phases) is less than 1.1 . }
#'  \item{rates}{The ILI/ARI rates that were used as input in the model.}
#'  \item{seasonRates}{The ILI/ARI rates for the entire season, if available.}
#'  \item{weights}{The set of observation weights used (usually a vector of ones).}
#'  \item{logSE}{The log standard error of the rates if available; \code{NULL} otherwise.}
#' }
#' In addition, if the object is also of class `FluJointHMM', it also contains the following elements:
#' \describe{
#'  \item{isolates}{The numbers of isolates that were used as input in the model.}
#'  \item{muIsol}{A vector with the fitted mean number of isolates per week}
#' }
#'
#' @export
FluHMM <- function(rates, seasonRates=rates, isolates=NULL, weights=NULL, logSE=NULL,
            K=3, initConv=TRUE, maxit=95000) {
  if (sum(is.na(rates))>0)
    stop("No missing values allowed in 'rates' vector...\n")
  if (!is.null(weights) && (!is.numeric(weights) || length(rates)!=length(weights)))
    stop("Argument 'logSE' must be a numeric vector of same length as 'rates'")
  if (is.null(weights)) weights <- rep(1, length(rates))
  dat <- list(
    nweeks=length(rates),
    rate=rates,
    SMAX = max(rates),
    SR = diff(range(rates)),
    SSD = sqrt((length(rates)-1)*var(rates)/qchisq(0.025, length(rates)-1)),
    K = K,
    weights = weights
  )
  if (is.null(logSE)) {
    fluModel <- list(fluModelChunkA_noErr, fluModelChunkB)
  } else {
    fluModel <- list(fluModelChunkA_err, fluModelChunkB)
    dat$logSE.rate <- logSE
  }
  .Object <- list()
  class(.Object) <- "FluHMM"
  if (!is.null(isolates)) {
    fluModel <- c(fluModel, fluModelChunkIsol)[c(1,3,2)]
    dat$nweeksIsol <- length(isolates)
    dat$nIsol <- isolates
    dat$SMAXisol <- max(isolates)
    dat$SRisol <- diff(range(isolates))
    dat$SSDisol <- sqrt((length(isolates)-1)*var(isolates)/qchisq(0.025, length(isolates)-1))
    class(.Object) <- c("FluJointHMM", "FluHMM")
  }
  fluModel <- do.call(paste, c(fluModel, sep="\n"))
  if (is.null(isolates)) {
    cat("\nFitting a FluHMM object.\n\nRunning initial 5000 iterations, please wait...\n")
  } else {
    cat("\nFitting a FluJointHMM object.\n\nRunning initial 5000 iterations, please wait...\n")
  }
  ini <- function(chain) {
    ini.values <- list(
      sigma0 = c(0.1, 0.1),
      muPre = 0
    )
    if (!is.null(logSE)) ini.values$trueRate <- rates
    if (!is.null(isolates)) {
      ini.values$sigmaIsol0 <- c(0.1, 0.1)
      ini.values$muPreIsol <- 0
    }
    return(ini.values)
  }
  t <- unname(system.time({
    .Object$model <- jags.model(file = textConnection(fluModel), data = dat, inits=ini,
                        n.chains = 6, n.adapt = 0, quiet=TRUE)
  })[3])
  .Object$elapsedTime <- t
  .Object$rates <- rates
  if (!is.null(isolates)) {
    .Object$isolates <- isolates
  }
  .Object$K <- K
  .Object$seasonRates <- seasonRates
  .Object$weights <- weights
  .Object$logSE <- logSE
  .Object$initConv <- FALSE
  update(.Object, iter=5000, enlarge=FALSE)
  cat("Initial sampling complete.\n")
  if (initConv==TRUE) {
    autoInitConv(.Object, iter=5000, maxit=maxit)
  }
  return(.Object)
}








fluModelChunkA_noErr <- 'model {
  for (i in 1:nweeks) {
    rate[i] ~ dnorm(mu[i], tau[state[i]]*weights[i])
  }'


fluModelChunkA_err <- 'data {
  for (i in 1:nweeks) {
    lograte[i] <- log(rate[i])
  }
}
model {
  for (i in 1:nweeks) {
    logTrueRate[i] <- log(trueRate[i])
    lograte[i] ~ dnorm(logTrueRate[i], tauT[i])
    tauT[i] <- pow(logSE.rate[i], -2)
  }

  for (i in 1:nweeks) {
    trueRate[i] ~ dnorm(mu[i], tau[state[i]]*weights[i])T(0,)
  }'


fluModelChunkB <- '

  # ILI rate part

  sigma0[1] ~ dunif(0.1, SSD)
  sigma0[2] ~ dunif(0.1, SSD)
  sigma[1:2] <- sort(sigma0)

  tau[1] <- pow(sigma[1], -2)
  for (i in 2:4) {
    tau[i] <- pow(sigma[2], -2)
  }
  tau[5] <- pow(sigma[1], -2)
  mu[1] <- muPre
  for (i in 2:nweeks) {
    mu[i] <- max(0, equals(state[i],1)*(mu[i-1] + beta[1]) +
        equals(state[i],2)*(mu[i-1] + beta[2]) +
        equals(state[i],3)*(mu[i-1]) +
        equals(state[i],4)*(mu[i-1] + beta[3]) +
        equals(state[i],5)*(mu[i-1]*beta[4]))
  }

  muPre ~ dnorm(0, 0.001)T(0,SMAX)

  binc ~ dnorm(0, 0.001)T(0,SR*2)
  bprop ~ dbeta(2,2)

  beta[1] ~ dnorm(0, 0.001)T(0,SR)
  beta[2] <- max(beta[1]*2, binc*bprop)
  beta[3] <- min(-beta[1]*2, binc*(bprop-1))
  beta[4] ~ dbeta(0.5,0.5)T(0.5,1)


  #Hidden Markov layer definition

  for (i in 1:K) {
    state[i] <- 1
  }
  for (i in (K+1):nweeks) {
    state[i] ~ dcat(
      ifelse((state[i-1]==2) && (state[max(i-3, 1)]!=2),
        c(0,1,0,0,0),
        ifelse((state[i-1]==4) && (state[max(i-3, 1)]!=4),
          c(0,0,0,1,0),
          Pmat[state[i-1], ]
        )
      )
    )
  }

  #Hyperparameters of the hidden layer
  P12 ~ dbeta(0.5, 0.5)
  P23 ~ dbeta(0.5, 0.5)
  P24 ~ dbeta(0.5, 0.5)
  P34 ~ dbeta(0.5, 0.5)
  P45 ~ dbeta(0.5, 0.5)

  Pmat[1,1] <- 1 - P12
  Pmat[1,2] <- P12
  Pmat[1,3] <- 0
  Pmat[1,4] <- 0
  Pmat[1,5] <- 0

  Pmat[2,1] <- 0
  Pmat[2,2] <- 1 - P23 - P24
  Pmat[2,3] <- P23
  Pmat[2,4] <- P24
  Pmat[2,5] <- 0

  Pmat[3,1] <- 0
  Pmat[3,2] <- 0
  Pmat[3,3] <- 1 - P34
  Pmat[3,4] <- P34
  Pmat[3,5] <- 0

  Pmat[4,1] <- 0
  Pmat[4,2] <- 0
  Pmat[4,3] <- 0
  Pmat[4,4] <- 1 - P45
  Pmat[4,5] <- P45

  Pmat[5,1] <- 0
  Pmat[5,2] <- 0
  Pmat[5,3] <- 0
  Pmat[5,4] <- 0
  Pmat[5,5] <- 1

}'


fluModelChunkIsol <- '

  # Isolates part

  for (i in 1:nweeksIsol) {
    nIsol[i] ~ dpois(nTrueIsol[i])
    nTrueIsol[i] ~ dnorm(muIsol[i], tauIsol[state[i]])T(0,)
  }

  sigmaIsol0[1] ~ dunif(0.1, SSDisol)
  sigmaIsol0[2] ~ dunif(0.1, SSDisol)
  sigmaIsol[1:2] <- sort(sigmaIsol0)

  tauIsol[1] <- pow(sigmaIsol[1], -2)
  for (i in 2:4) {
    tauIsol[i] <- pow(sigmaIsol[2], -2)
  }
  tauIsol[5] <- pow(sigmaIsol[1], -2)
  muIsol[1] <- muPreIsol
  for (i in 2:nweeksIsol) {
    muIsol[i] <- max(0, equals(state[i],1)*(muIsol[i-1] + betaIsol[1]) +
        equals(state[i],2)*(muIsol[i-1] + betaIsol[2]) +
        equals(state[i],3)*(muIsol[i-1]) +
        equals(state[i],4)*(muIsol[i-1] + betaIsol[3]) +
        equals(state[i],5)*(muIsol[i-1]*betaIsol[4]))
  }

  muPreIsol ~ dnorm(0, 0.001)T(0,SMAXisol)

  bincIsol ~ dnorm(0, 0.001)T(0,SRisol*2)
  bpropIsol ~ dbeta(2,2)

  betaIsol[1] ~ dnorm(0, 0.001)T(0,SRisol)
  betaIsol[2] <- max(betaIsol[1]*2, bincIsol*bpropIsol)
  betaIsol[3] <- min(-betaIsol[1]*2, bincIsol*(bpropIsol-1))
  betaIsol[4] ~ dbeta(0.5,0.5)T(0.5,1)

  '


