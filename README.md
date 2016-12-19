### FluHMM - Hidden Markov Model for influenza sentinel surveillance

This R package is meant to be used in the context of seasonal influenza ILI/ARI sentinel 
surveillance. The usual influenza surveillance season is from week 40 to week 20 of the
next year, and for this period every year, the sentinel surveillance system collects weekly
Influenza-Like Illness / Acute Respiratory Infection (ILI/ARI) rates from a network of 
(usually primary-care) physicians.

The package allows segmenting this set of weekly rates to five phases: (1) pre-epidemic, (2)
epidemic growth, (3) epidemic plateau, (4) epidemic decline, and (5) post-epidemic. This is done
by fitting a Bayesian Hidden Markov Model on the rates, and calculating the posterior
probabilities of each phase per week, using only the available data.

This way an alert can be raised when the posterior probability of the epidemic growth phase rises
beyond a certain threshold, signifying that influenza activity in the community is now rising.
This method uses no historical data other than the previous weeks of the current season. (In the
future it is planned to take account of historical seasons via appropriate priors).

The parametrization of the model and an evaluation of its performance are planned to be described 
in a peer-reviewed scientific article. For now, anyone interested can study the package code.

Besides an R installation, you need to have JAGS installed and the rjags package. These are the
only dependencies of this package.

**Installation**

Open R, and give:

      devtools::install_git("https://github.com/thlytras/FluHMM.git")

If you do not have the package "devtools", first install it from CRAN with:

      install.packages("devtools")

