# Bayesian-ARFIMA-Estimation
Code for the paper ARFIMA Bayesian Estimation using MCMC and ABC methods. The programs can be run via command line using the following arguments:

MCMC:
Rscript [true d value] [true phi value] [true theta value] [num sims] [repetition start id] [repetition end id]

ABC:
Rscript [true d value] [true phi value] [true theta value] [num sims] [repetition start id] [repetition end id]

The output location for the simulation files can be set using the *location* parameter at the top of the code.
Additional parameters such as the innovation variance, series length, and ABC quantiles can also be set at the top of the code.
