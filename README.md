# Bayesian-ARFIMA-Estimation
Code for the paper ARFIMA Bayesian Estimation using MCMC and ABC methods. The programs can be run via command line using the following arguments:

MCMC:
Rscript mcmc_simulation.R [true d value] [true phi value] [true theta value] [num sims] [repetition start id] [repetition end id]

ABC:
Rscript abc_simulation.R [true d value] [true phi value] [true theta value] [num sims] [repetition start id] [repetition end id]

The first 3 arguments correspond to the true parameter values used in the simulated ARFIMA series. the *num sims* argument is the number of MCMC or ABC simulations to perform. *repetition start id* and *repetition end id* represent the start and end id's of for a sequence of simulated ARFIMA series. 

**Example**
Rscript mcmc_simulation.R 0.2 0.5 0.2 5000 1 5

The above code will run the MCMC code for parameters $d=0.2, \phi=0.5, \theta=0.5$ with $5000$ simulations. This will be repeated for 5 simulated true ARFIMA series, with ID's 1,2,3,4,5.

The output location for the simulation files can be set using the *location* parameter at the top of the code.

Additional parameters such as the innovation variance, series length, and ABC quantiles can also be set at the top of the code in both files.

The output consists of a file for each simulation, model, and repetition. Each file contains the posterior draws at a given simulation step. MCMC output also includes the proposal values and their respective posterior values. 
