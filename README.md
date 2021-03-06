# bayespcaSim
### Simulation Study for the Bayes PCA method

The method (Bayes PCA) is compared with the sparse PCA specification of the model. Simulations are run with the R (>= 3.3.0) programming language. Bayesian PCA is implemented via the ```bayespca``` package by D. Vidotto (https://github.com/davidevdt/bayespca), while Sparse PCA is implemented with the ```elasticnet``` package by H. Zou (https://github.com/cran/elasticnet). 


### Run the simulations
To run the simulations: 
 1. Install the required R packages: 
     * ```devtools::install_github("cran/elasticnet")```
     * ```devtools::install_github("davidevdt/bayespca")```
 2. Launch ```main.R``` ; in this file, simulation parameters and plotting functions can be specified 
     * modify the simulation parameters by changing the values that appear before ```runSim()```
     * select the type of results you want to visualize: ```plotN = 1``` for Tucker congruence, ```plotN = 2``` for proportion of correct zeros/nonzeros, ```plotN = 3``` for the reconstruction errors
