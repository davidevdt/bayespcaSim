# main: Simulation Study ----------------------------------------------------------------------------------
# devtools::install_github("cran/elasticnet")
# devtools::install_github("davidevdt/bayespca")
library(elasticnet)
library(dplyr)
library(bayespca)
library(ggplot2	)


# Load simulation functions ----------------------------------------------------------------------------
source("https://raw.githubusercontent.com/davidevdt/bayespcaSim/master/loadFunctions.R")



# Select final plots ------------------------------------------------------------------------------------
plotN <- 2 									# 1. for Tucker congruence; 2. for % correct zero/nonzeros; 3. for Reconstruction Error

# Simulation Conditions ----------------------------------------------------------------------------------  
nsim <- 30
Icond <- c(25, 50, 100)						# Sample size conditions 
Jcond <- 50									# Number of variables   	
noiseCond <- c(0.05, 0.25)					# Prop. Noise conditions
sparsityCond <- c(0.5, 0.9)					# Prop. Sparsity conditions	

# Componenents' Parameters ------------------------------------------------------------------------------
D = 3										# Number of components   
varComp <- c(200, 100, 50)					# Variance of principal components

# Hyperparameters -------------------------------------------------------------------------
# Tuning parameters (spca)
beta <- c(0., 1e-04, 1e-03, 1e-02, 1e-01, 1, 2, 5, 10, 100)  				# Ridge grid penalty
numFolds <- 5 
sdRule <- TRUE																# Select parameters with S.E. rule


# Tuning parameter (bayesPCA - parameters for InverseGamma prior)
alphaIG <- c(  1, 2, 5, 10, 20 )
betaIG <- c(  1, 2, 4, 5, 10 )


# Hyperparameters (Stochastic Variable Selection)
SVS <- TRUE 							# If SVS == FALSE: use HPD intervals
propSpike <- 1e-4						# proportion of prior 'spike' variance 
priorInclusion <- rep(0.5, D) 			# prior inclusion probabilities 
beta1pi <- 1 							# For uniform Beta prior set beta1pi = 1
beta2pi <- 1 							# For uniform Beta prior set beta2pi = 1
threshold <- 0.50						# Probability threshold to mark elements of W as 0's 

# Other controls ----------------------------------------------------------------------------------   
maxiter <- 1e+05						 
tolerance <- 1e-02						# Convergence criterion -- Set to 1e-02 otherwise oracle elasticnet is too slow with high sparsity 
										# (decrease it for more precise results)
typeTuck <- 2 							# If typeTuck == 2: set to 0 elements with Pr(inclusion) < 0.5 
selType <- 2 							# If selType == 1: work with scaled observed data 
normalise <- FALSE
updatetau <- FALSE 
priorvar <- 'invgamma'
global.var <- FALSE  
useOrig <- TRUE							# If useOrig == FALSE : use weight matrix estimated without SVS 
origElbo <- TRUE						# If origElbo == FALSE : use Elbo computed without SVS
probHPDI <- 0.9							
 
# Run the simulations ----------------------------------------------------------------------------------  
simRes <- runSim( nsim, Icond, Jcond, noiseCond, sparsityCond,
					D, varComp, numFolds, 
					threshold, maxiter, tolerance, 
					typeTuck, selType, propSpike, 
					alphaIG, betaIG, beta, 
					SVS, normalise, beta1pi, beta2pi, 
					updatetau, priorvar, 
					priorInclusion, global.var, sdRule,
					useOrig, origElbo, probHPDI )

### RESULTS --------------------------------------------------------------------------------------------
globalResults <- simRes$globalResults 
globalResults$Method <- factor(globalResults$Method, levels = as.character(unique(globalResults$Method)) )
	

globalResults2 <- globalResults 
globalResults <- globalResults %>% filter(
	Method == "G(1,1)" |		
	Method == "G(10,5)" | 
	Method == "G(20,1)" | 
	Method == "G(5,1)" |  
	Method == "G(20,5)" |
	Method == "spca"
)
globalResults$I <- as.numeric(globalResults$I)
globalAvgMatrix <- simRes$globalAvgMatrix 
allConditions <- simRes$allConditions

				  
				  
				  
				  

# Plot the results ----------------------------------------------------------------------------------  
# Prepare Datasets for ggplot				
aggrRecErr <- groupIndex(globalResults)
	
# Labels for the conditions 
I.labs <- paste0("I = ", Icond)
names(I.labs) <- Icond
sp.labs <- paste0("Sparsity = ", sparsityCond)
names(sp.labs) <- sparsityCond
noise.labs <- paste0("Noise = ", noiseCond)
names(noise.labs) <- noiseCond




# Make plots -----------------------------------------------------------------------------------------------------------
if( plotN == 1 ){

	ggplot(globalResults, aes(x = Method, y = Tucker)) + 
		geom_boxplot(aes(fill = Method), alpha = 0.5) + 
		scale_fill_brewer(palette = "Set1") + 
		facet_grid(I~Noise+Sparsity, 
					labeller = labeller(I = I.labs, Noise = noise.labs, Sparsity = sp.labs )) + 
		ylab("Tucker Congruence") + 
		theme(strip.text.x = element_text(colour = "blue")) +
		theme(strip.text.y = element_text(colour = "blue")) +
		theme(legend.position = "none") +
		scale_x_discrete(position = "bottom") + 
		theme(axis.title.x=element_blank(),
		axis.ticks.x=element_line(size=1), 
		axis.text.x = element_text(face="bold"),
		legend.title = element_blank(), 
		legend.text = element_text(face="plain",size=10)	)+ 
		theme(legend.position = "right") +
		guides(fill = guide_legend(nrow = 7, override.aes = list(size = 0.5))) + 
		ggtitle("Tucker Congruence")+
		theme(plot.title = element_text(hjust = 0.5, face="bold"))

}else if( plotN == 2 ){

	ggplot(globalResults, aes(x = Method, y = PropCorrect)) + 
		geom_boxplot(aes(fill = Method), alpha = 0.5) + 
		scale_fill_brewer(palette = "Set1") + 
		facet_grid(I~Noise+Sparsity, 
					labeller = labeller(I = I.labs, Noise = noise.labs, Sparsity = sp.labs )) + 				
		ylab("Proportion of Correctly Identified Weights") + 
		theme(strip.text.x = element_text(colour = "blue")) +
		theme(strip.text.y = element_text(colour = "blue")) +
		theme(legend.position = "none") +
			scale_x_discrete(position = "bottom") + 
		theme(axis.title.x=element_blank(),
		axis.ticks.x=element_line(size=1), 
		axis.text.x = element_text(face="bold"),
		legend.title = element_blank(), 
		legend.text = element_text(face="plain",size=10)	)+ 
		theme(legend.position = "right") +
		guides(fill = guide_legend(nrow = 7, override.aes = list(size = 0.5))) + 
		ggtitle("Proportion of Correctly Identified Weights")+
		theme(plot.title = element_text(hjust = 0.5, face="bold"))
		
}else{

	ggplot(globalResults, aes(x = Method, y = RecErr)) + 
		geom_boxplot(aes(fill = Method), alpha = 0.5) + 
		scale_fill_brewer(palette = "Set1") + 
		# geom_line(data=aggrRecErr, aes(x=Method, y=RecErr, group=1), alpha = 0.7,col = "deepskyblue", size = 1) + 
		# geom_point(data=aggrRecErr, aes(x=Method, y=RecErr), alpha = 0.8, size = 2, color = "deepskyblue") + 
		# scale_colour_manual(values = c("deepskyblue"), labels=c( 'Average Reconstruction Error') ) +
		facet_grid(I~Noise+Sparsity, 
					labeller = labeller(I = I.labs, Noise = noise.labs, Sparsity = sp.labs ), scales = "free_y") + 				
		ylab("Rel. Reconstruction Error") + 
		theme(strip.text.x = element_text(colour = "blue")) +
		theme(strip.text.y = element_text(colour = "blue")) +
		theme(legend.position = "none") +
		scale_x_discrete(position = "bottom") + 
		theme(axis.title.x=element_blank(),
		axis.ticks.x=element_line(size=1), 
		axis.text.x = element_text(face="bold"),
		legend.title = element_blank(), 
		legend.text = element_text(face="plain",size=10)	)+ 
		theme(legend.position = "right") +
		guides(fill = guide_legend(nrow = 7, override.aes = list(size = 0.5))) + 
		ggtitle("Relative Reconstruction Error")+
		theme(plot.title = element_text(hjust = 0.5, face="bold"))
		
}
