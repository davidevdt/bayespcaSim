# main: Simulation Study ----------------------------------------------------------------------------------
# devtools::install_github("erichson/spca")
# devtools::install_github("davidevdt/bayespca")
library(sparsepca)
library(dplyr)
library(bayespca)
library(ggplot2	)



# Load simulation functions ----------------------------------------------------------------------------
source("https://raw.githubusercontent.com/davidevdt/bayespcaSim/master/loadFunctions.R")





# Select final plots ------------------------------------------------------------------------------------
plotN <- 2 									# 1. For Tucker congruence; 2. For % correct zero/nonzeros

# Simulation Conditions ----------------------------------------------------------------------------------  
nsim <- 30
Icond <- c(25, 50, 100)						# Sample size conditions 
Jcond <- 50									# Number of variables   	
noiseCond <- c(0.1, 0.3)					# Prop. Noise conditions
sparsityCond <- c(0.3, 0.8)					# Prop. Sparsity conditions	



# Componenents' Parameters ------------------------------------------------------------------------------
D = 3										# Number of components   
varComp <- c(200, 100, 50)					# Variance of principal components




# Hyperparameters -------------------------------------------------------------------------
# Tuning parameters (spca)
alpha <- c(0.,.0001, .001, .01, .1, .5, .9, 1., 2.)					# Lasso 
beta <- c(0., .0001, .001, .01, .1, .5, .9, 1)						# Ridge 
numFolds <- 5 
sdRule <- TRUE														# Select parameters with S.E. rule


# Tuning parameter (bayesPCA - parameters for InverseGamma prior)
alphaIG <- c(  1, 5, 10, 20, 50 )
betaIG <- c(  1, 5, 10, 20)



# Hyperparameters (Stochastic Variable Selection)
SVS <- TRUE 							# If SVS == FALSE: use HPD intervals
propSpike <- 1e-04						# proportion of prior 'spike' variance 
priorInclusion <- rep(0.5, D) 			# prior inclusion probabilities 
beta1pi <- 1 
beta2pi <- 1 
threshold <- 0.50						# Probability threshold to mark elements of W as 0's 




# Other controls ----------------------------------------------------------------------------------   
maxiter <- 1e+05						 
tolerance <- 1e-05		 
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
					alphaIG, betaIG, alpha, beta, 
					SVS, normalise, beta1pi, beta2pi, 
					updatetau, priorvar, 
					priorInclusion, global.var, sdRule,
					useOrig, origElbo, probHPDI )
				  
				  


				  
				  
				  
				  
				  
				  
				  
				  

# Plot the results ----------------------------------------------------------------------------------  
# Prepare Datasets for ggplot				
filterMethods <- TRUE

globalResults <- simRes$globalResults 
globalResults$Method <- factor(globalResults$Method, levels = as.character(unique(globalResults$Method)) )

if( filterMethods ){
	globalResults2 <- globalResults 
	globalResults <- globalResults %>% filter(
		Method == "IG(1,1)" | 
		Method == "IG(10,5)" | 
		Method == "IG(20,5)" | 
		Method == "IG(5,1)" |  
		Method == "IG(50,20)" |
		Method == "spca"
	)
}



globalAvgMatrix <- simRes$globalAvgMatrix 
allConditions <- simRes$allConditions



aggrElbos <- groupIndex(globalResults)
aggrElbos <- adjustAggregated(aggrElbos0b, normalize = TRUE )
		
	
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
		geom_line(data=aggrElbos, aes(x=Method, y=RecErr,group=1,color=Type), alpha = 0.8,col = "deepskyblue", size = 1) + 
		geom_point(data=aggrElbos, aes(x=Method, y=RecErr,color=Type), alpha = 0.8, size = 2) + 
		 scale_colour_manual(values = c("deepskyblue"), labels=c( 'Reconstruction Error') ) +
		facet_grid(I~Noise+Sparsity, 
					labeller = labeller(I = I.labs, Noise = noise.labs, Sparsity = sp.labs )) + 
		ylab("Tucker Congruence") + 
		theme(strip.text.x = element_text(face = "bold", colour = "blue")) +
		theme(strip.text.y = element_text(face = "bold", colour = "blue")) +
		theme(legend.position = "none") +
		scale_x_discrete(position = "top") + 
		theme(axis.title.x=element_blank(),
		axis.ticks.x=element_line(size=1), 
		axis.title = element_text(face = "bold"),
		legend.title = element_blank(), 
		legend.text = element_text(face="bold",size=10)	)+ 
		theme(legend.position = "bottom") +
		guides(fill = guide_legend(nrow = 1, override.aes = list(size = 1))) + 
		ggtitle("Tucker Congruence")+
		theme(plot.title = element_text(hjust = 0.5, face="bold"))

}else{

	ggplot(globalResults, aes(x = Method, y = PropCorrect)) + 
		geom_boxplot(aes(fill = Method), alpha = 0.5) + 
		scale_fill_brewer(palette = "Set1") + 
		geom_line(data=aggrElbos, aes(x=Method, y=RecErr,group=1,color=Type), alpha = 0.8,col = "deepskyblue", size = 1) + 
		geom_point(data=aggrElbos, aes(x=Method, y=RecErr,color=Type), alpha = 0.8, size = 2) + 
		 scale_colour_manual(values = c("deepskyblue"), labels=c( 'Reconstruction Error') ) +
		facet_grid(I~Noise+Sparsity, 
					labeller = labeller(I = I.labs, Noise = noise.labs, Sparsity = sp.labs )) + 				
		ylab("Proportion of Correctly Identified Weights") + 
		theme(strip.text.x = element_text(face = "bold", colour = "blue")) +
		theme(strip.text.y = element_text(face = "bold", colour = "blue")) +
		theme(legend.position = "none") +
		scale_x_discrete(position = "top") + 
		theme(axis.title.x=element_blank(),
		axis.ticks.x=element_line(size=1), 
		axis.title = element_text(face = "bold"),
		legend.title = element_blank(), 
		legend.text = element_text(face="bold",size=10)	)+ 
		theme(legend.position = "bottom") +
		guides(fill = guide_legend(nrow = 1, override.aes = list(size = 1))) + 
		ggtitle("Proportion of Correctly Identified Weights")+
		theme(plot.title = element_text(hjust = 0.5, face="bold")) 


}

				
