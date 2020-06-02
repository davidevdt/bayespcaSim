# runSim: read simulation parameters and create conditions  --------------------------------------------------------------
runSim <- function( nsim, Icond, Jcond, noiseCond, sparsityCond,
					D, varComp, numFolds, 
					threshold, maxiter, tolerance_elastic_net,
					tolerance_vbpca, 
					typeTuck, selType,
					alphaIG, betaIG, beta, 
					normalise, updatetau,
					global.var, sdRule,
					useOrig, origElbo, hpdi, probHPDI ){
				  					
				
	# Simulation conditions 
	allConditions <- expand.grid(I = Icond, J = Jcond, percNoise = noiseCond,
								 percSpars = sparsityCond)
	allConditions <- allConditions[order(allConditions$I),]
	row.names(allConditions) <- 1:nrow(allConditions)

	# spca and InverseGamma hyperparameters 
	spcaPars <- expand.grid(beta = beta)
	vbpcaPars <- expand.grid(alpha = alphaIG, beta = betaIG)

	# Zero/NonZero matrix 
	zeroMatList <- list()

	set.seed(1)
	for(sp in 1:length(sparsityCond)){
		zeroMatList[[as.character(sparsityCond[sp])]] <- list()
		for( j in 1:length(Jcond)  ){
				nameList <- as.character(Jcond[j])
				zeroMatList[[as.character(sparsityCond[sp])]][[nameList]] <- genZeroMat(Jcond[j], D, sparsityCond[sp])
		}
	}



	# Results (all conditions)
	globalResults <- data.frame()					  
	globalAvgMatrix <- list()


	startingTime <- proc.time()[3]
	
	# Run simulations 
	for( cond in 1:nrow(allConditions)  ){
		print( sprintf('######################################################') ) 
		print( sprintf('############## CONDITION: %1d of %1d  ###################', cond, nrow(allConditions) ) )
		print( sprintf('######################################################') ) 
		cat("\n")
		
		I <- allConditions[cond,]$I
		J <- allConditions[cond,]$J
		percNoise <- allConditions[cond,]$percNoise 
		percSpars <- allConditions[cond,]$percSpars
		zeroMat <- zeroMatList[[as.character(percSpars)]][[as.character(J)]]
		
		
		
		set.seed(42)
		####################  Run the simulations and get the results ######################### 
		res <- estimateModels( nsim, I, J, D, zeroMat, percNoise, varComp, 
					   vbpcaPars, spcaPars, numFolds, 
					   threshold, maxiter, tolerance_elastic_net,
					   tolerance_vbpca, 
					   typeTuck, selType,
					   normalise, updatetau,  
					   global.var, sdRule,
					   useOrig, origElbo, hpdi, probHPDI )
						 	 
		
		# Condition names  
		sampleSize <- rep( as.character(I), nrow( res[[1]] ) )
		nVariables <- rep( as.character(J), nrow( res[[1]] )  )
		sparsityLevel <- rep( as.character(percSpars), nrow( res[[1]] ) )
		noiseLevel <- rep( as.character(percNoise), nrow( res[[1]] ) )

		makeDf <- data.frame( "I" = sampleSize, "J" = nVariables,
							  "Sparsity" = sparsityLevel, "Noise" = noiseLevel, check.names = FALSE )
							  
		tmpDf <- cbind(makeDf, res[[1]])

		globalResults <- rbind(globalResults, tmpDf)
		globalAvgMatrix[[cond]] <- res[[2]]
			

	}
	print( sprintf('Total Simulation Time: %4.2e', proc.time()[3]-startingTime  ) )
	
	return(list(
		globalResults = globalResults, 
		globalAvgMatrix = globalAvgMatrix, 
		zeroMatList = zeroMatList, 
		allConditions = allConditions
	))

}
