# runModels: estimate PCA models with input condition --------------------------------------------------------------
# Data are generated during the simulation iterations 
estimateModels <- function( nsim, I, J, D, zeroMat, percNoise, varComp, 
					   vbpcaPars, spcaPars, numFolds, 
					   threshold, maxiter, tolerance, 
					   typeTuck = 1, selType = 2,
					   propSpike, SVS, normalise, 
					   beta1pi, beta2pi, 
					   updatetau, priorvar, priorInclusion, global.var,
					   sdRule, useOrig = TRUE, origElbo = TRUE, probHPDI){ 
					   	  
	# Store the results 
	Tucker <- matrix(0, nrow = (nrow(vbpcaPars) + 1), ncol = nsim ) 
	PercAll <- matrix(0, nrow = (nrow(vbpcaPars) + 1), ncol = nsim ) 
	PercZeros <- matrix(0, nrow = (nrow(vbpcaPars) + 1), ncol = nsim ) 
	PercOnes <- matrix(0, nrow = (nrow(vbpcaPars) + 1), ncol = nsim ) 
	ELBOS <- matrix(0, nrow = (nrow(vbpcaPars) + 1), ncol = nsim )
	RECERR <- matrix(0, nrow = (nrow(vbpcaPars) + 1), ncol = nsim )
	MODEL <- matrix(0, nrow = (nrow(vbpcaPars) + 1), ncol = nsim )
	avgMatrix <- rep(list( matrix(0, J, D) ), (nrow(vbpcaPars) + 1) )  
	
	ind0 <- which(zeroMat == 0)
	ind1 <- which(zeroMat != 0)

					  
	for( b in 1:nsim ){
	
		print( sprintf('##### Working on replication: %1d of %1d #####', b, nsim) )
	
		######################################## Data Generation ########################################
		# --> 1. Generate loading matrix 
		set.seed(387 + b)	
		gl <- genLoading(zeroMat, maxIt = 1e+04)
		
		P <- gl[[1]]
		if(gl$converged == 0){
			stop("Loading matrix not created.")
		}
		
		# --> 2. Generate data 
		set.seed(388 + b)
		dg <- datageneration( I, P, percNoise, varComp )
		Xobs <- dg[[1]]
		Wmat <- dg[[3]]		
		
		

		
		# Select observed/scaled data 
		if( selType == 1 ){
			
			# Scale data 
			Xscaled <- scaleX(Xobs, center = TRUE, scalingFactor = 0 )
			attr(Xscaled, "center") <- NULL
			attr(Xscaled, "scaling") <- NULL		
		
			Xsel <- Xscaled 
		
		}else{
			Xsel <- Xobs
		}
		
		
		
		######################################## spca ########################################
		print( sprintf('### 1. spca model ###') )
		# Tuning
		print( sprintf('Tuning...') )	
		
		if( sdRule == FALSE ){
			minMse <- Inf 
			selPar <- NULL
			
			
			for( i in 1:nrow(spcaPars)  ){
			
				set.seed(71)
				newMSE <- EigenCV(X = Xsel, alpha = spcaPars[i, ]$alpha, 
								  beta = spcaPars[i, ]$beta, nFolds = numFolds,
								  D = D, maxIt = maxiter, tol = tolerance, scaleDat = FALSE,
								  verbose = FALSE)$MSE
				if(newMSE < minMse){
					minMse <- newMSE 
					selPar <- i 
				}	
			}
			cat("Selected tuning set: ", selPar, "\n" )
			alphaSel <- spcaPars[selPar, "alpha"]
			betaSel <- spcaPars[selPar, "beta"]
		}else{
			mse <- rep(0, nrow(spcaPars))
			minMse <- Inf
			sdMse <- 0 
			selPar <- NULL
			
			
			for( i in 1:nrow(spcaPars)  ){
			
				set.seed(71)
				eCV <- EigenCV(X = Xsel, alpha = spcaPars[i, ]$alpha, 
								  beta = spcaPars[i, ]$beta, nFolds = numFolds,
								  D = D, maxIt = maxiter, tol = tolerance, scaleDat = FALSE,
								  verbose = FALSE)
				mse[i] <- eCV$MSE
				
				if(mse[i] < minMse){
					minMse <- mse[i]
					sdMse <- eCV$sdMSE
				}	
								
			}
			
			indx <- which(mse <= minMse + sdMse )
			selPar <- max(indx)
			
			cat("Selected tuning set: ", selPar, "\n" )
			alphaSel <- spcaPars[selPar, "alpha"]
			betaSel <- spcaPars[selPar, "beta"]		
		
		
		}
		
		# Estimation 
		print( sprintf('Estimation') )
		set.seed(71)
		mod <- sparsepca::spca(X = Xsel, k = D, alpha = alphaSel, beta = betaSel,
							   center = FALSE, scale = FALSE, max_iter = maxiter, 
							   tol = tolerance, verbose = FALSE )


		Tucker[1,b] <- TuckerCoef(Wmat, mod$loadings)$tucker_value  
		PercZeros[1,b] <- sum( mod$loadings[ind0] == 0 ) / length(ind0) 
		PercOnes[1,b] <- sum( mod$loadings[ind1] != 0 ) / length(ind1)
		PercAll[1,b] <- sum( ( (mod$loadings != 0 ) * 1 ) == zeroMat ) / (J * D)
		ELBOS[1,b] <- NA 
		
		RECERR[1,b] <- recErr(Xsel, mod$loadings, mod$transform)
		
		MODEL[1, b] <- "spca"
		avgMatrix[[1]] <- avgMatrix[[1]] + ( (mod$loadings != 0 ) * 1 ) 
		
		
		
		
		
		######################################## vbpca  ########################################
		print( sprintf('### 2. vbpca model: ###') )		
		
		if( SVS ){
		
			for( i in 1:nrow(vbpcaPars) ){
				print( sprintf('   par %1d of %1d', i, nrow(vbpcaPars) ) )
				
				alphaInvGamma <- vbpcaPars[i, 1]
				betaInvGamma <- vbpcaPars[i, 2]
				
				ctrl <- bayespca::vbpca_control(center = FALSE, scale. = FALSE,
								   svdStart = TRUE, normalise = normalise, 
								   seed = 71, plot.lowerbound = FALSE,
								   hpdi = FALSE, probHPDI = 0.9, 
								   alphatau = alphaInvGamma, betatau = betaInvGamma,
								   beta1pi = beta1pi, beta2pi = beta2pi, 
								   v0 = propSpike )	


				set.seed(71)
				mod <- bayespca::vbpca( X = Xsel, D = D, maxIter = maxiter, 
									 tolerance = tolerance, verbose = FALSE, 
									 tau = 1, updatetau = updatetau, 
									 priorvar = priorvar, SVS = SVS, 
									 priorInclusion = priorInclusion, 
									 global.var = global.var, 
									 control = ctrl, suppressWarnings = FALSE )							   
				
				
				probZeroMat <- mod[[9]]				
				
				
				if( useOrig ){
					tuckMat <- mod[[1]]
				}else{
				
					ctrl <- bayespca::vbpca_control(center = FALSE, scale. = FALSE,
								   svdStart = TRUE, normalise = normalise, 
								   seed = 71, plot.lowerbound = FALSE,
								   hpdi = FALSE, probHPDI = 0.9, 
								   alphatau = alphaInvGamma, betatau = betaInvGamma,
								   beta1pi = beta1pi, beta2pi = beta2pi, 
								   v0 = propSpike )	


					set.seed(71)
					modB <- bayespca::vbpca( X = Xsel, D = D, maxIter = maxiter, 
									 tolerance = tolerance, verbose = FALSE, 
									 tau = 1, updatetau = updatetau, 
									 priorvar = priorvar, SVS = FALSE, 
									 priorInclusion = priorInclusion, 
									 global.var = global.var, 
									 control = ctrl, suppressWarnings = FALSE )		

					tuckMat <- modB[[1]]
				
				}
				
				if( typeTuck == 2 ){
					tuckMat[probZeroMat <= threshold] <- 0  
				}
				
				Tucker[(i+1),b] <- TuckerCoef(Wmat, tuckMat)$tucker_value  
				PercZeros[(i+1),b] <- sum( probZeroMat[ind0] <= threshold ) / length(ind0) 
				PercOnes[(i+1),b] <- sum( probZeroMat[ind1] > threshold ) / length(ind1)
				PercAll[(i+1),b] <- sum( ( (probZeroMat > threshold ) * 1 ) == zeroMat ) / (J * D)
				
				
				
				RECERR[(i+1),b] <- recErr(Xsel, tuckMat, mod$P)
				
				
				if( !origElbo & !useOrig ){
					ELBOS[(i+1),b] <- modB$elbo 				
				}else{
					ELBOS[(i+1),b] <- mod$elbo 
				}

				MODEL[(i+1), b] <- paste0("IG(",alphaInvGamma,",",betaInvGamma,")")		
				avgMatrix[[i+1]] <- avgMatrix[[i+1]] + ( (probZeroMat > threshold ) * 1 ) 
			
			}
			
		}else{
		
			for( i in 1:nrow(vbpcaPars) ){
				print( sprintf('   par %1d of %1d', i, nrow(vbpcaPars) ) )
				
				alphaInvGamma <- vbpcaPars[i, 1]
				betaInvGamma <- vbpcaPars[i, 2]
				
				ctrl <- bayespca::vbpca_control(center = FALSE, scale. = FALSE,
								   svdStart = TRUE, normalise = normalise, 
								   seed = 71, plot.lowerbound = FALSE,
								   hpdi = TRUE, probHPDI = probHPDI, 
								   alphatau = alphaInvGamma, betatau = betaInvGamma,
								   beta1pi = beta1pi, beta2pi = beta2pi, 
								   v0 = propSpike )	


				set.seed(71)
				mod <- bayespca::vbpca( X = Xsel, D = D, maxIter = maxiter, 
									 tolerance = tolerance, verbose = FALSE, 
									 tau = 1, updatetau = updatetau, 
									 priorvar = priorvar, SVS = SVS, 
									 priorInclusion = priorInclusion, 
									 global.var = global.var, 
									 control = ctrl, suppressWarnings = FALSE )							   
				
				tuckMat <- mod[[1]]
				estZeroMat <- zeroHPDI( mod[[5]] )

				
				if( typeTuck == 2 ){
					tuckMat[estZeroMat == 0] <- 0  
				}
				
				Tucker[(i+1),b] <- TuckerCoef(Wmat, tuckMat)$tucker_value  
				PercZeros[(i+1),b] <- sum( estZeroMat == 0 ) / length(ind0) 
				PercOnes[(i+1),b] <- sum( estZeroMat != 0 ) / length(ind1)
				PercAll[(i+1),b] <- sum( ( (estZeroMat != 0 ) * 1 ) == zeroMat ) / (J * D)
				
				RECERR[(i+1),b] <- recErr(Xsel, tuckMat, mod$P)

				ELBOS[(i+1),b] <- mod$elbo 
				MODEL[(i+1), b] <- paste0("IG(",alphaInvGamma,",",betaInvGamma,")")		
				avgMatrix[[i+1]] <- avgMatrix[[i+1]] + ( ( estZeroMat != 0 ) * 1 ) 
			
			}	
		
		
		}
	
	}
					  
	avgMatrix <- lapply( avgMatrix, function(x) x / nsim )	
	names(avgMatrix) <- MODEL[,1]

	retDataFrame <- data.frame( 
						'Method' = c( t(MODEL) ), 
						'Tucker' = c( t(Tucker) ), 
						'PropCorrect' = c( t(PercAll) ), 
						'PercZeros' = c( t(PercZeros) ), 
						'PercOnes' = c( t(PercOnes) ),
						'RecErr' = c( t(RECERR) ),
						'ELBO' = c( t(ELBOS)  )	
					)
					  
						
	return( list(retDataFrame, avgMatrix ) )



}





