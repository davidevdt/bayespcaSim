# runModels: estimate PCA models with input condition --------------------------------------------------------------
# Data are generated during the simulation iterations 
estimateModels <- function( nsim, I, J, D, zeroMat, percNoise, varComp, 
							   vbpcaPars, spcaPars, numFolds, 
							   threshold, maxiter, tolerance_elastic_net,
							   tolerance_vbpca, 
							   typeTuck, selType,
							   normalise, updatetau,  
							   global.var, sdRule,
							   useOrig, origElbo, hpdi, probHPDI ){
					   	  
	# Store the results 
	Tucker <- matrix(0, nrow = (nrow(vbpcaPars) + 1), ncol = nsim ) 
	PercAll <- matrix(0, nrow = (nrow(vbpcaPars) + 1), ncol = nsim ) 
	PercZeros <- matrix(0, nrow = (nrow(vbpcaPars) + 1), ncol = nsim ) 
	PercOnes <- matrix(0, nrow = (nrow(vbpcaPars) + 1), ncol = nsim ) 
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
		
		
		
		######################################## spca [sparsepca] ########################################		--> Not run 
		# print( sprintf('### 1. spca model [sparsepca] ###') )
		# Tuning
		# print( sprintf('Tuning...') )	
		
		# if( sdRule == FALSE ){
		# 	minMse <- Inf 
		# 	selPar <- NULL
			
			
		#	for( i in 1:nrow(spcaPars)  ){
		#	
		#		set.seed(71)
		#		newMSE <- EigenCV(X = Xsel, alpha = spcaPars[i, ]$alpha, 
		#						  beta = spcaPars[i, ]$beta, nFolds = numFolds,
		#						  D = D, maxIt = maxiter, tol = tolerance_elastic_net, scaleDat = FALSE,
		#						  verbose = FALSE)$MSE
		#		if(newMSE < minMse){
		#			minMse <- newMSE 
		#			selPar <- i 
		#		}	
		#	}
		#	cat("Selected tuning set: ", selPar, "\n" )
		#	alphaSel <- spcaPars[selPar, "alpha"]
		#	betaSel <- spcaPars[selPar, "beta"]
		# }else{
		#	mse <- rep(0, nrow(spcaPars))
		#	minMse <- Inf
		#	sdMse <- 0 
		#	selPar <- NULL
		#	
		#	
		#	for( i in 1:nrow(spcaPars)  ){
		#	
		#		set.seed(71)
		#		eCV <- EigenCV(X = Xsel, alpha = spcaPars[i, ]$alpha, 
		#						  beta = spcaPars[i, ]$beta, nFolds = numFolds,
		#						  D = D, maxIt = maxiter, tol = tolerance_elastic_net, scaleDat = FALSE,
		#						  verbose = FALSE)
		#		mse[i] <- eCV$MSE
		#		
		#		if(mse[i] < minMse){
		#			minMse <- mse[i]
		#			sdMse <- eCV$sdMSE
		#		}	
		#						
		#	}
		#	
		#	indx <- which(mse <= minMse + sdMse )
		#	selPar <- max(indx)
		#	
		#	cat("Selected tuning set: ", selPar, "\n" )
		#	alphaSel <- spcaPars[selPar, "alpha"]
		#	betaSel <- spcaPars[selPar, "beta"]		
		#
		#
		# }
		#
		# Estimation 
		# print( sprintf('Estimation') )
		# set.seed(71)
		# mod <- sparsepca::spca(X = Xsel, k = D, alpha = alphaSel, beta = betaSel,
		#					   center = FALSE, scale = FALSE, max_iter = maxiter, 
		#					   tol = tolerance_elastic_net, verbose = FALSE )
        #
        #
		# Tucker[1,b] <- TuckerCoef(Wmat, mod$loadings)$tucker_value  
		# PercZeros[1,b] <- sum( mod$loadings[ind0] == 0 ) / length(ind0) 
		# PercOnes[1,b] <- sum( mod$loadings[ind1] != 0 ) / length(ind1)
		# PercAll[1,b] <- sum( ( (mod$loadings != 0 ) * 1 ) == zeroMat ) / (J * D)
		#
		# RECERR[1,b] <- recErr(Xsel, mod$loadings, mod$transform)
		# 
		# MODEL[1, b] <- "spca1"
		# avgMatrix[[1]] <- avgMatrix[[1]] + ( (mod$loadings != 0 ) * 1 ) 
		
		
		
		
		######################################## spca [elasticnet] ########################################
		# Note: oracle method; CV implemented to find the optimal Ridge parameter
		
		print( sprintf('### 1. spca model [elasticnet] ###') )	
		
		
		n_loadings <- apply(zeroMat, 2, sum) 
		
		# Tuning
		print( sprintf('Tuning...') )	
		
		beta_grid <- unique( spcaPars$beta )		
		
		if( sdRule == FALSE ){
			minMse <- Inf 
			selPar <- NULL
					
					
			for( i in 1:length(beta_grid) ){
					
				set.seed(71)
				newMSE <- EigenCVel_net(X = Xsel, alpha = n_loadings, 
								  beta = beta_grid[i], nFolds = numFolds,
								  D = D, maxIt = maxiter, tol = tolerance_elastic_net, scaleDat = FALSE,
								  verbose = FALSE)$MSE
				if(newMSE < minMse){
					minMse <- newMSE 
					selPar <- i 
				}	
			}
			cat("Selected tuning set: ", selPar, "\n" )
			betaSel <- beta_grid[selPar]

		}else{
			mse <- rep(0, length(beta_grid))
			minMse <- Inf
			sdMse <- 0 
			selPar <- NULL
					
					
			for( i in 1:length(beta_grid) ){
					
				set.seed(71)
				eCV <- EigenCVel_net(X = Xsel, alpha = n_loadings, 
								  beta = beta_grid[i], nFolds = numFolds,
								  D = D, maxIt = maxiter, tol = tolerance_elastic_net, scaleDat = FALSE,
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
			betaSel <- beta_grid[selPar]		
				
				
		}
		
			
		# Estimation 
		print( sprintf('Estimation') )
		set.seed(71)

		mod <- elasticnet::spca(x = Xsel, K = D, para = n_loadings, lambda = betaSel,
								sparse = "varnum", max.iter = maxiter, 
							   eps.conv = tolerance_elastic_net, trace = FALSE, 
							   type = "predictor", use.corr = FALSE)	

		XTX <- t(Xsel) %*% Xsel
		sVd <- svd(XTX %*% mod$loadings)
		P <- sVd$u %*% sVd$v		

		Tucker[1,b] <- TuckerCoef(Wmat, mod$loadings)$tucker_value  
		PercZeros[1,b] <- sum( mod$loadings[ind0] == 0 ) / length(ind0) 
		PercOnes[1,b] <- sum( mod$loadings[ind1] != 0 ) / length(ind1)
		PercAll[1,b] <- sum( ( (mod$loadings != 0 ) * 1 ) == zeroMat ) / (J * D)
		RECERR[1,b] <- recErr(Xsel, mod$loadings, P)
		MODEL[1, b] <- "spca"
		avgMatrix[[1]] <- avgMatrix[[1]] + ( (mod$loadings != 0 ) * 1 ) 
			


	
		
		######################################## vbpca  ########################################
		print( sprintf('### 2. vbpca model: ###') )		
		
		for( i in 1:nrow(vbpcaPars) ){
			print( sprintf('   par %1d of %1d', i, nrow(vbpcaPars) ) )
				
			alphaGamma <- vbpcaPars[i, 1]
			betaGamma <- vbpcaPars[i, 2]
				
			set.seed(71)
			mod <- bayespca::vbpca( X = Xsel, D = D, maxIter = maxiter,
									center = FALSE, scalecorrection = -1, 
									svdStart = TRUE, normalise = normalise,
									seed = 71, plot.lowerbound = FALSE, 
									hpdi = hpdi, probHPDI = probHPDI,
									alphatau = alphaGamma, betatau = betaGamma,
									tolerance = tolerance_vbpca, verbose = FALSE, 
									tau = 1, updatetau = updatetau,									
									 global.var = global.var, 
									 control = ctrl, suppressWarnings = TRUE )							
			
			tuckMat <- mod[[1]]
			
			if( hpdi ){
				estZeroMat <- zeroHPDI( mod[[5]] )
			}else{
				estZeroMat <- matrix(1, nrow(tuckMat), ncol(tuckMat)) 
				estZeroMat[plotheatmap(mod, matrix_type = "W", bound_tau = threshold)$W == 0] <- 0 			
			}
				
			if( typeTuck == 2 ){
				tuckMat[estZeroMat == 0] <- 0
			}
				
			Tucker[(i+1),b] <- TuckerCoef(Wmat, tuckMat)$tucker_value  
			PercZeros[(i+1),b] <- sum( estZeroMat == 0 ) / length(ind0) 
			PercOnes[(i+1),b] <- sum( estZeroMat != 0 ) / length(ind1)
			PercAll[(i+1),b] <- sum( ( (estZeroMat != 0 ) * 1 ) == zeroMat ) / (J * D)
				
			RECERR[(i+1),b] <- recErr(Xsel, tuckMat, mod$P)

			MODEL[(i+1), b] <- paste0("IG(",alphaGamma,",",betaGamma,")")		
			avgMatrix[[i+1]] <- avgMatrix[[i+1]] + ( ( estZeroMat != 0 ) * 1 ) 
			
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
						'RecErr' = c( t(RECERR) )
					)
					  
						
	return( list(retDataFrame, avgMatrix ) )


}





