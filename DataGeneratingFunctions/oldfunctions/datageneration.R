# Data Generating Function (version 1)
datageneration <- function( I, loadMat, noiseProp, varComps = NULL, maxVar = 50 ){

	J <- nrow(loadMat)
	D <- ncol(loadMat)

	if( is.null(varComps) ){
		varComps <- sort(sample(1:maxVar, D, rep = F) , decreasing = TRUE )
	}else{
		if( length(varComps) != D ){
			stop("varComps must be a D-dimensional vector.")
		}
	}
	
	# Generate components matrix 
	Z <- MASS::mvrnorm(n = I, mu = rep(0, D), Sigma = diag(varComps), empirical = TRUE )
	
	# Generate signal matrix (de-noised)
	Xsignal <- Z %*% t(loadMat) 
	SStrue <- var(as.vector(Xsignal))
	# True weight matrix 
	TrueW <- svd(Xsignal)[[3]][ ,1:D]
	
	# Add noise 
	E <- matrix( rnorm(I*J, 0, 1), I, J )
	correctFactor <- sqrt( (SStrue*noiseProp) / ( var(as.vector(E))*(1-noiseProp) ) ) 
	X <- Xsignal + E*correctFactor
	SSX <- var( as.vector( X ) )
	
	row.names(X) <- 1:I
	colnames(X) <- paste("var", 1:J, sep="")
	
	errorRatio <- 1 - (SStrue/SSX)
	
	# Check noise ratio 
	# ei <- eigen( t(X) %*% X )[[1]]
	# cei <- cumsum(ei) / sum(ei)
	
	return( list( 
				 X = X, 
				 Xsignal = Xsignal, 
				 TrueW = TrueW, 
				 errorRatio = errorRatio, 
				 scores = Z
				 # eigenvals = ei, 
				 # cumEigenvals = cei 
				 )  )
	
}	



#################################
# 			2nd Version
#################################

# Data Generating Function (version 2)
datageneration2 <- function( I, loadMat, variances ){

	J <- nrow(loadMat)
	D <- ncol(loadMat)

	# Variance-covariance matrix 
	Sigma <- loadMat %*% diag(variances) %*% t(loadMat)

	
	# Generate data matrix 
	X <- MASS::mvrnorm(n = I, mu = rep(0, J), Sigma = Sigma, empirical = TRUE ) 

	# Check noise ratio 
	ei <- eigen( t(X) %*% X )[[1]]
	cei <- cumsum(ei) / sum(ei)
	
	return( list( 
				 X = X, 
				 Sigma = Sigma,
				 eigenvals = ei, 
				 cumEigenvals = cei 
				 )  )	
}	
