# Data Generating Function (version 1)
datageneration <- function( I, structMat, variances ){

    D <- ncol(structMat) 
    J <- nrow(structMat) 

    P <- matrix(rnorm(J * D), J, D)
    P[structMat == 0] <- 0
    P <- cbind(P, matrix(rnorm(J*(J-D)), J, J-D))

    result <- genLoading(P)
    if(result$converged == 1){
        P <- result$P
        Sigma <- P %*% diag(variances) %*% t(P)
        X <- MASS::mvrnorm(I, mu = rep(0, J), Sigma, empirical=F)
        # X <- MASS::mvrnorm(I, mu = rep(0, J), Sigma, empirical=T)
        return( list(X = X, P = P, Sigma = Sigma) )
    } else {
        stop("Data generation failed\n")
    }	
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
