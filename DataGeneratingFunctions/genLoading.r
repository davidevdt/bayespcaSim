# Function that generates an orthonormal loading matrix starting from a baseline zero/nonzero structure matrix. 
genLoading <- function( structMat, maxIt = 1000 ){
	
	it <- 0
	
	J <- nrow(structMat)
	D <- ncol(structMat)
	
	retMat <- matrix( rnorm(J*D), J, D )
	# retMat <- matrix( runif(J*D, -1, 1), J, D )
	retMat[structMat == 0] <- 0
	
	checkMat <- t(retMat) %*% retMat
	
	# while the matrix is not orthonormal, orthogonalize its columns 
	while( ( all.equal( checkMat, diag(D)) != TRUE ) && it < maxIt ){
	# while( ( all.equal(round(checkMat, 30), diag(D)) != TRUE ) && it < maxIt ){
	
		for( d2 in 2:D  ){
			for( d1 in 1:(d2-1) ){
			
				inter <- intersect( which(retMat[,d1] != 0), which(retMat[,d2] != 0)  )
				
				u <- retMat[inter, d1]
				v <- retMat[inter, d2]
				orthV <- v - as.numeric( (t(u) %*% v) / (t(u) %*% u) ) * u 

				retMat[inter, d2] <- orthV
			
			}
		}
			
		retMat <- Normalise( retMat )
		checkMat <- t(retMat) %*% retMat
		
		it <- it + 1 
				
	}
	
	if( it < maxIt ){
		return(list(P = retMat, converged=1))
	}else{
		return(list(P = retMat, converged=0))	
	}
	
}