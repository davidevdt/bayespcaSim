# Function that generates an orthonormal loading matrix starting from a baseline zero/nonzero structure matrix. 
genLoading <- function( P, maxIt = 1000 ){
	it <- 0
	
	J <- nrow(P)
	D <- ncol(P)
	
	checkMat <- t(P) %*% P
	
	# while the matrix is not orthonormal, orthogonalize its columns 
	while( ( all.equal( checkMat, diag(J)) != TRUE ) && it < maxIt ){
	
		for( d2 in 2:D  ){
			for( d1 in 1:(d2-1) ){			
				inter <- intersect( which(P[,d1] != 0), which(P[,d2] != 0)  )
				
				u <- P[inter, d1]
				v <- P[inter, d2]
				orthV <- v - as.numeric( (t(u) %*% v) / (t(u) %*% u) ) * u 

				P[inter, d2] <- orthV			
			}
		}
			
		P <- Normalise( P )
		checkMat <- t(P) %*% P
		
		it <- it + 1 
				
	}
	
	if( it < maxIt ){
		return(list(P = P, converged=1))
	}else{
		return(list(P = P, converged=0))	
	}	
}