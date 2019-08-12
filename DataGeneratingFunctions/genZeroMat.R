# Function that generates a matrix with random zero/nonzero elements. 
genZeroMat <- function(J, D, sparsityprop){

	numZeros <- ceiling(sparsityprop * J * D)
	retMat <- matrix(1, J, D)
	
	zeroElem <- sample(1:(J*D), numZeros, rep = FALSE )
	retMat[zeroElem] <- 0 
	
	return( retMat )
	
	
}