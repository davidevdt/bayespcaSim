# Function that generates a matrix with random zero/nonzero elements. 
genZeroMat <- function(originalMat, sparsityprop){
		
	J <- nrow(originalMat)
	D <- ncol(originalMat)
	numZeros <- round(sparsityprop * J)
	
	for(d in 1:D) {
		originalMat[sample(1:J, numZeros), d] <- 0
	}
	
	return(originalMat)
	
}
