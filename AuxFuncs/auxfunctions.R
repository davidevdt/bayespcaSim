# Function that allows centering and scaling of the data 
scaleX <- function(X, center = FALSE, scalingFactor = -1){
	Xnew <- X 
	if( center ){
		means <- apply(X, 2, mean)
		Xnew <- t( t(X) - means) 
		attr(Xnew, "center") <- means 
		rm(means)
	}
	
	if(scalingFactor >= 0){
	
		sds <- apply(Xnew, 2, function(x){
										 sqrt( sum(x^2, na.rm = TRUE) / (length(na.omit(x))- scalingFactor) ) 	
										 })
		Xnew <- t( t(Xnew) / sds )
		attr(Xnew, "scaling") <- sds
		rm(sds)
	}
	
	return( Xnew )

}





# Make variance of unimportant components 
varFunc <- function(varComps, J, error, scaleFirst = 2){

	D <- length(varComps) 
	otherVars <- exp( seq(log(1e-05), log(mean(varComps)/scaleFirst),
					  length.out = J - D ) )[(J-D):1]
	scaleFactor <- (error*sum(varComps)) / ((1-error)*sum(otherVars))

	return( c(varComps, otherVars*scaleFactor)  )
				

}


## Function that normalises the columns of a matrix 
Normalise <- function( M ){
	return(apply( M, 2, function(x)  x / sqrt(sum(x^2)) ))
}


## Function that calculates Tucker Congruence between two matrices 
TuckerCoef <- function(MatrixA, MatrixB){
  
  nrow_data <- dim(MatrixA)[1]
  ncol_data <- dim(MatrixA)[2]
  INDIC_Mat <- gtools::permutations(ncol_data, ncol_data)
  ncol_INDIC <- dim(INDIC_Mat)[1]
  TUCK <- array(NA, dim = c(ncol_INDIC, ncol_data))
  tucker_values <- array()
  tuckerr <- array()
  for(i in 1: ncol_INDIC) {
    MatrixB_perm <- MatrixB[, INDIC_Mat[i,]]
    teller <- 1

    for (r in 1: ncol_data){
      vec1 <- MatrixA[, r]
      vec2 <- MatrixB_perm[, r]
      cp <- t(vec1) %*% vec2
      var1 <- t(vec1) %*% vec1
      var2 <- t(vec2) %*% vec2

      if (var1 > 0 & var2 > 0){
        tuckerr[teller] <- psych::tr(cp)/sqrt(psych::tr(var1)*psych::tr(var2))
        teller <- teller + 1
      } else if (var2 == 0){
        tuckerr[teller] <- 0
        teller <- teller + 1
      }
    }

    tucker_values[i] <- mean(abs(tuckerr))
    TUCK[i,] <- tuckerr
  }

  k <- which(tucker_values == max(tucker_values))
  k <- k[1]

  perm <- INDIC_Mat[k,]
  tucker_value <- max(tucker_values)
  tucker_vector <- TUCK[k, ]

  return_tucker <- list()
  return_tucker$perm <- perm
  return_tucker$tucker_value <- tucker_value
  return_tucker$tucker_vector <- tucker_vector
  return(return_tucker)
}





# Function that takes the  data frame with sim. results and creates a new one for True/False positive rate 
# makeTNRateData <- function(res){
	
#	NRows <- nrow(res)
#	PosRate <- res %>% select(I,J,Sparsity,Noise,Method,PercOnes,ELBO)
#	NegRate <- res %>% select(I,J,Sparsity,Noise,Method,PercZeros,ELBO)
	
#	Weights <- c(rep("Non-zero Weights", NRows) , rep("Zero Weights", NRows) )
	
#	colnames(PosRate)[6] <- "PropCorrect"
#	colnames(NegRate)[6] <- "PropCorrect"
	
#	retData <- rbind(PosRate, NegRate)
#	retData <- data.frame(retData, "Weights" = Weights)
	
#	retData

#} 





# Function the returns average, rescaled ELBO's given a dataframe results 
groupELBOs <- function(res, crit = c("Method", "I", "J", "Sparsity", "Noise")){

	crit2 <- paste0("res$",crit)
	lst <- vector("list", length(crit))
	for( i in 1:length(lst) ){
		lst[[i]] <- eval(parse(text = crit2[i]))
	}

	ElboInd <- length(crit) + 1 
	
	
	ret <- aggregate(res$ELBO, by = lst, mean)
	rescaled <- (1 - ret[!is.na(ret[,ElboInd]),ElboInd] / min(ret[!is.na(ret[,ElboInd]),ElboInd]) ) 
	reRescaled <- rescaled + (1 - max(rescaled)) 
	
	ret[!is.na(ret[,ElboInd]), ElboInd] <-  reRescaled
	colnames(ret) <- c(crit, "ELBO" )
	ret 
} 




# Function the returns average, rescaled ELBO's given a dataframe results 
groupIndex <- function(res, crit = c("Method", "I", "J", "Sparsity", "Noise"), posNegRate = FALSE, ELBO = FALSE){

	if( posNegRate == TRUE ){
		crit <- c(crit, "Weights")
	}

	crit2 <- paste0("res$",crit)
	lst <- vector("list", length(crit))
	for( i in 1:length(lst) ){
		lst[[i]] <- eval(parse(text = crit2[i]))
	}

	
	Ind <- length(crit) + 1 
	
	if( ELBO ){
		ret <- aggregate(res$ELBO, by = lst, mean)

		rescaled <- (1 - ret[!is.na(ret[,Ind]),Ind] / min(ret[!is.na(ret[,Ind]),Ind]) ) 
		reRescaled <- rescaled + (1 - max(rescaled)) 
	
		ret[!is.na(ret[,Ind]), Ind] <-  reRescaled
		
		colnames(ret) <- c(crit, "ELBO" )
	}else{
		ret <- aggregate(res$RecErr, by = lst, mean)
		
		# rescaled <- 1 - (ret[,Ind] / max(ret[,Ind]))
		reRescaled <-  (ret[,Ind] / max(ret[,Ind])) 
		# reRescaled <- 1 - rescaled 
		# reRescaled <- ret[,Ind] / max( ret[ ,Ind]  )
		
		ret[,Ind] <- reRescaled
		colnames(ret) <- c(crit, "RecErr")
	}


	
	ret 
} 



### Adjust "aggregated" data 
adjustAggregated <- function(res, critMatrix = NULL, ELBO = FALSE, normalize = FALSE, a = 0.5, b = 1 ){

	if( ELBO ){
		IndFilt <- which( res$Method != "spca"  )
	}else{
		IndFilt <- 1:nrow(res)
	}
	
	
	IndMeth <- which(colnames(res) == "Method")
	Ind <- ncol(res)
	Ind0 <- (1:(Ind-1))[-IndMeth]
	
	
	ELBs <- res[IndFilt, Ind]

	if( is.null(critMatrix) ){
		aggr <- unique(res[,-c(IndMeth, Ind)])
		
		if(!normalize){
			for( j in 1:nrow(aggr) ){	
			
				tmp <- t(apply(res[,Ind0], 1, function(x) x == aggr[j,]) )
				selInd <- which(apply(tmp, 1, sum) == length(Ind0) ) 
				
				
				mx <- max(ELBs[IndFilt %in% selInd])
				ELBs[IndFilt %in% selInd] <- 1 + (ELBs[IndFilt %in% selInd] - mx) 
				
			}
		}else{
			for( j in 1:nrow(aggr) ){	
			
				tmp <- t(apply(res[,Ind0], 1, function(x) x == aggr[j,]) )
				selInd <- which(apply(tmp, 1, sum) == length(Ind0) ) 
				
				
				mx <- max(ELBs[IndFilt %in% selInd])
				mn <- min(ELBs[IndFilt %in% selInd])
				ELBs[IndFilt %in% selInd] <- a + ( ( (ELBs[IndFilt %in% selInd] - mn)*(b - a) ) / (mx - mn) )	
		
			}
		}		
	}else{
		aggr <- as.numeric(unique( critMatrix  ))
		
		if( !normalize ){
			for( j in 1:length(aggr) ){	
			
				tmp <- t(apply(critMatrix, 2, function(x) x == aggr[j]) )
				selInd <- which(tmp == TRUE ) 
				
				
				mx <- max(ELBs[IndFilt %in% selInd])
				ELBs[IndFilt %in% selInd] <- 1 + (ELBs[IndFilt %in% selInd] - mx) 
				
			}
		}else{
			for( j in 1:length(aggr) ){	
			
				tmp <- t(apply(critMatrix, 2, function(x) x == aggr[j]) )
				selInd <- which(tmp == TRUE ) 
				
				
				mx <- max(ELBs[IndFilt %in% selInd])
				mn <- min(ELBs[IndFilt %in% selInd])
				ELBs[IndFilt %in% selInd] <- a + ( ( (ELBs[IndFilt %in% selInd] - mn)*(b - a) ) / (mx - mn) )	
				
			}		
		}
	}
		
	res[IndFilt,Ind] <- ELBs
	res
	
}


# Function that calculates reconstruction error 
recErr <- function(Xorig, W, P){
	
	recErr <- sum( (Xorig - ( Xorig %*% W %*% t(P) ) )^2)	
	recErr

}




# Count whether there is the 0 value between two vectors 
zeroInVector <- function( vec ){

	ret <- 1 
	
	if( vec[1] <= 0 & vec[2] >= 0  ){
		ret <- 0
	}

	ret 
	
}





# Find zeros in HPD intervals 
zeroHPDI <- function( lst ){

	D <- length(lst)
	J <- nrow(lst[[1]])
	retMat <- matrix(-1, J, D)
	
	for( d in 1:D ){
	
		retMat[,d] <- apply(lst[[d]], 1, zeroInVector)
	
	}
	
	retMat

}




