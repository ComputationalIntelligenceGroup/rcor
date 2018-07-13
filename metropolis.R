require(compiler)
enableJIT(3)


# sample correlations matrices via polar parametrization of the cholesky factor
rPolarUnif <- function(N = 1,p = 10, ...){
 sam <- array(dim = c(p,p,N), data = 0)
 for (i in 1:N){
     temp <- rcoef_polar(p, ...)
     sam[,,i] <- temp %*% t(temp) 
 }
 return(sam)
}




#' Sample coefficient matrix with polar parametrization
#'
#' Sample a polar parametrization of the Cholesky factor in the LL' decomposition of
#' the concentration matrix.
#'
#' @rdname polar
#'
#' @param p positive integer, the dimension of the square matrix
#' @param method String, method to sample from `sin(x)^k`: `numeric` or
#' `recursive`
#'
#' @details The mcoef factor is sampled such that mcoefmcoef' is uniformly
#'distributed among the correlation matrices.
#'
#' @return `mcoef` a lower triangular matrix such that mcoefmcoef' has 1 on diagonal.
#
rcoef_polar <- function(p = 100,
                         method = 'numeric') 
{	

	mcoef <- matrix(nrow = p,ncol = p,data = 0)
        diag(mcoef) <- 1
	theta <- matrix( data = pi/2, nrow = p, ncol = p)

	diag(mcoef) <- 1

	for (icol in 1:(p - 1)) {
		for (irow in (icol + 1):p) {
				theta[irow, icol] <- .rsin(n = 1, k = p - icol, method = method)
				mcoef[irow, icol] <- cos(theta[irow, icol])
		}
		if (icol >= 2) {
			for (j in 1:(icol - 1)) {
				mcoef[icol:p, icol] <- mcoef[icol:p, icol] * sin(theta[icol:p, j])
			}
		}	
	}
	mcoef[p, p] <- prod(sin(theta[p, 1:(p - 1)]))
	return(mcoef)
}



#' @rdname polar
#'
#' @param k exponent for `sin^k`
#' @param x value where `sin^k` will be calculated
.sin_k <- function(x, k) {
	return (sin(x)^k)
}

#' @rdname polar
#'
#' @importFrom stats integrate
.sin_k_cum <- function(x, k, method = "numeric") {
	
	if (x <= 0) {
		return (0)
	} 
	
	if (x >= pi) {
		return (1)
	}

	if (method == "numeric") {
		const <- integrate(.sin_k, lower = 0, upper = pi, k = k)$value
		return(integrate(.sin_k, lower = 0, upper = x, k = k)$value / const)
	}
	const <- .sin_int(pi, k)
	return(.sin_int(x, k) / const)
}

#' @rdname polar
#' 
#' @details recursive computation 
.sin_int <- function(x, k = 1){
	if (length(x)>1){
		return(sapply(x,.sin_int,k))
	}
	if (x <= 0) {
		x<-0
	}
	if (x > pi) {
		x<-pi
	}
	if (k < 0) return(0)
	if (k == 0) {
		return(x)
	} else if (k == 1) {
		return(1-cos(x))
	} else {
		return(-(1/k)*cos(x)*.sin_k(x, k-1) + ((k-1)/k)*.sin_int(x, k-2)  )
	}
}

#' @rdname polar
#'
#' @param n Number of samples to generate
#'
#' @importFrom stats uniroot
.rsin <- function(n, k = 1, method = 'numeric') {
	.sin_k_inv_unif <- function(x, u) {
		return (.sin_k_cum(x, k, method) - u)
	}

	.sin_k_invsampl <- function(u) {
		return(uniroot(.sin_k_inv_unif, u = u, interval = c(0, pi), 
					   extendInt = "upX")$root)
	}

	return(sapply(runif(n), .sin_k_invsampl))
}


#' sampling on sphere proportionally to a power of the first coordinate
#' 
#' @param N sample size
#' @param p dimension
#' @param i exponent of the density
#' @param h heating phase size
#' 
#'  Metropolis-Hasting algorithm to sample in the n dimensional semi-sphere (x_1>0)  
#' @export
sphereSample <- function(N = 1, p, i = 1, h = 100, eps = 0.01,returnAll = FALSE){
    Tot <- h + N #total number of iteration of MH
    Sample <- matrix(nrow = Tot, ncol = p) #obj initialization 
    Sample[1,] <- rnorm(n = p, mean = 0, sd = 1) #first point
    Sample[1,1]<-abs(Sample[1,1]) #absolute value first component (has to be positive)
    Sample[1,] <- Sample[1,]/sqrt(sum(Sample[1, ] ^ 2)) #normalization 
    for (j in 2:Tot){
        prop <- Sample[j-1,] + rnorm(n = p, mean = 0, sd = eps) # perturbate previous sample
        prop <- prop / sqrt(sum(prop ^ 2)) #normalize propozed
        if ((prop[1]>0) && (log(runif(1)) <= i*log((prop[1])) - i*log(Sample[j - 1 , 1]))){
            Sample[j,] <- prop
        }else{
            Sample[j,] <- Sample[j - 1 , ]
        }
    }
    if (!returnAll){
        Sample <- Sample[(h + 1): Tot, ]
    }

    return(Sample)
}


#'
#' @export
sampleU <- function(N = 1,p=10, h = 100, eps = 0.1){
    U <- array(dim=c(p,p,N),data = 0)
    U[p,p,1:N] <- 1
    for (i in 1:(p-1)){
        su <- sphereSample(N = N,p = p - i + 1 , i = i ,h = h, eps = eps)
        U[i,i:p,1:N] <- t(su)
    }
    return(U)
}


#' rCholUnif
#' 
#' @param  N 
#' @param return.minvector logical, if TRUE the minimimal vector representation is returned (useful to plot in the elliptope)
#' @param ... additional parameters
#' 
#' @export
rCholUnif <- function(N=1, p = 10, return.minvector = FALSE, ... ){
    sU <- sampleU(N = N, p = p, ... )
    vsC <- apply(sU,MARGIN = 3, function(U) return(U%*%t(U)) )
    sC <- array(data=vsC, dim=dim(sU))
    if (return.minvector){
        mv <- apply(sC, MARGIN = 3, function(m){
                        return(m[upper.tri(m)])
})
        return(t(mv))
    }else{
        return(sC)     
    }
}


