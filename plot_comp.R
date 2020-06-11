library("ggplot2")

ronion <- function(N, p) {
	sample <- array(dim = c(p, p, N))
	for (i in 1:N){
		sample[, , i] <- clusterGeneration::genPositiveDefMat(dim = p, rangeVar = c(1,1),
																													covMethod  ="onion")$Sigma
	}
	return(sample)
}
rvine <- function(N, p) {
	sample <- array(dim = c(p, p, N))
	for (i in 1:N){
		sample[, , i] <- clusterGeneration::genPositiveDefMat(dim = p, rangeVar = c(1,1),
																													covMethod  ="c-vine")$Sigma
	}
	return(sample)
}

rpolar <- function(N, p) {
	sample <- array(dim = c(p, p, N))
	for (i in 1:N){
		sample[, , i] <- randcorr::randcorr(p = p)
	}
	return(sample)
}

rmh <- function(N, p, ...) {
	return(gmat::chol_mh(N = N, p = p, ...))
}
f_sample <- c("polar" = rpolar, "chol" = rmh,  "onion" = ronion, "vine" = rvine)
method <- names(f_sample)

sub_sample <- function(x, bounds){
	
	p <- ncol(x)
  for (i in 1:p) {
    x <- x[x[, i] > bounds[1] & x[, i] < bounds[2], ]
  }
	return(x)
}

plot_elliptope <- function(x, bounds = NULL,
                         sort.fun = function(v){ return(sum(v^2))}, col = "red",
                         pallette =  colorRampPalette(c("red","blue")), file=NULL,... ){
  N <- nrow(x)
  if (!is.null(bounds)){
    x <- sub_sample(x = x, bounds  = bounds)
  }
  if (is.function(sort.fun)){
    idx <- sort(apply(x, MARGIN = 1, FUN = sort.fun), index.return  = TRUE)$ix
  }else{
    idx <- 1:N
  }
  if (is.function(pallette)){
    col = pallette(N)
  }else{
    col = col
  }
  rgl::plot3d(x, col = col,
              xlab = "x", ylab="y", zlab = "z",...)
  if (is.character(file)){
    rgl.postscript(file,fmt = "pdf")
  }
  
}

vol_elliptope <- function(n){
	if (n %% 2 == 0){
		num <- 2^((3 * n ^ 2 - 4 * n) / 4) * (gamma(n/2)) ^ n * prod(gamma(2 * (1 : ((n - 2) / 2))))
		den <- (gamma(n))^(n - 1)
		return(pi ^ ((n * (n - 2)) / 4) * num / den)
	}else{
		num <- prod(gamma(2 * (1 : ((n - 1) / 2))))
		den <- 2 ^ (((n - 1) ^ 2) / 4) * (gamma((n + 1) / 2 ))^(n - 1)
		return(pi ^ ((n ^ 2 - 1) / 4) * num/den)
	}
}

#### compute qube volume inside elliptope in 3 dimensions 
N <- 10000
p <- 3
f_sample <- c("polar" = rpolar, "chol" = rmh,  "onion" = ronion, "vine" = rvine)
method <- names(f_sample)
exp_volume <- N * 0.2^3 / vol_elliptope(n = 3)
rep <- 50

volume <- numeric(length(method))
names(volume) <- method
for (m in method) {
	for (r in 1:rep) {
		sample <- f_sample[[m]](N, p)
		sample <- gmat::vectorize(sample)

		region <- sub_sample(x = sample, bounds = c(-0.1, 0.1))
	
		volume[m] <- volume[m] + nrow(region)
	}
	volume[m] <- volume[m] / rep
}
vol <- list("exp" = exp_volume, "res" = volume)
saveRDS(vol, file = "vol.rds")

### plotting elliptope with different methods
N <- 10000; p <- 3

pdf(file="elliptope.pdf")
sample <- rmh(N = N, p = p, h = 1000, eps = 0.5)
sample <- gmat::vectorize(sample)
plot_elliptope(x = sample)
pairs(x = sample, lwd = 1 , pch  =20, cex  = 0.3, asp = 1, labels = c("x", "y", "z"))

sample <- rpolar(N = N, p = p)
sample <- gmat::vectorize(sample)
plot_elliptope(x = sample)
pairs(x = sample, lwd = 1 , pch  =20, cex  = 0.3, asp = 1, labels = c("x", "y", "z"))

sample <- ronion(N = N, p = p)
sample <- gmat::vectorize(sample)
plot_elliptope(x = sample)
pairs(x = sample, lwd = 1 , pch  =20, cex  = 0.3, asp = 1, labels = c("x", "y", "z"))

sample <- gmat::chol_iid(N = N, p = p)
sample <- gmat::vectorize(sample)
plot_elliptope(x = sample)
pairs(x = sample, lwd = 1 , pch  =20, cex  = 0.3, asp = 1, labels = c("x", "y", "z"))
dev.off()
