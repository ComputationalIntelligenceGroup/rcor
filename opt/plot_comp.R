library("ggplot2")
source("plot_utils.R")

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

rmh <- function(N, p) {
	return(gmat::chol_mh(N = N, p = p, h = 1000, eps = 0.1))
}

f_sample <- c("polar" = rpolar, "chol" = rmh,  "onion" = ronion, "vine" = rvine)
method <- names(f_sample)

#### compute eigenvalues (100 matrices of dimension 50)
wd <- getwd()
dir.create(paste0(wd, "/plot"), showWarnings = FALSE)
eigen_list <- list()
for (m in method) {
	sample <- readRDS(paste0("res_comp/", m, ".rds"))
	
	eigen_list[[m]] <- apply(sample, MARGIN = 3, function(M) return(eigen(M)$values))
}
df <- dplyr::bind_rows(eigen_list)
df <- reshape2::melt(df, measure.vars = 1:ncol(df), variable.name = "method")
pl <- ggplot(df, aes(value, group = method, color = method)) + 
	geom_density(size = 1) + 
	theme(text = element_text(size = 20), legend.position = "bottom") +
	ylab("Estimated density") +
	xlab("Eigenvalue")
ggsave(filename = "densities_eigen.pdf", plot = pl, device = "pdf", path = "plot",
			 width = 7, height = 4)

#### compute qube volume inside elliptope in 3 dimensions 
N <- 10000
p <- 3
f_sample <- c("polar" = rpolar, "chol" = rmh,  "onion" = ronion, "vine" = rvine)
method <- names(f_sample)
volume <- matrix(data = 0, nrow = length(method), ncol = 1, 
									dimnames = list("method" = method, 
																	"case" = c("-0.1 -- 0.1")))
exp_volume <- N * 0.2^3 / vol_elliptope(n = 3)
rep <- 50
for (m in method) {
	for (r in 1:rep) {
		sample <- f_sample[[m]](N, p)
		sample <- gmat::vectorize(sample)

		region <- sub_sample(x = sample, bounds = c(-0.1, 0.1))
	
		volume[m, 1] <- volume[m, 1] + nrow(region)
	}
	volume[m, ] <- volume[m, ] / rep
}
vol <- list("exp" = exp_volume, "res" = volume)
saveRDS(vol, file = "vol.rds")

sample <- rmh(N = 10000, p = 3)
sample <- gmat::vectorize(sample)
plot_elliptope(x = sample)
