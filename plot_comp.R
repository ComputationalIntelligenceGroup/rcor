library("ggplot2")
source("plot_utils.R")
devtools::install_github("irenecrsn/gmat")

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

rmh <- function(N, p) {
	return(gmat::chol_mh(N = N, p = p, h = 1000, eps = 0.1))
}

f_sample <- c("polar" = gmat::chol_polar, "chol" = rmh,  "onion" = ronion, "vine" = rvine)
method <- names(f_sample)

#### compute eigenvalues (100 matrices of dimension 50)
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

