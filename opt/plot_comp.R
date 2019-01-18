library("ggplot2")
source("plot_utils.R")

f_sample <- c("polar" = rpolar, "chol" = rmh,  "onion" = ronion, "vine" = rvine)
method <- names(f_sample)

#### compute eigenvalues (100 matrices of dimension 50)
wd <- getwd()
dir.create(paste0(wd, "/plot"), showWarnings = FALSE)
eigen_list <- list()
N <- 100
p <- 50
for (m in method) {
	sample <- f_sample[[m]](N, p)
	
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
exp_volume <- N * 0.2^3 / vol_elliptope(n = 3)
rep <- 50

volume <- matrix(data = 0, nrow = length(method), ncol = 1, 
								 dimnames = list("method" = method, 
								 								"case" = c("-0.1 -- 0.1")))
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

### plotting elliptope with different methods
N <- 10000; p <- 3

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

