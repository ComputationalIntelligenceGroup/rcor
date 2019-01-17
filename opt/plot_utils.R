library("RColorBrewer")
library("reshape2")
library("ggplot2")

plot_acceptance_iter <- function(N, p, i, eps, h) {

	wd <- getwd()
	dir.create(paste0(wd, "/plot"), showWarnings = FALSE)

	sam <- gmat::mh_sphere(N = N, k = p - i + 1, i = i, eps = eps, h = h)
	ars <- sapply(X = 2:N, FUN = function(n) {
		return(computeAcceptanceRatios(sam[1:n, ]))
	})
		
	df <- data.frame(x = 2:N, y = ars)

	pl <- ggplot(df, aes(x = x, y = y)) +
		geom_line() +
		theme(text = element_text(size = 20), legend.position = "bottom") +
		xlab("Iteration") +
		ylab("Acceptance ratio") +
		ylim(0, 1) + 
		ggtitle(paste0("Row number = ", i, ". Variance = ", eps))
		
	ggsave(filename = "acceptance_iter.pdf", plot = pl, device = "pdf", path = "plot",
				 width = 7, height = 5)
}

plot_averages <- function(N, p, i, eps, file_name, n_paths = 5, plot_title) {
	
	avgs <- matrix(nrow = n_paths, ncol = N)
	for (path in 1:n_paths) {
		sam <- gmat::mh_sphere(N = N, k = p - i + 1, i = i, eps = eps)
		avgs[path, ] <- sapply(1:N, function(n) {return (mean(sam[1:n, 2]))})
	}

	wd <- getwd() 
	dir.create(paste0(wd, "/plot/"), showWarnings = FALSE)
	
	palette <- colorRampPalette(colors = brewer.pal(n = min(9, n_paths), 
																									name = "Dark2"))
	colors <- palette(n_paths)
	
	df <- melt(avgs, varnames = c("Repetition", "Iteration"))
	df$Repetition <- as.factor(df$Repetition)
	
	pl <- ggplot(df, aes(x = Iteration, y = value, group = Repetition, color = Repetition)) +
		geom_line() +
		theme(text = element_text(size = 20), legend.position = "bottom") +
		scale_color_manual(values = colors) +
		xlab("Iteration") +
		ylab("Empirical average") +
		ylim(-0.5, 0.5) +
		ggtitle(plot_title)
	
	ggsave(filename = file_name, plot = pl, device = "pdf", path = "plot",
				 width = 7, height = 5)
}

plot_iterations <- function(N, p, i, eps, h) {
	
	wd <- getwd()
	dir.create(paste0(wd, "/plot"), showWarnings = FALSE)
	
	pdf(file = "plot/iterations.pdf", width = 7, height = 5)
	
	for (j in 1:length(eps)) {
		
		for (k in 1:length(i)) {
			sam <- gmat::mh_sphere(N = N, k = p - i[k] + 1, i = i[k], eps = eps[j], h = h)
		
			df <- data.frame(x = 1:N, y = sam[, 2])

			pl <- ggplot(df, aes(x = x, y = y)) +
				geom_line() +
				theme(text = element_text(size = 20)) +
				xlab("Iteration") +
				ylab("Sample value") +
				ggtitle(paste0("Row number = ", i[k], ". Variance = ", eps[j]))
		
			print(pl)
		}
	}
	
	dev.off()
}

plot_pointwise_3D_second_row <- function(N, eps, h) {
	
	p <- 3
	i <- 2
	
	wd <- getwd()
	dir.create(paste0(wd, "/plot"), showWarnings = FALSE)
	
	pdf(file = "plot/pointwise_3D_second_row.pdf", width = 7, height = 5)
	
	for (j in 1:length(eps)) {
		
		sam <- gmat::mh_sphere(N = N, k = p - i + 1, i = i, eps = eps[j], h = h)
			
		df <- data.frame(x = sam[, 1], y = sam[, 2], Iteration = 1:N)
			
		pl <- ggplot(df, aes(x = x, y = y)) +
			geom_point(aes(colour = Iteration)) +
			theme(text = element_text(size = 20)) +
			xlab("First coordinate") +
			ylab("Second coordinate") +
			xlim(0, 1) + 
			ylim(-1, 1) +
			ggtitle(paste0("Variance = ", eps[j]))
			
		print(pl)
	}
	
	dev.off()
}

plot_pointwise_3D_first_row <- function(N, eps, h) {
	
	p <- 3
	i <- 3
	
	wd <- getwd()
	dir.create(paste0(wd, "/plot"), showWarnings = FALSE)
	
	pdf(file = "plot/pointwise_3D_first_row.pdf", width = 7, height = 5)
	
	for (j in 1:length(eps)) {
		
		sam <- gmat::mh_sphere(N = N, k = p - i + 1, i = i, eps = eps[j], h = h)
		
		df <- data.frame(x = sam[, 1], y = sam[, 2], Iteration = 1:N)
		
		pl <- ggplot(df, aes(x = x, y = y)) +
			geom_point(aes(colour = Iteration)) +
			theme(text = element_text(size = 20)) +
			xlab("First coordinate") +
			ylab("Second coordinate") +
			xlim(0, 1) + 
			ylim(-1, 1) +
			ggtitle(paste0("Variance = ", eps[j]))
		
		print(pl)
	}
	
	dev.off()
}

##' sub sample 
##' 
##' @param x data row=observation, column = variables
##' @param bounds vector of bounding,
##'               bounds[1] contains the inferior limits
##'               bounds[2] contains the upper limits
sub_sample <- function(x, bounds){
  for (i in 1:ncol(x)) {
    x <- x[x[, i] > bounds[1] & x[, i] < bounds[2], ]
  }
  return(x)
}

plot_elliptope <- function(x, bounds = NULL,
                         sort.fun = function(v){ return(sum(v^2))}, col = "red",
                         pallette =  colorRampPalette(c("red","blue")), file=NULL,... ){
  N <- dim(x)[1]
  if (!is.null(bounds)){
    x <- sub_sample(x = x, bounds  = bounds)
  }
  if (is.function(sort.fun)){
    idx <- sort(apply(x, MARGIN = 1,FUN = sort.fun), index.return  = TRUE)$ix
  }else{
    idx <- 1:N
  }
  if (is.function(pallette)){
    col = pallette(N)
  }else{
    col = col
  }
  rgl::plot3d(x[,], col = col,
              xlab = "v1", ylab="v2", zlab = "v3",...)
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

