library("RColorBrewer")
library("reshape2")
library("ggplot2")

computeAcceptanceRatios <- function(x){
  N <- dim(x)[1]
  aR <- sapply(2:N, function(i){
    return(x[i] == x[i-1])
  })
  aRval <- 1 - length(aR[1:(N-1)][aR[1:(N-1)]])/(N - 1)
  return(aRval)
}

#### plotting acceptance ratio vs eps and i
plot_acceptance_3d <- function(N, p, i, eps, h) {

  ars <- array(dim = c(length(eps), length(i)), dimnames = list(eps = eps, i = i))
  for (k in 1:length(i)) {
    for (j in 1:length(eps)) {
    	sam <- sphereSample(N = N, p = p - i[k] + 1, i = i[k], eps = eps[j], h = h)
    	ars[j, k] <- computeAcceptanceRatios(sam)
    }
  }
  
  wd <- getwd()
  dir.create(paste0(wd, "/plot/"), showWarnings = FALSE)
  
  df <- melt(ars)

  pl <- ggplot(df, aes(x = eps, y = i, z = value)) +
  	geom_contour(aes(colour = ..level..)) +
  	scale_color_gradient(low = "blue", high = "red") +
  	theme(text = element_text(size = 20)) +
  	xlab("Perturbation variance") +
  	ylab("Row number") +
  	ggtitle("Contour plot of the acceptance ratio surface")
  
  ggsave(filename = "acceptance_3d.pdf", plot = pl, device = "pdf", path = "plot",
  			 width = 9, height = 4)
}

#### plotting acceptance ratio vs i
plot_acceptance_k <- function(N, p, i, eps, h) {
  eps_len <- length(eps)
  
  ars <- array(dim = c(length(i), length(eps)), dimnames = list(i = i, eps = eps))
  for (k in 1:length(i)) {
  	for (j in 1:eps_len) {
  		sam <- sphereSample(N = N, p = p - i[k] + 1, i = i[k], eps = eps[j], h = h)
  		ars[k, j] <- computeAcceptanceRatios(sam)
  	}
  }
  
  wd <- getwd() 
  dir.create(paste0(wd, "/plot/"), showWarnings = FALSE)
	
  palette <- colorRampPalette(colors = c("black","red"))
  colors <- palette(eps_len)
  
  df <- melt(ars)
  df$eps <- as.factor(df$eps)
  
  pl <- ggplot(df, aes(x = i, y = value, group = eps, color = eps)) +
  	geom_line() +
  	geom_point() +
  	theme(text = element_text(size = 20), legend.position = "bottom") +
  	scale_color_manual(values = colors) +
  	xlab("Row number") +
  	ylab("Acceptance ratio") +
  	ggtitle("")
  

  ggsave(filename = "acceptance_k.pdf", plot = pl, device = "pdf", path = "plot",
  			 width = 7, height = 5)
}

#### plotting acceptance ratio vs i
plot_acceptance_eps <- function(N, p, i, eps, h) {

	ars <- array(dim = c(length(i), length(eps)), dimnames = list(i = i, eps = eps))
	for (k in 1:length(i)) {
		for (j in 1:length(eps)) {
			sam <- sphereSample(N = N, p = p - i[k] + 1, i = i[k], eps = eps[j], h = h)
			ars[k, j] <- computeAcceptanceRatios(sam)
		}
	}
	
	wd <- getwd() 
	dir.create(paste0(wd, "/plot/"), showWarnings = FALSE)
	
	palette <- colorRampPalette(colors = c("black","red"))
	colors <- palette(length(i))
	
	df <- melt(ars)
	df$i <- as.factor(df$i)

	pl <- ggplot(df, aes(x = eps, y = value, group = i, color = i)) +
		geom_line() +
		geom_point() +
		theme(text = element_text(size = 20), legend.position = "bottom") +
		scale_color_manual(values = colors) +
		xlab("Perturbation variance") +
		ylab("Acceptance ratio") +
		ggtitle("")
	
	
	ggsave(filename = "acceptance_eps.pdf", plot = pl, device = "pdf", path = "plot",
				 width = 7, height = 5)
}


plot_time <- function(p, method, fname, dir_name = "res", log = FALSE) {
	
	res <- array(dim = c(length(p), length(method)), dimnames = list(p = p, method = method))
	
	for (i in 1:length(p)) {
		for (j in 1:length(method)) {
			time <- readRDS(file = paste0(dir_name,"/t_",method[j],"_",p[i], ".rds"))
			res[i, j] <- as.double(time, unit = "secs")
		}
	}
	
	wd <- getwd()
	dir.create(paste0(wd, "/plot/"), showWarnings = FALSE)
	
	palette <- colorRampPalette(colors = c("black","red"))
	colors <- palette(length(method))
	
	df <- melt(res)
	df$method <- as.factor(df$method)
	
	pl <- ggplot(df, aes(x = p, y = value, color = method, group = method)) +
		geom_line() +
		geom_point() +
		theme(text = element_text(size = 20), legend.position = "bottom") +
		xlab("Number of nodes") +
		ylab("Time in seconds") +
		scale_color_manual(values = colors)
	
	if (log == TRUE) {
		pl <- pl + scale_y_log10()
	}
		
	ggsave(filename = fname, plot = pl, device = "pdf", path = "plot/", width = 7,
				 height = 4)
}
