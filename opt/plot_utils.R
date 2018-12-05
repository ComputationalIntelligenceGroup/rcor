library("RColorBrewer")
library("reshape2")
library("ggplot2")

plot_acceptance_iter <- function(N, p, i, eps, h) {

	wd <- getwd()
	dir.create(paste0(wd, "/plot"), showWarnings = FALSE)

	sam <- gmat::mh_row(N = N, p = p - i + 1, i = i, eps = eps, h = h)
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
		sam <- gmat::mh_row(N = N, p = p - i + 1, i = i, eps = eps)
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
			sam <- gmat::mh_row(N = N, p = p - i[k] + 1, i = i[k], eps = eps[j], h = h)
		
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
		
		sam <- gmat::mh_row(N = N, p = p - i + 1, i = i, eps = eps[j], h = h)
			
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
		
		sam <- gmat::mh_row(N = N, p = p - i + 1, i = i, eps = eps[j], h = h)
		
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


plot_hist_dens_circle <- function(k = 10, N = 20000, eps = 0.1, h = 1000, title 
                                  = paste0("k = ", k, ", eps = ", eps)) {
  sample <- gmat::mh_row(N = N, p = 2, i = k, eps = eps, h = h, returnAll = F)
  circularDensityPlot(k = k, newPlot = T, lwd = 3, makeCircle = T, 
                      xaxt='n', ann=FALSE, yaxt = 'n', bty="n" , col = "black",lty  = 6 )
  title(title)
  circularHistogram(sample, makeCircle = F, newPlot = F, lwd = 1, breaks = "FD")
  ## plotting some references
  lines(c(0,0),c(-2,2), lwd=2)
  #text(x = -0.2, y = 0, "(0,0)")
  lines(c(0,1),c(0,0),lwd = 2)
}


#' plot density on circle
#' 
#' @example 
#'   circularDensityPlot(newPlot = T, makeCircle = T, k=10, norm=T)
#'   sapply(9:1, circularDensityPlot, norm=TRUE)

circularDensityPlot <- function(k = 1, fdens = function(th) return(cos(th)^k),  
                                norm = TRUE,  from = -pi/2, to = pi/2, by=0.03,
                                makeCircle=FALSE, newPlot=FALSE, col="red",...){
  gth <- seq(from = from , to = to, by = by)
  c <- cos(gth)
  s <- sin(gth)
  if (norm){
    if (k %% 2 == 0) { # even k
      l <- 1:(k/2)
      const <- pi*prod((2*l - 1)/(2*l))
    } else {
      if (k == 1) {
        const <- 2
      } else {
        l <- 1:((k - 1)/2)
        const <- 2*prod((2*l)/(2*l + 1))
      }
    }
  }else{
    const <- 1
  }
  f <- sapply(gth, fdens) / const
  if (newPlot){
    plot(c+f*c,s+f*s,col=col,type="l",asp=1,...)
  }else{
    lines(c+f*c,s+f*s,col=col,type="l",asp=1,...)
  }
  if (makeCircle){
    lines(c,s,type="l",asp=1, ...)
  }
}


#' plot a unit circle
plotCircle <- function(from = -pi/2, to = pi/2, by = 0.01, newPlot = TRUE, type = "l", ... ){
  gth <- seq(from,to,by)
  if (newPlot){
    plot(cos(gth),sin(gth), type = type, asp=1, ...)
  }else{
    lines(cos(gth),sin(gth),type = type, asp=1,...)
  }
  
}


#' plot a circular histogram
circularHistogram <- function(x, y = NULL, makeCircle = TRUE, newPlot = FALSE,
                              breaks = "Sturges" , ...){
  if (!is.null(dim(x))){
    if (length(dim(x)) > 1){
      if (dim(x)[2] == 1){
        x <- as.numeric(x)
      }else{
        x <- atan2(x[,2],x[,1])
      }
    }else{
      x <- as.numeric(x)
    }
  }else if (!is.null(y)){ # we interpret x,y as cartesian coordinates
    x <- atan2(y, x)
  } else{ #otherwise x is intepreted as angular values
    #do nada 
  }
  H <- hist(x , plot = F, breaks = breaks)
  cbreaks <- cos(H$breaks)
  sbreaks <- sin(H$breaks)
  gth <- seq(from = H$breaks[1] , to = H$breaks[length(H$breaks)], by = 0.01)
  c <- cos(gth)
  s <- sin(gth)
  if (newPlot){
    if (makeCircle){
      plot(c,s,type="l",asp=1, ...)
    }else{
      plot(c,s,type="l",asp=1, col="white", ...)
    }    
  }else{
    if (makeCircle){
      lines(c,s,type="l",asp=1, ...)
    }
  }
  for (j in 1:(length(H$breaks)-1)){
    lines(c(cbreaks[j], cbreaks[j]*(1+H$density[j])),
          c(sbreaks[j],sbreaks[j]*(1+H$density[j])) , ... )
    lines(c(cbreaks[j+1], cbreaks[j+1]*(1+H$density[j])),
          c(sbreaks[j+1],sbreaks[j+1]*(1+H$density[j])) , ... )
    lines(cbreaks[c(j,j+1)]*(1+H$density[j]),sbreaks[c(j,j+1)]*(1+H$density[j]), ...)
    
  }
  
}

##' sub sample 
##' 
##' @param x data row=observation, column = variables
##' @param bounds data.frame of bounding, same column as x and two row
##'               bounds[1,] contains the inferior limits
##'               bounds[2,] contains the upper limits
subSample <- function(x, bounds){
  for (i in 1:(dim(x)[2])){
    x <- x[x[,i]>= bounds[1,i],]
    x <-  x[x[,i]<= bounds[2,i],]
  }
  # x <- t(apply(x, MARGIN = 1, function(t){
  #   if (all(t>= bounds[1,]) && all(t<=bounds[2,]) ){
  #     return(t)
  #   }else{
  #     return(rep(NA,length(t)))
  #   }
  # }) )
  # x <- x[!is.na(x[,1]),]
  return(x)
}

### 3D plot ellittope
plotSample3D <- function(x, bounds = NULL,
                         sort.fun = function(v){ return(sum(v^2))}, col = "red",
                         pallette =  colorRampPalette(c("red","blue")), file=NULL,... ){
  N <- dim(x)[1]
  if (!is.null(bounds)){
    x <- subSample(x = x, bounds  = bounds)
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


volumeQube <- function(bounds){
	return(prod(apply(bounds, MARGIN = 2, function(b) return(abs(b[1]-b[2])))))
} 
volumeElliptope <- function(n){
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

plotSample3D <- function(x, bounds = NULL,
                         sort.fun = function(v){ return(sum(v^2))}, col = "red",
                         pallette =  colorRampPalette(c("red","blue")), file=NULL,... ){
  N <- dim(x)[1]
  if (!is.null(bounds)){
    x <- subSample(x = x, bounds  = bounds)
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


