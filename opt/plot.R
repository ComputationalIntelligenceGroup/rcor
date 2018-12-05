source('plot_utils.R')
library("rgl")
library("RColorBrewer")

N <- 5000; h <- 1000; p <- 1000;

eps <- c(0.0001, 0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1)
i <- c(10, 25, 50, 100, 250, 400, 600, 750, 850, 999)

#### plotting acceptance ratio vs iterations
plot_acceptance_iter(N = N, p = 3, i = 2, eps = 0.1, h = h)

#### plotting behaviour of empirical average coordinate-wise
## seems coordinate-independent, justification?
plot_averages(N = N, p = p, i = i[1], eps = eps[1],	n_paths = 6, file_name = "averages_down.pdf",
							plot_title = paste0("Row number = ", i[1], ". Variance = ", eps[1]))
plot_averages(N = N, p = p, i = i[round(length(i)/2)], eps = eps[round(length(eps)/2)],	n_paths = 6, file_name = "averages_middle.pdf",
							plot_title = paste0("Row number = ", i[round(length(i)/2)], ". Variance = ", eps[round(length(eps)/2)]))
plot_averages(N = N, p = p, i = i[length(i)], eps = eps[length(eps)],	n_paths = 6, file_name = "averages_up.pdf",
							plot_title = paste0("Row number = ", i[length(i)], ". Variance = ", eps[length(eps)]))


#### plotting behaviour of sample coordinate-wise
plot_iterations(N = N, p = p, i = i, eps = eps, h = h)

#### plotting behaviour of sample pointwise
plot_pointwise_3D_second_row(N = N, eps = eps, h = h)

### uniform sampling with Cholesky decomposition
N <- 10000
s_unif <- gmat::chol_mh(N = N, p = 3, h = 1000, eps = 0.5)
s_unif <- gmat::vectorize(s_unif)
plot3d(x = s_unif[,1], y = s_unif[,2], z = s_unif[,3], xlab = "x", ylab = "y", zlab = "z") 
pairs(x = s_unif, lwd = 1 , pch  =20, cex  = 0.3, asp = 1, labels = c("x", "y", "z"))

### sample iid coefficients 
N <- 10000
s_iidc <- gmat::chol_iid(N = N, p = 3)
s_iidc <- gmat::vectorize(s_iidc)
plot3d(x = s_iidc[1,], y = s_iidc[2,], z = s_iidc[3,], xlab = "x", ylab = "y", zlab = "z") 
pairs(x = t(s_iidc), lwd = 1 , pch  =20, cex  = 0.3, asp = 1, labels = c("x", "y", "z"))

### different densities k=1:kmax
kmax <- 20
col <- brewer.pal(9, "OrRd")
pall2 <- colorRampPalette(col[-c(1,2)])
col <- pall2(kmax) #useless if kmax <= 12

wd <- getwd()
dir.create(paste0(wd, "/plot"), showWarnings = FALSE)
pdf(file = "plot/rwm_densities.pdf")   
circularDensityPlot(k=kmax, col = col[kmax], newPlot =  TRUE, makeCircle = TRUE,
                    lwd = 2,xaxt='n', ann=FALSE, yaxt = 'n', bty="n",ylim=c(-1.5,1.5) )
for (k in  seq(from=kmax-2, to = 1, by = -2)){
  circularDensityPlot(k = k, col = col[k], lwd=2, makeCircle = F  )
}
lines(c(0,0),c(-2,2), lwd=2)
lines(c(0,1),c(0,0),lwd = 2)
legend(-1.2, 1, legend=paste("k =",c(1,kmax)),
       col=col[c(1,kmax)], lty=1, cex=0.8)
dev.off()

#### plotting real density and histogram
wd <- getwd()
dir.create(paste0(wd, "/plot"), showWarnings = FALSE)
pdf(file = paste("plot/hist_dens_circle.pdf"))
k <- c(1, 3, 5, 10, 20, 30)
N <- 20000
eps <- c(0.01, 0.1, 1, 10, 100, 1000)
for (k1 in k) {
  for (eps1 in eps) {
    plot_hist_dens_circle(k = k1, N = N, eps = eps1)  
  }
}
dev.off()

## 3D plots
##generate pallette
col <- brewer.pal(9, "Greys")
pallette <- colorRampPalette(col)
N <- 10000
dag <- gmat::rgraph(p = 4, d = 0.6, dag = TRUE)
plot(dag)
sam <- gmat::chol_mh(N = N, dag = dag)
sam <- gmat::vectorize(sam)
plotSample3D(sam[,c(2,3,5)],pallette = pallette, sort.fun = function(t) return(t[1]))

plot(sam[,c(2,3)],asp=1)
hist(sam[,1],freq = F)

sort(sam[,1],index.return =T) ->sorted
index <- 1:N
index <- sorted$ix
plot3d(sam[index,], col = pallette(N),
        xlab = "v1", ylab="v2", zlab = "v3")

########## unif iid
sam_iid <- gmat::rgbn_iid(N = N, p = 3, dag = dag, return.minvector = TRUE)
plot3d(x = sam_iid[,1], y = sam_iid[,2], z = sam_iid[,3]) 
plotSample3D(sam_iid,pallette = pallette, sort.fun = function(t) return(t[1]))

col <- brewer.pal(9, "Greens")
pallette <- colorRampPalette(col)
sam <- gmat::mh_row(N = N, p = 3, i = 8, eps = 0.1, h = 10000)
plotSample3D(sam[,]*1.005, pallette = pallette, sort.fun = function(t) return(t[1]),
             type="p",aspect=F, size=8, add=T, point_antialias = T)

### non saturated model
N <- 10000
dag2 <- gmat::rgraph(p = 4, d = 0.6, dag = TRUE)
plot(dag2)

sam <- gmat::chol_mh(N = N, dag = dag2)
sam <- gmat::vectorize(sam)
sam[1,]
plot(x = sam[,1], y = sam[,2], asp=1)
pairs(x = sam[,], lwd = 1 , pch  =20, cex  = 0.3)

sam_iid <- gmat::chol_iid(N = N, p = 3, dag = dag2)
sam_iid <- gmat::vectorize()
dim(sam_iid)
plot(x = sam_iid[,1], y = sam_iid[,2], asp=1) 


