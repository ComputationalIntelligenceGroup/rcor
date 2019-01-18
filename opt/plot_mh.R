source('plot_utils.R')

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


