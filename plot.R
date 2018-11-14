source('plot_utils.R')
devtools::install_github("irenecrsn/gmat", ref = "rchol")

N <- 5000; h <- 1000; p <- 1000;

#### plotting acceptance ratio vs i/eps
eps <- c(0.0001, 0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1, 5, 10, 100)
i <- c(10, 25, 50, 100, 250, 400, 600, 750, 850, 999)
plot_acceptance_k(N = N, p = p, i = i, eps = eps, h = h)
plot_acceptance_eps(N = N, p = p, i = i, eps = eps, h = h)
plot_acceptance_3d(N = N, p = p, i = i, eps = eps, h = h)

## plotting time
p <- seq(from = 10, to = 100, by = 10)
method <- c("chol", "c-vine", "onion", "polar")
dir_name <- "res/"
plot_time(p = p, method = method, dir_name = "res", fname = "time.pdf", log = FALSE)
plot_time(p = p, method = method, dir_name = "res", fname = "time_log.pdf", log = TRUE)


