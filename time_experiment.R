dir.create("res", showWarnings = FALSE)
N <- 5000 
Ps <- c(10*(1:10))
results <- array(dim=c(length(Ps), 5), dimnames = list(Ps, c("method","p","elapsed_time", "N", "avg_time")))

for (i in 1:length(Ps)){
    p <- Ps[i]
    results[i, "p"] <- p
    t1 <- Sys.time()
    invisible(gmat::chol_mh(N = N , p = p))
    t2 <- Sys.time()
    results[i,"elapsed_time"] <- as.numeric(difftime(t2,t1,unit = "sec"))
    results[i,"N"] <- N
    results[i,"avg_time"] <- results[i, "elapsed_time"]/N 
    results[i,"method"] <- 1
    saveRDS(difftime(t2,t1,unit = "sec"), file =paste0("res/t_chol_", p,".rds") )
}


write.csv(results, file = paste0("res/t_chol_",  format(Sys.time(), "%b%d%H-%M-%S%Y.csv")), quote = TRUE, row.names  = FALSE)


for (i in 1:length(Ps)){
    p <- Ps[i]
    results[i, "p"] <- p
    t1 <- Sys.time()
    invisible(sapply(1:N, FUN = function(n){return(randcorr::randcorr(p = p))},
	simplify = FALSE, USE.NAMES = FALSE))
    t2 <- Sys.time()
    results[i,"elapsed_time"] <- as.numeric(difftime(t2,t1,unit = "sec"))
    results[i,"N"] <- N
    results[i,"avg_time"] <- results[i, "elapsed_time"]/N 
    results[i,"method"] <- 2
    saveRDS(difftime(t2,t1,unit = "sec"), file =paste0("res/t_polar_", p,".rds") )
}


write.csv(results, file = paste0("res/t_polar_",  format(Sys.time(), "%b%d%H-%M-%S%Y.csv")), quote = TRUE, row.names  = FALSE)

for (i in 1:length(Ps)){
    p <- Ps[i]
    results[i, "p"] <- p
    t1 <- Sys.time()
	invisible(sapply(1:N, clusterGeneration::genPositiveDefMat, dim = p,
	covMethod = "onion", rangeVar = c(1,1), simplify = FALSE, USE.NAMES =
	FALSE))
    t2 <- Sys.time()
    results[i,"elapsed_time"] <- as.numeric(difftime(t2,t1,unit = "sec"))
    results[i,"N"] <- N
    results[i,"avg_time"] <- results[i, "elapsed_time"]/N 
    results[i,"method"] <- 3
    saveRDS(difftime(t2,t1,unit = "sec"), file =paste0("res/t_onion_", p,".rds") )
}


write.csv(results, file = paste0("res/t_onion_",  format(Sys.time(), "%b%d%H-%M-%S%Y.csv")), quote = TRUE, row.names  = FALSE)


for (i in 1:length(Ps)){
    p <- Ps[i]
    results[i, "p"] <- p
    t1 <- Sys.time()
	invisible(sapply(1:N, clusterGeneration::genPositiveDefMat, dim = p,
	covMethod = "c-vine", rangeVar = c(1,1), simplify = FALSE, USE.NAMES =
	FALSE))
    t2 <- Sys.time()
    results[i,"elapsed_time"] <- as.numeric(difftime(t2,t1,unit = "sec"))
    results[i,"N"] <- N
    results[i,"avg_time"] <- results[i, "elapsed_time"]/N 
    results[i,"method"] <- 3
    saveRDS(difftime(t2,t1,unit = "sec"), file =paste0("res/t_c-vine_", p,".rds") )
}


write.csv(results, file = paste0("res/t_c-vine_",  format(Sys.time(), "%b%d%H-%M-%S%Y.csv")), quote = TRUE, row.names  = FALSE)
