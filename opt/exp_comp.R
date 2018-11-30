devtools::install_github("irenecrsn/gmat")

N <- 100
p <- 50

sample_polar <- gmat::chol_polar(N = N, p = p)
sample_mh <- gmat::chol_mh(N = N, p = p, h=1000, eps = 0.1)
sample_onion <- array(dim=dim(sample_mh))
for (i in 1:N){
 sample_onion[,,i] <- clusterGeneration::genPositiveDefMat(dim = p, rangeVar = c(1,1), 
                                    covMethod  ="onion")$Sigma
}
sample_vine <- array(dim = dim(sample_mh))
for (i in 1:N){
 sample_vine[,,i] <- clusterGeneration::genPositiveDefMat(dim = p, rangeVar = c(1,1),
                                   covMethod  ="c-vine")$Sigma
}

wd <- getwd()
dir_name <- "res_comp"
dir.create(paste0(wd, "/",dir_name), showWarnings = FALSE)
saveRDS(sample_polar, file = paste0(dir_name,"/polar.rds"))
saveRDS(sample_mh, file = paste0(dir_name,"/chol.rds"))
saveRDS(sample_onion, file = paste0(dir_name,"/onion.rds"))
saveRDS(sample_vine, file = paste0(dir_name,"/vine.rds"))

