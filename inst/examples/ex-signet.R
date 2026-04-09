library(SIGNET)

## random seed for reproducibility
set.seed(314159)

## adjacency matrix
p <- 50
mydat <- matrix(rnorm(1000 * p), ncol = p)
myadj_mat <- cor(mydat)

## implementation of signet
myres <- signet(adj_mat = myadj_mat, Kmax = 10, method = "sil-gap")
