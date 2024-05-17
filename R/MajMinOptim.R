#' @name MajMinOptim

#' @title majorization-minimization optimization
#' @description Finding the optimized majorization-minimization centers
#' @param  X matrix of data  (dim 1: samples (must be equal to dim 1 of X), dim 2: attributes (must be equal to dim 2 of X))
#' @param  Z  is a n by k matrix where for all i and j, zi,j is abinary variable that is equal to 1 if the case i is assigned to cluster j and zero otherwise. (dim 1: samples (must be equal to dim 1 of X), dim 2: attributes (must be equal to dim 2 of X))
#' @param M  initial seleted centroids (randomly or another method)
#' @param eps a threshold value assumed as 0.0001
#' @param lambda a threshold value assumed as 0.5
#' @import stats
#' @import graphics
#' @import MASS
#' @return The optimized majorization-minimization centers.
#' @examples
#' {
#' X=rbind(matrix(rnorm(1000*2 ,4,.1),1000,2),matrix(rnorm(1000*2, 3, 0.2),1000,2))
#' M <- X[sample(nrow(X), 2),]
#' distsToCenters <- Euclid(X, M)
#' clusters <- apply(distsToCenters, 1, which.min)
#' Z <- matrix(0, nrow = NROW(X), ncol = 1)
#' for(i in 1:NROW(X))
#' if (clusters[[i]] == 1)
#'     Z[i,]=clusters[[i]]
#' Z=cbind(Z, 1-Z)
#' MajMinOptim(X,Z,M ,eps=1e-4, lambda=.5)
#' }
#' @export
MajMinOptim <- function(X,Z, M, eps, lambda) {
  D <- as.matrix(dist(Z%*%M))
  dim(Z)
  D <- pmax(D, eps)
  w_matrix <- matrix(0, nrow = NROW(X), ncol = NROW(X))
  w_matrix[cbind(seq(1, NROW(X) - 1),
                 seq(2, NROW(X)))] <- 1
  w_matrix[cbind(seq(2, NROW(X)),
                 seq(1, NROW(X) - 1))] <- 1
  W=w_matrix
  V <- -W/D
  diag(V) <- -rowSums(V)
  ZX <- crossprod(Z, X)
  ZVZ <- crossprod(Z, crossprod(V, Z))
  n <- nrow(X)
  G <- crossprod(Z,Z%*%M)-ZX+n*lambda*crossprod(ZVZ,M)
  Mopt <- solve(crossprod(Z)+n*lambda*ZVZ, ZX)
 return(Mopt)
}

