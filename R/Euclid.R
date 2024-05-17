#' @name Euclid

#' @title Euclidian distance
#' @description Calculates the Euclidian distance between points. This function can use in \code{kmeans} function to do the clustering procedure using the Euclidian distance.
#' @param  x matrix of data  (dim 1: samples (must be equal to dim 1 of X), dim 2: attributes (must be equal to dim 2 of X))
#' @param mu  initial seleted centroids (randomly or another method).
#' @import stats
#' @import graphics
#' @import MASS
#' @return Euclidian distance between two points.
#' @examples
#' {
#' X=rbind(matrix(rnorm(1000*2 ,4,.1),1000,2),matrix(rnorm(1000*2, 3, 0.2),1000,2))
#' M <- X[sample(nrow(X), 2),]
#' Euclid(X,M)
#' }
#' @export
#'
Euclid <- function(x, mu) {
  distanceMatrix <- matrix(NA, nrow=dim(x)[1], ncol=dim(mu)[1])
  for(i in 1:nrow(mu)) {
    distanceMatrix[,i] <- sqrt(rowSums(t(t(x)-mu[i,])^2))
  }
  distanceMatrix
}
