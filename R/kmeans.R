#' @name kmeans

#' @title k-means function
#' @description k-means algorithm in clustering. This function export the clustered results based on one replication of the k-means method
#' @param  x matrix of data  (dim 1: samples (must be equal to dim 1 of X), dim 2: attributes (must be equal to dim 2 of X))
#' @param centers  initial seleted centroids (randomly or another method)
#' @param nItter Number of itteration function
#' @import stats
#' @import MASS
#' @import graphics
#' @return clustered results based on k-means methods.
#' @examples
#' {
#' X=rbind(matrix(rnorm(1000*2 ,4,.1),1000,2),matrix(rnorm(1000*2, 3, 0.2),1000,2))
#' M <- X[sample(nrow(X), 2),]
#' kmeans(X,M, 4)
#' }
#' @export
kmeans <- function(x, centers, nItter=4) {
  clusterHistory <- vector(nItter, mode="list")
  centerHistory <- vector(nItter, mode="list")
  for(i in 1:nItter) {
    distsToCenters <- Euclid(x, centers)
    clusters <- apply(distsToCenters, 1, which.min)
    centers <- apply(x, 2, tapply, clusters, mean)
    dis2 <- apply(x, 2, tapply, clusters, dist)
    dis2=as.matrix(dis2)
    clusterHistory[[i]] <- clusters
    centerHistory[[i]] <- centers
  }
  list(clusters=clusterHistory, centers=centerHistory )
}

