
#' @name clusters_km

#' @title clustering results of the k-mean algorithm
#' @description clusters data into two clusters. This functionis uses the \code{kmeans} function to cluster the data and exports the clustering results as well as the sum of square (SS) of clustering using the Euclidian distance.
#' @param  x matrix of data  (dim 1: samples (must be equal to dim 1 of X), dim 2: attributes (must be equal to dim 2 of X))
#' @param  k number of clusters ( this version considers 2 clusters )
#' @import stats
#' @import MASS
#' @import graphics
#' @return sum of square (SS) of clustring
#' @examples
#' {
#' X=rbind(matrix(rnorm(1000*2 ,4,.1),1000,2),matrix(rnorm(1000*2, 3, 0.2),1000,2))
#' M<- X[sample(nrow(X), 2),]
#' clusters_km(X,2)
#' }
#' @export
#'


clusters_km <- function(x, k=2) {
  mat=x
  centers <- mat[sample(nrow(mat), k),]
  theResult=kmeans(mat, centers, 4)
  mat1=cbind(mat,theResult$clusters[[4]])
  lst <- setNames(lapply(split(1:nrow(mat1), mat1[,3]), function(i) mat1[i,]), c("CL1", "CL2" ))
  CL1=lst$CL1[,-3]
  CL2=lst$CL2[,-3]
  CNT1=apply(CL1,2, mean)
  CNT2=apply(CL2,2, mean)
  ss1<- matrix(NA, nrow=nrow(CL1), ncol=1)
  for (i in 1:nrow(CL1)){
    ss1 [[i]]= (CL1[i,1]-CNT1[1])^2+(CL1[i,2]-CNT1[2])^2
  }
  ss2<- matrix(NA, nrow=nrow(CL2), ncol=1)
  for (i in 1:nrow(CL2)){
    ss2 [[i]]= (CL2[i,1]-CNT2[1])^2+(CL2[i,2]-CNT2[2])^2
  }
  ss=sum(ss1)+sum(ss2)
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  par(mfrow=c(1,2))
  for(i in 1:4) {
    plot(mat, col=theResult$clusters[[i]], main=paste("itteration:", i), xlab="x", ylab="y")
    points(theResult$centers[[i]], lwd = 8, pch=1, col=c(2,6))
  }
  list(clusters=theResult$clusters[[4]], centers=theResult$centers[[4]] )
}

