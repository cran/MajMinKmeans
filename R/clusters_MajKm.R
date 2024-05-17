#' @name clusters_MajKm

#' @title clustering results of the majorized k-mean algorithm
#' @description clusters data into two clusters with a majorization k-means This functionis use a hybrid of the k-means and the majorizaion-minimazation method to cluster the data and exports the clustering results as well as the sum of square (SS) of clustering
#' @param  X matrix of data  (dim 1: samples (must be equal to dim 1 of X), dim 2: attributes (must be equal to dim 2 of X))
#' @param  k number of clusters ( this version considers 2 clusters )
#' @param  La the tunnung parameter
#' @import stats
#' @import MASS
#' @import graphics
#' @return sum of square (SS) of clustring and the 'delta' (difference of two successive majorization function).
#' @examples
#' {
#' X=rbind(matrix(rnorm(1000*2 ,4,.1),1000,2),matrix(rnorm(1000*2, 3, 0.2),1000,2))
#' M <- X[sample(nrow(X), 2),]
#' clusters_MajKm(X,2, 0.5)
#' }
#' @export

clusters_MajKm <- function(X, k=2, La) {
  old <- options("digits")
  M <- X[sample(nrow(X), 2),]
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  par(mfrow=c(1,2))
  M1=M
  plot(X)
  points(M1, lwd = 7, col=2)
  distsToCenters <- Euclid(X, M)
  clusters <- apply(distsToCenters, 1, which.min)
  Z <- matrix(0, nrow = NROW(X), ncol = 1)
  for(i in 1:NROW(X))
    if (clusters[[i]] == 1)
      Z[i,]=clusters[[i]]
  Z=cbind(Z, 1-Z)
  w_matrix <- matrix(0, nrow = NROW(X), ncol = NROW(X))
  w_matrix[cbind(seq(1, NROW(X) - 1),
                 seq(2, NROW(X)))] <- 1
  w_matrix[cbind(seq(2, NROW(X)),
                 seq(1, NROW(X) - 1))] <- 1
  W=w_matrix
  lambda=0.5
  n <- nrow(X)
  eps <- 1e-4
  k <- 1
  d <- 1
  g2=1
  while (k == 1 || d > 0.1) {
    D <- as.matrix(dist(Z%*%M))
    dim(Z)
    D <- pmax(D, eps)
    dim(W)
    V <- -W/D
    diag(V) <- -rowSums(V)
    ZX <- crossprod(Z, X)
    ZVZ <- crossprod(Z, crossprod(V, Z))
    G <- crossprod(Z,Z%*%M)-ZX+n*lambda*crossprod(ZVZ,M)
    M <- solve(crossprod(Z)+n*lambda*ZVZ, ZX)
    points(M, lwd = 7, col=5)
    distsToCenters <- Euclid(X, M)
    clusters <- apply(distsToCenters, 1, which.min)
    Z <- matrix(0, nrow = NROW(X), ncol = 1)
    for(i in 1:NROW(X))
      if (clusters[[i]] == 1)
        Z[i,]=clusters[[i]]
    Z=cbind(Z, 1-Z)
    z=Z
    if (NCOL(z)==2){
      mm=z[,1]+2*z[,2]
    }
    mat2=cbind(X,mm)
    M <- apply(X, 2, tapply, mm, mean)
    k <- k+1
    points(M, lwd = 7, col=8)

    #g <- sum(G^2)/n
    d=abs(g2-sum(G^2)/n)
    g2 <- sum(G^2)/n

  }
  newM=M
  points(newM, lwd = 7, col=6)
  z=Z
  if (NCOL(z)==2){
    mm=z[,1]+2*z[,2]
  }
  mat2=cbind(X,mm)
  centers_1 <- apply(X, 2, tapply, mm, mean)
  plot(X)
  points(mat2, col=mat2[,3])
  points(centers_1, lwd = 5, pch=c(2,2),col=c(5,4))
  lst1 <- setNames(lapply(split(1:nrow(mat2), mat2[,3]),
                          function(i) mat2[i,]), c("CL_1", "CL_2" ))

  CL_1=lst1$CL_1[,-3]
  CL_2=lst1$CL_2[,-3]
  CNT_1=apply(CL_1,2, mean)
  CNT_2=apply(CL_2,2, mean)
  CNT_1
  CNT_2
  s1<- matrix(NA, nrow=nrow(CL_1), ncol=1)
  for (i in 1:nrow(CL_1)){
    s1 [[i]]= (CL_1[i,1]-CNT_1[1])^2+(CL_1[i,2]-CNT_1[2])^2
  }
  s2<- matrix(NA, nrow=nrow(CL_2), ncol=1)
  for (i in 1:nrow(CL_2)){
    s2 [[i]]= (CL_2[i,1]-CNT_2[1])^2+(CL_2[i,2]-CNT_2[2])^2
  }
  s=sum(s1)+sum(s2)
  options(max.print=100)
  on.exit(options(old))
  list(clusters=clusters, centers=M )


}

