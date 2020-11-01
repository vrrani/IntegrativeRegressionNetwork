library(SNFtool)
library(doParallel)
setwd('~/SNF')
source('SimSilhouette.R')
source('getAffyMat.R')

GFLasso <- read.csv("~/Lung/Beta-GFLasso-Lung.csv",header = TRUE)
Lasso <- read.csv("~/Lung/Beta-Lasso-Lung.csv",header = TRUE)
SGL <- read.csv("~/Lung/Beta-SGL-Lung.csv",header = TRUE)
SIOL <- read.csv("~/Lung/Beta-SIOL-Lung - 0.1.csv",header = TRUE)

GFLasso <- GFLasso[,-1]
Lasso <- Lasso[,-1]
SGL <- SGL[,-1]
SIOL <- SIOL[,-1]

#input <- list( (GFLasso), (Lasso), (SGL), (SIOL) )
input <- list( t(GFLasso), t(Lasso), t(SGL),t(SIOL) )
names( input ) <- c( "GFLasso", "Lasso", "SGL", "SIOL")

findCorrelation <- function (W0,W1,W2,W3,W4) {
  
  WLst <- list( W0, W1, W2, W3, W4)
  finalCorrelation <- (SIM(W0,W1)+SIM(W0,W2)+SIM(W0,W3)+SIM(W0,W4)+1-SIM(W1,W2)+1-SIM(W1,W3)+1-SIM(W1,W4)+1-SIM(W2,W3)+1-SIM(W2,W4)+1-SIM(W3,W4))/ 10
  
  corMat <- matrix(0,5,5)
  for( i in 1:4 ) {
    for( j in (i+1):5 ) {
      corMat[i,j] <- SIM( WLst[[i]], WLst[[j]])
    }
  }
  
  fused <- sum( abs( corMat[1,2:5] ))
  single <- 6-sum( abs( corMat[2:5,2:5]) )
  return( ( fused + single ) / 10.0 )
  # return( 0.5*(fused/4+single/6) )   
}

SIM <- function(x,y) {
    ux <- x[upper.tri(x)]
    uy <- y[upper.tri(y)]
    return( cor(ux,uy) )
}

distMatrix <- list()
for( name in names(input) ) {
  X <- as.matrix(input[[name]])
  SX <- standardNormalization( X ) # this function is to normalize data, not to normalize distance        
  distMatrix[[name]] <- dist2( SX, SX )
}

runSNF <- function( distMatrix, K, alpha ) {
library(SNFtool)
    iter <- 20
    WS <- list()
    for( name in names(distMatrix) ) {
        D <- distMatrix[[name]]
        WS[[name]] <- getAffyMat(D, K, alpha)
    }
    W <- SNF( WS, K, iter )
    Correlation <- findCorrelation(W, WS[[1]], WS[[2]], WS[[3]], WS[[4]])
    return(Correlation)
}

# runSNF( distMatrix, 4, 0.4 )

retMat <- matrix( 0, 20, 6)
for( K in 2:20 ) {
    registerDoParallel(core=6)
    result <- foreach(alpha=c(0.3,0.4,0.5,0.6,0.7,0.8)) %dopar% runSNF(distMatrix, K, alpha)
    print( unlist(result) )
    retMat[K,] <- unlist(result)
}
which(retMat == max(data.matrix(retMat)), arr.ind = TRUE)
max(data.matrix(retMat))
