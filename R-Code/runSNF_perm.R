library(SNFtool)
library(doParallel)
library(permute)
library(igraph)
setwd('~/SNF')
source('SimSilhouette.R')
source('getAffyMat.R')

GFLasso <- read.csv("C:/Users/USER/Documents/Research/BK21/BarPlot & VennDiagram/BetaFiles/Breast/Beta-GFLasso-Breast.csv",header = TRUE)
Lasso <- read.csv("C:/Users/USER/Documents/Research/BK21/BarPlot & VennDiagram/BetaFiles/Breast/Beta-Lasso-Breast.csv",header = TRUE)
SGL <- read.csv("C:/Users/USER/Documents/Research/BK21/BarPlot & VennDiagram/BetaFiles/Breast/Beta-SGL-Breast.csv",header = TRUE)
SIOL <- read.csv("C:/Users/USER/Documents/Research/BK21/BarPlot & VennDiagram/BetaFiles/Breast/Beta-SIOL-Breast - 0.1.csv",header = TRUE)

GFLasso <- GFLasso[,-1]
Lasso <- Lasso[,-1]
SGL <- SGL[,-1]
SIOL <- SIOL[,-1]

input <- list( GFLasso, Lasso, SGL,SIOL )
names( input ) <- c( "GFLasso", "Lasso", "SGL", "SIOL")


permuteMat <- function( Y ) {
    for( i in 1:nrow(Y) ) {
        Y[i,] <- Y[i,shuffle(ncol(Y))]
    }
    return(Y)
}


getNumEdges <- function( G, interval ) {
    fq <- vector()
    for( t in interval ) {
        fq <- cbind( fq, sum( G[upper.tri(G)]>t ) )
    }  
    return(fq)
}

getMaxComponent <- function( G, interval ) {
     fq <- vector()
     for( t in interval ){
     
      WG <- graph.adjacency(G>t, mode="undirected")      
      c <- clusters( WG ) 
      maxc <- max( c$csize )
       fq <- cbind( fq, maxc ) 
     }  
     return(fq)
 }
 
runSNF <- function( ind, input, needPermute ) {
  alpha <- 0.4 
  K <- 15
  iter <- 20
  library(SNFtool)
  library(permute)

  print( sprintf("%d-%d", ind, needPermute) )
  distMatrix <- list()
  for( name in names(input) ) {
    X <- as.matrix(input[[name]])
   
    if( needPermute ) {
        X <- permuteMat( X )
    }
    SX <- standardNormalization( X ) 
    distMatrix[[name]] <- dist2( SX, SX )
  }
  
  WS <- list()
  
  for( name in names(distMatrix) ) {
    D <- distMatrix[[name]]
    WS[[name]] <- getAffyMat(D, K, alpha)
  }
  
  WF <- SNF( WS, K, iter )
  WS[["fused"]] <- WF

  return(WS)
}

registerDoParallel(core=4)
maxIter <- 100
ptime <- system.time({
    Wset <- foreach( i=0:maxIter ) %dopar% runSNF(i,input,i!=0)
})

print( ptime )

# Code for findinga cut off
finalcutoff<- list()

for( method in c("fused", "GFLasso", "Lasso", "SGL", "SIOL") ) {
    finalEdgescount<- NULL	
    finalComponentscount <- NULL
    print( method )
    maxW <- 0    
    for( i in 1:101 ){
        X <- Wset[[i]][[method]]
        maxW = max( maxW, max(X[upper.tri(X)]) )
    }
    iv <- seq( 0, maxW, maxW/40.0 )

    for( i in 2:101 ){
        getEdgescount <- getNumEdges(Wset[[i]][[method]], iv)
        finalEdgescount <- rbind(finalEdgescount, getEdgescount)

        getComponentscount <- getMaxComponent(Wset[[i]][[method]], iv)
        finalComponentscount <- rbind(finalComponentscount, getComponentscount)
    }
    numerator1 <- colMeans(finalEdgescount)
    denom1 <- getNumEdges(Wset[[1]][[method]], iv)

    numerator2 <- colMeans(finalComponentscount)
    denom2 <- getMaxComponent(Wset[[1]][[method]], iv)

    rat <- 0.5*(numerator1/denom1 + numerator2/denom2)
    print(rat)
    cutoff <- iv[which(rat==min(rat, na.rm=TRUE))]
    if(cutoff == 0) {
        nextmin <-sort(rat)[2]
        cutoff <- iv[which(rat== nextmin)] 
    }
    finalcutoff[[method]] <- cutoff 
}

finalcutoff

