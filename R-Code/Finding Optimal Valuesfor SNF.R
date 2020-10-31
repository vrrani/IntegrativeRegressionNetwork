library(SNFtool)
library(analogue)
library(doParallel)
library(Hmisc)

setwd('~/SNF')
source('SimSilhouette.R')
source('getAffyMat.R')


GFLasso <- read.csv("D:/Research/BK21/BarPlot & VennDiagram/BetaFiles/Breast/Beta-GFLasso-Breast.csv",header = TRUE)
Lasso <- read.csv("D:/Research/BK21/BarPlot & VennDiagram/BetaFiles/Breast/Beta-Lasso-Breast.csv",header = TRUE)
SGL <- read.csv("D:/Research/BK21/BarPlot & VennDiagram/BetaFiles/Breast/Beta-SGL-Breast.csv",header = TRUE)
SIOL <- read.csv("D:/Research/BK21/BarPlot & VennDiagram/BetaFiles/Breast/Beta-SIOL-Breast - 0.1.csv",header = TRUE)

GFLasso <- GFLasso[,-1]
Lasso <- Lasso[,-1]
SGL <- SGL[,-1]
SIOL <- SIOL[,-1]

input <- list( GFLasso, Lasso, SGL,SIOL )
names( input ) <- c( "GFLasso", "Lasso", "SGL", "SIOL")

findCorrelation <- function (W0,W1,W2,W3,W4) {
   finalCorrelation <- (SIM(W0,W1)+SIM(W0,W2)+SIM(W0,W3)+SIM(W0,W4)+1-SIM(W1,W2)+1-SIM(W1,W3)+1-SIM(W1,W4)+1-SIM(W2,W3)+1-SIM(W2,W4)+1-SIM(W3,W4))/ 10
   return(finalCorrelation)   
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
    
iter <- 20
 for(K in 2:20) {
  for(alpha in seq(0.3,0.8, by = 0.1))  {
      WS <- list()
      for( name in names(distMatrix) ) {
          D <- distMatrix[[name]]
          WS[[name]] <- getAffyMat(D, K, alpha)
      }

      W <- SNF( WS, K, iter )
      Correlation <- findCorrelation(W, WS[[1]], WS[[2]], WS[[3]], WS[[4]])
      result <- sprintf('%d\t%f\t%f', K, alpha, Correlation)
      print(result)
      write(result, file = "SNF_20150321.txt", append = TRUE)
  }
 }
    
        
