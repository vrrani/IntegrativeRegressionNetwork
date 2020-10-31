library(SNFtool)
library(cluster)
source('SimSilhouette.R')
source('getAffyMat.R')
setwd('~/SNF')

cancer <- c( "Breast", "Colon", "GBM" )

write("cancer\tK\talpah\tC\tsi", file = "SNF_2015022.txt", append = FALSE)

for( cname in cancer ) {
    GFL <- read.csv( sprintf('%s/Beta-GFLasso-%s.csv', cname, cname), row.names=1 )
    Lasso <- read.csv( sprintf('%s/Beta-Lasso-%s.csv', cname, cname), row.names=1 )
    SGL <- read.csv( sprintf('%s/Beta-SGL-%s.csv', cname, cname), row.names=1 )

    input <- list( GFL, Lasso, SGL )
    names( input ) <- c( "GFL", "Lasso", "SGL" )

    distMatrix <- list()

    for( name in names(input) ) {
        X <- t(as.matrix(input[[name]]))
        SX <- standardNormalization( X ) # this function is to normalize data, not to normalize distance
        distMatrix[[name]] <- dist2( X, X )
    }


    iter <- 20

    for( K in 2:20 ) {
      for( alpha in seq(0.3,0.8,by=0.1) ) {
        for( C in 2:20 ) {
          
          WS <- list()
          
          for( name in names(distMatrix) ) {
              D <- distMatrix[[name]]
              WS[[name]] <- getAffyMat(D, K, alpha)
          }
          
          W <- SNF( WS, K, iter )
          group <- spectralClustering( W, C )
          si <- SimSilhouette( W, group )

          result <- sprintf('%s\t%d\t%f\t%d\t%f', cname, K, alpha, C, mean(si) )
          print(result)
          write(result, file = "SNF_2015022.txt", append = TRUE)      
        }
      }
    }
}
