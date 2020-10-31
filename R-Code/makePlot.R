maxIter <- 100

for( method in c("fused", "GFLasso", "Lasso", "SGL", "SIOL") ) {
    maxX <- 0

    for( i in 1:length( Wset ) ) {
        X <- Wset[[i]][[method]]
        maxX <- max( maxX, max(X[upper.tri(X)]) )
    }

    print( sprintf("%s %f", method, maxX) )

    pdf( sprintf('figure_Breast_%s.pdf',method),width=9,height=13)
    par( mfrow=c(2,1) )

    iv <- seq( 0, maxX, maxX / 40.0 )
    ivAxis <- seq( 0, maxX, maxX / 20.0 )
    plot( 0, type='n', xlim=c(0,maxX), ylim=c(0,log(2*10^5,10)), xaxt="n", ylab="log10(num. edges)", xlab="weight", title = method )
    axis(1, at = sprintf("%.4f",iv) , las=2, cex.axis=0.8)

    for( i in 2:maxIter ) {
        lines( iv, log( getNumEdges( Wset[[i]][[method]], iv ), 10 ), col='gray' )
    }
    lines( iv, log( getNumEdges(Wset[[1]][[method]], iv),10), lwd=1.5, col='red')

    plot( 0, type='n', xlim=c(maxX,0), ylim=c(600,0), xaxt="n", ylab="num. component", xlab="weight", title = method  )
    axis(1, at = sprintf("%.4f",iv) , las=2, cex.axis=0.8)

    for( i in 2:maxIter ) {
        lines( iv, getMaxComponent( Wset[[i]][[method]], iv ), col='gray' )
    }
    lines( iv, getMaxComponent(Wset[[1]][[method]], iv), lwd=1.5, col='red')
    dev.off()
}
