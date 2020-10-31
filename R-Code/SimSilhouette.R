# D: distance n \times n matrix
# C: vector for lable of each point
# reference: http://ac.els-cdn.com/0377042787901257/1-s2.0-0377042787901257-main.pdf?_tid=1cdd6fd4-b9a4-11e4-8dec-00000aacb360&acdnat=1424507703_0ddf35794da4a78e4c09da2ad6a85f03

SimSilhouette <- function( D, C ) {
    stopifnot( nrow(D) == ncol(D) )
    stopifnot( nrow(D) == length(C) )

    clusters <- unique(C)
    N <- length(C) # number of points
    maxC <- max(C) 
    S <- rep(0, N)
    A <- rep(0, N)
   
    for( i in 1:N ) {
        inSameC <- which(C == C[i])
        inSameC <- inSameC[inSameC!=i] # discard itself

        avgSimInSameC <- mean( D[i,inSameC] ) # this is a(i)
        if( is.nan(avgSimInSameC) ) {
            avgSimInSameC <- 0.0
        }

        avgSimInDiffC <- rep(-Inf,maxC)

        for( c in clusters ) {
            if( c == C[i] ) {
                next
            }
            inDiffC <- which(C == c)
            avgSimInDiffC[c] <- mean( D[i,inDiffC] ) 
        }
        maxSimInDiffC <- max( avgSimInDiffC ) # this is b(i)
       
        # print( c( i, avgSimInSameC, maxSimInDiffC ) )
        if( avgSimInSameC > maxSimInDiffC ) {
            S[i] <- 1 - maxSimInDiffC / avgSimInSameC
        }
        else if( avgSimInSameC < maxSimInDiffC ) {
            S[i] <- avgSimInSameC / maxSimInDiffC - 1
        }
    }

    return(S)
}
