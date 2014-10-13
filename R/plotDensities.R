## Plot of the densities
densitiesPlot <- function(matrix.x, matrix.y, quantiles,
                          sampleNames, bw, main, xlab, xlim,
                          ylim, lty = 1, lwd =1, mean = FALSE,
                          col, from, to, bty = "o"){
    meanSample <- apply(matrix.x,1,mean)
    plot(matrix.x[,1] , matrix.x[,2], 
         main = main, 
         ylab = "Density", 
         xlab = xlab,
         col = "white",
         xlim = xlim, 
         ylim = ylim, bty = bty
         )
    ##rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "black")

    ## Will plot individual densities        
    if (!mean){
        ##Wanna make sure colors are not randomly drawn
        o <- order(col, decreasing=FALSE)
        col <- col[o]
        matrix.x <- matrix.x[,o]
        matrix.y <- matrix.y[,o]
        
        
        for (j in 1:ncol(matrix.x)){
            lines(matrix.x[,j], matrix.y[,j],
                  col= col[j], lwd = lwd, lty=lty)
        }
    } else {
        ## Make sure to remove missing values
        non.na.indices <- which(!is.na(col))
        col  <- col[non.na.indices]
        quantiles  <- quantiles[, non.na.indices]
        
        unique.col <- unique(col)
        n <- length(unique.col)
        for (i in 1:n){
            indices <- which(col==unique.col[[i]])
            current.color <- col[indices][1]
            if (length(indices==1)){
                currentMean <- quantiles[,indices]
            } else {
                currentMean <- apply(quantiles[,indices],1,mean)
            }
            lines(density(currentMean, bw = bw, from= from, to = to),
                  lwd = lwd, lty=lty, col=current.color)
        }
    }
}

## To add selected density
addHoverDensity <- function(selectedSamples = c(), sampleNames,
                            matrix.x, matrix.y, col){
    n <- length(selectedSamples)
    if (n >= 1){
        for (kk in 1:n){
            lines(matrix.x[,selectedSamples[kk]],matrix.y[,selectedSamples[kk]],
                  col="black",lwd = 4, lty=3)
        }
        lines(matrix.x[,selectedSamples[n]],matrix.y[,selectedSamples[n]],
              col="black",lwd = 3, lty=1)
    }
}
