plotQC <- function(unmethQuantiles, methQuantiles, 
                   sampleNames, col, bty="o"){
    
    slideNames <- substr(sampleNames,1,10)
    
    med <- as.integer(nrow(methQuantiles[[3]])/2)
    x <- log2(unlist(unmethQuantiles[[3]][med,]))
    y <- log2(unlist(methQuantiles[[3]][med,]))
    
    range.x <- max(x)-min(x)
    range.y <- max(y)-min(y)
    xlim1 <- min(x)-0.2*range.x
    xlim2 <- max(x)+0.2*range.x
    ylim1 <- min(y)-0.2*range.y
    ylim2 <- max(y)+0.2*range.y
    
    plot(x, y, xlim = c(xlim1, xlim2), ylim = c(ylim1, ylim2),
         cex = 1, pch = 20, col = col, main = "QC Plot", bty=bty)
    grid()
}

addHoverQC <- function(y, selectedSamples = c(),
                       unmethQuantiles, methQuantiles){
    med <- as.integer(nrow(methQuantiles[[3]])/2)
    mediansU <- unlist(unmethQuantiles[[3]][med,])
    mediansM <- unlist(methQuantiles[[3]][med,])
    
    ## To put a circle around the last entry: 
    n <- length(selectedSamples)
    if (n>=1){
        points(log2(mediansU[selectedSamples[n]]), 
               log2(mediansM[selectedSamples[n]]),
               col = "black", cex=3, pch=1, lwd=2)
    }
    ## Make the points black:
    points(log2(mediansU[selectedSamples]), 
           log2(mediansM[selectedSamples]),
           col = "black", cex = 1, pch = 17)
}
