plotPCA <- function(pca, pc1, pc2, col, covariates, selectedCov, bty="o"){
    xMin <- min(pca[,as.numeric(pc1)])
    xMax <- max(pca[,as.numeric(pc1)])
    xRange <- xMax - xMin
    xlim <- c(xMin-0.05*xRange, xMax+0.20*xRange)
    xlab <- paste("PC",as.numeric(pc1), " scores", sep="")
    ylab <- paste("PC",as.numeric(pc2), " scores", sep="")
    
    plot(pca[,as.numeric(pc1)], pca[,as.numeric(pc2)],
         col = col, pch = 18, cex = 2, xlab = xlab,
         ylab = ylab, xlim = xlim,
         main = "Principal component analysis (PCA)",
         cex.main = 1.5, cex.lab = 1.5, bty = bty)
    uColor <- unique(col)
    uCov   <- unique(covariates[,match(selectedCov, colnames(covariates))])
    
    legend("bottomright", legend = uCov, pch = 18, col = uColor,
           cex = 1.5, title = selectedCov, bty = "n")
    grid()
}
