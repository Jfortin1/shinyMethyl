## Return gender covariate if provided in the covariates
returnGivenGender <- function(covariates){
    givenGender <- NULL
    
    if (!is.null(covariates)){
        possibilities <- c("gender","Gender","sex","Sex","GENDER","SEX")
        sum <- sum(possibilities %in% colnames(covariates))
        if (sum!=0){
            goodColumn <- possibilities[possibilities %in% colnames(covariates)][1]
            givenGender <- as.character(covariates[,goodColumn])
            givenGender <- substr(toupper(givenGender),1,1)
            givenGender <- as.character(givenGender)
        }}
    
    return(givenGender)
}

## Predict gender with X and Y chr intensities
returnPredictedGender <- function(cnQuantiles, cutoff = (-0.3)){
    ## It assumes that cnQuantiles is a list containing $X and $Y
    X <- cnQuantiles$X
    Y <- cnQuantiles$Y
    med <- floor(nrow(X)/2)
    
    x = X[med,]
    y = Y[med,]
    diff <- log2(y)-log2(x)
    n <- length(diff)
    
    predictedGender = rep("M", n)
    predictedGender[which(diff < cutoff)] <- "F"     
    names(predictedGender) <- colnames(X)   
    
    predictedGender <- as.character(predictedGender)
    return(predictedGender)
}

plotPredictedGender <- function(cnQuantiles,cutoff =(-0.3),
                                color = NULL, legend=TRUE, bty="o"){
    X <- cnQuantiles$X
    Y <- cnQuantiles$Y
    med <- floor(nrow(X)/2)
    x <- X[med,]
    y <- Y[med,]
    diff <- log2(y)-log2(x)
    n <- length(diff)
    
    if (is.null(color)){
        color <-  rep("lightskyblue", n)
        color[which(diff < cutoff)] <- "orange"
    }
    
    plot(diff, 
         jitter(rep(0,n),factor=1.5), 
         ylim = c(-1,1), 
         pch=18, 
         cex=2,
         col= color, 
         yaxt="n", 
         xlab="median CN(Y)  - median CN(X)", 
         ylab="", bty=bty
         )            
    abline(v=as.numeric(cutoff),lty=3,lwd=2)
    
    
    if (legend){
        legend("topright",
               c("Predicted male","Predicted female"),
               cex=2, 
               pch=18, 
               col=c("lightskyblue","orange"),
               bty="n"
               )
    }
    
}

plotDiscrepancyGenders <- function(cutoff = (-0.3), covariates,
                                   cnQuantiles, bty="o"){
    predictedGender <- returnPredictedGender(cutoff = cutoff,
                                             cnQuantiles = cnQuantiles)
    givenGender     <- returnGivenGender(covariates)
    ## In the case a gender information was not provided:
    if (is.null(givenGender)){
        plotPredictedGender(cnQuantiles = cnQuantiles, 
                            cutoff = cutoff, 
                            color = NULL,
                            legend = TRUE, bty = bty)
    } else {
    ## In the case a gender information WAS provided:
        X <- cnQuantiles$X
        Y <- cnQuantiles$Y
        med <- floor(nrow(X)/2)
        x = X[med,]
        y = Y[med,]
        diff <- log2(y)-log2(x)
        n <- length(diff)
        
        color <- predictedGender
        color[color == "M"] <- "lightskyblue"
        color[color == "F"] <- "orange"
        
        ## For the samples with different predicted gender:
        non.matching <- predictedGender  != givenGender
        color[non.matching] <- "black"
        
                                        # For the samples with missing values:
        na.positions <- is.na(givenGender)
        color[na.positions] <- "grey"
        
        plotPredictedGender(cnQuantiles = cnQuantiles, 
                            cutoff = cutoff, 
                            color = color,
                            legend = FALSE, bty=bty)
        
        ## To add the sample names in the plot for the discrepancy samples:
        ## if (sum(color=="black" & color!="grey")>=1){
        ## indices <- which(color=="black" & color != "grey")
        ## text(diff[indices],-0.2, names(X)[indices])
        
        ## To add the legend: 
        legend.colors <- c("lightskyblue","orange")
        legend.names  <- c("Predicted male","Predicted female")
        
        ## In the case there are unmatching samples:
        if (sum(color=="black" & color!="grey")>=1){
            legend.colors <- c(legend.colors, "black")
            legend.names  <- c(legend.names, "Unmatching samples")
        }
        
        ## In the case there are missing values in the provided sex covariate:
        if (sum(color=="grey")>=1){
            legend.colors <- c(legend.colors, "grey")
            legend.names  <- c(legend.names, "Sex not provided by user")
        }
        
        legend( "topright",
               legend = legend.names,
               col = legend.colors,
               cex = 1.5, 
               pch = 18,
               bty="n")
    }
}

