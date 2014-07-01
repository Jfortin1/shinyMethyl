### To plot the internal controls
plotInternalControls <- function(y, col, main, sampleNames, bty="o"){
    n <- length(y)
    slideNames <- substr(sampleNames, 1, 10)
    plot(1:n, y, pch=20, cex=0.7, main = main, xlab = "Sample",
         ylab = "Control intensity", ylim = c(0, 1.3*max(y)),
         col = col, bty = bty)
    grid()
    abline(h=0, lty=3)
}

### To add the selected array to the internal controls plot
addHoverPoints <- function(y, sampleNames, selectedSamples=c()){
    ## To put a circle around the last entry: 
    n <- length(selectedSamples)
    if (n>=1){
        points(selectedSamples[n], y[selectedSamples[n]], pch = 1,
               cex = 3, lwd = 3)
        points(selectedSamples[1:n], y[selectedSamples[1:n]],
               pch = 17, cex = 1, lwd = 3)
    }
}

returnControlStat <- function(controlType, greenControls,
                              redControls, controlNames){
    index <- match(controlType, controlNames)
    greenControls.current <- greenControls[[index]]
    redControls.current <- redControls[[index]]
    
    if (controlType=="BISULFITE CONVERSION I"){
        controlStat <- (apply(greenControls.current[1:3,], 2, mean) +
                        apply(redControls.current[7:9,], 2, mean))/2 
    } else {
        if (controlType=="BISULFITE CONVERSION II"){
            controlStat <- apply(redControls.current, 2, mean)
        } else {
            if (controlType=="EXTENSION"){
                controlStat <- (apply(greenControls.current[3:4,], 2, mean) +
                                apply(redControls.current[1:2,], 2, mean))/2 
            } else {
                if (controlType=="NEGATIVE"){
                    controlStat <- (apply(greenControls.current, 2, mean) +
                                    apply(redControls.current, 2, mean))/2 
                } else {
                    if (controlType=="HYBRIDIZATION"){
                        controlStat <- apply(greenControls.current, 2, mean)
                    } else {controlStat <- apply(greenControls.current, 2, mean)
                        }
                }
            }
        }
    }
    return(controlStat)
}
