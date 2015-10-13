plotDesign450k <- function(sampleNames, covariates, legendTitle){
    boxColors <- buildColors(sampleNames, covariates)
    nLevels <- length(unique(boxColors))
    naPos <- boxColors=="NA"
    boxColors <- factor(boxColors)
    oldLevels <- levels(boxColors)
    levels(boxColors) <- 1:(nLevels)
    boxColors <- as.character(boxColors)
    if (sum(naPos)>0){
        boxColors[naPos] <- "white"
    }
    
    numberOfChips <- length(unique(substr(sampleNames,1,10)))
    nrows <- ceiling(numberOfChips/8)
    par(pin=c(width=5,heigth=5),mai=c(0,0,1,0))
    plot(0:(nrows+1), seq(0,nrows,nrows/(nrows+1)),
         type="n",axes=FALSE,xlab="",ylab="",
         main = "Illumina HumanMethylation 450K Array Design")
    
    for (i in 1:nrows){
        plotRow(leftMargin= 0, plateytop=(i+0.75)-1, plateybottom=i-1,
                color = boxColors[(96*(i-1)+1):(96*i)])
    }
    legendColors <- 1: (nLevels)
    legendPch <- rep(15, nLevels)
    if (sum(naPos)>0){
        legendColors[oldLevels=="NA"] <- "black"
        legendPch[oldLevels=="NA"] <- 0
    }
}

plotLegendDesign450k <- function(sampleNames, covariates, legendTitle){
    boxColors <- buildColors(sampleNames, covariates)
    nLevels <- length(unique(boxColors))
    naPos <- boxColors=="NA"
    boxColors <- factor(boxColors)
    oldLevels <- levels(boxColors)
    levels(boxColors) <- 1:(nLevels)
    boxColors <- as.character(boxColors)
    if (sum(naPos)>0){
        boxColors[naPos] <- "white"
    }
    
    legendColors <- 1: (nLevels)
    legendPch <- rep(15, nLevels)
    
    if ("NA" %in% oldLevels){
        legendColors[oldLevels=="NA"] <- "white"
    }
    
    plot(0,xaxt='n',yaxt='n',bty='n',pch='',ylab='',xlab='')
    legend("topright", pch= legendPch, col = legendColors,
           oldLevels, cex=1.5, title= legendTitle)  
}

## To build the color palette for the rectangles
buildColors <- function(sampleNames, colCovariate){
    colCovariate[is.na(colCovariate)] <- "NA"
    colCovariate <- as.character(colCovariate)
    chipInfo <- substr(sampleNames,1,10)
    rowInfo  <- substr(sampleNames,14,14)
    columnInfo <- substr(sampleNames,17,17)
    numberOfChips <- length(unique(chipInfo))
    boxColors <- rep("NA",numberOfChips*12)
    
    for (i in 1:numberOfChips){
        for (columnIndex in 1:2){
            for (rowIndex in 1:6){
                indices <- intersect(which(chipInfo == unique(chipInfo)[i]),
                                     which(columnInfo == columnIndex))
                indices <- intersect(indices, which(rowInfo == rowIndex))
                if (!length(indices)==0){
                    boxColors[12*(i-1)+1+6*(columnIndex-1)+(rowIndex-1)] <- colCovariate[indices]
                }
            }
        }
    }
    
    return(boxColors)
}


### To plot a row of plates
plotRow <- function(leftMargin, plateytop, plateybottom, color){
    unit <- (plateytop - plateybottom)/6
    rightMargin = 1+23*unit
    platePosition <- seq(leftMargin,rightMargin, (rightMargin-leftMargin)/(7))
    for (j in 1:length(platePosition)){
        plotPlate(plateybottom=plateybottom, plateytop=plateytop,
                  platePosition[j], color[(12*(j-1)+1):(12*j)])
    }
}


## To plot a single array
plotPlate <- function(plateybottom, plateytop, platexleft, color){
    platexright <- platexleft + (plateytop - plateybottom)/3
    rect(ybottom = plateybottom, ytop = plateytop, 
         xleft= platexleft, xright = platexright, 
         lwd=2, col=0)
    
    startsColumn1 <- sort(seq(plateybottom, plateytop,
                              (plateytop-plateybottom)/6), decreasing=TRUE)[-1]
    endsColumn1 <-  sort(seq(plateybottom, plateytop,
                             (plateytop-plateybottom)/6), decreasing=TRUE)[-7]
    
    ## For plotting the left column:
    for (k in 1:length(startsColumn1)){
        rect(ybottom = startsColumn1[k], 
             ytop=endsColumn1[k], 
             xleft = platexleft, 
             xright = (platexright+platexleft)/2, 
             col=color[k],
             lwd=3)
    }
    
    ## For plotting the right column:
    for (v in 1:length(startsColumn1)){
        rect(ybottom = startsColumn1[v], 
             ytop=endsColumn1[v],
             xleft = (platexleft+platexright)/2, 
             xright = platexright, 
             col = color[6+v],
             lwd=3)
    }
}
