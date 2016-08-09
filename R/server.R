#### Created by Jean-Philippe Fortin
#### March 28, 2014

server.shinyMethyl <- function(shinyMethylSet1, shinyMethylSet2=NULL){
    function(input, output, session) { 
        betaQuantiles   <-  getBeta(shinyMethylSet1)
        mQuantiles      <-  getM(shinyMethylSet1)
        methQuantiles   <-  getMeth(shinyMethylSet1)
        unmethQuantiles <-  getUnmeth(shinyMethylSet1)
        cnQuantiles     <-  getCN(shinyMethylSet1)
        greenControls   <-  getGreenControls(shinyMethylSet1)
        redControls     <-  getRedControls(shinyMethylSet1)
        covariates      <<- pData(shinyMethylSet1)
        pca             <-  getPCA(shinyMethylSet1)$scores
        sampleNames     <-  sampleNames(shinyMethylSet1)
        slideNames      <-  substr(sampleNames,1,10)
        arrayNames      <-  substr(sampleNames,12,17)
        plateNames      <-  substr(sampleNames,1,6)
        controlNames    <-  names(greenControls)
        
        method <- shinyMethylSet1@originObject
        
        
        ## In the case covariates is empty:
        if (ncol(covariates)==0){
            covariates <- data.frame(slide = slideNames, plate = plateNames)
            rownames(covariates) <- sampleNames
            covariates <<- covariates
        }
        
        ## Global variables:
        ## that will be accessed by the knitr report
        mouse.click.indices <<- c()
        colorSet     <<- "Set1" # Default color set
        sampleColors <<- as.numeric(as.factor(plateNames)) ## Default sample colors
        genderCutoff <<- -0.4 # Default cutoff for gender prediction
        current.control.type <<- "BISULFITE CONVERSION I"
        current.probe.type   <<- "I Green"
        current.density.type <<- "M-value"
        
        
        
#########################################################
        
###########        Colors 
        
#########################################################
        
        
        
        
                                        # To change the global color scheme:
	setColor <- reactive({
            colorSet <<- input$colorChoice
	})
        
	set.palette <- function(n, name){
            ## The name of the colors are part of the RColorBrewer package
            if (name != "Pault" & name != "Rainbow"){
                palette(brewer.pal(n = n, name = name ))
                
                ## Custom palette
            } else if (name == "Pault"){
                colors <- c("#332288","#88CCEE","#44AA99","#117733","#999933",
                            "#DDCC77","#661100","#CC6677","#882255","#AA4499")
                colors <- c("#4477AA","#CC6677","#DDCC77","#117733","#88CCEE",
                            "#AA4499","#44AA99","#999933","#882255","#661100",
                            "#6699CC", "#AA4466")
                palette(colors)
            } else if (name == "Rainbow"){
                colors <- c("#781C81", "#3F56A7","#4B91C0","#5FAA9F","#91BD61",
                            "#D8AF3D","#E77C30","#D92120")
                palette(colors)
            }
	}
        ## To choose the colors according to the phenotype:
	sampleColors <- reactive({
            return(as.numeric(as.factor(covariates[,match(input$phenotype,
                                                          colnames(covariates))])))
	})

#########################################################
###########        Computation of the densities
#########################################################

	## Return x and y matrices of the densities for raw data
	returnDensityMatrix <- reactive({
            index <- match(input$probeType,c("I Green","I Red","II","X","Y"))
            if (input$mOrBeta=="Beta-value"){
                bw <- input$bandwidth 
                quantiles <- betaQuantiles[[index]]
            } else {
                bw <- input$bandwidth2
                quantiles <- mQuantiles[[index]]
            }
            
            matrix <- apply(quantiles, 2, function(x){
                d <- density(x, bw = bw, n =512)
                c(d$x,d$y)
            })
            matrix.x <- matrix[1:512,]
            matrix.y <- matrix[513:(2*512),]
            return(list(matrix.x = matrix.x, matrix.y = matrix.y))
	})
        
    ## Return x and y matrices of the densities for normalized data
	returnDensityMatrixNorm <- reactive({
            if (is.null(shinyMethylSet2)){
                return(NULL)
            } else {
                index <- match(input$probeType,
                               c("I Green","I Red","II","X","Y"))
            }
            if (input$mOrBeta=="Beta-value"){
                bw <- input$bandwidth 
                quantiles <- shinyMethylSet2@betaQuantiles[[index]]
            } else {
                bw <- input$bandwidth2
                quantiles <- shinyMethylSet2@mQuantiles[[index]]
            }
            matrix <- apply(quantiles, 2, function(x){
                d <- density(x, bw = bw, n =512)
                c(d$x,d$y)
            })
            matrix.x <- matrix[1:512,]
            matrix.y <- matrix[513:(2*512),]
            return(list(matrix.x = matrix.x, matrix.y = matrix.y))	
	})

#########################################################
###########         Mouse clicks section
#########################################################

    ## Check if the selected sample is already in the list or not:
	check.mouse.clicks <- function(clickedIndex, mouse.click.indices){
           
            if (length(mouse.click.indices)!=0){
                if (clickedIndex %in% mouse.click.indices){
                    mouse.click.indices <- mouse.click.indices[mouse.click.indices!=clickedIndex]
                } else {
                    mouse.click.indices <- c(mouse.click.indices, clickedIndex)
                }
            } else {
                mouse.click.indices <- c(mouse.click.indices, clickedIndex)
            }

            return(mouse.click.indices)
	}
    
    ## To update the list of mouse clicks from internal controls plot:
	updateMouseClicks <- reactive({
            
        mouse.x <- input$controlsHover$x
        mouse.y <- input$controlsHover$y
        if (!is.null(mouse.x) & !is.null(mouse.y)){
            y <- reactiveControlStat()  
            n <- length(y)
            xDiff <- ((1:n)-rep(mouse.x, n)) / n
            yDiff <- (y - mouse.y)/(1.3*max(y))
            clickedIndex <- which.min(xDiff^2 + yDiff^2)
            names(clickedIndex) <- shinyMethylSet1@sampleNames[clickedIndex]
            
         
            if (current.control.type == as.character(input$controlType)){
                mouse.click.indices <<- check.mouse.clicks(clickedIndex, mouse.click.indices)
            }
            current.control.type <<- as.character(input$controlType)
        }
        mouse.click.indices
	})
    ## To update the list of mouse clicks from quality control plot:
	updateMouseClicksQC <- reactive({
        mouse.x <- input$qualityHover$x
        mouse.y <- input$qualityHover$y
        
        if (!is.null(mouse.x) & !is.null(mouse.y)){
            med <- as.integer(nrow(methQuantiles[[3]])/2)
    	    mediansU <- unlist(unmethQuantiles[[3]][med,])        
    	    mediansM <- unlist(methQuantiles[[3]][med,]) 
                
            x <- log2(mediansU)
            y <- log2(mediansM)
            n <- length(y)
            
            range.x <- max(x) - min(x)
            range.y <- max(y) - min(y)
            
            xDiff <- (x - mouse.x) / (1.4*range.x)
            yDiff <- (y - mouse.y)/  (1.4*range.y)
            clickedIndex <- which.min(xDiff^2 + yDiff^2)
            names(clickedIndex) <- shinyMethylSet1@sampleNames[clickedIndex]
            mouse.click.indices <<- check.mouse.clicks(clickedIndex,
                                                       mouse.click.indices)
        }
        mouse.click.indices
	})
        
    ## To update the list of mouse clicks from densities plot:
	updateMouseClicksDensities <- reactive({
        mouse.x <- input$densitiesHover$x
        mouse.y <- input$densitiesHover$y
        #print(mouse.x)
        #print(mouse.y)
        
        if (!is.null(mouse.x) & !is.null(mouse.y)){
        # If Beta values are selected
            if (input$mOrBeta=="Beta-value") {
                if (input$probeType == "II"){
                    ylim <- c(0,6)
                } else {
                    ylim <- c(0,10)
                }
                xlim  = c(-0.2,1.2)
            } else {
                xlim  = c(-8,8)
                ylim = c(0,0.35)
            }
            
            range.x <- xlim[2]-xlim[1]
            range.y <- ylim[2]-ylim[1]
            
            if (!is.null(mouse.x)){
                density.matrix <- returnDensityMatrix()
                x.matrix <- density.matrix[[1]]
                y.matrix <- density.matrix[[2]]
                x.matrix <- ((x.matrix - mouse.x)/range.x)^2
                y.matrix <- ((y.matrix - mouse.y)/range.y)^2
                clickedIndex <- arrayInd(which.min(x.matrix+y.matrix),
                                         c(512,ncol(x.matrix)))[2]
                names(clickedIndex) <- shinyMethylSet1@sampleNames[clickedIndex]
                if (current.probe.type == as.character(input$probeType) &
                    current.density.type == as.character(input$mOrBeta)){
                    mouse.click.indices <<- check.mouse.clicks(clickedIndex,
                                                               mouse.click.indices)
                }
                current.probe.type <<- as.character(input$probeType)
                current.density.type <<- as.character(input$mOrBeta)
            }
        }
        mouse.click.indices
	})
        
    ## To update the list of mouse clicks from normalized densities plot:
	updateMouseClicksDensitiesNorm <- reactive({
        mouse.x <- input$normHover$x
        mouse.y <- input$normHover$y
        
        if (!is.null(mouse.x) & !is.null(mouse.y)){
            ## If Beta values are selected
            if (input$mOrBeta=="Beta-value") {
                if (input$probeType == "II"){
                    ylim <- c(0,6)
                } else {
                    ylim <- c(0,10)
                }
                xlim  = c(-0.2,1.2)
            } else {
                xlim  = c(-8,8)
                ylim = c(0,0.35)
            }
            
            range.x <- xlim[2]-xlim[1]
            range.y <- ylim[2]-ylim[1]
            
            if (!is.null(mouse.x)){
                density.matrix <- returnDensityMatrixNorm()
                x.matrix <- density.matrix[[1]]
                y.matrix <- density.matrix[[2]]
                x.matrix <- ((x.matrix - mouse.x)/range.x)^2
                y.matrix <- ((y.matrix - mouse.y)/range.y)^2
                clickedIndex <- arrayInd(which.min(x.matrix+y.matrix),
                                         c(512,ncol(x.matrix)))[2]
                names(clickedIndex) <- shinyMethylSet2@sampleNames[clickedIndex]
                mouse.click.indices <<- check.mouse.clicks(clickedIndex,
                                                           mouse.click.indices)
            }
        }
        mouse.click.indices
	})
        
    ## To update the list from all plots:
	reactive.mouse.click.indices <- reactive({
            updateMouseClicks()
            updateMouseClicksQC()
            updateMouseClicksDensities()
            updateMouseClicksDensitiesNorm()
            return(mouse.click.indices)
	})
        
    ## Print the list of selected samples:
	output$cumulativeListPrint <- renderPrint({
            names <- names(reactive.mouse.click.indices())
            n <- length(names)
            if (n >=1){
                cat("Selected samples: \n \n")
                for (ii in 1:n){
                    cat(paste0(names[ii], "\n"))
                }
            }
	})
	
	
	cumulativeList <- reactive({
            names <- names(reactive.mouse.click.indices())
            return(names)
	})
	
        ## To write the selected samples into a csv file:
	output$selectedSamples <- downloadHandler(
            
            filename <- "selectedSamples.csv",
            content <- function(con){
                write.csv(cumulativeList(),con)
            }
    )
        
        
        
        
#########################################################
###########         Densities plots
#########################################################

        ## Densities plot for raw data  
        output$rawDensities <- renderPlot({
            
            set.palette(n=8, name=setColor())
            colors <- sampleColors()
            lwd <- as.numeric(input$lwd)
            lty <- as.numeric(input$lty)
            index = match(input$probeType,c("I Green","I Red","II","X","Y"))
            
            
            
            ## Plot specifications                  
            if (input$mOrBeta=="Beta-value"){
                xlim <- c(-0.2,1.2)
                if (input$probeType == "II"){
                    ylim <- c(0,6)
                } else {
                    ylim <- c(0,10)
                }
                from = -4; to = 4;
                main = "BETA-VALUE DENSITIES"
                xlab = "Beta-values"
                bw <- input$bandwidth
                quantiles <- betaQuantiles[[index]]
            } else {
                xlim <- c(-8,8)
                ylim <- c(0, 0.35)
                from = -10; to = 10;
                main = "M-VALUE DENSITIES" 
                xlab = "M-values"
                quantiles <- mQuantiles[[index]]
                bw <- input$bandwidth2
            }
            
            
            densitiesPlot(matrix.x = returnDensityMatrix()[[1]],
                          matrix.y = returnDensityMatrix()[[2]],
                          quantiles = quantiles,
                          sampleNames = sampleNames,
                          main = main, xlab = xlab,
                          xlim = xlim, ylim = ylim, col = colors,
                          mean = input$mean, 
                          lwd = lwd, lty = lty,
                          from = from, to = to, bw = bw)
            
            if (!input$mean){
                addHoverDensity(selectedSamples = reactive.mouse.click.indices(), 
                                sampleNames = sampleNames,
                                matrix.x = returnDensityMatrix()[[1]],
                                matrix.y = returnDensityMatrix()[[2]],
                                col = colors)   
            }
            
            
            
            
            ## To draw the lines:              
            if (input$mOrBeta=="Beta-value"){
        	abline(v=0,lty=3,lwd=3)
                abline(v=1,lty=3,lwd=3)
            } else {
        	abline(v=0,lty=3,lwd=3)
            }         
        })
                

        ## Density plot for normalized data:
        output$normDensities <- renderPlot({
            if (is.null(shinyMethylSet2)){
                return(NULL)
            }	else {
                
                set.palette(n=8, name=setColor())
                colors <- sampleColors()
                lwd <- as.numeric(input$lwd)
                lty <- as.numeric(input$lty)
                index <- match(input$probeType,
                               c("I Green","I Red","II","X","Y"))
                
                ## Plot specifications                  
                if (input$mOrBeta=="Beta-value"){
                    xlim <- c(-0.2,1.2)
                    if (input$probeType == "II"){
                        ylim <- c(0,6)
                    } else {
                        ylim <- c(0,10)
                    }
                    from = -4; to = 4;
                    main = "Normalized data (Beta values)"
                    xlab = "Beta-values"
                    bw <- input$bandwidth
                    quantiles =  shinyMethylSet2@betaQuantiles[[index]]
		} else {
                    xlim <- c(-8,8)
                    ylim <- c(0, 0.35)
                    from = -10; to = 10;
                    main = "Normalized data (M values)" 
                    xlab = "M-values"
                    quantiles =  shinyMethylSet2@mQuantiles[[index]]
                    bw <- input$bandwidth2
		}
                
		densitiesPlot(matrix.x = returnDensityMatrixNorm()[[1]],
                              matrix.y = returnDensityMatrixNorm()[[2]],
                              quantiles = quantiles,
                              sampleNames = sampleNames,
                              main = main, xlab = xlab,
                              xlim = xlim, ylim = ylim, col = colors,
                              mean = input$mean, 
                              lwd = lwd, lty = lty,
                              from = from, to = to, bw=bw)
                
                if (!input$mean){
                    addHoverDensity(selectedSamples = reactive.mouse.click.indices(), 
                                    sampleNames = sampleNames,
                                    matrix.x = returnDensityMatrixNorm()[[1]],
                                    matrix.y = returnDensityMatrixNorm()[[2]],
                                    col = colors)  
                }
                
                ## To draw the lines:              
                if (input$mOrBeta=="Beta-value"){
                    abline(v=0,lty=3,lwd=3)
                    abline(v=1,lty=3,lwd=3)
                } else {
                    abline(v=0,lty=3,lwd=3)
                }  
            }
        })
        


#########################################################
###########         Control probes and QC plot
#########################################################

        ## Return a summary measure of the selected control probes:
        reactiveControlStat <- reactive({
            controlStat <- returnControlStat(input$controlType, 
                                             greenControls = greenControls,
                                             redControls = redControls, 
                                             controlNames = controlNames
                                             )
            return(controlStat)
        })
        
        ## Internal controls plots
        output$internalControls <- renderPlot({
            colors <- sampleColors()
            set.palette(n=8, name=setColor())
            controlStat <- reactiveControlStat()
            plotInternalControls(y = controlStat,
                                 col = colors,
                                 main = input$controlType,
                                 sampleNames = sampleNames)
            addHoverPoints(y = controlStat, 
                           selectedSamples = reactive.mouse.click.indices(), 
                           sampleNames = sampleNames)
        })

        ## Quality control plot
        output$medianChannels <- renderPlot({
            
            set.palette(n=8, name=setColor())
            colors <- sampleColors()
            plotQC(unmethQuantiles = unmethQuantiles, 
                   methQuantiles = methQuantiles,
                   sampleNames = sampleNames,
                   col = colors)
            
            controlStat <- reactiveControlStat()
            addHoverQC(y = controlStat, 
                       selectedSamples = reactive.mouse.click.indices(), 
                       unmethQuantiles = unmethQuantiles,
                       methQuantiles = methQuantiles)
        })
        
        
        
#########################################################
###########         Sex plot
#########################################################


	## To change the gender cutoff for prediction:
	setGenderCutoff <- reactive({
            if (!is.null(input$genderCutoff)){
                genderCutoff <<- input$genderCutoff$x
            }
	})
	output$diffPrint <- renderPrint({
            if (ncol(data())!=3){
                cat("No gender was provided in the phenotype data")
            } else {
                diff.genders <- diffGenders()
                diff.genders <- diff.genders[!is.na(diff.genders)]
                if (!is.null(diffGenders())){
                    n <- length(diffGenders())
                    for (i in 1:n){
                        cat(paste0(diffGenders()[i]),"\n")
                    }
                } else {
                    cat("The provided gender agrees with the predicted gender for all samples")
                }
            }
	})
        
        
        
        
        
        ## Sex plot
        output$genderClustering <- renderPlot({
            
            setGenderCutoff()
            plotPredictedGender(cutoff = genderCutoff,
                                cnQuantiles = cnQuantiles)
            plotDiscrepancyGenders(cutoff = genderCutoff,
                                   cnQuantiles = cnQuantiles,
                                   covariates = covariates)
        })



        data <- reactive({
            setGenderCutoff()
            predictedGender <- returnPredictedGender(cutoff = genderCutoff,
                                                     cnQuantiles = cnQuantiles)
            givenGender     <- returnGivenGender(covariates)
            dataToReturn <- data.frame(predictedGender = predictedGender)
            if (!is.null(givenGender)){
                non.matching.samples    <- predictedGender  != givenGender
                na.samples <- is.na(givenGender)
                n <- length(predictedGender)
                agree <- rep("YES",n)
                agree[non.matching.samples] <- "NO"
                agree[na.samples] <- NA
                dataToReturn$givenGender <- givenGender
                dataToReturn$agree <- agree
            }
            rownames(dataToReturn) <- sampleNames
            return(dataToReturn)
        })
        
        
        output$downloadClusters <- downloadHandler(
            filename <- "predictedGender.csv",
            content <- function(con){
                write.csv(data(),con)
        })
        
        
        diffGenders <- reactive({
            
            non.matching.samples <- c()
            diffs <- data()
            if (ncol(diffs)==3){
 		non.matching.samples <- shinyMethylSet1@sampleNames[diffs$agree=="NO"]
 		non.matching.samples <- non.matching.samples[complete.cases(non.matching.samples)]
            }
            return(non.matching.samples)
            
        })
        output$diffPrint <- renderPrint({
            if (ncol(data())!=3){
		cat("No gender was provided in the phenotype data")
            } else {
		diff.genders <- diffGenders()
		diff.genders <- diff.genders[!is.na(diff.genders)]
		if (!is.null(diffGenders())){
                    n <- length(diffGenders())
                    for (i in 1:n){
                        cat(paste0(diffGenders()[i]),"\n")
                    }
		} else {
                    cat("The provided gender agrees with the predicted gender for all samples")
		}
            }
        })
        
        
        
        
        
        
        ## Densities plot for gender X
        output$densitiesGenderX <- renderPlot({
            setGenderCutoff()
            lwd <- as.numeric(input$lwd)
            lty <- as.numeric(input$lty)
            bw <- input$bandwidth
            index <- match("X",c("I Green","I Red","II","X","Y"))
            matrix <- apply(betaQuantiles[[index]], 2, function(x){
                d <- density(x, bw = bw, n =512)
                c(d$x,d$y)
            })
            matrix.x <- matrix[1:512,]
            matrix.y <- matrix[513:(2*512),]
            
            predictedGender <- returnPredictedGender(cutoff = genderCutoff,
                                                     cnQuantiles = cnQuantiles)
            colors <- rep("lightskyblue", length(predictedGender))
            colors[predictedGender=="F"] <- "orange"
            
            
            densitiesPlot(matrix.x = matrix.x,
                          matrix.y = matrix.y,
                          quantiles = betaQuantiles[[index]],
                          sampleNames = sampleNames,
                          main = "X Chromosome", 
                          xlab = "Beta-values",
                          xlim  = c(-0.2,1.2),
                          ylim = c(0,7),
                          bw  = input$bandwidth,
                          lwd = lwd,
                          lty = lty, 
                          mean = FALSE,
                          col = colors, 
                          from=-4, 
                          to =4)
            
            abline(v=0,lty=3,lwd=3)
	    abline(v=1,lty=3,lwd=3)
            
        })
        
        
#########################################################

###########        PCA plot

#########################################################

        output$pcaPlot <- renderPlot({
            set.palette(n=8, name=setColor())
            ##palette(brewer.pal(8,setColor()))
            if (method!="Merging"){
                plotPCA(pca = pca, 
                        pc1 = input$pc1,
                        pc2 = input$pc2,
                        col = sampleColors(),
                        covariates = covariates,
                        selectedCov = input$phenotype
                        )
            }
        })
        ## Print the summary of the PCA regression: 
        output$modelPrint <- renderPrint({
            if (method!="Merging"){
		y <- shinyMethylSet1@pca$scores[,as.numeric(input$pcToExplore)]
                cov <- covariates[match(rownames(shinyMethylSet1@pca$scores),
                                        rownames(covariates)),]
		x <- (as.factor(cov[,match(input$covToRegress,colnames(cov))]))
                model <- lm(y~ x)
                return(summary(model))
            }
            else {
                return("Merging was used to create the shinyMethylSet1; no PCA was performed. ")
            }
        })
        
        
        
        
#########################################################
###########        Design plot
#########################################################

        output$arrayDesign <- renderPlot({
            set.palette(n=8, name=setColor())
            color <- covariates[,match(input$phenotype,colnames(covariates))]
            plotDesign450k(as.character(sampleNames), covariates = color ,
                           legendTitle = input$phenotype)
        })
        output$arrayDesignLegend <- renderPlot({
            set.palette(n=8, name=setColor())
            color <- covariates[,match(input$phenotype,colnames(covariates))]
            plotLegendDesign450k(as.character(sampleNames), covariates =color ,
                                 legendTitle = input$phenotype)
        })

        


#########################################################
###########       Probe Bias Plot
#########################################################



## Densities plot for raw data  
        output$probeBiasPlot <- renderPlot({
            set.palette(n=8, name=setColor())
            colors <- sampleColors()
            lwd <- as.numeric(input$lwd)
            lty <- as.numeric(input$lty)
            index = match(input$probeType,c("I Green","I Red","II","X","Y"))
            selectedSample <- as.character(input$selectedSampleBias)
            sampleIndex <- match(selectedSample, sampleNames)
            indexIGreen <- 1
            indexIRed   <- 2 
            indexII     <- 3
            
            bw <- input$bandwidth 
            quantilesIGreen <- betaQuantiles[[indexIGreen]][,sampleIndex]
            quantilesIRed <- betaQuantiles[[indexIRed]][,sampleIndex]
            quantilesII <- betaQuantiles[[indexII]][,sampleIndex]
            quantilesIRed <- as.vector(quantilesIRed)
            quantilesIGreen <- as.vector(quantilesIGreen)
            quantilesII <- as.vector(quantilesII)
            
            
            ## Plot specifications                  
            xlim <- c(-0.2,1.2)
            if (input$probeType == "II"){
   		ylim <- c(0,6)
            } else {
   		ylim <- c(0,10)
            }
            from = -4; to = 4;
            main = "Probe Type Differences - Raw Data"
            xlab = "Beta-values"
            bw <- input$bandwidth
            
            plot(density(quantilesIGreen, bw=bw), 
                 main = main, 
                 ylab = "Density", 
                 xlab = xlab,
                 col = "olivedrab",
                 xlim = xlim, 
                 ylim = ylim, lty=lty, lwd=6)   
            lines(density(quantilesII, bw=bw), col="black", lty=lty, lwd=6)
            lines(density(quantilesIRed, bw=bw), col="brown2", lty=lty, lwd=6)
            
            
            ## To draw the lines:              
            abline(v=0,lty=3,lwd=3)
            abline(v=1,lty=3,lwd=3)
            
            legend(x=0.4, y=8, 
                   c("Type I Red","Type I Green","Type II"), 
                   col=c("brown2","olivedrab","black"), lwd=6, lty=1, bty="n")
            
        })
        
        ## Densities plot for raw data  
        output$probeBiasPlotNorm <- renderPlot({
            
            if (is.null(shinyMethylSet2)){
                return(NULL)
            }	else {
                
                set.palette(n=8, name=setColor())
                colors <- sampleColors()
                lwd <- as.numeric(input$lwd)
                lty <- as.numeric(input$lty)
                index = match(input$probeType,c("I Green","I Red","II","X","Y"))
                
                selectedSample <- as.character(input$selectedSampleBias)
                sampleIndex <- match(selectedSample, sampleNames)
                
                indexIGreen <- 1
                indexIRed   <- 2 
                indexII     <- 3
                
                bw = input$bandwidth 
                quantilesIGreen <- shinyMethylSet2@betaQuantiles[[indexIGreen]][,sampleIndex]
                quantilesIRed <- shinyMethylSet2@betaQuantiles[[indexIRed]][,sampleIndex]
                quantilesII <- shinyMethylSet2@betaQuantiles[[indexII]][,sampleIndex]
                quantilesIRed <- as.vector(quantilesIRed)
                quantilesIGreen <- as.vector(quantilesIGreen)
                quantilesII <- as.vector(quantilesII)
                
                
                ## Plot specifications                  
                xlim <- c(-0.2,1.2)
                if (input$probeType == "II"){
                    ylim <- c(0,6)
                } else {
                    ylim <- c(0,10)
                }
                from = -4; to = 4;
                main = "Probe Type Differences - Normalized Data"
                xlab = "Beta-values"
                bw <- input$bandwidth
                
                plot(density(quantilesIGreen, bw=bw), 
                     main = main, 
                     ylab = "Density", 
                     xlab = xlab,
                     col = "olivedrab",
                     xlim = xlim, 
                     ylim = ylim, lty=lty, lwd=6)   
                lines(density(quantilesII, bw=bw), col="black", lty=lty, lwd=6)
                lines(density(quantilesIRed, bw=bw), col="brown2", lty=lty, lwd=6)
                
                legend(x=0.4, y=8, 
                       c("Type I Red","Type I Green","Type II"), 
                       col=c("brown2","olivedrab","black"), lwd=6, lty=1, bty="n")
                
                                        # To draw the lines:              
                abline(v=0,lty=3,lwd=3)
                abline(v=1,lty=3,lwd=3)
                
            }
            
        })
        
        
        ## End 
    }    
}
