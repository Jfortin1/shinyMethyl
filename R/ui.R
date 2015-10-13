#### Created by Jean-Philippe Fortin
#### March 28, 2014

ui.shinyMethyl <- function(shinyMethylSet1, shinyMethylSet2=NULL){
    betaQuantiles   <- getBeta(shinyMethylSet1)
    mQuantiles      <- getM(shinyMethylSet1)
    methQuantiles   <- getMeth(shinyMethylSet1)
    unmethQuantiles <- getUnmeth(shinyMethylSet1)
    cnQuantiles     <- getCN(shinyMethylSet1)
    greenControls   <- getGreenControls(shinyMethylSet1)
    redControls     <- getRedControls(shinyMethylSet1)
    covariates      <<-pData(shinyMethylSet1)
    pca             <- getPCA(shinyMethylSet1)$scores
    sampleNames     <- sampleNames(shinyMethylSet1)
    slideNames      <- substr(sampleNames,1,10)
    arrayNames      <- substr(sampleNames,12,17)
    plateNames      <- substr(sampleNames,1,6)
    controlNames    <- names(greenControls)
    slides <- unique(slideNames)
    method <- shinyMethylSet1@originObject
    sampleColors <<- as.numeric(as.factor(plateNames))
    
    ## In the case covariates is empty:
    if (ncol(covariates)==0){
	covariates <- data.frame(slide = slideNames, plate = plateNames)
	rownames(covariates) <- sampleNames
	covariates <<- covariates
    }
    
    ## shinyUI(pageWithSidebar(
    pageWithSidebar(
###########################  ---  Header ------------ ##############
        #headerPanel(
        #    HTML("<p style=\"color:#000000;font-family:\"Times New Roman\",Georgia,Serif\">
        # shiny<span style=\"color:#E56717\">M</span><span style=\"color:#000000\">ethyl</span></p>")
        #    ),
        headerPanel("shinyMethyl"),
###########################  ---  Sidebar ---------- ###############
        sidebarPanel(
            wellPanel(
		HTML("<p><span style=\"color:#336666;font-size:16px\">
			      Color choice:</span></p>"),
                
                selectInput("colorChoice", "Color set:",
                            list("Pault","Rainbow","Set1","Set2","Set3",
                                 "Paired","Dark2","Accent"),
                            multiple=FALSE, 
                            selected="Set1"),
                HTML("<p><span style=\"color:#336666;font-size:16px\">
			      Quality control exploration:</span></p>"),
                
                selectInput("mOrBeta", "Methylation measure:",
                            list("Beta-value","M-value"),
                            multiple=FALSE, 
                            selected="M-value"),
                selectInput("probeType", "Choose a probe type for the density curves:", 
                            choices = c("I Green","I Red","II","X","Y"),
                            selected="I Green"),
                selectInput("controlType", "Choose a control type:", 
                            choices = controlNames,selected=controlNames[1]),
                if ("plate" %in% colnames(covariates)){
                    choices <- colnames(covariates)
                    selectInput("phenotype", "Choose a phenotype:",
                                choices,
                                multiple=FALSE, 
                                selected ="plate")
		} else {
                    choices <- colnames(covariates)
                    selectInput("phenotype", "Choose a phenotype:",
                                choices,
                                multiple=FALSE)},          	
                
                checkboxInput("mean","Average density by phenotypic level"),
                selectInput("lty", "Density line type (lty):",
                            list(1,2,3,4,5,6),
                            multiple=FALSE, 
                            selected=1),
                selectInput("lwd", "Density line width (lwd):",
                            list(0.5,1,1.5,2,3,4,5,6,7,8,9,10),
                            multiple=FALSE, 
                            selected=1)	
                ),
##########  -- Sliders       ---######       
            wellPanel(
                sliderInput(inputId = "bandwidth",
                            label = "Bandwidth Beta-value",
                            min = 0, max = 0.05, step = 0.001, value = 0.02),
                sliderInput(inputId = "bandwidth2",
                            label = "Bandwidth M-value",
                            min = 0, max = 0.5, step = 0.001, value = 0.35)
                )
            ),
        mainPanel(
            tabsetPanel(
###########################  ---  Home   ------------ 
                tabPanel("Home",
                         HTML("<br>
  		<p style=\"width:500px;text-align:justify\"><span style=\"color:#000000;font-size:16px\">
  		Welcome to
  		 <span style=\"font-weight:bold\">shiny</span><span style=\"color:#E56717;font-weight:bold\">M</span><span style=\"font-weight:bold\">ethyl</span>, 
  		 an interactive visualization R package for 
  		exploration and quality control of methylation data. The current version
  		is designed for Illumina Human Methylation 450K arrays. 
 <br><br>For more information, please visit <span style=\"font-style:italic\">shinyMethyl.com</span><br><br><br></p>")),

###########################  ---  Quality control --- 
                tabPanel("Quality control",
                         ## Densities plot:
                         HTML('<table border=0 width="100%"><tr bgcolor="#f5f5f5"><td>'),
                         div(style="width:100%;max-width:600px;",
                             ## Internal controls:
                             plotOutput("internalControls", clickId = "controlsHover")), 
                         HTML('</td><td>'),
                         ## Fast quality control plot:
                         plotOutput("medianChannels", clickId = "qualityHover"),
                         HTML('</td></tr></table>'),
                         plotOutput("rawDensities",clickId = "densitiesHover"),
                         conditionalPanel(condition= "!is.null(shinyMethylSet2)",
                                          plotOutput("normDensities",
                                                     clickId = "normHover")),
                         verbatimTextOutput(outputId = "cumulativeListPrint"),
                         downloadLink("selectedSamples","selectedSamples.csv")
                         ),
######################   ----   Array Design  -------- 
                tabPanel("Array Design",
                         HTML("<br>
  		<p style=\"width:800px;text-align:justify\"><span style=\"color:#000000;font-size:16px\">
  		The Illumina 450k samples are divided into slides. A slide contains 12 samples (6 by 2 grid) an a plate contains 8 slides (96 samples). The plot below shows the allocation of the samples to the plates and the coloring allows the user to judge if the design is well-balanced for different phenotype covariates.
  		</span></p><br><br>") ,
                         wellPanel(
                             div(style="max-height:800px;",
                                 plotOutput("arrayDesign")
                                 ),
                             div(style="max-height:800px;",
                                 plotOutput("arrayDesignLegend")
                                 )
                             )
                         ),
######################   ----   Gender clustering  --------  
                tabPanel("Gender clustering",
                         HTML("<br>
  		<p style=\"width:800px;text-align:justify\"><span style=\"color:#000000;font-size:16px\">
  		By comparing the median total intensity of the Y-chromosome-mapped probes to the median total intensity of the X-chromosome-mapped probes, where the total intensity is the sum of the methylated and unmethylated signals, it is possible to predict the gender of the sample by looking at the two distinct clusters of intensities. See the minfi function <span style=\"font-style:italic\">getSex()</span>.  		</span></p><br><br>") ,
                         HTML('<table border=0 width="100%"><tr bgcolor="#f5f5f5"><td>'),
                         div(style="width:100%;max-width:600px;",
                             plotOutput("genderClustering", clickId = "genderCutoff")), 
                         HTML('</td><td>'),
                         plotOutput("densitiesGenderX"),
                         HTML('</td></tr></table>'),
                         sidebarPanel(
                             HTML("
<p style=\"color:#000000;font-size:17px\">A. Click on the plot to choose the cutoff</span></p>
"),
                             HTML("
<p style=\"color:#000000;font-size:17px\">B. Save the predicted gender in a csv file:</span></p>
"),   
                             downloadLink("downloadClusters","predictedGender.csv")
                             ),
                         sidebarPanel(
                             HTML("
<p style=\"color:#000000;font-size:17px\">List of samples whose predicted gender do not agree on the gender provided in the phenotype data:</span></p>
"),
                             verbatimTextOutput(outputId = "diffPrint")
                             )
                         ),
######################   ----   PCA --------  
                tabPanel("PCA",
                         plotOutput("pcaPlot"),
                         HTML("
<p style=\"color:#000000;font-size:17px\">A. Choose two principal components to visualize: </span></p>
"),
                         selectInput("pc1", "PC in x:",
                                     seq(1,ncol(betaQuantiles[[1]]),1),
                                     multiple=FALSE, 
                                     selected=1),
                         selectInput("pc2", "PC in y:",
                                     seq(1,ncol(betaQuantiles[[1]]),1),
                                     multiple=FALSE, 
                                     selected=2),
                         HTML("
<p style=\"color:#000000;font-size:17px\">B. Choose a principal component to explore below:</span></p>
"),   
                         selectInput("pcToExplore", "PC:",
                                     seq(1,ncol(betaQuantiles[[1]]),1),
                                     multiple=FALSE, 
                                     selected=1),
                         HTML("
<p style=\"color:#000000;font-size:17px\">C. Choose a covariate to regress against the chosen PC:</span></p>
"),   
                         if (exists("covariates")){
                             choices <- colnames(covariates)
                             selectInput("covToRegress", "Covariate:",
                                         choices,
                                         multiple= FALSE)
                         } else {
                             selectInput("covToRegress", "Covariate:",
                                         list("Batch"),
                                         multiple= FALSE, 
                                         selected="Batch")},
                         HTML("
<p style=\"color:#000000;font-size:17px\">Association with the PC:</span></p>
"), 
                         verbatimTextOutput(outputId = "modelPrint")
                         ),
######################   ----   Type I/TypeII Bias --------  
		tabPanel("Probe Type Diff",
                         selectInput("selectedSampleBias", "Sample:",
                                     sampleNames,
                                     multiple= FALSE),
                         plotOutput("probeBiasPlot"),
                         conditionalPanel(condition= "!is.null(shinyMethylSet2)", plotOutput("probeBiasPlotNorm"))
                         ),
######################   ----   Reproducible report  --------  
                ## tabPanel("Reproducible Report",
                ##      HTML("<br>
                ## <p style=\"width:500px;text-align:justify\"><span style=\"color:#000000;font-size:16px\">
                                        #   <span style = \"font-weight:bold\">Reproducible Report</span><br><br>
                                        #   <span style=\"font-style:italic\">shinyMethyl</span> visualizations are reproducible in the sense that every plot can be reproduced with R functions included in the package, with the tuning parameters specified by the user.
                ## <br><br>
                ## </p>"),  
                
                ## numericInput("gender.cutoff.markdown", "Enter the prediction threshold used in the Gender clustering panel:", -0.3, min=-5, max=5, step=0.1),
                
                ## actionButton('create.report',"Create HTML Report"),
		
                ##   verbatimTextOutput(outputId = "reportPrint")
		
		
                
                                        ##  ),
                
######################   ----   About  --------  
                
############################################################
                
		tabPanel("About",
                         HTML("<br>
  		<p style=\"width:500px;text-align:justify\"><span style=\"color:#000000;font-size:16px\">
  		<span style = \"font-weight:bold\">About</span><br><br>
  		<span style=\"font-style:italic\">shinyMethyl</span> is a complementary tool to <span style=\"font-style:italic\">minfi</span> for visualizing methylation data from Illumina 450K arrays. 
  		<br><br>
          <span style = \"font-weight:bold\">Acknowledgements</span><br><br>
           The <span style=\"font-style:italic\">shinyMethyl</span> application is based on the package <span style=\"font-style:italic\">minfi</span>, and a large part of the source code is inspired by the work done by the <span style=\"font-style:italic\">minfi</span>'s authors. The gender clustering panel and the lower QC plot of the quality control panel are based on algorithms available in <span style=\"font-style:italic\">minfi</span> at the Bioconductor project.
          <br><br> 
          <span style=\"font-style:italic\">shinyMethyl</span> is currently developed at the Johns Hopkins Department of Biostatistics. Many thanks to Elizabeth M. Sweeney, John Muschelli and Leonardo Collado Torres for their help and for providing precious feedbacks. 
          </br></br>
          <span style = \"font-weight:bold\">Authors:</span>  Jean-Philippe Fortin and Kasper Daniel Hansen
         <br><br>
  		</span></p>")
                         )
                )
            )
        )
}
