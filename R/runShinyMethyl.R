#' @title Run the interactive shinyMethyl session
#' 
#' @description Function to run the interactive shinyMethyl session
#'     from a shinyMethylSet object. 
#' 
#' @param shinyMethylSet1 \code{shinyMethylSet} that must be extracted from
#'     an \code{RGChannelSet} object. 
#' @param shinyMethylSet2 Optional \code{shinyMethylSet} that must be
#'     extracted from a \code{GenomicRatioSet}. 
#' 
#' @return No value returned. Instead the shinyMethyl interactive
#'     session is launched.
#' 
#' @author Jean-Philippe Fortin
#' 
#' @seealso \code{\link{shinyMethylSet}}
#' 
#' @examples 
#' \dontrun{
#' if (interactive()){
#' library(minfi)
#' baseDir <- system.file("extdata", package = "minfiData")
#' targets <- read.metharray.sheet(baseDir)
#' targets$Sample_Plate <- substr(targets$Slide,1,7)
#' RGSet <- read.metharray.exp(targets=targets)
#' summarized.data <- shinySummarize(RGSet)
#' runShinyMethyl(summarized.data)
#' }
#' }
#' @export
runShinyMethyl <- function(shinyMethylSet1,
                           shinyMethylSet2=NULL
){
    directory <- system.file(package="shinyMethyl", "shinyMethyl")
    if (shinyMethylSet1@originObject != "RGChannelSet") {
        stop("First shinyMethylSet object must be created from an RGChannelSet")
    }

    shinyMethylSet1 <- shinyMethyl::orderByName(shinyMethylSet1)
    
    # If a second shinyMethylSet is provided (from GenomicRatioSet):
    if (!is.null(shinyMethylSet2)) {
        if (shinyMethylSet2@originObject != "GenomicRatioSet") {
            stop("Second shinyMethylSet object must be created from a GenomicRatioSet")
        }
        ## Need to make sure both shinyMethylSet's are compatible
        if (all.equal(shinyMethylSet1@sampleNames,
                      shinyMethylSet2@sampleNames)!=TRUE){
            stop("The two shinyMethylSet objects are not compatible.",
                 " Both shinyMethylSet objects must be created from ",
                 "the same samples.")
        }
        shinyMethylSet2 <- shinyMethyl::orderByName(shinyMethylSet2)
    } 

    ui <- ui.shinyMethyl(shinyMethylSet1, shinyMethylSet2)
    server <- server.shinyMethyl(shinyMethylSet1, shinyMethylSet2)
    shinyMethyl.app <- list(ui=ui, server=server)
    runApp(shinyMethyl.app)
}

