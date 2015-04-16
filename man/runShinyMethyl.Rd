\name{runShinyMethyl}
\alias{runShinyMethyl}
\title{
Run the interactive shinyMethyl session
}
\description{
  Function to run the interactive shinyMethyl session from a shinyMethylSet. 
}
\usage{
runShinyMethyl(shinyMethylSet1, shinyMethylSet2 = NULL)
}
\arguments{
  \item{shinyMethylSet1}{
    \code{shinyMethylSet} that must be extracted from an \code{RGChannelSet}
  }
  \item{shinyMethylSet2}{
    optional \code{shinyMethylSet} that must be extracted from an
    \code{GenomicRatioSet} 
  }
}

\value{
No value returned. Instead the shinyMethyl interactive session is launched.
}
\author{
Jean-Philippe Fortin
}
\seealso{
\code{\link{shinyMethylSet}}
}
\examples{
\dontrun{
if (require(minfiData)){
  baseDir <- system.file("extdata", package = "minfiData")
  targets <- read.450k.sheet(baseDir)
  targets$Sample_Plate <- substr(targets$Slide,1,7)
  RGSet <- read.450k.exp(base = baseDir, targets = targets)

  summarized.data <- shinySummarize(RGSet)
  runShinyMethyl(summarized.data)
}
}
}
