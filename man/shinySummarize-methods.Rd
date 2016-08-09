\name{shinySummarize-methods}
\docType{methods}
\alias{shinySummarize}
\alias{shinySummarize-methods}
\alias{shinySummarize,RGChannelSet-method}
\alias{shinySummarize,GenomicRatioSet-method}
\title{Summarizing a large 450K experiment into a shinyMethylSet}
\description{
 Summarizing methylation data from a RGChannelSet or GenomicRatioSet into a shinyMethylSet needed to launch the interactive interface of shinyMethyl. That works for both the 450k and EPIC arrays. 
}

\usage{
\S4method{shinySummarize}{RGChannelSet}(object)
\S4method{shinySummarize}{GenomicRatioSet}(object)
}

\arguments{
  \item{object}{Either a \code{RGChannelSet} or a
    \code{GenomicRatioSet}.} 
}


\details{
  Either a \code{RGChannelSet} is transformed into a \code{shinyMethylSet}, or a \code{GenomicRatioSet} is transformed to a \code{shinyMethylSet}
}

\author{Jean-Philippe Fortin \email{jfortin@jhsph.edu}}

\seealso{
  \code{\link{shinyMethylSet}} for the output object and
  \code{\link{RGChannelSet}} or \code{\link{GenomicRatioSet}} for the input object. 
}
\examples{
if (require(minfiData)){
  library(minfiData)
  summarized.data <- shinySummarize(RGsetEx)
}
}



\keyword{methods}

