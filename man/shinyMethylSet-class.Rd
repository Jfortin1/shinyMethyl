\name{shinyMethylSet-class}
\Rdversion{1.1}
\docType{class}
\alias{shinyMethylSet-class}
\alias{shinyMethylSet}
%%generics
\alias{getMeth}
\alias{getUnmeth}
\alias{getBeta}
\alias{getM}
\alias{getCN}
\alias{getGreenControls}
\alias{getRedControls}
\alias{getPCA}
\alias{pData}
\alias{sampleNames}
\alias{orderByName}
\alias{combine}
%%methods
\alias{combine,shinyMethylSet,shinyMethylSet-method}
\alias{getMeth,shinyMethylSet-method}
\alias{getUnmeth,shinyMethylSet-method}
\alias{getBeta,shinyMethylSet-method}
\alias{getM,shinyMethylSet-method}
\alias{getCN,shinyMethylSet-method}
\alias{show,shinyMethylSet-method}
\alias{pData,shinyMethylSet-method}
\alias{sampleNames,shinyMethylSet-method}

\title{shinyMethylSet instances}

\description{
  This class holds summarized data from Illumina methylation
  microarrays for interactive visualization purpose.
}

\usage{
## Constructor


shinyMethylSet(sampleNames = new("character"), 
               phenotype = new("data.frame"),  
               mQuantiles = new(vector("list",5)),
               betaQuantiles = new(vector("list",5)),
               methQuantiles = new(vector("list",5)),
               unmethQuantiles = new(vector("list",5)),
               cnQuantiles = new(vector("list",5)),
               greenControls = new(vector("list",12)),
               redControls = new(vector("list",12)),
               pca = new("list"),
               originObject = new("character"),
               array = new("character"))

## Data extraction / Accessors


\S4method{getMeth}{shinyMethylSet}(object)
\S4method{getUnmeth}{shinyMethylSet}(object)
\S4method{getBeta}{shinyMethylSet}(object)
\S4method{getM}{shinyMethylSet}(object)
\S4method{getCN}{shinyMethylSet}(object)
\S4method{pData}{shinyMethylSet}(object)
\S4method{sampleNames}{shinyMethylSet}(object)
}

\arguments{
  \item{object}{A \code{shinyMethylSet}.}
  \item{sampleNames}{A character vector}
  \item{phenotype}{A \code{data.frame} object.}
  \item{methQuantiles}{A list of 5 matrices containing several quantiles
    for the methylation values (between zero and infinity) separated by
    probe tpye and autosomal/sex probes. Each row is a quantile and each
    column is a sample.}
  \item{unmethQuantiles}{Similar to \code{Meth}}
  \item{betaQuantiles}{Similar to \code{Meth} but for Beta-values}
  \item{mQuantiles}{Similar to \code{Meth} but for M-values}
  \item{cnQuantiles}{Similar to \code{Meth} but for Copy Number}
  \item{greenControls}{A list containing the matrices of different raw
    control probes intensities in the green channel} 
  \item{redControls}{Similar to \code{greenControls} but for the red channel}
  \item{pca}{List containing the PCA scores for the 20,000 most variable
    CpGs and the percentages of variance explained} 
  \item{originObject}{Name of the minfi object from which the data were
    extracted} 
  \item{array}{Name of the Illumina array platform: 450k or EPIC} 
}

\details{
  These class is a representation of a \code{Meth} matrix and a
  \code{Unmeth} matrix linked to a \code{pData} data frame.

  A creationMethod slot is present to indicate from which object type
  the \code{shinyMethylSet} has been created: either a
  \code{RGChannelSet} or \code{GenomicRatioSet}
  
  A \code{MethylSet} stores \code{meth} and \code{Unmeth}. }


\section{Constructor}{
Instances are constructed using the \code{shinyMethylSet} function with the
arguments outlined above.
}

\section{Accessors}{  
  In the following code, \code{object} is a \code{shinyMethylSet}.
  
  \describe{
    \item{\code{getMeth(object)}, \code{getUnmeth(object)}}{Get the
      Meth or Unmeth matrix.}
    \item{\code{getBeta(object)}}{get Beta, see details.}
    \item{\code{getM(object)}}{get M-values, see details.}
    \item{\code{getCN(object)}}{get copy number values which are defined
      as the sum of the methylation and unmethylation channel.}
    \item{\code{sampleNames(object)}}{Get
      the sampleNames}
    \item{\code{combine(object1,object2)}}{combine \code{shinyMethylSet} objects}
  }
}

\author{Jean-Philippe Fortin \email{jfortin@jhsph.edu}}

\seealso{
  Objects of this class are typically created by using the function
  \code{\link{shinySummarize}} on a \code{\link{RGChannelSet}} or
  \code{\link{GenomicRatioSet}}.}

\examples{
  showClass("shinyMethylSet")
}
