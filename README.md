shinyMethyl
===========

Authors: [Jean-Philippe Fortin](mailto:zerbino@ebi.ac.uk) and [Kasper Daniel Hansen](mailto:khansen@jhsph.edu)

Welcome to shinyMethyl, an interactive R application based on the shiny package for exploration of DNA methylation data from Illumina 450K arrays. Please find an online demo of shinyMethyl [here](http://spark.rstudio.com/jfortin/shinyMethyl/)

Installation
------------
First, you will need to install at least the following packages from Bioconductor

```{r}
source("http://bioconductor.org/biocLite.R")
biocLite("minfi")
biocLite("minfiData")
```
and from CRAN
```{r}
install.packages("devtools")
install.packages("matrixStats")
install.packages("RColorBrewer")
```
To install the development version of shinyMethyl:
```{r}
library(devtools)
install_github("shiny", "rstudio")
install_github("shinyMethyl", "jfortin1", quick=TRUE)
install_github("shinyMethylData", "jfortin1")
```

Vignette
------------
You can find the vignette for shinyMethyl at [https://github.com/Jfortin1/shinyMethyl/blob/master/vignettes/shinyMethyl.pdf](https://github.com/Jfortin1/shinyMethyl/blob/master/vignettes/shinyMethyl.pdf)

Quick example
------------
After installation, you can launch shinyMethyl with an example dataset from TCGA with the following code:
```{r}
library(shinyMethyl)
library(shinyMethylData)
runShinyMethyl(tcga.summary.raw, tcga.summary.norm)
```

# Citation


To cite package __shinyMethyl__ in publications use:

Jean-Philippe Fortin and Kasper Daniel Hansen (2014). shinyMethyl: 
Interactive visualization of 450k methylation data. R package 
https://github.com/Jfortin1/shinyMethyl

A BibTeX entry for LaTeX users is

@Manual{,
    title = {shinyMethyl: 
Interactive visualization of 450k methylation data},
    author = {Jean-Philippe Fortin and Kasper Daniel Hansen},
    year = {2014},
    note = {R package version},
    url = {https://github.com/Jfortin1/shinyMethyl},
}

Updates
------------

- The package is now build using S4 classes; the current functions of shinyMethyl are no longer in used. The function shinySummarize() applied to a RGChannelSet object is now used to produce a shinyMethylSet, which is passed to runShinyMethyl() to launch a shinyMethyl session. Please see the vignette.
- Fixed returnPCScores() dependency on annotation package. Thanks to Maarten van Iterson 
- extractFromTargets450k() is temporarily non-available. Please use extractFromRGSet450k() instead. 
- Bug fixed: M-value densities are now visible
- With the update of minfi in Bioconductor 2.13, IlluminaHumanMethylation450kannotation.ilmn.v1.2 is no longer supported. ShinyMethyl is now updated to be annotation free (shinyMethyl does not really need an annotation). Thanks to Kathleen Fisch for pointing that out.
- Within the predictGender panel, it is now possible to download a .csv file containing not only the predicted gender, but also the actual gender if it was provided in the phenotype data, it will be also included in the file. The function also handles missing data. Thanks to Brent Pedersen


