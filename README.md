shinyMethyl
===========

[![Build Status](https://travis-ci.org/Jfortin1/shinyMethyl.svg?branch=master)](https://travis-ci.org/Jfortin1/shinyMethyl)

Authors: [Jean-Philippe Fortin](mailto:fortin946@gmail.com) and [Kasper Daniel Hansen](mailto:khansen@jhsph.edu)

Welcome to `shinyMethyl`, an interactive R application based on the `shiny` package for exploration of DNA methylation data from Illumina arrays (450k and EPIC arrays). `shinyMethyl` is part of the [Bioconductor project](http://www.bioconductor.org/packages/devel/bioc/html/shinyMethyl.html).

The `shinyMethyl` paper can be found [here](http://f1000research.com/articles/3-175/v2)


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
install.packages("httpuv")
install.packages("devtools")
install.packages("matrixStats")
install.packages("RColorBrewer")
```
To install the development version of shinyMethyl:
```{r}
library(devtools)
install_github("rstudio/shiny")
install_github("jfortin1/shinyMethyl")
install_github("jfortin1/shinyMethylData")
```

Vignette
------------
You can find the vignette for `shinyMethyl` at [https://github.com/Jfortin1/shinyMethyl/blob/master/vignettes/shinyMethyl.pdf](https://github.com/Jfortin1/shinyMethyl/blob/master/vignettes/shinyMethyl.pdf)

Quick example
------------
After installation, you can launch `shinyMethyl` with an example dataset from TCGA with the following code:
```{r}
library(shinyMethyl)
library(shinyMethylData)
runShinyMethyl(summary.tcga.raw, summary.tcga.norm)
```

# Citation


To cite package __shinyMethyl__ in publications use:

Fortin J, Fertig EJ and Hansen KD (2014). “shinyMethyl: interactive
quality control of Illumina 450k DNA methylation arrays in R.”, F1000Research

A BibTeX entry for LaTeX users is

 @Article{shinymethyl,
    author = {Jean-Philippe Fortin and Elana J. Fertig and Kasper D. Hansen},
    title = {{shinyMethyl: interactive quality control of Illumina 450k DNA methylation arrays in R}},
    journal = {F1000Research},
    year = {2014},
    volume = {3},
    pages = {175},
    doi = {10.12688/f1000research.4680.2},
    pubmed = {25285208}
  }
