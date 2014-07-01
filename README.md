shinyMethyl
===========

Welcome to the shinyMethyl, an interactive R application based on the shiny package for exploration of DNA methylation data.

Online demo: http://spark.rstudio.com/jfortin/shinyMethyl/

Some updates: 

- Fixed returnPCScores() dependency on annotation package. Thanks to Maarten van Iterson 
- extractFromTargets450k() is temporarily non-available. Please use extractFromRGSet450k() instead. 
- Bug fixed: M-value densities are now visible
- With the update of minfi in Bioconductor 2.13, IlluminaHumanMethylation450kannotation.ilmn.v1.2 is no longer supported. ShinyMethyl is now updated to be annotation free (shinyMethyl does not really need an annotation). Thanks to Kathleen Fisch for pointing that out.
- Within the predictGender panel, it is now possible to download a .csv file containing not only the predicted gender, but also the actual gender if it was provided in the phenotype data, it will be also included in the file. The function also handles missing data. Thanks to Brent Pedersen


