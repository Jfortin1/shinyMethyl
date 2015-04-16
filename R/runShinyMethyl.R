runShinyMethyl <- function(shinyMethylSet1, shinyMethylSet2 = NULL) {
	directory <- system.file(package = "shinyMethyl", "shinyMethyl")
	if (shinyMethylSet1@originObject != "RGChannelSet") {
		stop("First shinyMethylSet object must be created from an RGChannelSet")
	}



	#shinyMethylSet1   <<- shinyMethyl::orderByName(shinyMethylSet1)
	shinyMethylSet1 <- shinyMethyl::orderByName(shinyMethylSet1)
	if (!is.null(shinyMethylSet2)){
		shinyMethylSet2 <- shinyMethyl::orderByName(shinyMethylSet2)
	}

	## If a second shinyMethylSet is provided (from GenomicRatioSet):
	if (!is.null(shinyMethylSet2)) {
		if (shinyMethylSet2@originObject != "GenomicRatioSet") {
			stop("Second shinyMethylSet object must be created from a GenomicRatioSet")
		}
		## Need to make sure both shinyMethylSet's are compatible
		if (all.equal(shinyMethylSet1@sampleNames, shinyMethylSet2@sampleNames)!=TRUE) {
			stop("The two shinyMethylSet objects are not compatible. Both shinyMethylSet objects must be created from the same samples.")
		}
		#shinyMethylSet2 <<- shinyMethyl::orderByName(shinyMethylSet2)
	} else {
		#shinyMethylSet2 <<- NULL
		shinyMethylSet2 <- NULL
	}

	#source(paste0(directory, "/", "ui.R"))
	#source(paste0(directory, "/", "server.R"))

	shinyMethyl.app <- list(ui = ui.shinyMethyl(shinyMethylSet1, shinyMethylSet2), server = server.shinyMethyl(shinyMethylSet1, 
		shinyMethylSet2))
	#runApp(directory)
	runApp(shinyMethyl.app)
}
