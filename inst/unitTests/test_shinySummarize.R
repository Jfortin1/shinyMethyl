test_shinySummarize <- function(){
	stopifnot(require(minfiData))
	stopifnot(require(minfi))
    stopifnot(require(digest))
    stopifnot(require(shinyMethyl))
    load(file.path(path.package("shinyMethyl"), "unitTests", "testDigests.rda"))
    shinyMethylSet <- shinySummarize(RGsetEx)

    beta = getBeta(shinyMethylSet)
	m    = getM(shinyMethylSet)
	meth = getMeth(shinyMethylSet)
	unmeth = getUnmeth(shinyMethylSet)
	cn = getCN(shinyMethylSet)
	green <- getGreenControls(shinyMethylSet)
	red <- getRedControls(shinyMethylSet)
	pca <- shinyMethylSet@pca

	all.matrices <- c(beta,m,meth,unmeth,cn, green, red, pca)
	currentDigests <- lapply(all.matrices, minfi:::.digestMatrix)
	for (i in 1:length(currentDigests)){
		checkEquals(testDigests[[i]], currentDigests[[i]])
	}

	checkEquals(names(currentDigests), names(testDigests))
}