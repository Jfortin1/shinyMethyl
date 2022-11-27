library(minfiData)
library(minfi)
library(shinyMethyl)
   
test_that("shinySummarize function", {
    
    data(RGsetEx, package="minfiData")
    shinyMethylSet <- shinySummarize(RGsetEx)
    beta <- getBeta(shinyMethylSet)
    m    <- getM(shinyMethylSet)
    meth <- getMeth(shinyMethylSet)
    unmeth <- getUnmeth(shinyMethylSet)
    cn <- getCN(shinyMethylSet)
    green <- getGreenControls(shinyMethylSet)
    red <- getRedControls(shinyMethylSet)
    pca <- shinyMethylSet@pca
    all.matrices <- c(beta=beta,
                      m=m,
                      meth=meth,
                      unmeth=unmeth,
                      cn=cn,
                      green=green,
                      red=red,
                      pca=pca)
    currentDigests <- lapply(all.matrices, minfi:::.digestMatrix)

    for (i in 1:length(currentDigests)){
        file <- file.path("objects",names(currentDigests)[i])
        expect_equal_to_reference(currentDigests[[i]],file=file)
    }
})
