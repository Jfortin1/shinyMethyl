library(minfiData)
library(shinyMethyl)
library(digest)
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
testDigests  <- lapply(all.matrices, minfi:::.digestMatrix)

save(testDigests, file = "../unitTests/testDigests.rda")

gc()
sessionInfo()                         

rm(list = ls())




