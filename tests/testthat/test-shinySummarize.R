#library(minfiData)
#library(minfi)
#library(shinyMethyl)
   
# test_that("shinySummarize function", {
    
#     skip_on_os(os="windows")

#     data(RGsetEx, package="minfiData")
#     shinyMethylSet <- shinySummarize(RGsetEx)
#     beta <- getBeta(shinyMethylSet)
#     m    <- getM(shinyMethylSet)
#     meth <- getMeth(shinyMethylSet)
#     unmeth <- getUnmeth(shinyMethylSet)
#     cn <- getCN(shinyMethylSet)
#     green <- getGreenControls(shinyMethylSet)
#     red <- getRedControls(shinyMethylSet)
#     all.matrices <- c(beta=beta,
#                       m=m,
#                       meth=meth,
#                       unmeth=unmeth,
#                       cn=cn,
#                       green=green,
#                       red=red)
#     currentDigests <- lapply(all.matrices, minfi:::.digestMatrix)

#     for (i in seq_along(currentDigests)){
#         name <- paste0(names(currentDigests)[i], ".rds")
#         name <- gsub(" ", "_", name)
#         file <- file.path("objects", name)
#         expect_equal_to_reference(currentDigests[[i]],
#                                   file=file)
#     }
# })


