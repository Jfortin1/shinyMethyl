#' @export
setGeneric("shinySummarize",
           function(object, ...) standardGeneric("shinySummarize"))



#' @importFrom minfi getAnnotation
#' @importFrom minfi getRed getGreen
#' @importFrom minfi getControlAddress
#' @importFrom minfi preprocessRaw getMeth getUnmeth
#' @importFrom minfi mapToGenome
#' @importFrom BiocGenerics updateObject
#' @importFrom stats prcomp
#' @importFrom MatrixGenerics rowVars colQuantiles
#' @export
setMethod("shinySummarize",
          signature(object = "RGChannelSet"),
          function(object){

    .createIndices <- function(object, betaMatrix) {
        ann  <- getAnnotation(object)
        chr  <- ann$chr
        probeType  <- paste0(ann$Type,ann$Color)
        probeNames <- rownames(ann)
        autosomal  <- chr %in% paste0("chr", 1:22)
        indices <- list(IGrn = probeNames[probeType == "IGrn" & autosomal],
                        IRed = probeNames[probeType == "IRed" & autosomal],
                        II = probeNames[probeType == "II" & autosomal],
                        X = probeNames[chr == "chrX"],
                        Y = probeNames[chr == "chrY"])
        for (i in seq_along(indices)){
            indices[[i]] <- which(rownames(betaMatrix) %in% indices[[i]])
        }
        indices
    }
        
    controlType <- c("BISULFITE CONVERSION I",
                     "BISULFITE CONVERSION II",
                     "EXTENSION",
                     "HYBRIDIZATION",
                     "NEGATIVE",
                     "NON-POLYMORPHIC",
                     "NORM_A",
                     "NORM_C",
                     "NORM_G",
                     "NORM_T", 
                     "SPECIFICITY I",
                     "SPECIFICITY II",
                     "TARGET REMOVAL",
                     "STAINING")

    object <- updateObject(object)
    cat("[shinySummarize] Extracting Red and Green channels \n")
    r <- minfi::getRed(object)
    g <- minfi::getGreen(object)

    ## Extraction of the controls
    greenControls <- redControls <- vector("list",length(controlType))
    names(greenControls) <- names(redControls) <- controlType

    for (i in 1:length(controlType)){
        ctrlAddress <- minfi::getControlAddress(object,
                                                controlType=controlType[i])
        redControls[[i]]   <- r[ctrlAddress,,drop=FALSE]
        greenControls[[i]] <- g[ctrlAddress,,drop=FALSE]
    }

    rm(r)
    rm(g)
    gc(verbose=FALSE)

    cat("[shinySummarize] Raw preprocessing \n")
    mset <- minfi::preprocessRaw(object)
    methMatrix <- minfi::getMeth(mset) 
    unmethMatrix <- minfi::getUnmeth(mset)
    rm(mset)
    betaMatrix <- methMatrix/(methMatrix+unmethMatrix+100)
    mMatrix    <- log2((methMatrix+1)/(unmethMatrix+1))
    cnMatrix   <- log2(methMatrix+unmethMatrix)

    cat("[shinySummarize] Mapping to genome \n")
    gmSet <- minfi::mapToGenome(object)

    probe.indices <- .createIndices(gmSet, betaMatrix)
    autosomal <- unlist(probe.indices[1:3])

    rm(gmSet)
    gc(verbose=FALSE)
    probs <- seq(0,1,length.out=500)
    betaQuantiles <- vector("list", length(probe.indices)) 
    names(betaQuantiles) <- names(probe.indices)
    mQuantiles <- methQuantiles <- unmethQuantiles <- cnQuantiles <- betaQuantiles

    cat("[shinySummarize] Computing quantiles \n") 
    for (i in seq_along(probe.indices)){
        betaQuantiles[[i]] <- t(colQuantiles(betaMatrix[probe.indices[[i]],],
                                             probs=probs, na.rm=TRUE))
        mQuantiles[[i]] <- t(colQuantiles(mMatrix[probe.indices[[i]],],
                                          probs=probs, na.rm=TRUE))
        methQuantiles[[i]] <- t(colQuantiles(methMatrix[probe.indices[[i]],],
                                             probs=probs, na.rm=TRUE))
        unmethQuantiles[[i]] <- t(colQuantiles(unmethMatrix[probe.indices[[i]],],
                                               probs=probs,na.rm=TRUE))
        cnQuantiles[[i]] <- t(colQuantiles(cnMatrix[probe.indices[[i]],],
                                           probs=probs, na.rm=TRUE))
        
        colnames(betaQuantiles[[i]])   <- colnames(object)
        colnames(mQuantiles[[i]])      <- colnames(object)
        colnames(methQuantiles[[i]])   <- colnames(object)
        colnames(unmethQuantiles[[i]]) <- colnames(object)
        colnames(cnQuantiles[[i]])     <- colnames(object)
    }

    rm(methMatrix)
    rm(unmethMatrix)
    rm(mMatrix)
    rm(cnMatrix)
    gc(verbose=FALSE)
    
  
    cat("[shinySummarize] Computing principal components \n")
    numPositions <- 20000
    autMatrix <- betaMatrix[autosomal,,drop=FALSE]
    rm(betaMatrix)
    gc(verbose=FALSE)
    o <- order(-rowVars(autMatrix))[seq_len(numPositions)]
    pca <- prcomp(t(autMatrix[o,]))
    pca  <-  list(scores=pca$x,
                  percs=(pca$sdev^2)/(sum(pca$sdev^2))*100)
    names(pca$percs) <- colnames(object)
    
    
    
    object <-  shinyMethylSet(sampleNames = colnames(object),
                              phenotype   = as.data.frame(minfi::pData(object)),
                              mQuantiles  = mQuantiles,
                              betaQuantiles = betaQuantiles,
                              methQuantiles = methQuantiles,
                              unmethQuantiles = unmethQuantiles,
                              cnQuantiles = cnQuantiles ,
                              greenControls = greenControls ,
                              redControls = redControls ,
                              pca = pca,
                              originObject = "RGChannelSet",
                              array = object@annotation[["array"]])
    return(object)
})



#' @importFrom minfi getBeta getM getCN
#' @export
setMethod("shinySummarize",
          signature(object = "GenomicRatioSet"),
          function(object){
    object <- updateObject(object)

    .createIndices <- function(object, betaMatrix) {
        ann  <- getAnnotation(object)
        chr  <- ann$chr
        probeType  <- paste0(ann$Type,ann$Color)
        probeNames <- rownames(ann)
        autosomal  <- chr %in% paste0("chr", 1:22)
        indices <- list(IGrn = probeNames[probeType == "IGrn" & autosomal],
                        IRed = probeNames[probeType == "IRed" & autosomal],
                        II = probeNames[probeType == "II" & autosomal],
                        X = probeNames[chr == "chrX"],
                        Y = probeNames[chr == "chrY"])
        for (i in seq_along(indices)){
            indices[[i]] <- which(rownames(betaMatrix) %in% indices[[i]])
        }
        indices
    }
    
    
    cat("[shinySummarize] Computing methylation values \n")
    betaMatrix <- minfi::getBeta(object)
    mMatrix    <- minfi::getM(object)
    cnMatrix   <- minfi::getCN(object)
    
    probe.indices <- .createIndices(object, betaMatrix)
    autosomal <- unlist(probe.indices[1:3])
    probs <- seq(0,1, length.out = 500)
    
    betaQuantiles <- vector("list", length(probe.indices)) 
    names(betaQuantiles) <- names(probe.indices)
    mQuantiles <- cnQuantiles <-  betaQuantiles 
    
    cat("[shinySummarize] Computing quantiles \n")   
    for (i in 1:length(probe.indices)){
        betaQuantiles[[i]] <- t(colQuantiles(betaMatrix[probe.indices[[i]],],
                                             probs=probs, na.rm=TRUE))
        mQuantiles[[i]] <- t(colQuantiles(mMatrix[probe.indices[[i]],],
                                          probs=probs, na.rm=TRUE))
        cnQuantiles[[i]] <- t(colQuantiles(cnMatrix[probe.indices[[i]],],
                                           probs=probs, na.rm=TRUE))
        names(betaQuantiles[[i]]) <- colnames(object) 
        names(mQuantiles[[i]])    <- colnames(object)
        names(cnQuantiles[[i]])   <- colnames(object)
    }
    
    cat("[shinySummarize] Computing principal components \n")
    ## To compute the principal components:
    numPositions <- 20000
    autMatrix <- betaMatrix[autosomal,,drop=FALSE]
    rm(betaMatrix)
    o <- order(-rowVars(autMatrix))[seq_len(numPositions)]
    pca <- prcomp(t(autMatrix[o,]))
    pca <- list(scores=pca$x,
                percs=(pca$sdev^2)/(sum(pca$sdev^2))*100)
    names(pca$percs) <- colnames(object)

    object <- shinyMethylSet(sampleNames = colnames(object),
                             phenotype   = as.data.frame(minfi::pData(object)),
                             mQuantiles  = mQuantiles,
                             betaQuantiles = betaQuantiles,
                             methQuantiles = list(NULL),
                             unmethQuantiles = list(NULL),
                             cnQuantiles = cnQuantiles,
                             greenControls = list(NULL),
                             redControls = list(NULL),
                             pca = pca,
                             originObject = "GenomicRatioSet",
                             array = object@annotation[["array"]])
    return(object)
})
