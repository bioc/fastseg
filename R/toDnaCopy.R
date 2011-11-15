#' Function to create a DNAcopy object for plot functions.
#' 
#' @param segData The results of the segmentation.
#' @param chrom The vector of the chromosomes from the original data.
#' @param maploc A vector with the physical positions of the original data.
#' @param genomdat A matrix with the original data.
#' @param sampleNames The sample names of the original data.
#' @return An DNAcopy equivalent object.
#' @author Andreas Mitterecker
#' @export
#' @examples 
#' library(fastseg)
#' 
#' #####################################################################
#' ### the data
#' #####################################################################
#' data(coriell)
#' head(coriell)
#' 
#' samplenames <- colnames(coriell)[4:5]
#' data <- as.matrix(coriell[4:5])
#' data[is.na(data)] <- median(data, na.rm=TRUE)
#' chrom <- coriell$Chromosome
#' maploc <- coriell$Position
#' 
#' 
#' ###########################################################
#' ## GRanges 
#' ###########################################################
#' 
#' library("GenomicRanges")
#' 
#' ## with both individuals
#' gr <- GRanges(seqnames=chrom,
#'         ranges=IRanges(maploc, end=maploc))
#' elementMetadata(gr) <- data
#' colnames(elementMetadata(gr)) <- samplenames
#' res <- fastseg(gr)
#' 
#' segres <- toDNAcopyObj(
#'         segData     = res, 
#'         chrom       = as.character(seqnames(gr)), 
#'         maploc      = as.numeric(start(gr)), 
#'         genomdat    = data, 
#'         sampleNames = samplenames)
#' 
#' ## with one individual
#' gr2 <- gr
#' data2 <- as.matrix(data[, 1])
#' colnames(data2) <- "sample1"
#' elementMetadata(gr2) <- data2
#' res <- fastseg(gr2)
#' 
#' segres <- toDNAcopyObj(
#'         segData     = res, 
#'         chrom       = as.character(seqnames(gr)), 
#'         maploc      = as.numeric(start(gr)), 
#'         genomdat    = as.matrix(data2), 
#'         sampleNames = unique(elementMetadata(res)$ID))
#' 
#' 
#' ###########################################################
#' ## vector
#' ###########################################################
#' data2 <- data[, 1]
#' res <- fastseg(data2)
#' segres <- toDNAcopyObj(
#'         segData     = res, 
#'         chrom       = rep(1, length(data2)), 
#'         maploc      = 1:length(data2), 
#'         genomdat    = as.matrix(data2), 
#'         sampleNames = "sample1")
#' 
#' 
#' ###########################################################
#' ## matrix
#' ###########################################################
#' data2 <- data[1:400, ]
#' res <- fastseg(data2)
#' segres <- toDNAcopyObj(
#'         segData     = res, 
#'         chrom       = rep(1, nrow(data2)), 
#'         maploc      = 1:nrow(data2), 
#'         genomdat    = as.matrix(data2), 
#'         sampleNames = colnames(data2))
#' 
#' 
#' ###########################################################
#' ## Expression set object
#' ###########################################################
#' library(oligo)
#' eSet <- new("ExpressionSet")
#' assayData(eSet) <- list(intensity=data)
#' 
#' featureData(eSet) <- new("AnnotatedDataFrame", 
#'         data=data.frame(
#'                 chrom = chrom,
#'                 start = maploc, 
#'                 end   = maploc))
#' phenoData(eSet) <- new("AnnotatedDataFrame", 
#'         data=data.frame(samples=samplenames))
#' sampleNames(eSet) <- samplenames
#' res <- fastseg(eSet)
#' 
#' segres <- toDNAcopyObj(
#'         segData     = res, 
#'         chrom       = rep(1, nrow(data)), 
#'         maploc      = maploc, 
#'         genomdat    = as.matrix(data), 
#'         sampleNames = colnames(data))
#' 
#' 
#' #####################################################################
#' ### plot the segments
#' #####################################################################
#' 
#' library(DNAcopy)
#' plot(segres)
toDNAcopyObj <- function(segData, chrom, maploc, genomdat, sampleNames) {

    if(!all(sampleNames %in% elementMetadata(segData)$ID)) {
        stop("The sample names of the segments do not fit these of the original data!")
    }

    genomdat <- as.matrix(genomdat)
    colnames(genomdat) <- sampleNames
    
    data.type <- "logratio"
    zzz <- data.frame(chrom=I(chrom), maploc=maploc, genomdat)
    attr(zzz, "data.type") <- data.type
    class(zzz) <- c("CNA", "data.frame")

    segDataTmp <- as.data.frame(segData)

	colnames(segDataTmp)[1:3] <- c("chrom", "loc.start", "loc.end")
	#segDataTmp$chrom <- as.numeric(segDataTmp$chrom)
	segDataTmp$ID <- as.character(segDataTmp$ID)
	
	
    segres <- list()
    segres$data <- zzz
    segres$output <- segDataTmp[, c(6, 1:3, 7:8)]
    #segres$segRows <- segDataTmp[, 9:10]
	segres$segRows <- segDataTmp[, 2:3]
    segres$call <- "unknown"    
    class(segres) <- "DNAcopy"

    return(segres)
}
