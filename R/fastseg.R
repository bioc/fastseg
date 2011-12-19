# Copyright (C) 2011 Klambauer Guenter 
# <klambauer@bioinf.jku.at>

#' Detection of breakpoints using a fast segmentation algorithm based
#' on the cyber t-test.
#' 
#' @param x Values to be segmented either in the format of a sorted GRanges 
#' object, ExpressionSet object, matrix or vector.
#' @param type Parameter that sets the type of test. If set to 1 a test of 
#' the left against the right window is performend. If set to 2 the segment 
#' is also tested against the global mean. (Default = 1).
#' @param alpha A value between 0 and 1 is interpreted as the ratio of
#' initial breakpoints. An integer greater than one is interpreted as number
#' of desired breakpoints. Increasing this parameter leads to more segments. 
#' (Default = 0.05)
#' @param segMedianT A numeric vector of length two with the thresholds of 
#' segments' median values that are considered as significant. Only segments
#' with a median above the first or below the second value are kept in a final
#' merging step. (Default = "missing").
#' @param minSeg The minimal segment length. (Default = 4). 
#' @param eps Minimal difference between consecutive values. Only consecutive 
#' values with a minimium difference of "eps" are tested. This makes the 
#' segmentation algorithm even faster. If all values should be tested "eps" can
#' be set to zero. (Default = 0).
#' @param delta Segment extension parameter. If delta consecutive extensions
#' of the left and the right segment do not lead to a better p-value the testing
#' is stopped. (Default = 5).
#' @param maxInt Maximal length of the left and the right segment. (Default =
#' 40).
#' @param squashing  The degree of squashing of the input values. If set to zero 
#' no squashing is performed. (Default = 0).
#' @param cyberWeight The nu parameter of the cyber t-test. Can be interpreted
#' as the weight of the global variance. The higher the value the more small 
#' segments with high variance will be significant. (Default = 10).
#' @param segPlot Whether a standardized plot of the result should be plotted.
#' Set to FALSE if plotting to file is desired. (Default = TRUE).  
#' @param ... Further arguments passed to the plot function.
#' @examples 
#' x <- rnorm(n=500,sd=0.5)
#' x[150:200] <- rnorm(n=51,mean=3,sd=0.5)
#' fastseg(x)
#' @return A data frame containing the segments.
#' @author Guenter Klambauer \email{klambauer@@bioinf.jku.at}
#' @noRd
segmentGeneral <- function(x, type = 2, alpha = 0.05, segMedianT, minSeg = 4, 
		eps=0, delta = 5, maxInt = 10, squashing = 0, cyberWeight = 10, 
        segPlot = TRUE, ...) {
    
    if (missing("segMedianT")) {
        segMedianT <- c()
        segMedianT[1] <- median(x, na.rm=TRUE)+1.5*sd(x, na.rm=TRUE)
        segMedianT[2] <- median(x, na.rm=TRUE)-1.5*sd(x, na.rm=TRUE)
    } else {
        if (length(segMedianT)==1){
            segMedianT <- c(abs(segMedianT), -abs(segMedianT))
        }
    }
    
	if (any(is.na(x))) {
		x[is.na(x)] <- median(x, na.rm=TRUE)
	}
	
	if (missing("eps")) {
		eps <- quantile(abs(diff(x)), probs=0.75)
	}
	
	if (type==1) {
		res <-  .Call("segment", x, as.double(eps), as.integer(delta), 
				as.integer(maxInt), as.integer(minSeg),
				as.integer(squashing), as.double(cyberWeight))
	} else if (type==2) {
		res <- .Call("segmentCyberT", x, as.double(eps), as.integer(3), 
				as.integer(delta), as.integer(maxInt), as.integer(minSeg), 
				as.integer(squashing), as.double(cyberWeight), 0.15)
	}
	
	if (alpha >= 1) {
		alpha <- as.integer(alpha)
		brkptsInit <- sort(order(res$stat,decreasing=TRUE)[1:alpha])
	} else if (alpha < 1 & alpha > 0) {
		pValT <- quantile(res$stat, probs=1-alpha)
		brkptsInit <- which(res$stat > pValT)
	} else {
		stop(paste("Alpha must be either between 0 and 1 or an integer",
						"greater than 1."))
	}
    
	brkptsInit <- unique(c(0, brkptsInit, length(x)))
	nbrOfBrkpts <- length(brkptsInit)-1
	
	med <- vector(length=nbrOfBrkpts)
	m <- vector(length=nbrOfBrkpts)
	start <- vector(length=nbrOfBrkpts)
	end <- vector(length=nbrOfBrkpts)
    
	for (i in seq_len(nbrOfBrkpts)){
		y <- x[(brkptsInit[i]+1):brkptsInit[i+1]]
		m[i] <- mean(y)
		med[i] <- median(y)
		start[i] <- brkptsInit[i]+1
		end[i] <- brkptsInit[i+1]
	}
	
	df <- data.frame("start"=start, "end"=end, "mean"=m, "median"=med)
	
	if (all(segMedianT==0)) {
		
		#message("No merging of segments.")
		ir <- IRanges::IRanges(df$start, df$end)
		ir <- ir[which(width(ir)>=minSeg)]
		
		
		irAll <- IRanges::IRanges(1, length(x))
		segsFinal <- IRanges::as.data.frame(IRanges::sort(
						c(ir, IRanges::setdiff(irAll, ir))))
		
		if (segPlot) plot(x, pch=15, cex=0.5, ...)
		nbrOfSegs <- nrow(segsFinal)
		med <- vector(length=nbrOfSegs)
		m <- vector(length=nbrOfSegs)
		start <- vector(length=nbrOfSegs)
		end <- vector(length=nbrOfSegs)
		for (i in seq_len(nbrOfSegs)) {
			y <- x[segsFinal$start[i]:segsFinal$end[i]]
			m[i] <- mean(y)
			med[i] <- median(y)
			
			if (segPlot) lines(c(segsFinal$start[i], segsFinal$end[i]), 
						c(med[i], med[i]), lwd=5, col="green")
		}
		df2 <- data.frame("start"=segsFinal$start, "end"=segsFinal$end, 
				"mean"=m, "median"=med)
		
		l <- list(df2, df)
		names(l) <- c("finalSegments", "unmergedSegments")
		return(l)
		
	} else {
		dfAmp <- df[which(df$median > segMedianT[1]), ]
		irAmp <- IRanges::IRanges(dfAmp$start, dfAmp$end)
		irAmp <- IRanges::reduce(irAmp)
		
		dfLoss <- df[which(df$median < segMedianT[2]), ]
		irLoss <- IRanges::IRanges(dfLoss$start, dfLoss$end)
		irLoss <- IRanges::reduce(irLoss)
		
		ir <- IRanges::sort(c(irAmp, irLoss))
        
		ir <- ir[which(width(ir)>=minSeg)]
		
		rm(irAmp, irLoss, dfAmp, dfLoss)    
		
		irAll <- IRanges::IRanges(1, length(x))
		segsFinal <- IRanges::as.data.frame(IRanges::sort(
						c(ir, IRanges::setdiff(irAll, ir))))
		
		if (segPlot) plot(x, pch=15, cex=0.5, ...)
		nbrOfSegs <- nrow(segsFinal)
		med <- vector(length=nbrOfSegs)
		m <- vector(length=nbrOfSegs)
		start <- vector(length=nbrOfSegs)
		end <- vector(length=nbrOfSegs)
		for (i in seq_len(nbrOfSegs)) {
			y <- x[segsFinal$start[i]:segsFinal$end[i]]
			m[i] <- mean(y)
			med[i] <- median(y)
			
			if (segPlot) lines(c(segsFinal$start[i], segsFinal$end[i]), 
						c(med[i], med[i]), lwd=5, col="green")
		}
		
		df2 <- data.frame("start"=segsFinal$start, "end"=segsFinal$end, 
				"mean"=m, "median"=med)
		
		l <- list(df2, df)
		names(l) <- c("finalSegments", "unmergedSegments")
		return(l)
	}
}



#' Detection of breakpoints using a fast segmentation algorithm based
#' on the cyber t-test.
#' 
#' @param x Values to be segmented either in the format of a sorted GRanges 
#' object, ExpressionSet object, matrix or vector.
#' @param type Parameter that sets the type of test. If set to 1 a test of 
#' the left against the right window is performend. If set to 2 the segment 
#' is also tested against the global mean. (Default = 1).
#' @param alpha A value between 0 and 1 is interpreted as the ratio of
#' initial breakpoints. An integer greater than one is interpreted as number
#' of desired breakpoints. Increasing this parameter leads to more segments. 
#' (Default = 0.05)
#' @param segMedianT A numeric vector of length two with the thresholds of 
#' segments' median values that are considered as significant. Only segments
#' with a median above the first or below the second value are kept in a final
#' merging step. If missing the algorithm will try to find a reasonable value
#' by using z-scores. (Default "missing".)
#' @param minSeg The minimal segment length. (Default = 4). 
#' @param eps Minimal distance between consecutive values. Only consecutive 
#' values with a minimium distance of "eps" are tested. This makes the 
#' segmentation algorithm even faster. If all values should be tested "eps" can
#' be set to zero. If missing the algorithm will try to find a reasonable value
#' by using quantiles. (Default "missing".)
#' @param delta Segment extension parameter. If delta consecutive extensions
#' of the left and the right segment do not lead to a better p-value the testing
#' is stopped. (Default = 5).
#' @param maxInt Maximal length of the left and the right segment. (Default =
#' 40).
#' @param squashing  The degree of squashing of the input values. If set to zero 
#' no squashing is performed. (Default = 0).
#' @param cyberWeight The nu parameter of the cyber t-test. Can be interpreted
#' as the weight of the global variance. The higher the value the more small 
#' segments with high variance will be significant. (Default = 10). 
#' @param ... Further arguments passed to the plot function.
#' @return A data frame containing the segments.
#' @author Guenter Klambauer \email{klambauer@@bioinf.jku.at}
#' @export
#' @useDynLib fastseg
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
#' ## with one individual
#' gr2 <- gr
#' data2 <- as.matrix(data[, 1])
#' colnames(data2) <- "sample1"
#' elementMetadata(gr2) <- data2
#' res <- fastseg(gr2)
#' 
#' 
#' ###########################################################
#' ## vector
#' ###########################################################
#' data2 <- data[, 1]
#' res <- fastseg(data2)
#' 
#' 
#' 
#' ###########################################################
#' ## matrix
#' ###########################################################
#' data2 <- data[1:400, ]
#' res <- fastseg(data2)
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

fastseg <- function(x, type = 1, alpha = 0.1, segMedianT, minSeg = 4, 
        eps = 0, delta = 5, maxInt = 40, squashing = 0, cyberWeight = 10, ...) {
    
    if (inherits(x, "ExpressionSet")) {
        
        if (!("intensity" %in% names(assayData(x)))) {
            stop("ExpressionSet needs to have an assayData slot named intensity!")
        } 
        if(!all(colnames(fData(x))[1:3] == c("chrom", "start", "end"))) {
            stop("The first 3 colnames of featureData need to be names: chrom, start, end")
        }
        if(length(sampleNames(x)) != ncol(assayData(x)$intensity)) {
            stop("sampleNames must be assigned and have a correct dimension!")
        }
        
        y <- split(x, featureData(x)$chrom)
        nbrOfSeq <- length(y)
        
        res02 <- list()
        for (seq in seq_len(nbrOfSeq)) {
            x <- y[[seq]]
            
            res <- list()
            for (sampleIdx in seq_len(ncol(assayData(x)$intensity))) {
                z01 <- assayData(x)$intensity[, sampleIdx]
                sample <- sampleNames(x)[sampleIdx]
                
                resTmp <- segmentGeneral(z01, type, alpha, segMedianT, minSeg, 
                        eps, delta, maxInt, squashing, cyberWeight, 
                        segPlot = FALSE, ...)$finalSegments
                resTmp$sample <- sample
                res[[sampleIdx]] <- resTmp
            }
            res <- do.call("rbind", res)
            res$num.mark <- res$end - res$start
            
            chrom <- rep(featureData(x)$chrom[1], nrow(res))
            start <- featureData(x)$start[res$start]
            end <- featureData(x)$end[res$end]
            
            resX <- data.frame(
                    ID = res$sample, 
                    chrom = chrom, 
                    loc.start = start, 
                    loc.end = end, 
                    num.mark = res$num.mark, 
                    seg.mean = res$mean, 
                    startRow = res$start, 
                    endRow = res$end,stringsAsFactors=FALSE)
            
            res02[[seq]] <- resX
            
        }
        res03 <- do.call("rbind", res02)
        
        finalRes <- GRanges(seqnames = Rle(res03$chrom),
                ranges   = IRanges(start = res03$loc.start, end = res03$loc.end),
                ID = res03$ID, 
                num.mark  = res03$num.mark,
                seg.mean  = res03$seg.mean,
                startRow  = res03$startRow,
                endRow    = res03$endRow)
                
    } else if (inherits(x, "GRanges")) {
#        if (!all(lapply(IRanges::elementMetadata(x), mode) == "numeric")) {
#            stop("All elementMetadata of GRanges object needs to be numeric!")
#		}
#        
        y <- split(x, as.character(seqnames(x)))
		
        nbrOfSeq <- length(y)
     
        res02 <- list()
        for (seq in seq_len(nbrOfSeq)) {
            x <- y[[seq]]
            
            res <- list()
            for (sampleIdx in seq_len(ncol(elementMetadata(x)))) {
                z01 <- elementMetadata(x)[[sampleIdx]]
                sample <- names(elementMetadata(x))[sampleIdx]
                resTmp <- segmentGeneral(z01, type, alpha, segMedianT, minSeg, 
                        eps, delta, maxInt, squashing, cyberWeight, 
                        segPlot = FALSE, ...)$finalSegments
                resTmp$sample <- sample
                res[[sampleIdx]] <- resTmp
            }
			
            res <- do.call("rbind", res)

            res$num.mark <- res$end - res$start

            chrom <- as.character(seqnames(x)[1])
            start <- start(x)[res$start]
            end <- start(x)[res$end] + width(x)[res$end]-1
            
            resX <- data.frame(
                    ID = res$sample, 
                    chrom = chrom, 
                    loc.start = start, 
                    loc.end = end, 
                    num.mark = res$num.mark, 
                    seg.mean = res$mean, 
                    startRow = res$start, 
                    endRow = res$end,stringsAsFactors=FALSE)
            resX <- resX[order(resX$chrom,resX$loc.start), ]
            res02[[seq]] <- resX
            
        }
		
        res03 <- do.call("rbind", res02)

        finalRes <- GRanges(seqnames = Rle(res03$chrom),
                ranges   = IRanges(start = res03$loc.start, end = res03$loc.end),
                ID = res03$ID, 
                num.mark  = res03$num.mark,
                seg.mean  = res03$seg.mean,
                startRow  = res03$startRow,
                endRow    = res03$endRow)

    } else if (is.matrix(x)) {
        nbrOfSamples <- ncol(x)
        samples <- colnames(x)
        
        segsTmp <- list()
        for (i in seq_len(nbrOfSamples)) {
            segsTmp[[i]] <- segmentGeneral(x[, i], type, alpha, segMedianT, minSeg, 
                    eps, delta, maxInt, squashing, cyberWeight, 
                    segPlot = FALSE, ...)$finalSegments
            segsTmp[[i]]$sample <- samples[i]
        }

        res02 <- do.call("rbind", segsTmp)

        finalRes <- GRanges(seqnames = Rle(rep(1, nrow(res02))),
                ranges   = IRanges(start = res02$start, end = res02$end),
                ID = res02$sample, 
                num.mark  = res02$end - res02$start + 1,
                seg.mean  = res02$mean,
                startRow  = res02$start,
                endRow    = res02$end)
        
    } else if (is.vector(x)) {
        
        res02 <- segmentGeneral(x, type, alpha, segMedianT, minSeg, 
                eps, delta, maxInt, squashing, cyberWeight, 
                segPlot = FALSE, ...)$finalSegments
       
        finalRes <- GRanges(seqnames = Rle(rep(1, nrow(res02))),
                ranges   = IRanges(start = res02$start, end = res02$end),
                ID = rep("sample1", nrow(res02)), 
                num.mark  = res02$end - res02$start + 1,
                seg.mean  = res02$mean,
                startRow  = res02$start,
                endRow    = res02$end)
        
    } else {
        
        stop("GRanges object, ExpressionSet object, vector or matrix as input expected!")    
        
    }
    
    finalRes <- finalRes[order(
                    as.character(elementMetadata(finalRes)$ID), 
                    (as.character(seqnames(finalRes))), 
                    as.numeric(start(finalRes))), ]
    
    return(finalRes)
    
}
