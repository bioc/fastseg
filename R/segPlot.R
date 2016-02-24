# Plot method by Venkatraman E. Seshan, Adam Olshen (DNAcopy, v1.29)
# modified for fastseg input.



#' Plots the data from a copy number array experiment (aCGH, ROMA etc.)
#' along with the results of segmenting it into regions of equal copy
#' numbers.
#' @param x The object that was segmented by fastseg.
#' @param res The result of fastseg.
#' @param plot.type the type of plot. (Default = "s").
#' @param altcol logical flag to indicate if chromosomes should be
#'   plotted in alternating colors in the whole genome plot. (Default = TRUE).
#' @param sbyc.layout \code{layout} settings for the multifigure grid layout
#'    for the `samplebychrom' type.  It should be specified as a vector of
#'    two integers which are the number of rows and columns.  The default
#'    values are chosen based on the number of chromosomes to produce a
#'    near square graph.   For normal genome it is 4x6 (24 chromosomes)
#'    plotted by rows. (Default = NULL).
#' @param cbys.layout \code{layout} settings for the multifigure grid layout
#'   for the `chrombysample' type.  As above it should be specified as
#'    number of rows and columns and the default chosen based on the
#'    number of samples. (Default = NULL).
#' @param cbys.nchrom the number of chromosomes per page in the layout.
#'(Default = 1).
#' @param include.means logical flag to indicate whether segment means
#'   are to be drawn. (Default = TRUE).
#' @param zeroline logical flag to indicate whether a horizontal line at
#'    y=0 is to be drawn. (Default = TRUE).
#' @param pt.pch the plotting character used for plotting the log-ratio
#'    values. (Default = ".")
#' @param pt.cex the size of plotting character used for the log-ratio
#'    values (Default = 3).
#' @param pt.cols the color list for the points. The colors alternate
#'    between chromosomes. (Default = c("green","black").
#' @param segcol the color of the lines indicating the segment means.
#' (Default = "red").
#' @param zlcol the color of the zeroline. (Default = "grey").
#' @param ylim this argument is present to override the default limits
#'    which is the range of symmetrized log-ratios. (Default = NULL).
#' @param lwd line weight of lines for segment mean and zeroline. (Default = 3).
#' @param ... other arguments which will be passed to \code{plot}
#'   commands.
#' @examples
#' data(coriell)
#' head(coriell)
#' samplenames <- colnames(coriell)[4:5]
#' data <- as.matrix(coriell[4:5])
#' chrom <- coriell$Chromosome
#' maploc <- coriell$Position
#' library("GenomicRanges")
#' gr <- GRanges(seqnames=chrom,
#' 		ranges=IRanges(maploc, end=maploc))
#' mcols(gr) <- data
#' colnames(mcols(gr)) <- samplenames
#' res <- fastseg(gr)
#' segPlot(gr,res)
#' @return A plot of the values and segments. 
#' @author klambaue
#' @export
segPlot <- function(x,res,
		plot.type="chrombysample", 
		altcol=TRUE,sbyc.layout=NULL,
		cbys.nchrom=1, cbys.layout=NULL, include.means=TRUE, zeroline=TRUE,
		pt.pch = NULL, pt.cex=NULL, pt.cols = NULL, segcol=NULL, 
		zlcol =NULL, ylim=NULL,lwd=NULL,...){
	
	if (is.vector(x)){
		xdat <- as.matrix(x,ncol=1)
		nsample <- 1
		chrom <- as.character(seqnames(res))
		if (length(unique(chrom))!=1){
			stop("Data points and result indicate different number of chromosomes.")	
		}
		chrom <- rep(unique(chrom),length(x))
		maploc <- 1:length(x)
		xmaploc=FALSE 
	} else if (is.matrix(x)){
		if (is.null(colnames(x))){
			colnames(x) <- paste("Sample",1:ncol(x),sep="_")
		}
		xdat <- x
		nsample <- ncol(x)
		chrom <- as.character(seqnames(res))
		if (length(unique(chrom))!=1){
			stop("Data points and result indicate different number of chromosomes.")	
		}
		chrom <- rep(unique(chrom),nrow(x))
		maploc <- 1:nrow(x)
		xmaploc=FALSE 
	} else if (inherits(x, "GRanges")){
		#cat("..")	
		xdat <- do.call("cbind",values(x)@listData)
		nsample <- ncol(values(x))
		chrom <- as.character(seqnames(x))
	
		#chrom <- rep(unique(chrom),length(x))
		maploc <- (  start(ranges(x)) + end(ranges(x))  )/2
		xmaploc=TRUE
	}else{
		stop("x must be of type GRanges, vector or matrix!")
		
	}	
	
	if(missing(ylim)) {
		uylim <- max(abs(xdat[which(is.finite(xdat))]), na.rm=TRUE)
		ylim <- c(-uylim, uylim)
	}
	
	
	xres <- data.frame("ID" = values(res)$ID,
	"num.mark" = (-values(res)$startRow+values(res)$endRow)+1,
	"chrom" = as.character(seqnames(res)),
	"seg.mean" = values(res)$seg.mean,stringsAsFactors=FALSE)
	xres <- xres[order(values(res)$ID,as.character(seqnames(res)),
				values(res)$startRow), ]
	
	if(dev.cur() <= 1) dev.new()
	int.dev <- dev.interactive()
	#plot.type <- match.arg(plot.type)
	op <- par(no.readonly = TRUE)
	parask <- par("ask")
	if (int.dev & !parask & nsample>1) par(ask = TRUE)
	sampleid <- colnames(xdat)
	if (is.null(sampleid)){sampleid <- paste("sample",1:nsample,sep="")}
	#browser()
	chrom0 <- chrom
	uchrom <- unique(chrom0)
	nchrom <- length(uchrom)
	if (xmaploc) {
		maploc0 <- as.numeric(maploc)
		if(length(uchrom)>1 & max(maploc0[chrom0==uchrom[1]]) > min(maploc0[chrom0==uchrom[2]])) {
			plen <- max(maploc0[chrom0==uchrom[1]])
			for(i in 2:nchrom) {
				maploc0[chrom0==uchrom[i]] <- plen + maploc0[chrom0==uchrom[i]]
				plen <- max(maploc0[chrom0==uchrom[i]])
			}
		}
	}
	if (missing(pt.pch)) pt.pch <- "."
	if (missing(pt.cex)) {
		if (pt.pch==".") { pt.cex <- 3}
		else {pt.cex <- 1}
	}
	wcol0 <- rep(1, length(chrom0))
	if (altcol) {
		j <- 0
		for (i in uchrom) {
			j <- (j+1) %% 2
			wcol0[chrom0==i] <- 1+j
		}
	}
	
	if (missing(pt.cols)) pt.cols <- c("black","green")
	if (missing(segcol)) segcol <- "red"
	if (missing(zlcol)) zlcol <- "grey"
	if (missing(lwd)) lwd <- 3
	if (plot.type == "chrombysample" | plot.type == "c" ) {
		#cat("Setting multi-figure configuration\n")
	#browser()	
	par(mar = c(0, 4, 0, 2), oma = c(4, 0, 4, 0), mgp = c(2, 0.7, 0))
		if (missing(cbys.layout)) {
			nrow <- ncol <- ceiling(sqrt(nsample))
			if (nrow*ncol - nsample > 0) {
				nrow <- nrow - 1
				ncol <- ncol + 1
			}
			if (nrow*ncol - nsample >= nrow) ncol <- ncol - 1
			cbys.layout <- c(nrow, ncol)
		}
		
		lmat0 <- lmat1 <- c(1:nsample, rep(-cbys.nchrom*nsample,
						max(prod(cbys.layout) - nsample,1)))
		for(i in 1:(cbys.nchrom-1)) {
			lmat1 <- c(lmat1,lmat0+nsample*i)
		}
		lmat1[lmat1<0] <- 0
		lmat <- matrix(lmat1[1:prod(cbys.layout)], nrow = cbys.layout[1], 
				ncol = cbys.nchrom*cbys.layout[2], byrow = FALSE)
		#lmat <- matrix(lmat1, nrow = cbys.layout[1], ncol = cbys.layout[2], byrow = FALSE)
		#cat("lmat: ",lmat,"\n")
		#cat("cbys: ",cbys.layout,"\n")
		
		layout(lmat)
	}
	if (plot.type == "samplebychrom"| plot.type == "s" ) {
		#cat("Setting multi-figure configuration\n")
		par(mar = c(4, 4, 4, 2), oma = c(0, 0, 2, 0), mgp = c(2, 0.7, 0))
		if (missing(sbyc.layout)) {
			nrow <- ncol <- ceiling(sqrt(nchrom))
			if (nrow*ncol - nchrom > 0) {
				nrow <- nrow - 1
				ncol <- ncol + 1
			}
			if (nrow*ncol - nchrom > ncol) ncol <- ncol - 1
			sbyc.layout <- c(nrow, ncol)
		}
		lmat <- matrix(c(1:nchrom, rep(0,prod(sbyc.layout)-nchrom)),
				nrow = sbyc.layout[1], ncol = sbyc.layout[2], byrow=TRUE)
		layout(lmat)
	}
	if (plot.type == "chrombysample" | plot.type == "c") {
		atchrom <- 0.5/cbys.nchrom
		for (ichrom in uchrom) {
			if (xmaploc) maploc1 <- maploc0[chrom0==ichrom]
			for (isamp in 1:nsample) {
				genomdat <- xdat[chrom0==ichrom, isamp]
				ina <- which(is.finite(genomdat))
				genomdat <- genomdat[ina]
				if (xmaploc) maploc <- maploc1[ina]
				ii <- cumsum(c(0, 
								xres$num.mark[xres$ID == sampleid[isamp] & xres$chrom==ichrom]))
				mm <- xres$seg.mean[xres$ID == sampleid[isamp] & xres$chrom==ichrom]
				kk <- length(ii)
				zz <- cbind(ii[-kk] + 1, ii[-1])
				if (xmaploc) {
					plot(maploc, genomdat, pch = pt.pch, cex=pt.cex, 
							xaxt="n", ylim = ylim, ylab = sampleid[isamp])
				} else {
					plot(genomdat, pch = pt.pch, cex=pt.cex, xaxt="n",
							ylim = ylim, ylab = sampleid[isamp])
				}
				if(zeroline) abline(h=0, col=zlcol, lwd=lwd)
				if (isamp%%cbys.layout[1] == 0) {
					axis(1, outer=TRUE)
					title(xlab="Index")
				}
				if (include.means) {
					for (i in 1:(kk - 1)) {
						if (xmaploc) { 
							lines(maploc[zz[i, ]], rep(mm[i], 2), col = segcol, lwd=lwd)
						} else {
							lines(zz[i, ], rep(mm[i], 2), col = segcol, lwd=lwd)
						}
					}
				}
			}
			mtext(paste("Chromosome",ichrom), side = 3, line = 1, 
					at = atchrom, outer=TRUE, font=2)
			atchrom <- atchrom + 1/cbys.nchrom
			atchrom <- atchrom - floor(atchrom)
		}
	} else {
		for (isamp in 1:nsample)
		{
			#browser()
			xres <- xres[order(match(xres$chrom,uchrom)), ]
			genomdat <- xdat[, isamp]
			ina <- which(is.finite(genomdat))
			genomdat <- genomdat[ina]
			wcol <- wcol0[ina]
			chrom <- chrom0[ina]
			if (xmaploc) maploc <- maploc0[ina]
			ii <- cumsum(c(0, xres$num.mark[xres$ID == sampleid[isamp]]))
			mm <- xres$seg.mean[xres$ID == sampleid[isamp]]
			kk <- length(ii)
			zz <- cbind(ii[-kk] + 1, ii[-1])
			if(missing(ylim)) ylim <- range(c(genomdat, -genomdat))
			if (plot.type=="whole" | plot.type == "w")
			{	#browser()
				if (xmaploc) {
					plot(maploc, genomdat, pch = pt.pch, cex=pt.cex, 
							col=pt.cols[wcol], main = sampleid[isamp], 
							ylab = "", ylim = ylim)
					if(zeroline) abline(h=0, col=zlcol, lwd=lwd)
				} else {
					plot(genomdat, pch = pt.pch, cex=pt.cex, 
							col=pt.cols[wcol], main = sampleid[isamp], 
							ylab = "", ylim = ylim)
					if(zeroline) abline(h=0, col=zlcol, lwd=lwd)
				}
				if (include.means) {
					for (i in 1:(kk - 1))
					{
						if (xmaploc) { 
							lines(maploc[zz[i, ]], rep(mm[i], 2),
									col = segcol, lwd=lwd)
						} else {
							lines(zz[i, ], rep(mm[i], 2), col = segcol, lwd=lwd)
						}
					}
				}
			}
			if (plot.type=="samplebychrom" | plot.type == "s")
			{
				#browser()
				
				cc <- xres$chrom[xres$ID == sampleid[isamp]]
				for (ichrom in uchrom)
				{
					if (xmaploc) {
						plot(maploc[chrom == ichrom], 
								genomdat[chrom == ichrom], pch = pt.pch, 
								cex=pt.cex, xlab="maploc", ylab = "", 
								main = paste("Chromosome", ichrom), ylim = ylim)
					} else {
						plot(genomdat[chrom == ichrom], pch = pt.pch, 
								cex=pt.cex, ylab = "", 
								main = paste("Chromosome", ichrom), ylim = ylim)
					}
					if(zeroline) abline(h=0, col=zlcol, lwd=lwd)
					if (include.means) {
						jj <- which(cc==ichrom)
						jj0 <- min(jj)
						for (i in jj)
						{
							if (xmaploc) {
								lines(maploc[zz[i, ]], rep(mm[i], 2), 
										col = segcol, lwd=lwd)
							} else {
								lines(1+zz[i, ]-zz[jj0,1], rep(mm[i], 2), 
										col = segcol, lwd=lwd)
							}
						}
					}
				}
				mtext(sampleid[isamp], side = 3, line = 0, outer = TRUE, font=2)
			}
			if (plot.type=="plateau" | plot.type == "p")
			{
				omm <- order(mm)
				ozz <- zz[omm,]
				ina <- unlist(apply(ozz, 1, function(ii) ii[1]:ii[2]))
				plot(genomdat[ina], pch = pt.pch, cex=pt.cex, 
						main = sampleid[isamp], ylab = "", ylim = ylim)
				if(zeroline) abline(h=0, col=zlcol, lwd=lwd)
				if (include.means) {
					ii <- cumsum(c(0, 
									xres$num.mark[xres$ID == sampleid[isamp]][omm]))
					smm <- mm[omm]
					zz <- cbind(ii[-kk] + 1, ii[-1])
					for (i in 1:(kk-1)) lines(zz[i, ], 
								rep(smm[i], 2), col = segcol, lwd=lwd)
				}
			}
		}
	}
}

