\name{toDNAcopyObj}
\alias{toDNAcopyObj}
\title{Function to create a DNAcopy object for plot functions.}
\usage{
  toDNAcopyObj(segData, chrom, maploc, genomdat,
    sampleNames)
}
\arguments{
  \item{segData}{The results of the segmentation.}

  \item{chrom}{The vector of the chromosomes from the
  original data.}

  \item{maploc}{A vector with the physical positions of the
  original data.}

  \item{genomdat}{A matrix with the original data.}

  \item{sampleNames}{The sample names of the original
  data.}
}
\value{
  An DNAcopy equivalent object.
}
\description{
  Function to create a DNAcopy object for plot functions.
}
\examples{
library(fastseg)

#####################################################################
### the data
#####################################################################
data(coriell)
head(coriell)

samplenames <- colnames(coriell)[4:5]
data <- as.matrix(coriell[4:5])
data[is.na(data)] <- median(data, na.rm=TRUE)
chrom <- coriell$Chromosome
maploc <- coriell$Position


###########################################################
## GRanges
###########################################################

library("GenomicRanges")

## with both individuals
gr <- GRanges(seqnames=chrom,
        ranges=IRanges(maploc, end=maploc))
elementMetadata(gr) <- data
colnames(elementMetadata(gr)) <- samplenames
res <- fastseg(gr)

segres <- toDNAcopyObj(
        segData     = res,
        chrom       = as.character(seqnames(gr)),
        maploc      = as.numeric(start(gr)),
        genomdat    = data,
        sampleNames = samplenames)

## with one individual
gr2 <- gr
data2 <- as.matrix(data[, 1])
colnames(data2) <- "sample1"
elementMetadata(gr2) <- data2
res <- fastseg(gr2)

segres <- toDNAcopyObj(
        segData     = res,
        chrom       = as.character(seqnames(gr)),
        maploc      = as.numeric(start(gr)),
        genomdat    = as.matrix(data2),
        sampleNames = unique(elementMetadata(res)$ID))


###########################################################
## vector
###########################################################
data2 <- data[, 1]
res <- fastseg(data2)
segres <- toDNAcopyObj(
        segData     = res,
        chrom       = rep(1, length(data2)),
        maploc      = 1:length(data2),
        genomdat    = as.matrix(data2),
        sampleNames = "sample1")


###########################################################
## matrix
###########################################################
data2 <- data[1:400, ]
res <- fastseg(data2)
segres <- toDNAcopyObj(
        segData     = res,
        chrom       = rep(1, nrow(data2)),
        maploc      = 1:nrow(data2),
        genomdat    = as.matrix(data2),
        sampleNames = colnames(data2))



#####################################################################
### plot the segments
#####################################################################

library(DNAcopy)
plot(segres)
}
\author{
  Andreas Mitterecker
}

