\name{fastseg}
\alias{fastseg}
\title{Detection of breakpoints using a fast segmentation algorithm based
on the cyber t-test.}
\usage{
  fastseg(x, type = 1, alpha = 0.05, segMedianT, minSeg = 4,
    eps = 0, delta = 5, maxInt = 40, squashing = 0,
    cyberWeight = 10)
}
\arguments{
  \item{x}{Values to be segmented either in the format of a
  sorted GRanges object, ExpressionSet object, matrix or
  vector.}

  \item{type}{Parameter that sets the type of test. If set
  to 1 a test of the left against the right window is
  performend. If set to 2 the segment is also tested
  against the global mean. (Default = 1).}

  \item{alpha}{A value between 0 and 1 is interpreted as
  the ratio of initial breakpoints. An integer greater than
  one is interpreted as number of desired breakpoints.
  Increasing this parameter leads to more segments.
  (Default = 0.1)}

  \item{segMedianT}{A numeric vector of length two with the
  thresholds of segments' median values that are considered
  as significant. Only segments with a median above the
  first or below the second value are kept in a final
  merging step. If missing the algorithm will try to find a
  reasonable value by using z-scores. (Default "missing".)}

  \item{minSeg}{The minimal segment length. (Default = 4).}

  \item{eps}{Minimal distance between consecutive values.
  Only consecutive values with a minimium distance of "eps"
  are tested. This makes the segmentation algorithm even
  faster. If all values should be tested "eps" can be set
  to zero. If missing the algorithm will try to find a
  reasonable value by using quantiles. (Default = 0.)}

  \item{delta}{Segment extension parameter. If delta
  consecutive extensions of the left and the right segment
  do not lead to a better p-value the testing is stopped.
  (Default = 5).}

  \item{maxInt}{Maximal length of the left and the right
  segment. (Default = 40).}

  \item{squashing}{The degree of squashing of the input
  values. If set to zero no squashing is performed.
  (Default = 0).}

  \item{cyberWeight}{The nu parameter of the cyber t-test.
  Can be interpreted as the weight of the global variance.
  The higher the value the more small segments with high
  variance will be significant. (Default = 10).}

  \item{...}{Further arguments passed to the plot
  function.}
}
\value{
  A data frame containing the segments.
}
\description{
  Detection of breakpoints using a fast segmentation
  algorithm based on the cyber t-test.
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
mcols(gr) <- data
colnames(mcols(gr)) <- samplenames
res <- fastseg(gr)

## with one individual
gr2 <- gr
data2 <- as.matrix(data[, 1])
colnames(data2) <- "sample1"
mcols(gr2) <- data2
res <- fastseg(gr2)


###########################################################
## vector
###########################################################
data2 <- data[, 1]
res <- fastseg(data2)



###########################################################
## matrix
###########################################################
data2 <- data[1:400, ]
res <- fastseg(data2)
}
\author{
  Guenter Klambauer \email{klambauer@bioinf.jku.at}
}

