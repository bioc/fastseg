\name{segPlot}
\alias{segPlot}
\title{Plots the data from a copy number array experiment (aCGH, ROMA etc.)
along with the results of segmenting it into regions of equal copy
numbers.}
\usage{
  segPlot(x, res, plot.type = "chrombysample",
    altcol = TRUE, sbyc.layout = NULL, cbys.nchrom = 1,
    cbys.layout = NULL, include.means = TRUE,
    zeroline = TRUE, pt.pch = NULL, pt.cex = NULL,
    pt.cols = NULL, segcol = NULL, zlcol = NULL,
    ylim = NULL, lwd = NULL, ...)
}
\arguments{
  \item{x}{The object that was segmented by fastseg.}

  \item{res}{The result of fastseg.}

  \item{plot.type}{the type of plot. (Default = "s").}

  \item{altcol}{logical flag to indicate if chromosomes
  should be plotted in alternating colors in the whole
  genome plot. (Default = TRUE).}

  \item{sbyc.layout}{\code{layout} settings for the
  multifigure grid layout for the `samplebychrom' type.  It
  should be specified as a vector of two integers which are
  the number of rows and columns.  The default values are
  chosen based on the number of chromosomes to produce a
  near square graph.  For normal genome it is 4x6 (24
  chromosomes) plotted by rows. (Default = NULL).}

  \item{cbys.layout}{\code{layout} settings for the
  multifigure grid layout for the `chrombysample' type.  As
  above it should be specified as number of rows and
  columns and the default chosen based on the number of
  samples. (Default = NULL).}

  \item{cbys.nchrom}{the number of chromosomes per page in
  the layout. (Default = 1).}

  \item{include.means}{logical flag to indicate whether
  segment means are to be drawn. (Default = TRUE).}

  \item{zeroline}{logical flag to indicate whether a
  horizontal line at y=0 is to be drawn. (Default = TRUE).}

  \item{pt.pch}{the plotting character used for plotting
  the log-ratio values. (Default = ".")}

  \item{pt.cex}{the size of plotting character used for the
  log-ratio values (Default = 3).}

  \item{pt.cols}{the color list for the points. The colors
  alternate between chromosomes. If missing the point
  colors are black and green.}

  \item{segcol}{the color of the lines indicating the
  segment means. (Default = "red").}

  \item{zlcol}{the color of the zeroline. (Default =
  "grey").}

  \item{ylim}{this argument is present to override the
  default limits which is the range of symmetrized
  log-ratios. (Default = NULL).}

  \item{lwd}{line weight of lines for segment mean and
  zeroline. (Default = 3).}

  \item{...}{other arguments which will be passed to
  \code{plot} commands.}
}
\value{
  A plot of the values and segments.
}
\description{
  Plots the data from a copy number array experiment (aCGH,
  ROMA etc.) along with the results of segmenting it into
  regions of equal copy numbers.
}
\examples{
data(coriell)
head(coriell)
samplenames <- colnames(coriell)[4:5]
data <- as.matrix(coriell[4:5])
chrom <- coriell$Chromosome
maploc <- coriell$Position
library("GenomicRanges")
gr <- GRanges(seqnames=chrom,
		ranges=IRanges(maploc, end=maploc))
elementMetadata(gr) <- data
colnames(elementMetadata(gr)) <- samplenames
res <- fastseg(gr)
segPlot(gr,res)
}
\author{
  klambaue
}

