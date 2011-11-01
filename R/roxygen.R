#' Array CGH data set of Coriell cell lines
#'
#' These are two data array CGH studies sets of Corriel cell lines taken 
#' from the reference below.
#' @name coriell
#' @docType data
#' @references \url{http://www.nature.com/ng/journal/v29/n3/suppinfo/ng754\_S1.html}
#' @usage data(coriell)
#' @format A data frame containing five variables: first is clone name, 
#' second is clone chromosome, third is clone position, fourth and fifth are 
#' log2ratio for two cell lines.
#' @keywords data
#' @references Snijders et al., Assembly of microarrays for genome-wide 
#' measurement of DNA copy number, Nature Genetics, 2001 
NULL

#' Example data set for fastseg
#' 
#' The data is a small subset of copy number calls which were
#' produced by the cn.farms algorithm from an Affymetrix SNP microarray 
#' experiment of a HapMap sample.
#' 
#' @name fastsegData
#' @docType data
#' @usage data(fastsegData)
#' @format A simple vector with a copy number call as produced by the cn.farms algorithm.
#' @references \url{http://nar.oxfordjournals.org/content/early/2011/04/12/nar.gkr197.abstract}
#' Clevert et al., cn.FARMS: a latent variable model to detect copy number 
#' variations in microarray data with a low false discovery rate, NAR, 2011
#' @keywords data
NULL

#' @importFrom graphics lines plot
#' @importFrom stats, median, quantile, sd
#' @importFrom IRanges, IRanges, as.data.frame, sort, setdiff, reduce

