#' Large-Scale Parallelled Protein Sequence Similarity Calculation based on Smith-Waterman Local Alignment
#'
#' Large-Scale Parallelled Protein Sequence Similarity Calculation based on Smith-Waterman Local Alignment
#'
#' This function implemented the parallelled version for calculating 
#' protein sequence similarity based on Smith-Waterman local alignment.
#' 
#' @param x A character vector, as the input protein sequence.
#' @param propmat A matrix containing the properties for the amino acids. 
#'        Each row represent one amino acid type, each column represents one property.
#'        Note that the one-letter row names must be provided for we need them to seek 
#'        the properties for each AA type.
#' @param pc Integer. Use the first pc principal components as the scales.
#'        Must be no greater than the number of AA properties provided.
#' @param lag The lag parameter. Must be less than the amino acids.
#' @param scale Logical. Should we auto-scale the property matrix 
#'        (\code{propmat}) before PCA? Default is \code{TRUE}.
#' @param silent Logical. Whether we print the standard deviation, 
#'        proportion of variance and the cumulative proportion of 
#'        the selected principal components or not.
#'        Default is \code{TRUE}.
#'
#' @return A length \code{lag * p^2} named vector, 
#'         \code{p} is the number of scales (principal components) selected.
#' 
#' @note Note
#' 
#' @keywords smith-waterman sequence local alignment parallel similarity parSeqSim
#'
#' @aliases parSeqSim
#' 
#' @author Xiao Nan <\url{http://www.road2stat.com}>
#' 
#' @seealso See \code{\link{parGOSim}} for paralleled protein similarity
#' calculation based on Gene Ontology (GO) similarity.
#' 
#' @export parSeqSim
#' 
#' @examples
#' x = readFASTA(system.file('protseq/P00750.fasta', package = 'protr'))[[1]]
#' AAidxmat = t(na.omit(as.matrix(AAindex[, 7:26])))
#' scales = parSeqSim(x, propmat = AAidxmat, pc = 5, lag = 7, silent = FALSE)
#' 

parSeqSim = function (x, propmat, pc, lag, scale = TRUE, silent = TRUE) {

  return(NULL)
  
}
