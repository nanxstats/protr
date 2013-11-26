#' Parallellized Protein Sequence Similarity Calculation based on Gene Ontology (GO) Similarity
#'
#' Parallellized Protein Sequence Similarity Calculation based on Gene Ontology (GO) Similarity
#'
#' This function implemented the parallellized version for calculating 
#' protein sequence similarity based on Gene Ontology (GO) similarity.
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
#' @keywords GO Gene Ontology parallel similarity parGOSim
#'
#' @aliases parGOSim
#' 
#' @author Nan Xiao <\url{http://www.road2stat.com}>
#' 
#' @seealso See \code{\link{parSeqSim}} for paralleled protein similarity
#' calculation based on Smith-Waterman local alignment.
#' 
#' @export parGOSim
#' 
#' @examples
#' x = readFASTA(system.file('protseq/P00750.fasta', package = 'protr'))[[1]]
#' scales = parGOSim(x, propmat = AAidxmat, pc = 5, lag = 7, silent = FALSE)
#' 

parGOSim = function (x, propmat, pc, lag, scale = TRUE, silent = TRUE) {
  
  return(NULL)
  
}
