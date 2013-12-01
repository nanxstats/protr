#' Parallellized Protein Sequence Similarity Calculation based on Sequence Alignment
#'
#' Parallellized Protein Sequence Similarity Calculation based on Sequence Alignment
#'
#' This function implemented the parallellized version for calculating 
#' protein sequence similarity based on sequence alignment.
#' 
#' @param x A character vector, as the input protein sequence.
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
#' @author Nan Xiao <\url{http://www.road2stat.com}>
#' 
#' @seealso See \code{\link{parGOSim}} for paralleled protein similarity
#' calculation based on Gene Ontology (GO) similarity.
#' 
#' @export parSeqSim
#' 
#' @examples
#' parSeqSim(x)
#' 

parSeqSim = function (x) {

  return(NULL)
  
}
