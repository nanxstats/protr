#' Parallellized Protein Sequence Similarity Calculation based on Gene Ontology (GO) Similarity
#'
#' Parallellized Protein Sequence Similarity Calculation based on Gene Ontology (GO) Similarity
#'
#' This function implemented the parallellized version for calculating 
#' protein sequence similarity based on Gene Ontology (GO) similarity.
#' 
#' @param x A character vector, as the input protein sequence.
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
#' parGOSim(x)
#' 

parGOSim = function (x) {
  
  return(NULL)
  
}
