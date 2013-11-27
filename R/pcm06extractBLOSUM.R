#' Generalized BLOSUM Matrix-Derived Descriptors
#' 
#' Generalized BLOSUM Matrix-Derived Descriptors
#' 
#' This function calculates the generalized BLOSUM matrix-derived descriptors
#' (Dim: \code{20 + (n * lambda)}, 
#' \code{n} is the number of properties selected, default is 80).
#' 
#' @param x A character vector, as the input protein sequence.
#'
#' @return A length \code{20 + n * lambda} named vector, 
#'         \code{n} is the number of properties selected.
#' 
#' @note Note
#' 
#' @keywords extract BLOSUM extractBLOSUM PCM
#'
#' @aliases extractBLOSUM
#' 
#' @author Nan Xiao <\url{http://www.road2stat.com}>
#' 
#' @export extractBLOSUM
#' 
#' @references
#' TBD
#' 
#' @examples
#' x = readFASTA(system.file('protseq/P00750.fasta', package = 'protr'))[[1]]
#' # extractBLOSUM(x)
#' 

extractBLOSUM = function (x) {
  
  if (protcheck(x) == FALSE) stop('x has unrecognized amino acid type')
  
  return(NULL)
  
}

