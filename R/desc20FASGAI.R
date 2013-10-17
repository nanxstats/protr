#' Factor Analysis Scales of Generalized AA Information (FASGAI) Descriptor
#'
#' Factor Analysis Scales of Generalized AA Information (FASGAI) Descriptor
#' 
#' This function calculates the Factor Analysis Scales 
#' of Generalized AA Information (FASGAI) descriptors.
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
#' @keywords extract FASGAI extractFASGAI PCM
#'
#' @aliases extractFASGAI
#' 
#' @author Xiao Nan <\url{http://www.road2stat.com}>
#' 
#' @export extractFASGAI
#' 
#' @references
#' TBD
#' 
#' @examples
#' x = readFASTA(system.file('protseq/P00750.fasta', package = 'protr'))[[1]]
#' extractFASGAI(x)
#' 

extractFASGAI = function (x) {
  
  if (protcheck(x) == FALSE) stop('x has unrecognized amino acid type')
  
  return(NULL)
  
}

