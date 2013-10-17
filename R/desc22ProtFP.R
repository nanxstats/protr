#' Protein Fingerprint (ProtFP)
#'
#' Protein Fingerprint (ProtFP)
#' 
#' This function calculates the Protein Fingerprint (ProtFP)
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
#' @keywords extract fingerprint ProtFP extractProtFP PCM
#'
#' @aliases extractProtFP
#' 
#' @author Xiao Nan <\url{http://www.road2stat.com}>
#' 
#' @export extractProtFP
#' 
#' @references
#' TBD
#' 
#' @examples
#' x = readFASTA(system.file('protseq/P00750.fasta', package = 'protr'))[[1]]
#' extractProtFP(x)
#' 

extractProtFP = function (x) {
  
  if (protcheck(x) == FALSE) stop('x has unrecognized amino acid type')
  
  return(NULL)
  
}

