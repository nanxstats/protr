#' T-scales Descriptor
#'
#' T-scales Descriptor
#' 
#' This function calculates the T-scales descriptor
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
#' @keywords extract T-scales extractTscales PCM
#'
#' @aliases extractTscales
#' 
#' @author Xiao Nan <\url{http://www.road2stat.com}>
#' 
#' @seealso See \code{\link{extractZscales}}, \code{\link{extractSTscales}} 
#' and \code{\link{extractVHSE}} for Z-scales, ST-scales and VHSE descriptors.
#' See \code{\link{extractScales}} for 
#' generalized scales-based descriptors.
#' 
#' @export extractTscales
#' 
#' @references
#' TBD
#' 
#' @examples
#' x = readFASTA(system.file('protseq/P00750.fasta', package = 'protr'))[[1]]
#' extractTscales(x)
#' 

extractTscales = function (x) {
  
  if (protcheck(x) == FALSE) stop('x has unrecognized amino acid type')
  
  return(NULL)
  
}

