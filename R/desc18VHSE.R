#' Vectors of Hydrophobic, Steric, and Electronic (VHSE) Descriptor
#'
#' Vectors of Hydrophobic, Steric, and Electronic (VHSE) Descriptor
#' 
#' This function calculates the Vectors of Hydrophobic, Steric, 
#' and Electronic (VHSE) descriptor.
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
#' @keywords extract VHSE extractVHSE PCM
#'
#' @aliases extractVHSE
#' 
#' @author Xiao Nan <\url{http://www.road2stat.com}>
#' 
#' @seealso See \code{\link{extractZscales}}, \code{\link{extractTscales}} 
#' and \code{\link{extractSTscales}} for Z-scales, T-scales and ST-scales descriptors.
#' See \code{\link{extractScales}} for 
#' generalized scales-based descriptors.
#' 
#' @export extractVHSE
#' 
#' @references
#' TBD
#' 
#' @examples
#' x = readFASTA(system.file('protseq/P00750.fasta', package = 'protr'))[[1]]
#' extractVHSE(x)
#' 

extractVHSE = function (x) {
  
  if (protcheck(x) == FALSE) stop('x has unrecognized amino acid type')
  
  return(NULL)
  
}

