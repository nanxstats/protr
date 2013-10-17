#' Generalized Scales-Based Descriptors
#'
#' Generalized Scales-Based Descriptors
#'
#' This function calculates the generalized scales-based descriptors.
#' This function implements the core computation needed for 
#' Z-scales T-scales, ST-scales and VHSE descriptors 
#' in the protr package.
#' 
#' @param x A character vector, as the input protein sequence.
#' @param propmat A matrix containing the properties for the amino acids. 
#'        Each row represent one amino acid type, each column represents one property.
#'        Note that the one-letter row names must be provided for we need them to seek 
#'        the properties for each AA type.
#' @param pc Integer. Use the first pc principle components as the scales.
#'        Must be no greater than the number of AA properties provided.
#' @param lag The lag parameter. Must be less than the amino acids.
#' @param scale Logical. Should we auto-scale the property matrix 
#'        (\code{propmat}) before PCA?
#'
#' @return A length \code{lag * p^2} named vector, 
#'         \code{p} is the number of scales (principle components) selected.
#' 
#' @note Note
#' 
#' @keywords extract scales extractScales PCM
#'
#' @aliases extractScales
#' 
#' @author Xiao Nan <\url{http://www.road2stat.com}>
#' 
#' @seealso See \code{\link{extractZscales}}, \code{\link{extractTscales}}, 
#' \code{\link{extractSTscales}} and \code{\link{extractVHSE}}
#' for Z-scales, T-scales, ST-scales and VHSE descriptors.
#' 
#' @export extractScales
#' 
#' @examples
#' x = readFASTA(system.file('protseq/P00750.fasta', package = 'protr'))[[1]]
#' extractScales(x)
#' 

extractScales = function (x, propmat, pc, lag, scale = TRUE) {
  
  if (protcheck(x) == FALSE) stop('x has unrecognized amino acid type')
  
  if (pc > ncol(x)) stop('PCs must be less than variables')
  
  
  
  print(paste('First', pc, 'principle components explained', var, '% of the total variance'))
  
  return(NULL)
  
}

# x = readFASTA(system.file('protseq/P00750.fasta', package = 'protr'))[[1]]
# pc = 5
# lag = 5
# 
# require(protr)
# data(AAindex)
# propmat = as.matrix(AAindex[, 7:26])
# propmat = na.omit(propmat)
# propmat = t(propmat)
# 
# if (scale = TRUE) propmat = scale(propmat)
#
# prop.pr = prcomp(propmat)
# prop.pred = predict(prop.pr)
# prop.pred[, 1:pc]
# 
# accmat = matrix(0, pc, nchar(x))
# x1 = strsplit(x, '')[[1]]
# 
# for (i in 1:nchar(x)) {
#   accmat[, i] = prop.pred[x1[i], 1:pc]
# }
# 
# result = acc(accmat, lag)
