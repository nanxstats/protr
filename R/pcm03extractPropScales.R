#' Generalized AA-Properties Based Scales Descriptors
#'
#' Generalized AA-Properties Based Scales Descriptors
#'
#' This function calculates the generalized amino acid properties based scales descriptors.
#' Users could provide customized amino acid descriptor matrices.
#' 
#' @param x A character vector, as the input protein sequence.
#' @param descmat A matrix containing the properties for the amino acids. 
#'        Each row represent one amino acid type, each column represents one property.
#'        Note that the one-letter row names must be provided for we need them to seek 
#'        the properties for each AA type.
#' @param pc Integer. Use the first pc principal components as the scales.
#'        Must be no greater than the number of AA properties provided.
#' @param lag The lag parameter. Must be less than the amino acids.
#' @param scale Logical. Should we auto-scale the property matrix 
#'        (\code{descmat}) before PCA? Default is \code{TRUE}.
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
#' @keywords extract descriptor scales extractPropScales PCM
#'
#' @aliases extractPropScales
#' 
#' @author Xiao Nan <\url{http://www.road2stat.com}>
#' 
#' @seealso See \code{\link{extractScales}} for generalized scales-based descriptors.
#' 
#' @export extractPropScales
#' 
#' @examples
#' x = readFASTA(system.file('protseq/P00750.fasta', package = 'protr'))[[1]]
#' # data(AAindex)
#' # AAidxmat = t(na.omit(as.matrix(AAindex[, 7:26])))
#' # scales = extractPropScales(x, descmat = AAidxmat, pc = 5, lag = 7, silent = FALSE)
#' 

extractPropScales = function (x, descmat, pc, lag, scale = TRUE, silent = TRUE) {
  
  if (protcheck(x) == FALSE) stop('x has unrecognized amino acid types')
  
  pc = min(pc, ncol(descmat), nrow(descmat))
  
  prop.pr = prcomp(descmat, scale = scale)
  prop.pred = predict(prop.pr)
  
  accmat = matrix(0, pc, nchar(x))
  x.split = strsplit(x, '')[[1]]
  
  for (i in 1:nchar(x)) {
    accmat[, i] = prop.pred[x.split[i], 1:pc]
  }
  
  result = acc(accmat, lag)
  
  if (!silent) {
    cat('Summary of the first', pc,'principal components:\n')
    print(summary(prop.pr)$importance[, 1:pc])
  }
  
  return(result)
  
}
