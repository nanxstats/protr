#' Scales-Based Descriptors derived by Principal Components Analysis
#'
#' Scales-Based Descriptors derived by Principal Components Analysis
#'
#' This function calculates scales-based descriptors
#' derived by Principal Components Analysis (PCA).
#' Users could provide customized amino acid property matrices.
#' This function implements the core computation procedure needed for
#' the scales-based descriptors derived by AA-Properties (AAindex)
#' and scales-based descriptors derived by 20+ classes of 2D and 3D
#' molecular descriptors (Topological, WHIM, VHSE, etc.) in the protr package.
#'
#' @param x A character vector, as the input protein sequence.
#' @param propmat A matrix containing the properties for the amino acids.
#'        Each row represent one amino acid type, each column represents
#'        one property. Note that the one-letter row names must be provided
#'        for we need them to seek the properties for each AA type.
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
#' @keywords extract scales PCA PCM
#'
#' @aliases extractScales
#'
#' @author Nan Xiao <\url{http://nanx.me}>
#'
#' @seealso See \code{\link{extractDescScales}} scales descriptors based on
#' 20+ classes of molecular descriptors, and \code{\link{extractProtFP}}
#' for amino acid property based scales descriptors (protein fingerprint).
#'
#' @export extractScales
#'
#' @examples
#' x = readFASTA(system.file('protseq/P00750.fasta', package = 'protr'))[[1]]
#' data(AAindex)
#' AAidxmat = t(na.omit(as.matrix(AAindex[, 7:26])))
#' scales = extractScales(x, propmat = AAidxmat, pc = 5, lag = 7, silent = FALSE)

extractScales = function (x, propmat, pc, lag, scale = TRUE, silent = TRUE) {

  if (protcheck(x) == FALSE) stop('x has unrecognized amino acid types')

  pc = min(pc, ncol(propmat), nrow(propmat))

  prop.pr = prcomp(propmat, scale = scale)
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
