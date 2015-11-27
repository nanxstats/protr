#' Scales-Based Descriptors derived by Multidimensional Scaling
#'
#' Scales-Based Descriptors derived by Multidimensional Scaling
#'
#' This function calculates scales-based descriptors
#' derived by Multidimensional Scaling (MDS).
#' Users could provide customized amino acid property matrices.
#'
#' @param x A character vector, as the input protein sequence.
#' @param propmat A matrix containing the properties for the amino acids.
#'        Each row represent one amino acid type, each column represents one property.
#'        Note that the one-letter row names must be provided for we need them to seek
#'        the properties for each AA type.
#' @param k Integer. The maximum dimension of the space which the data
#'        are to be represented in.
#'        Must be no greater than the number of AA properties provided.
#' @param lag The lag parameter. Must be less than the amino acids.
#' @param scale Logical. Should we auto-scale the property matrix
#'        (\code{propmat}) before doing MDS? Default is \code{TRUE}.
#' @param silent Logical. Whether we print the \code{k} eigenvalues
#'        computed during the scaling process or not.
#'        Default is \code{TRUE}.
#'
#' @return A length \code{lag * p^2} named vector,
#'         \code{p} is the number of scales (dimensionality) selected.
#'
#' @keywords extract MDS PCM
#'
#' @aliases extractMDSScales
#'
#' @author Nan Xiao <\url{http://nanx.me}>
#'
#' @seealso See \code{\link{extractScales}} for scales-based
#' descriptors derived by Principal Components Analysis.
#'
#' @export extractMDSScales
#'
#' @references
#' Venkatarajan, M. S., & Braun, W. (2001).
#' New quantitative descriptors of amino acids based on multidimensional scaling
#' of a large number of physical-chemical properties.
#' Molecular modeling annual, 7(12), 445--453.
#'
#' @examples
#' x = readFASTA(system.file('protseq/P00750.fasta', package = 'protr'))[[1]]
#' data(AATopo)
#' tprops = AATopo[, c(37:41, 43:47)]  # select a set of topological descriptors
#' mds = extractMDSScales(x, propmat = tprops, k = 5, lag = 7, silent = FALSE)

extractMDSScales = function (x, propmat, k, lag, scale = TRUE, silent = TRUE) {

  if (protcheck(x) == FALSE) stop('x has unrecognized amino acid type')

  k = min(k, ncol(propmat) - 1, nrow(propmat) - 1)

  if (scale) propmat = scale(propmat)

  d = dist(propmat)  # euclidean distances between the rows
  mds = cmdscale(d, k = k, eig = TRUE)

  accmat = matrix(0, k, nchar(x))
  x.split = strsplit(x, '')[[1]]

  for (i in 1:nchar(x)) {
    accmat[, i] = mds$points[x.split[i], 1:k]
  }

  result = acc(accmat, lag)

  if (!silent) {
    cat('Eigenvalues computed during the scaling process:\n')
    print(mds$eig)
  }

  return(result)

}
