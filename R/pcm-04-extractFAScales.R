#' Scales-Based Descriptors derived by Factor Analysis
#'
#' Scales-Based Descriptors derived by Factor Analysis
#'
#' This function calculates scales-based descriptors
#' derived by Factor Analysis (FA).
#' Users could provide customized amino acid property matrices.
#'
#' @param x A character vector, as the input protein sequence.
#' @param propmat A matrix containing the properties for the amino acids.
#'        Each row represent one amino acid type, each column represents one property.
#'        Note that the one-letter row names must be provided for we need them to seek
#'        the properties for each AA type.
#' @param factors Integer. The number of factors to be fitted.
#'        Must be no greater than the number of AA properties provided.
#' @param scores Type of scores to produce. The default is \code{"regression"},
#'        which gives Thompson's scores, \code{"Bartlett"} given Bartlett's weighted
#'        least-squares scores.
#' @param lag The lag parameter. Must be less than the amino acids number
#'            in the protein sequence.
#' @param scale Logical. Should we auto-scale the property matrix
#'        (\code{propmat}) before doing Factor Analysis? Default is \code{TRUE}.
#' @param silent Logical. Whether we print the SS loadings,
#'        proportion of variance and the cumulative proportion of
#'        the selected factors or not.
#'        Default is \code{TRUE}.
#'
#' @return A length \code{lag * p^2} named vector,
#'         \code{p} is the number of scales (factors) selected.
#'
#' @keywords extract Factor PCM
#'
#' @aliases extractFAScales
#'
#' @author Nan Xiao <\url{http://nanx.me}>
#'
#' @export extractFAScales
#'
#' @references
#' Atchley, W. R., Zhao, J., Fernandes, A. D., & Druke, T. (2005).
#' Solving the protein sequence metric problem.
#' Proceedings of the National Academy of Sciences of the United States of America,
#' 102(18), 6395-6400.
#'
#' @examples
#' x = readFASTA(system.file('protseq/P00750.fasta', package = 'protr'))[[1]]
#' data(AATopo)
#' tprops = AATopo[, c(37:41, 43:47)]  # select a set of topological descriptors
#' fa = extractFAScales(x, propmat = tprops, factors = 5, lag = 7, silent = FALSE)

extractFAScales = function (x, propmat, factors, scores = 'regression', lag, scale = TRUE, silent = TRUE) {

  if (protcheck(x) == FALSE) stop('x has unrecognized amino acid type')

  factors = min(factors, ncol(propmat), nrow(propmat))

  if (scale) propmat = scale(propmat)

  prop.fa = factanal(propmat, factors = factors, scores = scores)
  prop.scores = prop.fa$scores

  accmat = matrix(0, factors, nchar(x))
  x.split = strsplit(x, '')[[1]]

  for (i in 1:nchar(x)) {
    accmat[, i] = prop.scores[x.split[i], 1:factors]
  }

  result = acc(accmat, lag)

  if (!silent) {
    cat('Summary of the factor analysis result:\n')
    print(prop.fa)
  }

  return(result)

}
