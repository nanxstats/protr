.protcheckgap = function (x) {

  AADict = c('A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I',
             'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V',
             '-')

  return(all(strsplit(x, split = '')[[1]] %in% AADict))

}

#' Scales-Based Descriptors derived by Principal Components Analysis
#' (with Gap Support)
#'
#' Scales-Based Descriptors derived by Principal Components Analysis
#' (with Gap Support)
#'
#' This function calculates scales-based descriptors
#' derived by Principal Components Analysis (PCA), with gap support.
#' Users could provide customized amino acid property matrices.
#' This function implements the core computation procedure needed for
#' the scales-based descriptors derived by AA-Properties (AAindex)
#' and scales-based descriptors derived by 20+ classes of 2D and 3D
#' molecular descriptors (Topological, WHIM, VHSE, etc.) in the protr package.
#'
#' @param x A character vector, as the input protein sequence.
#'        Use '\code{-}' to represent gaps in the sequence.
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
#' @keywords extract scales PCA PCM gap
#'
#' @aliases extractScalesGap
#'
#' @author Nan Xiao <\url{http://nanx.me}>
#'
#' @seealso See \code{\link{extractProtFPGap}} for amino acid property based
#' scales descriptors (protein fingerprint) with gap support.
#'
#' @export extractScalesGap
#'
#' @examples
#' # amino acid sequence with gaps
#' x = readFASTA(system.file('protseq/align.fasta', package = 'protr'))$`IXI_235`
#' data(AAindex)
#' AAidxmat = t(na.omit(as.matrix(AAindex[, 7:26])))
#' scales = extractScalesGap(x, propmat = AAidxmat, pc = 5, lag = 7, silent = FALSE)

extractScalesGap = function (x, propmat, pc, lag,
                             scale = TRUE, silent = TRUE) {

  if (.protcheckgap(x) == FALSE) stop('x has unrecognized amino acid types. Note: use "-" to represent gaps.')

  gapmat = t(matrix(rep(0L, ncol(propmat))))
  row.names(gapmat) = '-'
  propmat = rbind(propmat, gapmat)

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
