#' Scales-Based Descriptors with 20+ classes of Molecular Descriptors
#'
#' Scales-Based Descriptors with 20+ classes of Molecular Descriptors
#'
#' This function calculates the scales-based descriptors with
#' molecular descriptors sets calculated by
#' Dragon, Discovery Studio and MOE.
#' Users could specify which molecular descriptors to select from one of these
#' deseriptor sets by specify the numerical or character index of the
#' molecular descriptors in the descriptor set.
#'
#' @param x A character vector, as the input protein sequence.
#' @param propmat The matrix containing the descriptor set for the amino acids,
#'        which could be chosen from
#'        \code{AAMOE2D}, \code{AAMOE3D}, \code{AACPSA},
#'        \code{AADescAll}, \code{AA2DACOR}, \code{AA3DMoRSE},
#'        \code{AAACF}, \code{AABurden}, \code{AAConn},
#'        \code{AAConst}, \code{AAEdgeAdj}, \code{AAEigIdx},
#'        \code{AAFGC}, \code{AAGeom}, \code{AAGETAWAY},
#'        \code{AAInfo}, \code{AAMolProp}, \code{AARandic},
#'        \code{AARDF}, \code{AATopo}, \code{AATopoChg},
#'        \code{AAWalk}, \code{AAWHIM}.
#' @param index Integer vector or character vector. Specify which molecular descriptors
#'        to select from one of these deseriptor sets by specify the
#'        numerical or character index of the molecular descriptors in the descriptor set.
#'        Default is \code{NULL}, means selecting all the molecular descriptors
#'        in this descriptor set.
#' @param pc Integer. The maximum dimension of the space which the data
#'        are to be represented in.
#'        Must be no greater than the number of AA properties provided.
#' @param lag The lag parameter. Must be less than the amino acids.
#' @param scale Logical. Should we auto-scale the property matrix
#'        (\code{propmat}) before doing MDS? Default is \code{TRUE}.
#' @param silent Logical. Whether we print the standard deviation,
#'        proportion of variance and the cumulative proportion of
#'        the selected principal components or not.
#'        Default is \code{TRUE}.
#'
#' @return A length \code{lag * p^2} named vector,
#'         \code{p} is the number of scales selected.
#'
#' @keywords extract scales PCM
#'
#' @aliases extractDescScales
#'
#' @author Nan Xiao <\url{http://nanx.me}>
#'
#' @export extractDescScales
#'
#' @examples
#' x = readFASTA(system.file('protseq/P00750.fasta', package = 'protr'))[[1]]
#' descscales = extractDescScales(x, propmat = 'AATopo', index = c(37:41, 43:47),
#'                                pc = 5, lag = 7, silent = FALSE)

extractDescScales = function (x, propmat, index = NULL, pc, lag, scale = TRUE, silent = TRUE) {

  propmat = get(propmat)
  if (!is.null(index)) propmat = propmat[, index]

  result = extractScales(x = x, propmat = propmat, pc = pc, lag = lag, scale = scale, silent = silent)

  return(result)

}
