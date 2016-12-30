#' CTD Descriptors - Composition (with Customized Amino Acid Classification Support)
#'
#' CTD Descriptors - Composition (with Customized Amino Acid Classification Support)
#'
#' This function calculates the Composition descriptor of the
#' CTD descriptors, with customized amino acid classification support.
#'
#' @param x A character vector, as the input protein sequence.
#' @param aagroup1 A named list which contains the first group of customized
#' amino acid classification. See example below.
#' @param aagroup2 A named list which contains the second group of customized
#' amino acid classification. See example below.
#' @param aagroup3 A named list which contains the third group of customized
#' amino acid classification. See example below.
#'
#' @return A length \code{k * 3} named vector, \code{k} is the number of
#' amino acid properties used.
#'
#' @keywords extract CTD Composition
#'
#' @aliases extractCTDCClass
#'
#' @author Nan Xiao <\url{http://nanx.me}>
#'
#' @seealso See \code{\link{extractCTDTClass}} and \code{\link{extractCTDDClass}}
#'          for Transition and Distribution of the CTD descriptors with
#'          customized amino acid classification support.
#'
#' @export extractCTDCClass
#'
#' @note For this descriptor type, users need to intelligently evaluate
#' the underlying details of the descriptors provided, instead of using
#' this function with their data blindly. It would be wise to use some
#' negative and positive control comparisons where relevant to help guide
#' interpretation of the results.
#'
#' @references
#' Inna Dubchak, Ilya Muchink, Stephen R. Holbrook and Sung-Hou Kim.
#' Prediction of protein folding class using global description of
#' amino acid sequence. \emph{Proceedings of the National Academy of Sciences}.
#' USA, 1995, 92, 8700-8704.
#'
#' Inna Dubchak, Ilya Muchink, Christopher Mayor, Igor Dralyuk and Sung-Hou Kim.
#' Recognition of a Protein Fold in the Context of the SCOP classification.
#' \emph{Proteins: Structure, Function and Genetics}, 1999, 35, 401-407.
#'
#' @examples
#' x = readFASTA(system.file('protseq/P00750.fasta', package = 'protr'))[[1]]
#'
#' # using five customized amino acid property classification
#' group1 = list(hydrophobicity  = c('R', 'K', 'E', 'D', 'Q', 'N'),
#'               normwaalsvolume = c('G', 'A', 'S', 'T', 'P', 'D', 'C'),
#'               polarizability  = c('G', 'A', 'S', 'D', 'T'),
#'               secondarystruct = c('E', 'A', 'L', 'M', 'Q', 'K', 'R', 'H'),
#'               solventaccess   = c('A', 'L', 'F', 'C', 'G', 'I', 'V', 'W'))
#'
#' group2 = list(hydrophobicity  = c('G', 'A', 'S', 'T', 'P', 'H', 'Y'),
#'               normwaalsvolume = c('N', 'V', 'E', 'Q', 'I', 'L'),
#'               polarizability  = c('C', 'P', 'N', 'V', 'E', 'Q', 'I', 'L'),
#'               secondarystruct = c('V', 'I', 'Y', 'C', 'W', 'F', 'T'),
#'               solventaccess   = c('R', 'K', 'Q', 'E', 'N', 'D'))
#'
#' group3 = list(hydrophobicity  = c('C', 'L', 'V', 'I', 'M', 'F', 'W'),
#'               normwaalsvolume = c('M', 'H', 'K', 'F', 'R', 'Y', 'W'),
#'               polarizability  = c('K', 'M', 'H', 'F', 'R', 'Y', 'W'),
#'               secondarystruct = c('G', 'N', 'P', 'S', 'D'),
#'               solventaccess   = c('M', 'S', 'P', 'T', 'H', 'Y'))
#'
#' extractCTDCClass(x, aagroup1 = group1, aagroup2 = group2, aagroup3 = group3)

extractCTDCClass = function (x, aagroup1, aagroup2, aagroup3) {

  if (protcheck(x) == FALSE) stop('x has unrecognized amino acid type')

  if ((length(aagroup1) != length(aagroup2) |
       length(aagroup1) != length(aagroup3)) |
      (length(aagroup2) != length(aagroup3)))
    stop('The three groups must have the same property numbers')

  xSplitted = strsplit(x, split = '')[[1]]
  n  = nchar(x)

  propnum = length(aagroup1)

  G = vector('list', propnum)
  for (i in 1L:propnum) G[[i]] = rep(NA, n)

  # Get groups for each property & each amino acid

  for (i in 1L:propnum) {
    try(G[[i]][which(xSplitted %in% aagroup1[[i]])] <- 'G1')
    try(G[[i]][which(xSplitted %in% aagroup2[[i]])] <- 'G2')
    try(G[[i]][which(xSplitted %in% aagroup3[[i]])] <- 'G3')
  }

  G = lapply(G, as.factor)

  CTDC = unlist(lapply(G, summary))/n
  names(CTDC) = paste('prop', rep(1L:propnum, each = 3L),
                      '.', names(CTDC), sep = '')

  return(CTDC)

}
