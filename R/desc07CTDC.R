#' CTD Descriptors - Composition
#'
#' CTD Descriptors - Composition
#'
#' This function calculates the Composition descriptor of the
#' CTD descriptors (Dim: 21).
#'
#' @param x A character vector, as the input protein sequence.
#'
#' @return A length 21 named vector
#'
#' @keywords extract CTD Composition
#'
#' @aliases extractCTDC
#'
#' @author Nan Xiao <\url{http://nanx.me}>
#'
#' @seealso See \code{\link{extractCTDT}} and \code{\link{extractCTDD}}
#'          for Transition and Distribution of the CTD descriptors.
#'
#' @export extractCTDC
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
#' extractCTDC(x)

extractCTDC = function (x) {

  if (protcheck(x) == FALSE) stop('x has unrecognized amino acid type')

  group1 = list(hydrophobicity  = c('R', 'K', 'E', 'D', 'Q', 'N'),
                normwaalsvolume = c('G', 'A', 'S', 'T', 'P', 'D', 'C'),
                polarity        = c('L', 'I', 'F', 'W', 'C', 'M', 'V', 'Y'),
                polarizability  = c('G', 'A', 'S', 'D', 'T'),
                charge          = c('K', 'R'),
                secondarystruct = c('E', 'A', 'L', 'M', 'Q', 'K', 'R', 'H'),
                solventaccess   = c('A', 'L', 'F', 'C', 'G', 'I', 'V', 'W'))

  group2 = list(hydrophobicity  = c('G', 'A', 'S', 'T', 'P', 'H', 'Y'),
                normwaalsvolume = c('N', 'V', 'E', 'Q', 'I', 'L'),
                polarity        = c('P', 'A', 'T', 'G', 'S'),
                polarizability  = c('C', 'P', 'N', 'V', 'E', 'Q', 'I', 'L'),
                charge          = c('A', 'N', 'C', 'Q', 'G', 'H', 'I', 'L',
                                    'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'),
                secondarystruct = c('V', 'I', 'Y', 'C', 'W', 'F', 'T'),
                solventaccess   = c('R', 'K', 'Q', 'E', 'N', 'D'))

  group3 = list(hydrophobicity  = c('C', 'L', 'V', 'I', 'M', 'F', 'W'),
                normwaalsvolume = c('M', 'H', 'K', 'F', 'R', 'Y', 'W'),
                polarity        = c('H', 'Q', 'R', 'K', 'N', 'E', 'D'),
                polarizability  = c('K', 'M', 'H', 'F', 'R', 'Y', 'W'),
                charge          = c('D', 'E'),
                secondarystruct = c('G', 'N', 'P', 'S', 'D'),
                solventaccess   = c('M', 'S', 'P', 'T', 'H', 'Y'))

  xSplitted = strsplit(x, split = '')[[1]]
  n  = nchar(x)

  # Get groups for each property & each amino acid

  g1 = lapply(group1, function(g) length(which(xSplitted %in% g)))
  names(g1) = paste(names(g1), 'Group1', sep = '.')
  g2 = lapply(group2, function(g) length(which(xSplitted %in% g)))
  names(g2) = paste(names(g2), 'Group2', sep = '.')
  g3 = lapply(group3, function(g) length(which(xSplitted %in% g)))
  names(g3) = paste(names(g3), 'Group3', sep = '.')
  CTDC = unlist(c(g1, g2, g3))/n
  # This reorders the groups from
  # p1.g1, p2.g1, p3.g1 ...
  # to
  # p1.g1, p1.g2, p1.g3 ...
  # This reordering is not really important
  # just to make the results look pretty
  ids = unlist(lapply(1:7, function(x) x + c(0, 7, 14)))

  return(CTDC[ids])

}
