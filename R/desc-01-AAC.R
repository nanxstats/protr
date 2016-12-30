#' Amino Acid Composition Descriptor
#'
#' Amino Acid Composition Descriptor
#'
#' This function calculates the Amino Acid Composition descriptor (Dim: 20).
#'
#' @param x A character vector, as the input protein sequence.
#'
#' @return A length 20 named vector
#'
#' @keywords extract Composition
#'
#' @aliases extractAAC
#'
#' @author Nan Xiao <\url{http://nanx.me}>
#'
#' @seealso See \code{\link{extractDC}} and \code{\link{extractTC}}
#'          for Dipeptide Composition and Tripeptide Composition descriptors.
#'
#' @export extractAAC
#'
#' @references
#' M. Bhasin, G. P. S. Raghava.
#' Classification of Nuclear Receptors Based on
#' Amino Acid Composition and Dipeptide Composition.
#' \emph{Journal of Biological Chemistry}, 2004, 279, 23262.
#'
#' @examples
#' x = readFASTA(system.file('protseq/P00750.fasta', package = 'protr'))[[1]]
#' extractAAC(x)

extractAAC = function (x) {

  if (protcheck(x) == FALSE) stop('x has unrecognized amino acid type')

  # The 20 Amino Acid Abbrevation Dictionary is from
  # http://en.wikipedia.org/wiki/Amino_acid#Table_of_standard_amino_acid_abbreviations_and_properties

  AADict = c('A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I',
             'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V')

  AAC = summary(factor(strsplit(x, split = '')[[1]], levels = AADict),
                maxsum = 21)/nchar(x)

  return(AAC)

}
