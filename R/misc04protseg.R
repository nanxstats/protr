#' Protein Sequence Segmentation
#'
#' Protein Sequence Segmentation
#'
#' This function extracts the segmentations from the protein sequence.
#'
#' @param x A character vector, as the input protein sequence.
#'
#' @param aa A character, the amino acid type. One of
#'           \code{'A'}, \code{'R'}, \code{'N'}, \code{'D'},
#'           \code{'C'}, \code{'E'}, \code{'Q'}, \code{'G'},
#'           \code{'H'}, \code{'I'}, \code{'L'}, \code{'K'},
#'           \code{'M'}, \code{'F'}, \code{'P'}, \code{'S'},
#'           \code{'T'}, \code{'W'}, \code{'Y'}, \code{'V'}.
#'
#' @param k A positive integer, specifys the window size (half of the window),
#'          default is 7.
#'
#' @return A named list, each component contains one of the
#'         segmentations (a character string), names of the list components
#'         are the positions of the specified amino acid in the sequence.
#'
#' @keywords segmentation
#'
#' @aliases protseg
#'
#' @author Nan Xiao <\url{http://nanx.me}>
#'
#' @export protseg
#'
#' @examples
#' x = readFASTA(system.file('protseq/P00750.fasta', package = 'protr'))[[1]]
#' protseg(x, aa = 'R', k = 5)

protseg = function (x, aa = c('A', 'R', 'N', 'D', 'C',
                              'E', 'Q', 'G', 'H', 'I',
                              'L', 'K', 'M', 'F', 'P',
                              'S', 'T', 'W', 'Y', 'V'), k = 7) {

  aa = match.arg(aa)

  xSplitted = strsplit(x, split = '')[[1]]
  n = nchar(x)

  CenterIdx = which(xSplitted == aa)

  if (length(CenterIdx) < 0.5) stop(paste('Did not find AA', aa, 'in the sequence'))

  CenterIdx = CenterIdx[(CenterIdx - k) > 0.5 & (CenterIdx + k) < (n + 0.5)]

  if (length(CenterIdx) < 0.5) stop(paste('Segmentation does not exist for AA', aa, 'and step', k))

  Segments = vector('list', length(CenterIdx))

  for (i in 1:length(CenterIdx)) {
    Segments[[i]] = paste(xSplitted[(CenterIdx[i] - k):(CenterIdx[i] + k)],
                          collapse = '')
  }

  names(Segments) = as.character(CenterIdx)

  return(Segments)

}
