#' Conjoint Triad Descriptor (with Customized Amino Acid Classification Support)
#'
#' Conjoint Triad Descriptor (with Customized Amino Acid Classification Support)
#'
#' This function calculates the Conjoint Triad descriptor, with customized
#' amino acid classification support.
#'
#' @param x A character vector, as the input protein sequence.
#' @param aaclass A list containing the customized amino acid
#' classification. See example below.
#'
#' @return A length \code{k^3} named vector, where \code{k} is the number
#' of customized classes of the amino acids.
#'
#' @keywords extract Conjoint Triad
#'
#' @aliases extractCTriadClass
#'
#' @author Nan Xiao <\url{http://nanx.me}>
#'
#' @export extractCTriadClass
#'
#' @note For this descriptor type, users need to intelligently evaluate
#' the underlying details of the descriptors provided, instead of using
#' this function with their data blindly. It would be wise to use some
#' negative and positive control comparisons where relevant to help guide
#' interpretation of the results.
#'
#' @references
#' J.W. Shen, J. Zhang, X.M. Luo, W.L. Zhu,
#' K.Q. Yu, K.X. Chen, Y.X. Li, H.L. Jiang.
#' Predicting Protein-protein Interactions Based Only on Sequences Information.
#' \emph{Proceedings of the National Academy of Sciences}. 007, 104, 4337--4341.
#'
#' @examples
#' x = readFASTA(system.file('protseq/P00750.fasta', package = 'protr'))[[1]]
#'
#' # using customized amino acid classification (normalized van der Waals volume)
#' newclass = list(c('G', 'A', 'S', 'T', 'P', 'D', 'C'),
#'                 c('N', 'V', 'E', 'Q', 'I', 'L'),
#'                 c('M', 'H', 'K', 'F', 'R', 'Y', 'W'))
#'
#' extractCTriadClass(x, aaclass = newclass)

extractCTriadClass = function (x, aaclass) {

  if (protcheck(x) == FALSE) stop('x has unrecognized amino acid type')

  nclass = length(aaclass)

  vspace = expand.grid(1L:nclass, 1L:nclass, 1L:nclass)
  veclen = nrow(vspace)
  CTDict = vector('list', veclen)

  for (i in 1L:veclen) {
    tmp = as.vector(outer(aaclass[[vspace[i, 1L]]],
                          aaclass[[vspace[i, 2L]]], paste, sep = ''))
    CTDict[[i]] = as.vector(outer(tmp,
                                  aaclass[[vspace[i, 3L]]], paste, sep = ''))
  }

  CTDict = unlist(CTDict)

  classlen = sapply(aaclass, length)
  CTIndexeach = rep(NA, veclen)

  for (i in 1L:veclen) CTIndexeach[i] = classlen[vspace[i, 1]] *
    classlen[vspace[i, 2]] * classlen[vspace[i, 3]]

  CTIndex = rep(1L:veclen, CTIndexeach)

  xSplitted = strsplit(x, split = '')[[1L]]
  n  = nchar(x)
  CTAll = summary(factor(paste(paste(xSplitted[-c(n, n-1L)],
                                     xSplitted[-c(1L, n)], sep = ''),
                               xSplitted[-c(1L, 2L)], sep = ''),
                         levels = CTDict), maxsum = length(CTDict)+1L)

  MatchedIndex = which(CTAll != 0L)
  MatchedNames = names(CTAll[MatchedIndex])
  MatchedTimes = as.integer(CTAll[MatchedIndex])
  CTAll = rep(MatchedNames, times = MatchedTimes)

  CT = rep(0L, veclen)

  for (i in 1L:length(MatchedNames)) {
    idx = CTIndex[which(CTDict == MatchedNames[i])]
    CT[idx] = CT[idx] + MatchedTimes[i]
  }

  CT = (CT - min(CT))/max(CT)

  names(CT) = paste0('VS', paste0(paste0(vspace[, 1], vspace[, 2]),
                                  vspace[, 3]))

  return(CT)

}
