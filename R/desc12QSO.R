#' Quasi-Sequence-Order (QSO) Descriptor
#'
#' Quasi-Sequence-Order (QSO) Descriptor
#'
#' This function calculates the Quasi-Sequence-Order (QSO) descriptor
#' (Dim: \code{20 + 20 + (2 * nlag)}, default is 100).
#'
#' @param x A character vector, as the input protein sequence.
#'
#' @param nlag The maximum lag, defualt is 30.
#'
#' @param w The weighting factor, default is 0.1.
#'
#' @return A length \code{20 + 20 + (2 * nlag)} named vector
#'
#' @keywords extract QSO
#'
#' @aliases extractQSO
#'
#' @author Nan Xiao <\url{http://nanx.me}>
#'
#' @seealso See \code{\link{extractSOCN}} for sequence-order-coupling numbers.
#'
#' @export extractQSO
#'
#' @references
#' Kuo-Chen Chou. Prediction of Protein Subcellar Locations by
#' Incorporating Quasi-Sequence-Order Effect.
#' \emph{Biochemical and Biophysical Research Communications},
#' 2000, 278, 477-483.
#'
#' Kuo-Chen Chou and Yu-Dong Cai. Prediction of Protein Sucellular Locations by
#' GO-FunD-PseAA Predictor.
#' \emph{Biochemical and Biophysical Research Communications},
#' 2004, 320, 1236-1239.
#'
#' Gisbert Schneider and Paul Wrede. The Rational Design of
#' Amino Acid Sequences by Artifical Neural Networks and Simulated
#' Molecular Evolution: Do Novo Design of an Idealized Leader Cleavge Site.
#' \emph{Biophys Journal}, 1994, 66, 335-344.
#'
#' @examples
#' x = readFASTA(system.file('protseq/P00750.fasta', package = 'protr'))[[1]]
#' extractQSO(x)

extractQSO = function (x, nlag = 30, w = 0.1) {

  if (protcheck(x) == FALSE) stop('x has unrecognized amino acid type')

  N = nchar(x)
  if (N <= nlag) stop('Length of the protein sequence must be greater than "nlag"')

  DistMat1 = read.csv(system.file('sysdata/Schneider-Wrede.csv', package = 'protr'), header = TRUE)
  DistMat2 = read.csv(system.file('sysdata/Grantham.csv', package = 'protr'), header = TRUE)
  row.names(DistMat1) = as.character(DistMat1[, 1])
  DistMat1 = DistMat1[, -1]
  row.names(DistMat2) = as.character(DistMat2[, 1])
  DistMat2 = DistMat2[, -1]

  xSplitted = strsplit(x, split = '')[[1]]

  # Compute Schneider.tau_d

  tau1 = vector('list', nlag)

  for (d in 1:nlag) {
    for (i in 1:(N - d)) {
      tau1[[d]][i] = (DistMat1[xSplitted[i], xSplitted[i + d]])^2
    }
  }

  tau1 = sapply(tau1, sum)

  # Compute Grantham.tau_d

  tau2 = vector('list', nlag)

  for (d in 1:nlag) {
    for (i in 1:(N - d)) {
      tau2[[d]][i] = (DistMat2[xSplitted[i], xSplitted[i + d]])^2
    }
  }

  tau2 = sapply(tau2, sum)

  # Compute fr

  AADict = c('A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I',
             'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V')

  fr = summary(factor(xSplitted, levels = AADict), maxsum = 21)

  # Compute Schneider.Xr Grantham.Xr Schneider.Xd Grantham.Xd

  Xr1 = fr/(1 + (w * sum(tau1)))
  names(Xr1) = paste('Schneider.Xr.', names(Xr1), sep = '')

  Xr2 = fr/(1 + (w * sum(tau2)))
  names(Xr2) = paste('Grantham.Xr.', names(Xr2), sep = '')

  Xd1 = (w * tau1)/(1 + (w * sum(tau1)))
  names(Xd1) = paste('Schneider.Xd.', 1:nlag, sep = '')

  Xd2 = (w * tau2)/(1 + (w * sum(tau2)))
  names(Xd2) = paste('Grantham.Xd.', 1:nlag, sep = '')

  QSO = c(Xr1, Xr2, Xd1, Xd2)

  return(QSO)

}
