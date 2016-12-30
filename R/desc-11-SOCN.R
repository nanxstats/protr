#' Sequence-Order-Coupling Numbers
#'
#' Sequence-Order-Coupling Numbers
#'
#' This function calculates the Sequence-Order-Coupling Numbers
#' (Dim: \code{nlag * 2}, default is 60).
#'
#' @param x A character vector, as the input protein sequence.
#'
#' @param nlag The maximum lag, defualt is 30.
#'
#' @return A length \code{nlag * 2} named vector
#'
#' @keywords extract SOCN Order Coupling
#'
#' @aliases extractSOCN
#'
#' @author Nan Xiao <\url{http://nanx.me}>
#'
#' @seealso See \code{\link{extractQSO}} for quasi-sequence-order descriptors.
#'
#' @export extractSOCN
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
#' extractSOCN(x)

extractSOCN = function (x, nlag = 30) {

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

  names(tau1) = paste('Schneider.lag', 1:nlag, sep = '')

  # Compute Grantham.tau_d

  tau2 = vector('list', nlag)

  for (d in 1:nlag) {
    for (i in 1:(N - d)) {
      tau2[[d]][i] = (DistMat2[xSplitted[i], xSplitted[i + d]])^2
    }
  }

  tau2 = sapply(tau2, sum)

  names(tau2) = paste('Grantham.lag', 1:nlag, sep = '')

  tau = c(tau1, tau2)

  return(tau)

}
