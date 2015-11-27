#' Auto Cross Covariance (ACC) for Generating Scales-Based Descriptors of the Same Length
#'
#' Auto Cross Covariance (ACC) for Generating Scales-Based Descriptors of the Same Length
#'
#' This function calculates the auto covariance and auto cross covariance
#' for generating scale-based descriptors of the same length.
#'
#' @param mat A \code{p * n} matrix. Each row represents one scale
#' (total \code{p} scales), each column represents one amino acid position
#' (total \code{n} amino acids).
#'
#' @param lag The lag parameter. Must be less than the amino acids.
#'
#' @return A length \code{lag * p^2} named vector, the element names are
#'         constructed by: the scales index (crossed scales index) and lag index.
#'
#' @note To know more details about auto cross covariance, see the references.
#'
#' @keywords acc covariance
#'
#' @aliases acc
#'
#' @author Nan Xiao <\url{http://nanx.me}>
#'
#' @seealso See \code{\link{extractScales}} for scales-based descriptors.
#' For more details, see \code{\link{extractDescScales}}
#' and \code{\link{extractProtFP}}.
#'
#' @export acc
#'
#' @references
#' Wold, S., Jonsson, J., Sj\"{o}rstr\"{o}m, M., Sandberg, M., & R\"{a}nnar, S. (1993).
#' DNA and peptide sequences and chemical processes multivariately modelled
#' by principal component analysis and partial least-squares projections
#' to latent structures.
#' \emph{Analytica chimica acta}, 277(2), 239--253.
#'
#' Sj\"{o}str\"{o}m, M., R\"{a}nnar, S., & Wieslander, A. (1995).
#' Polypeptide sequence property relationships in \emph{Escherichia coli}
#' based on auto cross covariances.
#' \emph{Chemometrics and intelligent laboratory systems}, 29(2), 295--305.
#'
#' @examples
#' p = 8    # p is the scales number
#' n = 200  # n is the amino acid number
#' lag = 7  # the lag paramter
#' mat = matrix(rnorm(p * n), nrow = p, ncol = n)
#' acc(mat, lag)

acc = function (mat, lag) {

  p = nrow(mat)
  n = ncol(mat)

  if (lag > n) stop('lag must be smaller than the amino acids in the sequence')

  # auto covariance: p elements

  acc1 = matrix(0, nrow = p, ncol = lag)

  for (j in 1:p) {
    for (i in 1:lag) {
      acc1[j, i] = sum(mat[j, 1:(n - i)] * mat[j, ((1:(n - i)) + i)]) / (n - i)
    }
  }

  acc1 = as.vector(acc1)
  names(acc1) = as.vector(outer(paste0('scl', 1:p), paste0('lag', 1:lag), paste, sep = '.'))

  # auto cross covariance: p^2 - p elements

  acc2 = matrix(0, nrow = p^2 - p, ncol = lag)

  idx = cbind(combn(1:p, 2), combn(1:p, 2)[2:1, ])

  for (j in 1:ncol(idx)) {
    for (i in 1:lag) {
      acc2[j, i] = sum(mat[idx[1, j], 1:(n - i)] * mat[idx[2, j], ((1:(n - i)) + i)]) / (n - i)
    }
  }

  acc2 = as.vector(acc2)
  names(acc2) = as.vector(outer(paste0('scl', paste(idx[1, ], idx[2, ], sep = '.')), paste0('lag', 1:lag), paste, sep = '.'))

  ACC = c(acc1, acc2)

  return(ACC)

}
