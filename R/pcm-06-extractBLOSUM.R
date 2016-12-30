#' BLOSUM and PAM Matrix-Derived Descriptors
#'
#' BLOSUM and PAM Matrix-Derived Descriptors
#'
#' This function calculates BLOSUM matrix-derived descriptors.
#' For users' convenience, \code{protr} provides the
#' BLOSUM45, BLOSUM50, BLOSUM62, BLOSUM80, BLOSUM100,
#' PAM30, PAM40, PAM70, PAM120, and PAM250 matrices
#' for the 20 amino acids to select.
#'
#' @param x A character vector, as the input protein sequence.
#' @param submat Substitution matrix for the 20 amino acids. Should be one of
#'        \code{AABLOSUM45}, \code{AABLOSUM50}, \code{AABLOSUM62},
#'        \code{AABLOSUM80}, \code{AABLOSUM100}, \code{AAPAM30},
#'        \code{AAPAM40}, \code{AAPAM70}, \code{AAPAM120}, \code{AAPAM250}.
#'        Default is \code{'AABLOSUM62'}.
#' @param k Integer. The number of selected scales (i.e. the first \code{k} scales)
#'        derived by the substitution matrix.
#'        This could be selected according to the printed relative importance values.
#' @param lag The lag parameter. Must be less than the amino acids.
#' @param scale Logical. Should we auto-scale the substitution matrix
#'        (\code{submat}) before doing eigen decomposition? Default is \code{TRUE}.
#' @param silent Logical. Whether we print the relative importance of each scales
#'        (diagnal value of the eigen decomposition result matrix B) or not.
#'        Default is \code{TRUE}.
#' @return  A length \code{lag * p^2} named vector,
#'         \code{p} is the number of scales selected.
#'
#' @keywords extract BLOSUM PAM PCM
#'
#' @aliases extractBLOSUM
#'
#' @author Nan Xiao <\url{http://nanx.me}>
#'
#' @export extractBLOSUM
#'
#' @references
#' Georgiev, A. G. (2009).
#' Interpretable numerical descriptors of amino acid space.
#' Journal of Computational Biology, 16(5), 703--723.
#'
#' @examples
#' x = readFASTA(system.file('protseq/P00750.fasta', package = 'protr'))[[1]]
#' blosum = extractBLOSUM(x, submat = 'AABLOSUM62', k = 5, lag = 7, scale = TRUE, silent = FALSE)

extractBLOSUM = function (x, submat = 'AABLOSUM62', k, lag, scale = TRUE, silent = TRUE) {

  if (protcheck(x) == FALSE) stop('x has unrecognized amino acid type')

  k = min(k, 20)

  submat = get(submat)
  if (scale) submat = scale(submat)

  eig = eigen(submat)
  A = eig$vectors
  B = eig$values
  rownames(A) = rownames(submat)
  # the equation: submat == A %*% diag(B) %*% t(A)

  accmat = matrix(0, k, nchar(x))
  x.split = strsplit(x, '')[[1]]

  for (i in 1:nchar(x)) {
    accmat[, i] = A[x.split[i], 1:k]
  }

  result = acc(accmat, lag)

  if (!silent) {
    cat('Relative importance of all the possible 20 scales:\n')
    print(B)
  }

  return(result)

}
