#' Profile-based protein representation derived by PSSM
#' (Position-Specific Scoring Matrix) and auto cross covariance
#'
#' Profile-based protein representation derived by PSSM
#' (Position-Specific Scoring Matrix) and auto cross covariance
#'
#' This function calculates the feature vector based on the PSSM
#' by running PSI-Blast and auto cross covariance tranformation.
#'
#' @param pssmmat The PSSM computed by \code{\link{extractPSSM}}.
#' @param lag The lag parameter. Must be less than the number of amino acids
#' in the sequence (i.e. the number of columns in the PSSM matrix).
#'
#' @return A length \code{lag * 20^2} named numeric vector,
#' the element names are derived by the amino acid name abbreviation
#' (crossed amino acid name abbreviation) and lag index.
#'
#' @seealso \link{extractPSSM} \link{extractPSSMFeature}
#'
#' @keywords extract PSSM Blast Alignment
#'
#' @aliases extractPSSMAcc
#'
#' @author Nan Xiao <\url{http://nanx.me}>
#'
#' @export extractPSSMAcc
#'
#' @references
#' Wold, S., Jonsson, J., Sj\"{o}rstr\"{o}m, M., Sandberg,
#' M., & R\"{a}nnar, S. (1993).
#' DNA and peptide sequences and chemical processes multivariately modelled
#' by principal component analysis and partial least-squares projections
#' to latent structures.
#' \emph{Analytica chimica acta}, 277(2), 239--253.
#'
#' @examples
#' if (Sys.which('makeblastdb') == '' | Sys.which('psiblast') == '') {
#'   cat('Could not find makeblastdb or psiblast. Please install NCBI Blast+ first.')
#' } else {
#'   x = readFASTA(system.file('protseq/P00750.fasta', package = 'protr'))[[1]]
#'   dbpath = tempfile('tempdb', fileext = '.fasta')
#'   invisible(file.copy(from = system.file('protseq/Plasminogen.fasta',
#'                                          package = 'protr'), to = dbpath))
#'   pssmmat = extractPSSM(seq = x, database.path = dbpath)
#'   pssmacc = extractPSSMAcc(pssmmat, lag = 3)
#'   tail(pssmacc)
#' }

extractPSSMAcc = function(pssmmat, lag) {

  # Normalize PSSM scores to (0, 1)
  mat = 1/(1 + exp(pssmmat))

  accpssm = function (mat, lag) {

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

    AADict = c('A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I',
               'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V')
    names(acc1) = as.vector(outer(AADict, paste0('lag', 1:lag),
                                  paste, sep = '.'))

    # auto cross covariance: p^2 - p elements

    acc2 = matrix(0, nrow = p^2 - p, ncol = lag)

    idx = cbind(combn(1:p, 2), combn(1:p, 2)[2:1, ])

    for (j in 1:ncol(idx)) {
      for (i in 1:lag) {
        acc2[j, i] = sum(mat[idx[1, j], 1:(n - i)] * mat[idx[2, j], ((1:(n - i)) + i)]) / (n - i)
      }
    }

    acc2 = as.vector(acc2)
    names(acc2) = as.vector(outer(paste(AADict[idx[1, ]], AADict[idx[2, ]],
                                        sep = '.'),
                                  paste0('lag', 1:lag), paste, sep = '.'))

    ACC = c(acc1, acc2)

    return(ACC)

  }

  result = accpssm(mat, lag)

  return(result)

}
