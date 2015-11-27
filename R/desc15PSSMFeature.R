#' Profile-based protein representation derived by PSSM
#' (Position-Specific Scoring Matrix)
#'
#' Profile-based protein representation derived by PSSM
#' (Position-Specific Scoring Matrix)
#'
#' This function calculates the profile-based protein representation
#' derived by PSSM. The feature vector is based on the PSSM computed by
#' \code{\link{extractPSSM}}. For a given sequence,
#' The PSSM feature represents the log-likelihood of the substitution of the
#' 20 types of amino acids at that position in the sequence.
#' Each PSSM feature value in the vector represents the degree of conservation
#' of a given amino acid type. The value is normalized to
#' interval (0, 1) by the transformation 1/(1+e^(-x)).
#'
#' @param pssmmat The PSSM computed by \code{\link{extractPSSM}}.
#'
#' @return A numeric vector which has \code{20 x N} named elements,
#' where \code{N} is the size of the window (number of rows of the PSSM).
#'
#' @seealso \link{extractPSSM} \link{extractPSSMAcc}
#'
#' @keywords extract PSSM Blast Alignment
#'
#' @aliases extractPSSMFeature
#'
#' @author Nan Xiao <\url{http://nanx.me}>
#'
#' @export extractPSSMFeature
#'
#' @references
#' Ye, Xugang, Guoli Wang, and Stephen F. Altschul.
#' "An assessment of substitution scores for protein profile-profile comparison."
#' \emph{Bioinformatics} 27.24 (2011): 3356--3363.
#'
#' Rangwala, Huzefa, and George Karypis.
#' "Profile-based direct kernels for remote homology detection and fold recognition."
#' \emph{Bioinformatics} 21.23 (2005): 4239--4247.
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
#'   pssmfeature = extractPSSMFeature(pssmmat)
#'   head(pssmfeature)
#' }

extractPSSMFeature = function(pssmmat) {

  # Normalize PSSM scores to (0, 1)
  result = as.vector(1/(1 + exp(pssmmat)))

  AADict = c('A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I',
             'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V')
  names(result) = paste0('PSSM_', paste0(rep(1L:ncol(pssmmat), each = 20L),
                                         '_', AADict))

  return(result)

}
