#' Feature Extraction for Protein Sequences
#'
#' Feature Extraction for Protein Sequences
#' 
#' This function takes a protein sequence and descriptor type(s)
#' as arguments and returns the calculated descriptors.
#' This is basically a wrapper function for the \code{extractXXX()}
#' functions in this package.
#' 
#' @param x      A character vector, as the input protein sequence.
#' 
#' @param type   a character vector specifying the descriptor types
#'               to be calculated. See below for details.
#'               
#' @param \dots  Additional options used by specific descriptor 
#'               as described below.
#' 
#' @return The result named vector
#'
#' The elements that could be used in the \code{type} argument:
#' \describe{
#' \item{\code{'AAC'}}{Amino Acid Composition (Amino Acid Composition Descriptor)}
#' \item{\code{'DC'}}{Dipeptide Composition (Amino Acid Composition Descriptor)}
#' \item{\code{'TC'}}{Tripeptide Composition (Amino Acid Composition Descriptor)}
#' \item{\code{'Moreau-Broto'}}{Normalized Moreau-Broto Autocorrelation (Autocorrelation Descriptor)}
#' \item{\code{'Moran'}}{Moran Autocorrelation (Autocorrelation Descriptor)}
#' \item{\code{'Geary'}}{Geary Autocorrelation (Autocorrelation Descriptor)}
#' \item{\code{'CTD-C'}}{Composition (CTD Descriptor)}
#' \item{\code{'CTD-T'}}{Transition (CTD Descriptor)}
#' \item{\code{'CTD-D'}}{Distribution (CTD Descriptor)}
#' \item{\code{'C-Triad'}}{Conjoint Triad}
#' \item{\code{'SOCN'}}{Sequence Order Coupling Number (Quasi-sequence Order Descriptor)}
#' \item{\code{'QSO'}}{Quasi-sequence Order Descriptors (Quasi-sequence Order Descriptor)}
#' \item{\code{'PAAC'}}{Pseudo Amino Acid Composition (Pseudo Amino Acid Composition Descriptor)}
#' \item{\code{'APAAC'}}{Pseudo Amino Acid Composition - Amphiphilic Pseudo Amino Acid Composition}}
#' 
#' @keywords rdpi extract feature extraction descriptor
#'
#' @aliases extract
#' 
#' @author Xiao Nan <\url{http://www.road2stat.com}>
#' 
#' @seealso See \code{\link{extractAAC}} for Amino Acid Composition, etc.
#' 
#' @export extract
#' 
#' @examples
#' # Load protein sequences
#' # prot  = readFASTA(system.file('AAseq/P00750.fasta', package = 'rdpi'))[[1]]
#' # prots = 
#' # extract single feature type of single protein sequence
#' # extract(prot, 'AAC')
#' # extract multiple feature type of single protein sequence
#' # extract(prot, type = c('AAC', 'PAAC', 'APAAC'))
#' # extract single feature type of multiple protein sequence
#' # lapply(prots, extract, 'AAC')
#' # extract multiple feature type of multiple protein sequences
#' # lapply(prots, extract, type = c('AAC', 'PAAC', 'APAAC'))
#' 

extract = function (x, type, ...) {

  type = c('AAC', 'DC', 'TC', 
           'Moreau-Broto', 'Moran', 'Geary', 
           'CTD-C', 'CTD-T', 'CTD-D',
           'C-Triad', 
           'SOCN', 'QSO', 
           'PAAC', 'APAAC')
  
}

