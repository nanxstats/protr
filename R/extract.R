#' Feature Extraction for Protein Sequences
#'
#' Feature Extraction for Protein Sequences
#' 
#' This function takes a protein sequence and descriptor type 
#' as arguments and returns the result as a data frame.
#' 
#' @param x      A character vector, as the input protein sequence. 
#'               Length \code{n}.
#' 
#' @param type   a vector of size n * 1 describing the chunklets:
#'               -1 in the i'th place says that point i doesn't belong to any chunklet;
#'               integer j in place i says that point i belongs to chunklet j.
#'               The chunklets indexes should be 1:(number of chunklets).
#' @param \dots  Additional options used by specific descriptor 
#'               as described below.
#' 
#' @return The result data frame
#'
#' The three returned argument are just different forms of the same output.
#' If one is interested in a Mahalanobis metric over the original data space, 
#' the first argument is all she/he needs. If a transformation into another
#' space (where one can use the Euclidean metric) is preferred, the second
#' returned argument is sufficient. Using A and B is equivalent in the 
#' following sense:
#' 
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
#' \item{\code{'QSOD'}}{Quasi-sequence Order Descriptors (Quasi-sequence Order Descriptor)}
#' \item{\code{'PAAC'}}{Pseudo Amino Acid Composition (Pseudo Amino Acid Composition Descriptor)}
#' \item{\code{'APAAC'}}{Pseudo Amino Acid Composition - Amphiphilic Pseudo Amino Acid Composition}}
#' 
#' @keywords rdpi extract feature extraction descriptors
#'
#' @aliases extract
#' 
#' @note Note that any different sets of instances (chunklets),
#'       e.g. {1, 3, 7} and {4, 6}, might belong to the 
#'       same class and might belong to different classes.
#' 
#' @author Xiao Nan <\url{http://www.road2stat.com}>
#' 
#' @seealso See \code{\link{extractAAC}} for Amino Acid Composition, etc.
#' 
#' @export extract
#' 
#' @references
#' Aharon Bar-Hillel, Tomer Hertz, Noam Shental, and Daphna Weinshall (2003).
#' Learning Distance Functions using Equivalence Relations.
#' \emph{Proceedings of 20th International Conference on
#' Machine Learning (ICML2003)}.
#' 
#' @examples
#' A06852 = readFASTA(system.file('AAseq/A06852.fasta', package = 'rdpi'))
#' # extract(A06852, 'AAC')
#' 

extract = function (x, type, ...) {

  type = c('AAC', 'DC', 'TC', 
           'Moreau-Broto', 'Moran', 'Geary', 
           'CTD-C', 'CTD-T', 'CTD-D',
           'C-Triad', 
           'SOCN', 'QSOD', 
           'PAAC', 'APAAC')
  
}
